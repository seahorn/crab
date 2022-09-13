/*******************************************************************************
 * Array smashing domain
 *
 * Word-level assumption: all array accesses access the same number of
 * bytes. This is important since this domain always smashes an array
 * so all accesses should have same size.  The domain checks for that.
 *
 * This domain also assumes that no region or references are passed to
 * it but it doesn't check for it.
 *
 * Ghost variables:
 *
 * The array smashing domain creates a ghost scalar variable for each
 * smashed array. Very importantly, if the array smashing domain is
 * called by another domain D with ghost variables created by D then
 * these ghost variables cannot change during the analysis. Otherwise,
 * the array smashing domain might be unsound.
 *
 * On the other hand, the array smashing domain DOES NOT BREAK
 * MODULARITY of the base domain. That is, the base domain can create
 * ghost variables based on the ghost variables created by the array
 * smashing domain.
 *
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/constant.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {

namespace domains {

// All arrays are `smashed` into a single summarized variable.
template <typename BaseNumDomain>
class array_smashing final
    : public abstract_domain_api<array_smashing<BaseNumDomain>> {

public:
  using number_t = typename BaseNumDomain::number_t;
  using varname_t = typename BaseNumDomain::varname_t;

private:
  using array_smashing_t = array_smashing<BaseNumDomain>;
  using abstract_domain_t = abstract_domain_api<array_smashing_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using base_dom_t = BaseNumDomain;
  using interval_t = ikos::interval<number_t>;

private:
  using variable_factory_t = typename varname_t::variable_factory_t;
  using bytes_t = constant<uint64_t>;
  using last_access_env_t = ikos::separate_domain<variable_t, bytes_t>;

  // map an array to the number of *bytes* last being accessed
  last_access_env_t m_last_access_env;
  // the base domain
  base_dom_t m_base_dom;

  // set last access to array_var to sz bytes
  void set_size(const variable_t &array_var, uint64_t sz) {
    m_last_access_env.set(array_var, sz);
  }

  // get the size of the last access to array_var
  bytes_t get_size(const variable_t &array_var) {
    return m_last_access_env.at(array_var);
  }

  // return true if the last access to array_var is sz
  bool equal_size(const variable_t &array_var, uint64_t sz) {
    if (const bytes_t *val = m_last_access_env.find(array_var)) {
      return (val->is_constant() && (val->get_constant() == sz));
    }
    return false;
  }

  variable_t mk_scalar_var(const variable_t &array_var,
                           uint64_t size /*bytes*/) {
    auto array_ty = array_var.get_type();
    assert(array_ty.is_bool_array() || array_ty.is_integer_array() ||
           array_ty.is_real_array());

    variable_type ghost_ty(INT_TYPE, size * 8);
    if (array_ty.is_bool_array()) {
      ghost_ty = variable_type(BOOL_TYPE);
    } else if (array_ty.is_real_array()) {
      ghost_ty = variable_type(REAL_TYPE);
    }

    // The variable factory does caching so it returns always the same
    // variable for the same array variable.
    auto &vfac =
        const_cast<varname_t *>(&(array_var.name()))->get_var_factory();
    return variable_t(vfac.get(array_var.name(), ".smashed"), ghost_ty);
  }

  // Convert elem_size to an uint64_t
  uint64_t check_and_get_elem_size(const linear_expression_t &elem_size) const {
    auto to_interval = [this](const linear_expression_t &e) {
      interval_t r(e.constant());
      for (auto kv : e) {
        interval_t c(kv.first);
        r += c * m_base_dom.at(kv.second);
      }
      return r;
    };
    interval_t i_elem_size = to_interval(elem_size);
    if (boost::optional<number_t> n_bytes = i_elem_size.singleton()) {
      if (static_cast<int64_t>(*n_bytes) > 0) {
        return (uint64_t) static_cast<int64_t>(*n_bytes);
      }
    }
    CRAB_ERROR("array domains expect constant array element sizes ",
               "between 1 and ", std::numeric_limits<uint64_t>::max(),
               ". Found ", elem_size);
  }

  // The internal representation contains ghost variables to represent
  // the contents of arrays. They shouldn't be exposed outside via
  // linear constraints.
  linear_constraint_system_t
  filter_ghost_vars(linear_constraint_system_t &&csts) const {
    linear_constraint_system_t res;
    for (auto const &cst : csts) {
      if (std::all_of(
              cst.expression().variables_begin(),
              cst.expression().variables_end(), [](const variable_t &v) {
                return (v.name().str().find(".smashed") == std::string::npos);
              })) {
        res += cst;
      }
    }
    return res;
  }

  // lhs and rhs only involve scalar variables
  void do_update(base_dom_t &dom, const variable_t &lhs,
		 const linear_expression_t &rhs, bool weak) {
    auto ty = lhs.get_type();
    if (ty.is_bool()) {
      if (rhs.is_constant()) {
        if (rhs.constant() >= number_t(1)) {
	  if (!weak) {
	    dom.assign_bool_cst(lhs, linear_constraint_t::get_true());
	  } else {
	    dom.weak_assign_bool_cst(lhs, linear_constraint_t::get_true());
	  }
        } else {
	  if (!weak) {
	    dom.assign_bool_cst(lhs, linear_constraint_t::get_false());
	  } else {
	    dom.weak_assign_bool_cst(lhs, linear_constraint_t::get_false());
	  }
        }
      } else if (auto rhs_v = rhs.get_variable()) {
	if (!weak) {
	  dom.assign_bool_var(lhs, (*rhs_v), false);
	} else {
	  dom.weak_assign_bool_var(lhs, (*rhs_v), false);
	}
      }
    } else {
      assert(ty.is_integer() || ty.is_real());
      if (!weak) {
	dom.assign(lhs, rhs);
      } else {
	dom.weak_assign(lhs, rhs);
      }
    }
  }

  void do_strong_update(const variable_t &lhs, const linear_expression_t &rhs) {
    do_update(m_base_dom, lhs, rhs, false /*!weak*/);
  }

  void do_weak_update(const variable_t &lhs, const linear_expression_t &rhs) {
    do_update(m_base_dom, lhs, rhs, true /*weak*/);
  }

  array_smashing(last_access_env_t &&last_access_env, base_dom_t &&base_dom)
      : m_last_access_env(std::move(last_access_env)),
        m_base_dom(std::move(base_dom)) {}

public:
  array_smashing() {}

  array_smashing make_top() const override { return array_smashing(); }

  array_smashing make_bottom() const override {
    array_smashing res;
    res.set_to_bottom();
    return res;
  }

  void set_to_top() override {
    m_last_access_env = last_access_env_t::top();
    m_base_dom.set_to_top();
  }

  void set_to_bottom() override {
    m_last_access_env.set_to_bottom();
    m_base_dom.set_to_bottom();
  }

  array_smashing(const array_smashing_t &other)
      : m_last_access_env(other.m_last_access_env),
        m_base_dom(other.m_base_dom) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  array_smashing(array_smashing_t &&other) = default;

  array_smashing_t &operator=(const array_smashing_t &other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      m_last_access_env = other.m_last_access_env;
      m_base_dom = other.m_base_dom;
    }
    return *this;
  }

  array_smashing_t &operator=(array_smashing_t &&other) = default;

  bool is_bottom() const override { return m_base_dom.is_bottom(); }

  bool is_top() const override { return m_base_dom.is_top(); }

  bool operator<=(const array_smashing_t &other) const override {
    return (m_base_dom <= other.m_base_dom);
  }

  void operator|=(const array_smashing_t &other) override {
    m_last_access_env = m_last_access_env | other.m_last_access_env;
    m_base_dom |= other.m_base_dom;
  }

  array_smashing_t operator|(const array_smashing_t &other) const override {
    return array_smashing_t(m_last_access_env | other.m_last_access_env,
                            m_base_dom | other.m_base_dom);
  }

  void operator&=(const array_smashing_t &other) override {
    m_last_access_env = m_last_access_env & other.m_last_access_env;
    m_base_dom &= other.m_base_dom;
  }
  
  array_smashing_t operator&(const array_smashing_t &other) const override {
    return array_smashing_t(m_last_access_env & other.m_last_access_env,
                            m_base_dom & other.m_base_dom);
  }

  array_smashing_t operator||(const array_smashing_t &other) const override {
    return array_smashing_t(m_last_access_env || other.m_last_access_env,
                            m_base_dom || other.m_base_dom);
  }

  array_smashing_t
  widening_thresholds(const array_smashing_t &other,
                      const thresholds<number_t> &ts) const override {
    return array_smashing_t(
        m_last_access_env.widening_thresholds(other.m_last_access_env, ts),
        m_base_dom.widening_thresholds(other.m_base_dom, ts));
  }

  array_smashing_t operator&&(const array_smashing_t &other) const override {
    return array_smashing_t(m_last_access_env && other.m_last_access_env,
                            m_base_dom && other.m_base_dom);
  }

  virtual interval_t operator[](const variable_t &v) override {
    return m_base_dom[v];
  }

  virtual interval_t at(const variable_t &v) const override {
    return m_base_dom.at(v);
  }

  void forget(const variable_vector_t &variables) override {
    variable_vector_t remove_vars;
    for (auto v : variables) {
      if (v.get_type().is_array()) {
        bytes_t size = get_size(v);
        if (size.is_constant()) {
          remove_vars.push_back(mk_scalar_var(v, size.get_constant()));
          m_last_access_env -= v;
        }
      } else {
        remove_vars.push_back(v);
      }
    }
    // No array variables should go to m_base_dom
    m_base_dom.forget(remove_vars);
  }

  void project(const variable_vector_t &variables) override {
    variable_vector_t keep_vars, keep_array_vars;
    for (auto v : variables) {
      if (v.get_type().is_array()) {
        bytes_t size = get_size(v);
        if (size.is_constant()) {
          variable_t gvar(mk_scalar_var(v, size.get_constant()));
          keep_vars.push_back(gvar);
          keep_array_vars.push_back(v);
        }
      } else {
        keep_vars.push_back(v);
      }
    }
    // No array variables should go to m_base_dom
    m_base_dom.project(keep_vars);
    m_last_access_env.project(keep_array_vars);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (var.get_type() != new_var.get_type()) {
      CRAB_ERROR(domain_name(), "::expand must preserve the same type");
    }

    if (var.get_type().is_array()) {
      bytes_t size = get_size(var);
      if (size.is_constant()) {
        variable_t gvar(mk_scalar_var(var, size.get_constant()));
        variable_t gnew_var(mk_scalar_var(new_var, size.get_constant()));
        // No array variables should go to m_base_dom
        m_base_dom.expand(gvar, gnew_var);
        set_size(new_var, size.get_constant());
      }
    } else {
      m_base_dom.expand(var, new_var);
    }
  }

  void normalize() override { m_base_dom.normalize(); }

  void minimize() override { m_base_dom.minimize(); }

  void operator+=(const linear_constraint_system_t &csts) override {
    m_base_dom += csts;
  }

  bool entails(const linear_constraint_t &cst) const override {
    return m_base_dom.entails(cst);
  }
  
  void operator-=(const variable_t &var) override {
    if (var.get_type().is_array()) {
      bytes_t size = get_size(var);
      if (size.is_constant()) {
        variable_t gvar(mk_scalar_var(var, size.get_constant()));
        // No array variable should go to m_base_dom
        m_base_dom -= gvar;
        m_last_access_env -= var;
      }
    } else {
      m_base_dom -= var;
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_base_dom.assign(x, e);

    CRAB_LOG("smashing", crab::outs()
                             << "apply " << x << " := " << e << *this << "\n";);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    m_base_dom.weak_assign(x, e);

    CRAB_LOG("smashing", crab::outs()
	     << "weak_assign(" << x << "," << e << ")" <<  *this << "\n";);
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    m_base_dom.select(lhs, cond, e1, e2);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const array_smashing_t &inv) override {
    m_base_dom.backward_assign(x, e, inv.m_base_dom);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const array_smashing_t &inv) override {
    m_base_dom.backward_apply(op, x, y, z, inv.m_base_dom);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const array_smashing_t &inv) override {
    m_base_dom.backward_apply(op, x, y, z, inv.m_base_dom);
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    m_base_dom.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_base_dom.apply(op, x, y, k);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << k << *this << "\n";);
  }

  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    m_base_dom.assign_bool_cst(lhs, rhs);
  }

  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) override {
    m_base_dom.assign_bool_ref_cst(lhs, rhs);
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    m_base_dom.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  virtual void weak_assign_bool_cst(const variable_t &lhs,
				    const linear_constraint_t &rhs) override {
    m_base_dom.weak_assign_bool_cst(lhs, rhs);
  }

  virtual void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
				    bool is_not_rhs) override {
    m_base_dom.weak_assign_bool_var(lhs, rhs, is_not_rhs);
  }
  
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    m_base_dom.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    m_base_dom.assume_bool(v, is_negated);
  }

  virtual void select_bool(const variable_t &lhs, const variable_t &cond,
                           const variable_t &b1,
                           const variable_t &b2) override {
    m_base_dom.select_bool(lhs, cond, b1, b2);
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const array_smashing_t &inv) override {
    m_base_dom.backward_assign_bool_cst(lhs, rhs, inv.m_base_dom);
  }

  virtual void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const array_smashing_t &inv) override {
    m_base_dom.backward_assign_bool_ref_cst(lhs, rhs, inv.m_base_dom);
  }

  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const array_smashing_t &inv) override {
    m_base_dom.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv.m_base_dom);
  }

  virtual void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const array_smashing_t &inv) override {
    m_base_dom.backward_apply_binary_bool(op, x, y, z, inv.m_base_dom);
  }

  // array_operators_api

  // All the array elements are initialized to val
  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size /*bytes*/,
                          const linear_expression_t & /*lb_idx*/,
                          const linear_expression_t & /*ub_idx*/,
                          const linear_expression_t &val) override {

    uint64_t size = check_and_get_elem_size(elem_size);
    set_size(a, size);
    variable_t scalar_var(mk_scalar_var(a, size));

    auto ty = scalar_var.get_type();
    if (ty.is_bool()) {
      if (val.is_constant()) {
        if (val.constant() >= number_t(1)) {
          m_base_dom.assign_bool_cst(scalar_var,
                                     linear_constraint_t::get_true());
        } else {
          m_base_dom.assign_bool_cst(scalar_var,
                                     linear_constraint_t::get_false());
        }
      } else if (auto var = val.get_variable()) {
        m_base_dom.assign_bool_var(scalar_var, (*var), false);
      }
    } else {
      assert(ty.is_integer() || ty.is_real());
      m_base_dom.assign(scalar_var, val);
    }
    CRAB_LOG("smashing", crab::outs() << "forall i:: " << a << "[i]==" << val
                                      << " -- " << *this << "\n";);
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size /*bytes*/,
                          const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.load");
    crab::ScopedCrabStats __st__(domain_name() + ".load");

    uint64_t size = check_and_get_elem_size(elem_size);
    if (equal_size(a, size)) {
      variable_t scalar_var(mk_scalar_var(a, size));

      // We need to be careful when assigning a summarized variable a
      // into a non-summarized variable lhs. Simply m_base_dom.assign(lhs,a)
      // is not sound.
      auto &vfac = const_cast<varname_t *>(&(a.name()))->get_var_factory();
      variable_t copy_scalar_var(vfac.get(scalar_var.name(), ".copy"),
                                 scalar_var.get_type());
      m_base_dom.expand(scalar_var, copy_scalar_var);
      auto ty = scalar_var.get_type();
      if (ty.is_bool()) {
        m_base_dom.assign_bool_var(lhs, copy_scalar_var, false);
      } else {
        assert(ty.is_integer() || ty.is_real());
        m_base_dom.assign(lhs, copy_scalar_var);
      }
      m_base_dom -= copy_scalar_var;
    } else {
      m_base_dom -= lhs;
    }

    CRAB_LOG("smashing", crab::outs() << lhs << ":=" << a << "[" << i
                                      << "]  -- " << *this << "\n";);
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.store");
    crab::ScopedCrabStats __st__(domain_name() + ".store");

    uint64_t size = check_and_get_elem_size(elem_size);
    if (is_strong_update) {
      set_size(a, size);
    }

    if (equal_size(a, size)) {
      variable_t scalar_var(mk_scalar_var(a, size));
      if (is_strong_update) {
        do_strong_update(scalar_var, val);
      } else {
        do_weak_update(scalar_var, val);
      }
    }

    CRAB_LOG("smashing", crab::outs() << a << "[" << i << "]:=" << val << " -- "
                                      << *this << "\n";);
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.store");
    crab::ScopedCrabStats __st__(domain_name() + ".store");

    uint64_t size = check_and_get_elem_size(elem_size);
    if (equal_size(a, size)) {
      variable_t scalar_var(mk_scalar_var(a, size));
      do_weak_update(scalar_var, val);
    }
    CRAB_LOG("smashing", crab::outs() << a << "[" << i << ".." << j << "]:="
                                      << val << " -- " << *this << "\n";);
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {

    bytes_t size = get_size(rhs);
    if (size.is_constant()) {
      set_size(lhs, size.get_constant());
      variable_t scalar_lhs(mk_scalar_var(lhs, size.get_constant()));
      variable_t scalar_rhs(mk_scalar_var(rhs, size.get_constant()));

      auto ty = scalar_lhs.get_type();
      if (ty.is_bool()) {
        m_base_dom.assign_bool_var(scalar_lhs, scalar_rhs, false);
      } else {
        assert(ty.is_integer() || ty.is_real());
        m_base_dom.assign(scalar_lhs, scalar_rhs);
      }
    }
  }

  // backward array operations
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const array_smashing_t &invariant) override {
    CRAB_WARN("backward_array_init in array smashing domain not implemented");
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const array_smashing_t &invariant) override {
    CRAB_WARN("backward_array_load in array smashing domain not implemented");
    this->operator-=(lhs);
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const array_smashing_t &invariant) override {
    CRAB_WARN("backward_array_store in array smashing domain not implemented");
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const array_smashing_t &invariant) override {
    CRAB_WARN(
        "backward_array_store_range in array smashing domain not implemented");
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const array_smashing_t &invariant) override {
    CRAB_WARN("backward_array_assign in array smashing domain not implemented");
  }

  /// array_smashing is a functor domain that implements all
  /// operations except region/reference operations.
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(array_smashing_t)

  linear_constraint_system_t to_linear_constraint_system() const override {
    return filter_ghost_vars(
        std::move(m_base_dom.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;

    auto disj_csts = m_base_dom.to_disjunctive_linear_constraint_system();
    for (auto &csts : disj_csts) {
      auto filtered_csts = filter_ghost_vars(std::move(csts));
      if (!filtered_csts.is_true()) {
        res += filtered_csts;
      }
    }
    return res;
  }

  /* Deprecated: do not use them */
  base_dom_t &get_content_domain() { return m_base_dom; }
  const base_dom_t &get_content_domain() const { return m_base_dom; }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    m_base_dom.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const array_smashing_t &invariant) override {
    m_base_dom.backward_intrinsic(name, inputs, outputs, invariant.m_base_dom);
  }
  /* end intrinsics operations */

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename expects vectors same sizes");
    }

    variable_vector_t old_vars, new_vars;
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      if (from[i].get_type() != to[i].get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve the same type");
      }
      const variable_t &old_var = from[i];
      const variable_t &new_var = to[i];

      if (old_var.get_type().is_array()) {
        bytes_t size = get_size(old_var);
        if (size.is_constant()) {
          set_size(new_var, size.get_constant());
          variable_t gold_var(mk_scalar_var(old_var, size.get_constant()));
          variable_t gnew_var(mk_scalar_var(new_var, size.get_constant()));
          old_vars.push_back(gold_var);
          new_vars.push_back(gnew_var);
        }
      } else {
        old_vars.push_back(old_var);
        new_vars.push_back(new_var);
      }
    }
    // No array variables should go to m_base_dom
    m_base_dom.rename(old_vars, new_vars);
  }

  void write(crab_os &o) const override { o << m_base_dom; }

  std::string domain_name() const override {
    std::string name("ArraySmashing(" + m_base_dom.domain_name() + ")");
    return name;
  }

}; // end array_smashing

template <typename BaseDomain>
struct abstract_domain_traits<array_smashing<BaseDomain>> {
  using number_t = typename BaseDomain::number_t;
  using varname_t = typename BaseDomain::varname_t;
};

// template <typename BaseDom>
// class checker_domain_traits<array_smashing<BaseDom>> {
// public:
//   using this_type = array_smashing<BaseDom>;
//   using linear_constraint_t = typename this_type::linear_constraint_t;
//   using disjunctive_linear_constraint_system_t =
//       typename this_type::disjunctive_linear_constraint_system_t;
//   static bool entail(this_type &lhs,
//                      const disjunctive_linear_constraint_system_t &rhs) {
//     BaseDom &lhs_dom = lhs.get_content_domain();
//     return checker_domain_traits<BaseDom>::entail(lhs_dom, rhs);
//   }
//   static bool entail(const disjunctive_linear_constraint_system_t &lhs,
//                      this_type &rhs) {
//     BaseDom &rhs_dom = rhs.get_content_domain();
//     return checker_domain_traits<BaseDom>::entail(lhs, rhs_dom);
//   }
//   static bool intersect(this_type &inv, const linear_constraint_t &cst) {
//     BaseDom &dom = inv.get_content_domain();
//     return checker_domain_traits<BaseDom>::intersect(dom, cst);
//   }
// };

} // namespace domains
} // namespace crab
