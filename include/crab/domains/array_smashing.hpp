/*******************************************************************************
 * Array smashing domain
 *
 * Word-level assumption: for any array load, it assumes that the
 * number of read bytes is equal to the number of bytes written during
 * the last array write. If this assumption is too strong then use the
 * array_adaptive domain.
 *
 * For efficiency reasons, the domain does not generate a ghost
 * variable for the summarized variable (i.e., the variable that
 * represents the smashed array). This means that array variables are
 * passed to the bool/numerical domain. In particular, the smashing
 * domain passes an array variable to the base domain in the following
 * operations: assign_bool_cst, assign_bool_var, assign, and expand
 * for the modeling of array operations.
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {

namespace domains {

// Abstract domain to reason about summarized variables. All
// array elements are `smashed` into a single variable.
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
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using base_dom_t = BaseNumDomain;
  using interval_t = ikos::interval<number_t>;

private:
  using bound_t = ikos::bound<number_t>;

  // Contain scalar and summarized array variables.
  //
  // XXX: We need to be careful in methods such as
  // to_linear_constraint_system and
  // to_disjunctive_linear_constraint_system. These methods
  // convert the internal representation to linear constraints
  // which should not contain array variables.
  base_dom_t m_base_dom;

  array_smashing(base_dom_t &&base_dom) : m_base_dom(std::move(base_dom)) {}

  void do_strong_update(base_dom_t &dom, const variable_t &a,
                        const linear_expression_t &rhs) {
    auto ty = a.get_type();
    if (ty.is_bool_array()) {
      if (rhs.is_constant()) {
        if (rhs.constant() >= number_t(1)) {
          dom.assign_bool_cst(a, linear_constraint_t::get_true());
        } else {
          dom.assign_bool_cst(a, linear_constraint_t::get_false());
        }
      } else if (auto rhs_v = rhs.get_variable()) {
        dom.assign_bool_var(a, (*rhs_v), false);
      }
    } else if (ty.is_integer_array() || ty.is_real_array()) {
      dom.assign(a, rhs);
    } else {
      /* unreachable */
    }
  }

  void do_strong_update(const variable_t &a, const linear_expression_t &rhs) {
    do_strong_update(m_base_dom, a, rhs);
  }

  // We perform the strong update on a copy of *this. Then, we join
  // the copy of *this with *this.
  void do_weak_update(const variable_t &a, const linear_expression_t &rhs) {
    base_dom_t other(m_base_dom);
    do_strong_update(other, a, rhs);
    m_base_dom |= other;
  }

  // The internal representation contains variables of array
  // type and add them as dimensions in the underlying numerical
  // domain. This is OK but it shouldn't be exposed outside via
  // linear constraints.
  linear_constraint_system_t
  filter_noninteger_vars(linear_constraint_system_t &&csts) const {
    linear_constraint_system_t res;
    for (auto const &cst : csts) {
      if (std::all_of(
              cst.expression().variables_begin(),
              cst.expression().variables_end(), [](const variable_t &v) {
                return v.get_type().is_integer() || v.get_type().is_bool();
              })) {
        res += cst;
      }
    }
    return res;
  }

public:
  array_smashing() { m_base_dom.set_to_top(); }

  array_smashing make_top() const override {
    base_dom_t base_dom;
    array_smashing out(base_dom.make_top());
    return out;
  }

  array_smashing make_bottom() const override {
    base_dom_t base_dom;
    array_smashing out(base_dom.make_bottom());
    return out;
  }

  void set_to_top() override {
    base_dom_t base_dom;
    array_smashing abs(base_dom.make_top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    base_dom_t base_dom;
    array_smashing abs(base_dom.make_bottom());
    std::swap(*this, abs);
  }

  array_smashing(const array_smashing_t &other) : m_base_dom(other.m_base_dom) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  array_smashing(const array_smashing_t &&other)
      : m_base_dom(std::move(other.m_base_dom)) {}

  array_smashing_t &operator=(const array_smashing_t &other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      m_base_dom = other.m_base_dom;
    }
    return *this;
  }

  array_smashing_t &operator=(const array_smashing_t &&other) {
    if (this != &other) {
      m_base_dom = std::move(other.m_base_dom);
    }
    return *this;
  }

  bool is_bottom() const override { return (m_base_dom.is_bottom()); }

  bool is_top() const override { return (m_base_dom.is_top()); }

  bool operator<=(const array_smashing_t &other) const override {
    return (m_base_dom <= other.m_base_dom);
  }

  void operator|=(const array_smashing_t &other) override {
    m_base_dom |= other.m_base_dom;
  }

  array_smashing_t operator|(const array_smashing_t &other) const override {
    return array_smashing_t(m_base_dom | other.m_base_dom);
  }

  array_smashing_t operator&(const array_smashing_t &other) const override {
    return array_smashing_t(m_base_dom & other.m_base_dom);
  }

  array_smashing_t operator||(const array_smashing_t &other) const override {
    return array_smashing_t(m_base_dom || other.m_base_dom);
  }

  array_smashing_t widening_thresholds(
      const array_smashing_t &other,
      const iterators::thresholds<number_t> &ts) const override {
    return array_smashing_t(m_base_dom.widening_thresholds(other.m_base_dom, ts));
  }

  array_smashing_t operator&&(const array_smashing_t &other) const override {
    return array_smashing_t(m_base_dom && other.m_base_dom);
  }

  virtual interval_t operator[](const variable_t &v) override {
    return m_base_dom[v];
  }

  virtual interval_t at(const variable_t &v) const override {
    return m_base_dom.at(v);
  }

  void forget(const variable_vector_t &variables) override {
    m_base_dom.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    m_base_dom.project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (var.get_type() != new_var.get_type()) {
      CRAB_ERROR(domain_name(), "::expand must preserve the same type");
    }
    m_base_dom.expand(var, new_var);
  }

  void normalize() override { m_base_dom.normalize(); }

  void minimize() override { m_base_dom.minimize(); }

  void operator+=(const linear_constraint_system_t &csts) override {
    m_base_dom += csts;
  }

  void operator-=(const variable_t &var) override { m_base_dom -= var; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_base_dom.assign(x, e);

    CRAB_LOG("smashing", crab::outs()
                             << "apply " << x << " := " << e << *this << "\n";);
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
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
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

  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    m_base_dom.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    m_base_dom.assume_bool(v, is_negated);
  }

  virtual void select_bool(const variable_t &lhs, const variable_t &cond,
			   const variable_t &b1, const variable_t &b2) override {
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
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t & /*lb_idx*/,
                          const linear_expression_t & /*ub_idx*/,
                          const linear_expression_t &val) override {
    auto ty = a.get_type();
    if (ty.is_bool_array()) {
      if (val.is_constant()) {
        if (val.constant() >= number_t(1)) {
          m_base_dom.assign_bool_cst(a, linear_constraint_t::get_true());
        } else {
          m_base_dom.assign_bool_cst(a, linear_constraint_t::get_false());
        }
      } else if (auto var = val.get_variable()) {
        m_base_dom.assign_bool_var(a, (*var), false);
      }
    } else if (ty.is_integer_array() || ty.is_real_array()) {
      m_base_dom.assign(a, val);
    }
    CRAB_LOG("smashing", crab::outs() << "forall i:: " << a << "[i]==" << val
                                      << " -- " << *this << "\n";);
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.load");
    crab::ScopedCrabStats __st__(domain_name() + ".load");

    // We need to be careful when assigning a summarized variable a
    // into a non-summarized variable lhs. Simply m_base_dom.assign(lhs,a)
    // is not sound.
    auto &vfac = const_cast<varname_t *>(&(a.name()))->get_var_factory();
    variable_t a_prime(vfac.get());
    m_base_dom.expand(a, a_prime);
    auto ty = a.get_type();
    if (ty.is_bool_array()) {
      m_base_dom.assign_bool_var(lhs, a_prime, false);
    } else if (ty.is_integer_array() || ty.is_real_array()) {
      m_base_dom.assign(lhs, a_prime);
    }
    m_base_dom -= a_prime;

    CRAB_LOG("smashing", crab::outs() << lhs << ":=" << a << "[" << i
                                      << "]  -- " << *this << "\n";);
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t & /*elem_size*/,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.store");
    crab::ScopedCrabStats __st__(domain_name() + ".store");

    if (is_strong_update) {
      do_strong_update(a, val);
    } else {
      do_weak_update(a, val);
    }

    CRAB_LOG("smashing", crab::outs() << a << "[" << i << "]:=" << val << " -- "
                                      << *this << "\n";);
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t & /*elem_size*/,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.store");
    crab::ScopedCrabStats __st__(domain_name() + ".store");
    do_weak_update(a, val);
    CRAB_LOG("smashing", crab::outs() << a << "[" << i << ".." << j << "]:="
                                      << val << " -- " << *this << "\n";);
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    auto ty = lhs.get_type();
    if (ty.is_bool_array()) {
      m_base_dom.assign_bool_var(lhs, rhs, false);
    } else if (ty.is_integer_array() || ty.is_real_array()) {
      m_base_dom.assign(lhs, rhs);
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
    return filter_noninteger_vars(
        std::move(m_base_dom.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;

    auto disj_csts = m_base_dom.to_disjunctive_linear_constraint_system();
    for (auto &csts : disj_csts) {
      auto filtered_csts = filter_noninteger_vars(std::move(csts));
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
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
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
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      if (from[i].get_type() != to[i].get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve the same type");
      }
    }
    m_base_dom.rename(from, to);
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

template <typename BaseDom>
class checker_domain_traits<array_smashing<BaseDom>> {
public:
  using this_type = array_smashing<BaseDom>;
  using linear_constraint_t = typename this_type::linear_constraint_t;
  using disjunctive_linear_constraint_system_t =
      typename this_type::disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    BaseDom &lhs_dom = lhs.get_content_domain();
    return checker_domain_traits<BaseDom>::entail(lhs_dom, rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    BaseDom &rhs_dom = rhs.get_content_domain();
    return checker_domain_traits<BaseDom>::entail(lhs, rhs_dom);
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    BaseDom &lhs_dom = lhs.get_content_domain();
    return checker_domain_traits<BaseDom>::entail(lhs_dom, rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    BaseDom &dom = inv.get_content_domain();
    return checker_domain_traits<BaseDom>::intersect(dom, cst);
  }
};

} // namespace domains
} // namespace crab
