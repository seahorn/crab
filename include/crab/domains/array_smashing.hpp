/*******************************************************************************
 * Array smashing domain
 *
 * Assume all array accesses are aligned wrt to the size of the array
 * element (e.g., if the size of the array element is 4 bytes then all
 * array accesses must be multiple of 4). Note that this assumption
 * does not hold in general, so the client must ensure that.
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
template <typename NumDomain>
class array_smashing final : public abstract_domain<array_smashing<NumDomain>> {

public:
  using number_t = typename NumDomain::number_t;
  using varname_t = typename NumDomain::varname_t;

private:
  using array_smashing_t = array_smashing<NumDomain>;
  using abstract_domain_t = abstract_domain<array_smashing_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using content_domain_t = NumDomain;
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
  NumDomain _inv;

  array_smashing(NumDomain &&inv) : _inv(std::move(inv)) {}

  void do_strong_update(NumDomain &dom, const variable_t &a,
                        const linear_expression_t &rhs) {
    switch (a.get_type()) {
    case ARR_BOOL_TYPE:
      if (rhs.is_constant()) {
        if (rhs.constant() >= number_t(1)) {
          dom.assign_bool_cst(a, linear_constraint_t::get_true());
        } else {
          dom.assign_bool_cst(a, linear_constraint_t::get_false());
        }
      } else if (auto rhs_v = rhs.get_variable()) {
        dom.assign_bool_var(a, (*rhs_v), false);
      }
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      dom.assign(a, rhs);
      break;
    default:; /* unreachable */
    }
  }

  void do_strong_update(const variable_t &a, const linear_expression_t &rhs) {
    do_strong_update(_inv, a, rhs);
  }

  // We perform the strong update on a copy of *this. Then, we join
  // the copy of *this with *this.
  void do_weak_update(const variable_t &a, const linear_expression_t &rhs) {
    NumDomain other(_inv);
    do_strong_update(other, a, rhs);
    _inv |= other;
  }

  // The internal representation contains variables of array
  // type and add them as dimensions in the underlying numerical
  // domain. This is OK but it shouldn't be exposed outside via
  // linear constraints.
  linear_constraint_system_t
  filter_noninteger_vars(linear_constraint_system_t &&csts) const {
    linear_constraint_system_t res;
    for (auto const &cst : csts) {
      if (std::all_of(cst.expression().variables_begin(),
                      cst.expression().variables_end(),
                      [](const variable_t &v) {
                        return v.is_int_type() || v.is_bool_type();
                      })) {
        res += cst;
      }
    }
    return res;
  }

public:
  array_smashing() { _inv.set_to_top(); }

  array_smashing make_top() const override {
    NumDomain inv;
    array_smashing out(inv.make_top());
    return out;
  }

  array_smashing make_bottom() const override {
    NumDomain inv;
    array_smashing out(inv.make_bottom());
    return out;
  }

  void set_to_top() override {
    NumDomain inv;
    array_smashing abs(inv.make_top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    NumDomain inv;
    array_smashing abs(inv.make_bottom());
    std::swap(*this, abs);
  }

  array_smashing(const array_smashing_t &other) : _inv(other._inv) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  array_smashing(const array_smashing_t &&other)
      : _inv(std::move(other._inv)) {}

  array_smashing_t &operator=(const array_smashing_t &other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      _inv = other._inv;
    }
    return *this;
  }

  array_smashing_t &operator=(const array_smashing_t &&other) {
    if (this != &other) {
      _inv = std::move(other._inv);
    }
    return *this;
  }

  bool is_bottom() const override { return (_inv.is_bottom()); }

  bool is_top() const override { return (_inv.is_top()); }

  bool operator<=(const array_smashing_t &other) const override {
    return (_inv <= other._inv);
  }

  void operator|=(const array_smashing_t &other) override {
    _inv |= other._inv;
  }

  array_smashing_t operator|(const array_smashing_t &other) const override {
    return array_smashing_t(_inv | other._inv);
  }

  array_smashing_t operator&(const array_smashing_t &other) const override {
    return array_smashing_t(_inv & other._inv);
  }

  array_smashing_t operator||(const array_smashing_t &other) const override {
    return array_smashing_t(_inv || other._inv);
  }

  array_smashing_t widening_thresholds(
      const array_smashing_t &other,
      const iterators::thresholds<number_t> &ts) const override {
    return array_smashing_t(_inv.widening_thresholds(other._inv, ts));
  }

  array_smashing_t operator&&(const array_smashing_t &other) const override {
    return array_smashing_t(_inv && other._inv);
  }

  virtual interval_t operator[](const variable_t &v) override {
    return _inv[v];
  }

  void forget(const variable_vector_t &variables) override {
    _inv.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    _inv.project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    CRAB_WARN("array smashing expand not implemented");
  }

  void normalize() override { _inv.normalize(); }

  void minimize() override { _inv.minimize(); }

  void operator+=(const linear_constraint_system_t &csts) override {
    _inv += csts;
  }

  void operator-=(const variable_t &var) override { _inv -= var; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    _inv.assign(x, e);

    CRAB_LOG("smashing", crab::outs()
                             << "apply " << x << " := " << e << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const array_smashing_t &inv) override {
    _inv.backward_assign(x, e, inv._inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const array_smashing_t &inv) override {
    _inv.backward_apply(op, x, y, z, inv._inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const array_smashing_t &inv) override {
    _inv.backward_apply(op, x, y, z, inv._inv);
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    _inv.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    _inv.apply(op, x, y, k);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << k << *this << "\n";);
  }

  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    _inv.assign_bool_cst(lhs, rhs);
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    _inv.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    _inv.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    _inv.assume_bool(v, is_negated);
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const array_smashing_t &inv) override {
    _inv.backward_assign_bool_cst(lhs, rhs, inv._inv);
  }

  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const array_smashing_t &inv) override {
    _inv.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv._inv);
  }

  virtual void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const array_smashing_t &inv) override {
    _inv.backward_apply_binary_bool(op, x, y, z, inv._inv);
  }

  // array_operators_api

  // All the array elements are initialized to val
  virtual void array_init(const variable_t &a,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t & /*lb_idx*/,
                          const linear_expression_t & /*ub_idx*/,
                          const linear_expression_t &val) override {
    switch (a.get_type()) {
    case ARR_BOOL_TYPE: {
      if (val.is_constant()) {
        if (val.constant() >= number_t(1)) {
          _inv.assign_bool_cst(a, linear_constraint_t::get_true());
        } else {
          _inv.assign_bool_cst(a, linear_constraint_t::get_false());
        }
      } else if (auto var = val.get_variable()) {
        _inv.assign_bool_var(a, (*var), false);
      }
      break;
    }
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(a, val);
      break;
    default:;
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
    // into a non-summarized variable lhs. Simply _inv.assign(lhs,a)
    // is not sound.
    auto &vfac = const_cast<varname_t *>(&(a.name()))->get_var_factory();
    variable_t a_prime(vfac.get());
    _inv.expand(a, a_prime);
    switch (a.get_type()) {
    case ARR_BOOL_TYPE:
      _inv.assign_bool_var(lhs, a_prime, false);
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(lhs, a_prime);
      break;
    default:;
    }
    _inv -= a_prime;

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
    switch (lhs.get_type()) {
    case ARR_BOOL_TYPE:
      _inv.assign_bool_var(lhs, rhs, false);
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(lhs, rhs);
      break;
    default:;
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

  // reference operations
  void region_init(const memory_region &reg) override {}
  void ref_make(const variable_t &ref, const memory_region &reg) override {}
  void ref_load(const variable_t &ref, const memory_region &reg,
                const variable_t &res) override {}
  void ref_store(const variable_t &ref, const memory_region &reg,
                 const linear_expression_t &val) override {}
  void ref_gep(const variable_t &ref1, const memory_region &reg1,
               const variable_t &ref2, const memory_region &reg2,
               const linear_expression_t &offset) override {}
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const memory_region &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {}
  void ref_store_to_array(const variable_t &ref, const memory_region &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {}
  void ref_assume(const reference_constraint_t &cst) override {}

  linear_constraint_system_t to_linear_constraint_system() const override {
    return filter_noninteger_vars(
        std::move(_inv.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;

    auto disj_csts = _inv.to_disjunctive_linear_constraint_system();
    for (auto &csts : disj_csts) {
      auto filtered_csts = filter_noninteger_vars(std::move(csts));
      if (!filtered_csts.is_true()) {
        res += filtered_csts;
      }
    }
    return res;
  }

  NumDomain &get_content_domain() { return _inv; }

  NumDomain get_content_domain() const { return _inv; }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    _inv.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const array_smashing_t &invariant) override {
    _inv.backward_intrinsic(name, inputs, outputs, invariant._inv);
  }
  /* end intrinsics operations */

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    _inv.rename(from, to);
  }

  void write(crab_os &o) const override { o << _inv; }

  std::string domain_name() const override {
    std::string name("ArraySmashing(" + _inv.domain_name() + ")");
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
