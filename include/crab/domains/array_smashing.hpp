/*******************************************************************************
 * Array smashing domain
 *
 * Assume all array accesses are aligned wrt to the size of the array
 * element (e.g., if the size of the array element is 4 bytes then all
 * array accesses must be multiple of 4). Note that this assumption
 * does not hold in general, so the client must ensure that.
 ******************************************************************************/

#pragma once

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

namespace crab {

namespace domains {

// Abstract domain to reason about summarized variables. All
// array elements are `smashed` into a single variable.
template <typename NumDomain>
class array_smashing final : public abstract_domain<array_smashing<NumDomain>> {

public:
  typedef typename NumDomain::number_t number_t;
  typedef typename NumDomain::varname_t varname_t;

private:
  typedef array_smashing<NumDomain> array_smashing_t;
  typedef abstract_domain<array_smashing_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef crab::pointer_constraint<variable_t> ptr_cst_t;
  typedef NumDomain content_domain_t;
  typedef interval<number_t> interval_t;

private:
  typedef bound<number_t> bound_t;

  // Contain scalar and summarized array variables.
  //
  // XXX: We need to be careful in methods such as
  // to_linear_constraint_system and
  // to_disjunctive_linear_constraint_system. These methods
  // convert the internal representation to linear constraints
  // which should not contain array variables.
  NumDomain _inv;

  array_smashing(NumDomain &&inv) : _inv(std::move(inv)) {}

  void do_strong_update(variable_t a, linear_expression_t rhs) {
    switch (a.get_type()) {
    case ARR_BOOL_TYPE:
      if (rhs.is_constant()) {
        if (rhs.constant() >= number_t(1)) {
          _inv.assign_bool_cst(a, linear_constraint_t::get_true());
        } else {
          _inv.assign_bool_cst(a, linear_constraint_t::get_false());
        }
      } else if (auto rhs_v = rhs.get_variable()) {
        _inv.assign_bool_var(a, (*rhs_v), false);
      }
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(a, rhs);
      break;
    case ARR_PTR_TYPE:
      if (rhs.is_constant() && rhs.constant() == number_t(0)) {
        _inv.pointer_mk_null(a);
      } else if (auto rhs_v = rhs.get_variable()) {
        _inv.pointer_assign(a, (*rhs_v), number_t(0));
      }
      break;
    default:; /* unreachable */
    }
  }

  // We perform the strong update on a copy of *this using a_new
  // while adding the assignment a_new = a_old on *this (if
  // a_old is defined). Then, we join the copy of *this with *this.
  void do_weak_update(variable_t a_new, boost::optional<variable_t> a_old,
                      linear_expression_t rhs) {
    NumDomain other(_inv);
    switch (a_new.get_type()) {
    case ARR_BOOL_TYPE:
      if (rhs.is_constant()) {
        if (rhs.constant() >= number_t(1)) {
          other.assign_bool_cst(a_new, linear_constraint_t::get_true());
          if (a_old)
            _inv.assign_bool_var(a_new, *a_old, false);
        } else {
          other.assign_bool_cst(a_new, linear_constraint_t::get_false());
          if (a_old)
            _inv.assign_bool_var(a_new, *a_old, false);
        }
      } else if (auto rhs_v = rhs.get_variable()) {
        other.assign_bool_var(a_new, (*rhs_v), false);
        if (a_old)
          _inv.assign_bool_var(a_new, *a_old, false);
      }
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      other.assign(a_new, rhs);
      if (a_old)
        _inv.assign(a_new, *a_old);
      break;
    case ARR_PTR_TYPE:
      if (rhs.is_constant() && rhs.constant() == number_t(0)) {
        other.pointer_mk_null(a_new);
        if (a_old)
          _inv.pointer_assign(a_new, *a_old, number_t(0));
      } else if (auto rhs_v = rhs.get_variable()) {
        other.pointer_assign(a_new, (*rhs_v), number_t(0));
        if (a_old)
          _inv.pointer_assign(a_new, *a_old, number_t(0));
      }
      break;
    default:; /* unreachable */
    }

    // join
    _inv |= other;
  }

  // The internal representation contains variables of array
  // type and add them as dimensions in the underlying numerical
  // domain. This is OK but it shouldn't be exposed outside via
  // linear constraints.
  linear_constraint_system_t
  filter_noninteger_vars(linear_constraint_system_t &&csts) {
    linear_constraint_system_t res;
    for (auto &cst : csts) {
      auto vars = cst.variables();
      if (std::all_of(vars.begin(), vars.end(), [](const variable_t &v) {
            return v.is_int_type() || v.is_bool_type();
          })) {
        res += cst;
      }
    }
    return res;
  }

public:
  array_smashing() : _inv(NumDomain::top()) {}

  void set_to_top() {
    array_smashing abs(NumDomain::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    array_smashing abs(NumDomain::bottom());
    std::swap(*this, abs);
  }

  array_smashing(const array_smashing_t &other) : _inv(other._inv) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  array_smashing(const array_smashing_t &&other)
      : _inv(std::move(other._inv)) {}

  array_smashing_t &operator=(const array_smashing_t &other) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
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

  bool is_bottom() { return (_inv.is_bottom()); }

  bool is_top() { return (_inv.is_top()); }

  bool operator<=(array_smashing_t other) { return (_inv <= other._inv); }

  void operator|=(array_smashing_t other) { _inv |= other._inv; }

  array_smashing_t operator|(array_smashing_t other) {
    return array_smashing_t(_inv | other._inv);
  }

  array_smashing_t operator&(array_smashing_t other) {
    return array_smashing_t(_inv & other._inv);
  }

  array_smashing_t operator||(array_smashing_t other) {
    return array_smashing_t(_inv || other._inv);
  }

  array_smashing_t
  widening_thresholds(array_smashing_t other,
                      const iterators::thresholds<number_t> &ts) {
    return array_smashing_t(_inv.widening_thresholds(other._inv, ts));
  }

  array_smashing_t operator&&(array_smashing_t other) {
    return array_smashing_t(_inv && other._inv);
  }

  interval_t operator[](const variable_t &v) { return _inv[v]; }

  void forget(const variable_vector_t &variables) { _inv.forget(variables); }

  void project(const variable_vector_t &variables) { _inv.project(variables); }

  void expand(variable_t var, variable_t new_var) {
    CRAB_WARN("array smashing expand not implemented");
  }

  void normalize() { _inv.normalize(); }

  void minimize() { _inv.minimize(); }

  void operator+=(linear_constraint_system_t csts) { _inv += csts; }

  void operator-=(variable_t var) { _inv -= var; }

  void assign(variable_t x, linear_expression_t e) {
    _inv.assign(x, e);

    CRAB_LOG("smashing", crab::outs()
                             << "apply " << x << " := " << e << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t z) {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       array_smashing_t inv) {
    _inv.backward_assign(x, e, inv.get_content_domain());
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      array_smashing_t inv) {
    _inv.backward_apply(op, x, y, z, inv.get_content_domain());
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      array_smashing_t inv) {
    _inv.backward_apply(op, x, y, z, inv.get_content_domain());
  }

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    _inv.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    _inv.apply(op, x, y, z);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << z << *this << "\n";);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    _inv.apply(op, x, y, k);

    CRAB_LOG("smashing", crab::outs() << "apply " << x << " := " << y << " "
                                      << op << " " << k << *this << "\n";);
  }

  // boolean operators
  virtual void assign_bool_cst(variable_t lhs,
                               linear_constraint_t rhs) override {
    _inv.assign_bool_cst(lhs, rhs);
  }

  virtual void assign_bool_var(variable_t lhs, variable_t rhs,
                               bool is_not_rhs) override {
    _inv.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  virtual void apply_binary_bool(bool_operation_t op, variable_t x,
                                 variable_t y, variable_t z) override {
    _inv.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(variable_t v, bool is_negated) override {
    _inv.assume_bool(v, is_negated);
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                        array_smashing_t inv) {
    _inv.backward_assign_bool_cst(lhs, rhs, inv.get_content_domain());
  }

  virtual void backward_assign_bool_var(variable_t lhs, variable_t rhs,
                                        bool is_not_rhs, array_smashing_t inv) {
    _inv.backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                  inv.get_content_domain());
  }

  virtual void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                          variable_t y, variable_t z,
                                          array_smashing_t inv) {
    _inv.backward_apply_binary_bool(op, x, y, z, inv.get_content_domain());
  }

  // pointer_operators_api
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    _inv.pointer_load(lhs, rhs);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    _inv.pointer_store(lhs, rhs);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    _inv.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    _inv.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    _inv.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    _inv.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(ptr_cst_t cst) override {
    _inv.pointer_assume(cst);
  }

  virtual void pointer_assert(ptr_cst_t cst) override {
    _inv.pointer_assert(cst);
  }

  // array_operators_api

  // All the array elements are initialized to val
  virtual void array_init(variable_t a, linear_expression_t /*elem_size*/,
                          linear_expression_t /*lb_idx*/,
                          linear_expression_t /*ub_idx*/,
                          linear_expression_t val) override {
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
    case ARR_PTR_TYPE:
      if (val.is_constant() && val.constant() == number_t(0)) {
        _inv.pointer_mk_null(a);
      } else if (auto var = val.get_variable()) {
        _inv.pointer_assign(a, (*var), number_t(0));
      }
      break;
    default:;
    }
    CRAB_LOG("smashing", crab::outs() << "forall i:: " << a << "[i]==" << val
                                      << " -- " << *this << "\n";);
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t /*elem_size*/,
                          linear_expression_t i) override {
    crab::CrabStats::count(getDomainName() + ".count.load");
    crab::ScopedCrabStats __st__(getDomainName() + ".load");

    // We need to be careful when assigning a summarized variable a
    // into a non-summarized variable lhs. Simply _inv.assign(lhs,
    // a) is not sound.
    /* ask for a temp var */
    variable_t a_prime(a.name().get_var_factory().get());
    _inv.expand(a, a_prime);
    switch (a.get_type()) {
    case ARR_BOOL_TYPE:
      _inv.assign_bool_var(lhs, a_prime, false);
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(lhs, a_prime);
      break;
    case ARR_PTR_TYPE:
      _inv.pointer_assign(lhs, a_prime, number_t(0));
      break;
    default:;
    }
    _inv -= a_prime;

    CRAB_LOG("smashing", crab::outs() << lhs << ":=" << a << "[" << i
                                      << "]  -- " << *this << "\n";);
  }

  virtual void array_store(variable_t a, linear_expression_t /*elem_size*/,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    crab::CrabStats::count(getDomainName() + ".count.store");
    crab::ScopedCrabStats __st__(getDomainName() + ".store");

    if (is_strong_update) {
      do_strong_update(a, val);
    } else {
      do_weak_update(a, boost::none, val);
    }

    CRAB_LOG("smashing", crab::outs() << a << "[" << i << "]:=" << val << " -- "
                                      << *this << "\n";);
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t /*elem_size*/,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    crab::CrabStats::count(getDomainName() + ".count.store");
    crab::ScopedCrabStats __st__(getDomainName() + ".store");

    if (is_strong_update) {
      do_strong_update(a_new, val);
      // a_old is unused if strong update
    } else {
      do_weak_update(a_new, a_old, val);
    }

    CRAB_LOG("smashing", crab::outs() << a_new << ":= " << a_old << "[" << i
                                      << "<-" << val << "]"
                                      << " -- " << *this << "\n";);
  }

  virtual void array_store_range(variable_t a,
                                 linear_expression_t /*elem_size*/,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.store");
    crab::ScopedCrabStats __st__(getDomainName() + ".store");
    do_weak_update(a, boost::none, val);
    CRAB_LOG("smashing", crab::outs() << a << "[" << i << ".." << j << "]:="
                                      << val << " -- " << *this << "\n";);
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t /*elem_size*/,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.store");
    crab::ScopedCrabStats __st__(getDomainName() + ".store");
    do_weak_update(a_new, a_old, val);
    CRAB_LOG("smashing", crab::outs()
                             << a_new << ":= " << a_old << "[" << i << ".." << j
                             << "<-" << val << "] -- " << *this << "\n";);
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    switch (lhs.get_type()) {
    case ARR_BOOL_TYPE:
      _inv.assign_bool_var(lhs, rhs, false);
      break;
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
      _inv.assign(lhs, rhs);
      break;
    case ARR_PTR_TYPE:
      _inv.pointer_assign(lhs, rhs, number_t(0));
      break;
    default:;
    }
  }

  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           array_smashing_t invariant) {
    CRAB_WARN("backward_array_init in array smashing domain not implemented");
  }
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           array_smashing_t invariant) {
    CRAB_WARN("backward_array_load in array smashing domain not implemented");
    this->operator-=(lhs);
  }
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, array_smashing_t invariant) {
    CRAB_WARN("backward_array_store in array smashing domain not implemented");
  }
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, array_smashing_t invariant) {
    CRAB_WARN("backward_array_store in array smashing domain not implemented");
  }
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  array_smashing_t invariant) {
    CRAB_WARN(
        "backward_array_store_range in array smashing domain not implemented");
  }
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  array_smashing_t invariant) {
    CRAB_WARN(
        "backward_array_store_range in array smashing domain not implemented");
  }
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             array_smashing_t invariant) {
    CRAB_WARN("backward_array_assign in array smashing domain not implemented");
  }

  linear_constraint_system_t to_linear_constraint_system() {
    return filter_noninteger_vars(
        std::move(_inv.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
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
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  array_smashing_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    _inv.rename(from, to);
  }
  
  void write(crab_os &o) { o << _inv; }

  static std::string getDomainName() {
    std::string name("ArraySmashing(" + NumDomain::getDomainName() + ")");
    return name;
  }

}; // end array_smashing

template <typename BaseDomain>
struct abstract_domain_traits<array_smashing<BaseDomain>> {
  typedef typename BaseDomain::number_t number_t;
  typedef typename BaseDomain::varname_t varname_t;
};

template <typename BaseDom>
class checker_domain_traits<array_smashing<BaseDom>> {
public:
  typedef array_smashing<BaseDom> this_type;
  typedef typename this_type::linear_constraint_t linear_constraint_t;
  typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;

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
