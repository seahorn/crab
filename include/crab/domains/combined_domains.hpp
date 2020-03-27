#pragma once

/**************************************************************************
 * Combination of domains:
 *
 * (1) reduced product of two arbitrary domains with only lattice
 *     operations.
 *
 * (2) reduced product of two arbitrary domains with all operations.
 *
 * The reduction in (1) and (2) is simply done by making bottom the
 * abstract state if one of them is bottom.
 *
 * (3) reduced product of two numerical domains. The reduction is done
 *     via a special push operation that must be defined by the
 *     domains. This is often more precise than (2).
 *
 * (4) reduced product of an arbitrary numerical domain and congruences.
 *
 * (5) reduced product of an arbitrary numerical domain and nullity.
 **************************************************************************/

#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>

#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/nullity.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

// Provided by Ikos
// Reduced product of two arbitrary domains with only lattice
// operations.
template <typename Domain1, typename Domain2>
class basic_domain_product2 : public writeable {

public:
  typedef basic_domain_product2<Domain1, Domain2> basic_domain_product2_t;
  typedef Domain1 first_type;
  typedef Domain2 second_type;

private:
  bool _is_bottom;
  Domain1 _first;
  Domain2 _second;

  void canonicalize() {
    if (!this->_is_bottom) {
      this->_is_bottom = this->_first.is_bottom() || this->_second.is_bottom();
      if (this->_is_bottom) {
        this->_first = Domain1::bottom();
        this->_second = Domain2::bottom();
      }
    }
  }

public:
  basic_domain_product2()
      : _is_bottom(false), _first(Domain1::top()), _second(Domain2::top()) {}

  basic_domain_product2(Domain1 &&first, Domain2 &&second,
                        bool &&apply_reduction = true)
      : _is_bottom(false), _first(std::move(first)),
        _second(std::move(second)) {
    if (apply_reduction) {
      this->canonicalize();
    }
  }

  basic_domain_product2(const basic_domain_product2_t &other)
      : _is_bottom(other._is_bottom), _first(other._first),
        _second(other._second) {}

  basic_domain_product2(const basic_domain_product2_t &&other)
      : _is_bottom(std::move(other._is_bottom)),
        _first(std::move(other._first)), _second(std::move(other._second)) {}

  basic_domain_product2_t &operator=(const basic_domain_product2_t &other) {
    if (this != &other) {
      this->_is_bottom = other._is_bottom;
      this->_first = other._first;
      this->_second = other._second;
    }
    return *this;
  }

  basic_domain_product2_t &operator=(const basic_domain_product2_t &&other) {
    if (this != &other) {
      this->_is_bottom = std::move(other._is_bottom);
      this->_first = std::move(other._first);
      this->_second = std::move(other._second);
    }
    return *this;
  }

  static basic_domain_product2_t top() {
    return basic_domain_product2_t(Domain1::top(), Domain2::top());
  }

  static basic_domain_product2_t bottom() {
    return basic_domain_product2_t(Domain1::bottom(), Domain2::bottom());
  }

  bool is_bottom() {
    this->canonicalize();
    return this->_is_bottom;
  }

  bool is_top() { return (this->_first.is_top() && this->_second.is_top()); }

  Domain1 &first(bool apply_reduction = true) {
    if (apply_reduction) {
      this->canonicalize();
    }
    return this->_first;
  }

  Domain2 &second(bool apply_reduction = true) {
    if (apply_reduction) {
      this->canonicalize();
    }
    return this->_second;
  }

  bool operator<=(basic_domain_product2_t other) {
    if (this->is_bottom()) {
      return true;
    } else if (other.is_bottom()) {
      return false;
    } else {
      return (this->_first <= other._first) && (this->_second <= other._second);
    }
  }

  bool operator==(basic_domain_product2_t other) {
    return (this->operator<=(other) && other.operator<=(*this));
  }

  void operator|=(basic_domain_product2_t other) {
    if (this->is_bottom()) {
      *this = other;
    } else if (other.is_bottom()) {
      return;
    } else {
      this->_first |= other._first;
      this->_second |= other._second;
    }
  }

  basic_domain_product2_t operator|(basic_domain_product2_t other) {
    if (this->is_bottom()) {
      return other;
    } else if (other.is_bottom()) {
      return *this;
    } else {
      return basic_domain_product2_t(this->_first | other._first,
                                     this->_second | other._second);
    }
  }

  basic_domain_product2_t operator||(basic_domain_product2_t other) {
    return basic_domain_product2_t(this->_first || other._first,
                                   this->_second || other._second,
                                   false /* do not apply reduction */);
  }

  basic_domain_product2_t operator&(basic_domain_product2_t other) {
    if (this->is_bottom() || other.is_bottom()) {
      return bottom();
    } else {
      return basic_domain_product2_t(this->_first & other._first,
                                     this->_second & other._second);
    }
  }

  basic_domain_product2_t operator&&(basic_domain_product2_t other) {
    if (this->is_bottom() || other.is_bottom()) {
      return bottom();
    } else {
      return basic_domain_product2_t(this->_first && other._first,
                                     this->_second && other._second);
    }
  }

  void write(crab::crab_os &o) {
    if (this->is_bottom()) {
      o << "_|_";
    } else {
      o << "(" << this->_first << ", " << this->_second << ")";
    }
  }

  static std::string getDomainName() {
    std::string name = "Product(" + Domain1::getDomainName() + "," +
                       Domain2::getDomainName() + ")";
    return name;
  }
}; // class basic_domain_product2

// Provided by Ikos
// Reduced product of two arbitrary domains with all operations.
template <typename Number, typename VariableName, typename Domain1,
          typename Domain2>
class domain_product2 final
    : public abstract_domain<
          domain_product2<Number, VariableName, Domain1, Domain2>> {
public:
  typedef domain_product2<Number, VariableName, Domain1, Domain2>
      domain_product2_t;
  typedef abstract_domain<domain_product2_t> abstract_domain_t;
  typedef Domain1 first_type;
  typedef Domain2 second_type;

  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;

private:
  typedef basic_domain_product2<Domain1, Domain2> basic_domain_product2_t;

  basic_domain_product2_t _product;

  domain_product2(basic_domain_product2_t &&product)
      : _product(std::move(product)) {}

  void reduce() {
    if (this->_product.first().is_bottom() ||
        this->_product.second().is_bottom()) {
      _product = basic_domain_product2_t::bottom();
    }
  }

public:
  void set_to_top() {
    domain_product2_t abs(basic_domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    domain_product2_t abs(basic_domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  domain_product2() : _product() {}

  domain_product2(Domain1 first, Domain2 second)
      : _product(basic_domain_product2_t(first, second)) {}

  domain_product2(const domain_product2_t &other) : _product(other._product) {}

  domain_product2(const domain_product2_t &&other)
      : _product(std::move(other._product)) {}

  domain_product2_t &operator=(const domain_product2_t &other) {
    if (this != &other)
      this->_product = other._product;
    return *this;
  }

  domain_product2_t &operator=(const domain_product2_t &&other) {
    if (this != &other)
      this->_product = std::move(other._product);
    return *this;
  }

  bool is_bottom() { return this->_product.is_bottom(); }

  bool is_top() { return this->_product.is_top(); }

  Domain1 &first() { return this->_product.first(); }

  Domain2 &second() { return this->_product.second(); }

  bool operator<=(domain_product2_t other) {
    return (this->_product <= other._product);
  }

  bool operator==(domain_product2_t other) {
    return (this->_product == other._product);
  }

  void operator|=(domain_product2_t other) { this->_product |= other._product; }

  domain_product2_t operator|(domain_product2_t other) {
    return domain_product2_t(this->_product | other._product);
  }

  domain_product2_t operator&(domain_product2_t other) {
    return domain_product2_t(this->_product & other._product);
  }

  domain_product2_t operator||(domain_product2_t other) {
    return domain_product2_t(this->_product || other._product);
  }

  domain_product2_t
  widening_thresholds(domain_product2_t other,
                      const iterators::thresholds<number_t> &ts) {
    bool apply_reduction = false;
    return domain_product2_t(basic_domain_product2_t(
        std::move(this->_product.first(apply_reduction)
                      .widening_thresholds(other._product.first(), ts)),
        std::move(this->_product.second(apply_reduction)
                      .widening_thresholds(other._product.second(), ts)),
        std::move(apply_reduction)));
  }

  domain_product2_t operator&&(domain_product2_t other) {
    return domain_product2_t(this->_product && other._product);
  }

  void assign(variable_t x, linear_expression_t e) {
    this->_product.first().assign(x, e);
    this->_product.second().assign(x, e);
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.first().apply(op, x, y, z);
    this->_product.second().apply(op, x, y, z);
    this->reduce();
  }

  void apply(operation_t op, variable_t x, variable_t y, Number k) {
    this->_product.first().apply(op, x, y, k);
    this->_product.second().apply(op, x, y, k);
    this->reduce();
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       domain_product2_t invariant) {
    this->_product.first().backward_assign(x, e, invariant.first());
    this->_product.second().backward_assign(x, e, invariant.second());
    this->reduce();
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, Number k,
                      domain_product2_t invariant) {
    this->_product.first().backward_apply(op, x, y, k, invariant.first());
    this->_product.second().backward_apply(op, x, y, k, invariant.second());
    this->reduce();
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      domain_product2_t invariant) {
    this->_product.first().backward_apply(op, x, y, z, invariant.first());
    this->_product.second().backward_apply(op, x, y, z, invariant.second());
    this->reduce();
  }

  void operator+=(linear_constraint_system_t csts) {
    this->_product.first() += csts;
    this->_product.second() += csts;
    this->reduce();
  }

  void operator-=(variable_t v) {
    this->_product.first() -= v;
    this->_product.second() -= v;
  }

  // cast operators

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    this->_product.first().apply(op, dst, src);
    this->_product.second().apply(op, dst, src);
    this->reduce();
  }

  // bitwise operators

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.first().apply(op, x, y, z);
    this->_product.second().apply(op, x, y, z);
    this->reduce();
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, Number k) {
    this->_product.first().apply(op, x, y, k);
    this->_product.second().apply(op, x, y, k);
    this->reduce();
  }

  // array operators

  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {
    this->_product.first().array_init(a, elem_size, lb_idx, ub_idx, val);
    this->_product.second().array_init(a, elem_size, lb_idx, ub_idx, val);
    this->reduce();
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {

    this->_product.first().array_load(lhs, a, elem_size, i);
    this->_product.second().array_load(lhs, a, elem_size, i);
    this->reduce();
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    this->_product.first().array_store(a, elem_size, i, val, is_strong_update);
    this->_product.second().array_store(a, elem_size, i, val, is_strong_update);
    this->reduce();
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t elem_size, linear_expression_t i,
                           linear_expression_t val,
                           bool is_strong_update) override {
    this->_product.first().array_store(a_new, a_old, elem_size, i, val,
                                       is_strong_update);
    this->_product.second().array_store(a_new, a_old, elem_size, i, val,
                                        is_strong_update);
    this->reduce();
  }

  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    this->_product.first().array_store_range(a, elem_size, i, j, val);
    this->_product.second().array_store_range(a, elem_size, i, j, val);
    this->reduce();
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    this->_product.first().array_store_range(a_new, a_old, elem_size, i, j,
                                             val);
    this->_product.second().array_store_range(a_new, a_old, elem_size, i, j,
                                              val);
    this->reduce();
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    this->_product.first().array_assign(lhs, rhs);
    this->_product.second().array_assign(lhs, rhs);
    this->reduce();
  }

  // backward array operations

  virtual void backward_array_init(variable_t a, linear_expression_t elem_size,
                                   linear_expression_t lb_idx,
                                   linear_expression_t ub_idx,
                                   linear_expression_t val,
                                   domain_product2_t invariant) override {
    this->_product.first().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                               val, invariant.first());
    this->_product.second().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                                val, invariant.second());
    this->reduce();
  }

  virtual void backward_array_load(variable_t lhs, variable_t a,
                                   linear_expression_t elem_size,
                                   linear_expression_t i,
                                   domain_product2_t invariant) override {

    this->_product.first().backward_array_load(lhs, a, elem_size, i,
                                               invariant.first());
    this->_product.second().backward_array_load(lhs, a, elem_size, i,
                                                invariant.second());
    this->reduce();
  }

  virtual void backward_array_store(variable_t a, linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    domain_product2_t invariant) override {
    this->_product.first().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.first());
    this->_product.second().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.second());
    this->reduce();
  }

  virtual void backward_array_store(variable_t a_new, variable_t a_old,
                                    linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    domain_product2_t invariant) override {
    this->_product.first().backward_array_store(
        a_new, a_old, elem_size, i, val, is_strong_update, invariant.first());
    this->_product.second().backward_array_store(
        a_new, a_old, elem_size, i, val, is_strong_update, invariant.second());
    this->reduce();
  }

  virtual void
  backward_array_store_range(variable_t a, linear_expression_t elem_size,
                             linear_expression_t i, linear_expression_t j,
                             linear_expression_t val,
                             domain_product2_t invariant) override {
    this->_product.first().backward_array_store_range(a, elem_size, i, j, val,
                                                      invariant.first());
    this->_product.second().backward_array_store_range(a, elem_size, i, j, val,
                                                       invariant.second());
    this->reduce();
  }

  virtual void backward_array_store_range(
      variable_t a_new, variable_t a_old, linear_expression_t elem_size,
      linear_expression_t i, linear_expression_t j, linear_expression_t val,
      domain_product2_t invariant) override {
    this->_product.first().backward_array_store_range(
        a_new, a_old, elem_size, i, j, val, invariant.first());
    this->_product.second().backward_array_store_range(
        a_new, a_old, elem_size, i, j, val, invariant.second());
    this->reduce();
  }

  virtual void backward_array_assign(variable_t lhs, variable_t rhs,
                                     domain_product2_t invariant) override {
    this->_product.first().backward_array_assign(lhs, rhs, invariant.first());
    this->_product.second().backward_array_assign(lhs, rhs, invariant.second());
    this->reduce();
  }

  // pointer operators
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    this->_product.first().pointer_load(lhs, rhs);
    this->_product.second().pointer_load(lhs, rhs);
    this->reduce();
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    this->_product.first().pointer_store(lhs, rhs);
    this->_product.second().pointer_store(lhs, rhs);
    this->reduce();
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    this->_product.first().pointer_assign(lhs, rhs, offset);
    this->_product.second().pointer_assign(lhs, rhs, offset);
    this->reduce();
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    this->_product.first().pointer_mk_obj(lhs, address);
    this->_product.second().pointer_mk_obj(lhs, address);
    this->reduce();
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    this->_product.first().pointer_function(lhs, func);
    this->_product.second().pointer_function(lhs, func);
    this->reduce();
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    this->_product.first().pointer_mk_null(lhs);
    this->_product.second().pointer_mk_null(lhs);
    this->reduce();
  }

  virtual void pointer_assume(pointer_constraint_t cst) override {
    this->_product.first().pointer_assume(cst);
    this->_product.second().pointer_assume(cst);
    this->reduce();
  }

  virtual void pointer_assert(pointer_constraint_t cst) override {
    this->_product.first().pointer_assert(cst);
    this->_product.second().pointer_assert(cst);
    this->reduce();
  }

  // boolean operators
  virtual void assign_bool_cst(variable_t lhs,
                               linear_constraint_t rhs) override {
    this->_product.first().assign_bool_cst(lhs, rhs);
    this->_product.second().assign_bool_cst(lhs, rhs);
    this->reduce();
  }

  virtual void assign_bool_var(variable_t lhs, variable_t rhs,
                               bool is_not_rhs) override {
    this->_product.first().assign_bool_var(lhs, rhs, is_not_rhs);
    this->_product.second().assign_bool_var(lhs, rhs, is_not_rhs);
    this->reduce();
  }

  virtual void apply_binary_bool(bool_operation_t op, variable_t x,
                                 variable_t y, variable_t z) override {
    this->_product.first().apply_binary_bool(op, x, y, z);
    this->_product.second().apply_binary_bool(op, x, y, z);
    this->reduce();
  }

  virtual void assume_bool(variable_t v, bool is_negated) override {
    this->_product.first().assume_bool(v, is_negated);
    this->_product.second().assume_bool(v, is_negated);
    this->reduce();
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                        domain_product2_t inv) {
    this->_product.first().backward_assign_bool_cst(lhs, rhs, inv.first());
    this->_product.second().backward_assign_bool_cst(lhs, rhs, inv.second());
    this->reduce();
  }

  virtual void backward_assign_bool_var(variable_t lhs, variable_t rhs,
                                        bool is_not_rhs,
                                        domain_product2_t inv) {
    this->_product.first().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                    inv.first());
    this->_product.second().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                     inv.second());
    this->reduce();
  }

  virtual void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                          variable_t y, variable_t z,
                                          domain_product2_t inv) {
    this->_product.first().backward_apply_binary_bool(op, x, y, z, inv.first());
    this->_product.second().backward_apply_binary_bool(op, x, y, z,
                                                       inv.second());
    this->reduce();
  }

  virtual void forget(const variable_vector_t &variables) {
    this->_product.first().forget(variables);
    this->_product.second().forget(variables);
  }

  virtual void project(const variable_vector_t &variables) {
    this->_product.first().project(variables);
    this->_product.second().project(variables);
  }

  virtual void expand(variable_t var, variable_t new_var) {
    this->_product.first().expand(var, new_var);
    this->_product.second().expand(var, new_var);
  }

  virtual void normalize() {
    this->_product.first().normalize();
    this->_product.second().normalize();
  }

  virtual void minimize() {
    this->_product.first().minimize();
    this->_product.second().minimize();
  }

  virtual linear_constraint_system_t to_linear_constraint_system() {
    linear_constraint_system_t csts;
    // XXX: We might add redundant constraints.
    csts += this->_product.first().to_linear_constraint_system();
    csts += this->_product.second().to_linear_constraint_system();
    return csts;
  }

  virtual disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    disjunctive_linear_constraint_system_t csts;
    // XXX: We might add redundant constraints.
    csts += this->_product.first().to_disjunctive_linear_constraint_system();
    csts += this->_product.second().to_disjunctive_linear_constraint_system();
    return csts;
  }

  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) override {
    this->_product.first().rename(from, to);
    this->_product.second().rename(from, to);
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    this->_product.first().intrinsic(name, inputs, outputs);
    this->_product.second().intrinsic(name, inputs, outputs);    
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  domain_product2_t invariant) override {
    this->_product.first().backward_intrinsic(name, inputs, outputs, invariant.first());
    this->_product.second().backward_intrinsic(name, inputs, outputs, invariant.second());        
  }
  /* end intrinsics operations */
  
  void write(crab::crab_os &o) { this->_product.write(o); }

  static std::string getDomainName() {
    return basic_domain_product2_t::getDomainName();
  }

}; // class domain_product2

namespace reduced_product_impl {
class default_params {
public:
  enum { left_propagate_equalities = 1 };
  enum { right_propagate_equalities = 1 };
  enum { left_propagate_inequalities = 1 };
  enum { right_propagate_inequalities = 1 };
  enum { left_propagate_intervals = 1 };
  enum { right_propagate_intervals = 1 };
  enum { disable_reduction = 0 };
  enum { apply_reduction_only_add_constraint = 0 };
};

class term_dbm_params {
public:
  enum { left_propagate_equalities = 1 };
  enum { right_propagate_equalities = 1 };
  enum { left_propagate_inequalities = 0 };
  enum { right_propagate_inequalities = 0 };
  enum { left_propagate_intervals = 0 };
  enum { right_propagate_intervals = 0 };
  enum { disable_reduction = 0 };
  enum { apply_reduction_only_add_constraint = 1 };
};
} // namespace reduced_product_impl

// This domain is similar to domain_product2 but it uses a more
// precise reduction operation.
template <typename Domain1, typename Domain2,
          class Params = reduced_product_impl::default_params>
class reduced_numerical_domain_product2 final
    : public abstract_domain<
          reduced_numerical_domain_product2<Domain1, Domain2, Params>> {

public:
  typedef reduced_numerical_domain_product2<Domain1, Domain2, Params>
      reduced_numerical_domain_product2_t;
  typedef abstract_domain<reduced_numerical_domain_product2_t>
      abstract_domain_t;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  // Assume that Domain1 and Domain2 have the same types for
  // number_t and varname_t
  typedef typename Domain1::number_t number_t;
  typedef typename Domain1::varname_t varname_t;
  typedef interval<number_t> interval_t;

private:
  typedef patricia_tree_set<variable_t> variable_set_t;
  typedef domain_product2<number_t, varname_t, Domain1, Domain2>
      domain_product2_t;

  domain_product2_t _product;

  reduced_numerical_domain_product2(const domain_product2_t &product)
      : _product(product) {}

  linear_constraint_system_t to_linear_constraints(variable_t v,
                                                   interval_t i) const {
    linear_constraint_system_t csts;
    if (i.lb().is_finite() && i.ub().is_finite()) {
      auto lb = *(i.lb().number());
      auto ub = *(i.ub().number());
      if (lb == ub) {
        csts += (v == lb);
      } else {
        csts += (v >= lb);
        csts += (v <= ub);
      }
    } else if (i.lb().is_finite()) {
      auto lb = *(i.lb().number());
      csts += (v >= lb);
    } else if (i.ub().is_finite()) {
      auto ub = *(i.ub().number());
      csts += (v <= ub);
    }
    return csts;
  }

  void reduce_variable(const variable_t &v) {
    crab::CrabStats::count(getDomainName() + ".count.reduce");
    crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

    if (!is_bottom() && !Params::disable_reduction) {

      // We just propagate from one domain to another.  We could
      // propagate in the other direction ... and repeat it
      // computing a fixpoint+narrowing of descending iterations.

      Domain1 &inv1 = _product.first();
      Domain2 &inv2 = _product.second();

      //////
      // propagate interval constraints between domains
      //////
      if (Params::left_propagate_intervals) {
        interval_t i1 = inv1[v];
        if (!i1.is_top())
          inv2 += to_linear_constraints(v, i1);
      }
      if (Params::right_propagate_intervals) {
        interval_t i2 = inv2[v];
        if (!i2.is_top())
          inv1 += to_linear_constraints(v, i2);
      }

      //////
      // propagate other constraints expressed by the domains
      //////
      if ((Params::left_propagate_equalities ||
           Params::left_propagate_inequalities) &&
          (Params::right_propagate_equalities ||
           Params::right_propagate_inequalities)) {
        linear_constraint_system_t csts1, csts2, filtered_csts1, filtered_csts2;
        bool propagate_only_equalities;

        propagate_only_equalities = !Params::left_propagate_inequalities;
        crab::domains::reduced_domain_traits<Domain1>::extract(
            inv1, v, csts1, propagate_only_equalities);

        propagate_only_equalities = !Params::right_propagate_inequalities;
        crab::domains::reduced_domain_traits<Domain2>::extract(
            inv2, v, csts2, propagate_only_equalities);

        // filter out those redundant constraints (i.e.,
        // constraints that the other domain already knows about)

        for (auto &c1 : csts1) {
          if (std::find_if(csts2.begin(), csts2.end(),
                           [c1](const linear_constraint_t &c2) {
                             return c2.equal(c1);
                           }) == csts2.end()) {
            filtered_csts1 += c1;
          }
        }

        {
          std::string k(getDomainName() + ".count.reduce.equalities_from_" +
                        _product.first().getDomainName());
          crab::CrabStats::uset(k, crab::CrabStats::get(k) +
                                       filtered_csts1.size());
        }
        inv2 += filtered_csts1;

        for (auto &c2 : csts2) {
          if (std::find_if(csts1.begin(), csts1.end(),
                           [c2](const linear_constraint_t &c1) {
                             return c1.equal(c2);
                           }) == csts1.end()) {
            filtered_csts2 += c2;
          }
        }
        {
          std::string k(getDomainName() + ".count.reduce.equalities_from_" +
                        _product.second().getDomainName());
          crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts2.size());
        }
        inv1 += filtered_csts2;
      } else if (Params::left_propagate_equalities ||
                 Params::left_propagate_inequalities) {
        linear_constraint_system_t csts1;
        const bool propagate_only_equalities =
            !Params::left_propagate_inequalities;
        crab::domains::reduced_domain_traits<Domain1>::extract(
            inv1, v, csts1, propagate_only_equalities);
        std::string k(getDomainName() + ".count.reduce.equalities_from_" +
                      _product.first().getDomainName());
        crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts1.size());
        inv2 += csts1;
      } else if (Params::right_propagate_equalities ||
                 Params::right_propagate_inequalities) {
        linear_constraint_system_t csts2;
        const bool propagate_only_equalities =
            !Params::right_propagate_inequalities;
        crab::domains::reduced_domain_traits<Domain2>::extract(
            inv2, v, csts2, propagate_only_equalities);
        std::string k(getDomainName() + ".count.reduce.equalities_from_" +
                      _product.second().getDomainName());
        crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts2.size());
        inv1 += csts2;
      }
    }
  }

public:
  void set_to_top() {
    reduced_numerical_domain_product2_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    reduced_numerical_domain_product2_t abs(domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  reduced_numerical_domain_product2() : _product() {}

  reduced_numerical_domain_product2(
      const reduced_numerical_domain_product2_t &other)
      : _product(other._product) {}

  reduced_numerical_domain_product2_t &
  operator=(const reduced_numerical_domain_product2_t &other) {
    if (this != &other)
      this->_product = other._product;

    return *this;
  }

  bool is_bottom() { return this->_product.is_bottom(); }

  bool is_top() { return this->_product.is_top(); }

  Domain1 &first() { return this->_product.first(); }

  Domain2 &second() { return this->_product.second(); }

  bool operator<=(reduced_numerical_domain_product2_t other) {
    return this->_product <= other._product;
  }

  void operator|=(reduced_numerical_domain_product2_t other) {
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";);
    this->_product |= other._product;
    CRAB_LOG("reduced-dom", crab::outs() << *this << "\n----------------\n";);
  }

  reduced_numerical_domain_product2_t
  operator|(reduced_numerical_domain_product2_t other) {
    reduced_numerical_domain_product2_t res(this->_product | other._product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";
             crab::outs() << res << "\n================\n");
    return res;
  }

  reduced_numerical_domain_product2_t
  operator&(reduced_numerical_domain_product2_t other) {
    reduced_numerical_domain_product2_t res(this->_product & other._product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ MEET ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator||(reduced_numerical_domain_product2_t other) {
    reduced_numerical_domain_product2_t res(this->_product || other._product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  widening_thresholds(reduced_numerical_domain_product2_t other,
                      const iterators::thresholds<number_t> &ts) {
    reduced_numerical_domain_product2_t res(
        this->_product.widening_thresholds(other._product, ts));
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator&&(reduced_numerical_domain_product2_t other) {
    reduced_numerical_domain_product2_t res(this->_product && other._product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ NARROWING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  void set(variable_t v, interval_t x) {
    this->_product.first().set(v, x);
    this->_product.second().set(v, x);
  }

  interval_t operator[](variable_t v) {
    // We can choose either first or second domain
    return this->_product.second()[v];
  }

  void operator+=(linear_constraint_system_t csts) {
    this->_product += csts;
    for (auto v : csts.variables()) {
      reduce_variable(v);
    }
    CRAB_LOG("combined-domain", crab::outs() << "Added constraints " << csts
                                             << "=" << *this << "\n");
  }

  void operator-=(variable_t v) { this->_product -= v; }

  void assign(variable_t x, linear_expression_t e) {
    this->_product.assign(x, e);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain", crab::outs()
                                    << x << ":=" << e << "=" << *this << "\n");
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << z << "=" << *this << "\n");
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << k << "=" << *this << "\n");
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       reduced_numerical_domain_product2_t invariant) {
    this->_product.backward_assign(x, e, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      for (auto v : e.variables())
        this->reduce_variable(v);
    }
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t k,
                      reduced_numerical_domain_product2_t invariant) {
    this->_product.backward_apply(op, x, y, k, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      this->reduce_variable(y);
    }
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      reduced_numerical_domain_product2_t invariant) {
    this->_product.backward_apply(op, x, y, z, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      this->reduce_variable(y);
      this->reduce_variable(z);
    }
  }

  // cast operators

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    this->_product.apply(op, dst, src);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(dst);
    }
  }

  // bitwise operators

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
  }

  /*
     Begin unimplemented operations

     reduced_numerical_domain_product2 implements only standard
     abstract operations of a numerical domain.  The
     implementation of boolean, array, or pointer operations is
     empty because they should never be called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(crab::domains::bool_operation_t op, variable_t x,
                         variable_t y, variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                reduced_numerical_domain_product2_t invariant) {
  }
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                reduced_numerical_domain_product2_t invariant) {
  }
  void
  backward_apply_binary_bool(crab::domains::bool_operation_t op, variable_t x,
                             variable_t y, variable_t z,
                             reduced_numerical_domain_product2_t invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t val) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t val) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           reduced_numerical_domain_product2_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           reduced_numerical_domain_product2_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            reduced_numerical_domain_product2_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            reduced_numerical_domain_product2_t invariant) {}
  void
  backward_array_store_range(variable_t a, linear_expression_t elem_size,
                             linear_expression_t i, linear_expression_t j,
                             linear_expression_t val,
                             reduced_numerical_domain_product2_t invariant) {}
  void backward_array_store_range(
      variable_t a_new, variable_t a_old, linear_expression_t elem_size,
      linear_expression_t i, linear_expression_t j, linear_expression_t val,
      reduced_numerical_domain_product2_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             reduced_numerical_domain_product2_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs) {}
  void pointer_store(variable_t lhs, variable_t rhs) {}
  void pointer_assign(variable_t lhs, variable_t rhs,
                      linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    this->_product.rename(from, to);
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    this->_product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  reduced_numerical_domain_product2_t invariant) override {
    this->_product.backward_intrinsic(name, inputs, outputs, invariant._product);
  }
  /* end intrinsics operations */
  
  void forget(const variable_vector_t &variables) {
    this->_product.forget(variables);
  }

  void project(const variable_vector_t &variables) {
    this->_product.project(variables);
  }

  void expand(variable_t var, variable_t new_var) {
    this->_product.expand(var, new_var);
  }

  void normalize() { this->_product.normalize(); }

  void minimize() { this->_product.minimize(); }

  void write(crab_os &o) { this->_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    return this->_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    return this->_product.to_disjunctive_linear_constraint_system();
  }

  static std::string getDomainName() {
    std::string name = "ReducedProduct(" + Domain1::getDomainName() + "," +
                       Domain2::getDomainName() + ")";
    return name;
  }

}; // class reduced_numerical_domain_product2

/*
 *  The reduce operator based on "Static Analysis of Arithmetical
 *  Congruences" by P. Granger published in International Journal of
 *  Computer Mathematics, 1989.
 */
template <typename Number> class interval_congruence {
public:
  typedef interval_congruence<Number> interval_congruence_t;

private:
  typedef interval<Number> interval_t;
  typedef bound<Number> bound_t;
  typedef congruence<Number> congruence_t;

private:
  interval_t _first;
  congruence_t _second;

private:
  interval_congruence(bool is_bottom)
      : _first(is_bottom ? interval_t::bottom() : interval_t::top()),
        _second(is_bottom ? congruence_t::bottom() : congruence_t::top()) {}

public:
  static interval_congruence_t top() { return interval_congruence(false); }

  static interval_congruence_t bottom() { return interval_congruence(true); }

private:
  inline Number abs(Number x) { return x < 0 ? -x : x; }

  // operator % can return a negative number
  // mod(a, b) always returns a positive number
  inline Number mod(Number a, Number b) {
    Number m = a % b;
    if (m < 0)
      return m + b;
    else
      return m;
  }

  // R(c,a) is the least element of c greater or equal than a
  inline Number R(congruence_t c, Number a) {
    Number m = c.get_modulo();
    Number p = c.get_remainder();
    return a + mod(p - a, abs(m));
  }

  // L(c,a) is the greatest element of c smaller or equal than a
  inline Number L(congruence_t c, Number a) {
    Number m = c.get_modulo();
    Number p = c.get_remainder();
    return a - mod(a - p, abs(m));
  }

public:
  interval_congruence(Number n)
      : _first(interval_t(n)), _second(congruence_t(n)) {}

  interval_congruence(interval_t i, congruence_t c) : _first(i), _second(c) {
    this->reduce();
  }

  interval_congruence(interval_t i) : _first(i), _second(congruence_t::top()) {
    this->reduce();
  }

  interval_congruence(congruence_t c) : _first(interval_t::top()), _second(c) {
    this->reduce();
  }

  interval_congruence(const interval_congruence &other)
      : _first(other._first), _second(other._second) {}

  interval_congruence_t &operator=(interval_congruence_t other) {
    this->_first = other._first;
    this->_second = other._second;
    return *this;
  }

  bool is_bottom() {
    return this->_first.is_bottom() || this->_second.is_bottom();
  }

  bool is_top() { return this->_first.is_top() && this->_second.is_top(); }

  interval_t &first() { return this->_first; }

  congruence_t &second() { return this->_second; }

  /*
     Let (i,c) be a pair of interval and congruence these are the
     main rules described by Granger:

     if (c.is_bottom() || i.is_bottom()) (bottom(), bottom());
     if (c = 0Z+a and a \notin i)        (bottom(), bottom());
     if (c = 0Z+a)                       ([a,a]   , c);
     if (i=[a,b] and R(c,a) > L(c,b))    (bottom(), bottom());
     if (i=[a,b])                        ([R(c,a), L(c,b)], c);
     if (i=[a,+oo])                      ([R(c,a), +oo], c);
     if (i=[-oo,b])                      ([-oo, L(c,b)], c);
     otherwise                           (i,c)
   */

  void reduce() {
    interval_t &i = first();
    congruence_t &c = second();

    if (i.is_bottom() || c.is_bottom()) {
      i = interval_t::bottom();
      c = congruence_t::bottom();
    }

    // congruence is top and interval is a singleton
    if (c.is_top()) {
      boost::optional<Number> n = i.singleton();
      if (n) {
        c = congruence_t(*n);
      }
      return;
    }

    Number modulo = c.get_modulo();
    if (modulo == 0) {
      // congruence is a singleton so we refine the interval
      interval_t a(c.get_remainder());
      if (!(a <= i)) {
        i = interval_t::bottom();
        c = congruence_t::bottom();
      } else {
        i = a;
      }
    } else {
      // refine lower and upper bounds of the interval using
      // congruences
      bound_t lb = i.lb();
      bound_t ub = i.ub();

      if (lb.is_finite() && ub.is_finite()) {
        Number x = R(c, *(lb.number()));
        Number y = L(c, *(ub.number()));
        if (x > y) {
          i = interval_t::bottom();
          c = congruence_t::bottom();
        } else if (x == y) {
          i = interval_t(x);
          c = congruence_t(x);
        } else {
          i = interval_t(bound_t(x), bound_t(y));
        }
      } else if (lb.is_finite()) {
        Number x = R(c, *(lb.number()));
        i = interval_t(bound_t(x), bound_t::plus_infinity());
      } else if (ub.is_finite()) {
        Number y = L(c, *(ub.number()));
        i = interval_t(bound_t::minus_infinity(), bound_t(y));
      } else {
        // interval is top
      }
    }
  }

  void write(crab_os &o) const {
    o << "(" << this->_first << ", " << this->_second << ")";
  }

public:
  interval_congruence_t operator+(interval_congruence_t x) {
    return interval_congruence_t(this->_first.operator+(x.first()),
                                 this->_second.operator+(x.second()));
  }

  interval_congruence_t operator-(interval_congruence_t x) {
    return interval_congruence_t(this->_first.operator-(x.first()),
                                 this->_second.operator-(x.second()));
  }

  interval_congruence_t operator*(interval_congruence_t x) {
    return interval_congruence_t(this->_first.operator*(x.first()),
                                 this->_second.operator*(x.second()));
  }

  interval_congruence_t operator/(interval_congruence_t x) {
    return interval_congruence_t(this->_first.operator/(x.first()),
                                 this->_second.operator/(x.second()));
  }

  interval_congruence_t operator|(interval_congruence_t other) {
    return interval_congruence_t(this->_first | other._first,
                                 this->_second | other._second);
  }

  interval_congruence_t operator&(interval_congruence_t other) {
    return interval_congruence_t(this->_first & other._first,
                                 this->_second & other._second);
  }

public:
  // division and remainder operations

  interval_congruence_t SDiv(interval_congruence_t x) {
    return interval_congruence_t(this->_first.SDiv(x.first()),
                                 this->_second.SDiv(x.second()));
  }

  interval_congruence_t UDiv(interval_congruence_t x) {
    return interval_congruence_t(this->_first.UDiv(x.first()),
                                 this->_second.UDiv(x.second()));
  }

  interval_congruence_t SRem(interval_congruence_t x) {
    return interval_congruence_t(this->_first.SRem(x.first()),
                                 this->_second.SRem(x.second()));
  }

  interval_congruence_t URem(interval_congruence_t x) {
    return interval_congruence_t(this->_first.URem(x.first()),
                                 this->_second.URem(x.second()));
  }

  // bitwise operations

  interval_congruence_t Trunc(unsigned width) {
    return interval_congruence_t(this->_first.Trunc(width),
                                 this->_second.Trunc(width));
  }

  interval_congruence_t ZExt(unsigned width) {
    return interval_congruence_t(this->_first.ZExt(width),
                                 this->_second.ZExt(width));
  }

  interval_congruence_t SExt(unsigned width) {
    return interval_congruence_t(this->_first.SExt(width),
                                 this->_second.SExt(width));
  }

  interval_congruence_t And(interval_congruence_t x) {
    return interval_congruence_t(this->_first.And(x.first()),
                                 this->_second.And(x.second()));
  }

  interval_congruence_t Or(interval_congruence_t x) {
    return interval_congruence_t(this->_first.Or(x.first()),
                                 this->_second.Or(x.second()));
  }

  interval_congruence_t Xor(interval_congruence_t x) {
    return interval_congruence_t(this->_first.Xor(x.first()),
                                 this->_second.Xor(x.second()));
  }

  interval_congruence_t Shl(interval_congruence_t x) {
    return interval_congruence_t(this->_first.Shl(x.first()),
                                 this->_second.Shl(x.second()));
  }

  interval_congruence_t LShr(interval_congruence_t x) {
    return interval_congruence_t(this->_first.LShr(x.first()),
                                 this->_second.LShr(x.second()));
  }

  interval_congruence_t AShr(interval_congruence_t x) {
    return interval_congruence_t(this->_first.AShr(x.first()),
                                 this->_second.AShr(x.second()));
  }
};

template <typename Number>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const interval_congruence<Number> &v) {
  v.write(o);
  return o;
}

// Reduced product of a numerical domain with congruences.
template <typename NumAbsDom>
class numerical_congruence_domain final
    : public abstract_domain<numerical_congruence_domain<NumAbsDom>> {

  typedef numerical_congruence_domain<NumAbsDom> rnc_domain_t;
  typedef abstract_domain<rnc_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef typename NumAbsDom::number_t number_t;
  typedef typename NumAbsDom::varname_t varname_t;

  typedef congruence_domain<number_t, varname_t> congruence_domain_t;
  typedef interval_congruence<number_t> interval_congruence_t;
  typedef interval<number_t> interval_t;

private:
  typedef patricia_tree_set<variable_t> variable_set_t;
  typedef domain_product2<number_t, varname_t, NumAbsDom, congruence_domain_t>
      domain_product2_t;

  domain_product2_t _product;

  numerical_congruence_domain(const domain_product2_t &product)
      : _product(product) {}

  void reduce_variable(const variable_t &v) {
    crab::CrabStats::count(getDomainName() + ".count.reduce");
    crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

    if (is_bottom())
      return;

    auto i = this->_product.first()[v]; // project on intervals
    auto c = this->_product.second()[v];
    interval_congruence_t val(i, c);

    if (val.is_bottom()) {
      set_to_bottom();
    } else {
      if (val.first() != i) {
        // FIXME: this is imprecise for relational domains
        this->_product.first().set(v, val.first());
      }

      if (val.second() != c)
        this->_product.second().set(v, val.second());
    }
  }

  void reduce_variables(variable_set_t variables) {
    for (typename variable_set_t::iterator it = variables.begin();
         !is_bottom() && it != variables.end(); ++it) {
      this->reduce_variable((*it));
    }
  }

public:
  void set_to_top() {
    rnc_domain_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    rnc_domain_t abs(domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  numerical_congruence_domain() : _product() {}

  numerical_congruence_domain(const rnc_domain_t &other)
      : _product(other._product) {}

  rnc_domain_t &operator=(const rnc_domain_t &other) {
    if (this != &other)
      this->_product = other._product;

    return *this;
  }

  bool is_bottom() { return this->_product.is_bottom(); }

  bool is_top() { return this->_product.is_top(); }

  NumAbsDom &first() { return this->_product.first(); }

  congruence_domain_t &second() { return this->_product.second(); }

  bool operator<=(rnc_domain_t other) {
    return this->_product <= other._product;
  }

  bool operator==(rnc_domain_t other) {
    return this->_product == other._product;
  }

  void operator|=(rnc_domain_t other) { this->_product |= other._product; }

  rnc_domain_t operator|(rnc_domain_t other) {
    return rnc_domain_t(this->_product | other._product);
  }

  rnc_domain_t operator&(rnc_domain_t other) {
    return rnc_domain_t(this->_product & other._product);
  }

  rnc_domain_t operator||(rnc_domain_t other) {
    return rnc_domain_t(this->_product || other._product);
  }

  rnc_domain_t widening_thresholds(rnc_domain_t other,
                                   const iterators::thresholds<number_t> &ts) {
    return rnc_domain_t(this->_product.widening_thresholds(other._product, ts));
  }

  rnc_domain_t operator&&(rnc_domain_t other) {
    return rnc_domain_t(this->_product && other._product);
  }

  // pre: x is already reduced
  void set(variable_t v, interval_congruence_t x) {
    this->_product.first().set(v, x.first());
    this->_product.second().set(v, x.second());
  }

  interval_congruence_t get(variable_t v) {
    return interval_congruence_t(this->_product.first()[v],
                                 this->_product.second()[v]);
  }

  interval_t operator[](variable_t v) {
    interval_congruence_t x = get(v);
    return x.first();
  }

  void operator+=(linear_constraint_system_t csts) {
    this->_product += csts;
    this->reduce_variables(csts.variables());
  }

  void operator-=(variable_t v) { this->_product -= v; }

  void assign(variable_t x, linear_expression_t e) {
    this->_product.assign(x, e);
    this->reduce_variable(x);
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
    this->reduce_variable(x);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
    this->reduce_variable(x);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       rnc_domain_t invariant) {
    this->_product.backward_assign(x, e, invariant._product);
    // reduce the variables in the right-hand side
    for (auto v : e.variables())
      this->reduce_variable(v);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t k,
                      rnc_domain_t invariant) {
    this->_product.backward_apply(op, x, y, k, invariant._product);
    // reduce the variables in the right-hand side
    this->reduce_variable(y);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      rnc_domain_t invariant) {
    this->_product.backward_apply(op, x, y, z, invariant._product);
    // reduce the variables in the right-hand side
    this->reduce_variable(y);
    this->reduce_variable(z);
  }

  // cast operators

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    this->_product.apply(op, dst, src);
    this->reduce_variable(dst);
  }

  // bitwise operators

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
    this->reduce_variable(x);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
    this->reduce_variable(x);
  }

  /*
     Begin unimplemented operations

     numerical_congruence_domain implements only standard abstract
     operations of a numerical domain.  The implementation of
     boolean, array, or pointer operations is empty because they
     should never be called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(crab::domains::bool_operation_t op, variable_t x,
                         variable_t y, variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                rnc_domain_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                rnc_domain_t invariant) {}
  void backward_apply_binary_bool(crab::domains::bool_operation_t op,
                                  variable_t x, variable_t y, variable_t z,
                                  rnc_domain_t invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t val) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t val) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           rnc_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           rnc_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, rnc_domain_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, rnc_domain_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t val,
                                  rnc_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t val,
                                  rnc_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             rnc_domain_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs) {}
  void pointer_store(variable_t lhs, variable_t rhs) {}
  void pointer_assign(variable_t lhs, variable_t rhs,
                      linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void forget(const variable_vector_t &variables) {
    this->_product.forget(variables);
  }

  void project(const variable_vector_t &variables) {
    this->_product.project(variables);
  }

  void expand(variable_t var, variable_t new_var) {
    this->_product.expand(var, new_var);
  }

  void normalize() { this->_product.normalize(); }

  void minimize() { this->_product.minimize(); }

  void write(crab_os &o) { this->_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    return this->_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    return this->_product.to_disjunctive_linear_constraint_system();
  }

  static std::string getDomainName() {
    return domain_product2_t::getDomainName();
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    this->_product.rename(from, to);    
  }


  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    this->_product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			   rnc_domain_t invariant) override {
    this->_product.backward_intrinsic(name, inputs, outputs, invariant._product);
  }
  /* end intrinsics operations */
  
}; // class numerical_congruence_domain

// Reduced product of a numerical domain with nullity
// Reduction is done by as in domain_product2.
template <typename NumAbsDom>
class numerical_nullity_domain final
    : public abstract_domain<numerical_nullity_domain<NumAbsDom>> {

  typedef numerical_nullity_domain<NumAbsDom> nn_domain_t;
  typedef abstract_domain<nn_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef typename NumAbsDom::number_t number_t;
  typedef typename NumAbsDom::varname_t varname_t;
  typedef interval<number_t> interval_t;

private:
  typedef nullity_domain<number_t, varname_t> nullity_domain_t;
  typedef domain_product2<number_t, varname_t, NumAbsDom, nullity_domain_t>
      domain_product2_t;

  domain_product2_t _product;

  numerical_nullity_domain(const domain_product2_t &product)
      : _product(product) {}

public:
  void set_to_top() {
    nn_domain_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    nn_domain_t abs(domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  numerical_nullity_domain() : _product() {}

  numerical_nullity_domain(const nn_domain_t &other)
      : _product(other._product) {}

  nn_domain_t &operator=(const nn_domain_t &other) {
    if (this != &other)
      this->_product = other._product;

    return *this;
  }

  bool is_bottom() { return this->_product.is_bottom(); }

  bool is_top() { return this->_product.is_top(); }

  NumAbsDom &first() { return this->_product.first(); }

  nullity_domain_t &second() { return this->_product.second(); }

  bool operator<=(nn_domain_t other) {
    return this->_product <= other._product;
  }

  bool operator==(nn_domain_t other) {
    return this->_product == other._product;
  }

  void operator|=(nn_domain_t other) { this->_product |= other._product; }

  nn_domain_t operator|(nn_domain_t other) {
    return nn_domain_t(this->_product | other._product);
  }

  nn_domain_t operator&(nn_domain_t other) {
    return nn_domain_t(this->_product & other._product);
  }

  nn_domain_t operator||(nn_domain_t other) {
    return nn_domain_t(this->_product || other._product);
  }

  nn_domain_t widening_thresholds(nn_domain_t other,
                                  const iterators::thresholds<number_t> &ts) {
    return nn_domain_t(this->_product.widening_thresholds(other._product, ts));
  }

  nn_domain_t operator&&(nn_domain_t other) {
    return nn_domain_t(this->_product && other._product);
  }

  // numerical_domains

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
  }

  void assign(variable_t x, linear_expression_t e) {
    this->_product.assign(x, e);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       nn_domain_t invariant) {
    this->_product.backward_assign(x, e, invariant._product);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t k,
                      nn_domain_t invariant) {
    this->_product.backward_apply(op, x, y, k, invariant._product);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      nn_domain_t invariant) {
    this->_product.backward_apply(op, x, y, z, invariant._product);
  }

  void operator+=(linear_constraint_system_t csts) { this->_product += csts; }

  void operator-=(variable_t v) { this->_product -= v; }

  void set(variable_t v, interval_t intv) {
    this->_product.first().set(v, intv);
  }

  interval_t operator[](variable_t v) { return this->_product.first()[v]; }

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {
    this->_product.assign_bool_cst(lhs, rhs);
  }

  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {
    this->_product.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {
    this->_product.apply_binary_bool(op, x, y, z);
  }

  void assume_bool(variable_t v, bool is_negated) {
    this->_product.assume_bool(v, is_negated);
  }

  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                nn_domain_t invariant) {
    this->_product.backward_assign_bool_cst(lhs, rhs, invariant._product);
  }

  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                nn_domain_t invariant) {
    this->_product.backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                            invariant._product);
  }

  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  nn_domain_t invariant) {
    this->_product.backward_apply_binary_bool(op, x, y, z, invariant._product);
  }

  // cast operators

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    this->_product.apply(op, dst, src);
  }

  // bitwise operators

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    this->_product.apply(op, x, y, z);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    this->_product.apply(op, x, y, k);
  }

  // array operators

  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {
    this->_product.array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {
    this->_product.array_load(lhs, a, elem_size, i);
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    this->_product.array_store(a, elem_size, i, val, is_strong_update);
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t elem_size, linear_expression_t i,
                           linear_expression_t val,
                           bool is_strong_update) override {
    this->_product.array_store(a_new, a_old, elem_size, i, val,
                               is_strong_update);
  }

  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    this->_product.array_store_range(a, elem_size, i, j, val);
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    this->_product.array_store_range(a_new, a_old, elem_size, i, j, val);
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    this->_product.array_assign(lhs, rhs);
  }

  // backward array operators

  virtual void backward_array_init(variable_t a, linear_expression_t elem_size,
                                   linear_expression_t lb_idx,
                                   linear_expression_t ub_idx,
                                   linear_expression_t val,
                                   nn_domain_t invariant) override {
    this->_product.backward_array_init(a, elem_size, lb_idx, ub_idx, val,
                                       invariant._product);
  }

  virtual void backward_array_load(variable_t lhs, variable_t a,
                                   linear_expression_t elem_size,
                                   linear_expression_t i,
                                   nn_domain_t invariant) override {
    this->_product.backward_array_load(lhs, a, elem_size, i,
                                       invariant._product);
  }

  virtual void backward_array_store(variable_t a, linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    nn_domain_t invariant) override {
    this->_product.backward_array_store(a, elem_size, i, val, is_strong_update,
                                        invariant._product);
  }

  virtual void backward_array_store(variable_t a_new, variable_t a_old,
                                    linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    nn_domain_t invariant) override {
    this->_product.backward_array_store(a_new, a_old, elem_size, i, val,
                                        is_strong_update, invariant._product);
  }

  virtual void backward_array_store_range(variable_t a,
                                          linear_expression_t elem_size,
                                          linear_expression_t i,
                                          linear_expression_t j,
                                          linear_expression_t val,
                                          nn_domain_t invariant) override {
    this->_product.backward_array_store_range(a, elem_size, i, j, val,
                                              invariant._product);
  }

  virtual void backward_array_store_range(variable_t a_new, variable_t a_old,
                                          linear_expression_t elem_size,
                                          linear_expression_t i,
                                          linear_expression_t j,
                                          linear_expression_t val,
                                          nn_domain_t invariant) override {
    this->_product.backward_array_store_range(a_new, a_old, elem_size, i, j,
                                              val, invariant._product);
  }

  virtual void backward_array_assign(variable_t lhs, variable_t rhs,
                                     nn_domain_t invariant) override {
    this->_product.backward_array_assign(lhs, rhs, invariant._product);
  }

  // pointer operators
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    this->_product.pointer_load(lhs, rhs);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    this->_product.pointer_store(lhs, rhs);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    this->_product.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    this->_product.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    this->_product.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    this->_product.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(pointer_constraint_t cst) override {
    this->_product.pointer_assume(cst);
  }

  virtual void pointer_assert(pointer_constraint_t cst) override {
    this->_product.pointer_assert(cst);
  }

  void write(crab_os &o) { this->_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    return this->_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    return this->_product.to_disjunctive_linear_constraint_system();
  }

  static std::string getDomainName() {
    return domain_product2_t::getDomainName();
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    this->_product.rename(from, to);
  }

  /* begin intrinsics operations */    
  virtual void intrinsic(std::string name,
			 const variable_vector_t &inputs,
			 const variable_vector_t &outputs) override {
    this->_product.intrinsic(name, inputs, outputs);
  }

  virtual void backward_intrinsic(std::string name,
				  const variable_vector_t &inputs,
				  const variable_vector_t &outputs,
				  nn_domain_t invariant) override {
    this->_product.backward_intrinsic(name, inputs, outputs, invariant._product);
  }
  /* end intrinsics operations */
  
  void forget(const variable_vector_t &variables) {
    this->_product.forget(variables);
  }

  void project(const variable_vector_t &variables) {
    this->_product.project(variables);
  }

  void expand(variable_t var, variable_t new_var) {
    this->_product.expand(var, new_var);
  }

  void normalize() { this->_product.normalize(); }

  void minimize() { this->_product.minimize(); }

}; // class numerical_nullity_domain

template <typename Number, typename VariableName, typename Domain1,
          typename Domain2>
struct abstract_domain_traits<
    domain_product2<Number, VariableName, Domain1, Domain2>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

template <typename Domain1, typename Domain2, class Params>
struct abstract_domain_traits<
    reduced_numerical_domain_product2<Domain1, Domain2, Params>> {
  typedef typename Domain1::number_t number_t;
  typedef typename Domain1::varname_t varname_t;
};

template <typename NumAbsDom>
struct abstract_domain_traits<numerical_congruence_domain<NumAbsDom>> {
  typedef typename NumAbsDom::number_t number_t;
  typedef typename NumAbsDom::varname_t varname_t;
};

template <typename Domain>
struct abstract_domain_traits<numerical_nullity_domain<Domain>> {
  typedef typename Domain::number_t number_t;
  typedef typename Domain::varname_t varname_t;
};

template <typename Dom>
struct array_graph_domain_helper_traits<numerical_nullity_domain<Dom>> {
  typedef numerical_nullity_domain<Dom> num_null_domain_t;
  typedef typename num_null_domain_t::linear_constraint_t linear_constraint_t;
  typedef typename Dom::variable_t variable_t;

  static bool is_unsat(num_null_domain_t &inv, linear_constraint_t cst) {
    return array_graph_domain_helper_traits<Dom>::is_unsat(inv.first(), cst);
  }

  static void active_variables(num_null_domain_t &inv,
                               std::vector<variable_t> &out) {
    array_graph_domain_helper_traits<Dom>::active_variables(inv.first(), out);
  }
};

} // end namespace domains
} // namespace crab
