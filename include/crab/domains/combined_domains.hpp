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
 **************************************************************************/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

// Provided by Ikos
// Reduced product of two arbitrary domains with only lattice
// operations.
template <typename Domain1, typename Domain2>
class basic_domain_product2 {

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
      // we don't apply normalization when widening
      canonicalize();
    }
  }

  basic_domain_product2(const basic_domain_product2_t &other) = default;
  basic_domain_product2(basic_domain_product2_t &&other) = default;
  basic_domain_product2_t &operator=(const basic_domain_product2_t &other)  = default;
  basic_domain_product2_t &operator=(basic_domain_product2_t &&other) = default;

  static basic_domain_product2_t top() {
    return basic_domain_product2_t(Domain1::top(), Domain2::top());
  }

  static basic_domain_product2_t bottom() {
    return basic_domain_product2_t(Domain1::bottom(), Domain2::bottom());
  }

  bool is_bottom() const {
    //this->canonicalize();
    // return this->_is_bottom;
    if (_is_bottom) {
      return true;
    } else {
      return _first.is_bottom() || _second.is_bottom();
    }
  }

  bool is_top() const {
    return (this->_first.is_top() && this->_second.is_top());
  }

  Domain1 &first(bool apply_reduction = true) {
    if (apply_reduction) {
      // if called during widening we don't normalize
      this->canonicalize();
    }
    return this->_first;
  }

  Domain2 &second(bool apply_reduction = true) {
    if (apply_reduction) {
      // if called during widening we don't normalize      
      this->canonicalize();
    }
    return this->_second;
  }

  const Domain1 &first() const {
    return this->_first;
  }

  const Domain2 &second() const {
    return this->_second;
  }
  
  bool operator<=(const basic_domain_product2_t &other) const {
    if (this->is_bottom()) {
      return true;
    } else if (other.is_bottom()) {
      return false;
    } else {
      return (this->_first <= other._first) && (this->_second <= other._second);
    }
  }

  bool operator==(const basic_domain_product2_t &other) const {
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

  void write(crab::crab_os &o) const {
    if (this->is_bottom()) {
      o << "_|_";
    } else {
      o << "(" << this->_first << ", " << this->_second << ")";
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, basic_domain_product2_t &dom) {
    dom.write(o);
    return o;
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
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::interval_t;
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
  void set_to_top() override {
    domain_product2_t abs(basic_domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    domain_product2_t abs(basic_domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  domain_product2() : _product() {}

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

  bool is_bottom() const override { return this->_product.is_bottom(); }

  bool is_top() const override { return this->_product.is_top(); }

  Domain1 &first() { return this->_product.first(); }
  const Domain1 &first() const { return this->_product.first(); }

  Domain2 &second() { return this->_product.second(); }
  const Domain2 &second() const { return this->_product.second(); }

  bool operator<=(const domain_product2_t &other) const override {
    return (this->_product <= other._product);
  }

  bool operator==(domain_product2_t other) {
    return (this->_product == other._product);
  }

  void operator|=(domain_product2_t other) override {
    this->_product |= other._product;
  }

  domain_product2_t operator|(domain_product2_t other) override {
    return domain_product2_t(this->_product | other._product);
  }

  domain_product2_t operator&(domain_product2_t other) override {
    return domain_product2_t(this->_product & other._product);
  }

  domain_product2_t operator||(domain_product2_t other) override {
    return domain_product2_t(this->_product || other._product);
  }

  domain_product2_t
  widening_thresholds(domain_product2_t other,
                      const iterators::thresholds<number_t> &ts) override {
    bool apply_reduction = false;
    return domain_product2_t(basic_domain_product2_t(
        std::move(this->_product.first(apply_reduction)
                      .widening_thresholds(other._product.first(), ts)),
        std::move(this->_product.second(apply_reduction)
                      .widening_thresholds(other._product.second(), ts)),
        std::move(apply_reduction)));
  }

  domain_product2_t operator&&(domain_product2_t other) override {
    return domain_product2_t(this->_product && other._product);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    this->_product.first().assign(x, e);
    this->_product.second().assign(x, e);
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.first().apply(op, x, y, z);
    this->_product.second().apply(op, x, y, z);
    this->reduce();
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, Number k) override {
    this->_product.first().apply(op, x, y, k);
    this->_product.second().apply(op, x, y, k);
    this->reduce();
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       domain_product2_t invariant) override {
    this->_product.first().backward_assign(x, e, invariant.first());
    this->_product.second().backward_assign(x, e, invariant.second());
    this->reduce();
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, Number k,
                      domain_product2_t invariant) override {
    this->_product.first().backward_apply(op, x, y, k, invariant.first());
    this->_product.second().backward_apply(op, x, y, k, invariant.second());
    this->reduce();
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, const variable_t &z,
                      domain_product2_t invariant) override {
    this->_product.first().backward_apply(op, x, y, z, invariant.first());
    this->_product.second().backward_apply(op, x, y, z, invariant.second());
    this->reduce();
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    this->_product.first() += csts;
    this->_product.second() += csts;
    this->reduce();
  }

  void operator-=(const variable_t &v) override {
    this->_product.first() -= v;
    this->_product.second() -= v;
  }

  // cast operators

  void apply(int_conv_operation_t op,
	     const variable_t &dst, const variable_t &src) override {
    this->_product.first().apply(op, dst, src);
    this->_product.second().apply(op, dst, src);
    this->reduce();
  }

  // bitwise operators

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.first().apply(op, x, y, z);
    this->_product.second().apply(op, x, y, z);
    this->reduce();
  }

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, Number k) override {
    this->_product.first().apply(op, x, y, k);
    this->_product.second().apply(op, x, y, k);
    this->reduce();
  }

  // array operators

  virtual void array_init(const variable_t &a, const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    this->_product.first().array_init(a, elem_size, lb_idx, ub_idx, val);
    this->_product.second().array_init(a, elem_size, lb_idx, ub_idx, val);
    this->reduce();
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) override {

    this->_product.first().array_load(lhs, a, elem_size, i);
    this->_product.second().array_load(lhs, a, elem_size, i);
    this->reduce();
  }

  virtual void array_store(const variable_t &a, const linear_expression_t &elem_size,
                           const linear_expression_t &i, const linear_expression_t &val,
                           bool is_strong_update) override {
    this->_product.first().array_store(a, elem_size, i, val, is_strong_update);
    this->_product.second().array_store(a, elem_size, i, val, is_strong_update);
    this->reduce();
  }

  virtual void array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                                 const linear_expression_t &i, const linear_expression_t &j,
                                 const linear_expression_t &val) override {
    this->_product.first().array_store_range(a, elem_size, i, j, val);
    this->_product.second().array_store_range(a, elem_size, i, j, val);
    this->reduce();
  }

  virtual void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    this->_product.first().array_assign(lhs, rhs);
    this->_product.second().array_assign(lhs, rhs);
    this->reduce();
  }

  // backward array operations

  virtual void backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                                   const linear_expression_t &lb_idx,
                                   const linear_expression_t &ub_idx,
                                   const linear_expression_t &val,
                                   domain_product2_t invariant) override {
    this->_product.first().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                               val, invariant.first());
    this->_product.second().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                                val, invariant.second());
    this->reduce();
  }

  virtual void backward_array_load(const variable_t &lhs, const variable_t &a,
                                   const linear_expression_t &elem_size,
                                   const linear_expression_t &i,
                                   domain_product2_t invariant) override {

    this->_product.first().backward_array_load(lhs, a, elem_size, i,
                                               invariant.first());
    this->_product.second().backward_array_load(lhs, a, elem_size, i,
                                                invariant.second());
    this->reduce();
  }

  virtual void backward_array_store(const variable_t &a, const linear_expression_t &elem_size,
                                    const linear_expression_t &i,
                                    const linear_expression_t &val,
                                    bool is_strong_update,
                                    domain_product2_t invariant) override {
    this->_product.first().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.first());
    this->_product.second().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.second());
    this->reduce();
  }

  virtual void
  backward_array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                             const linear_expression_t &i, const linear_expression_t &j,
                             const linear_expression_t &val,
                             domain_product2_t invariant) override {
    this->_product.first().backward_array_store_range(a, elem_size, i, j, val,
                                                      invariant.first());
    this->_product.second().backward_array_store_range(a, elem_size, i, j, val,
                                                       invariant.second());
    this->reduce();
  }

  virtual void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                                     domain_product2_t invariant) override {
    this->_product.first().backward_array_assign(lhs, rhs, invariant.first());
    this->_product.second().backward_array_assign(lhs, rhs, invariant.second());
    this->reduce();
  }

  // reference  operators
  virtual void region_init(const memory_region &reg) override {
    this->_product.first().region_init(reg);
    this->_product.second().region_init(reg);
    this->reduce();
  }
  
  virtual void ref_make(const variable_t &ref, const memory_region &reg) override {
    this->_product.first().ref_make(ref, reg);
    this->_product.second().ref_make(ref, reg);
    this->reduce();
  }

  virtual void ref_load(const variable_t &ref, const memory_region &reg,
			const variable_t &res) override {
    this->_product.first().ref_load(ref, reg, res);
    this->_product.second().ref_load(ref, reg, res);
    this->reduce();
  }

  virtual void ref_store(const variable_t &ref, const memory_region &reg,
			 const linear_expression_t &val) override {
    this->_product.first().ref_store(ref, reg, val);
    this->_product.second().ref_store(ref, reg, val);
    this->reduce();
  }

  virtual void ref_gep(const variable_t &ref1, const memory_region &reg1,
		       const variable_t &ref2, const memory_region &reg2,
		       const linear_expression_t &offset) override {
    this->_product.first().ref_gep(ref1, reg1, ref2, reg2, offset);
    this->_product.second().ref_gep(ref1, reg1, ref2, reg2, offset);
    this->reduce();
  }

  virtual void ref_load_from_array(const variable_t &lhs,
				   const variable_t &ref, const memory_region &region,
				   const linear_expression_t &index,
				   const linear_expression_t &elem_size) override {
    this->_product.first().ref_load_from_array(lhs, ref, region, index, elem_size);
    this->_product.second().ref_load_from_array(lhs, ref, region, index, elem_size);
    this->reduce();
  }

  virtual void ref_store_to_array(const variable_t &ref, const memory_region &region,
				  const linear_expression_t &index,
				  const linear_expression_t &elem_size,
				  const linear_expression_t &val) override {
    this->_product.first().ref_store_to_array(ref, region, index, elem_size, val);
    this->_product.second().ref_store_to_array(ref, region, index, elem_size, val);
    this->reduce();
  }

  virtual void ref_assume(const reference_constraint_t &cst) override {
    this->_product.first().ref_assume(cst);
    this->_product.second().ref_assume(cst);
    this->reduce();
  }


  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    this->_product.first().assign_bool_cst(lhs, rhs);
    this->_product.second().assign_bool_cst(lhs, rhs);
    this->reduce();
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    this->_product.first().assign_bool_var(lhs, rhs, is_not_rhs);
    this->_product.second().assign_bool_var(lhs, rhs, is_not_rhs);
    this->reduce();
  }

  virtual void apply_binary_bool(bool_operation_t op,
				 const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.first().apply_binary_bool(op, x, y, z);
    this->_product.second().apply_binary_bool(op, x, y, z);
    this->reduce();
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    this->_product.first().assume_bool(v, is_negated);
    this->_product.second().assume_bool(v, is_negated);
    this->reduce();
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs,
                                        domain_product2_t inv) override {
    this->_product.first().backward_assign_bool_cst(lhs, rhs, inv.first());
    this->_product.second().backward_assign_bool_cst(lhs, rhs, inv.second());
    this->reduce();
  }

  virtual void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                        bool is_not_rhs,
                                        domain_product2_t inv) override {
    this->_product.first().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                    inv.first());
    this->_product.second().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                     inv.second());
    this->reduce();
  }

  virtual void backward_apply_binary_bool(bool_operation_t op,
					  const variable_t &x,
                                          const variable_t &y, const variable_t &z,
                                          domain_product2_t inv) override {
    this->_product.first().backward_apply_binary_bool(op, x, y, z, inv.first());
    this->_product.second().backward_apply_binary_bool(op, x, y, z,
                                                       inv.second());
    this->reduce();
  }

  virtual void forget(const variable_vector_t &variables) override {
    this->_product.first().forget(variables);
    this->_product.second().forget(variables);
  }
 
  virtual void project(const variable_vector_t &variables) override {
    this->_product.first().project(variables);
    this->_product.second().project(variables);
  }

  virtual void expand(const variable_t &var, const variable_t &new_var) override {
    this->_product.first().expand(var, new_var);
    this->_product.second().expand(var, new_var);
  }

  virtual void normalize() override {
    this->_product.first().normalize();
    this->_product.second().normalize();
  }

  virtual void minimize() override {
    this->_product.first().minimize();
    this->_product.second().minimize();
  }

  virtual interval_t operator[](const variable_t &v) override {
    return this->_product.first()[v] & this->_product.second()[v];
  }
  
  virtual linear_constraint_system_t to_linear_constraint_system() const override {
    linear_constraint_system_t csts;
    // XXX: We might add redundant constraints.
    csts += this->_product.first().to_linear_constraint_system();
    csts += this->_product.second().to_linear_constraint_system();
    return csts;
  }

  virtual disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
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
  
  void write(crab::crab_os &o) const override { this->_product.write(o); }

  std::string domain_name() const override {
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
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::interval_t;
  typedef typename Domain1::number_t number_t;
  typedef typename Domain1::varname_t varname_t;

  static_assert(std::is_same<number_t, typename Domain2::number_t>::value,
		"Domain1 and Domain2 must have same type for number_t");
  static_assert(std::is_same<varname_t, typename Domain2::varname_t>::value,
		"Domain1 and Domain2 must have same type for varname_t");
  
private:
  typedef domain_product2<number_t, varname_t, Domain1, Domain2>
      domain_product2_t;

  domain_product2_t _product;

  reduced_numerical_domain_product2(const domain_product2_t &product)
      : _product(product) {}

  linear_constraint_system_t to_linear_constraints(const variable_t &v,
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
    crab::CrabStats::count(domain_name() + ".count.reduce");
    crab::ScopedCrabStats __st__(domain_name() + ".reduce");

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
          std::string k(domain_name() + ".count.reduce.equalities_from_" +
                        _product.first().domain_name());
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
          std::string k(domain_name() + ".count.reduce.equalities_from_" +
                        _product.second().domain_name());
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
        std::string k(domain_name() + ".count.reduce.equalities_from_" +
                      _product.first().domain_name());
        crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts1.size());
        inv2 += csts1;
      } else if (Params::right_propagate_equalities ||
                 Params::right_propagate_inequalities) {
        linear_constraint_system_t csts2;
        const bool propagate_only_equalities =
            !Params::right_propagate_inequalities;
        crab::domains::reduced_domain_traits<Domain2>::extract(
            inv2, v, csts2, propagate_only_equalities);
        std::string k(domain_name() + ".count.reduce.equalities_from_" +
                      _product.second().domain_name());
        crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts2.size());
        inv1 += csts2;
      }
    }
  }

public:
  void set_to_top() override {
    reduced_numerical_domain_product2_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
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

  bool is_bottom() const override { return this->_product.is_bottom(); }

  bool is_top() const override { return this->_product.is_top(); }

  //Domain1 &first() { return this->_product.first(); }
  //Domain2 &second() { return this->_product.second(); }

  bool operator<=(const reduced_numerical_domain_product2_t &other) const override {
    return this->_product <= other._product;
  }

  void operator|=(reduced_numerical_domain_product2_t other) override {
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";);
    this->_product |= other._product;
    CRAB_LOG("reduced-dom", crab::outs() << *this << "\n----------------\n";);
  }

  reduced_numerical_domain_product2_t
  operator|(reduced_numerical_domain_product2_t other) override {
    reduced_numerical_domain_product2_t res(this->_product | other._product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";
             crab::outs() << res << "\n================\n");
    return res;
  }

  reduced_numerical_domain_product2_t
  operator&(reduced_numerical_domain_product2_t other) override {
    reduced_numerical_domain_product2_t res(this->_product & other._product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ MEET ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator||(reduced_numerical_domain_product2_t other) override {
    reduced_numerical_domain_product2_t res(this->_product || other._product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  widening_thresholds(reduced_numerical_domain_product2_t other,
                      const iterators::thresholds<number_t> &ts) override {
    reduced_numerical_domain_product2_t res(
        this->_product.widening_thresholds(other._product, ts));
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator&&(reduced_numerical_domain_product2_t other) override {
    reduced_numerical_domain_product2_t res(this->_product && other._product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ NARROWING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  void set(const variable_t &v, interval_t x) {
    this->_product.first().set(v, x);
    this->_product.second().set(v, x);
  }

  interval_t operator[](const variable_t &v) override {
    return this->_product.first()[v] & this->_product.second()[v];
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    this->_product += csts;
    for (auto const& cst: csts) {
      for (auto const&v: cst.variables()) {
	reduce_variable(v);
      }
    }
    CRAB_LOG("combined-domain", crab::outs() << "Added constraints " << csts
                                             << "=" << *this << "\n");
  }

  void operator-=(const variable_t &v) override { this->_product -= v; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    this->_product.assign(x, e);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain", crab::outs()
                                    << x << ":=" << e << "=" << *this << "\n");
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << z << "=" << *this << "\n");
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    this->_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << k << "=" << *this << "\n");
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       reduced_numerical_domain_product2_t invariant) override {
    this->_product.backward_assign(x, e, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      for (auto const&v : e.variables())
        this->reduce_variable(v);
    }
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, number_t k,
                      reduced_numerical_domain_product2_t invariant) override {
    this->_product.backward_apply(op, x, y, k, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      this->reduce_variable(y);
    }
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, const variable_t &z,
                      reduced_numerical_domain_product2_t invariant) override {
    this->_product.backward_apply(op, x, y, z, invariant._product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      this->reduce_variable(y);
      this->reduce_variable(z);
    }
  }

  // cast operators

  void apply(int_conv_operation_t op, const variable_t &dst, const variable_t &src) override {
    this->_product.apply(op, dst, src);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(dst);
    }
  }

  // bitwise operators

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
  }

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    this->_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      this->reduce_variable(x);
    }
  }

  /*
     Begin unimplemented operations

     reduced_numerical_domain_product2 implements only standard
     abstract operations of a numerical domain.  
  */

  // boolean operations
  void assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs) override {}
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs) override {}
  void apply_binary_bool(bool_operation_t op,
			 const variable_t &x, const variable_t &y, const variable_t &z) override {}
  void assume_bool(const variable_t &v, bool is_negated) override {}
  // backward boolean operations
  void backward_assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs,
                                reduced_numerical_domain_product2_t invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs,
                                reduced_numerical_domain_product2_t invariant) override {}
  void
  backward_apply_binary_bool(bool_operation_t op,
			     const variable_t &x, const variable_t &y, const variable_t &z,
                             reduced_numerical_domain_product2_t invariant) override {}
  // array operations
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  void array_load(const variable_t &lhs, const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {}
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {}
  void array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                         const linear_expression_t &i, const linear_expression_t &j,
                         const linear_expression_t &val) override {}
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {}
  // backward array operations
  void backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx, const linear_expression_t &val,
                           reduced_numerical_domain_product2_t invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size, const linear_expression_t &i,
                           reduced_numerical_domain_product2_t invariant) override {}
  void backward_array_store(const variable_t &a, const linear_expression_t &elem_size,
                            const linear_expression_t &i, const linear_expression_t &v,
                            bool is_strong_update,
                            reduced_numerical_domain_product2_t invariant) override {}
  void
  backward_array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                             const linear_expression_t &i, const linear_expression_t &j,
                             const linear_expression_t &val,
                             reduced_numerical_domain_product2_t invariant) override {}
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             reduced_numerical_domain_product2_t invariant) override {}
  // reference operations
  void region_init(const memory_region &reg) override {}         
  void ref_make(const variable_t &ref, const memory_region &reg) override {}
  void ref_load(const variable_t &ref, const memory_region &reg, const variable_t &res) override {}
  void ref_store(const variable_t &ref, const memory_region &reg,
		 const linear_expression_t &val) override {}
  void ref_gep(const variable_t &ref1, const memory_region &reg1,
	       const variable_t &ref2, const memory_region &reg2,
	       const linear_expression_t &offset) override {}
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref, const memory_region &region,
			   const linear_expression_t &index,
			   const linear_expression_t &elem_size) override {}
  void ref_store_to_array(const variable_t &ref, const memory_region &region,
			  const linear_expression_t &index, const linear_expression_t &elem_size,
			  const linear_expression_t &val) override {}
  void ref_assume(const reference_constraint_t &cst) override {}
  /* End unimplemented operations */

  void rename(const variable_vector_t &from, const variable_vector_t &to) override {
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
  
  void forget(const variable_vector_t &variables) override {
    this->_product.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    this->_product.project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    this->_product.expand(var, new_var);
  }

  void normalize() override { this->_product.normalize(); }

  void minimize() override { this->_product.minimize(); }

  void write(crab_os &o) const override { this->_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() const override {
    return this->_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return this->_product.to_disjunctive_linear_constraint_system();
  }

  std::string domain_name() const override {
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
  typedef ikos::interval<Number> interval_t;
  typedef ikos::congruence<Number> congruence_t;
  typedef ikos::bound<Number> bound_t;

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
  const interval_t &first() const { return this->_first; }  

  congruence_t &second() { return this->_second; }
  const congruence_t &second() const { return this->_second; }  

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
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::interval_t;
  typedef typename NumAbsDom::number_t number_t;
  typedef typename NumAbsDom::varname_t varname_t;

  typedef ikos::congruence_domain<number_t, varname_t> congruence_domain_t;
  typedef interval_congruence<number_t> interval_congruence_t;

private:
  typedef domain_product2<number_t, varname_t, NumAbsDom, congruence_domain_t>
      domain_product2_t;

  domain_product2_t _product;

  numerical_congruence_domain(const domain_product2_t &product)
      : _product(product) {}

  void reduce_variable(const variable_t &v) {
    crab::CrabStats::count(domain_name() + ".count.reduce");
    crab::ScopedCrabStats __st__(domain_name() + ".reduce");

    if (is_bottom())
      return;

    auto i = this->_product.first()[v]; // project on intervals
    auto c = this->_product.second().to_congruence(v);
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

public:
  void set_to_top() override {
    rnc_domain_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
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

  bool is_bottom() const override { return this->_product.is_bottom(); }

  bool is_top() const override { return this->_product.is_top(); }

  //NumAbsDom &first() { return this->_product.first(); }
  //congruence_domain_t &second() { return this->_product.second(); }

  bool operator<=(const rnc_domain_t &other) const override {
    return this->_product <= other._product;
  }

  bool operator==(rnc_domain_t other) {
    return this->_product == other._product;
  }

  void operator|=(rnc_domain_t other) override { this->_product |= other._product; }

  rnc_domain_t operator|(rnc_domain_t other) override {
    return rnc_domain_t(this->_product | other._product);
  }

  rnc_domain_t operator&(rnc_domain_t other) override {
    return rnc_domain_t(this->_product & other._product);
  }

  rnc_domain_t operator||(rnc_domain_t other) override {
    return rnc_domain_t(this->_product || other._product);
  }

  rnc_domain_t widening_thresholds(rnc_domain_t other,
                                   const iterators::thresholds<number_t> &ts) override {
    return rnc_domain_t(this->_product.widening_thresholds(other._product, ts));
  }

  rnc_domain_t operator&&(rnc_domain_t other) override {
    return rnc_domain_t(this->_product && other._product);
  }

  // pre: x is already reduced
  void set(const variable_t &v, interval_congruence_t x) {
    this->_product.first().set(v, x.first());
    this->_product.second().set(v, x.second());
  }

  interval_congruence_t get(const variable_t &v) {
    return interval_congruence_t(this->_product.first()[v],
                                 this->_product.second().to_congruence(v));
  }

  interval_t operator[](const variable_t &v) override {
    interval_congruence_t x = get(v);
    return x.first();
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    this->_product += csts;

    if (!is_bottom()) {
      for (auto const& cst: csts) {
	for (auto const& v: cst.variables()) {
	  this->reduce_variable(v);
	  if (is_bottom()) {
	    return;
	  }
	}
      }
    }
  }

  void operator-=(const variable_t &v) override { this->_product -= v; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    this->_product.assign(x, e);
    this->reduce_variable(x);
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.apply(op, x, y, z);
    this->reduce_variable(x);
  }

  void apply(arith_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    this->_product.apply(op, x, y, k);
    this->reduce_variable(x);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       rnc_domain_t invariant) override {
    this->_product.backward_assign(x, e, invariant._product);
    // reduce the variables in the right-hand side
    for (auto const&v : e.variables())
      this->reduce_variable(v);
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, number_t k,
                      rnc_domain_t invariant) override {
    this->_product.backward_apply(op, x, y, k, invariant._product);
    // reduce the variables in the right-hand side
    this->reduce_variable(y);
  }

  void backward_apply(arith_operation_t op,
		      const variable_t &x, const variable_t &y, const variable_t &z,
                      rnc_domain_t invariant) override {
    this->_product.backward_apply(op, x, y, z, invariant._product);
    // reduce the variables in the right-hand side
    this->reduce_variable(y);
    this->reduce_variable(z);
  }

  // cast operators

  void apply(int_conv_operation_t op,
	     const variable_t &dst, const variable_t &src) override {
    this->_product.apply(op, dst, src);
    this->reduce_variable(dst);
  }

  // bitwise operators

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    this->_product.apply(op, x, y, z);
    this->reduce_variable(x);
  }

  void apply(bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    this->_product.apply(op, x, y, k);
    this->reduce_variable(x);
  }

  /*
     Begin unimplemented operations

     numerical_congruence_domain implements only standard abstract
     operations of a numerical domain.
  */

  // boolean operations
  void assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs) override {}
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs) override {}
  void apply_binary_bool(bool_operation_t op,
			 const variable_t &x, const variable_t &y, const variable_t &z) override {}
  void assume_bool(const variable_t &v, bool is_negated) override {}
  // backward boolean operations
  void backward_assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs,
                                rnc_domain_t invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs,
                                rnc_domain_t invariant) override {}
  void backward_apply_binary_bool(bool_operation_t op,
                                  const variable_t &x, const variable_t &y, const variable_t &z,
                                  rnc_domain_t invariant) override {}
  // array operations
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  void array_load(const variable_t &lhs, const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {}
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {}
  void array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                         const linear_expression_t &i, const linear_expression_t &j,
                         const linear_expression_t &val) override {}
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {}
  // backward array operations
  void backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
			   const linear_expression_t &val,
                           rnc_domain_t invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
			   const linear_expression_t &i,
                           rnc_domain_t invariant) override {}
  void backward_array_store(const variable_t &a, const linear_expression_t &elem_size,
                            const linear_expression_t &i, const linear_expression_t &v,
                            bool is_strong_update, rnc_domain_t invariant) override {}
  void backward_array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                                  const linear_expression_t &i, const linear_expression_t &j,
                                  const linear_expression_t &val,
                                  rnc_domain_t invariant) override {}
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             rnc_domain_t invariant) override {}
  // reference operations
  void region_init(const memory_region &reg) override {}         
  void ref_make(const variable_t &ref, const memory_region &reg) override {}
  void ref_load(const variable_t &ref, const memory_region &reg, const variable_t &res) override {}
  void ref_store(const variable_t &ref, const memory_region &reg,
		 const linear_expression_t &val) override {}
  void ref_gep(const variable_t &ref1, const memory_region &reg1,
	       const variable_t &ref2, const memory_region &reg2,
	       const linear_expression_t &offset) override {}
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref, const memory_region &region,
			   const linear_expression_t &index,
			   const linear_expression_t &elem_size) override {}
  void ref_store_to_array(const variable_t &ref, const memory_region &region,
			  const linear_expression_t &index, const linear_expression_t &elem_size,
			  const linear_expression_t &val) override {}
  void ref_assume(const reference_constraint_t &cst) override {}
  /* End unimplemented operations */

  void forget(const variable_vector_t &variables) override {
    this->_product.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    this->_product.project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    this->_product.expand(var, new_var);
  }

  void normalize() override { this->_product.normalize(); }

  void minimize() override { this->_product.minimize(); }

  void write(crab_os &o) const override { this->_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() const override {
    return this->_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return this->_product.to_disjunctive_linear_constraint_system();
  }

  std::string domain_name() const override {
    return domain_product2_t::getDomainName();
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) override {
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


} // end namespace domains
} // namespace crab
