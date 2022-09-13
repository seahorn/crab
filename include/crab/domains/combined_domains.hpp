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
 **************************************************************************/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/lattice_domain.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>
#include <algorithm>

namespace crab {
namespace domains {
  
// Reduced product of two arbitrary domains with only lattice
// operations.
template <typename Domain1, typename Domain2>
class basic_domain_product2:
    public lattice_domain_api<basic_domain_product2<Domain1, Domain2>> {

public:
  using basic_domain_product2_t = basic_domain_product2<Domain1, Domain2>;
  using first_type = Domain1;
  using second_type = Domain2;

private:
  bool m_is_bottom;
  Domain1 m_first;
  Domain2 m_second;

  void canonicalize() {
    if (!m_is_bottom) {
      m_is_bottom = m_first.is_bottom() || m_second.is_bottom();
      if (m_is_bottom) {
        m_first.set_to_bottom();
        m_second.set_to_bottom();
      }
    }
  }

public:
  basic_domain_product2() : m_is_bottom(false) {
    m_first.set_to_top();
    m_second.set_to_top();
  }

  basic_domain_product2(Domain1 &&first, Domain2 &&second,
                        bool &&apply_reduction = true)
      : m_is_bottom(false), m_first(std::move(first)),
        m_second(std::move(second)) {
    if (apply_reduction) {
      // we don't apply normalization when widening
      canonicalize();
    }
  }

  basic_domain_product2(const basic_domain_product2_t &other) = default;
  basic_domain_product2(basic_domain_product2_t &&other) = default;
  basic_domain_product2_t &
  operator=(const basic_domain_product2_t &other) = default;
  basic_domain_product2_t &operator=(basic_domain_product2_t &&other) = default;

  basic_domain_product2_t make_top() const override {
    Domain1 dom1;
    Domain2 dom2;
    dom1.set_to_top();
    dom2.set_to_top();
    return basic_domain_product2_t(std::move(dom1), std::move(dom2));
  }

  basic_domain_product2_t make_bottom() const override {
    Domain1 dom1;
    Domain2 dom2;
    dom1.set_to_bottom();
    dom2.set_to_bottom();
    return basic_domain_product2_t(std::move(dom1), std::move(dom2));
  }

  void set_to_top() override {
    m_is_bottom = false;
    m_first.set_to_top();
    m_second.set_to_top();
  }

  void set_to_bottom() override {
    m_is_bottom = true;
    m_first.set_to_bottom();
    m_second.set_to_bottom();
  }

  bool is_bottom() const override {
    // canonicalize();
    // return m_is_bottom;
    if (m_is_bottom) {
      return true;
    } else {
      return m_first.is_bottom() || m_second.is_bottom();
    }
  }

  bool is_top() const override {
    return (m_first.is_top() && m_second.is_top());
  }

  Domain1 &first() {
    canonicalize();
    return m_first;
  }

  Domain2 &second() {
    canonicalize();
    return m_second;
  }

  const Domain1 &first() const { return m_first; }

  const Domain2 &second() const { return m_second; }

  bool operator<=(const basic_domain_product2_t &other) const override {
    if (is_bottom()) {
      return true;
    } else if (other.is_bottom()) {
      return false;
    } else {
      return (m_first <= other.m_first) && (m_second <= other.m_second);
    }
  }

  bool operator==(const basic_domain_product2_t &other) const {
    return (operator<=(other) && other.operator<=(*this));
  }

  void operator|=(const basic_domain_product2_t &other) override {
    if (is_bottom()) {
      *this = other;
    } else if (other.is_bottom()) {
      return;
    } else {
      m_first |= other.m_first;
      m_second |= other.m_second;
    }
  }

  basic_domain_product2_t
  operator|(const basic_domain_product2_t &other) const override {
    if (is_bottom()) {
      return other;
    } else if (other.is_bottom()) {
      return *this;
    } else {
      return basic_domain_product2_t(m_first | other.m_first,
                                     m_second | other.m_second);
    }
  }

  basic_domain_product2_t
  operator||(const basic_domain_product2_t &other) const override {
    return basic_domain_product2_t(m_first || other.m_first,
                                   m_second || other.m_second,
                                   false /* do not apply reduction */);
  }

  void operator&=(const basic_domain_product2_t &other) {
    if (is_bottom()) {
      // do nothing
    } else if (other.is_bottom()) {
      *this = other;
    } else {
      m_first &= other.m_first;
      m_second &= other.m_second;
    }
  }
  
  basic_domain_product2_t
  operator&(const basic_domain_product2_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      return basic_domain_product2_t(m_first & other.m_first,
                                     m_second & other.m_second);
    }
  }

  basic_domain_product2_t
  operator&&(const basic_domain_product2_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      return basic_domain_product2_t(m_first && other.m_first,
                                     m_second && other.m_second);
    }
  }

  void write(crab::crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else {
      o << "(";
      m_first.write(o);
      o << ", ";
      m_second.write(o);
      o << ")";
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   basic_domain_product2_t &dom) {
    dom.write(o);
    return o;
  }

  std::string domain_name() const override {
    std::string name =
        "Product(" + m_first.domain_name() + "," + m_second.domain_name() + ")";
    return name;
  }

}; // class basic_domain_product2

// Reduced product of two arbitrary domains with all operations.
template <typename Number, typename VariableName, typename Domain1,
          typename Domain2>
class reduced_domain_product2 final
    : public abstract_domain_api<
          reduced_domain_product2<Number, VariableName, Domain1, Domain2>> {
public:
  using reduced_domain_product2_t =
      reduced_domain_product2<Number, VariableName, Domain1, Domain2>;
  using abstract_domain_t = abstract_domain_api<reduced_domain_product2_t>;
  using first_type = Domain1;
  using second_type = Domain2;

  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using number_t = Number;
  using varname_t = VariableName;

private:
  using basic_domain_product2_t = basic_domain_product2<Domain1, Domain2>;

  basic_domain_product2_t m_product;

  reduced_domain_product2(basic_domain_product2_t &&product)
      : m_product(std::move(product)) {}

  // Reduce operation
  void reduce() {
    if (m_product.first().is_bottom() ||
        m_product.second().is_bottom()) {
      m_product.set_to_bottom();
    }
  }

public:
  reduced_domain_product2_t make_top() const override {
    basic_domain_product2_t dom_prod;
    return reduced_domain_product2_t(dom_prod.make_top());
  }

  reduced_domain_product2_t make_bottom() const override {
    basic_domain_product2_t dom_prod;
    return reduced_domain_product2_t(dom_prod.make_bottom());
  }

  void set_to_top() override {
    basic_domain_product2_t dom_prod;
    reduced_domain_product2_t dom(dom_prod.make_top());
    std::swap(*this, dom);
  }

  void set_to_bottom() override {
    basic_domain_product2_t dom_prod;
    reduced_domain_product2_t dom(dom_prod.make_bottom());
    std::swap(*this, dom);
  }

  reduced_domain_product2() : m_product() {}

  reduced_domain_product2(Domain1 &&val1, Domain2 &&val2)
    : m_product(std::move(val1), std::move(val2)) {}
  
  reduced_domain_product2(const reduced_domain_product2_t &other)
    : m_product(other.m_product) {}

  reduced_domain_product2(const reduced_domain_product2_t &&other)
      : m_product(std::move(other.m_product)) {}

  reduced_domain_product2_t &operator=(const reduced_domain_product2_t &other) {
    if (this != &other)
      m_product = other.m_product;
    return *this;
  }

  reduced_domain_product2_t &operator=(const reduced_domain_product2_t &&other) {
    if (this != &other)
      m_product = std::move(other.m_product);
    return *this;
  }

  bool is_bottom() const override { return m_product.is_bottom(); }

  bool is_top() const override { return m_product.is_top(); }

  Domain1 &first() { return m_product.first(); }
  const Domain1 &first() const { return m_product.first(); }

  Domain2 &second() { return m_product.second(); }
  const Domain2 &second() const { return m_product.second(); }

  bool operator<=(const reduced_domain_product2_t &other) const override {
    return (m_product <= other.m_product);
  }

  bool operator==(const reduced_domain_product2_t &other) const {
    return (m_product == other.m_product);
  }

  void operator|=(const reduced_domain_product2_t &other) override {
    m_product |= other.m_product;
  }

  reduced_domain_product2_t operator|(const reduced_domain_product2_t &other) const override {
    return reduced_domain_product2_t(m_product | other.m_product);
  }

  void operator&=(const reduced_domain_product2_t &other) override {
    m_product &= other.m_product;
  }
  
  reduced_domain_product2_t operator&(const reduced_domain_product2_t &other) const override {
    return reduced_domain_product2_t(m_product & other.m_product);
  }

  reduced_domain_product2_t operator||(const reduced_domain_product2_t &other) const override {
    return reduced_domain_product2_t(m_product || other.m_product);
  }

  reduced_domain_product2_t widening_thresholds(
      const reduced_domain_product2_t &other,
      const thresholds<number_t> &ts) const override {
    bool apply_reduction = false;
    return reduced_domain_product2_t(basic_domain_product2_t(
        std::move(m_product.first().widening_thresholds(
            other.m_product.first(), ts)),
        std::move(m_product.second().widening_thresholds(
            other.m_product.second(), ts)),
        std::move(apply_reduction)));
  }

  reduced_domain_product2_t operator&&(const reduced_domain_product2_t &other) const override {
    return reduced_domain_product2_t(m_product && other.m_product);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.first().assign(x, e);
    m_product.second().assign(x, e);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.first().weak_assign(x, e);
    m_product.second().weak_assign(x, e);
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.first().apply(op, x, y, z);
    m_product.second().apply(op, x, y, z);
    reduce();
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             Number k) override {
    m_product.first().apply(op, x, y, k);
    m_product.second().apply(op, x, y, k);
    reduce();
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    m_product.first().select(lhs, cond, e1, e2);
    m_product.second().select(lhs, cond, e1, e2);
    reduce();
  }  

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_assign(x, e, invariant.first());
    m_product.second().backward_assign(x, e, invariant.second());
    reduce();
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, Number k,
                      const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_apply(op, x, y, k, invariant.first());
    m_product.second().backward_apply(op, x, y, k, invariant.second());
    reduce();
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_apply(op, x, y, z, invariant.first());
    m_product.second().backward_apply(op, x, y, z, invariant.second());
    reduce();
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    m_product.first() += csts;
    m_product.second() += csts;
    reduce();
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (m_product.first().entails(cst)) {
      return true;
    }
    if (m_product.second().entails(cst)) {
      return true;
    }
    return false;
  }
  
  void operator-=(const variable_t &v) override {
    m_product.first() -= v;
    m_product.second() -= v;
  }

  // cast operators

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    m_product.first().apply(op, dst, src);
    m_product.second().apply(op, dst, src);
    reduce();
  }

  // bitwise operators

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.first().apply(op, x, y, z);
    m_product.second().apply(op, x, y, z);
    reduce();
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             Number k) override {
    m_product.first().apply(op, x, y, k);
    m_product.second().apply(op, x, y, k);
    reduce();
  }

  // array operators

  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    m_product.first().array_init(a, elem_size, lb_idx, ub_idx, val);
    m_product.second().array_init(a, elem_size, lb_idx, ub_idx, val);
    reduce();
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) override {

    m_product.first().array_load(lhs, a, elem_size, i);
    m_product.second().array_load(lhs, a, elem_size, i);
    reduce();
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    m_product.first().array_store(a, elem_size, i, val, is_strong_update);
    m_product.second().array_store(a, elem_size, i, val, is_strong_update);
    reduce();
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &val) override {
    m_product.first().array_store_range(a, elem_size, i, j, val);
    m_product.second().array_store_range(a, elem_size, i, j, val);
    reduce();
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    m_product.first().array_assign(lhs, rhs);
    m_product.second().array_assign(lhs, rhs);
    reduce();
  }

  // backward array operations

  virtual void
  backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                      const linear_expression_t &lb_idx,
                      const linear_expression_t &ub_idx,
                      const linear_expression_t &val,
                      const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                               val, invariant.first());
    m_product.second().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                                val, invariant.second());
    reduce();
  }

  virtual void
  backward_array_load(const variable_t &lhs, const variable_t &a,
                      const linear_expression_t &elem_size,
                      const linear_expression_t &i,
                      const reduced_domain_product2_t &invariant) override {

    m_product.first().backward_array_load(lhs, a, elem_size, i,
                                               invariant.first());
    m_product.second().backward_array_load(lhs, a, elem_size, i,
                                                invariant.second());
    reduce();
  }

  virtual void backward_array_store(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &val,
      bool is_strong_update, const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.first());
    m_product.second().backward_array_store(
        a, elem_size, i, val, is_strong_update, invariant.second());
    reduce();
  }

  virtual void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &val,
      const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_array_store_range(a, elem_size, i, j, val,
                                                      invariant.first());
    m_product.second().backward_array_store_range(a, elem_size, i, j, val,
                                                       invariant.second());
    reduce();
  }

  virtual void
  backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                        const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_array_assign(lhs, rhs, invariant.first());
    m_product.second().backward_array_assign(lhs, rhs, invariant.second());
    reduce();
  }

  // region/reference operators
  virtual void region_init(const variable_t &reg) override {
    m_product.first().region_init(reg);
    m_product.second().region_init(reg);
    reduce();
  }

  virtual void region_copy(const variable_t &lhs_reg,
                           const variable_t &rhs_reg) override {
    m_product.first().region_copy(lhs_reg, rhs_reg);
    m_product.second().region_copy(lhs_reg, rhs_reg);
    reduce();
  }

  virtual void region_cast(const variable_t &src_reg,
                           const variable_t &dst_reg) override {
    m_product.first().region_cast(src_reg, dst_reg);
    m_product.second().region_cast(src_reg, dst_reg);
    reduce();
  }
  
  virtual void ref_make(const variable_t &ref, const variable_t &reg,
			const variable_or_constant_t &size,
			const allocation_site &as) override {
    m_product.first().ref_make(ref, reg, size, as);
    m_product.second().ref_make(ref, reg, size, as);
    reduce();
  }

  virtual void ref_free(const variable_t &reg, const variable_t &ref) override {
    m_product.first().ref_free(reg, ref);
    m_product.second().ref_free(reg, ref);
    reduce();
  }
  
  virtual void ref_load(const variable_t &ref, const variable_t &reg,
                        const variable_t &res) override {
    m_product.first().ref_load(ref, reg, res);
    m_product.second().ref_load(ref, reg, res);
    reduce();
  }

  virtual void ref_store(const variable_t &ref, const variable_t &reg,
                         const variable_or_constant_t &val) override {
    m_product.first().ref_store(ref, reg, val);
    m_product.second().ref_store(ref, reg, val);
    reduce();
  }

  virtual void ref_gep(const variable_t &ref1, const variable_t &reg1,
                       const variable_t &ref2, const variable_t &reg2,
                       const linear_expression_t &offset) override {
    m_product.first().ref_gep(ref1, reg1, ref2, reg2, offset);
    m_product.second().ref_gep(ref1, reg1, ref2, reg2, offset);
    reduce();
  }

  virtual void ref_assume(const reference_constraint_t &cst) override {
    m_product.first().ref_assume(cst);
    m_product.second().ref_assume(cst);
    reduce();
  }

  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    m_product.first().ref_to_int(reg, ref_var, int_var);
    m_product.second().ref_to_int(reg, ref_var, int_var);
    reduce();
  }

  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    m_product.first().int_to_ref(int_var, reg, ref_var);
    m_product.second().int_to_ref(int_var, reg, ref_var);
    reduce();
  }
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
		  const variable_t &cond,
		  const variable_or_constant_t &ref1,
		  const boost::optional<variable_t> &rgn1,
		  const variable_or_constant_t &ref2,
		  const boost::optional<variable_t> &rgn2) override {
    m_product.first().select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
    m_product.second().select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
    reduce();
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    return m_product.first().is_null_ref(ref) & m_product.second().is_null_ref(ref);
  }

  bool get_allocation_sites(const variable_t &ref,
			    std::vector<allocation_site> &out) override {
    std::vector<allocation_site> s1, s2;
    bool b1 = m_product.first().get_allocation_sites(ref, s1);
    bool b2 = m_product.first().get_allocation_sites(ref, s2);
    if (b1 && b2) {
      std::sort(s1.begin(), s1.end());
      std::sort(s2.begin(), s2.end());
      std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
			    std::back_inserter(out));
      return true;
    } else if (b1) {
      out.assign(s1.begin(), s1.end());
      return true;
    } else if (b2) {
      out.assign(s2.begin(), s2.end());
      return true;
    }
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
		std::vector<uint64_t> &out) override {
    std::vector<uint64_t> s1, s2;
    bool b1 = m_product.first().get_tags(rgn, ref, s1);
    bool b2 = m_product.first().get_tags(rgn, ref, s2);
    if (b1 && b2) {
      std::sort(s1.begin(), s1.end());
      std::sort(s2.begin(), s2.end());
      std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
			    std::back_inserter(out));
      return true;
    } else if (b1) {
      out.assign(s1.begin(), s1.end());
      return true;
    } else if (b2) {
      out.assign(s2.begin(), s2.end());
      return true;
    }
    return false;
  }
  
  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    m_product.first().assign_bool_cst(lhs, rhs);
    m_product.second().assign_bool_cst(lhs, rhs);
    reduce();
  }

  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) override {
    m_product.first().assign_bool_ref_cst(lhs, rhs);
    m_product.second().assign_bool_ref_cst(lhs, rhs);
    reduce();
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    m_product.first().assign_bool_var(lhs, rhs, is_not_rhs);
    m_product.second().assign_bool_var(lhs, rhs, is_not_rhs);
    reduce();
  }

  virtual void weak_assign_bool_cst(const variable_t &lhs,
				    const linear_constraint_t &rhs) override {
    m_product.first().weak_assign_bool_cst(lhs, rhs);
    m_product.second().weak_assign_bool_cst(lhs, rhs);
    reduce();
  }

  virtual void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
				    bool is_not_rhs) override {
    m_product.first().weak_assign_bool_var(lhs, rhs, is_not_rhs);
    m_product.second().weak_assign_bool_var(lhs, rhs, is_not_rhs);
    reduce();
  }
  
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    m_product.first().apply_binary_bool(op, x, y, z);
    m_product.second().apply_binary_bool(op, x, y, z);
    reduce();
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    m_product.first().assume_bool(v, is_negated);
    m_product.second().assume_bool(v, is_negated);
    reduce();
  }

  virtual void select_bool(const variable_t &lhs, const variable_t &cond,
			   const variable_t &b1, const variable_t &b2) override {
    m_product.first().select_bool(lhs, cond, b1, b2);
    m_product.second().select_bool(lhs, cond, b1, b2);
    reduce();
  }
  
  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const reduced_domain_product2_t &inv) override {
    m_product.first().backward_assign_bool_cst(lhs, rhs, inv.first());
    m_product.second().backward_assign_bool_cst(lhs, rhs, inv.second());
    reduce();
  }

  virtual void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const reduced_domain_product2_t &inv) override {
    m_product.first().backward_assign_bool_ref_cst(lhs, rhs, inv.first());
    m_product.second().backward_assign_bool_ref_cst(lhs, rhs,
                                                         inv.second());
    reduce();
  }

  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const reduced_domain_product2_t &inv) override {
    m_product.first().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                    inv.first());
    m_product.second().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                     inv.second());
    reduce();
  }

  virtual void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const reduced_domain_product2_t &inv) override {
    m_product.first().backward_apply_binary_bool(op, x, y, z, inv.first());
    m_product.second().backward_apply_binary_bool(op, x, y, z,
                                                       inv.second());
    reduce();
  }

  virtual void forget(const variable_vector_t &variables) override {
    m_product.first().forget(variables);
    m_product.second().forget(variables);
  }

  virtual void project(const variable_vector_t &variables) override {
    m_product.first().project(variables);
    m_product.second().project(variables);
  }

  virtual void expand(const variable_t &var,
                      const variable_t &new_var) override {
    m_product.first().expand(var, new_var);
    m_product.second().expand(var, new_var);
  }

  virtual void normalize() override {
    m_product.first().normalize();
    m_product.second().normalize();
  }

  virtual void minimize() override {
    m_product.first().minimize();
    m_product.second().minimize();
  }

  virtual interval_t operator[](const variable_t &v) override {
    return m_product.first()[v] & m_product.second()[v];
  }

  virtual interval_t at(const variable_t &v) const override {
    return m_product.first().at(v) & m_product.second().at(v);
  }  

  virtual linear_constraint_system_t
  to_linear_constraint_system() const override {
    linear_constraint_system_t csts;
    // XXX: We might add redundant constraints.
    csts += m_product.first().to_linear_constraint_system();
    csts += m_product.second().to_linear_constraint_system();
    return csts;
  }

  virtual disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t csts;
    // XXX: We might add redundant constraints.
    csts += m_product.first().to_disjunctive_linear_constraint_system();
    csts += m_product.second().to_disjunctive_linear_constraint_system();
    return csts;
  }

  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) override {
    m_product.first().rename(from, to);
    m_product.second().rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    m_product.first().intrinsic(name, inputs, outputs);
    m_product.second().intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const reduced_domain_product2_t &invariant) override {
    m_product.first().backward_intrinsic(name, inputs, outputs,
					 invariant.first());
    m_product.second().backward_intrinsic(name, inputs, outputs,
					  invariant.second());
  }
  /* end intrinsics operations */

  void write(crab::crab_os &o) const override { m_product.write(o); }

  std::string domain_name() const override {
    return m_product.domain_name();
  }

}; // class reduced_domain_product2

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

// This domain is similar to reduced_domain_product2 but it combines two
// numerical domains and it defines a more precise, customizable
// reduction operation.
template <typename Domain1, typename Domain2,
          class Params = reduced_product_impl::default_params>
class reduced_numerical_domain_product2 final
    : public abstract_domain_api<
          reduced_numerical_domain_product2<Domain1, Domain2, Params>> {

public:
  using reduced_numerical_domain_product2_t =
      reduced_numerical_domain_product2<Domain1, Domain2, Params>;
  using abstract_domain_t =
      abstract_domain_api<reduced_numerical_domain_product2_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using number_t = typename Domain1::number_t;
  using varname_t = typename Domain1::varname_t;

  static_assert(std::is_same<number_t, typename Domain2::number_t>::value,
                "Domain1 and Domain2 must have same type for number_t");
  static_assert(std::is_same<varname_t, typename Domain2::varname_t>::value,
                "Domain1 and Domain2 must have same type for varname_t");

private:
  using reduced_domain_product2_t =
      reduced_domain_product2<number_t, varname_t, Domain1, Domain2>;

  reduced_domain_product2_t m_product;

  reduced_numerical_domain_product2(const reduced_domain_product2_t &product)
      : m_product(product) {}

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

      Domain1 &inv1 = m_product.first();
      Domain2 &inv2 = m_product.second();

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
                        m_product.first().domain_name());
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
                        m_product.second().domain_name());
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
                      m_product.first().domain_name());
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
                      m_product.second().domain_name());
        crab::CrabStats::uset(k, crab::CrabStats::get(k) + csts2.size());
        inv1 += csts2;
      }
    }
  }

public:
  reduced_numerical_domain_product2_t make_top() const override {
    reduced_domain_product2_t dom_prod;
    return reduced_numerical_domain_product2_t(dom_prod.make_top());
  }

  reduced_numerical_domain_product2_t make_bottom() const override {
    reduced_domain_product2_t dom_prod;
    return reduced_numerical_domain_product2_t(dom_prod.make_bottom());
  }

  void set_to_top() override {
    reduced_domain_product2_t dom_prod;
    reduced_numerical_domain_product2_t abs(dom_prod.make_top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reduced_domain_product2_t dom_prod;
    reduced_numerical_domain_product2_t abs(dom_prod.make_bottom());
    std::swap(*this, abs);
  }

  reduced_numerical_domain_product2()
    : m_product() {}

  reduced_numerical_domain_product2(Domain1 &&val1, Domain2 &&val2)
    : m_product(std::move(val1), std::move(val2)) {}
    
  reduced_numerical_domain_product2(
      const reduced_numerical_domain_product2_t &other)
    : m_product(other.m_product) {}

  reduced_numerical_domain_product2(
      reduced_numerical_domain_product2_t &&other)
    : m_product(std::move(other.m_product)) {}

  
  reduced_numerical_domain_product2_t &
  operator=(const reduced_numerical_domain_product2_t &other) {
    if (this != &other) {
      m_product = other.m_product;
    }
    return *this;
  }

  reduced_numerical_domain_product2_t &
  operator=(reduced_numerical_domain_product2_t &&other) {
    if (this != &other) {
      m_product = std::move(other.m_product);
    }
    return *this;
  }
  
  bool is_bottom() const override { return m_product.is_bottom(); }

  bool is_top() const override { return m_product.is_top(); }

  bool
  operator<=(const reduced_numerical_domain_product2_t &other) const override {
    return m_product <= other.m_product;
  }

  void operator|=(const reduced_numerical_domain_product2_t &other) override {
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";);
    m_product |= other.m_product;
    CRAB_LOG("reduced-dom", crab::outs() << *this << "\n----------------\n";);
  }

  reduced_numerical_domain_product2_t
  operator|(const reduced_numerical_domain_product2_t &other) const override {
    reduced_numerical_domain_product2_t res(m_product | other.m_product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ JOIN ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";
             crab::outs() << res << "\n================\n");
    return res;
  }

  void operator&=(const reduced_numerical_domain_product2_t &other) override {
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ MEET ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------";);
    m_product &= other.m_product;
    CRAB_LOG("reduced-dom", crab::outs() << *this << "\n----------------\n";);
  }
  
  reduced_numerical_domain_product2_t
  operator&(const reduced_numerical_domain_product2_t &other) const override {
    reduced_numerical_domain_product2_t res(m_product & other.m_product);
    CRAB_LOG("combined-domain", crab::outs()
                                    << "============ MEET ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator||(const reduced_numerical_domain_product2_t &other) const override {
    reduced_numerical_domain_product2_t res(m_product || other.m_product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t widening_thresholds(
      const reduced_numerical_domain_product2_t &other,
      const thresholds<number_t> &ts) const override {
    reduced_numerical_domain_product2_t res(
        m_product.widening_thresholds(other.m_product, ts));
    CRAB_LOG("combined-domain",
             crab::outs() << "============ WIDENING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  reduced_numerical_domain_product2_t
  operator&&(const reduced_numerical_domain_product2_t &other) const override {
    reduced_numerical_domain_product2_t res(m_product && other.m_product);
    CRAB_LOG("combined-domain",
             crab::outs() << "============ NARROWING ==================";
             crab::outs() << *this << "\n----------------";
             crab::outs() << other << "\n----------------\n";);
    return res;
  }

  void set(const variable_t &v, interval_t x) {
    m_product.first().set(v, x);
    m_product.second().set(v, x);
  }

  interval_t operator[](const variable_t &v) override {
    return m_product.first()[v] & m_product.second()[v];
  }

  interval_t at(const variable_t &v) const override {
    return m_product.first().at(v) & m_product.second().at(v);
  }  

  void operator+=(const linear_constraint_system_t &csts) override {
    m_product += csts;
    if (!is_bottom()) {
      for (auto const &cst : csts) {
	for (auto const &v : cst.variables()) {
	  reduce_variable(v);
	  if (is_bottom()) {
	    return;
	  }
	}
      }
    }
    CRAB_LOG("combined-domain", crab::outs() << "Added constraints " << csts
                                             << "=" << *this << "\n");
  }


  bool entails(const linear_constraint_t &cst) const override {
    if (m_product.first().entails(cst)) {
      return true;
    }
    if (m_product.second().entails(cst)) {
      return true;
    }
    return false;
  }
  
  void operator-=(const variable_t &v) override { m_product -= v; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.assign(x, e);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
    CRAB_LOG("combined-domain", crab::outs()
                                    << x << ":=" << e << "=" << *this << "\n");
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.weak_assign(x, e);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
    CRAB_LOG("combined-domain", crab::outs()
	     << "weak_assign(" << x << "," << e << ")=" << *this << "\n");
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << z << "=" << *this << "\n");
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
    CRAB_LOG("combined-domain",
             crab::outs() << x << ":=" << y << op << k << "=" << *this << "\n");
  }

  void backward_assign(
      const variable_t &x, const linear_expression_t &e,
      const reduced_numerical_domain_product2_t &invariant) override {
    m_product.backward_assign(x, e, invariant.m_product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      for (auto const &v : e.variables())
        reduce_variable(v);
    }
  }

  void backward_apply(
      arith_operation_t op, const variable_t &x, const variable_t &y,
      number_t k,
      const reduced_numerical_domain_product2_t &invariant) override {
    m_product.backward_apply(op, x, y, k, invariant.m_product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      reduce_variable(y);
    }
  }

  void backward_apply(
      arith_operation_t op, const variable_t &x, const variable_t &y,
      const variable_t &z,
      const reduced_numerical_domain_product2_t &invariant) override {
    m_product.backward_apply(op, x, y, z, invariant.m_product);
    if (!Params::apply_reduction_only_add_constraint) {
      // reduce the variables in the right-hand side
      reduce_variable(y);
      reduce_variable(z);
    }
  }

  // cast operators

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    m_product.apply(op, dst, src);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(dst);
    }
  }

  // bitwise operators

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(x);
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    m_product.select(lhs, cond, e1, e2);
    if (!Params::apply_reduction_only_add_constraint) {
      reduce_variable(lhs);
    }
  }
  
  /// reduced_numerical_domain_product2 implements only standard
  /// abstract operations of a numerical domain so it is intended to be
  /// used as a leaf domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(reduced_numerical_domain_product2_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(reduced_numerical_domain_product2_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(
      reduced_numerical_domain_product2_t)

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    m_product.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    m_product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(
      std::string name,
      const variable_or_constant_vector_t &inputs,
      const variable_vector_t &outputs,
      const reduced_numerical_domain_product2_t &invariant) override {
    m_product.backward_intrinsic(name, inputs, outputs,
				 invariant.m_product);
  }
  /* end intrinsics operations */

  void forget(const variable_vector_t &variables) override {
    m_product.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    m_product.project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    m_product.expand(var, new_var);
  }

  void normalize() override { m_product.normalize(); }

  void minimize() override { m_product.minimize(); }

  void write(crab_os &o) const override { m_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() const override {
    return m_product.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return m_product.to_disjunctive_linear_constraint_system();
  }

  std::string domain_name() const override {
    const Domain1 &dom1 = m_product.first();
    const Domain2 &dom2 = m_product.second();

    std::string name =
        "ReducedProduct(" + dom1.domain_name() + "," + dom2.domain_name() + ")";
    return name;
  }

}; // class reduced_numerical_domain_product2


template <typename Number, typename VariableName, typename Domain1,
          typename Domain2>
struct abstract_domain_traits<
    reduced_domain_product2<Number, VariableName, Domain1, Domain2>> {
  using number_t = Number;
  using varname_t = VariableName;
};

template <typename Domain1, typename Domain2, class Params>
struct abstract_domain_traits<
    reduced_numerical_domain_product2<Domain1, Domain2, Params>> {
  using number_t = typename Domain1::number_t;
  using varname_t = typename Domain1::varname_t;
};

} // end namespace domains
} // namespace crab
