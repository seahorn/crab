#pragma once

/**
 *  Reduced product of sign and constant domains.
 */

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/constant_domain.hpp>
#include <crab/domains/sign_domain.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
class sign_constant_domain final
    : public abstract_domain_api<sign_constant_domain<Number, VariableName>> {

  using signed_constant_domain_t = sign_constant_domain<Number, VariableName>;
  using abstract_domain_t = abstract_domain_api<signed_constant_domain_t>;

public:
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

  using constant_domain_t = constant_domain<number_t, varname_t>;
  using sign_domain_t = sign_domain<number_t, varname_t>;
  using constant_t = typename constant_domain_t::constant_t;
  using sign_t = typename sign_domain_t::sign_t;

private:
  using reduced_domain_product2_t =
      reduced_domain_product2<number_t, varname_t, sign_domain_t, constant_domain_t>;

  reduced_domain_product2_t m_product;

  sign_constant_domain(reduced_domain_product2_t &&product)
      : m_product(std::move(product)) {}

  void reduce_variable(const variable_t &v) {
    crab::CrabStats::count(domain_name() + ".count.reduce");
    crab::ScopedCrabStats __st__(domain_name() + ".reduce");

    if (is_bottom()) {
      return;
    }

    sign_domain_t &sign_dom = m_product.first();
    constant_domain_t &constant_dom = m_product.second();

    sign_t s = sign_dom.get_sign(v);
    if (s.equal_zero()) {
      constant_dom.set_constant(v,
                                constant_dom.get_constant(v) & constant_t::zero());
    }
    constant_t c = constant_dom.get_constant(v);
    if (!c.is_constant()) {
      return;
    }

    sign_dom.set_sign(v, s & sign_t(c.get_constant()));
  }

public:
  signed_constant_domain_t make_top() const override {
    reduced_domain_product2_t dom_prod;
    return signed_constant_domain_t(dom_prod.make_top());
  }

  signed_constant_domain_t make_bottom() const override {
    reduced_domain_product2_t dom_prod;
    return signed_constant_domain_t(dom_prod.make_bottom());
  }

  void set_to_top() override {
    reduced_domain_product2_t dom_prod;
    signed_constant_domain_t abs(dom_prod.make_top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reduced_domain_product2_t dom_prod;
    signed_constant_domain_t abs(dom_prod.make_bottom());
    std::swap(*this, abs);
  }

  sign_constant_domain() : m_product() {}

  sign_constant_domain(const signed_constant_domain_t &other) = default;

  sign_constant_domain(signed_constant_domain_t &&other) = default;

  signed_constant_domain_t &
  operator=(const signed_constant_domain_t &other) = default;

  signed_constant_domain_t &
  operator=(signed_constant_domain_t &&other) = default;

  bool is_bottom() const override { return m_product.is_bottom(); }

  bool is_top() const override { return m_product.is_top(); }

  bool operator<=(const signed_constant_domain_t &other) const override {
    return m_product <= other.m_product;
  }

  void operator|=(const signed_constant_domain_t &other) override {
    m_product |= other.m_product;
  }

  signed_constant_domain_t
  operator|(const signed_constant_domain_t &other) const override {
    return signed_constant_domain_t(m_product | other.m_product);
  }

  void operator&=(const signed_constant_domain_t &other) override {
    m_product &= other.m_product;
  }

  signed_constant_domain_t
  operator&(const signed_constant_domain_t &other) const override {
    return signed_constant_domain_t(m_product & other.m_product);
  }

  signed_constant_domain_t
  operator||(const signed_constant_domain_t &other) const override {
    return signed_constant_domain_t(m_product || other.m_product);
  }

  signed_constant_domain_t widening_thresholds(
      const signed_constant_domain_t &other,
      const thresholds<number_t> &ts) const override {
    return signed_constant_domain_t(
        m_product.widening_thresholds(other.m_product, ts));
  }

  signed_constant_domain_t
  operator&&(const signed_constant_domain_t &other) const override {
    return signed_constant_domain_t(m_product && other.m_product);
  }

  interval_t operator[](const variable_t &v) override {
    return m_product.first()[v] & m_product.second()[v];
  }

  interval_t at(const variable_t &v) const override {
    return m_product.first().at(v) & m_product.second().at(v);
  }
  
  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom())
      return;

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
  }

  DEFAULT_ENTAILS(signed_constant_domain_t)
  
  void operator-=(const variable_t &v) override { m_product -= v; }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.assign(x, e);
    reduce_variable(x);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.weak_assign(x, e);
    reduce_variable(x);
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    reduce_variable(x);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    reduce_variable(x);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const signed_constant_domain_t &invariant) override {
    m_product.backward_assign(x, e, invariant.m_product);
    // reduce the variables in the right-hand side
    for (auto const &v : e.variables()) {
      reduce_variable(v);
    }
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const signed_constant_domain_t &invariant) override {
    m_product.backward_apply(op, x, y, k, invariant.m_product);
    // reduce the variables in the right-hand side
    reduce_variable(y);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const signed_constant_domain_t &invariant) override {
    m_product.backward_apply(op, x, y, z, invariant.m_product);
    // reduce the variables in the right-hand side
    reduce_variable(y);
    reduce_variable(z);
  }

  // cast operators

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    m_product.apply(op, dst, src);
    reduce_variable(dst);
  }

  // bitwise operators

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    reduce_variable(x);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    reduce_variable(x);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    m_product.select(lhs, cond, e1, e2);
    reduce_variable(lhs);
  }

  /// sign_constant_domain implements only standard abstract
  /// operations of a numerical domain so it is intended to be used as
  /// a leaf domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(signed_constant_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(signed_constant_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(signed_constant_domain_t)

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

  std::string domain_name() const override { return m_product.domain_name(); }

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

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const signed_constant_domain_t &invariant) override {
    m_product.backward_intrinsic(name, inputs, outputs, invariant.m_product);
  }
  /* end intrinsics operations */

}; // class sign_constant_domain
} // namespace domains
} // namespace crab

namespace crab {
namespace domains {
template <typename Number, typename VariableName>
struct abstract_domain_traits<sign_constant_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // namespace domains
} // namespace crab
