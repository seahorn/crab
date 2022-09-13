#pragma once

/**
 * Reduced product of a numerical domain and the congruence domain.
 **/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/interval_congruence.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {

// Reduced product of a numerical domain with interval x congruences.
template <typename NumAbsDom>
class numerical_congruence_domain final
    : public abstract_domain_api<numerical_congruence_domain<NumAbsDom>> {

  using rnc_domain_t = numerical_congruence_domain<NumAbsDom>;
  using abstract_domain_t = abstract_domain_api<rnc_domain_t>;

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
  using number_t = typename NumAbsDom::number_t;
  using varname_t = typename NumAbsDom::varname_t;

  using congruence_domain_t = ikos::congruence_domain<number_t, varname_t>;
  using interval_congruence_t = interval_congruence<number_t>;

private:
  using reduced_domain_product2_t =
      reduced_domain_product2<number_t, varname_t, NumAbsDom, congruence_domain_t>;

  reduced_domain_product2_t m_product;

  numerical_congruence_domain(const reduced_domain_product2_t &product)
      : m_product(product) {}

  void reduce_variable(const variable_t &v) {
    crab::CrabStats::count(domain_name() + ".count.reduce");
    crab::ScopedCrabStats __st__(domain_name() + ".reduce");

    if (is_bottom()) {
      return;
    }

    auto i = m_product.first()[v]; // project on intervals
    auto c = m_product.second().to_congruence(v);
    interval_congruence_t val(std::move(i), std::move(c));

    if (val.is_bottom()) {
      set_to_bottom();
    } else {
      if (val.first() != i) {
        // FIXME: method set is not part of the abstract_domain API so
        // it might not compile.
        m_product.first().set(v, val.first());
      }

      if (val.second() != c) {
        // FIXME: method set is not part of the abstract_domain API so
        // it might not compile.
        m_product.second().set(v, val.second());
      }
    }
  }

public:
  rnc_domain_t make_top() const override {
    reduced_domain_product2_t dom_prod;
    return rnc_domain_t(dom_prod.make_top());
  }

  rnc_domain_t make_bottom() const override {
    reduced_domain_product2_t dom_prod;
    return rnc_domain_t(dom_prod.make_bottom());
  }

  void set_to_top() override {
    reduced_domain_product2_t dom_prod;
    rnc_domain_t abs(dom_prod.make_top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reduced_domain_product2_t dom_prod;
    rnc_domain_t abs(dom_prod.make_bottom());
    std::swap(*this, abs);
  }

  numerical_congruence_domain() : m_product() {}

  numerical_congruence_domain(const rnc_domain_t &other)
      : m_product(other.m_product) {}

  rnc_domain_t &operator=(const rnc_domain_t &other) {
    if (this != &other)
      m_product = other.m_product;

    return *this;
  }

  bool is_bottom() const override { return m_product.is_bottom(); }

  bool is_top() const override { return m_product.is_top(); }

  bool operator<=(const rnc_domain_t &other) const override {
    return m_product <= other.m_product;
  }

  bool operator==(const rnc_domain_t &other) const {
    return m_product == other.m_product;
  }

  void operator|=(const rnc_domain_t &other) override {
    m_product |= other.m_product;
  }

  rnc_domain_t operator|(const rnc_domain_t &other) const override {
    return rnc_domain_t(m_product | other.m_product);
  }

  void operator&=(const rnc_domain_t &other) override {
    m_product &= other.m_product;
  }
  
  rnc_domain_t operator&(const rnc_domain_t &other) const override {
    return rnc_domain_t(m_product & other.m_product);
  }

  rnc_domain_t operator||(const rnc_domain_t &other) const override {
    return rnc_domain_t(m_product || other.m_product);
  }

  rnc_domain_t widening_thresholds(
      const rnc_domain_t &other,
      const thresholds<number_t> &ts) const override {
    return rnc_domain_t(m_product.widening_thresholds(other.m_product, ts));
  }

  rnc_domain_t operator&&(const rnc_domain_t &other) const override {
    return rnc_domain_t(m_product && other.m_product);
  }

  // pre: x is already reduced
  void set(const variable_t &v, interval_congruence_t x) {
    m_product.first().set(v, x.first());
    m_product.second().set(v, x.second());
  }

  interval_t operator[](const variable_t &v) override {
    return m_product.first()[v];
  }

  interval_t at(const variable_t &v) const override {
    return m_product.first().at(v);
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
  }

  bool entails(const linear_constraint_t &cst) const override {
    return m_product.entails(cst);
  }
  
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
                       const rnc_domain_t &invariant) override {
    m_product.backward_assign(x, e, invariant.m_product);
    // reduce the variables in the right-hand side
    for (auto const &v : e.variables())
      reduce_variable(v);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const rnc_domain_t &invariant) override {
    m_product.backward_apply(op, x, y, k, invariant.m_product);
    // reduce the variables in the right-hand side
    reduce_variable(y);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const rnc_domain_t &invariant) override {
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

  /// numerical_congruence_domain implements only standard abstract
  /// operations of a numerical domain so it is intended to be used as
  /// a leaf domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(rnc_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(rnc_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(rnc_domain_t)

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
                          const rnc_domain_t &invariant) override {
    m_product.backward_intrinsic(name, inputs, outputs, invariant.m_product);
  }
  /* end intrinsics operations */

}; // class numerical_congruence_domain

template <typename NumAbsDom>
struct abstract_domain_traits<numerical_congruence_domain<NumAbsDom>> {
  using number_t = typename NumAbsDom::number_t;
  using varname_t = typename NumAbsDom::varname_t;
};

} // end namespace domains
} // namespace crab
