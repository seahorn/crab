#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <string>

namespace crab {
namespace domains {

/*
 * Lookahead widening domain based on the paper "Lookahead Widening"
 * by D. Gopan and T. Reps published in CAV 2006.
 */
template <class Dom>
class lookahead_widening_domain
    : public abstract_domain_api<lookahead_widening_domain<Dom>> {
public:
  using this_type = lookahead_widening_domain<Dom>;
  using abstract_domain_api_t = abstract_domain_api<this_type>;
  using typename abstract_domain_api_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_api_t::interval_t;
  using typename abstract_domain_api_t::linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_t;
  using typename abstract_domain_api_t::linear_expression_t;
  using typename abstract_domain_api_t::number_t;
  using typename abstract_domain_api_t::reference_constraint_t;
  using typename abstract_domain_api_t::variable_or_constant_t;
  using typename abstract_domain_api_t::variable_or_constant_vector_t;
  using typename abstract_domain_api_t::variable_t;
  using typename abstract_domain_api_t::variable_vector_t;
  using typename abstract_domain_api_t::varname_t;

private:
  using product_domain_t =
      reduced_domain_product2<number_t, varname_t, Dom, Dom>;

  /**
   * The first value is used to keep the analysis within the current
   * loop phase: this value is used to decide which branch to take and
   * is never widened. The second value is used to compute the
   * solution for the current phase: both widening and narrowing are
   * applied to it. When the second value stabilizes, it is promoted
   * into the first value.

   * Invariants:
   * if m_is_bottom is false then
   *    m_product.first() is less or equal than m_product.second()
   *    m_product.first() is not bottom
  **/

  // TODO(OPTIMIZATION): keep only one abstract value if first and
  // second are the same.
  product_domain_t m_product;
  bool m_is_bottom;

  lookahead_widening_domain(product_domain_t &&product)
      : m_product(std::move(product)), m_is_bottom(false) {
    if (m_product.first().is_bottom()) {
      set_to_bottom();
    }
  }

public:
  lookahead_widening_domain(bool is_bottom = false)
      : m_product(product_domain_t()), m_is_bottom(is_bottom) {}

  lookahead_widening_domain(const this_type &o) = default;
  lookahead_widening_domain(this_type &&o) = default;
  this_type &operator=(const this_type &o) = default;
  this_type &operator=(this_type &&o) = default;

  void set_to_top() override {
    m_product.set_to_top();
    m_is_bottom = false;
  }

  void set_to_bottom() override {
    m_product.set_to_top();
    m_is_bottom = true;
  }

  this_type make_bottom() const override {
    this_type res(true);
    return res;
  }

  this_type make_top() const override {
    this_type res(false);
    return res;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override { return !m_is_bottom && m_product.is_top(); }

  bool operator<=(const this_type &other) const override {
    auto compare = [](const Dom &d1, const Dom &d2) -> unsigned {
      // return 0 if d1 == d2
      // return 1 if d1 < d2
      // return 2 otherwise
      bool v1 = d1 <= d2;
      bool v2 = d2 <= d1;
      if (v1 && v2) {
        return 0;
      } else if (v1 && !v2) {
        return 1;
      } else {
        return 2;
      }
    };

    if (is_bottom() || other.is_top()) {
      return true;
    } else if (is_top() || other.is_bottom()) {
      return false;
    } else {
      // Lexicographical ordering
      unsigned res = compare(m_product.first(), other.m_product.first());
      return (res == 1 ||
              (res == 0 && (m_product.second() <= other.m_product.second())));
    }
  }

  void operator|=(const this_type &other) override {
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      m_product |= other.m_product;
    }
  }

  this_type operator|(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      return this_type(m_product | other.m_product);
    }
  }

  void operator&=(const this_type &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      m_product &= other.m_product;
    }
  }

  this_type operator&(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      return this_type(m_product & other.m_product);
    }
  }

  this_type operator||(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      // As defined in Section 3.2, page 7.
      //
      // There is an alternative definition in Section 4 for
      // "accumulating" analyzers which is not the case of Crab.
      if (other.m_product.second() <= m_product.second()) {
        // end of phase: second stabilizes so it is promoted to first
        Dom first = other.m_product.second();
        Dom second = other.m_product.second();
        product_domain_t product(std::move(first), std::move(second));
        return this_type(std::move(product));
      } else {
        Dom first = m_product.first() | other.m_product.first();
        Dom second = m_product.second() || other.m_product.second();
        product_domain_t product(std::move(first), std::move(second));
        return this_type(std::move(product));
      }
    }
  }

  this_type widening_thresholds(const this_type &other,
                                const thresholds<number_t> &ts) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      if (other.m_product.second() <= m_product.second()) {
        // end of phase: second stabilizes so it is promoted to first
        Dom first = other.m_product.second();
        Dom second = other.m_product.second();
        product_domain_t product(std::move(first), std::move(second));
        return this_type(std::move(product));
      } else {
        Dom first = m_product.first() | other.m_product.first();
        Dom second = m_product.second().widening_thresholds(
            other.m_product.second(), ts);
        product_domain_t product(std::move(first), std::move(second));
        return this_type(std::move(product));
      }
    }
  }

  this_type operator&&(const this_type &other) const override {
    // Be careful with infinite descending chains
    return (*this & other);
  }

  void operator-=(const variable_t &var) override { m_product -= var; }

  interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    return m_product.first()[v];
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    return m_product.first().at(v);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (!is_bottom()) {
      m_product.first() += csts;
      if (m_product.first().is_bottom()) {
        set_to_bottom();
        return;
      }
      m_product.second() += csts;
      // second cannot be bottom because first is not bottom and
      // second is always an overapproximation of first.
      assert(!m_product.second().is_bottom());
    }
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    } else if (cst.is_tautology()) {
      return true;
    } else if (cst.is_contradiction()) {
      return false;
    }

    return m_product.first().entails(cst);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      m_product.assign(x, e);
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      m_product.weak_assign(x, e);
    }
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      m_product.apply(op, x, y, z);
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_product.apply(op, x, y, z);
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      m_product.apply(op, dst, src);
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_product.apply(op, x, y, z);
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      m_product.apply(op, x, y, k);
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond, 
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    if (!is_bottom()) {
      m_product.select(lhs, cond, e1, e2);
    }
  }
  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  // TODO: we should implement these operations
  BOOL_OPERATIONS_NOT_IMPLEMENTED(this_type)
  // TODO: we should implement these operations
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(this_type)
  // TODO: we should implement these operations
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(this_type)

  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      return linear_constraint_system_t(linear_constraint_t::get_false());
    } else if (is_top()) {
      return linear_constraint_system_t(linear_constraint_t::get_true());
    } else {
      return m_product.first().to_linear_constraint_system();
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    if (is_bottom()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else {
      return m_product.first().to_disjunctive_linear_constraint_system();
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (!(is_bottom() || is_top())) {
      m_product.forget(variables);
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      m_product.project(variables);
    }
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (!is_bottom()) {
      m_product.expand(var, new_var);
    }
  }

  void normalize() override {}
  void minimize() override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      m_product.rename(from, to);
    }
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    if (!is_bottom()) {
      m_product.intrinsic(name, inputs, outputs);
    }
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic not implemented");
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      o << m_product.first();
    }
  }

  friend crab_os &operator<<(crab_os &o, const this_type &absval) {
    absval.write(o);
    return o;
  }

  std::string domain_name() const override {
    Dom absval;
    return "LookaheadWideningDomain(" + absval.domain_name() + ")";
  }
};

template <typename Domain>
struct abstract_domain_traits<lookahead_widening_domain<Domain>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

} // end namespace domains
} // end namespace crab
