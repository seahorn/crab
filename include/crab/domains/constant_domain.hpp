#pragma once

/* Classical constant propagation domain */

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/constant.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
/**
 *  Each variable is mapped to an element in this lattice:
 *
 *            top
 *             |
 * ...,-3,-2,-1,0,1,2,3,...
 *             |
 *           bottom
 *
 * top means that it might not be a constant value
 **/
template <typename Number, typename VariableName>
class constant_domain final : public crab::domains::abstract_domain_api<
                                  constant_domain<Number, VariableName>> {
public:
  using constant_domain_t = constant_domain<Number, VariableName>;
  using abstract_domain_t =
      crab::domains::abstract_domain_api<constant_domain_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using constant_t = constant<Number>;
  using number_t = Number;
  using varname_t = VariableName;

private:
  using separate_domain_t = ikos::separate_domain<variable_t, constant_t>;

private:
  separate_domain_t m_env;

  constant_domain(separate_domain_t &&env) : m_env(std::move(env)) {}

  constant_t eval(const linear_expression_t &expr) const {
    assert(!is_bottom());
    constant_t r(expr.constant());
    for (auto const &kv : expr) {
      constant_t c(kv.first);
      r = r.Add(c.Mul(m_env.at(kv.second)));
      if (r.is_top()) {
        break;
      }
    }
    return r;
  }

  constant_t compute_residual(const linear_constraint_t &cst,
                              const variable_t &pivot) const {
    constant_t residual(cst.constant());
    for (auto const &kv : cst) {
      constant_t c(kv.first);
      const variable_t &v = kv.second;
      if (!(v == pivot)) {
        residual = residual.Sub(c.Mul(m_env.at(v)));
        if (residual.is_top()) {
          break;
        }
      }
    }
    return residual;
  }

  void propagate(const linear_constraint_t &cst) {
    if (is_bottom()) {
      return;
    }

    if (cst.is_inequality() || cst.is_strict_inequality() ||
        cst.is_disequation()) {
      constant_t e = eval(cst.expression());
      if (e.is_constant()) {
        if (cst.is_inequality()) {
          if (!(e.get_constant() <= number_t(0))) {
            set_to_bottom();
            return;
          }
        } else if (cst.is_disequation()) {
          if (!(e.get_constant() != number_t(0))) {
            set_to_bottom();
            return;
          }
        } else if (cst.is_strict_inequality()) {
          if (!(e.get_constant() < number_t(0))) {
            set_to_bottom();
            return;
          }
        }
      }
    } else if (cst.is_equality()) {
      for (auto kv : cst) {
        number_t c = kv.first;
        const variable_t &pivot = kv.second;
        constant_t new_c = compute_residual(cst, pivot).SDiv(c);
        if (!new_c.is_top()) {
          m_env.set(pivot, m_env.at(pivot) & new_c);
        }
      }
    }
  }

  void solve_constraints(const linear_constraint_system_t &csts) {
    for (auto const &c : csts) {
      if (is_bottom()) {
        return;
      }
      if (c.is_tautology()) {
        continue;
      }
      if (c.is_contradiction()) {
        set_to_bottom();
        return;
      }
      propagate(c);
    }
  }

public:
  constant_domain_t make_top() const override {
    return constant_domain_t(separate_domain_t::top());
  }

  constant_domain_t make_bottom() const override {
    return constant_domain_t(separate_domain_t::bottom());
  }

  void set_to_top() override {
    constant_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    constant_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  constant_domain() : m_env(separate_domain_t::top()) {}

  constant_domain(const constant_domain_t &e) : m_env(e.m_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  constant_domain(constant_domain_t &&e) : m_env(std::move(e.m_env)) {}

  constant_domain_t &operator=(const constant_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_env = o.m_env;
    }
    return *this;
  }

  constant_domain_t &operator=(constant_domain_t &&o) {
    if (this != &o) {
      m_env = std::move(o.m_env);
    }
    return *this;
  }

  constant_t get_constant(const variable_t &v) const { return m_env.at(v); }

  void set_constant(const variable_t &v, constant_t c) { m_env.set(v, c); }

  bool is_bottom() const override { return m_env.is_bottom(); }

  bool is_top() const override { return m_env.is_top(); }

  bool operator<=(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (m_env <= o.m_env);
  }

  void operator|=(const constant_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    CRAB_LOG("constant-domain",
             crab::outs() << "Join " << m_env << " and " << o.m_env << "\n";);
    m_env = m_env | o.m_env;
    CRAB_LOG("constant-domain", crab::outs() << "Res=" << m_env << "\n";);
  }

  constant_domain_t operator|(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    return (m_env | o.m_env);
  }

  void operator&=(const constant_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    CRAB_LOG("constant-domain",
             crab::outs() << "Meet " << m_env << " and " << o.m_env << "\n";);
    m_env = m_env & o.m_env;
    CRAB_LOG("constant-domain", crab::outs() << "Res=" << m_env << "\n";);
  }
  
  constant_domain_t operator&(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    return (m_env & o.m_env);
  }

  constant_domain_t operator||(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return (m_env || o.m_env);
  }

  constant_domain_t
  widening_thresholds(const constant_domain_t &o,
                      const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return m_env.widening_thresholds(o.m_env, ts);
  }

  constant_domain_t operator&&(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (m_env && o.m_env);
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    if (!is_bottom()) {
      m_env -= v;
    }
  }

  interval_t operator[](const variable_t &v) override { return at(v); }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      constant_t c = m_env.at(v);
      if (c.is_bottom()) {
	return interval_t::bottom();
      } else if (c.is_top()) {
	return interval_t::top();
      } else {
	assert(c.is_constant());
	return interval_t(c.get_constant());
      }
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    solve_constraints(csts);
  }

  DEFAULT_ENTAILS(constant_domain_t)
  
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    if (!is_bottom()) {
      if (boost::optional<variable_t> v = e.get_variable()) {
	m_env.set(x, m_env.at(*v));
      } else {
	m_env.set(x, eval(e));
      }
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.weak_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".weak_assign");
    if (!is_bottom()) {    
      if (boost::optional<variable_t> v = e.get_variable()) {
	m_env.join(x, m_env.at(*v));
      } else {
	m_env.join(x, eval(e));
      }
    }
  }
  
  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      constant_t yc = m_env.at(y);
      constant_t zc = m_env.at(z);
      constant_t xc = constant_t::top();
      switch (op) {
      case crab::domains::OP_ADDITION:
        xc = yc.Add(zc);
        break;
      case crab::domains::OP_SUBTRACTION:
        xc = yc.Sub(zc);
        break;
      case crab::domains::OP_MULTIPLICATION:
        xc = yc.Mul(zc);
        break;
      case crab::domains::OP_SDIV:
        xc = yc.SDiv(zc);
        break;
      case crab::domains::OP_SREM:
        xc = yc.SRem(zc);
        break;
      case crab::domains::OP_UDIV:
        xc = yc.UDiv(zc);
        break;
      case crab::domains::OP_UREM:
        xc = yc.URem(zc);
        break;
      }
      m_env.set(x, xc);
    }
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      constant_t yc = m_env.at(y);
      constant_t zc(k);
      constant_t xc = constant_t::top();
      switch (op) {
      case crab::domains::OP_ADDITION:
        xc = yc.Add(zc);
        break;
      case crab::domains::OP_SUBTRACTION:
        xc = yc.Sub(zc);
        break;
      case crab::domains::OP_MULTIPLICATION:
        xc = yc.Mul(zc);
        break;
      case crab::domains::OP_SDIV:
        xc = yc.SDiv(zc);
        break;
      case crab::domains::OP_SREM:
        xc = yc.SRem(zc);
        break;
      case crab::domains::OP_UDIV:
        xc = yc.UDiv(zc);
        break;
      case crab::domains::OP_UREM:
        xc = yc.URem(zc);
        break;
      }
      m_env.set(x, xc);
    }
  }

  // intrinsics operations
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const constant_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  // backward arithmetic operations
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const constant_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");
    // TODO
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const constant_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    // TODO
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const constant_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    // TODO
  }

  // cast operations
  void apply(crab::domains::int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<constant_domain_t>::apply(*this, op, dst, src);
  }

  // bitwise operations
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      constant_t yc = m_env.at(y);
      constant_t zc = m_env.at(z);
      constant_t xc = constant_t::top();
      switch (op) {
      case crab::domains::OP_AND:
        xc = yc.BitwiseAnd(zc);
        break;
      case crab::domains::OP_OR:
        xc = yc.BitwiseOr(zc);
        break;
      case crab::domains::OP_XOR:
        xc = yc.BitwiseXor(zc);
        break;
      case crab::domains::OP_SHL:
        xc = yc.BitwiseShl(zc);
        break;
      case crab::domains::OP_LSHR:
        xc = yc.BitwiseLShr(zc);
        break;
      case crab::domains::OP_ASHR: 
        xc = yc.BitwiseAShr(zc);
        break;
      }
      m_env.set(x, xc);
    }
  }

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      constant_t yc = m_env.at(y);
      constant_t zc(k);
      constant_t xc = constant_t::top();
      switch (op) {
      case crab::domains::OP_AND:
        xc = yc.BitwiseAnd(zc);
        break;
      case crab::domains::OP_OR:
        xc = yc.BitwiseOr(zc);
        break;
      case crab::domains::OP_XOR:
        xc = yc.BitwiseXor(zc);
        break;
      case crab::domains::OP_SHL:
        xc = yc.BitwiseShl(zc);
        break;
      case crab::domains::OP_LSHR:
        xc = yc.BitwiseLShr(zc);
        break;
      case crab::domains::OP_ASHR: 
        xc = yc.BitwiseAShr(zc);
        break;
      }
      m_env.set(x, xc);
    }
  }

  virtual void select(const variable_t &lhs, const linear_constraint_t &cond,
                      const linear_expression_t &e1,
                      const linear_expression_t &e2) override {
    crab::CrabStats::count(domain_name() + ".count.select");
    crab::ScopedCrabStats __st__(domain_name() + ".select");

    if (!is_bottom()) {
      constant_domain_t inv1(*this);
      inv1 += cond;
      if (inv1.is_bottom()) {
        assign(lhs, e2);
        return;
      }
      constant_domain_t inv2(*this);
      inv2 += cond.negate();
      if (inv2.is_bottom()) {
        assign(lhs, e1);
        return;
      }
      m_env.set(lhs, eval(e1) | eval(e2));
    }
  }

  /// constant_domain implements only standard abstract operations of
  /// a numerical domain so it is intended to be used as a leaf domain
  /// in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(constant_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(constant_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(constant_domain_t)

    
  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }
    for (auto const &var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");
    if (!is_bottom()) {
      m_env.project(variables);
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");
    if (!is_bottom()) {
      m_env.rename(from, to);
    }
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    m_env.set(new_x, m_env.at(x));
  }

  void normalize() override {}

  void minimize() override {}

  void write(crab::crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    m_env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;

    if (this->is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    for (auto it = m_env.begin(); it != m_env.end(); ++it) {
      const variable_t &v = it->first;
      const constant_t &c = it->second;
      if (c.is_constant()) {
        csts += linear_constraint_t(v == c.get_constant());
      }
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  std::string domain_name() const override { return "ConstantDomain"; }

}; // class constant_domain
} // namespace domains
} // namespace crab

namespace crab {
namespace domains {
template <typename Number, typename VariableName>
struct abstract_domain_traits<constant_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // namespace domains
} // namespace crab
