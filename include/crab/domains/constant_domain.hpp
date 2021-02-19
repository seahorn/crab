#pragma once

/* Classical constant propagation domain */

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
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
template <typename Number> class constant {
  boost::optional<Number> m_constant;
  bool m_is_bottom;
  using constant_t = constant<Number>;

public:
  constant(bool is_bottom) : m_constant(boost::none), m_is_bottom(is_bottom) {}

  constant(Number c) : m_constant(c), m_is_bottom(false) {}

  static constant_t bottom() { return constant(true); }

  static constant_t top() { return constant(false); }

  bool is_bottom() const { return m_is_bottom; }

  bool is_top() const { return (!is_bottom() && !m_constant); }

  bool is_constant() const { return m_constant != boost::none; }

  Number get_constant() const {
    assert(is_constant());
    return *m_constant;
  }

  bool operator<=(const constant_t &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (o.is_bottom() || is_top()) {
      return false;
    } else {
      assert(is_constant());
      assert(o.is_constant());
      return get_constant() == o.get_constant();
    }
  }

  bool operator==(const constant_t &o) const {
    return (m_is_bottom == o.m_is_bottom && m_constant == o.m_constant);
  }

  constant_t operator|(const constant_t &o) const {
    if (is_bottom() || o.is_top())
      return o;
    else if (is_top() || o.is_bottom())
      return *this;
    else {
      assert(is_constant());
      assert(o.is_constant());
      if (get_constant() == o.get_constant()) {
        return *this;
      } else {
        return constant_t::top();
      }
    }
  }

  constant_t operator||(const constant_t &o) const { return *this | o; }

  template <typename Thresholds>
  constant_t widening_thresholds(const constant_t &o,
                                 const Thresholds &ts /*unused*/) const {
    return *this | o;
  }

  constant_t operator&(const constant_t &o) const {
    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      assert(is_constant());
      assert(o.is_constant());
      if (get_constant() == o.get_constant()) {
        return *this;
      } else {
        return constant_t::bottom();
      }
    }
  }

  constant_t operator&&(const constant_t &o) const { return *this & o; }

  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      assert(is_constant());
      o << get_constant();
    }
  }

  friend inline crab_os &operator<<(crab_os &o, const constant_t &c) {
    c.write(o);
    return o;
  }
};

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
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using constant_t = constant<Number>;
  using number_t = Number;
  using varname_t = VariableName;

private:
  using interval_domain_t = ikos::interval_domain<number_t, varname_t>;
  using separate_domain_t = ikos::separate_domain<variable_t, constant_t>;

private:
  separate_domain_t m_env;

  constant_domain(separate_domain_t &&env) : m_env(std::move(env)) {}

  void solve_constraints(const linear_constraint_system_t &csts) {
    if (!is_bottom()) {
      interval_domain_t idom;
      for (auto const &c : csts) {
        if (c.is_inequality() && c.is_unsigned()) {
          // interval domain ignores these constraints
          continue;
        }
        if (c.is_tautology()) {
          continue;
        }
        if (c.is_contradiction()) {
          set_to_bottom();
          return;
        }

        for (auto v : c.variables()) {
          // add v to the interval domain
          idom.set(v, operator[](v));
          // forget v from the constant domain
          // operator-=(v);
        }
      }

      // add csts to the interval domain
      idom += csts;

      // translate back form the interval domain to the constant domain
      if (idom.is_bottom()) {
        set_to_bottom();
      } else {
        for (auto kv : idom) {
          set(kv.first, kv.second);
        }
      }
    }
  }

  interval_t operator[](const linear_expression_t &expr) const {
    if (is_bottom())
      return interval_t::bottom();
    interval_t r(expr.constant());
    for (auto kv : expr) {
      interval_t c(kv.first);
      r += c * to_interval(m_env[kv.second]);
    }
    return r;
  }

  constant_t to_constant(const interval_t &i) const {
    if (i.is_bottom()) {
      return constant_t::bottom();
    } else if (boost::optional<number_t> c = i.singleton()) {
      return constant_t(*c);
    } else {
      return constant_t::top();
    }
  }

  interval_t to_interval(const constant_t &c) const {
    if (c.is_bottom()) {
      return interval_t::bottom();
    } else if (c.is_top()) {
      return interval_t::top();
    } else {
      assert(c.is_constant());
      return interval_t(c.get_constant());
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

  constant_t get_constant(const variable_t &v) const {
    return m_env[v];
  }

  void set_constant(const variable_t &v, constant_t c) {
    m_env.set(v, c);
  }
  
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

  constant_domain_t widening_thresholds(
      const constant_domain_t &o,
      const crab::iterators::thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return m_env.widening_thresholds(o.m_env, ts);
  }

  constant_domain_t operator&&(const constant_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (m_env && o.m_env);
  }

  void set(const variable_t &v, interval_t i) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    m_env.set(v, to_constant(i));
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    m_env -= v;
  }

  interval_t operator[](const variable_t &v) override {
    return to_interval(m_env[v]);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    solve_constraints(csts);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    if (boost::optional<variable_t> v = e.get_variable()) {
      m_env.set(x, m_env[(*v)]);
    } else {
      interval_t ie = operator[](e);
      m_env.set(x, to_constant(ie));
    }
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = to_interval(m_env[y]);
    interval_t zi = to_interval(m_env[z]);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case crab::domains::OP_ADDITION:
      xi = yi + zi;
      break;
    case crab::domains::OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case crab::domains::OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case crab::domains::OP_SDIV:
      xi = yi / zi;
      break;
    case crab::domains::OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case crab::domains::OP_SREM:
      xi = yi.SRem(zi);
      break;
    case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    m_env.set(x, to_constant(xi));
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = to_interval(m_env[y]);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case crab::domains::OP_ADDITION:
      xi = yi + zi;
      break;
    case crab::domains::OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case crab::domains::OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case crab::domains::OP_SDIV:
      xi = yi / zi;
      break;
    case crab::domains::OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case crab::domains::OP_SREM:
      xi = yi.SRem(zi);
      break;
    case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    m_env.set(x, to_constant(xi));
  }

  // intrinsics operations
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
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
  void apply(crab::domains::int_conv_operation_t /*op*/, const variable_t &dst,
             const variable_t &src) override {
    // ignore the widths
    assign(dst, src);
  }

  // bitwise operations
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = to_interval(m_env[y]);
    interval_t zi = to_interval(m_env[z]);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case crab::domains::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case crab::domains::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case crab::domains::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case crab::domains::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case crab::domains::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case crab::domains::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    m_env.set(x, to_constant(xi));
  }

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = to_interval(m_env[y]);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case crab::domains::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case crab::domains::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case crab::domains::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case crab::domains::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case crab::domains::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case crab::domains::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    m_env.set(x, to_constant(xi));
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

      set(lhs, this->operator[](e1) | this->operator[](e2));
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

    m_env.project(variables);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    m_env.rename(from, to);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    m_env.set(new_x, m_env[x]);
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
