#pragma once

/**
 *  Classical sign domain.
 *
 *  This is actually the extended sign domain defined in "Tutorial on
 *  Static Inference of Numeric Invariants by Abstract Interpretation"
 *  by A. Mine 2017.
 */

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/sign.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
class sign_domain final : public crab::domains::abstract_domain_api<
                              sign_domain<Number, VariableName>> {
public:
  using sign_domain_t = sign_domain<Number, VariableName>;
  using abstract_domain_t = crab::domains::abstract_domain_api<sign_domain_t>;
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
  using sign_t = sign<Number>;
  using number_t = Number;
  using varname_t = VariableName;

private:
  using separate_domain_t = ikos::separate_domain<variable_t, sign_t>;

private:
  separate_domain_t m_env;

  sign_domain(separate_domain_t &&env) : m_env(std::move(env)) {}

  void solve_constraints(const linear_constraint_system_t &csts) {

    // Solve < and > constraints
    auto solve_strict_inequality = [this](const variable_t &v, sign_t rhs,
                                          bool is_less_than) {
      sign_t lhs = m_env.at(v);
      CRAB_LOG(
          "sign-domain", crab::outs() << v << "=" << lhs; if (is_less_than) {
            crab::outs() << " < ";
          } else { crab::outs() << " > "; } crab::outs() << rhs
                                                         << "\n";);

      if (rhs.equal_zero()) {
        m_env.set(v, m_env.at(v) &
                         (is_less_than ? sign_t::mk_less_than_zero()
                                       : sign_t::mk_greater_than_zero()));
        CRAB_LOG("sign-domain", crab::outs() << "\t"
                                             << "Refined " << v << "="
                                             << m_env.at(v) << "\n";);
        if (is_bottom())
          return;
      } else {
        if (!is_less_than) {
          std::swap(lhs, rhs);
        }
        if ((lhs.greater_than_zero() || lhs.greater_or_equal_than_zero()) &&
            (rhs.less_than_zero() || rhs.less_or_equal_than_zero())) {
          set_to_bottom();
          CRAB_LOG("sign-domain", crab::outs() << "\t"
                                               << "Refined _|_\n";);
          return;
        }
      }
    };

    // Solve <= and >= constraints
    auto solve_inequality = [this](const variable_t &v, sign_t rhs,
                                   bool is_less_equal) {
      sign_t lhs = m_env.at(v);
      CRAB_LOG(
          "sign-domain", crab::outs() << v << "=" << lhs; if (is_less_equal) {
            crab::outs() << " <= ";
          } else { crab::outs() << " >= "; } crab::outs() << rhs
                                                          << "\n";);

      if (rhs.equal_zero()) {
        m_env.set(v, m_env.at(v) &
                         (is_less_equal
                              ? sign_t::mk_less_or_equal_than_zero()
                              : sign_t::mk_greater_or_equal_than_zero()));
        CRAB_LOG("sign-domain", crab::outs() << "\t"
                                             << "Refined " << v << "="
                                             << m_env.at(v) << "\n";);
        if (is_bottom())
          return;
      } else {
        if (!is_less_equal) {
          std::swap(lhs, rhs);
        }
        if (lhs.greater_or_equal_than_zero() && rhs.less_than_zero()) {
          set_to_bottom();
          CRAB_LOG("sign-domain", crab::outs() << "\t"
                                               << "Refined _|_\n";);
          return;
        }
        if (lhs.greater_than_zero() &&
            (rhs.less_than_zero() || rhs.less_or_equal_than_zero())) {
          set_to_bottom();
          CRAB_LOG("sign-domain", crab::outs() << "\t"
                                               << "Refined _|_\n";);
          return;
        }
      }
    };

    if (!is_bottom()) {
      for (auto const &c : csts) {
        if (c.is_tautology()) {
          continue;
        }
        if (c.is_contradiction()) {
          set_to_bottom();
          return;
        }

        std::vector<std::pair<variable_t, sign_t>> less_than, less_equal;
        std::vector<std::pair<variable_t, sign_t>> greater_than, greater_equal;
        std::vector<std::pair<variable_t, sign_t>> equal, not_equal;

        extract_sign_constraints(c, less_than, less_equal, greater_than,
                                 greater_equal, equal, not_equal);

        CRAB_LOG("sign-domain", crab::outs() << "Adding " << c << "\n";);

        for (auto kv : less_than) {
          solve_strict_inequality(kv.first, kv.second, true /*less_than*/);
        }
        for (auto kv : greater_than) {
          solve_strict_inequality(kv.first, kv.second, false /*greater_than*/);
        }

        for (auto kv : less_equal) {
          solve_inequality(kv.first, kv.second, true /*less_equal_than*/);
        }

        for (auto kv : greater_equal) {
          solve_inequality(kv.first, kv.second, false /*greater_equal_than*/);
        }

        for (auto kv : equal) {
          sign_t lhs = m_env.at(kv.first);
          sign_t rhs = kv.second;
          CRAB_LOG("sign-domain", crab::outs() << kv.first << "=" << lhs
                                               << " == " << rhs << "\n";);
          m_env.set(kv.first, m_env.at(kv.first) & rhs);
          CRAB_LOG("sign-domain", crab::outs() << "\t"
                                               << "Refined " << kv.first << "="
                                               << m_env.at(kv.first) << "\n";);
          if (is_bottom())
            return;
        }

        for (auto kv : not_equal) {
          sign_t lhs = m_env.at(kv.first);
          sign_t rhs = kv.second;
          CRAB_LOG("sign-domain", crab::outs() << kv.first << "=" << lhs
                                               << " != " << rhs << "\n";);
          if (rhs.equal_zero()) {
            m_env.set(kv.first,
                      m_env.at(kv.first) & sign_t::mk_not_equal_zero());
            CRAB_LOG("sign-domain", crab::outs()
                                        << "\t"
                                        << "Refined " << kv.first << "="
                                        << m_env.at(kv.first) << "\n";);
            if (is_bottom())
              return;
          } else if (rhs.not_equal_zero()) {
            m_env.set(kv.first, m_env.at(kv.first) & sign_t::mk_equal_zero());
            CRAB_LOG("sign-domain", crab::outs()
                                        << "\t"
                                        << "Refined " << kv.first << "="
                                        << m_env.at(kv.first) << "\n";);
            if (is_bottom())
              return;
          }
        }
      }
    }
  }

  // Extract constraints of the form v OP sign where OP = {<, <=, >, >=, ==, !=
  // }
  void extract_sign_constraints(
      const linear_constraint_t &c,
      std::vector<std::pair<variable_t, sign_t>> &less_than,
      std::vector<std::pair<variable_t, sign_t>> &less_equal,
      std::vector<std::pair<variable_t, sign_t>> &greater_than,
      std::vector<std::pair<variable_t, sign_t>> &greater_equal,
      std::vector<std::pair<variable_t, sign_t>> &equal,
      std::vector<std::pair<variable_t, sign_t>> &not_equal) {
    auto e = c.expression();
    for (auto kv : e) {
      variable_t pivot = kv.second;
      sign_t res = compute_residual(e, pivot) / sign_t(kv.first);
      if (res.is_bottom()) {
        // this shouldn't happen
        continue;
      }
      if (!res.is_top()) {
        if (c.is_strict_inequality()) {
          if (kv.first < number_t(0)) {
            greater_than.push_back({pivot, res});
          } else {
            less_than.push_back({pivot, res});
          }
        } else if (c.is_inequality()) {
          if (kv.first < number_t(0)) {
            greater_equal.push_back({pivot, res});
          } else {
            less_equal.push_back({pivot, res});
          }
        } else if (c.is_equality()) {
          equal.push_back({pivot, res});
        } else if (c.is_disequation()) {
          not_equal.push_back({pivot, res});
        }
      }
    }
  }

  sign_t compute_residual(const linear_expression_t &e, variable_t pivot) {
    sign_t residual(-e.constant());
    for (auto kv : e) {
      const variable_t &v = kv.second;
      if (v.index() != pivot.index()) {
        residual = residual - (sign_t(kv.first) * m_env.at(v));
      }
    }
    return residual;
  }

  sign_t eval_expr(const linear_expression_t &expr) const {
    if (is_bottom())
      return sign_t::bottom();

    sign_t r(expr.constant());
    for (auto kv : expr) {
      sign_t c(kv.first);
      r = r + (c * m_env.at(kv.second));
    }
    return r;
  }

public:
  sign_domain_t make_top() const override {
    return sign_domain_t(separate_domain_t::top());
  }

  sign_domain_t make_bottom() const override {
    return sign_domain_t(separate_domain_t::bottom());
  }

  void set_to_top() override {
    sign_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    sign_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  sign_t get_sign(const variable_t &v) const { return m_env.at(v); }

  void set_sign(const variable_t &v, sign_t s) { m_env.set(v, s); }

  sign_domain() : m_env(separate_domain_t::top()) {}

  sign_domain(const sign_domain_t &e) : m_env(e.m_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  sign_domain(sign_domain_t &&e) : m_env(std::move(e.m_env)) {}

  sign_domain_t &operator=(const sign_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_env = o.m_env;
    }
    return *this;
  }

  sign_domain_t &operator=(sign_domain_t &&o) {
    if (this != &o) {
      m_env = std::move(o.m_env);
    }
    return *this;
  }

  bool is_bottom() const override { return m_env.is_bottom(); }

  bool is_top() const override { return m_env.is_top(); }

  bool operator<=(const sign_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (m_env <= o.m_env);
  }

  void operator|=(const sign_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    CRAB_LOG("sign-domain",
             crab::outs() << "Join " << m_env << " and " << o.m_env << "\n";);
    m_env = m_env | o.m_env;
    CRAB_LOG("sign-domain", crab::outs() << "Res=" << m_env << "\n";);
  }

  sign_domain_t operator|(const sign_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    return (m_env | o.m_env);
  }

  void operator&=(const sign_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    CRAB_LOG("sign-domain",
             crab::outs() << "Meet " << m_env << " and " << o.m_env << "\n";);
    m_env = m_env & o.m_env;
    CRAB_LOG("sign-domain", crab::outs() << "Res=" << m_env << "\n";);
  }
  
  sign_domain_t operator&(const sign_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    return (m_env & o.m_env);
  }

  sign_domain_t operator||(const sign_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return (m_env | o.m_env);
  }

  sign_domain_t
  widening_thresholds(const sign_domain_t &o,
                      const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return (m_env | o.m_env);
  }

  sign_domain_t operator&&(const sign_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (m_env & o.m_env);
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
      return m_env.at(v).to_interval();
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    solve_constraints(csts);
  }

  DEFAULT_ENTAILS(sign_domain_t)
  
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    CRAB_LOG("sign-domain", crab::outs() << x << " := " << e << "\n";);
    if (!is_bottom()) {    
      if (boost::optional<variable_t> v = e.get_variable()) {
	m_env.set(x, m_env.at(*v));
      } else {
	m_env.set(x, eval_expr(e));
      }
    }
    CRAB_LOG("sign-domain", crab::outs() << "RES=" << m_env.at(x) << "\n";);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.weak_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".weak_assign");
    CRAB_LOG("sign-domain", crab::outs() << "weak_assign(" << x << "," << e << ")\n";);
    if (!is_bottom()) {        
      if (boost::optional<variable_t> v = e.get_variable()) {
	m_env.join(x, m_env.at(*v));
      } else {
	m_env.join(x, eval_expr(e));
      }
    }
    CRAB_LOG("sign-domain", crab::outs() << "RES=" << m_env.at(x) << "\n";);
  }
  

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (!is_bottom()) {        
      sign_t yi = m_env.at(y);
      sign_t zi = m_env.at(z);
      sign_t xi = sign_t::bottom();
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
      }
      m_env.set(x, xi);
    }
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (!is_bottom()) {        
      sign_t yi = m_env.at(y);
      sign_t zi(k);
      sign_t xi = sign_t::bottom();
      
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
      }
      m_env.set(x, xi);
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
                          const sign_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  // backward arithmetic operations
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const sign_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");
    // TODO
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const sign_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    // TODO
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const sign_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    // TODO
  }

  // cast operations
  void apply(crab::domains::int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<sign_domain_t>::apply(*this, op, dst, src);
  }

  // bitwise operations
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (!is_bottom()) {        
      sign_t yi = m_env.at(y);
      sign_t zi = m_env.at(z);
      sign_t xi = sign_t::bottom();
      
      switch (op) {
      case crab::domains::OP_AND: 
	xi = yi.And(zi);
	break;
      case crab::domains::OP_OR: 
	xi = yi.Or(zi);
	break;
      case crab::domains::OP_XOR: 
	xi = yi.Xor(zi);
	break;
      case crab::domains::OP_SHL: 
	xi = yi.Shl(zi);
	break;
      case crab::domains::OP_LSHR: 
	xi = yi.LShr(zi);
	break;
      case crab::domains::OP_ASHR: 
	xi = yi.AShr(zi);
	break;
      }
      m_env.set(x, xi);
    }
  }

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    if (!is_bottom()) {            
      sign_t yi = m_env.at(y);
      sign_t zi(k);
      sign_t xi = sign_t::bottom();
      
      switch (op) {
      case crab::domains::OP_AND: 
	xi = yi.And(zi);
	break;
      case crab::domains::OP_OR: 
	xi = yi.Or(zi);
	break;
      case crab::domains::OP_XOR: 
	xi = yi.Xor(zi);
	break;
      case crab::domains::OP_SHL: 
	xi = yi.Shl(zi);
	break;
      case crab::domains::OP_LSHR: 
	xi = yi.LShr(zi);
	break;
      case crab::domains::OP_ASHR: 
	xi = yi.AShr(zi);
	break;
      }
      m_env.set(x, xi);
    }
  }

  virtual void select(const variable_t &lhs, const linear_constraint_t &cond,
                      const linear_expression_t &e1,
                      const linear_expression_t &e2) override {
    crab::CrabStats::count(domain_name() + ".count.select");
    crab::ScopedCrabStats __st__(domain_name() + ".select");

    if (!is_bottom()) {
      sign_domain_t inv1(*this);
      inv1 += cond;
      if (inv1.is_bottom()) {
        assign(lhs, e2);
        return;
      }

      sign_domain_t inv2(*this);
      inv2 += cond.negate();
      if (inv2.is_bottom()) {
        assign(lhs, e1);
        return;
      }

      m_env.set(lhs, eval_expr(e1) | eval_expr(e2));
    }
  }

  /// sign_domain implements only standard abstract operations of
  /// a numerical domain so it is intended to be used as a leaf domain
  /// in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(sign_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(sign_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(sign_domain_t)

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
      const sign_t &s = it->second;
      if (s.equal_zero()) {
        csts += linear_constraint_t(v == number_t(0));
      } else if (s.less_than_zero()) {
        csts += linear_constraint_t(v < number_t(0));
      } else if (s.greater_than_zero()) {
        csts += linear_constraint_t(v > number_t(0));
      } else if (s.less_or_equal_than_zero()) {
        csts += linear_constraint_t(v <= number_t(0));
      } else if (s.greater_or_equal_than_zero()) {
        csts += linear_constraint_t(v >= number_t(0));
      } else {
        // we cannot represent not_equal_zero as a linear constraint
        continue;
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

  std::string domain_name() const override { return "SignDomain"; }

}; // class sign_domain
} // namespace domains
} // namespace crab

namespace crab {
namespace domains {
template <typename Number, typename VariableName>
struct abstract_domain_traits<sign_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // namespace domains
} // namespace crab
