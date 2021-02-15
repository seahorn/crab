#pragma once

/*
   A simple flat 3-valued boolean domain and a reduced product of this
   flat bool domain with an arbitrary numerical domain
*/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {

namespace domains {

// A simple abstract domain for booleans
template <typename Number, typename VariableName>
class flat_boolean_domain final
    : public abstract_domain_api<flat_boolean_domain<Number, VariableName>> {
  using flat_boolean_domain_t = flat_boolean_domain<Number, VariableName>;

public:
  using abstract_domain_t = abstract_domain_api<flat_boolean_domain_t>;
  using number_t = Number;
  using varname_t = VariableName;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using bool_t = boolean_value;
  using separate_domain_t = ikos::separate_domain<variable_t, boolean_value>;
  using iterator = typename separate_domain_t::iterator;

private:
  separate_domain_t _env;

  flat_boolean_domain(separate_domain_t env) : _env(env) {}

public:
  flat_boolean_domain_t make_top() const override {
    return flat_boolean_domain_t(separate_domain_t::top());
  }

  flat_boolean_domain_t make_bottom() const override {
    return flat_boolean_domain_t(separate_domain_t::bottom());
  }

  void set_to_top() override {
    flat_boolean_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    flat_boolean_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_domain() : _env(separate_domain_t::top()) {}

  flat_boolean_domain(const flat_boolean_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  flat_boolean_domain_t &operator=(const flat_boolean_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o)
      _env = o._env;
    return *this;
  }

  iterator begin() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return _env.begin();
  }

  iterator end() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return _env.end();
  }

  bool is_bottom() const override { return _env.is_bottom(); }

  bool is_top() const override { return _env.is_top(); }

  bool operator<=(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (_env <= o._env);
  }

  flat_boolean_domain_t
  operator|(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    flat_boolean_domain_t res(_env | o._env);
    CRAB_LOG("flat-boolean", crab::outs() << "After join " << *this << " and "
                                          << o << "=" << res << "\n";);
    return res;
  }

  void operator|=(const flat_boolean_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    CRAB_LOG("flat-boolean",
             crab::outs() << "After join " << *this << " and " << o << "=");
    _env = _env | o._env;
    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  flat_boolean_domain_t
  operator&(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    flat_boolean_domain_t res(_env & o._env);
    CRAB_LOG("flat-boolean", crab::outs() << "After meet " << *this << " and "
                                          << o << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  operator||(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    flat_boolean_domain_t res(_env || o._env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  widening_thresholds(const flat_boolean_domain_t &o,
                      const iterators::thresholds<number_t> &) const override {

    flat_boolean_domain_t res(_env || o._env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  operator&&(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (_env && o._env);
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    if (!is_bottom())
      _env -= v;
  }

  // flat_boolean_domains_api

  // XXX: the flat boolean domain cannot reason about linear
  // constraints so we assign top to x.
  void assign_bool_cst(const variable_t &x,
                       const linear_constraint_t &cst) override {
    _env -= x;
    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << x << ":=" << bx << "\n");
  }

  void assign_bool_ref_cst(const variable_t &x,
                           const reference_constraint_t &cst) override {
    _env -= x;
    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << x << ":=" << bx << "\n");
  }

  void assign_bool_var(const variable_t &x, const variable_t &y,
                       bool is_not_y) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");
    _env.set(x, (is_not_y ? _env[y].Negate() : _env[y]));
    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << "After " << x << ":=";
             if (is_not_y) crab::outs() << "not(" << y << ")";
             else crab::outs() << y;
             crab::outs() << " --->" << x << "=" << bx << "\n");
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    switch (op) {
    case OP_BAND:
      _env.set(x, _env[y].And(_env[z]));
      break;
    case OP_BOR:
      _env.set(x, _env[y].Or(_env[z]));
      break;
    case OP_BXOR:
      _env.set(x, _env[y].Xor(_env[z]));
      break;
    default:
      CRAB_ERROR("Unknown boolean operator");
    }

    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << "After " << x << ":=" << y << " " << op << " " << z
                          << " --->" << x << "=" << bx << "\n");
  }

  void assume_bool(const variable_t &x, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_negated)
      _env.set(x, _env[x] & boolean_value::get_true());
    else
      _env.set(x, _env[x] & boolean_value::get_false());

    CRAB_LOG("flat-boolean", auto bx = _env[x];
             if (!is_negated) crab::outs()
             << "After assume(" << x << ") --> " << x << "=" << bx << "\n";
             else crab::outs() << "After assume(not(" << x << ")) --> " << x
                               << "=" << bx << "\n";);
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      const bool negate = true;

      flat_boolean_domain_t inv1(*this);
      inv1.assume_bool(cond, !negate);
      if (inv1.is_bottom()) {
	_env.set(lhs, _env[b2]);
	return;
      }
      
      
      flat_boolean_domain_t inv2(*this);
      inv2.assume_bool(cond, negate);
      if (inv2.is_bottom()) {
	_env.set(lhs, _env[b1]);
	return;
      }
      
      _env.set(lhs, _env[b1] | _env[b2]);
    }
  }
  
  // XXX: these methods are not actually part of boolean_operators
  // api but they are used by flat_boolean_numerical_domain and
  // domain_traits.

  void set_bool(const variable_t &x, boolean_value v) { _env.set(x, v); }

  boolean_value get_bool(const variable_t &x) { return _env[x]; }

  // backward boolean operators
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");
    if (is_bottom())
      return;

    /* nothing to do: flat_boolean_domain ignores this */
    _env -= lhs;
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");
    if (is_bottom())
      return;

    /* nothing to do: flat_boolean_domain ignores this */
    _env -= lhs;
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_var");

    if (is_bottom())
      return;
    /** TODO  **/
    /*
       assume(lhs == rhs);
       assume(lhs == not(rhs))
    */
    _env -= lhs;
    *this = *this & inv;
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply_binary_bool");

    if (is_bottom())
      return;

    /** TODO **/
    /*
       if x is true and op=AND then y=true and z=true
       if x is false and op=OR then y=false and z=false
    */
    _env -= x;
    *this = *this & inv;
  }

  /// flat_boolean_domain implements only boolean operations.  It is
  /// intended to be used as part of a reduced product with a
  /// numerical domain.
  NUMERICAL_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  
  // not part of the numerical_domains api but it should be
  void set(const variable_t &x, interval_t intv) {}

  std::string domain_name() const override { return "Boolean"; }

  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    _env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    if (is_bottom())
      return linear_constraint_t::get_false();

    if (is_top())
      return linear_constraint_t::get_true();

    linear_constraint_system_t res;
    for (auto kv : _env) {
      variable_t v(kv.first);
      if (kv.second.is_true()) {
        res += linear_constraint_t(v == number_t(1));
      } else if (kv.second.is_false()) {
        res += linear_constraint_t(v == number_t(0));
      } else {
        res += linear_constraint_t(v >= number_t(0));
        res += linear_constraint_t(v >= number_t(1));
      }
    }
    return res;
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

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    _env.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const flat_boolean_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    for (variable_t v : variables) {
      this->operator-=(v);
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    flat_boolean_domain_t res;
    for (variable_t v : variables) {
      res.set_bool(v, get_bool(v));
    }
    std::swap(*this, res);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set_bool(new_x, get_bool(x));
  }

  void normalize() override {}

  void minimize() override {}

}; // class flat_boolean_domain

template <typename Number, typename VariableName>
struct abstract_domain_traits<flat_boolean_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

// Simple reduced product of the flat boolean domain with an
// arbitrary numerical domain.
// The reduction happens in two situations:
//    (1) when bvar := linear_constraint from numerical to boolean
//    (2) when assume_bool(bvar) from boolean to numerical
// The step (2) is quite weak.
//
// For instance, code like this won't trigger any propagation so
// we won't know the value of x:
//
//   havoc(x);
//   b1 = (x >= 0);
//   y  = zext b1 to i32
//   assume (y >= 1);
//
// Instead, the code should be rewritten like this:
//   havoc(x);
//   b1 = (x >= 0);
//   assume_bool(b1);
//
template <typename NumDom>
class flat_boolean_numerical_domain final
    : public abstract_domain_api<flat_boolean_numerical_domain<NumDom>> {
  using bool_num_domain_t = flat_boolean_numerical_domain<NumDom>;
  using abstract_domain_t = abstract_domain_api<bool_num_domain_t>;

public:
  using number_t = typename NumDom::number_t;
  using varname_t = typename NumDom::varname_t;
  using bool_domain_t = flat_boolean_domain<number_t, varname_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using bound_t = ikos::bound<number_t>;

private:
  // This lattice is the dual of a discrete lattice where
  // elements are linear constraints.
  class lin_cst_set_domain {

    using lincons_domain_t = ikos::discrete_domain<linear_constraint_t>;
    lincons_domain_t m_lincons_set;

  public:
    using iterator = typename lincons_domain_t::iterator;

    lin_cst_set_domain(lincons_domain_t s) : m_lincons_set(s) {}
    lin_cst_set_domain()
        : m_lincons_set(lincons_domain_t::bottom()) /*top by default*/ {}
    lin_cst_set_domain(const lin_cst_set_domain &other)
        : m_lincons_set(other.m_lincons_set) {}

    static lin_cst_set_domain bottom() { return lincons_domain_t::top(); }
    static lin_cst_set_domain top() { return lincons_domain_t::bottom(); }

    bool is_top() const { return m_lincons_set.is_bottom(); }
    bool is_bottom() const { return m_lincons_set.is_top(); }

    bool operator<=(const lin_cst_set_domain &other) const {
      if (other.is_top() || is_bottom())
        return true;
      else
        return other.m_lincons_set <= m_lincons_set;
    }

    bool operator==(const lin_cst_set_domain &other) const {
      return (*this <= other && other <= *this);
    }

    void operator|=(const lin_cst_set_domain &other) {
      m_lincons_set = m_lincons_set & other.m_lincons_set;
    }

    lin_cst_set_domain operator|(const lin_cst_set_domain &other) const {
      return lin_cst_set_domain(m_lincons_set & other.m_lincons_set);
    }

    lin_cst_set_domain operator&(const lin_cst_set_domain &other) const {
      return lin_cst_set_domain(m_lincons_set | other.m_lincons_set);
    }

    lin_cst_set_domain operator||(const lin_cst_set_domain &other) const {
      return this->operator|(other);
    }

    lin_cst_set_domain operator&&(const lin_cst_set_domain &other) const {
      return this->operator&(other);
    }

    lin_cst_set_domain &operator+=(const linear_constraint_t &c) {
      m_lincons_set += c;
      return *this;
    }
    lin_cst_set_domain &operator-=(const linear_constraint_t &c) {
      m_lincons_set -= c;
      return *this;
    }

    std::size_t size() { return m_lincons_set.size(); }
    iterator begin() { return m_lincons_set.begin(); }
    iterator end() { return m_lincons_set.end(); }
    void write(crab::crab_os &o) {
      if (is_bottom())
        o << "_|_";
      else if (is_top())
        o << "top";
      else
        m_lincons_set.write(o);
    }

    friend crab::crab_os &operator<<(crab::crab_os &o,
                                     lin_cst_set_domain &dom) {
      dom.write(o);
      return o;
    }
  };

  using domain_product2_t =
      domain_product2<number_t, varname_t, bool_domain_t, NumDom>;

  // For performing reduction from the boolean domain to the
  // numerical one.
  // Map bool variables to sets of constraints such that if the
  // bool variable is true then the conjunction of the constraints
  // must be satisfiable.
  using var_lincons_map_t =
      ikos::separate_domain<variable_t, lin_cst_set_domain>;
  /** we need to keep track which constraints still hold at the
      time the reduction from boolean variables to numerical ones
      is done. For instance,
      a := x > y;
      // unchanged = {x,y}
      if (*)
        x := 0;
        // unchanged = {y}
      else
        // unchanged = {x,y}

      // unchanged = {y}
      assume(a);
      // we cannot say here that x>y holds since x might have been modified
  **/

  // Simple wrapper for performing a must-forward dataflow
  // analysis that keeps track whether a variable has been
  // modified from any path since the variable was defined up to
  // the current location.
  class invariance_domain {

  public:
    using variable_t = typename linear_constraint_t::variable_t;
    using var_domain_t = ikos::discrete_domain<variable_t>;
    using iterator = typename var_domain_t::iterator;

  private:
    var_domain_t m_varset;

  public:
    invariance_domain() : m_varset(var_domain_t::bottom()) /*top by default*/ {}
    invariance_domain(var_domain_t s) : m_varset(s) {}
    invariance_domain(variable_t v) : m_varset(var_domain_t::bottom()) {
      m_varset += v;
    }
    invariance_domain(const invariance_domain &other)
        : m_varset(other.m_varset) {}

    static invariance_domain bottom() { return var_domain_t::top(); }
    static invariance_domain top() { return var_domain_t::bottom(); }

    bool is_top() const { return m_varset.is_bottom(); }
    bool is_bottom() const { return m_varset.is_top(); }

    /*
             top  = everything might change
          {x,y,z} = nothing changes
              {x} = everything might change except x

                 top
                / | \
             {x} {y} {z}
             / \ / \ / \
           {x,y} {x,z} {y,z}
              \   |   /
               {x,y,z}
                  |
                 bot
     */
    bool operator<=(const invariance_domain &other) const {
      if (other.is_top() || is_bottom())
        return true;
      else
        return other.m_varset <= m_varset;
    }

    void operator|=(const invariance_domain &other) {
      m_varset = m_varset & other.m_varset;
    }

    invariance_domain operator|(const invariance_domain &other) const {
      return invariance_domain(m_varset & other.m_varset);
    }

    invariance_domain operator&(const invariance_domain &other) const {
      return invariance_domain(m_varset | other.m_varset);
    }

    invariance_domain operator||(const invariance_domain &other) const {
      return this->operator|(other);
    }

    invariance_domain operator&&(const invariance_domain &other) const {
      return this->operator&(other);
    }

    invariance_domain &operator+=(const variable_t &v) {
      m_varset += v;
      return *this;
    }
    invariance_domain &operator-=(const variable_t &v) {
      m_varset -= v;
      return *this;
    }

    bool operator[](const variable_t &v) {
      invariance_domain d(v);
      return (d <= *this);
    }
    std::size_t size() { return m_varset.size(); }
    iterator begin() { return m_varset.begin(); }
    iterator end() { return m_varset.end(); }
    void write(crab::crab_os &o) {
      if (is_bottom())
        o << "_|_";
      else if (is_top())
        o << "top";
      else
        m_varset.write(o);
    }
    friend crab::crab_os &operator<<(crab::crab_os &o, invariance_domain &dom) {
      dom.write(o);
      return o;
    }
  };

  domain_product2_t _product;
  var_lincons_map_t _var_to_csts;
  invariance_domain _unchanged_vars;

  flat_boolean_numerical_domain(domain_product2_t &&product,
                                var_lincons_map_t &&var_to_csts,
                                invariance_domain &&unchanged_vars)
      : _product(std::move(product)), _var_to_csts(std::move(var_to_csts)),
        _unchanged_vars(std::move(unchanged_vars)) {}

public:
  bool_num_domain_t make_top() const override {
    domain_product2_t prod;
    prod.set_to_top();
    return bool_num_domain_t(std::move(prod), var_lincons_map_t::top(),
                             invariance_domain::top());
  }

  bool_num_domain_t make_bottom() const override {
    domain_product2_t prod;
    prod.set_to_bottom();
    return bool_num_domain_t(std::move(prod), var_lincons_map_t::bottom(),
                             invariance_domain::bottom());
  }

  void set_to_top() override {
    domain_product2_t prod;
    prod.set_to_top();
    bool_num_domain_t abs(std::move(prod), var_lincons_map_t::top(),
                          invariance_domain::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    domain_product2_t prod;
    prod.set_to_bottom();
    bool_num_domain_t abs(std::move(prod), var_lincons_map_t::bottom(),
                          invariance_domain::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_numerical_domain()
      : _product(), _var_to_csts(), _unchanged_vars() {}

  flat_boolean_numerical_domain(const bool_num_domain_t &other)
      : _product(other._product), _var_to_csts(other._var_to_csts),
        _unchanged_vars(other._unchanged_vars) {}

  flat_boolean_numerical_domain(const bool_num_domain_t &&other)
      : _product(std::move(other._product)),
        _var_to_csts(std::move(other._var_to_csts)),
        _unchanged_vars(std::move(other._unchanged_vars)) {}

  bool_num_domain_t &operator=(const bool_num_domain_t &other) {
    if (this != &other) {
      _product = other._product;
      _var_to_csts = other._var_to_csts;
      _unchanged_vars = other._unchanged_vars;
    }
    return *this;
  }

  bool_num_domain_t &operator=(const bool_num_domain_t &&other) {
    if (this != &other) {
      _product = std::move(other._product);
      _var_to_csts = std::move(other._var_to_csts);
      _unchanged_vars = std::move(other._unchanged_vars);
    }
    return *this;
  }

  bool is_bottom() const override { return _product.is_bottom(); }

  bool is_top() const override { return _product.is_top(); }

  bool_domain_t &first() { return _product.first(); }

  NumDom &second() { return _product.second(); }

  bool operator<=(const bool_num_domain_t &other) const override {
    return _product <= other._product;
  }

  bool operator==(const bool_num_domain_t &other) const {
    return _product == other._product;
  }

  void operator|=(const bool_num_domain_t &other) override {
    _product |= other._product;
    _var_to_csts = _var_to_csts | other._var_to_csts;
    _unchanged_vars = _unchanged_vars | other._unchanged_vars;
  }

  bool_num_domain_t operator|(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(_product | other._product,
                             _var_to_csts | other._var_to_csts,
                             _unchanged_vars | other._unchanged_vars);
  }

  bool_num_domain_t operator&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(_product & other._product,
                             _var_to_csts & other._var_to_csts,
                             _unchanged_vars & other._unchanged_vars);
  }

  bool_num_domain_t operator||(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(_product || other._product,
                             _var_to_csts || other._var_to_csts,
                             _unchanged_vars || other._unchanged_vars);
  }

  bool_num_domain_t widening_thresholds(
      const bool_num_domain_t &other,
      const iterators::thresholds<number_t> &ts) const override {
    return bool_num_domain_t(_product.widening_thresholds(other._product, ts),
                             _var_to_csts || other._var_to_csts,
                             _unchanged_vars || other._unchanged_vars);
  }

  bool_num_domain_t operator&&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(_product && other._product,
                             _var_to_csts && other._var_to_csts,
                             _unchanged_vars && other._unchanged_vars);
  }

  // numerical_domains_api

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    _product.apply(op, x, y, z);
    _unchanged_vars -= variable_t(x);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    _product.apply(op, x, y, k);
    _unchanged_vars -= variable_t(x);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    _product.assign(x, e);
    _unchanged_vars -= variable_t(x);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    _product.select(lhs, cond, e1, e2);
    _unchanged_vars -= variable_t(lhs);
    
  }  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const bool_num_domain_t &invariant) override {
    _product.backward_assign(x, e, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const bool_num_domain_t &invariant) override {
    _product.backward_apply(op, x, y, z, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const bool_num_domain_t &invariant) override {
    _product.backward_apply(op, x, y, z, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraint");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraint");

    if (csts.is_true() || csts.is_false()) {
      return;
    }

    bool all_non_boolean = true;
    for (const linear_constraint_t &cst : csts) {
      if (std::any_of(
              cst.expression().variables_begin(),
              cst.expression().variables_end(),
              [](const variable_t &v) { return v.get_type().is_bool(); })) {
        all_non_boolean = false;
        break;
      }
    }

    if (all_non_boolean) {
      /// -- Common case: all constraints are non-boolean.
      _product.second() += csts;
    } else {
      /// -- We have both boolean and non-boolean constraints.

      // Normalization ensures that inequality pairs x <= k and x >= k
      // are replaced with a single equality x = k, where x always
      // appears as positive.
      linear_constraint_system_t norm_csts = csts.normalize();
      linear_constraint_system_t non_boolean_csts;
      for (const linear_constraint_t &cst : norm_csts) {
        if (cst.is_equality() && std::all_of(cst.expression().variables_begin(),
                                             cst.expression().variables_end(),
                                             [](const variable_t &v) {
                                               return v.get_type().is_bool();
                                             })) {
          // boolean component
          const linear_expression_t &exp = cst.expression();
          if (exp.size() == 1) {
            number_t coef = (*(exp.begin())).first;
            const variable_t &var = (*(exp.begin())).second;
            number_t k = exp.constant();
            if (coef == number_t(1)) {
              if (k == number_t(0)) {
                _product.first().assume_bool(var, true /*is_negated*/);
              } else if (k == number_t(1)) {
                _product.first().assume_bool(var, false /*is_negated*/);
              }
            }
          }
        } else {
          // numerical component
          non_boolean_csts += cst;
        }
      }
      _product.second() += non_boolean_csts;
    }

#if 0
    // update unchanged_vars
    for (auto const& cst: csts) {
      for (auto const&v : cst.variables()) {
      _unchanged_vars -= v;
      }
    }
#endif
  }

  void set(const variable_t &x, interval_t intv) {
    // domain_product2 does not define set method
    _product.second().set(x, intv); // only on the numerical domain
    _unchanged_vars -= variable_t(x);
  }

  interval_t operator[](const variable_t &v) override {
    // domain_product2 does not define [] method
    boolean_value bv = _product.first().get_bool(v);
    interval_t isecond = _product.second()[v];

    if (bv.is_bottom() || isecond.is_bottom())
      return interval_t::bottom();

    if (bv.is_true())
      return interval_t(number_t(1)) & isecond;
    else if (bv.is_false())
      return interval_t(number_t(0)) & isecond;
    else
      return isecond;
  }

  void operator-=(const variable_t &v) override {
    _product -= v;
    _var_to_csts -= v;
    _unchanged_vars -= variable_t(v);
  }

  // boolean_operators

  void assign_bool_cst(const variable_t &x,
                       const linear_constraint_t &cst) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    /// Reduction from the numerical domain to the flat boolean
    /// domain

    if (_product.is_bottom())
      return;

    _product.assign_bool_cst(x, cst);

    NumDom inv1(_product.second());
    inv1 += cst;
    if (inv1.is_bottom()) {
      // -- definitely false
      _product.first().set_bool(x, boolean_value::get_false());
    } else {
      NumDom inv2(_product.second());
      inv2 += cst.negate();
      if (inv2.is_bottom()) {
        // -- definitely true
        _product.first().set_bool(x, boolean_value::get_true());
      } else {
#if 0
	    // XXX: before we give up we convert into intervals and
	    // check again if the negated constraint is bottom.  This
	    // is useful because e.g., apron domains completely ignore
	    // disequations.
	    // 
	    // JN: I disable this code because AFIK all domains reason
	    // now about disequalities, included apron/elina domains.
	    interval_domain<number_t,varname_t> inv3;
	    inv3 += cst.negate();	    
	    for (auto c: _product.second().to_linear_constraint_system())
	    { inv3 += c;}
	    if (inv3.is_bottom()) {
	      // -- definitely true	  
	      _product.first().set_bool(x, boolean_value::get_true());
	    } else {
	      // -- inconclusive
	      _product.first().set_bool(x, boolean_value::top());
	    }
#else
        // -- inconclusive
        _product.first().set_bool(x, boolean_value::top());
#endif
      }
    }
    _var_to_csts.set(x, lin_cst_set_domain(cst));
    // We assume all variables in cst are unchanged unless the
    // opposite is proven
    for (auto const &v : cst.variables()) {
      _unchanged_vars += v;
    }

    CRAB_LOG("flat-boolean", auto bx = _product.first().get_bool(x);
             crab::outs() << "*** Reduction numerical --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);
  }

  void assign_bool_ref_cst(const variable_t &x,
                           const reference_constraint_t &cst) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (is_bottom())
      return;

    operator-=(x);
  }

  void assign_bool_var(const variable_t &x, const variable_t &y,
                       bool is_not_y) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (is_bottom())
      return;

    _product.assign_bool_var(x, y, is_not_y);

    if (!is_not_y)
      _var_to_csts.set(x, _var_to_csts[y]);
    else {
      auto csts = _var_to_csts[y];
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        _var_to_csts.set(x, lin_cst_set_domain(cst.negate()));
        return;
      }
      // we do not negate multiple conjunctions because it would
      // become a disjunction so we give up
      if (csts.size() > 1)
        _var_to_csts -= x;
    }

    CRAB_LOG("flat-boolean",
             crab::outs() << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (is_bottom())
      return;

    _product.apply_binary_bool(op, x, y, z);

    // // --- for reduction from boolean to the numerical domain
    // if (op == OP_BAND) {
    //   _var_to_csts.set
    //     (x, _var_to_csts [y] & _var_to_csts [z]);
    //   return;
    // }

    // // we almost lose precision with or and xor except if one of
    // // the operands is false
    // if (op == OP_BOR || op == OP_BXOR) {
    //   if (_product.first().get_bool(y).is_false()) {
    //     _var_to_csts.set(x, _var_to_csts [z]);
    //     return;
    //   }
    //   if (_product.first().get_bool(z).is_false()) {
    //     _var_to_csts.set(x, _var_to_csts [y]);
    //     return;
    //   }
    // }

    /// otherwise we give up
    _var_to_csts -= x;
  }

  void assume_bool(const variable_t &x, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (is_bottom())
      return;

    // return true if cst's variables haven't been modified
    auto unmodified_constraint = [this](const linear_constraint_t &cst) {
      typename invariance_domain::var_domain_t vars(
          cst.expression().variables_begin(), cst.expression().variables_end());
      return _unchanged_vars <= vars;
    };

    _product.assume_bool(x, is_negated);

    CRAB_LOG("flat-boolean",
             crab::outs() << "*** Reduction boolean --> numerical\n"
                          << "\tassume" << (is_negated ? "(not " : "(") << x
                          << ")\n"
                          << "\tINV=" << _product.second() << "\n"
                          << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);

    if (_var_to_csts[x].is_top() || _var_to_csts[x].is_bottom())
      return;

    // Perform reduction from the flat boolean domain to the
    // numerical domain.
    if (!is_negated) {
      for (auto cst : _var_to_csts[x]) {
        // -- we only apply reduction if we know that all the
        // constraint variables have not been modified since they
        // were added into _var_to_csts.
        if (unmodified_constraint(cst)) {
          _product.second() += cst;
          CRAB_LOG("flat-boolean", crab::outs() << "\t"
                                                << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean", crab::outs()
                                       << "\t"
                                       << "Cannot applied " << cst << "\n";);
        }
      }
    } else {
      // we only perform reduction if there is only one constraint
      auto csts = _var_to_csts[x];
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        if (unmodified_constraint(cst)) {
          _product.second() += cst.negate();
          CRAB_LOG("flat-boolean", crab::outs() << "\t"
                                                << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean", crab::outs()
                                       << "\t"
                                       << "Cannot apply " << cst << "\n";);
        }
      }
    }

    CRAB_LOG("flat-boolean",
             crab::outs() << "After reduction=" << _product.second() << "\n";);
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {

    if (!is_bottom()) {
      _product.select_bool(lhs, cond, b1, b2);
      // Maybe a bit too imprecise
      _var_to_csts -= lhs;
    }
  }

  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const bool_num_domain_t &inv) override {
    /** TODO **/
    /*
       if lhs is true than assume(rhs)
       if lhs is false then assume(not rhs)
    */
    /** TODO: this can be done better **/
    _var_to_csts -= lhs;
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const bool_num_domain_t &inv) override {
    /** TODO: this can be done better **/
    operator-=(lhs);
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const bool_num_domain_t &inv) override {
    _product.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv._product);
    /** TODO: this can be done better **/
    _var_to_csts -= lhs;
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const bool_num_domain_t &inv) override {
    _product.backward_apply_binary_bool(op, x, y, z, inv._product);
    /** TODO: this can be done better **/
    _var_to_csts -= x;
  }

  // cast_operators_api

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    auto get_bitwidth = [](const variable_t v) {
      auto ty = v.get_type();
      if (!(ty.is_integer() || ty.is_bool())) {
        CRAB_ERROR("unexpected types in cast operation");
      }
      return (ty.is_integer() ? ty.get_integer_bitwidth() : 1);
    };
    CRAB_LOG("flat-boolean", crab::outs()
                                 << src << ":" << get_bitwidth(src) << " " << op
                                 << " " << dst << ":" << get_bitwidth(dst)
                                 << " with " << *this << "\n");

    if (op == OP_TRUNC && (get_bitwidth(src) > 1 && get_bitwidth(dst) == 1)) {
      // -- int to bool:
      // assume that zero is false and non-zero is true
      interval_t i_src = _product.second()[src];
      interval_t zero = interval_t(number_t(0));
      if (i_src == zero) {
        _product.first().set_bool(dst, boolean_value::get_false());
      } else if (!(zero <= i_src)) {
        _product.first().set_bool(dst, boolean_value::get_true());
      } else {
        _product.first().set_bool(dst, boolean_value::top());
      }
    } else if ((op == OP_ZEXT || op == OP_SEXT) &&
               (get_bitwidth(src) == 1 && get_bitwidth(dst) > 1)) {
      // -- bool to int:
      // if OP_SEXT then true is -1 and false is zero
      // if OP_ZEXT then true is 1 and false is zero
      boolean_value b_src = _product.first().get_bool(src);
      if (b_src.is_true()) {
        // _product.second().assign(dst, linear_expression_t(op == OP_SEXT ? -1
        // : 1));
        _product.second().assign(dst, number_t(1));
      } else if (b_src.is_false()) {
        _product.second().assign(dst, number_t(0));
      } else {
	_product.second().apply(op, dst, src);
	
        // The flat boolean domain shouldn't know whether we try to
        // model integers faithfully or not (i.e., obeying
        // machine-arithmetic laws). This code should probably go to
        // domains who model machine arithmetic (e.g., wrapped
        // interval domain).
        //
        // if (op == OP_SEXT) {
        //   _product.second() += linear_constraint_t(variable_t(dst) >= -1);
        //   _product.second() += linear_constraint_t(variable_t(dst) <= 0);
        // } else {
        //   _product.second() += linear_constraint_t(variable_t(dst) >= 0);
        //   _product.second() += linear_constraint_t(variable_t(dst) <= 1);
        //}
      }
      _unchanged_vars -= variable_t(dst);
    } else {
      _product.apply(op, dst, src);
      _unchanged_vars -= variable_t(dst);
    }

    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  // bitwise_operators_api

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    _product.apply(op, x, y, z);
    _unchanged_vars -= variable_t(x);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    _product.apply(op, x, y, k);
    _unchanged_vars -= variable_t(x);
  }

  // array_operators_api

  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    _product.array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) override {
    _product.array_load(lhs, a, elem_size, i);
    if (lhs.get_type().is_integer() || lhs.get_type().is_real()) {
      _unchanged_vars -= variable_t(lhs);
    }
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    _product.array_store(a, elem_size, i, val, is_strong_update);
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &v) override {
    _product.array_store_range(a, elem_size, i, j, v);
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    _product.array_assign(lhs, rhs);
  }

  // backward array operations

  virtual void
  backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                      const linear_expression_t &lb_idx,
                      const linear_expression_t &ub_idx,
                      const linear_expression_t &val,
                      const bool_num_domain_t &invariant) override {
    _product.backward_array_init(a, elem_size, lb_idx, ub_idx, val,
                                 invariant._product);
  }

  virtual void
  backward_array_load(const variable_t &lhs, const variable_t &a,
                      const linear_expression_t &elem_size,
                      const linear_expression_t &i,
                      const bool_num_domain_t &invariant) override {
    _product.backward_array_load(lhs, a, elem_size, i, invariant._product);
    if (a.get_type().is_integer_array() || a.get_type().is_real_array()) {
      _unchanged_vars -= variable_t(lhs);
    }
  }

  virtual void backward_array_store(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &val,
      bool is_strong_update, const bool_num_domain_t &invariant) override {
    _product.backward_array_store(a, elem_size, i, val, is_strong_update,
                                  invariant._product);
  }

  virtual void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const bool_num_domain_t &invariant) override {
    _product.backward_array_store_range(a, elem_size, i, j, v,
                                        invariant._product);
  }

  virtual void
  backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                        const bool_num_domain_t &invariant) override {
    _product.backward_array_assign(lhs, rhs, invariant._product);
  }

  // region/reference api
  void region_init(const variable_t &reg) override {
    _product.region_init(reg);
  }

  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {
    _product.region_copy(lhs_reg, rhs_reg);
  }

  void ref_make(const variable_t &ref, const variable_t &reg,
		const allocation_site &as) override {
    _product.ref_make(ref, reg, as);
  }
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {
    _product.ref_load(ref, reg, res);
    if (res.get_type().is_integer() || res.get_type().is_real()) {
      _unchanged_vars -= variable_t(res);
    }    
  }
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {
    _product.ref_store(ref, reg, val);
  }
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {
    _product.ref_gep(ref1, reg1, ref2, reg2, offset);
  }
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const variable_t &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {
    _product.ref_load_from_array(lhs, ref, region, index, elem_size);
  }
  void ref_store_to_array(const variable_t &ref, const variable_t &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {
    _product.ref_store_to_array(ref, region, index, elem_size, val);
  }
  void ref_assume(const reference_constraint_t &cst) override {
    _product.ref_assume(cst);
  }
  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    _product.ref_to_int(reg, ref_var, int_var);
    _unchanged_vars -= variable_t(int_var);
  }
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    _product.int_to_ref(int_var, reg, ref_var);
  }
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
		  const variable_t &cond,
		  const variable_or_constant_t &ref1,
		  const boost::optional<variable_t> &rgn1,
		  const variable_or_constant_t &ref2,
		  const boost::optional<variable_t> &rgn2) override {
    _product.select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
  }
  boolean_value is_null_ref(const variable_t &ref) override {
    return _product.is_null_ref(ref);
  }
  bool get_allocation_sites(const variable_t &ref,
			    std::vector<allocation_site> &out) override {
    return _product.get_allocation_sites(ref, out);
  }
  
  void write(crab_os &o) const override { _product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() const override {
    linear_constraint_system_t res;
    res += _product.first().to_linear_constraint_system();
    res += _product.second().to_linear_constraint_system();
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;
    res += _product.first().to_disjunctive_linear_constraint_system();
    res += _product.second().to_disjunctive_linear_constraint_system();
    return res;
  }

  std::string domain_name() const override { return _product.domain_name(); }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    _product.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    _product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const bool_num_domain_t &invariant) override {
    _product.backward_intrinsic(name, inputs, outputs, invariant._product);
  }
  /* end intrinsics operations */

  void normalize() override { _product.normalize(); }

  void minimize() override { _product.minimize(); }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    _product.forget(variables);
    for (variable_t v : variables) {
      _var_to_csts -= v;
      _unchanged_vars -= v;
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    if (variables.empty()) {
      set_to_top();
      return;
    }

    _product.project(variables);

    var_lincons_map_t new_var_to_csts;
    invariance_domain new_unchanged_vars;
    for (variable_t v : variables) {
      new_var_to_csts.set(v, _var_to_csts[v]);
      if (_unchanged_vars[v]) {
        new_unchanged_vars += v;
      }
    }
    std::swap(_var_to_csts, new_var_to_csts);
    std::swap(_unchanged_vars, new_unchanged_vars);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    _product.expand(x, new_x);
    _var_to_csts.set(new_x, _var_to_csts[x]);
    if (_unchanged_vars[variable_t(x)]) {
      _unchanged_vars += variable_t(new_x);
    }
  }
}; // class flat_boolean_numerical_domain

template <typename Num>
struct abstract_domain_traits<flat_boolean_numerical_domain<Num>> {
  using number_t = typename Num::number_t;
  using varname_t = typename Num::varname_t;
};

template <typename NumDom>
class checker_domain_traits<flat_boolean_numerical_domain<NumDom>> {
public:
  using this_type = flat_boolean_numerical_domain<NumDom>;
  using linear_constraint_t = typename this_type::linear_constraint_t;
  using disjunctive_linear_constraint_system_t =
      typename this_type::disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    NumDom &lhs_dom = lhs.second();
    return checker_domain_traits<NumDom>::entail(lhs_dom, rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    NumDom &rhs_dom = rhs.second();
    return checker_domain_traits<NumDom>::entail(lhs, rhs_dom);
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    NumDom &lhs_dom = lhs.second();
    return checker_domain_traits<NumDom>::entail(lhs_dom, rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    NumDom &dom = inv.second();
    return checker_domain_traits<NumDom>::intersect(dom, cst);
  }
};

} // end namespace domains
} // end namespace crab
