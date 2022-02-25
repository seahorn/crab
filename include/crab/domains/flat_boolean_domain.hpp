#pragma once

/*
   A simple flat 3-valued boolean domain and a reduced product of this
   flat bool domain with an arbitrary domain.
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
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using bool_t = boolean_value;
  using separate_domain_t = ikos::separate_domain<variable_t, boolean_value>;
  using iterator = typename separate_domain_t::iterator;

private:
  separate_domain_t m_env;

  flat_boolean_domain(separate_domain_t env) : m_env(env) {}

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

  flat_boolean_domain() : m_env(separate_domain_t::top()) {}

  flat_boolean_domain(const flat_boolean_domain_t &e) : m_env(e.m_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  flat_boolean_domain_t &operator=(const flat_boolean_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o)
      m_env = o.m_env;
    return *this;
  }

  iterator begin() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return m_env.begin();
  }

  iterator end() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return m_env.end();
  }

  bool is_bottom() const override { return m_env.is_bottom(); }

  bool is_top() const override { return m_env.is_top(); }

  bool operator<=(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (m_env <= o.m_env);
  }

  flat_boolean_domain_t
  operator|(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    flat_boolean_domain_t res(m_env | o.m_env);
    CRAB_LOG("flat-boolean", crab::outs() << "After join " << *this << " and "
                                          << o << "=" << res << "\n";);
    return res;
  }

  void operator|=(const flat_boolean_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    CRAB_LOG("flat-boolean",
             crab::outs() << "After join " << *this << " and " << o << "=");
    m_env = m_env | o.m_env;
    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  flat_boolean_domain_t
  operator&(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    flat_boolean_domain_t res(m_env & o.m_env);
    CRAB_LOG("flat-boolean", crab::outs() << "After meet " << *this << " and "
                                          << o << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  operator||(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    flat_boolean_domain_t res(m_env || o.m_env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  widening_thresholds(const flat_boolean_domain_t &o,
                      const iterators::thresholds<number_t> &) const override {

    flat_boolean_domain_t res(m_env || o.m_env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  operator&&(const flat_boolean_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (m_env && o.m_env);
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    if (!is_bottom())
      m_env -= v;
  }

  // flat_boolean_domains_api

  // XXX: the flat boolean domain cannot reason about linear
  // constraints so we assign top to x.
  void assign_bool_cst(const variable_t &x,
                       const linear_constraint_t &cst) override {
    m_env -= x;
    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             crab::outs() << x << ":=" << bx << "\n");
  }

  void assign_bool_ref_cst(const variable_t &x,
                           const reference_constraint_t &cst) override {
    m_env -= x;
    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             crab::outs() << x << ":=" << bx << "\n");
  }

  void assign_bool_var(const variable_t &x, const variable_t &y,
                       bool is_not_y) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");
    m_env.set(x, (is_not_y ? m_env.at(y).Negate() : m_env.at(y)));
    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
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
      m_env.set(x, m_env.at(y).And(m_env.at(z)));
      break;
    case OP_BOR:
      m_env.set(x, m_env.at(y).Or(m_env.at(z)));
      break;
    case OP_BXOR:
      m_env.set(x, m_env.at(y).Xor(m_env.at(z)));
      break;
    default:
      CRAB_ERROR("Unknown boolean operator");
    }

    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             crab::outs() << "After " << x << ":=" << y << " " << op << " " << z
                          << " --->" << x << "=" << bx << "\n");
  }

  void assume_bool(const variable_t &x, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_negated)
      m_env.set(x, m_env.at(x) & boolean_value::get_true());
    else
      m_env.set(x, m_env.at(x) & boolean_value::get_false());

    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             if (!is_negated) crab::outs()
             << "After assume(" << x << ") --> " << x << "=" << bx << "\n";
             else crab::outs() << "After assume(not(" << x << ")) --> " << x
                               << "=" << bx << "\n";);
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      if (b1 == b2) {
	assign_bool_var(lhs, b1, false);
      } else {
	const bool negate = true;
	flat_boolean_domain_t inv1(*this);
	inv1.assume_bool(cond, !negate);
	if (inv1.is_bottom()) {
	  m_env.set(lhs, m_env.at(b2));
	  return;
	}
	
	
	flat_boolean_domain_t inv2(*this);
	inv2.assume_bool(cond, negate);
	if (inv2.is_bottom()) {
	  m_env.set(lhs, m_env.at(b1));
	  return;
	}
	
	m_env.set(lhs, m_env.at(b1) | m_env.at(b2));
      }
    }
  }
  
  // XXX: these methods are not actually part of boolean_operators
  // api but they are used by flat_boolean_numerical_domain and
  // domain_traits.

  void set_bool(const variable_t &x, boolean_value v) { m_env.set(x, v); }

  boolean_value get_bool(const variable_t &x) const { return m_env.at(x); }

  // backward boolean operators
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");
    if (is_bottom())
      return;

    /* nothing to do: flat_boolean_domain ignores this */
    m_env -= lhs;
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");
    if (is_bottom())
      return;

    /* nothing to do: flat_boolean_domain ignores this */
    m_env -= lhs;
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_var");

    if (is_bottom())
      return;
    /* TODO
       assume(lhs == rhs);
       assume(lhs == not(rhs))
    */
    m_env -= lhs;
    *this = *this & inv;
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const flat_boolean_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply_binary_bool");

    if (is_bottom())
      return;

    /* TODO
       if x is true and op=AND then y=true and z=true
       if x is false and op=OR then y=false and z=false
    */
    m_env -= x;
    *this = *this & inv;
  }

  /// flat_boolean_domain implements only boolean operations.  It is
  /// intended to be used as part of a reduced product with a
  /// another domain.
  NUMERICAL_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(flat_boolean_domain_t)
  
  // not part of the numerical_domains api but it should be
  void set(const variable_t &x, interval_t intv) {}

  std::string domain_name() const override { return "Boolean"; }

  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    m_env.write(o);
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
    for (auto kv : m_env) {
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

    m_env.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
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

// Simple reduced product of the flat boolean domain with an arbitrary
// abstract domain.
// 
// The reduction happens in two situations:
//    (1) when bvar := linear_constraint from numerical (non-boolean) to boolean
//    (2) when assume_bool(bvar) from boolean to numerical (non-boolean)
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
// Also, boolean binary operations are mostly ignored during the
// reduction. For instance, the code below won't trigger any reduction
// from booleans to numerical constraints.
//
//   b1 := (x>=0);
//   b2 := (y>=0);
//   b3 := b1 or b2;
//   assume(b3); // we will miss that either x or y is non-negative.
// 
// The reduction code also supports reference constraints:
//     bvar := reference_constraint 
//     ref_assume(bvar)
//
// However, note that if the flat_boolean_numerical_domain is part of
// the base domain used by the region_domain then the following
// abstract operations will never be called:
// 
// - assign_bool_ref_cst: because region domain map reference
//   variables to numerical ones.
// - all region/reference operations 
//
//
template <typename Dom>
class flat_boolean_numerical_domain final
    : public abstract_domain_api<flat_boolean_numerical_domain<Dom>> {
  using bool_num_domain_t = flat_boolean_numerical_domain<Dom>;
  using abstract_domain_t = abstract_domain_api<bool_num_domain_t>;

public:
  using number_t = typename Dom::number_t;
  using varname_t = typename Dom::varname_t;
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
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using bound_t = ikos::bound<number_t>;

private:
  

  // invariant_domain is an abstract domain used to perform a
  // must-forward dataflow analysis that keeps track whether a
  // variable has been modified from any path since the variable was
  // defined up to the current location.
  // 
  // We need to keep track which constraints still hold at the time
  // the reduction from boolean variables to numerical ones is
  // done. For instance,
  //   
  //  a := x > y;
  //  // unchanged = {x,y}
  //  if (*) {
  //    x := 0;
  //    // unchanged = {y}
  //  } else {
  //    // unchanged = {x,y}
  //  }
  //  // unchanged = {y}
  //  assume(a);
  //  // we cannot say a implies x>y since x might have been modified
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

    bool at(const variable_t &v) const{
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


  using reduced_domain_product2_t =
      reduced_domain_product2<number_t, varname_t, bool_domain_t, Dom>;
  
  // bool_to_lincons_env_t and var_refcons_map_t are abstract domains
  // used to map bool variables to sets of constraints such that if
  // the bool variable is true then the conjunction of the constraints
  // must be satisfiable.
  struct linear_constraint_compare {
    bool operator()(const linear_constraint_t &c1, const linear_constraint_t &c2) const {
      return c1.lexicographical_compare(c2);
    }
  };
  struct reference_constraint_compare {
    bool operator()(const reference_constraint_t &c1, const reference_constraint_t &c2) const {
      return c1.lexicographical_compare(c2);
    }
  };  
  using lincst_domain_t = dual_set_domain<linear_constraint_t, linear_constraint_compare>;
  using bool_to_lincons_env_t = ikos::separate_domain<variable_t, lincst_domain_t>;
  using refcst_domain_t = dual_set_domain<reference_constraint_t, reference_constraint_compare>;
  using bool_to_refcons_env_t = ikos::separate_domain<variable_t, refcst_domain_t>;

  // Reduced product of flat boolean domain and an arbitrary domain.
  reduced_domain_product2_t m_product;
  /** 
   * These three domains used to perform reduction from the boolean
   * domain to the other one.
   **/
  // Map from boolean variables to set of constraints such that if the
  // variable is evaluated to true then all the constraints hold.
  bool_to_lincons_env_t m_var_to_lincsts;
  bool_to_refcons_env_t m_var_to_refcsts;
  // Which variable has been unchanged since its last update
  invariance_domain m_unchanged_vars;

  /** Begin helpers to update m_var_to_lincsts and m_var_to_refcsts **/
  template<class BoolToCstEnv>
  void reduction_assign_bool_var(BoolToCstEnv &env,
				 const variable_t &x, const variable_t &y,
				 bool is_negated) {
    if (!is_negated) {
      env.set(x, env.at(y));
    } else {
      typename BoolToCstEnv::value_type csts = env.at(y);
      if (csts.size() == 1) {
	// cst is either linear_constraint or reference_constraint
        auto cst = *(csts.begin());
        env.set(x, typename BoolToCstEnv::value_type(cst.negate()));
        return;
      }
      // we do not negate multiple conjunctions because it would
      // become a disjunction so we give up
      if (csts.size() > 1) {
        env -= x;
      }
    }
  }

  // return true if cst's variables haven't been modified. As
  // side-effect, it adds cst into the non-boolean domain if the
  // returned value is true.
  bool add_if_unchanged(const linear_constraint_t &cst) {
    typename invariance_domain::var_domain_t vars(
          cst.expression().variables_begin(), cst.expression().variables_end());
    bool unchanged = m_unchanged_vars <= vars;
    if (unchanged) {
      m_product.second() += cst;
    }
    return unchanged;
  }

  // return true if cst's variables haven't been modified.  As
  // side-effect, it adds cst into the non-boolean domain if the
  // returned value is true.
  bool add_if_unchanged(const reference_constraint_t &cst) {
    auto cst_vars = cst.variables();
    typename invariance_domain::var_domain_t vars(
          cst_vars.begin(), cst_vars.end());
    bool unchanged = m_unchanged_vars <= vars;
    if (unchanged) {
      m_product.second().ref_assume(cst);
    }
    return unchanged;
  }
  
  template<class BoolToCstEnv>
  void reduction_assume_bool(BoolToCstEnv &env,
			     const variable_t &x, bool is_negated) {

    if (env.at(x).is_top() || env.at(x).is_bottom()) {
      return;
    }

    // Perform reduction from the flat boolean domain to the
    // other domain.
    if (!is_negated) {
      for (auto cst : env.at(x)) {
        // -- we only apply reduction if we know that all the
        // constraint variables have not been modified since they
        // were added into env.
        if (add_if_unchanged(cst)) {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Cannot applied " << cst << "\n";);
        }
      }
    } else {
      // we only perform reduction if there is only one constraint
      typename BoolToCstEnv::value_type csts = env.at(x);
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        if (add_if_unchanged(cst.negate())) {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Cannot apply " << cst << "\n";);
        }
      }
    }
  }

  // Return true if cond evaluates definitely to false
  bool eval_false(const variable_t &cond) const{
    bool_domain_t inv(m_product.first());
    inv.assume_bool(cond, false /*not negated*/);
    return inv.is_bottom();
  }

  // Return true if cond evaluates definitely to true
  bool eval_true(const variable_t &cond) const {
    bool_domain_t inv(m_product.first());
    inv.assume_bool(cond, true /*negated*/);
    return inv.is_bottom();
  }
  
  template<class BoolToCstEnv>
  void reduction_select_bool(BoolToCstEnv &env,
			     const variable_t &lhs, const variable_t &cond,
			     const variable_t &b1, const variable_t &b2) {

    if (eval_true(cond)) {
      env.set(lhs, env.at(b1)); // lhs:= select(cond, b1, b2) --> lhs := b1
    } else if (eval_false(cond)) {
      env.set(lhs, env.at(b2)); // lhs:= select(cond, b1, b2) --> lhs := b2
    } else {
      env -= lhs;
    } 
  }
  
  /** End helpers to update m_var_to_lincsts and m_var_to_refcsts **/


  flat_boolean_numerical_domain(reduced_domain_product2_t &&product,
                                bool_to_lincons_env_t &&var_to_lincsts,
				bool_to_refcons_env_t &&var_to_refcsts,
                                invariance_domain &&unchanged_vars)
      : m_product(std::move(product)),
	m_var_to_lincsts(std::move(var_to_lincsts)),
	m_var_to_refcsts(std::move(var_to_refcsts)),
        m_unchanged_vars(std::move(unchanged_vars)) {}

public:
  bool_num_domain_t make_top() const override {
    reduced_domain_product2_t prod;
    prod.set_to_top();
    return bool_num_domain_t(std::move(prod), bool_to_lincons_env_t::top(),
			     bool_to_refcons_env_t::top(),
                             invariance_domain::top());
  }

  bool_num_domain_t make_bottom() const override {
    reduced_domain_product2_t prod;
    prod.set_to_bottom();
    return bool_num_domain_t(std::move(prod), bool_to_lincons_env_t::bottom(),
			     bool_to_refcons_env_t::bottom(),
                             invariance_domain::bottom());
  }

  void set_to_top() override {
    reduced_domain_product2_t prod;
    prod.set_to_top();
    bool_num_domain_t abs(std::move(prod), bool_to_lincons_env_t::top(),
			  bool_to_refcons_env_t::top(),
                          invariance_domain::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reduced_domain_product2_t prod;
    prod.set_to_bottom();
    bool_num_domain_t abs(std::move(prod), bool_to_lincons_env_t::bottom(),
			  bool_to_refcons_env_t::bottom(),
                          invariance_domain::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_numerical_domain()
    : m_product(), m_var_to_lincsts(), m_var_to_refcsts(), m_unchanged_vars() {}

  flat_boolean_numerical_domain(const bool_num_domain_t &other)
      : m_product(other.m_product), m_var_to_lincsts(other.m_var_to_lincsts),
	m_var_to_refcsts(other.m_var_to_refcsts),
        m_unchanged_vars(other.m_unchanged_vars) {}

  flat_boolean_numerical_domain(const bool_num_domain_t &&other)
      : m_product(std::move(other.m_product)),
        m_var_to_lincsts(std::move(other.m_var_to_lincsts)),
        m_var_to_refcsts(std::move(other.m_var_to_refcsts)),	
        m_unchanged_vars(std::move(other.m_unchanged_vars)) {}

  bool_num_domain_t &operator=(const bool_num_domain_t &other) {
    if (this != &other) {
      m_product = other.m_product;
      m_var_to_lincsts = other.m_var_to_lincsts;
      m_var_to_refcsts = other.m_var_to_refcsts;      
      m_unchanged_vars = other.m_unchanged_vars;
    }
    return *this;
  }

  bool_num_domain_t &operator=(const bool_num_domain_t &&other) {
    if (this != &other) {
      m_product = std::move(other.m_product);
      m_var_to_lincsts = std::move(other.m_var_to_lincsts);
      m_var_to_refcsts = std::move(other.m_var_to_refcsts);      
      m_unchanged_vars = std::move(other.m_unchanged_vars);
    }
    return *this;
  }

  bool is_bottom() const override { return m_product.is_bottom(); }

  bool is_top() const override { return m_product.is_top(); }

  bool_domain_t &first() { return m_product.first(); }

  Dom &second() { return m_product.second(); }

  bool operator<=(const bool_num_domain_t &other) const override {
    return m_product <= other.m_product;
  }

  bool operator==(const bool_num_domain_t &other) const {
    return m_product == other.m_product;
  }

  void operator|=(const bool_num_domain_t &other) override {
    m_product |= other.m_product;
    m_var_to_lincsts = m_var_to_lincsts | other.m_var_to_lincsts;
    m_var_to_refcsts = m_var_to_refcsts | other.m_var_to_refcsts;    
    m_unchanged_vars = m_unchanged_vars | other.m_unchanged_vars;
  }

  bool_num_domain_t operator|(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product | other.m_product,
                             m_var_to_lincsts | other.m_var_to_lincsts,
                             m_var_to_refcsts | other.m_var_to_refcsts,
                             m_unchanged_vars | other.m_unchanged_vars);
  }

  bool_num_domain_t operator&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product & other.m_product,
                             m_var_to_lincsts & other.m_var_to_lincsts,
                             m_var_to_refcsts & other.m_var_to_refcsts,			     
                             m_unchanged_vars & other.m_unchanged_vars);
  }

  bool_num_domain_t operator||(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product || other.m_product,
                             m_var_to_lincsts || other.m_var_to_lincsts,
                             m_var_to_refcsts || other.m_var_to_refcsts,			     
                             m_unchanged_vars || other.m_unchanged_vars);
  }

  bool_num_domain_t widening_thresholds(
      const bool_num_domain_t &other,
      const iterators::thresholds<number_t> &ts) const override {
    return bool_num_domain_t(m_product.widening_thresholds(other.m_product, ts),
                             m_var_to_lincsts || other.m_var_to_lincsts,
                             m_var_to_refcsts || other.m_var_to_refcsts,			     
                             m_unchanged_vars || other.m_unchanged_vars);
  }

  bool_num_domain_t operator&&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product && other.m_product,
                             m_var_to_lincsts && other.m_var_to_lincsts,
                             m_var_to_refcsts && other.m_var_to_refcsts,			     
                             m_unchanged_vars && other.m_unchanged_vars);
  }

  // numerical_domains_api

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    m_unchanged_vars -= x;
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    m_unchanged_vars -= x;
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.assign(x, e);
    m_unchanged_vars -= x;
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    m_product.select(lhs, cond, e1, e2);
    m_unchanged_vars -= lhs;
    
  }  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const bool_num_domain_t &invariant) override {
    m_product.backward_assign(x, e, invariant.m_product);
    m_unchanged_vars -= x;
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const bool_num_domain_t &invariant) override {
    m_product.backward_apply(op, x, y, z, invariant.m_product);
    m_unchanged_vars -= x;
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const bool_num_domain_t &invariant) override {
    m_product.backward_apply(op, x, y, z, invariant.m_product);
    m_unchanged_vars -= x;
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraint");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraint");

    if (csts.is_true()) {
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
      m_product.second() += csts;
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
                m_product.first().assume_bool(var, true /*is_negated*/);
              } else if (k == number_t(1)) {
                m_product.first().assume_bool(var, false /*is_negated*/);
              }
            }
          }
        } else {
          // numerical component
          non_boolean_csts += cst;
        }
      }
      m_product.second() += non_boolean_csts;
    }

#if 0
    // update unchanged_vars
    for (auto const& cst: csts) {
      for (auto const&v : cst.variables()) {
      m_unchanged_vars -= v;
      }
    }
#endif
  }

  void set(const variable_t &x, interval_t intv) {
    // reduced_domain_product2 does not define set method
    m_product.second().set(x, intv); // only on the numerical domain
    m_unchanged_vars -= x;
  }

  interval_t operator[](const variable_t &v) override {
    // reduced_domain_product2 does not define [] method
    boolean_value bv = m_product.first().get_bool(v);
    interval_t isecond = m_product.second()[v];
    if (bv.is_bottom() || isecond.is_bottom()) {
      return interval_t::bottom();
    }     
    if (bv.is_true()) {
      return interval_t(number_t(1)) & isecond;
    } else if (bv.is_false()) {
      return interval_t(number_t(0)) & isecond;
    } else {
      return isecond;
    }
  }

  interval_t at(const variable_t &v) const override {
    // reduced_domain_product2 does not define [] method
    boolean_value bv = m_product.first().get_bool(v);
    interval_t isecond = m_product.second().at(v);
    if (bv.is_bottom() || isecond.is_bottom()) {
      return interval_t::bottom();
    }
    if (bv.is_true()) {
      return interval_t(number_t(1)) & isecond;
    } else if (bv.is_false()) {
      return interval_t(number_t(0)) & isecond;
    } else {
      return isecond;
    }
  }

  void operator-=(const variable_t &v) override {
    m_product -= v;
    m_var_to_lincsts -= v;
    m_var_to_refcsts -= v;    
    m_unchanged_vars -= v;
  }

  // boolean_operators

  void assign_bool_cst(const variable_t &x,
                       const linear_constraint_t &cst) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (m_product.is_bottom()) {
      return;
    }

    m_product.assign_bool_cst(x, cst);

    /** Reduction from the non-boolean domain to the flat boolean domain **/
    if (cst.is_tautology()) {
      m_product.first().set_bool(x, boolean_value::get_true());
    } else if (cst.is_contradiction()) {
      m_product.first().set_bool(x, boolean_value::get_false());	
    } else {
      Dom inv1(m_product.second());
      inv1 += cst;
      if (inv1.is_bottom()) {
	// -- definitely false
	m_product.first().set_bool(x, boolean_value::get_false());
      } else {
	// The call to entail is equivalent to:
	//   Dom inv2(m_product.second());
	//   inv2 += cst.negate();
	//   if (inv2.is_bottom()) { ...}
	// 
	// However, entail splits equalities into inequalities to
	// avoid disequalities after negation.
	if (crab::domains::checker_domain_traits<Dom>::
	    entail(m_product.second(), cst)) {
	  // -- definitely true
	  m_product.first().set_bool(x, boolean_value::get_true());
	} else {
	  // -- inconclusive
	  m_product.first().set_bool(x, boolean_value::top());
	}
      }
      m_var_to_lincsts.set(x, lincst_domain_t(cst));
      // We assume all variables in cst are unchanged unless the
      // opposite is proven
      for (auto const &v : cst.variables()) {
	m_unchanged_vars += v;
      }
    }

    CRAB_LOG("flat-boolean", auto bx = m_product.first().get_bool(x);
             crab::outs() << "*** Reduction non-boolean --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_var_to_lincsts << "\n"
                          << "\tINV=" << m_product << "\n";);
  }

  void assign_bool_ref_cst(const variable_t &x,
                           const reference_constraint_t &cst) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_ref_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_ref_cst");

    if (is_bottom()) {
      return;
    }

    m_product.assign_bool_ref_cst(x, cst);

    /** Reduction from the non-boolean domain to the flat boolean domain **/
    if (cst.is_tautology()) {
      m_product.first().set_bool(x, boolean_value::get_true());
    } else if (cst.is_contradiction()) {
      m_product.first().set_bool(x, boolean_value::get_false());	
    } else {
      Dom inv1(m_product.second());
      inv1.ref_assume(cst);
      if (inv1.is_bottom()) {
	// -- definitely false
	m_product.first().set_bool(x, boolean_value::get_false());
    } else {
	Dom inv2(m_product.second());
	inv2.ref_assume(cst.negate());
	if (inv2.is_bottom()) {
	  // -- definitely true
	  m_product.first().set_bool(x, boolean_value::get_true());
	} else {
	  // -- inconclusive
	  m_product.first().set_bool(x, boolean_value::top());
	}
      }
      m_var_to_refcsts.set(x, refcst_domain_t(cst));
      // We assume all variables in cst are unchanged unless the
      // opposite is proven
      for (auto const &v : cst.variables()) {
	m_unchanged_vars += v;
      }
    }

    CRAB_LOG("flat-boolean", auto bx = m_product.first().get_bool(x);
             crab::outs() << "*** Reduction non-boolean --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
                          << "\treference constraints for reduction=" << m_var_to_refcsts
                          << "\n";);
  }

				 
  void assign_bool_var(const variable_t &x, const variable_t &y,
                       bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (is_bottom()) {
      return;
    }

    m_product.assign_bool_var(x, y, is_negated);
    reduction_assign_bool_var(m_var_to_lincsts, x, y, is_negated);
    reduction_assign_bool_var(m_var_to_refcsts, x, y, is_negated);

    CRAB_LOG("flat-boolean",
             crab::outs() << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_var_to_lincsts << "\n"
                          << "\treference constraints for reduction=" << m_var_to_refcsts	     
                          << "\n";);
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (is_bottom()) {
      return;
    }

    m_product.apply_binary_bool(op, x, y, z);

    // // --- for reduction from boolean to the numerical domain
    // if (op == OP_BAND) {
    //   m_var_to_lincsts.set
    //     (x, m_var_to_lincsts [y] & m_var_to_lincsts [z]);
    //   return;
    // }

    // // we almost lose precision with or and xor except if one of
    // // the operands is false
    // if (op == OP_BOR || op == OP_BXOR) {
    //   if (m_product.first().get_bool(y).is_false()) {
    //     m_var_to_lincsts.set(x, m_var_to_lincsts [z]);
    //     return;
    //   }
    //   if (m_product.first().get_bool(z).is_false()) {
    //     m_var_to_lincsts.set(x, m_var_to_lincsts [y]);
    //     return;
    //   }
    // }

    /// otherwise we give up
    m_var_to_lincsts -= x;
    m_var_to_refcsts -= x;    
  }

  void assume_bool(const variable_t &x, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (is_bottom()) {
      return;
    }
    
    m_product.assume_bool(x, is_negated);

    CRAB_LOG("flat-boolean",
             crab::outs() << "*** Reduction boolean --> non-boolean\n"
                          << "\tassume" << (is_negated ? "(not " : "(") << x
                          << ")\n"
                          << "\tINV=" << m_product << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_var_to_lincsts << "\n"
                          << "\tref constraints for reduction=" << m_var_to_refcsts	     
                          << "\n";);

    if (is_bottom()) {
      return; 
    }
    
    reduction_assume_bool(m_var_to_lincsts, x, is_negated);
    reduction_assume_bool(m_var_to_refcsts, x, is_negated);    

    CRAB_LOG("flat-boolean",
             crab::outs() << "\tAfter reduction=" << m_product << "\n";);
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {


    
    if (!is_bottom()) {
      if (b1 == b2) {
	assign_bool_var(lhs, b1, false);
      } else {
	m_product.select_bool(lhs, cond, b1, b2);
	
	reduction_select_bool(m_var_to_lincsts, lhs, cond, b1, b2);
	reduction_select_bool(m_var_to_refcsts, lhs, cond, b1, b2);
      }
    }
  }

  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const bool_num_domain_t &inv) override {
    /* TODO
       if lhs is true than assume(rhs)
       if lhs is false then assume(not rhs)
    */
    /** TODO: this can be done better **/
    m_var_to_lincsts -= lhs;
    m_var_to_refcsts -= lhs;    
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
    m_product.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv.m_product);
    /** TODO: this can be done better **/
    m_var_to_lincsts -= lhs;
    m_var_to_refcsts -= lhs;        
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const bool_num_domain_t &inv) override {
    m_product.backward_apply_binary_bool(op, x, y, z, inv.m_product);
    /** TODO: this can be done better **/
    m_var_to_lincsts -= x;
    m_var_to_refcsts -= x;        
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
      interval_t i_src = m_product.second()[src];
      interval_t zero = interval_t(number_t(0));
      if (i_src == zero) {
        m_product.first().set_bool(dst, boolean_value::get_false());
      } else if (!(zero <= i_src)) {
        m_product.first().set_bool(dst, boolean_value::get_true());
      } else {
        m_product.first().set_bool(dst, boolean_value::top());
      }
    } else if ((op == OP_ZEXT || op == OP_SEXT) &&
               (get_bitwidth(src) == 1 && get_bitwidth(dst) > 1)) {
      // -- bool to int:
      // if OP_SEXT then true is -1 and false is zero
      // if OP_ZEXT then true is 1 and false is zero
      boolean_value b_src = m_product.first().get_bool(src);
      if (b_src.is_true()) {
        // m_product.second().assign(dst, linear_expression_t(op == OP_SEXT ? -1
        // : 1));
        m_product.second().assign(dst, number_t(1));
      } else if (b_src.is_false()) {
        m_product.second().assign(dst, number_t(0));
      } else {
	m_product.second().apply(op, dst, src);
	
        // The flat boolean domain shouldn't know whether we try to
        // model integers faithfully or not (i.e., obeying
        // machine-arithmetic laws). This code should probably go to
        // domains who model machine arithmetic (e.g., wrapped
        // interval domain).
        //
        // if (op == OP_SEXT) {
        //   m_product.second() += linear_constraint_t(variable_t(dst) >= -1);
        //   m_product.second() += linear_constraint_t(variable_t(dst) <= 0);
        // } else {
        //   m_product.second() += linear_constraint_t(variable_t(dst) >= 0);
        //   m_product.second() += linear_constraint_t(variable_t(dst) <= 1);
        //}
      }
      m_unchanged_vars -= dst;
    } else {
      m_product.apply(op, dst, src);
      m_unchanged_vars -= dst;
    }

    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  // bitwise_operators_api

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_product.apply(op, x, y, z);
    m_unchanged_vars -= x;
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_product.apply(op, x, y, k);
    m_unchanged_vars -= x;
  }

  // array_operators_api

  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    m_product.array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) override {
    m_product.array_load(lhs, a, elem_size, i);
    if (lhs.get_type().is_integer() || lhs.get_type().is_real()) {
      m_unchanged_vars -= lhs;
    }
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    m_product.array_store(a, elem_size, i, val, is_strong_update);
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &v) override {
    m_product.array_store_range(a, elem_size, i, j, v);
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    m_product.array_assign(lhs, rhs);
  }

  // backward array operations

  virtual void
  backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                      const linear_expression_t &lb_idx,
                      const linear_expression_t &ub_idx,
                      const linear_expression_t &val,
                      const bool_num_domain_t &invariant) override {
    m_product.backward_array_init(a, elem_size, lb_idx, ub_idx, val,
                                 invariant.m_product);
  }

  virtual void
  backward_array_load(const variable_t &lhs, const variable_t &a,
                      const linear_expression_t &elem_size,
                      const linear_expression_t &i,
                      const bool_num_domain_t &invariant) override {
    m_product.backward_array_load(lhs, a, elem_size, i, invariant.m_product);
    if (a.get_type().is_integer_array() || a.get_type().is_real_array()) {
      m_unchanged_vars -= lhs;
    }
  }

  virtual void backward_array_store(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &val,
      bool is_strong_update, const bool_num_domain_t &invariant) override {
    m_product.backward_array_store(a, elem_size, i, val, is_strong_update,
                                  invariant.m_product);
  }

  virtual void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const bool_num_domain_t &invariant) override {
    m_product.backward_array_store_range(a, elem_size, i, j, v,
                                        invariant.m_product);
  }

  virtual void
  backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                        const bool_num_domain_t &invariant) override {
    m_product.backward_array_assign(lhs, rhs, invariant.m_product);
  }

  // region/reference api
  void region_init(const variable_t &reg) override {
    m_product.region_init(reg);
  }

  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {
    m_product.region_copy(lhs_reg, rhs_reg);
  }

  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {
    m_product.region_cast(src_reg, dst_reg);
  }
  
  void ref_make(const variable_t &ref, const variable_t &reg,
		const variable_or_constant_t &size,
		const allocation_site &as) override {
    m_product.ref_make(ref, reg, size, as);
  }
  void ref_free(const variable_t &reg, const variable_t &ref)
    override {
    m_product.ref_free(reg, ref);
  }
  
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {
    m_product.ref_load(ref, reg, res);
    if (res.get_type().is_integer() || res.get_type().is_real()) {
      m_unchanged_vars -= res;
    }    
  }
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {
    m_product.ref_store(ref, reg, val);
  }
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {
    m_product.ref_gep(ref1, reg1, ref2, reg2, offset);
  }
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const variable_t &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {
    m_product.ref_load_from_array(lhs, ref, region, index, elem_size);
  }
  void ref_store_to_array(const variable_t &ref, const variable_t &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {
    m_product.ref_store_to_array(ref, region, index, elem_size, val);
  }
  void ref_assume(const reference_constraint_t &cst) override {
    m_product.ref_assume(cst);
  }
  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    m_product.ref_to_int(reg, ref_var, int_var);
    m_unchanged_vars -= int_var;
  }
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    m_product.int_to_ref(int_var, reg, ref_var);
  }
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
		  const variable_t &cond,
		  const variable_or_constant_t &ref1,
		  const boost::optional<variable_t> &rgn1,
		  const variable_or_constant_t &ref2,
		  const boost::optional<variable_t> &rgn2) override {
    m_product.select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
  }
  boolean_value is_null_ref(const variable_t &ref) override {
    return m_product.is_null_ref(ref);
  }
  bool get_allocation_sites(const variable_t &ref,
			    std::vector<allocation_site> &out) override {
    return m_product.get_allocation_sites(ref, out);
  }
  bool get_tags(const variable_t &rgn, const variable_t &ref,
		std::vector<uint64_t> &out) override {
    return m_product.get_tags(rgn, ref, out);
  }
  
  void write(crab_os &o) const override { m_product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() const override {
    linear_constraint_system_t res;
    res += m_product.first().to_linear_constraint_system();
    res += m_product.second().to_linear_constraint_system();
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;
    res += m_product.first().to_disjunctive_linear_constraint_system();
    res += m_product.second().to_disjunctive_linear_constraint_system();
    return res;
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
                          const bool_num_domain_t &invariant) override {
    m_product.backward_intrinsic(name, inputs, outputs, invariant.m_product);
  }
  /* end intrinsics operations */

  void normalize() override { m_product.normalize(); }

  void minimize() override { m_product.minimize(); }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    m_product.forget(variables);
    for (variable_t v : variables) {
      m_var_to_lincsts -= v;
      m_var_to_refcsts -= v;      
      m_unchanged_vars -= v;
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

    m_product.project(variables);

    bool_to_lincons_env_t new_var_to_lincsts;
    bool_to_refcons_env_t new_var_to_refcsts;    
    invariance_domain new_unchanged_vars;
    for (variable_t v : variables) {
      new_var_to_lincsts.set(v, m_var_to_lincsts.at(v));
      new_var_to_refcsts.set(v, m_var_to_refcsts.at(v));      
      
      if (m_unchanged_vars.at(v)) {
        new_unchanged_vars += v;
      }
    }
    std::swap(m_var_to_lincsts, new_var_to_lincsts);
    std::swap(m_var_to_refcsts, new_var_to_refcsts);    
    std::swap(m_unchanged_vars, new_unchanged_vars);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    m_product.expand(x, new_x);
    m_var_to_lincsts.set(new_x, m_var_to_lincsts.at(x));
    m_var_to_refcsts.set(new_x, m_var_to_refcsts.at(x));    
    if (m_unchanged_vars.at(x)) {
      m_unchanged_vars += new_x;
    }
  }
}; // class flat_boolean_numerical_domain

template <typename Num>
struct abstract_domain_traits<flat_boolean_numerical_domain<Num>> {
  using number_t = typename Num::number_t;
  using varname_t = typename Num::varname_t;
};

template <typename Dom>
class checker_domain_traits<flat_boolean_numerical_domain<Dom>> {
public:
  using this_type = flat_boolean_numerical_domain<Dom>;
  using linear_constraint_t = typename this_type::linear_constraint_t;
  using disjunctive_linear_constraint_system_t =
      typename this_type::disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    Dom &lhs_dom = lhs.second();
    return checker_domain_traits<Dom>::entail(lhs_dom, rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    Dom &rhs_dom = rhs.second();
    return checker_domain_traits<Dom>::entail(lhs, rhs_dom);
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    Dom &lhs_dom = lhs.second();
    return checker_domain_traits<Dom>::entail(lhs_dom, rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    Dom &dom = inv.second();
    return checker_domain_traits<Dom>::intersect(dom, cst);
  }
};

} // end namespace domains
} // end namespace crab
