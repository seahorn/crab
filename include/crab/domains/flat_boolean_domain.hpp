#pragma once

/**
 *  A simple flat 3-valued boolean domain and a reduced product of
 *  this flat bool domain with an arbitrary domain.
 **/

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
// A simple flat 3-valued boolean abstract domain
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

  void operator&=(const flat_boolean_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    CRAB_LOG("flat-boolean",
             crab::outs() << "After meet " << *this << " and " << o << "=");
    m_env = m_env & o.m_env;
    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
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
                      const thresholds<number_t> &) const override {

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

  void weak_assign_bool_cst(const variable_t &x,
			    const linear_constraint_t &cst) override {
    m_env -= x;
    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             crab::outs() << "weak_assign(" << x << "," << bx << ")\n");
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

  void weak_assign_bool_var(const variable_t &x, const variable_t &y,
			    bool is_not_y) override {
    crab::CrabStats::count(domain_name() + ".count.weak_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".weak_assign_bool_var");
    m_env.join(x, (is_not_y ? m_env.at(y).Negate() : m_env.at(y)));
    CRAB_LOG("flat-boolean", auto bx = m_env.at(x);
             crab::outs() << "After " << "weak_assign(" << x << ",";
             if (is_not_y) crab::outs() << "not(" << y << "))";
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
    /* TODO(backward)
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

    /* TODO(backward)
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
// REVISIT: this product domain should be specialized to a numerical
// domain so that we don't need to support operations over reference
// variables and constraints.
//
// The reduction happens in three situations:
//    (1) when bvar := linear_constraint from non-boolean domain to boolean.
//    (2) when bvar := reference_constraint from non-boolean domain to boolean.
//    (3) when assume_bool(bvar) from boolean to non-boolean domain.
// 
// The step (3) is quite weak. Very importantly, we don't trigger any
// reduction from the non-boolean domain to the boolean domain with
// assume(linear_constraint) or assume(reference_constraint).
// Therefore, code like the following won't trigger any propagation so
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
  
  using reduced_domain_product2_t =
      reduced_domain_product2<number_t, varname_t, bool_domain_t, Dom>;
  struct lincst_cmp {
    bool operator()(const linear_constraint_t &c1, const linear_constraint_t &c2) const {
      return c1.lexicographical_compare(c2);
    }
  };
  struct refcst_cmp {
    bool operator()(const reference_constraint_t &c1, const reference_constraint_t &c2) const {
      return c1.lexicographical_compare(c2);
    }
  };
  using lincst_set_t = dual_set_domain<set_domain<linear_constraint_t, lincst_cmp>>;
  using refcst_set_t = dual_set_domain<set_domain<reference_constraint_t, refcst_cmp>>;
  using bool_set_t = dual_set_domain<ikos::discrete_domain<variable_t>>;
  using invariance_domain_t = dual_set_domain<ikos::discrete_domain<variable_t>>;
  using bool_to_lincons_env_t = ikos::separate_domain<variable_t, lincst_set_t>;  
  using bool_to_refcons_env_t = ikos::separate_domain<variable_t, refcst_set_t>;
  using bool_to_bools_env_t = ikos::separate_domain<variable_t, bool_set_t>; 
  
  // Reduced product of flat boolean domain and an arbitrary domain.
  reduced_domain_product2_t m_product;
  /** 
   * These specialized subdomains (m_bool_to_X) are used to remember
   * simple boolean combinations (mostly and's) of constraints in a
   * very local way (i.e., involving few basic blocks within the same
   * function). This allows to trigger some reduction once some
   * Booleans becomes true or false.
   * 
   * REVISIT: a more systematic way would be to build these booleans
   * combinations of constraints as uninterpreted functions and
   * evaluate them later when Booleans are known to be true or false.
   **/
  // Map Boolean variables to set of constraints or other Booleans
  // such that if the variable is evaluated to true then the
  // constraints or the other Booleans are also true.
  bool_to_lincons_env_t m_bool_to_lincsts;
  bool_to_refcons_env_t m_bool_to_refcsts;
  bool_to_bools_env_t   m_bool_to_bools;
  /**
   * m_unchanged_vars keeps track of whether a variable has been
   * modified since the variable was used in a constraint through
   * statements b:= lin_cst or b:= ref_cst.  This allows us to decide
   * which constraints still hold at the time the reduction from
   * boolean variables to non-boolean ones is done. For instance,
   *   
   *  a := x > y;
   *  // unchanged = {x,y}
   *  if (*) {
   *   x := 0;
   *    // unchanged = {y}
   *  } else {
   *    // unchanged = {x,y}
   *  }
   *  // unchanged = {y}
   *  assume(a);
   *  // we cannot say a implies x>y since x might have been modified
   **/
  invariance_domain_t m_unchanged_vars;

  /** Begin helpers to update subdomains **/

  // REVISIT: this function is needed in operations such as forget,
  // project, rename, and expand. Unfortunately, they can destroy all
  // the sharing in separate_domain.
  template<typename BoolEnv>
  void transform_if(BoolEnv &env,
		    std::function<bool(const typename BoolEnv::mapped_type &value)> pred,
		    std::function<void(typename BoolEnv::mapped_type &value)> transform) {
    if (env.is_top() || env.is_bottom()) {
      return;
    }
    
    std::vector<typename BoolEnv::value_type> worklist;
    for (auto it=env.begin(), et=env.end(); it!=et; ++it) {
      if (pred(it->second)) {
	worklist.push_back(*it);
      }
    }
    while (!worklist.empty()) {
       auto pair = worklist.back();
      variable_t v(pair.first);
      typename BoolEnv::mapped_type value(pair.second);
      worklist.pop_back();
      transform(value);
      env -= v;
      env.set(v, value);
    }
  }

  
  template<class BoolToCstEnv>
  void propagate_assign_bool_var(BoolToCstEnv &env,
				 const variable_t &x, const variable_t &y,
				 bool is_negated) {
    if (!is_negated) {
      env.set(x, env.at(y));
    } else {
      typename BoolToCstEnv::mapped_type csts = env.at(y);
      if (csts.size() == 1) {
	// cst is either linear_constraint or reference_constraint
        auto cst = *(csts.begin());
        env.set(x, typename BoolToCstEnv::mapped_type(cst.negate()));
      } else if (csts.size() > 1) { 
	// we do not negate multiple conjunctions because it would
	// become a disjunction so we give up
        env -= x;
      }
    }
  }
  
  // return true if cst's variables haven't been modified. As
  // side-effect, it adds cst into the non-boolean domain if the
  // returned value is true.
  bool add_if_unchanged(const linear_constraint_t &cst) {
    typename invariance_domain_t::set_domain_t vars(
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
    typename invariance_domain_t::set_domain_t vars(
          cst_vars.begin(), cst_vars.end());
    bool unchanged = m_unchanged_vars <= vars;
    if (unchanged) {
      m_product.second().ref_assume(cst);
    }
    return unchanged;
  }

  // Return true if cond evaluates definitely to false
  bool eval_false(const variable_t &cond) const{
    return (m_product.first().get_bool(cond) == boolean_value::get_false());
  }

  // Return true if cond evaluates definitely to true
  bool eval_true(const variable_t &cond) const {
    return (m_product.first().get_bool(cond) == boolean_value::get_true());    
  }
  
  template<class BoolToCstEnv>
  void bwd_reduction_assume_bool(BoolToCstEnv &env,
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
      typename BoolToCstEnv::mapped_type csts = env.at(x);
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        if (add_if_unchanged(cst.negate())) {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Applied " << cst.negate() << "\n";);
        } else {
          CRAB_LOG("flat-boolean",
		   crab::outs() << "\t" << "Cannot apply " << cst.negate() << "\n";);
        }
      }
    }
  }

  /**
   * Given lhs:= select(cond, b1, b2) 
   * if cond is true then
   *     lhs := b1
   * if cond is false then
   *     lhs := b2
   **/
  void fwd_reduction_select_bool(const variable_t &lhs, const variable_t &cond,
				 const variable_t &b1, const variable_t &b2) {
    if (eval_true(cond)) {     
      m_bool_to_lincsts.set(lhs, m_bool_to_lincsts.at(b1));
      m_bool_to_refcsts.set(lhs, m_bool_to_refcsts.at(b1));
      m_bool_to_bools.set(lhs, m_bool_to_bools.at(b1) & bool_set_t(b1));
    } else if (eval_false(cond)) {
      m_bool_to_lincsts.set(lhs, m_bool_to_lincsts.at(b2));
      m_bool_to_refcsts.set(lhs, m_bool_to_refcsts.at(b2));
      m_bool_to_bools.set(lhs, m_bool_to_bools.at(b2) & bool_set_t(b2));
    } else {
      m_bool_to_lincsts -= lhs;
      m_bool_to_refcsts -= lhs;
      m_bool_to_bools -= lhs;
    } 
  }


  /**
   * Given lhs := ite(cond, b1, b2)
   *
   * if b1 is false then
   *    if lhs becomes true then not(cond) and b2 must be true
   *
   * if b2 is false then
   *    if lhs becomes true then cond and b1 must be true
   **/
  template<class BoolToCstEnv>
  void propagate_select_bool(BoolToCstEnv &env,
			     const boolean_value &b1_val, const boolean_value &b2_val,   
			     const variable_t &lhs, const variable_t &cond,
			     const variable_t &b1, const variable_t &b2) {
    // lhs := true false true
    if (b2_val.is_false()) {
      // if lhs becomes true later then cond and b1 must be true
      env.set(lhs, env.at(cond) & env.at(b1));
    } else if (b1_val.is_false()) {
      // if lhs becomes true later then !cond and b2 must be true
      auto cond_csts = env.at(cond);
      if (cond_csts.size() == 1) {
	auto cst = *(cond_csts.begin());
	env.set(lhs, typename BoolToCstEnv::mapped_type(cst.negate()) & env.at(b2));
      } else {
	// we lost the condition because we cannot negate without
	// introducing disjunctions.
	env.set(lhs, env.at(b2));
      }
    }
  }
  
  /** End helpers to update subdomains **/

  /**
   * Reduce from the boolean domain to the non-boolean domain.
   **/
  void reduce_bool_to_csts(const variable_t &x, bool is_negated) {
    if (!is_negated) {
      // We know that x is definitely true

      // backward reduction for m_bool_to_bools
      // if x is true then all Boolean variables in bool_vars must be true.
      auto bool_vars = m_bool_to_bools.at(x);
      for (auto const& v:  bool_vars) {
	m_product.first().assume_bool(v, is_negated);
	m_bool_to_lincsts.set(x, m_bool_to_lincsts.at(x) & m_bool_to_lincsts.at(v));
	m_bool_to_refcsts.set(x, m_bool_to_refcsts.at(x) & m_bool_to_refcsts.at(v));	
      }

      bwd_reduction_assume_bool(m_bool_to_lincsts, x, is_negated);
      bwd_reduction_assume_bool(m_bool_to_refcsts, x, is_negated);
      
    } else {
      /**
       * The kind of facts that we keep track in m_bool_to_lincsts
       * (resp m_bool_to_bools) is of the form "if b becomes true then
       * constraint C (resp. boolean b') MUST be true"
       *
       * Thus, if b is false then nothing can really be said
       * about C (resp b') because if b is false then it does not imply
       * that C or b' are false.
       *
       * TODO: Missing propagation.
       * 
       * However, if we know that b' is false and we have inferred the
       * fact "if b is true then b' is true" then we can conclude that
       * b must be false.
       **/
    }
  }
  

  /**
   * Reduction from the non-boolean domain to the flat boolean domain.
   */
  void reduce_num_cst_to_bool(const variable_t &x,
			      const linear_constraint_t &cst) {
    if (cst.is_tautology()) {
      m_product.first().set_bool(x, boolean_value::get_true());
    } else if (cst.is_contradiction()) {
      m_product.first().set_bool(x, boolean_value::get_false());	
    } else {
      if (m_product.second().entails(cst)) {
	// -- definitely true
	m_product.first().set_bool(x, boolean_value::get_true());
      } else if (m_product.second().entails(cst.negate())) {
	// -- definitely false
	m_product.first().set_bool(x, boolean_value::get_false());	
      } else {
	// -- inconclusive
	m_product.first().set_bool(x, boolean_value::top());
      }
      
      m_bool_to_lincsts.set(x, lincst_set_t(cst));
      // We assume all variables in cst are unchanged unless the
      // opposite is proven
      for (auto const &v : cst.variables()) {
	m_unchanged_vars += v;
      }
    }
    m_bool_to_bools -= x;
  }

  /**
   * Reduction from the non-boolean domain to the flat boolean domain.
   */  
  void reduce_ref_cst_to_bool(const variable_t &x,
			      const reference_constraint_t &cst) {
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
      m_bool_to_refcsts.set(x, refcst_set_t(cst));
      // We assume all variables in cst are unchanged unless the
      // opposite is proven
      for (auto const &v : cst.variables()) {
	m_unchanged_vars += v;
      }
    }
    m_bool_to_bools -= x;

  }

  
  flat_boolean_numerical_domain(reduced_domain_product2_t &&product,
                                bool_to_lincons_env_t &&bool_to_lincsts,
				bool_to_refcons_env_t &&bool_to_refcsts,
				bool_to_bools_env_t &&bool_to_bools,
                                invariance_domain_t &&unchanged_vars)
      : m_product(std::move(product)),
	m_bool_to_lincsts(std::move(bool_to_lincsts)),
	m_bool_to_refcsts(std::move(bool_to_refcsts)),
	m_bool_to_bools(std::move(bool_to_bools)),
        m_unchanged_vars(std::move(unchanged_vars)) {}

public:
  bool_num_domain_t make_top() const override {
    reduced_domain_product2_t prod;
    prod.set_to_top();
    return bool_num_domain_t(std::move(prod), bool_to_lincons_env_t::top(),
			     bool_to_refcons_env_t::top(),
			     bool_to_bools_env_t::top(),
                             invariance_domain_t::top());
  }

  bool_num_domain_t make_bottom() const override {
    reduced_domain_product2_t prod;
    prod.set_to_bottom();
    return bool_num_domain_t(std::move(prod), bool_to_lincons_env_t::bottom(),
			     bool_to_refcons_env_t::bottom(),
			     bool_to_bools_env_t::bottom(),
                             invariance_domain_t::bottom());
  }

  void set_to_top() override {
    reduced_domain_product2_t prod;
    prod.set_to_top();
    bool_num_domain_t abs(std::move(prod), bool_to_lincons_env_t::top(),
			  bool_to_refcons_env_t::top(),
			  bool_to_bools_env_t::top(),
                          invariance_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reduced_domain_product2_t prod;
    prod.set_to_bottom();
    bool_num_domain_t abs(std::move(prod), bool_to_lincons_env_t::bottom(),
			  bool_to_refcons_env_t::bottom(),
			  bool_to_bools_env_t::bottom(),
                          invariance_domain_t::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_numerical_domain()
    : m_product(), m_bool_to_lincsts(), m_bool_to_refcsts(), m_bool_to_bools(),
      m_unchanged_vars() {}

  flat_boolean_numerical_domain(const bool_num_domain_t &other)
      : m_product(other.m_product), m_bool_to_lincsts(other.m_bool_to_lincsts),
	m_bool_to_refcsts(other.m_bool_to_refcsts),
	m_bool_to_bools(other.m_bool_to_bools),
        m_unchanged_vars(other.m_unchanged_vars) {}

  flat_boolean_numerical_domain(const bool_num_domain_t &&other)
      : m_product(std::move(other.m_product)),
        m_bool_to_lincsts(std::move(other.m_bool_to_lincsts)),
        m_bool_to_refcsts(std::move(other.m_bool_to_refcsts)),
	m_bool_to_bools(std::move(other.m_bool_to_bools)),
        m_unchanged_vars(std::move(other.m_unchanged_vars)) {}

  bool_num_domain_t &operator=(const bool_num_domain_t &other) {
    if (this != &other) {
      m_product = other.m_product;
      m_bool_to_lincsts = other.m_bool_to_lincsts;
      m_bool_to_refcsts = other.m_bool_to_refcsts;
      m_bool_to_bools = other.m_bool_to_bools;
      m_unchanged_vars = other.m_unchanged_vars;
    }
    return *this;
  }

  bool_num_domain_t &operator=(const bool_num_domain_t &&other) {
    if (this != &other) {
      m_product = std::move(other.m_product);
      m_bool_to_lincsts = std::move(other.m_bool_to_lincsts);
      m_bool_to_refcsts = std::move(other.m_bool_to_refcsts);
      m_bool_to_bools = std::move(other.m_bool_to_bools);
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
    m_bool_to_lincsts = m_bool_to_lincsts | other.m_bool_to_lincsts;
    m_bool_to_refcsts = m_bool_to_refcsts | other.m_bool_to_refcsts;
    m_bool_to_bools = m_bool_to_bools | other.m_bool_to_bools;
    
    m_unchanged_vars = m_unchanged_vars | other.m_unchanged_vars;
    CRAB_LOG("flat-boolean",
             crab::outs() << "*** After join ***\n"
                          << "\tINV=" << m_product << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_bool_to_lincsts << "\n"
	                  << "\tref constraints for reduction=" << m_bool_to_refcsts << "\n"
	                  << "\tBooleans for reduction=" << m_bool_to_bools << "\n";);
    
  }

  bool_num_domain_t operator|(const bool_num_domain_t &other) const override {
    bool_num_domain_t res(m_product | other.m_product,
			     m_bool_to_lincsts | other.m_bool_to_lincsts,
			     m_bool_to_refcsts | other.m_bool_to_refcsts,
			     m_bool_to_bools | other.m_bool_to_bools,
			     m_unchanged_vars | other.m_unchanged_vars);
    CRAB_LOG("flat-boolean",
             crab::outs() << "*** After join ***\n"
                          << "\tINV=" << res.m_product << "\n"
                          << "\tunchanged vars=" << res.m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << res.m_bool_to_lincsts << "\n"
	                  << "\tref constraints for reduction=" << res.m_bool_to_refcsts << "\n"
	                  << "\tBooleans for reduction=" << res.m_bool_to_bools << "\n";);
    return res;
  }

  bool_num_domain_t operator&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product & other.m_product,
                             m_bool_to_lincsts & other.m_bool_to_lincsts,
                             m_bool_to_refcsts & other.m_bool_to_refcsts,
			     m_bool_to_bools & other.m_bool_to_bools,
                             m_unchanged_vars & other.m_unchanged_vars);
  }

  void operator&=(const bool_num_domain_t &other) override {
    m_product &= other.m_product;
    m_bool_to_lincsts = m_bool_to_lincsts & other.m_bool_to_lincsts;
    m_bool_to_refcsts = m_bool_to_refcsts & other.m_bool_to_refcsts;
    m_bool_to_bools = m_bool_to_bools & other.m_bool_to_bools;
    m_unchanged_vars = m_unchanged_vars & other.m_unchanged_vars;
  }
  
  bool_num_domain_t operator||(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product || other.m_product,
                             m_bool_to_lincsts || other.m_bool_to_lincsts,
                             m_bool_to_refcsts || other.m_bool_to_refcsts,
			     m_bool_to_bools || other.m_bool_to_bools,
                             m_unchanged_vars || other.m_unchanged_vars);
  }

  bool_num_domain_t widening_thresholds(
      const bool_num_domain_t &other,
      const thresholds<number_t> &ts) const override {
    return bool_num_domain_t(m_product.widening_thresholds(other.m_product, ts),
                             m_bool_to_lincsts || other.m_bool_to_lincsts,
                             m_bool_to_refcsts || other.m_bool_to_refcsts,
			     m_bool_to_bools || other.m_bool_to_bools,
                             m_unchanged_vars || other.m_unchanged_vars);
  }

  bool_num_domain_t operator&&(const bool_num_domain_t &other) const override {
    return bool_num_domain_t(m_product && other.m_product,
                             m_bool_to_lincsts && other.m_bool_to_lincsts,
                             m_bool_to_refcsts && other.m_bool_to_refcsts,
			     m_bool_to_bools && other.m_bool_to_bools,
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

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    m_product.weak_assign(x, e);
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

  bool entails(const linear_constraint_t &cst) const override {
    return m_product.second().entails(cst);
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
    if (v.get_type().is_bool()) {
      m_bool_to_lincsts -= v;
      m_bool_to_refcsts -= v;
      m_bool_to_bools -= v;
    } else {
      m_unchanged_vars -= v;
    }

    transform_if(m_bool_to_bools,
		 [&v](const bool_set_t &s) { return s.at(v);},
		 [&v](bool_set_t &s) { s -= v;});
    
    // We should also remove any constraint in
    // m_bool_to_lincsts/m_bool_to_refcsts that involves v.  We don't
    // do it because v is marked as possibly changed so those
    // constraints will not be considered anyway.
    
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
    reduce_num_cst_to_bool(x, cst);

    CRAB_LOG("flat-boolean", auto bx = m_product.first().get_bool(x);
             crab::outs() << "*** Reduction non-boolean --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_bool_to_lincsts << "\n"
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
    reduce_ref_cst_to_bool(x, cst);
    
    CRAB_LOG("flat-boolean", auto bx = m_product.first().get_bool(x);
             crab::outs() << "*** Reduction non-boolean --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
                          << "\treference constraints for reduction=" << m_bool_to_refcsts
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
    propagate_assign_bool_var(m_bool_to_lincsts, x, y, is_negated);
    propagate_assign_bool_var(m_bool_to_refcsts, x, y, is_negated);
    if (!is_negated) {
      m_bool_to_bools.set(x, m_bool_to_bools.at(y) & bool_set_t(y));
    } else {
      // TODO: we don't handle negative booleans in m_bool_to_bools.
      m_bool_to_bools -= x;
    }

    CRAB_LOG("flat-boolean",
             crab::outs() << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_bool_to_lincsts << "\n"
	                  << "\treference constraints for reduction=" << m_bool_to_refcsts << "\n"
	                  << "\tbooleans for reduction=" << m_bool_to_bools
                          << "\n";);
  }

  DEFAULT_WEAK_BOOL_ASSIGN(bool_num_domain_t)
  
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
    //   m_bool_to_lincsts.set
    //     (x, m_bool_to_lincsts [y] & m_bool_to_lincsts [z]);
    //   return;
    // }

    // // we almost lose precision with or and xor except if one of
    // // the operands is false
    // if (op == OP_BOR || op == OP_BXOR) {
    //   if (m_product.first().get_bool(y).is_false()) {
    //     m_bool_to_lincsts.set(x, m_bool_to_lincsts [z]);
    //     return;
    //   }
    //   if (m_product.first().get_bool(z).is_false()) {
    //     m_bool_to_lincsts.set(x, m_bool_to_lincsts [y]);
    //     return;
    //   }
    // }

    m_bool_to_lincsts -= x;
    m_bool_to_refcsts -= x;
    if (op == OP_BAND) {
      m_bool_to_bools.set(x, m_bool_to_bools.at(y) &
			     m_bool_to_bools.at(z) &
			     bool_set_t(y) & bool_set_t(z));
    } else {
      // TODO: we don't handle or/xor
      m_bool_to_bools -= x;
    }
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
	                  << "\tlinear constraints for reduction=" << m_bool_to_lincsts << "\n"
	                  << "\tref constraints for reduction=" << m_bool_to_refcsts << "\n"
	                  << "\tBooleans for reduction=" << m_bool_to_bools << "\n";);
    
    if (is_bottom()) {
      return; 
    }

    reduce_bool_to_csts(x, is_negated);
    
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
	fwd_reduction_select_bool(lhs, cond, b1, b2);
	auto val1 = m_product.first().get_bool(b1);
	auto val2 = m_product.first().get_bool(b2);
	propagate_select_bool(m_bool_to_lincsts, val1, val2, lhs, cond, b1, b2);
	propagate_select_bool(m_bool_to_refcsts, val1, val2, lhs, cond, b1, b2);
	if (val2.is_false()) {
	  m_bool_to_bools.set(lhs,
			      m_bool_to_bools.at(b1) & m_bool_to_bools.at(cond) &
			      bool_set_t(b1) & bool_set_t(cond));
	} else if (val1.is_false()) {
	  m_bool_to_bools.set(lhs, m_bool_to_bools.at(b2) & bool_set_t(b2));
	  // TODO: we don't handle negative booleans in
	  // m_bool_to_bools so we don't add not(cond)
	}
      }
    }

    CRAB_LOG("flat-boolean",
         crab::outs() << lhs << ":= select(" << cond << "," << b1 << "," << b2 << ")\n"
                          << "\tINV=" << m_product << "\n"
                          << "\tunchanged vars=" << m_unchanged_vars << "\n"
	                  << "\tlinear constraints for reduction=" << m_bool_to_lincsts << "\n"
	                  << "\tref constraints for reduction=" << m_bool_to_refcsts << "\n"
	                  << "\tBooleans for reduction=" << m_bool_to_bools << "\n";);	      
  }

  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const bool_num_domain_t &inv) override {
    /** TODO(backward): this can be done better **/
    m_bool_to_lincsts -= lhs;
    m_bool_to_refcsts -= lhs;
    m_bool_to_bools -= lhs;
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const bool_num_domain_t &inv) override {
    /** TODO(backward): this can be done better **/
    operator-=(lhs);
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const bool_num_domain_t &inv) override {
    m_product.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv.m_product);
    /** TODO(backward): this can be done better **/
    m_bool_to_lincsts -= lhs;
    m_bool_to_refcsts -= lhs;
    m_bool_to_bools -= lhs;
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const bool_num_domain_t &inv) override {
    m_product.backward_apply_binary_bool(op, x, y, z, inv.m_product);
    /** TODO(backward): this can be done better **/
    m_bool_to_lincsts -= x;
    m_bool_to_refcsts -= x;
    m_bool_to_bools -= x;
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
    
    for (const variable_t &v : variables) {
      if (v.get_type().is_bool()) {
	m_bool_to_lincsts -= v;
	m_bool_to_refcsts -= v;
	m_bool_to_bools -= v;
      } else {
	m_unchanged_vars -= v;
      }
    }

    transform_if(m_bool_to_bools,
		 [&variables](const bool_set_t &s) {
		   for (auto const& v: variables) {
		     if (s.at(v)) {
		       return true;
		     } 
		   }
		   return false;
		 },
		 [&variables](bool_set_t &s) {
		   for (auto const& v: variables) {
		     s -= v;
		   }
		 });

    // We should also remove any constraint in
    // m_bool_to_lincsts/m_bool_to_refcsts that involves any variable
    // in variables.  We don't do it because they are marked as
    // possibly changed so those constraints will not be considered
    // anyway.
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

    // REVISIT: we conservately throw away all the information in the
    // subdomains. We do it because to be precise here it can be kind
    // of expensive and our assumption is that the scope of these
    // subdomains is just few basic blocks and projection is mostly
    // used during the transformer of a call.
    m_bool_to_lincsts = std::move(bool_to_lincons_env_t::top());
    m_bool_to_refcsts = std::move(bool_to_refcons_env_t::top());
    m_bool_to_bools = std::move(bool_to_bools_env_t::top());
    m_unchanged_vars = std::move(invariance_domain_t::top());
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (is_bottom() || is_top()) {
      return;
    }
    
    assert(from.size() == to.size());
    m_product.rename(from, to);

    std::vector<variable_t> old_bools, new_bools;
    std::copy_if(from.begin(), from.end(), std::back_inserter(old_bools),
		 [](const variable_t &v) {
		   return v.get_type().is_bool();
		 });
    std::copy_if(to.begin(), to.end(), std::back_inserter(new_bools),
		 [](const variable_t &v) {
		   return v.get_type().is_bool();
		 });
    assert(old_bools.size() == new_bools.size());

    // REVISIT: not the most precise renaming but it should be okay.
    m_bool_to_lincsts.rename(old_bools, new_bools);
    m_bool_to_refcsts.rename(old_bools, new_bools);
    m_bool_to_bools = std::move(bool_to_bools_env_t::top());
    // Mark from's variables as possibly modified needed for
    // soundness of m_bool_to_lincsts and m_bool_to_refcsts.
    for (auto const&v: from) {
      if (!v.get_type().is_bool()) {
	m_unchanged_vars -= v;
      }
    }
  }
  
  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    m_product.expand(x, new_x);

    if (x.get_type().is_bool()) {
      m_bool_to_lincsts.set(new_x, m_bool_to_lincsts.at(x));
      m_bool_to_refcsts.set(new_x, m_bool_to_refcsts.at(x));
      // REVISIT: do nothing in m_bool_to_bools is not precise but sound.
    } else {
      if (m_unchanged_vars.at(x)) {
	m_unchanged_vars += new_x;
      }
    }
  }
}; // class flat_boolean_numerical_domain
  
template <typename Num>
struct abstract_domain_traits<flat_boolean_numerical_domain<Num>> {
  using number_t = typename Num::number_t;
  using varname_t = typename Num::varname_t;
};

// template <typename Dom>
// class checker_domain_traits<flat_boolean_numerical_domain<Dom>> {
// public:
//   using this_type = flat_boolean_numerical_domain<Dom>;
//   using linear_constraint_t = typename this_type::linear_constraint_t;
//   using disjunctive_linear_constraint_system_t =
//       typename this_type::disjunctive_linear_constraint_system_t;
//   static bool entail(this_type &lhs,
//                      const disjunctive_linear_constraint_system_t &rhs) {
//     Dom &lhs_dom = lhs.second();
//     return checker_domain_traits<Dom>::entail(lhs_dom, rhs);
//   }
//   static bool entail(const disjunctive_linear_constraint_system_t &lhs,
//                      this_type &rhs) {
//     Dom &rhs_dom = rhs.second();
//     return checker_domain_traits<Dom>::entail(lhs, rhs_dom);
//   }
//   static bool intersect(this_type &inv, const linear_constraint_t &cst) {
//     Dom &dom = inv.second();
//     return checker_domain_traits<Dom>::intersect(dom, cst);
//   }
// };

} // end namespace domains
} // end namespace crab
