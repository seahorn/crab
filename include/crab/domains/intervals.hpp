/*******************************************************************************
 *
 * Standard domain of intervals.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Contributors: Alexandre C. D. Wimmers (alexandre.c.wimmers@nasa.gov)
 *               Jorge A. Navas (jorge.navas@sri.com)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/stats.hpp>

namespace ikos {

template <typename Number, typename VariableName,
          std::size_t max_reduction_cycles = 10>
class interval_domain final
    : public crab::domains::abstract_domain_api<
          interval_domain<Number, VariableName, max_reduction_cycles>> {
public:
  using interval_domain_t =
      interval_domain<Number, VariableName, max_reduction_cycles>;
  using abstract_domain_t =
      crab::domains::abstract_domain_api<interval_domain_t>;
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
  using separate_domain_t = separate_domain<variable_t, interval_t>;
  using solver_t =
      linear_interval_solver<number_t, varname_t, separate_domain_t>;

public:
  using iterator = typename separate_domain_t::iterator;

private:
  separate_domain_t _env;

  interval_domain(separate_domain_t env) : _env(env) {}

  interval_t operator[](linear_expression_t expr) {
    interval_t r(expr.constant());
    for (typename linear_expression_t::const_iterator it = expr.begin();
         it != expr.end(); ++it) {
      interval_t c(it->first);
      r += c * this->_env.at(it->second);
    }
    return r;
  }
  void add(const linear_constraint_system_t &csts,
           std::size_t threshold = max_reduction_cycles) {
    if (!this->is_bottom()) {
      linear_constraint_system_t pp_csts;
      for (auto const &c : csts) {
	if (c.is_disequation()) {
	  // We try to convert a disequation into a strict inequality
	  crab::domains::constraint_simp_domain_traits<interval_domain_t>::
	    lower_disequality(*this, c, pp_csts);
	}
        pp_csts += c;
      }
      solver_t solver(pp_csts, threshold);
      solver.run(this->_env);
    }
  }

  interval_domain_t operator+(const linear_constraint_system_t &csts) {
    interval_domain_t e(this->_env);
    e += csts;
    return e;
  }

public:
  interval_domain_t make_top() const override {
    return interval_domain_t(separate_domain_t::top());
  }

  interval_domain_t make_bottom() const override {
    return interval_domain_t(separate_domain_t::bottom());
  }

  void set_to_top() override {
    interval_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    interval_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  interval_domain() : _env(separate_domain_t::top()) {}

  interval_domain(const interval_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  interval_domain_t &operator=(const interval_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }

  iterator begin() { return this->_env.begin(); }

  iterator end() { return this->_env.end(); }

  bool is_bottom() const override { return this->_env.is_bottom(); }

  bool is_top() const override { return this->_env.is_top(); }

  bool operator<=(const interval_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (this->_env <= e._env);
  }

  void operator|=(const interval_domain_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    this->_env = this->_env | e._env;
  }

  interval_domain_t operator|(const interval_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    return (this->_env | e._env);
  }

  void operator&=(const interval_domain_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    this->_env = this->_env & e._env;
  }
  
  interval_domain_t operator&(const interval_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    return (this->_env & e._env);
  }

  interval_domain_t operator||(const interval_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return (this->_env || e._env);
  }

  interval_domain_t widening_thresholds(
      const interval_domain_t &e,
      const crab::thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return this->_env.widening_thresholds(e._env, ts);
  }

  interval_domain_t operator&&(const interval_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return (this->_env && e._env);
  }

  void set(const variable_t &v, interval_t i) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    this->_env.set(v, i);
  }

  void set(const variable_t &v, number_t n) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    this->_env.set(v, interval_t(n));
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    this->_env -= v;
  }

  interval_t operator[](const variable_t &v) override { return at(v); }

  interval_t at(const variable_t &v) const override {
    return this->_env.at(v);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    this->add(csts);
  }

  virtual bool entails(const linear_constraint_t &cst) const override {	
    if (is_bottom()) {							
      return true;							
    } if (cst.is_tautology()) {						
      return true;							
    } if (cst.is_contradiction()) {					
      return false;							
    }

    // val is modified after the check
    auto entailmentFn = [](interval_domain_t &val, const linear_constraint_t &c) -> bool {	
      linear_constraint_t neg_c = c.negate();			        
      val += neg_c;						        
      return val.is_bottom();					        
    };

    // Get only relevant state wrt cst variables
    interval_domain_t val;
    for (auto const&v: cst.variables()) {
      val.set(v, at(v));
    }

    if (cst.is_equality()) {						
      linear_constraint_t pob1(cst.expression(),linear_constraint_t::INEQUALITY);
      interval_domain_t tmp(val);
      if (!entailmentFn(tmp, pob1)) {
	return false;
      }
      linear_constraint_t pob2(cst.expression() * number_t(-1), linear_constraint_t::INEQUALITY);   	 
      return entailmentFn(val, pob2);		
    } else {								
      return entailmentFn(val, cst);						
    }									
  }
  
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (boost::optional<variable_t> v = e.get_variable()) {
      this->_env.set(x, this->_env.at(*v));
    } else {
      interval_t r = e.constant();
      for (auto kv : e) {
        r += kv.first * this->_env.at(kv.second);
      }
      this->_env.set(x, r);
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.weak_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".weak_assign");

    if (boost::optional<variable_t> v = e.get_variable()) {
      this->_env.join(x, this->_env.at(*v));
    } else {
      interval_t r = e.constant();
      for (auto kv : e) {
        r += kv.first * this->_env.at(kv.second);
      }
      this->_env.join(x, r);
    }
  }

  
  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env.at(y);
    interval_t zi = this->_env.at(z);
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
    default:
      //case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    }
    this->_env.set(x, xi);
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env.at(y);
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
    default:
      //case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    }
    this->_env.set(x, xi);
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
                          const interval_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  // backward arithmetic operations
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const interval_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    crab::domains::BackwardAssignOps<interval_domain_t>::assign(*this, x, e,
                                                                inv);
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const interval_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const interval_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  // cast operations
  void apply(crab::domains::int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::domains::int_cast_domain_traits<interval_domain_t>::apply(*this, op, dst, src);
  }

  // bitwise operations
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env.at(y);
    interval_t zi = this->_env.at(z);
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
    default:
      //case crab::domains::OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    this->_env.set(x, xi);
  }

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env.at(y);
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
    default:
      //case crab::domains::OP_ASHR: 
      xi = yi.AShr(zi);
      break;
    }
    this->_env.set(x, xi);
  }

  virtual void select(const variable_t &lhs, const linear_constraint_t &cond,
		      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    crab::CrabStats::count(domain_name() + ".count.select");
    crab::ScopedCrabStats __st__(domain_name() + ".select");
    
    if (!is_bottom()) {
      interval_domain_t inv1(*this);
      inv1 += cond;
      if (inv1.is_bottom()) {
	assign(lhs, e2);
	return;
      }
      
      interval_domain_t inv2(*this);
      inv2 += cond.negate();
      if (inv2.is_bottom()) {
	assign(lhs, e1);
	return;
      }

      set(lhs, this->operator[](e1) | this->operator[](e2));
    }
  }
  
  /// interval_domain implements only standard abstract operations of
  /// a numerical domain so it is intended to be used as a leaf domain
  /// in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(interval_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(interval_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(interval_domain_t)

  
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

    _env.project(variables);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    _env.rename(from, to);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set(new_x, this->_env.at(x));
  }

  void normalize() override {}

  void minimize() override {}

  void write(crab::crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    this->_env.write(o);
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

    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      const variable_t &v = it->first;
      interval_t i = it->second;
      boost::optional<number_t> lb = i.lb().number();
      boost::optional<number_t> ub = i.ub().number();
      if (lb)
        csts += linear_constraint_t(v >= *lb);
      if (ub)
        csts += linear_constraint_t(v <= *ub);
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

  std::string domain_name() const override { return "Intervals"; }

}; // class interval_domain
} // namespace ikos

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
struct abstract_domain_traits<ikos::interval_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // namespace domains
} // namespace crab
