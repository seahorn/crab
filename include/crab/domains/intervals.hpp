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
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/domains/separate_domains.hpp>

namespace ikos {

template <typename Number, typename VariableName,
          std::size_t max_reduction_cycles = 10>
class interval_domain final
    : public crab::domains::abstract_domain<
          interval_domain<Number, VariableName, max_reduction_cycles>> {
public:
  typedef interval_domain<Number, VariableName, max_reduction_cycles>
      interval_domain_t;
  typedef crab::domains::abstract_domain<interval_domain_t> abstract_domain_t;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef interval<number_t> interval_t;

private:
  typedef separate_domain<variable_t, interval_t> separate_domain_t;
  typedef linear_interval_solver<number_t, varname_t, separate_domain_t>
      solver_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  interval_domain(separate_domain_t env) : _env(env) {}

public:
  void set_to_top() {
    interval_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    interval_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  interval_domain() : _env(separate_domain_t::top()) {}

  interval_domain(const interval_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  interval_domain_t &operator=(const interval_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }

  iterator begin() { return this->_env.begin(); }

  iterator end() { return this->_env.end(); }

  bool is_bottom() { return this->_env.is_bottom(); }

  bool is_top() { return this->_env.is_top(); }

  bool operator<=(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    return (this->_env <= e._env);
  }

  void operator|=(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    this->_env = this->_env | e._env;
  }

  interval_domain_t operator|(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    return (this->_env | e._env);
  }

  interval_domain_t operator&(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");
    return (this->_env & e._env);
  }

  interval_domain_t operator||(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    return (this->_env || e._env);
  }

  interval_domain_t
  widening_thresholds(interval_domain_t e,
                      const crab::iterators::thresholds<number_t> &ts) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    return this->_env.widening_thresholds(e._env, ts);
  }

  interval_domain_t operator&&(interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    return (this->_env && e._env);
  }

  void set(variable_t v, interval_t i) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, i);
  }

  void set(variable_t v, number_t n) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, interval_t(n));
  }

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
    this->_env -= v;
  }

  interval_t operator[](variable_t v) { return this->_env[v]; }

  interval_t operator[](linear_expression_t expr) {
    interval_t r(expr.constant());
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      interval_t c(it->first);
      r += c * this->_env[it->second];
    }
    return r;
  }

  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");
    this->add(csts);
  }

  void add(linear_constraint_system_t csts,
           std::size_t threshold = max_reduction_cycles) {
    if (!this->is_bottom()) {
      // XXX: filter out unsigned linear inequalities
      linear_constraint_system_t signed_csts;
      for (auto const &c : csts) {
        if (c.is_inequality() && c.is_unsigned()) {
          // CRAB_WARN("unsigned inequality skipped");
          continue;
        }
        signed_csts += c;
      }
      solver_t solver(signed_csts, threshold);
      solver.run(this->_env);
    }
  }

  interval_domain_t operator+(linear_constraint_system_t csts) {
    interval_domain_t e(this->_env);
    e += csts;
    return e;
  }

  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (boost::optional<variable_t> v = e.get_variable()) {
      this->_env.set(x, this->_env[(*v)]);
    } else {
      interval_t r = e.constant();
      for (typename linear_expression_t::iterator it = e.begin(); it != e.end();
           ++it) {
        r += it->first * this->_env[it->second];
      }
      this->_env.set(x, r);
    }
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi = this->_env[z];
    interval_t xi = interval_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi(k);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  // intrinsics operations
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  interval_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  
  // backward arithmetic operations
  void backward_assign(variable_t x, linear_expression_t e,
                       interval_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign");

    crab::domains::BackwardAssignOps<interval_domain_t>::assign(*this, x, e,
                                                                inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      interval_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      interval_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  /*
     Begin unimplemented operations

     interval_domain implements only standard abstract
     operations of a numerical domain.  The implementation of
     boolean, array, or pointer operations is empty because they
     should never be called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(crab::domains::bool_operation_t op, variable_t x,
                         variable_t y, variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                interval_domain_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                interval_domain_t invariant) {}
  void backward_apply_binary_bool(crab::domains::bool_operation_t op,
                                  variable_t x, variable_t y, variable_t z,
                                  interval_domain_t invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           interval_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           interval_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            interval_domain_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            interval_domain_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  interval_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  interval_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             interval_domain_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  // cast operations
  void apply(crab::domains::int_conv_operation_t /*op*/, variable_t dst,
             variable_t src) {
    // ignore the widths
    assign(dst, src);
  }

  // bitwise operations
  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi = this->_env[z];
    interval_t xi = interval_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi(k);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    separate_domain_t env;
    for (variable_t var : variables) {
      env.set(var, this->_env[var]);
    }
    std::swap(_env, env);
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    _env.rename(from, to);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set(new_x, this->_env[x]);
  }

  void normalize() {}

  void minimize() {}

  void write(crab::crab_os &o) {
    crab::CrabStats::count(getDomainName() + ".count.write");
    crab::ScopedCrabStats __st__(getDomainName() + ".write");

    this->_env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;

    if (this->is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      variable_t v = it->first;
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
  to_disjunctive_linear_constraint_system() {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  static std::string getDomainName() { return "Intervals"; }

}; // class interval_domain
} // namespace ikos

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
struct abstract_domain_traits<ikos::interval_domain<Number, VariableName>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

} // namespace domains
} // namespace crab
