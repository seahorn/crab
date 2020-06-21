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
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/stats.hpp>

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
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::interval_t;
  typedef Number number_t;
  typedef VariableName varname_t;

private:
  typedef separate_domain<variable_t, interval_t> separate_domain_t;
  typedef linear_interval_solver<number_t, varname_t, separate_domain_t>
      solver_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  interval_domain(separate_domain_t env) : _env(env) {}

  interval_t operator[](linear_expression_t expr) {
    interval_t r(expr.constant());
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      interval_t c(it->first);
      r += c * this->_env[it->second];
    }
    return r;
  }
  void add(const linear_constraint_system_t &csts,
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

  interval_domain_t operator+(const linear_constraint_system_t &csts) {
    interval_domain_t e(this->_env);
    e += csts;
    return e;
  }
  
public:
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

  bool operator<=(interval_domain_t e) override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return (this->_env <= e._env);
  }

  void operator|=(interval_domain_t e) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    this->_env = this->_env | e._env;
  }

  interval_domain_t operator|(interval_domain_t e) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    return (this->_env | e._env);
  }

  interval_domain_t operator&(interval_domain_t e) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    return (this->_env & e._env);
  }

  interval_domain_t operator||(interval_domain_t e) override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return (this->_env || e._env);
  }

  interval_domain_t
  widening_thresholds(interval_domain_t e,
                      const crab::iterators::thresholds<number_t> &ts) override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return this->_env.widening_thresholds(e._env, ts);
  }

  interval_domain_t operator&&(interval_domain_t e) override {
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

  interval_t operator[](const variable_t &v) override {
    return this->_env[v];
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    this->add(csts);
  }


  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (boost::optional<variable_t> v = e.get_variable()) {
      this->_env.set(x, this->_env[(*v)]);
    } else {
      interval_t r = e.constant();
      for (auto kv: e) {
        r += kv.first * this->_env[kv.second];
      }
      this->_env.set(x, r);
    }
  }

  void apply(crab::domains::arith_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi = this->_env[z];
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
    this->_env.set(x, xi);
  }

  void apply(crab::domains::arith_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env[y];
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
    this->_env.set(x, xi);
  }

  // intrinsics operations
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  interval_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());    
  }
  
  // backward arithmetic operations
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       interval_domain_t inv) override  {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    crab::domains::BackwardAssignOps<interval_domain_t>::assign(*this, x, e,
                                                                inv);
  }

  void backward_apply(crab::domains::arith_operation_t op,
		      const variable_t &x, const variable_t &y, number_t z,
                      interval_domain_t inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  void backward_apply(crab::domains::arith_operation_t op,
		      const variable_t &x, const variable_t &y, const variable_t &z,
                      interval_domain_t inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<interval_domain_t>::apply(*this, op, x, y,
                                                               z, inv);
  }

  /*
     Begin unimplemented operations

     interval_domain implements only standard abstract
     operations of a numerical domain. 
  */

  // boolean operations
  void assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs) override {}
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs) override {}
  void apply_binary_bool(crab::domains::bool_operation_t op,
			 const variable_t &x, const variable_t &y, const variable_t &z) override {}
  void assume_bool(const variable_t &v, bool is_negated) override {}
  // backward boolean operations
  void backward_assign_bool_cst(const variable_t &lhs, const linear_constraint_t &rhs,
                                interval_domain_t invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs, bool is_not_rhs,
                                interval_domain_t invariant) override {}
  void backward_apply_binary_bool(crab::domains::bool_operation_t op,
                                  const variable_t &x, const variable_t &y, const variable_t &z,
                                  interval_domain_t invariant) override {}
  // array operations
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  void array_load(const variable_t &lhs, const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    operator-=(lhs);
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {}
  void array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                         const linear_expression_t &i, const linear_expression_t &j,
                         const linear_expression_t &v) override {}
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {}
  // backward array operations
  void backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx, const linear_expression_t &val,
                           interval_domain_t invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size, const linear_expression_t &i,
                           interval_domain_t invariant) override {}
  void backward_array_store(const variable_t &a, const linear_expression_t &elem_size,
                            const linear_expression_t &i, const linear_expression_t &v,
                            bool is_strong_update,
                            interval_domain_t invariant) override {}
  void backward_array_store_range(const variable_t &a, const linear_expression_t &elem_size,
                                  const linear_expression_t &i, const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  interval_domain_t invariant) override {}
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             interval_domain_t invariant) override {}
  // reference operations
  void region_init(const crab::memory_region &reg) override {}          
  void ref_make(const variable_t &ref, const crab::memory_region &reg) override {}
  void ref_load(const variable_t &ref, const crab::memory_region &reg,
		const variable_t &res) override {}
  void ref_store(const variable_t &ref, const crab::memory_region &reg,
		 const linear_expression_t &val) override {}
  void ref_gep(const variable_t &ref1, const crab::memory_region &reg1,
	       const variable_t &ref2, const crab::memory_region &reg2,
	       const linear_expression_t &offset) override {}
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
			   const crab::memory_region &region,
			   const linear_expression_t &index,
			   const linear_expression_t &elem_size) override {}
  void ref_store_to_array(const variable_t &ref, const crab::memory_region &region,
			  const linear_expression_t &index, const linear_expression_t &elem_size,
			  const linear_expression_t &val) override {}
  void ref_assume(const reference_constraint_t &cst) override {}
  /* End unimplemented operations */

  // cast operations
  void apply(crab::domains::int_conv_operation_t /*op*/,
	     const variable_t &dst, const variable_t &src) override {
    // ignore the widths
    assign(dst, src);
  }

  // bitwise operations
  void apply(crab::domains::bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env[y];
    interval_t zi = this->_env[z];
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
    this->_env.set(x, xi);
  }

  void apply(crab::domains::bitwise_operation_t op,
	     const variable_t &x, const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    interval_t yi = this->_env[y];
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
    this->_env.set(x, xi);
  }

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    separate_domain_t env;
    for (variable_t var : variables) {
      env.set(var, this->_env[var]);
    }
    std::swap(_env, env);
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) override {
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

    set(new_x, this->_env[x]);
  }

  void normalize() override {}

  void minimize() override {}

  void write(crab::crab_os &o) override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    this->_env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() override {
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
  to_disjunctive_linear_constraint_system() override {
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
  typedef Number number_t;
  typedef VariableName varname_t;
};

} // namespace domains
} // namespace crab
