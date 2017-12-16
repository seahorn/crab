/*******************************************************************************
 *
 * Resolution of a system of linear constraints over the domain of intervals
 * is based on W. Harvey & P. J. Stuckey's paper: Improving linear constraint
 * propagation by changing constraint representation, in Constraints,
 * 8(2):173–207, 2003.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Contributors: Alexandre C. D. Wimmers (alexandre.c.wimmers@nasa.gov)
 *               Jorge Navas (jorge.navas@sri.com)
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

#include <vector>
#include <set>
#include <map>
#include <boost/optional.hpp>
#include <crab/common/types.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/bignums.hpp>
#include <crab/common/wrapint.hpp>
#include <crab/domains/linear_constraints.hpp>

namespace ikos {

  namespace linear_interval_solver_impl {

    template<typename Interval, typename Number>
    inline Interval trim_bound(Interval i, Number c);

    template<typename Interval, typename Number>
    inline Interval mk_interval(Number n, typename crab::wrapint::bitwidth_t /*bitwidth*/) {
      // default implementation ignores bitwidth
      return Interval(n);
    }
  } 

  template< typename Number, typename VariableName, typename IntervalCollection >
  class linear_interval_solver {
    
  public:
    typedef typename IntervalCollection::value_type Interval; 
    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    typedef typename variable_t::bitwidth_t bitwidth_t;
    
  private:
    typedef std::vector< linear_constraint_t > cst_table_t;
    typedef std::set< unsigned int > uint_set_t;
    typedef std::map< variable_t, uint_set_t > trigger_table_t;
    typedef typename linear_constraint_t::variable_set_t variable_set_t;

  private:
    class bottom_found { };

  private:
    std::size_t _max_cycles;
    std::size_t _max_op;
    bool _is_contradiction;
    bool _is_large_system;
    cst_table_t _cst_table;
    trigger_table_t _trigger_table;
    variable_set_t _refined_variables;
    std::size_t _op_count;
    
  private:
    static const std::size_t _large_system_cst_threshold = 3;
    // cost of one propagation cycle for a dense 3x3 system of constraints 
    static const std::size_t _large_system_op_threshold = 27; 

  private:
    void refine(variable_t v, Interval i, IntervalCollection& env) {
      Interval old_i = env[v];
      Interval new_i = old_i & i;
      if (new_i.is_bottom()) {
	throw bottom_found();
      }
      if (!(old_i == new_i)) {
	env.set(v, new_i);
	this->_refined_variables += v;
	++(this->_op_count);
      }
    }

    Interval compute_residual(linear_constraint_t cst, variable_t pivot, 
                              IntervalCollection& env) {
      bitwidth_t w = pivot.get_bitwidth();
      Interval residual= linear_interval_solver_impl::mk_interval<Interval>(cst.constant(), w);
      for (typename linear_constraint_t::iterator it = cst.begin(); 
           it != cst.end(); ++it) {
	variable_t v = it->second;
	if (!(v == pivot)) {
	  residual = residual - (linear_interval_solver_impl::mk_interval<Interval>(it->first, w) * env[v]);
	  ++(this->_op_count);
	}
      }
      return residual;
    }
    
    void propagate(linear_constraint_t cst, IntervalCollection& env) {
      for (typename linear_constraint_t::iterator it = cst.begin(); 
           it != cst.end(); ++it) {
	Number c = it->first;
	variable_t pivot = it->second;	
	Interval ic = linear_interval_solver_impl::mk_interval<Interval>(c, pivot.get_bitwidth());	
	Interval rhs = compute_residual(cst, pivot, env) / ic;
	if (cst.is_equality()) {
	  this->refine(pivot, rhs, env);
	} else if (cst.is_inequality()) {
	  if (c > 0) {
	    refine(pivot, rhs.lower_half_line(), env);
	  } else {
	    refine(pivot, rhs.upper_half_line(), env);
	  }
	} else if (cst.is_strict_inequality()) {
	  // do nothing
	} else {
	  // cst is a disequation
	  boost::optional<Number> c = rhs.singleton();
	  if (c) {
	    Interval old_i = env[pivot];
	    Interval new_i = linear_interval_solver_impl::trim_bound (old_i, *c);
	    if (new_i.is_bottom()) {
	      throw bottom_found();
	    }
	    if (!(old_i == new_i)) {
	      env.set(pivot, new_i);
	      this->_refined_variables += pivot;
	    }
	    ++(this->_op_count);
	  }
	}
      }
    }
    
    void solve_large_system(IntervalCollection& env) {
      this->_op_count = 0;
      this->_refined_variables.clear();
      for (typename cst_table_t::iterator it = this->_cst_table.begin(); 
           it != this->_cst_table.end(); ++it) {
	this->propagate(*it, env);
      }
      do {
	variable_set_t vars_to_process(this->_refined_variables);
	this->_refined_variables.clear();
	for (typename variable_set_t::iterator it = vars_to_process.begin(); 
               it != vars_to_process.end(); ++it) {
	  uint_set_t& csts = this->_trigger_table[*it];
	  for (typename uint_set_t::iterator cst_it = csts.begin(); 
                 cst_it != csts.end(); ++cst_it) {
	    this->propagate(this->_cst_table.at(*cst_it), env);
	  }
	}
      }
      while (!this->_refined_variables.empty() && 
             this->_op_count <= this->_max_op);
    }

    void solve_small_system(IntervalCollection& env) {
      std::size_t cycle = 0;
      do {
	++cycle;
	this->_refined_variables.clear();
	for (typename cst_table_t::iterator it = this->_cst_table.begin(); 
               it != this->_cst_table.end(); ++it) {
	  this->propagate(*it, env);
	}
      }
      while (!this->_refined_variables.empty() &&  cycle <= this->_max_cycles);
    }
    
    
  public:

    linear_interval_solver(linear_constraint_system_t csts, 
                           std::size_t max_cycles): 
        _max_cycles(max_cycles), 
        _is_contradiction(false), 
        _is_large_system(false), 
        _op_count(0) {

      std::size_t op_per_cycle = 0;
      for (typename linear_constraint_system_t::iterator it = csts.begin(); 
           it != csts.end(); ++it) {
	linear_constraint_t cst = *it;
	if (cst.is_contradiction()) {
	  this->_is_contradiction = true;
	  return;
	} else if (cst.is_tautology()) {
	  continue;
	} else {
	  std::size_t cst_size = cst.size();
	  if (cst.is_strict_inequality()) {
	    // convert e < c into {e <= c, e != c}
	    linear_constraint_t c1(cst.expression(), linear_constraint_t::kind_t::INEQUALITY);
	    linear_constraint_t c2(cst.expression(), linear_constraint_t::kind_t::DISEQUATION);
	    this->_cst_table.push_back(c1);
	    this->_cst_table.push_back(c2);
	    cst_size = c1.size() + c2.size();
	  } else {
	    this->_cst_table.push_back(cst);
	  }
	  // cost of one reduction step on the constraint in terms
	  // of accesses to the interval collection
	  op_per_cycle += cst_size * cst_size; 
	}
      }

      this->_is_large_system = (this->_cst_table.size() > 
                                _large_system_cst_threshold) || 
          (op_per_cycle > _large_system_op_threshold);
      
      if (!this->_is_contradiction && this->_is_large_system) {
	this->_max_op = op_per_cycle * max_cycles;
	for (unsigned int i = 0; i < this->_cst_table.size(); ++i) {
	  linear_constraint_t cst = this->_cst_table.at(i);
	  variable_set_t vars = cst.variables();
	  for (typename variable_set_t::iterator it = vars.begin(); 
                 it != vars.end(); ++it) {
	    this->_trigger_table[*it].insert(i);
	  }
	}
      }
    }
    
    void run(IntervalCollection& env) {
      if (this->_is_contradiction) {
        env.set_to_bottom();
      } else {
        try {
          if (this->_is_large_system) {
	    this->solve_large_system(env);
          } else {
            this->solve_small_system(env);
          }
        }
        catch (bottom_found& e) {
          env.set_to_bottom();
        }
      }
    }
    
  }; // class linear_interval_solver

} // namespace ikos

