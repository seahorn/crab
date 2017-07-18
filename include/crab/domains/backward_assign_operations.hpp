#ifndef GENERIC_BACKWARD_ASSIGNMENT_HPP
#define GENERIC_BACKWARD_ASSIGNMENT_HPP

/**
 * Implement generic backward assignments. Unless the assignment is
 * invertible the result is an over-approximation so we need to adapt
 * these operations in case we need under-approximations.
 **/


#include <crab/config.h>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>

namespace crab {
  namespace domains {

    template <class AbsDom>
    class BackwardAssignOps {
    public:
      typedef typename AbsDom::number_t number_t;
      typedef typename AbsDom::varname_t varname_t;
      typedef typename AbsDom::variable_t variable_t;
      typedef typename AbsDom::linear_constraint_t linear_constraint_t;
      typedef typename AbsDom::linear_expression_t linear_expression_t;

      /*
       * Backward x := e 
       *  case 1: if e is invertible
       *    x = y + k <--> y = x - k     
       *    x = y - k <--> y = x + k
       *    x = y * k <--> y = x / k  if (k != 0)
       *    x = y / k <--> y = x * k  if (k != 0)
       *  case 2: if x does not appear in e
       *   1) add constraint x = e
       *   2) forget x
       *  case 3: if x appears in e 
       *   TODO: (this works only for very expressive domains)
       *   1) add new variable x'
       *   1) add constraint x = e[x'/x]
       *   3) forget x 
       *   4) rename x' as x
       **/ 
      
      // x := e
      // pre: dom is not bottom
      static void assign(AbsDom& dom, varname_t x, linear_expression_t e,
			 AbsDom inv) {
	if (e.variables() >= variable_t(x)) {
	  CRAB_WARN("backwards x:=e when x appears in e not implemented");
	} else {
	  dom = dom & inv;
	  linear_expression_t x_e = e - variable_t(x);
	  dom += linear_constraint_t(x_e, linear_constraint_t::EQUALITY);
	  dom -= x;
	}
      }

      // x := y op k
      // pre: dom is not bottom
      static void apply(AbsDom& dom, operation_t op,
			varname_t x, varname_t y, number_t k,
			AbsDom inv) {
	switch(op) {
	  case OP_ADDITION: {
	    dom = dom & inv;	
	    dom.apply(OP_SUBTRACTION, y, x, k);
	    break;
	  }
   	  case OP_SUBTRACTION: {
	    dom = dom & inv;		    
	    dom.apply(OP_ADDITION, y, x, k);
	    break;
	  }
	  case OP_MULTIPLICATION: { 
	    if(k != 0) {
	      dom = dom & inv;		      
	      dom.apply(OP_DIVISION, y, x, k);
	    } else
	      CRAB_WARN("backwards x:= y * k is not invertible");
	    break;
	  }
	  case OP_DIVISION: { 
	    if(k != 0) {
	      dom = dom & inv;		      
	      dom.apply(OP_MULTIPLICATION, y, x, k);	      
	    } else
	      CRAB_WARN("backwards x:= y / k is not invertible");
	    break;
	  }
  	  default:;;  
	}
	return;
      }

      // x = y op z
      static void apply(AbsDom& dom, operation_t op,
			varname_t x, varname_t y, varname_t z,
			AbsDom inv) {
	if (x == y || x == z) {
	  CRAB_WARN("backwards x:=e when x appears in e not implemented");
	} else {
	  switch(op) {
	  case OP_ADDITION:
	    dom = dom & inv;	
	    dom += linear_constraint_t(variable_t(y) + variable_t(z) - variable_t(x),
				       linear_constraint_t::EQUALITY);
	    dom -= x;
	    break;
	  case OP_SUBTRACTION:
	    dom = dom & inv;		    
	    dom += linear_constraint_t(variable_t(y) - variable_t(z) - variable_t(x),
				       linear_constraint_t::EQUALITY);
	    dom -= x;
	    break;
	  case OP_MULTIPLICATION:
	    CRAB_WARN("backwards x = y * z not implemented");
	    break;	    
	  case OP_DIVISION:
	    CRAB_WARN("backwards x = y / z not implemented");
	    break;
	  }
	}
      }
    };

  } //end namespace domains
} // end namespace crab

#endif 
