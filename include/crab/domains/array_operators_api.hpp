/*******************************************************************************
 * Generic API for array operations
 ******************************************************************************/

#pragma once 

#include <crab/common/bignums.hpp>
#include <crab/domains/linear_constraints.hpp>

#include <vector>

namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class array_operators {

       protected:

        typedef ikos::linear_expression<Number, VariableName> _linear_exp_t;
	typedef ikos::variable<Number, VariableName> variable_t;
	
       public:

        virtual ~array_operators () { }

        // initialize all array contents in [lb_idx,ub_idx] to val
        // where elem_size is in bytes.
        virtual void array_init (variable_t a, _linear_exp_t elem_size,
				 _linear_exp_t lb_idx, _linear_exp_t ub_idx, 
				 _linear_exp_t val) {}

        // lhs := a[i] where elem_size is in bytes
        virtual void array_load (variable_t lhs,
				 variable_t a, _linear_exp_t elem_size,
				 _linear_exp_t i) {}  

        // a[i] := v where elem_size is in bytes
        virtual void array_store (variable_t a, _linear_exp_t elem_size,
				  _linear_exp_t i, _linear_exp_t v, 
                                  bool is_singleton) {}

        // a := b  (forall i :: a[i] := b[i])
        virtual void array_assign (variable_t lhs, variable_t rhs) {}
                            
      };

   } // namespace domains  
} // namespace crab
