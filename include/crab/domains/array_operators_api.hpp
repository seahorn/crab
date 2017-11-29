/*******************************************************************************
 * Generic API for array operations
 ******************************************************************************/

#pragma once 

#include <crab/common/bignums.hpp>
//#include <crab/domains/linear_constraints.hpp>

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

        // Pre: a is an array of type a_ty
        // assume all array contents in [lb_idx,ub_idx] are equal to val
        virtual void array_assume (variable_t a, variable_type a_ty, 
                                   _linear_exp_t lb_idx, _linear_exp_t ub_idx, 
                                   _linear_exp_t val) {}

        // Pre: a is an array of type a_ty
        // lhs := a[i] where bytes is the size of a elements in bytes
        virtual void array_load (variable_t lhs, variable_t a, variable_type a_ty,
				 _linear_exp_t i, ikos::z_number bytes) {}  

        // Pre: a is an array of type a_ty
        // a[i] := v where bytes is the size of a elements in bytes
        virtual void array_store (variable_t a, variable_type a_ty,
				  _linear_exp_t i, _linear_exp_t v, 
                                  ikos::z_number bytes, bool is_singleton) {}

        // Pre: lhs and rhs are arrays of the same type ty
        // a := b  (forall i :: a[i] := b[i])
        virtual void array_assign (variable_t lhs, variable_t rhs, variable_type ty) {}
                            
      };

   } // namespace domains  
} // namespace crab
