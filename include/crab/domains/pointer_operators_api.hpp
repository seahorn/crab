/*******************************************************************************
 * Generic API for pointer operations
 ******************************************************************************/

#pragma once

#include <crab/common/types.hpp>
#include <crab/domains/linear_constraints.hpp>


namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class pointer_operators {
       public:

        typedef ikos::linear_expression<Number, VariableName> lin_exp_t;
	typedef ikos::variable<Number, VariableName> variable_t;
	
        typedef pointer_constraint<variable_t> ptr_cst_t;
        
        virtual ~pointer_operators () {}
        
        // p := *q 
        virtual void pointer_load (variable_t lhs, variable_t rhs) {}
        // *p := q
        virtual void pointer_store (variable_t lhs, variable_t rhs) {} 
        // p := q + n
        virtual void pointer_assign (variable_t lhs, variable_t rhs, lin_exp_t offset) {}
        // p := &a;
        virtual void pointer_mk_obj (variable_t lhs, ikos::index_t address) {}
        // p := &func
        virtual void pointer_function (variable_t lhs, VariableName func) {}
        // p := null
        virtual void pointer_mk_null (variable_t lhs) {}
        // assume (cst);
        virtual void pointer_assume (ptr_cst_t cst) {}
        // assert (cst);
        virtual void pointer_assert (ptr_cst_t cst) {}
      };

   } // namespace domains  
} // namespace crab
