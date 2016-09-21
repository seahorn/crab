/*******************************************************************************
 * Generic API for pointer operations
 ******************************************************************************/

#ifndef POINTER_OPERATORS_API_HPP
#define POINTER_OPERATORS_API_HPP

#include <crab/common/types.hpp>
#include <crab/domains/linear_constraints.hpp>


namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class pointer_operators {
       public:

        typedef ikos::linear_expression<Number, VariableName> lin_exp_t;
        typedef pointer_constraint<VariableName> ptr_cst_t;
        
        virtual ~pointer_operators () {}
        
        // p := *q 
        virtual void pointer_load (VariableName lhs, VariableName rhs) {}
        // *p := q
        virtual void pointer_store (VariableName lhs, VariableName rhs) {} 
        // p := q + n
        virtual void pointer_assign (VariableName lhs, VariableName rhs, lin_exp_t offset) {}
        // p := &a;
        virtual void pointer_mk_obj (VariableName lhs, ikos::index_t address) {}
        // p := &func
        virtual void pointer_function (VariableName lhs, VariableName func) {}
        // p := null
        virtual void pointer_mk_null (VariableName lhs) {}
        // assume (cst);
        virtual void pointer_assume (ptr_cst_t cst) {}
        // assert (cst);
        virtual void pointer_assert (ptr_cst_t cst) {}
      };

   } // namespace domains  
} // namespace crab
#endif 
