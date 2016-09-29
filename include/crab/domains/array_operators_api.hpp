/*******************************************************************************
 * Generic API for array operations
 ******************************************************************************/

#ifndef ARRAY_OPERATORS_API_HPP
#define ARRAY_OPERATORS_API_HPP

#include <crab/common/bignums.hpp>
//#include <crab/domains/linear_constraints.hpp>

#include <vector>

namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class array_operators {

       protected:

        typedef ikos::linear_expression<Number, VariableName> _linear_exp_t;

       public:
        
        virtual ~array_operators () { }

        /* --- begin (temporary) operations for array smashing */
        // This need to be fixed
        // smashes a sequence of array stores 
        virtual void array_init (VariableName a, const std::vector<ikos::z_number>& vals) {}
        // assume all array contents are in [*lb_val, *ub_val]
        // if lb_val is null then -oo 
        // if ub_val is null then +oo
        // XXX: we do not use bound or interval because we don't want
        //      to introduce circular dependencies between header files
        virtual void array_assume (VariableName a, 
                                   boost::optional<Number> lb_val, boost::optional<Number> ub_val) {}
        /* end operations for array smashing */


        // Pre: a is an array of type a_ty
        // lhs := a[i] where bytes is the size of a elements in bytes
        virtual void array_load (VariableName lhs, VariableName a, variable_type a_ty, _linear_exp_t i, 
                                 ikos::z_number bytes) {}

        // Pre: a is an array of type a_ty
        // a[i] := v where bytes is the size of a elements in bytes
        virtual void array_store (VariableName a, variable_type a_ty, _linear_exp_t i, VariableName v, 
                                  ikos::z_number bytes, bool is_singleton) {}

        // Pre: lhs and rhs are arrays of the same type ty
        // a := b  (forall i :: a[i] := b[i])
        virtual void array_assign (VariableName lhs, VariableName rhs, variable_type ty) {}
                            
      };

   } // namespace domains  
} // namespace crab
#endif 
