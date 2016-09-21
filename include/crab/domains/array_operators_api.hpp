/*******************************************************************************
 * Generic API for array operations
 ******************************************************************************/

#ifndef ARRAY_OPERATORS_API_HPP
#define ARRAY_OPERATORS_API_HPP

#include <crab/common/bignums.hpp>
#include <crab/domains/linear_constraints.hpp>

#include <vector>

namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class array_operators {

       protected:

        typedef ikos::linear_expression<Number, VariableName> _linear_exp_t;

       public:
        
        virtual ~array_operators () { }

        // smashes a sequence of array stores 
        virtual void array_init (VariableName a, const std::vector<ikos::z_number>& vals) {}
        // assume all array contents are in [*lb_val, *ub_val]
        // if lb_val is null then -oo 
        // if ub_val is null then +oo
        // XXX: we do not use bound or interval because we don't want
        //      to introduce circular dependencies between header files
        virtual void array_assume (VariableName a, 
                                   boost::optional<Number> lb_val, boost::optional<Number> ub_val) {}

        virtual void array_load (VariableName lhs, VariableName a, VariableName i, 
                                 ikos::z_number bytes) {}
        virtual void array_store (VariableName a, VariableName i, _linear_exp_t v, 
                                  ikos::z_number bytes, bool is_singleton) {}
                                  
      };

   } // namespace domains  
} // namespace crab
#endif 
