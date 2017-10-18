/*******************************************************************************
 * Generic API for boolean operations
 ******************************************************************************/

#ifndef BOOLEAN_OPERATORS_API_HPP
#define BOOLEAN_OPERATORS_API_HPP

#include <crab/common/types.hpp>

namespace crab {

   namespace domains {

     typedef enum  { 
       OP_BAND, 
       OP_BOR, 
       OP_BXOR
    } bool_operation_t;

     inline crab::crab_os& operator<<(crab::crab_os&o, bool_operation_t op) {
       switch (op) {
         case OP_BAND: o << "&"; break;
         case OP_BOR : o << "|"; break;
         case OP_BXOR: o << "^"; break;
         default     : ;;
       }
       return o;
     }

     // XXX: we deliberately choose not to include a template
     // parameter Bool and have methods that take arguments of that
     // type.
     template<typename Number, typename VariableName>
     class boolean_operators {
     public:
       typedef ikos::linear_constraint<Number, VariableName> lin_cst_t;       
       
       virtual ~boolean_operators () {};
       
       virtual void assign_bool_cst(VariableName lhs, lin_cst_t rhs) {}
       virtual void assign_bool_var(VariableName lhs, VariableName rhs, bool is_not_rhs) {}
       virtual void apply_binary_bool(bool_operation_t op,
				      VariableName x,VariableName y,VariableName z) {}
       virtual void assume_bool(VariableName v, bool is_negated) {}
     };

     template<typename Number, typename VariableName, typename NumAbsDom>
     class backward_boolean_operators {
     public:
       typedef ikos::linear_constraint<Number, VariableName> lin_cst_t;       
       
       virtual void backward_assign_bool_cst(VariableName lhs, lin_cst_t rhs,
					     NumAbsDom invariant){}
       virtual void backward_assign_bool_var(VariableName lhs, VariableName rhs, bool is_not_rhs,
					     NumAbsDom invariant) {}
       virtual void backward_apply_binary_bool(bool_operation_t op,
					       VariableName x,VariableName y,VariableName z,
					       NumAbsDom invariant) {}
       virtual ~backward_boolean_operators () {};
     };
     
   } // namespace domains  
} // namespace crab

namespace crab {
  template<>
  inline boost::optional<domains::bool_operation_t> 
  conv_op (bool_binary_operation_t op) {     
    switch (op) {
    case BINOP_BAND: return domains::OP_BAND;
    case BINOP_BOR:  return domains::OP_BOR;
    case BINOP_BXOR: return domains::OP_BXOR;
    default: return boost::optional<domains::bool_operation_t> ();
    }
  }
}

#endif 
