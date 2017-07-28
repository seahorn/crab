/*******************************************************************************
 * Generic API for cast operations
 ******************************************************************************/

#ifndef CAST_OPERATORS_API_HPP
#define CAST_OPERATORS_API_HPP

#include <crab/common/types.hpp>
#include <crab/common/bignums.hpp>

namespace crab {
  namespace domains {

    typedef enum  { 
      OP_TRUNC, 
      OP_SEXT, 
      OP_ZEXT 
    } int_conv_operation_t;

    inline crab::crab_os& operator<<(crab::crab_os&o, int_conv_operation_t op) {
      switch (op) {
      case OP_TRUNC: o << "trunc"; break;
      case OP_SEXT : o << "sext"; break;
      case OP_ZEXT : o << "zext"; break;
      default:
	CRAB_ERROR("unreachable");
      }
      return o;
    }

    template<typename VariableName >
    class int_cast_operators {
    public:
      virtual void apply(int_conv_operation_t op,
			 VariableName dst, unsigned dst_width,
			 VariableName src, unsigned src_width)  = 0;

      virtual ~int_cast_operators() { }
    }; 
    
  } // end domains

  template<>
  inline boost::optional<domains::int_conv_operation_t> 
  conv_op (cast_operation_t op) {     
    switch (op) {
    case CAST_TRUNC: return domains::OP_TRUNC;
    case CAST_SEXT:  return domains::OP_SEXT;
    case CAST_ZEXT: return domains::OP_ZEXT;
    default: return boost::optional<domains::int_conv_operation_t> ();
    }
  }
  
} // end namespace crab
#endif 
