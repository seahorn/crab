/*******************************************************************************
 * Generic API for division and remainder operations.
 ******************************************************************************/

#ifndef IKOS_DIVISION_OPERATORS_API_HPP
#define IKOS_DIVISION_OPERATORS_API_HPP

namespace ikos {

  typedef enum { 
    OP_SDIV, 
    OP_UDIV, 
    OP_SREM, 
    OP_UREM
  } div_operation_t;

  template< typename Number, typename VariableName >
  class division_operators {

  public:
    virtual void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) = 0;
    virtual void apply(div_operation_t op, VariableName x, VariableName y, Number z) = 0;
    virtual ~division_operators() { }

  }; // class division_operators
  
} // namespace ikos

#endif // IKOS_DIVISION_OPERATORS_API_HPP
