/*******************************************************************************
 * Generic API for bitwise operations.
 ******************************************************************************/

#ifndef IKOS_BITWISE_OPERATORS_API_HPP
#define IKOS_BITWISE_OPERATORS_API_HPP

namespace ikos {

  typedef enum  { 
    OP_TRUNC, 
    OP_SEXT, 
    OP_ZEXT 
  } conv_operation_t;

  typedef enum  { 
    OP_AND, 
    OP_OR, 
    OP_XOR, 
    OP_SHL, 
    OP_LSHR, 
    OP_ASHR
  } bitwise_operation_t;

  template< typename Number, typename VariableName >
  class bitwise_operators {
   public:
    virtual void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) = 0;
    virtual void apply(conv_operation_t op, VariableName x, Number y, unsigned width) = 0;
    virtual void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) = 0;
    virtual void apply(bitwise_operation_t op, VariableName x, VariableName y, Number z) = 0;
    virtual ~bitwise_operators() { }

  }; // class bitwise_operators
  
} // namespace ikos

#endif // IKOS_BITWISE_OPERATORS_API_HPP
