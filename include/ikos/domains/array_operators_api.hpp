/*******************************************************************************
 * API for array operations.
 ******************************************************************************/

#ifndef IKOS_ARRAY_OPERATORS_API_HPP
#define IKOS_ARRAY_OPERATORS_API_HPP

namespace ikos {

  template< typename Number, typename VariableName >
  class array_operators 
  {
   public:

    typedef linear_expression< Number, VariableName > linear_expression_t;

    virtual void load (VariableName lhs, VariableName arr, VariableName idx) = 0;

    virtual void store (VariableName arr_out, VariableName arr_in, 
                        VariableName idx, linear_expression_t val) = 0;

    virtual ~array_operators() { }

  }; // class array_operators
  
} // namespace ikos

#endif // IKOS_ARRAY_OPERATORS_API_HPP
