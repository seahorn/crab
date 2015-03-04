/*******************************************************************************
 * Generic API for numerical domains.
 ******************************************************************************/


#ifndef IKOS_NUMERICAL_DOMAINS_API_HPP
#define IKOS_NUMERICAL_DOMAINS_API_HPP

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp>

namespace ikos {
  
  template< typename Number, typename VariableName >
  class numerical_domain {
    
  public:
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    
  public:
    virtual void apply(operation_t op, VariableName x, VariableName y, VariableName z) = 0; // x = y op z

    virtual void apply(operation_t op, VariableName x, VariableName y, Number k) = 0; // x = y op k

    virtual void assign(VariableName x, linear_expression_t e) = 0; // x = e

    virtual void operator+=(linear_constraint_system_t csts) = 0;

    virtual void operator-=(VariableName v) = 0;

    void operator+=(linear_constraint_t cst) {
      linear_constraint_system_t csts(cst);
      operator+=(csts);
    }

    virtual ~numerical_domain() { }
    
  }; // numerical_domain
  
} // namespace ikos

#endif // IKOS_NUMERICAL_DOMAINS_API_HPP
