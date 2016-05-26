/*******************************************************************************
 *
 * Generic API for numerical domains.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#ifndef IKOS_NUMERICAL_DOMAINS_API_HPP
#define IKOS_NUMERICAL_DOMAINS_API_HPP

#include <crab/common/types.hpp>
#include <crab/domains/linear_constraints.hpp>

namespace ikos {

  // Enumeration type for basic arithmetic operations
  typedef enum {
    OP_ADDITION,
    OP_SUBTRACTION,
    OP_MULTIPLICATION,
    OP_DIVISION
  } operation_t;

  inline crab::crab_os& operator<<(crab::crab_os&o, operation_t op) {
    switch (op) {
      case OP_ADDITION: o << "+"; break;
      case OP_SUBTRACTION: o << "-"; break;
      case OP_MULTIPLICATION: o << "*"; break;
    default: o << "/"; break;
    }
    return o;
  }
  
  template< typename Number, typename VariableName >
  class numerical_domain {
    
  public:
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    typedef variable< Number, VariableName > variable_t;
    typedef Number number_t;
    typedef VariableName varname_t;
    
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

namespace crab {
  using namespace ikos;

  template<>
  inline boost::optional<operation_t> 
  conv_op (binary_operation_t op) {     
    switch (op) {
      case BINOP_ADD: return OP_ADDITION;
      case BINOP_SUB: return OP_SUBTRACTION;
      case BINOP_MUL: return OP_MULTIPLICATION;
      case BINOP_SDIV: return OP_DIVISION;
      default: return boost::optional<operation_t> ();
    }
  }
}

#endif // IKOS_NUMERICAL_DOMAINS_API_HPP
