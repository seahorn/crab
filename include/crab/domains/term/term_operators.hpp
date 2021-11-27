#pragma once

#include <crab/domains/abstract_domain_operators.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
namespace term {

typedef enum {
  TERM_OP_ADD,
  TERM_OP_SUB,
  TERM_OP_MUL,
  TERM_OP_SDIV,
  TERM_OP_UDIV,
  TERM_OP_SREM,
  TERM_OP_UREM,
  TERM_OP_AND,
  TERM_OP_OR,
  TERM_OP_XOR,
  TERM_OP_SHL,
  TERM_OP_LSHR,
  TERM_OP_ASHR,
  TERM_OP_FUNCTION
} term_operator_t;

inline crab::crab_os &operator<<(crab::crab_os &o, term_operator_t op) {
  switch (op) {
  case TERM_OP_ADD:
    o << "+";
    break;
  case TERM_OP_SUB:
    o << "-";
    break;
  case TERM_OP_MUL:
    o << "*";
    break;
  case TERM_OP_SDIV:
    o << "/";
    break;
  case TERM_OP_UDIV:
    o << "/_u";
    break;
  case TERM_OP_SREM:
    o << "%";
    break;
  case TERM_OP_UREM:
    o << "%_u";
    break;
  case TERM_OP_AND:
    o << "&";
    break;
  case TERM_OP_OR:
    o << "|";
    break;
  case TERM_OP_XOR:
    o << "^";
    break;
  case TERM_OP_SHL:
    o << "<<";
    break;
  case TERM_OP_LSHR:
    o << ">>_l";
    break;
  case TERM_OP_ASHR:
    o << ">>_r";
    break;
  case TERM_OP_FUNCTION:
    o << "uf";
    break;
  default:
    CRAB_ERROR("unexpected binary operation ", op);
  }
  return o;
}

/* Convert between Crab operators and term domain uninterpreted functors */
term_operator_t conv2termop(arith_operation_t op);
term_operator_t conv2termop(bitwise_operation_t op);
boost::optional<arith_operation_t> conv2arith(term_operator_t op);
boost::optional<bitwise_operation_t> conv2bitwise(term_operator_t op);
} // end namespace term
} // end namespace domains
} // end namespace crab
