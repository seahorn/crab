#pragma once

#include <crab/support/os.hpp>

/* Types for abstract domain operations */

namespace crab {
namespace domains {

// Enumeration type for basic arithmetic operations
// Do not modify the order.
typedef enum {
  OP_ADDITION,
  OP_SUBTRACTION,
  OP_MULTIPLICATION,
  OP_SDIV,
  OP_UDIV,
  OP_SREM,
  OP_UREM
} arith_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, arith_operation_t op) {
  switch (op) {
  case OP_ADDITION:
    o << "+";
    break;
  case OP_SUBTRACTION:
    o << "-";
    break;
  case OP_MULTIPLICATION:
    o << "*";
    break;
  case OP_SDIV:
    o << "/";
    break;
  case OP_UDIV:
    o << "/_u";
    break;
  case OP_SREM:
    o << "%";
    break;
  default:
    o << "%_u";
    break;
  }
  return o;
}

// Enumeration type for bitwise operations
typedef enum {
  OP_AND,
  OP_OR,
  OP_XOR,
  OP_SHL,
  OP_LSHR,
  OP_ASHR
} bitwise_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, bitwise_operation_t op) {
  switch (op) {
  case OP_AND:
    o << "&";
    break;
  case OP_OR:
    o << "|";
    break;
  case OP_XOR:
    o << "^";
    break;
  case OP_SHL:
    o << "<<";
    break;
  case OP_LSHR:
    o << ">>_l";
    break;
  default:
    o << ">>_a";
    break;
  }
  return o;
}

// Enumeration type for cast operations
typedef enum { OP_TRUNC, OP_SEXT, OP_ZEXT } int_conv_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, int_conv_operation_t op) {
  switch (op) {
  case OP_TRUNC:
    o << "trunc";
    break;
  case OP_SEXT:
    o << "sext";
    break;
  default: /*OP_ZEXT*/
    o << "zext";
    break;
  }
  return o;
}

// Enumeration types for boolean operations
typedef enum { OP_BAND, OP_BOR, OP_BXOR } bool_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, bool_operation_t op) {
  switch (op) {
  case OP_BAND:
    o << "&";
    break;
  case OP_BOR:
    o << "|";
    break;
  default: /*OP_BXOR*/
    o << "^";
    break;
  }
  return o;
}
} // end namespace domains
} // end namespace crab
