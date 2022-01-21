#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace cfg {

// To group together statements over integers or real
typedef enum {
  BINOP_ADD,
  BINOP_SUB,
  BINOP_MUL,
  BINOP_SDIV,
  BINOP_UDIV,
  BINOP_SREM,
  BINOP_UREM,
  BINOP_AND,
  BINOP_OR,
  BINOP_XOR,
  BINOP_SHL,
  BINOP_LSHR,
  BINOP_ASHR,
} binary_operation_t;

typedef enum { BINOP_BAND, BINOP_BOR, BINOP_BXOR } bool_binary_operation_t;

typedef enum { CAST_TRUNC, CAST_SEXT, CAST_ZEXT } cast_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, binary_operation_t op) {
  switch (op) {
  case BINOP_ADD:
    o << "+";
    break;
  case BINOP_SUB:
    o << "-";
    break;
  case BINOP_MUL:
    o << "*";
    break;
  case BINOP_SDIV:
    o << "/";
    break;
  case BINOP_UDIV:
    o << "/_u";
    break;
  case BINOP_SREM:
    o << "%";
    break;
  case BINOP_UREM:
    o << "%_u";
    break;
  case BINOP_AND:
    o << "&";
    break;
  case BINOP_OR:
    o << "|";
    break;
  case BINOP_XOR:
    o << "^";
    break;
  case BINOP_SHL:
    o << "<<";
    break;
  case BINOP_LSHR:
    o << ">>_l";
    break;
  case BINOP_ASHR:
    o << ">>_r";
    break;
    // Clang complains switch already covers all values
    // default:
    // CRAB_ERROR("unexpected binary operation ", op);
  }
  return o;
}

inline crab::crab_os &operator<<(crab::crab_os &o, bool_binary_operation_t op) {
  switch (op) {
  case BINOP_BAND:
    o << "&";
    break;
  case BINOP_BOR:
    o << "|";
    break;
  case BINOP_BXOR:
    o << "^";
    break;
    // Clang complains switch already covers all values
    //default:
    //CRAB_ERROR("unexpected boolean binary operation ", op);
  }
  return o;
}

inline crab::crab_os &operator<<(crab::crab_os &o, cast_operation_t op) {
  switch (op) {
  case CAST_TRUNC:
    o << "trunc";
    break;
  case CAST_SEXT:
    o << "sext";
    break;
  case CAST_ZEXT:
    o << "zext";
    break;
    // Clang complains switch already covers all values
    //default:
    //CRAB_ERROR("unexpected cast operation", op);
  }
  return o;
}

} // end namespace cfg
} // end namespace crab
