#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace cfg {

enum binary_operation_t {
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
};

enum bool_binary_operation_t { 
  BINOP_BAND, 
  BINOP_BOR, 
  BINOP_BXOR 
};

enum cast_operation_t { 
  CAST_TRUNC, 
  CAST_SEXT, 
  CAST_ZEXT 
};

/**
 * The `operator<<` overloads below render an operation the way Crab prints
 * statements for humans, with terse symbols. The `op_name` overloads instead
 * give each operation a stable, spelled-out name, for serializers whose output
 * is matched on by other tools (see cfg_to_json.hpp).
 *
 * Keep the two in sync when adding an operation, but keep them *separate*:
 * respelling a display symbol should stay a cosmetic change, not one that
 * breaks consumers of a serialized CFG.
 *
 * None of the switches below has a `default` case, so that adding an operation
 * is caught at compile time by -Wswitch rather than at run time. The trailing
 * CRAB_ERROR is only for a value outside the enumerators, which cannot happen
 * without a cast.
 */

inline crab::crab_os &operator<<(crab::crab_os &o, binary_operation_t op) {
  switch (op) {
  case BINOP_ADD:
    return o << "+";
  case BINOP_SUB:
    return o << "-";
  case BINOP_MUL:
    return o << "*";
  case BINOP_SDIV:
    return o << "/";
  case BINOP_UDIV:
    return o << "/_u";
  case BINOP_SREM:
    return o << "%";
  case BINOP_UREM:
    return o << "%_u";
  case BINOP_AND:
    return o << "&";
  case BINOP_OR:
    return o << "|";
  case BINOP_XOR:
    return o << "^";
  case BINOP_SHL:
    return o << "<<";
  case BINOP_LSHR:
    return o << ">>_l";
  case BINOP_ASHR:
    return o << ">>_r";
  }
  CRAB_ERROR("unexpected binary operation");
}

constexpr const char *op_name(binary_operation_t op) {
  switch (op) {
  case BINOP_ADD:
    return "add";
  case BINOP_SUB:
    return "sub";
  case BINOP_MUL:
    return "mul";
  case BINOP_SDIV:
    return "sdiv";
  case BINOP_UDIV:
    return "udiv";
  case BINOP_SREM:
    return "srem";
  case BINOP_UREM:
    return "urem";
  case BINOP_AND:
    return "and";
  case BINOP_OR:
    return "or";
  case BINOP_XOR:
    return "xor";
  case BINOP_SHL:
    return "shl";
  case BINOP_LSHR:
    return "lshr";
  case BINOP_ASHR:
    return "ashr";
  }
  CRAB_ERROR("unexpected binary operation");
}

inline crab::crab_os &operator<<(crab::crab_os &o, bool_binary_operation_t op) {
  switch (op) {
  case BINOP_BAND:
    return o << "&";
  case BINOP_BOR:
    return o << "|";
  case BINOP_BXOR:
    return o << "^";
  }
  CRAB_ERROR("unexpected boolean binary operation");
}

constexpr const char *op_name(bool_binary_operation_t op) {
  switch (op) {
  case BINOP_BAND:
    return "and";
  case BINOP_BOR:
    return "or";
  case BINOP_BXOR:
    return "xor";
  }
  CRAB_ERROR("unexpected boolean binary operation");
}

inline crab::crab_os &operator<<(crab::crab_os &o, cast_operation_t op) {
  switch (op) {
  case CAST_TRUNC:
    return o << "trunc";
  case CAST_SEXT:
    return o << "sext";
  case CAST_ZEXT:
    return o << "zext";
  }
  CRAB_ERROR("unexpected cast operation");
}

constexpr const char *op_name(cast_operation_t op) {
  switch (op) {
  case CAST_TRUNC:
    return "trunc";
  case CAST_SEXT:
    return "sext";
  case CAST_ZEXT:
    return "zext";
  }
  CRAB_ERROR("unexpected cast operation");
}

} // end namespace cfg
} // end namespace crab
