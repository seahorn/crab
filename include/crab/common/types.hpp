#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <functional>
#include <iosfwd>
#include <memory>

/* Basic type definitions */

namespace ikos {
// Numerical type for indexed objects
typedef uint64_t index_t;
} // end namespace ikos

namespace crab {
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
  BINOP_FUNCTION
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
  case BINOP_FUNCTION:
    o << "uf";
    break;
  default:
    CRAB_ERROR("unexpected binary operation ", op);
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
  default:
    CRAB_ERROR("unexpected boolean binary operation ", op);
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
  default:
    CRAB_ERROR("unexpected cast operation", op);
  }
  return o;
}

template <typename T> inline boost::optional<T> conv_op(binary_operation_t op);

template <typename T>
inline boost::optional<T> conv_op(bool_binary_operation_t op);

template <typename T> inline boost::optional<T> conv_op(cast_operation_t op);

} // end namespace crab

namespace crab {
namespace domains {
using namespace ikos;
}
} // namespace crab


