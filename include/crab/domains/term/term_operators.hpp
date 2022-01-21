#pragma once

#include <crab/domains/abstract_domain_operators.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <boost/utility/string_ref.hpp>

namespace crab {
namespace domains {
namespace term {

class term_operator_t {
  uint32_t m_value;
  boost::string_ref m_name;
public:
  constexpr term_operator_t(uint32_t value)
    : m_value(value), m_name(boost::string_ref()) {}
  term_operator_t(uint32_t value, boost::string_ref name);
  ~term_operator_t() = default;
  constexpr operator uint32_t() const { return m_value; }
  constexpr uint32_t value() const { return m_value; }
  boost::string_ref name() const { return m_name; }
  constexpr static uint32_t first_nonreserved_value() { return 100;}
};

constexpr term_operator_t TERM_OP_ADD(0);
constexpr term_operator_t TERM_OP_SUB(1);
constexpr term_operator_t TERM_OP_MUL(2);
constexpr term_operator_t TERM_OP_SDIV(3);
constexpr term_operator_t TERM_OP_UDIV(4);
constexpr term_operator_t TERM_OP_SREM(5);
constexpr term_operator_t TERM_OP_UREM(6);
constexpr term_operator_t TERM_OP_NOT(7);
constexpr term_operator_t TERM_OP_BAND(8);
constexpr term_operator_t TERM_OP_BOR(9);
constexpr term_operator_t TERM_OP_BXOR(10);
constexpr term_operator_t TERM_OP_AND(11);
constexpr term_operator_t TERM_OP_OR(12);
constexpr term_operator_t TERM_OP_XOR(13);
constexpr term_operator_t TERM_OP_SHL(14);
constexpr term_operator_t TERM_OP_LSHR(15);
constexpr term_operator_t TERM_OP_ASHR(16);
constexpr term_operator_t TERM_OP_FUNCTION(17);  
  
inline crab::crab_os &operator<<(crab::crab_os &o, term_operator_t op) {

  // User-defined operator
  if (const char* name = op.name().data()) {
    o << name;
    return o;
  }

  // Crab predefined operators
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
  case TERM_OP_NOT:
    o << "not";
    break;
  case TERM_OP_BAND:
    o << "band";
    break;
  case TERM_OP_BOR:
    o << "bor";
    break;
  case TERM_OP_BXOR:
    o << "bxor";
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
term_operator_t conv2termop(bool_operation_t op);
boost::optional<arith_operation_t> conv2arith(term_operator_t op);
boost::optional<bitwise_operation_t> conv2bitwise(term_operator_t op);
boost::optional<bool_operation_t> conv2bool(term_operator_t op);
} // end namespace term
} // end namespace domains
} // end namespace crab
