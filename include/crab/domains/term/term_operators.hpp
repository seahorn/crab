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
  bool m_reserved;
public:
  // used for defining constexpr crab operators
  constexpr term_operator_t(uint32_t value, bool is_reserved = true)
    : m_value(value), m_reserved(is_reserved) {}
  // used for defining new user operators
  static term_operator_t make_operator(uint32_t value);  
  ~term_operator_t() = default;
  constexpr operator uint32_t() const { return m_value; }
  constexpr uint32_t value() const { return m_value; }
  constexpr bool is_reserved() const { return m_reserved;}
  constexpr static uint32_t first_nonreserved_value() { return 100;}  
  friend crab::crab_os &operator<<(crab::crab_os &o, term_operator_t op);
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
