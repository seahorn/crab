#include <crab/domains/term/term_operators.hpp>

namespace crab {
namespace domains {
namespace term {

term_operator_t::term_operator_t(uint32_t value, boost::string_ref name)
  : m_value(value), m_name(name) {
  if (m_value < first_nonreserved_value()) {
    CRAB_ERROR("First " + std::to_string(first_nonreserved_value()) +
	       "numbers are reserved in term_operator_t");
  }  
}

term::term_operator_t conv2termop(arith_operation_t op) {
  switch (op) {
  case OP_ADDITION:
    return term::TERM_OP_ADD;
  case OP_SUBTRACTION:
    return term::TERM_OP_SUB;
  case OP_MULTIPLICATION:
    return term::TERM_OP_MUL;
  case OP_SDIV:
    return term::TERM_OP_SDIV;
  case OP_UDIV:
    return term::TERM_OP_UDIV;
  case OP_SREM:
    return term::TERM_OP_SREM;
  default:
    return term::TERM_OP_UREM;
  }
}

boost::optional<arith_operation_t> conv2arith(term::term_operator_t op) {
  switch (op) {
  case term::TERM_OP_ADD:
    return OP_ADDITION;
  case term::TERM_OP_SUB:
    return OP_SUBTRACTION;
  case term::TERM_OP_MUL:
    return OP_MULTIPLICATION;
  case term::TERM_OP_SDIV:
    return OP_SDIV;
  case term::TERM_OP_UDIV:
    return OP_UDIV;
  case term::TERM_OP_SREM:
    return OP_SREM;
  case term::TERM_OP_UREM:
    return OP_UREM;
  default:
    return boost::optional<arith_operation_t>();
  }
}

term::term_operator_t conv2termop(bitwise_operation_t op) {
  switch (op) {
  case OP_AND:
    return term::TERM_OP_AND;
  case OP_OR:
    return term::TERM_OP_OR;
  case OP_XOR:
    return term::TERM_OP_XOR;
  case OP_SHL:
    return term::TERM_OP_SHL;
  case OP_LSHR:
    return term::TERM_OP_LSHR;
  default:
    return term::TERM_OP_ASHR;
  }
}

boost::optional<bitwise_operation_t> conv2bitwise(term::term_operator_t op) {
  switch (op) {
  case term::TERM_OP_AND:
    return OP_AND;
  case term::TERM_OP_OR:
    return OP_OR;
  case term::TERM_OP_XOR:
    return OP_XOR;
  case term::TERM_OP_SHL:
    return OP_SHL;
  case term::TERM_OP_LSHR:
    return OP_LSHR;
  case term::TERM_OP_ASHR:
    return OP_ASHR;
  default:
    return boost::optional<bitwise_operation_t>();
  }
}

term::term_operator_t conv2termop(bool_operation_t op) {
  switch (op) {
  case OP_BAND:
    return term::TERM_OP_BAND;
  case OP_BOR:
    return term::TERM_OP_BOR;
  default:
    return term::TERM_OP_BXOR;
  }
}

boost::optional<bool_operation_t> conv2bool(term::term_operator_t op) {
  switch (op) {
  case term::TERM_OP_BAND:
    return OP_BAND;
  case term::TERM_OP_BOR:
    return OP_BOR;
  case term::TERM_OP_BXOR:
    return OP_BXOR;
  default:
    return boost::optional<bool_operation_t>();
  }
}

} // namespace term
} // end namespace domains
} // end namespace crab
