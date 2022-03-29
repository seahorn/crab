#pragma once

#include <crab/domains/sign.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

template <typename Number>
sign<Number> sign<Number>::shiftOp(const sign<Number> &o) const {
  // The shift operation is shl, lshr, or ashr
  // zero is special
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else {
    if (equal_zero() || o.equal_zero()) {
      return *this;
    } else {
      return top();
    }
  }
}

template <typename Number>
sign<Number> sign<Number>::defaultOp(const sign<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number> sign<Number>::sign(sign_interval s) : m_sign(s) {}

template <typename Number>
sign<Number>::sign(bool is_bottom)
    : m_sign(is_bottom ? sign_interval::BOT : sign_interval::TOP) {}

template <typename Number>
sign<Number>::sign(Number c) : m_sign(sign_interval::TOP) {
  Number zero(get_zero());
  if (c == zero) {
    m_sign = sign_interval::EQZ;
  } else if (c < zero) {
    m_sign = sign_interval::LTZ;
  } else /* (c > zero) */ {
    m_sign = sign_interval::GTZ;
  }
}

/* MODIFY HERE if not integers */
template <typename Number>
sign<Number>
sign<Number>::from_interval(const ikos::interval<Number> &i) const {
  Number zero(get_zero());
  Number plus_one(get_plus_one());
  Number minus_one(get_minus_one());
  if (i.is_bottom()) {
    return sign<Number>::bottom();
  } else if (i.is_top()) {
    return sign<Number>::top();
  } else if (i <= ikos::interval<Number>(zero, zero)) {
    return sign<Number>(sign_interval::EQZ);
  } else if (i <=
             ikos::interval<Number>(bound_t::minus_infinity(), minus_one)) {
    return sign<Number>(sign_interval::LTZ);
  } else if (i <= ikos::interval<Number>(bound_t::minus_infinity(), zero)) {
    return sign<Number>(sign_interval::LEZ);
  } else if (i <= ikos::interval<Number>(plus_one, bound_t::plus_infinity())) {
    return sign<Number>(sign_interval::GTZ);
  } else if (i <= ikos::interval<Number>(zero, bound_t::plus_infinity())) {
    return sign<Number>(sign_interval::GEZ);
  } else {
    // unreachable
    return sign<Number>::top();
  }
}

/* MODIFY HERE if not integers */
template <typename Number>
ikos::interval<Number> sign<Number>::to_interval() const {
  Number zero(get_zero());
  Number plus_one(get_plus_one());
  Number minus_one(get_minus_one());
  if (is_bottom()) {
    return ikos::interval<Number>::bottom();
  } else if (is_top()) {
    return ikos::interval<Number>::top();
  } else if (equal_zero()) {
    return ikos::interval<Number>(zero);
  } else if (less_than_zero()) {
    return ikos::interval<Number>(bound_t::minus_infinity(), minus_one);
  } else if (greater_than_zero()) {
    return ikos::interval<Number>(plus_one, bound_t::plus_infinity());
  } else if (less_or_equal_than_zero()) {
    return ikos::interval<Number>(bound_t::minus_infinity(), zero);
  } else if (greater_or_equal_than_zero()) {
    return ikos::interval<Number>(zero, bound_t::plus_infinity());
  } else /* not_equal_zero*/ {
    // we cannot express not_equal_zero in an interval
    return ikos::interval<Number>::top();
  }
}

template <typename Number> sign<Number> sign<Number>::bottom() {
  return sign(true);
}

template <typename Number> sign<Number> sign<Number>::top() {
  return sign(false);
}

template <typename Number> sign<Number> sign<Number>::mk_equal_zero() {
  return sign<Number>(sign_interval::EQZ);
}

template <typename Number> sign<Number> sign<Number>::mk_less_than_zero() {
  return sign<Number>(sign_interval::LTZ);
}

template <typename Number> sign<Number> sign<Number>::mk_greater_than_zero() {
  return sign<Number>(sign_interval::GTZ);
}

template <typename Number>
sign<Number> sign<Number>::mk_less_or_equal_than_zero() {
  return sign<Number>(sign_interval::LEZ);
}

template <typename Number>
sign<Number> sign<Number>::mk_greater_or_equal_than_zero() {
  return sign<Number>(sign_interval::GEZ);
}

template <typename Number> sign<Number> sign<Number>::mk_not_equal_zero() {
  return sign<Number>(sign_interval::NEZ);
}

template <typename Number> bool sign<Number>::is_bottom() const {
  return m_sign == sign_interval::BOT;
}

template <typename Number> bool sign<Number>::is_top() const {
  return m_sign == sign_interval::TOP;
}

template <typename Number> bool sign<Number>::equal_zero() const {
  return m_sign == sign_interval::EQZ;
}

template <typename Number> bool sign<Number>::less_than_zero() const {
  return m_sign == sign_interval::LTZ;
}

template <typename Number> bool sign<Number>::greater_than_zero() const {
  return m_sign == sign_interval::GTZ;
}

template <typename Number> bool sign<Number>::less_or_equal_than_zero() const {
  return m_sign == sign_interval::LEZ;
}

template <typename Number>
bool sign<Number>::greater_or_equal_than_zero() const {
  return m_sign == sign_interval::GEZ;
}

template <typename Number> bool sign<Number>::not_equal_zero() const {
  return m_sign == sign_interval::NEZ;
}

template <typename Number>
bool sign<Number>::operator<=(const sign<Number> &o) const {
  if (is_bottom() || o.is_top()) {
    return true;
  } else if (o.is_bottom() || is_top()) {
    return false;
  } else {
    // operands are not either top or bottom
    if (m_sign == sign_interval::LTZ) {
      return o.m_sign == m_sign || o.m_sign == sign_interval::LEZ ||
             o.m_sign == sign_interval::NEZ;
    } else if (m_sign == sign_interval::GTZ) {
      return o.m_sign == m_sign || o.m_sign == sign_interval::GEZ ||
             o.m_sign == sign_interval::NEZ;
    } else if (m_sign == sign_interval::EQZ) {
      return o.m_sign == m_sign || o.m_sign == sign_interval::LEZ ||
             o.m_sign == sign_interval::GEZ;
    } else {
      /* m_sign == sign_interval::LEZ || m_sign == sign_interval::GEZ m_sign
       * == sign_interval::NEZ */
      return m_sign == o.m_sign;
    }
  }
}

template <typename Number>
bool sign<Number>::operator==(const sign<Number> &o) const {
  return (m_sign == o.m_sign);
}

template <typename Number>
sign<Number> sign<Number>::operator|(const sign<Number> &o) const {
  if (is_bottom() || o.is_top()) {
    return o;
  } else if (is_top() || o.is_bottom()) {
    return *this;
  } else {
    // operands are not either top or bottom
    if (m_sign == sign_interval::LTZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
        return *this;
      case sign_interval::GTZ:
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::NEZ);
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::LEZ);
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::TOP);
      default:
        CRAB_ERROR("sign::operator| unreachable 1");
      }
    } else if (m_sign == sign_interval::GTZ) {
      switch (o.m_sign) {
      case sign_interval::GTZ:
        return *this;
      case sign_interval::LTZ:
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::NEZ);
      case sign_interval::EQZ:
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::GEZ);
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::TOP);
      default:
        CRAB_ERROR("sign::operator| unreachable 2");
      }
    } else if (m_sign == sign_interval::EQZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::LEZ);
      case sign_interval::GTZ:
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::GEZ);
      case sign_interval::EQZ:
        return *this;
      case sign_interval::NEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator| unreachable 3");
      }
    } else if (m_sign == sign_interval::LEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return *this;
      case sign_interval::GTZ:
      case sign_interval::NEZ:
      case sign_interval::GEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator| unreachable 4");
      }
    } else if (m_sign == sign_interval::NEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::GTZ:
      case sign_interval::NEZ:
        return *this;
      case sign_interval::EQZ:
      case sign_interval::LEZ:
      case sign_interval::GEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator| unreachable 5");
      }
    } else if (m_sign == sign_interval::GEZ) {
      switch (o.m_sign) {
      case sign_interval::EQZ:
      case sign_interval::GTZ:
      case sign_interval::GEZ:
        return *this;
      case sign_interval::LTZ:
      case sign_interval::LEZ:
      case sign_interval::NEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator| unreachable 6");
      }
    }
  }
  CRAB_ERROR("sign::operator| unreachable 7");
}

template <typename Number>
sign<Number> sign<Number>::operator&(const sign<Number> &o) const {
  if (is_bottom() || o.is_top())
    return *this;
  else if (is_top() || o.is_bottom()) {
    return o;
  } else {
    // operands are not either top or bottom
    if (m_sign == sign_interval::LTZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
        return *this;
      case sign_interval::GTZ:
      case sign_interval::EQZ:
      case sign_interval::GEZ:
        return bottom();
      case sign_interval::LEZ:
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::LTZ);
      default:
        CRAB_ERROR("sign::operator& unreachable 1");
      }
    } else if (m_sign == sign_interval::GTZ) {
      switch (o.m_sign) {
      case sign_interval::GTZ:
        return *this;
      case sign_interval::LTZ:
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return bottom();
      case sign_interval::GEZ:
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::GTZ);
      default:
        CRAB_ERROR("sign::operator& unreachable 2");
      }
    } else if (m_sign == sign_interval::EQZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::GTZ:
      case sign_interval::NEZ:
        return bottom();
      case sign_interval::GEZ:
      case sign_interval::LEZ:
      case sign_interval::EQZ:
        return *this;
      default:
        CRAB_ERROR("sign::operator& unreachable 3");
      }
    } else if (m_sign == sign_interval::LEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return o;
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::LTZ);
      case sign_interval::GTZ:
        return bottom();
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::EQZ);
      default:
        CRAB_ERROR("sign::operator& unreachable 4");
      }
    } else if (m_sign == sign_interval::GEZ) {
      switch (o.m_sign) {
      case sign_interval::EQZ:
      case sign_interval::GTZ:
      case sign_interval::GEZ:
        return o;
      case sign_interval::LTZ:
        return bottom();
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::EQZ);
      case sign_interval::NEZ:
        return sign<Number>(sign_interval::GTZ);
      default:
        CRAB_ERROR("sign::operator& unreachable 5");
      }
    } else if (m_sign == sign_interval::NEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::GTZ:
      case sign_interval::NEZ:
        return o;
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::LTZ);
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::GTZ);
      case sign_interval::EQZ:
        return bottom();
      default:
        CRAB_ERROR("sign::operator& unreachable 6");
      }
    }
  }
  CRAB_ERROR("sign::operator& unreachable 7");
}

// addition
template <typename Number>
sign<Number> sign<Number>::operator+(const sign<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else if (is_top() || o.is_top()) {
    return top();
  } else {
    // operands are not either top or bottom
    if (m_sign == sign_interval::LTZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return *this;
      case sign_interval::GTZ:
      case sign_interval::NEZ:
      case sign_interval::GEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator+ unreachable 1");
      }
    } else if (m_sign == sign_interval::GTZ) {
      switch (o.m_sign) {
      case sign_interval::GTZ:
      case sign_interval::EQZ:
      case sign_interval::GEZ:
        return *this;
      case sign_interval::LTZ:
      case sign_interval::LEZ:
      case sign_interval::NEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator+ unreachable 2");
      }
    } else if (m_sign == sign_interval::LEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::EQZ:
      case sign_interval::LEZ:
        return *this;
      case sign_interval::GTZ:
      case sign_interval::NEZ:
      case sign_interval::GEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator+ unreachable 2");
      }
    } else if (m_sign == sign_interval::GEZ) {
      switch (o.m_sign) {
      case sign_interval::GTZ:
      case sign_interval::GEZ:
      case sign_interval::EQZ:
        return *this;
      case sign_interval::LTZ:
      case sign_interval::LEZ:
      case sign_interval::NEZ:
        return top();
      default:
        CRAB_ERROR("sign::operator+ unreachable 2");
      }
    } else if (m_sign == sign_interval::EQZ) {
      return o;
    } else /* (m_sign == sign_interval::NEZ)*/ {
      return top();
    }
  }
}

// subtraction
template <typename Number>
sign<Number> sign<Number>::operator-(const sign<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else if (is_top() || o.is_top()) {
    return top();
  } else {
    switch (o.m_sign) {
    case sign_interval::GTZ:
      return *this + sign<Number>(sign_interval::LTZ);
    case sign_interval::GEZ:
      return *this + sign<Number>(sign_interval::LEZ);
    case sign_interval::EQZ:
      return *this;
    case sign_interval::LTZ:
      return *this + sign<Number>(sign_interval::GTZ);
    case sign_interval::LEZ:
      return *this + sign<Number>(sign_interval::GEZ);
    default: /*sign_interval::NEZ*/
      return top();
    }
  }
}

// multiplication
template <typename Number>
sign<Number> sign<Number>::operator*(const sign<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else if (equal_zero() || o.equal_zero()) {
    return sign<Number>(sign_interval::EQZ);
  } else if (is_top() || o.is_top()) {
    return top();
  } else if (not_equal_zero() || o.not_equal_zero()) {
    return top();
  } else {
    // operands are not either top, bottom, zero, or non-zero
    if (m_sign == sign_interval::LTZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
        return sign<Number>(sign_interval::GTZ);
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::GEZ);
      case sign_interval::GTZ:
        return sign<Number>(sign_interval::LTZ);
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::LEZ);
      default:
        CRAB_ERROR("sign::operator* unreachable 1");
      }
    } else if (m_sign == sign_interval::GTZ) {
      return o;
    } else if (m_sign == sign_interval::LEZ) {
      switch (o.m_sign) {
      case sign_interval::LTZ:
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::GEZ);
      case sign_interval::GTZ:
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::LEZ);
      default:
        CRAB_ERROR("sign::operator* unreachable 2");
      }
    } else if (m_sign == sign_interval::GEZ) {
      switch (o.m_sign) {
      case sign_interval::GTZ:
      case sign_interval::GEZ:
        return sign<Number>(sign_interval::GEZ);
      case sign_interval::LTZ:
      case sign_interval::LEZ:
        return sign<Number>(sign_interval::LEZ);
      default:
        CRAB_ERROR("sign::operator* unreachable 3");
      }
    }
  }
  CRAB_ERROR("sign::operator* unreachable 4");
}

// signed division
template <typename Number>
sign<Number> sign<Number>::operator/(const sign<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else if (o.equal_zero()) {
    return bottom();
  } else if (equal_zero()) {
    return *this;
  } else if (is_top() || o.is_top()) {
    return top();
  } else if (not_equal_zero() || o.not_equal_zero()) {
    return top();
  } else {
    // Once we exclude top, bottom, zero, and non-zero
    // signed division is like multiplication
    return (*this) * o;
  }
}

// division and remainder operations

template <typename Number>
sign<Number> sign<Number>::UDiv(const sign<Number> &o) const {
  return defaultOp(o);
}

template <typename Number>
sign<Number> sign<Number>::SRem(const sign<Number> &o) const {
  return defaultOp(o);
}

template <typename Number>
sign<Number> sign<Number>::URem(const sign<Number> &o) const {
  return defaultOp(o);
}

// bitwise operations
template <typename Number>
sign<Number> sign<Number>::And(const sign<Number> &o) const {
  // zero is special
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else {
    if (equal_zero() || o.equal_zero()) {
      return sign<Number>(sign_interval::EQZ);
    } else {
      return top();
    }
  }
}

template <typename Number>
sign<Number> sign<Number>::Or(const sign<Number> &o) const {
  // zero is special
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  } else {
    if (equal_zero()) {
      return o;
    } else if (o.equal_zero()) {
      return *this;
    } else {
      return top();
    }
  }
}

template <typename Number>
sign<Number> sign<Number>::Xor(const sign<Number> &o) const {
  return Or(o);
}

template <typename Number>
sign<Number> sign<Number>::Shl(const sign<Number> &o) const {
  return shiftOp(o);
}

template <typename Number>
sign<Number> sign<Number>::LShr(const sign<Number> &o) const {
  return shiftOp(o);
}

template <typename Number>
sign<Number> sign<Number>::AShr(const sign<Number> &o) const {
  return shiftOp(o);
}

template <typename Number> void sign<Number>::write(crab::crab_os &o) const {
  switch (m_sign) {
  case sign_interval::BOT:
    o << "_|_";
    return;
  case sign_interval::LTZ:
    o << "[-oo,-1]";
    return;
  case sign_interval::GTZ:
    o << "[1,+oo]";
    return;
  case sign_interval::EQZ:
    o << "[0,0]";
    return;
  case sign_interval::NEZ:
    o << "[-oo,-1] U [1,+oo]";
    return;
  case sign_interval::GEZ:
    o << "[0,+oo]";
    return;
  case sign_interval::LEZ:
    o << "[-oo,0]";
    return;
  default: /*sign_interval::TOP*/
    o << "[-oo,+oo]";
  }
}

} // namespace domains
} // namespace crab
