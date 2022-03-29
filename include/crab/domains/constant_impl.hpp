#pragma once

#include <crab/domains/constant.hpp>

namespace crab {
namespace domains {

template <typename Number>
constant<Number>::constant(bool is_bottom)
    : m_constant(boost::none), m_is_bottom(is_bottom) {}

template <typename Number>
constant<Number>::constant(Number c) : m_constant(c), m_is_bottom(false) {}

template <typename Number> constant<Number> constant<Number>::bottom() {
  return constant(true);
}

template <typename Number> constant<Number> constant<Number>::top() {
  return constant(false);
}

template <typename Number> constant<Number> constant<Number>::zero() {
  return constant(Number(0));
}

template <typename Number> bool constant<Number>::is_bottom() const {
  return m_is_bottom;
}

template <typename Number> bool constant<Number>::is_top() const {
  return (!is_bottom() && !m_constant);
}

template <typename Number> bool constant<Number>::is_constant() const {
  return m_constant != boost::none;
}

template <typename Number> Number constant<Number>::get_constant() const {
  assert(is_constant());
  return *m_constant;
}

template <typename Number>
bool constant<Number>::operator<=(const constant<Number> &o) const {
  if (is_bottom() || o.is_top()) {
    return true;
  } else if (o.is_bottom() || is_top()) {
    return false;
  } else {
    assert(is_constant());
    assert(o.is_constant());
    return get_constant() == o.get_constant();
  }
}

template <typename Number>
bool constant<Number>::operator==(const constant<Number> &o) const {
  return (m_is_bottom == o.m_is_bottom && m_constant == o.m_constant);
}

template <typename Number>
constant<Number> constant<Number>::operator|(const constant<Number> &o) const {
  if (is_bottom() || o.is_top())
    return o;
  else if (is_top() || o.is_bottom())
    return *this;
  else {
    assert(is_constant());
    assert(o.is_constant());
    if (get_constant() == o.get_constant()) {
      return *this;
    } else {
      return constant<Number>::top();
    }
  }
}

template <typename Number>
constant<Number> constant<Number>::operator||(const constant<Number> &o) const {
  return *this | o;
}

template <typename Number>
constant<Number> constant<Number>::operator&(const constant<Number> &o) const {
  if (is_bottom() || o.is_top())
    return *this;
  else if (is_top() || o.is_bottom()) {
    return o;
  } else {
    assert(is_constant());
    assert(o.is_constant());
    if (get_constant() == o.get_constant()) {
      return *this;
    } else {
      return constant<Number>::bottom();
    }
  }
}

template <typename Number>
constant<Number> constant<Number>::operator&&(const constant<Number> &o) const {
  return *this & o;
}

template <typename Number>
constant<Number> constant<Number>::Add(const constant<Number> &o) const {
  if (is_constant() && o.is_constant()) {
    return constant<Number>(get_constant() + o.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
constant<Number> constant<Number>::Sub(const constant<Number> &o) const {
  if (is_constant() && o.is_constant()) {
    return constant<Number>(get_constant() - o.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
constant<Number> constant<Number>::Mul(const constant<Number> &o) const {
  if (is_constant() && o.is_constant()) {
    return constant<Number>(get_constant() * o.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
constant<Number> constant<Number>::SDiv(const constant<Number> &o) const {
  if (o.is_constant() && o.get_constant() == Number(0)) {
    return constant<Number>::bottom();
  }
  if (is_constant() && o.is_constant()) {
    return constant<Number>(get_constant() / o.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
constant<Number> constant<Number>::SRem(const constant<Number> &o) const {
  if (o.is_constant() && o.get_constant() == Number(0)) {
    return constant<Number>::bottom();
  }
  return constant<Number>::top();
}

template <typename Number>
constant<Number> constant<Number>::UDiv(const constant<Number> &o) const {
  if (o.is_constant() && o.get_constant() == Number(0)) {
    return constant<Number>::bottom();
  }
  return constant<Number>::top();
}

template <typename Number>
constant<Number> constant<Number>::URem(const constant<Number> &o) const {
  if (o.is_constant() && o.get_constant() == Number(0)) {
    return constant<Number>::bottom();
  }
  return constant<Number>::top();
}

template <typename Number>
constant<Number> constant<Number>::BitwiseAnd(const constant<Number> &o) const {
  return constant<Number>::top();
}
template <typename Number>
constant<Number> constant<Number>::BitwiseOr(const constant<Number> &o) const {
  return constant<Number>::top();
}
template <typename Number>
constant<Number> constant<Number>::BitwiseXor(const constant<Number> &o) const {
  return constant<Number>::top();
}
template <typename Number>
constant<Number> constant<Number>::BitwiseShl(const constant<Number> &o) const {
  return constant<Number>::top();
}
template <typename Number>
constant<Number>
constant<Number>::BitwiseLShr(const constant<Number> &o) const {
  return constant<Number>::top();
}
template <typename Number>
constant<Number>
constant<Number>::BitwiseAShr(const constant<Number> &o) const {
  return constant<Number>::top();
}

template <typename Number>
void constant<Number>::write(crab::crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else if (is_top()) {
    o << "top";
  } else {
    assert(is_constant());
    o << get_constant();
  }
}

} // namespace domains
} // namespace crab
