#pragma once

#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace ikos {

template <typename Number>
bound<Number>::bound(bool is_infinite, Number n)
    : _is_infinite(is_infinite), _n(n) {
  if (is_infinite) {
    if (n > 0)
      _n = 1;
    else
      _n = -1;
  }
}

template <typename Number>
bound<Number> bound<Number>::min(const bound<Number> &x,
                                 const bound<Number> &y) {
  return ((x <= y) ? x : y);
}

template <typename Number>
bound<Number> bound<Number>::min(const bound<Number> &x, const bound<Number> &y,
                                 const bound<Number> &z) {
  return min(x, min(y, z));
}

template <typename Number>
bound<Number> bound<Number>::min(const bound<Number> &x, const bound<Number> &y,
                                 const bound<Number> &z,
                                 const bound<Number> &t) {
  return min(x, min(y, z, t));
}

template <typename Number>
bound<Number> bound<Number>::max(const bound<Number> &x,
                                 const bound<Number> &y) {
  return ((x <= y) ? y : x);
}

template <typename Number>
bound<Number> bound<Number>::max(const bound<Number> &x, const bound<Number> &y,
                                 const bound<Number> &z) {
  return max(x, max(y, z));
}

template <typename Number>
bound<Number> bound<Number>::max(const bound<Number> &x, const bound<Number> &y,
                                 const bound<Number> &z,
                                 const bound<Number> &t) {
  return max(x, max(y, z, t));
}

template <typename Number> bound<Number> bound<Number>::plus_infinity() {
  return bound<Number>(true, 1);
}

template <typename Number> bound<Number> bound<Number>::minus_infinity() {
  return bound<Number>(true, -1);
}

template <typename Number>
bound<Number>::bound(int n) : _is_infinite(false), _n(n) {}

template <typename Number> bound<Number>::bound(std::string s) : _n(1) {
  if (s == "+oo") {
    _is_infinite = true;
  } else if (s == "-oo") {
    _is_infinite = true;
    _n = -1;
  } else {
    _is_infinite = false;
    _n = Number(s);
  }
}

template <typename Number>
bound<Number>::bound(Number n) : _is_infinite(false), _n(n) {}

template <typename Number> bool bound<Number>::is_infinite() const {
  return _is_infinite;
}

template <typename Number> bool bound<Number>::is_finite() const {
  return !_is_infinite;
}

template <typename Number> bool bound<Number>::is_plus_infinity() const {
  return (is_infinite() && _n > 0);
}

template <typename Number> bool bound<Number>::is_minus_infinity() const {
  return (is_infinite() && _n < 0);
}

template <typename Number> bound<Number> bound<Number>::operator-() const {
  return bound<Number>(_is_infinite, -_n);
}

template <typename Number>
bound<Number> bound<Number>::operator+(const bound<Number> &x) const {
  if (is_finite() && x.is_finite()) {
    return bound<Number>(_n + x._n);
  } else if (is_finite() && x.is_infinite()) {
    return x;
  } else if (is_infinite() && x.is_finite()) {
    return *this;
  } else if (_n == x._n) {
    return *this;
  } else {
    CRAB_ERROR("Bound: undefined operation -oo + +oo");
  }
}

template <typename Number>
bound<Number> &bound<Number>::operator+=(const bound<Number> &x) {
  return operator=(operator+(x));
}

template <typename Number>
bound<Number> bound<Number>::operator-(const bound<Number> &x) const {
  return operator+(x.operator-());
}

template <typename Number>
bound<Number> &bound<Number>::operator-=(const bound<Number> &x) {
  return operator=(operator-(x));
}

template <typename Number>
bound<Number> bound<Number>::operator*(const bound<Number> &x) const {
  if (x._n == 0)
    return x;
  else if (_n == 0)
    return *this;
  else
    return bound<Number>(_is_infinite || x._is_infinite, _n * x._n);
}

template <typename Number>
bound<Number> &bound<Number>::operator*=(const bound<Number> &x) {
  return operator=(operator*(x));
}

template <typename Number>
bound<Number> bound<Number>::operator/(const bound<Number> &x) const {
  if (x._n == 0) {
    CRAB_ERROR("Bound: division by zero");
  } else if (is_finite() && x.is_finite()) {
    return bound<Number>(false, _n / x._n);
  } else if (is_finite() && x.is_infinite()) {
    if (_n > 0) {
      return x;
    } else if (_n == 0) {
      return *this;
    } else {
      return x.operator-();
    }
  } else if (is_infinite() && x.is_finite()) {
    if (x._n > 0) {
      return *this;
    } else {
      return operator-();
    }
  } else {
    return bound<Number>(true, _n * x._n);
  }
}

template <typename Number>
bound<Number> &bound<Number>::operator/=(const bound<Number> &x) {
  return operator=(operator/(x));
}

template <typename Number>
bool bound<Number>::operator<(const bound<Number> &x) const {
  return !operator>=(x);
}

template <typename Number>
bool bound<Number>::operator>(const bound<Number> &x) const {
  return !operator<=(x);
}

template <typename Number>
bool bound<Number>::operator==(const bound<Number> &x) const {
  return (_is_infinite == x._is_infinite && _n == x._n);
}

template <typename Number>
bool bound<Number>::operator!=(const bound<Number> &x) const {
  return !operator==(x);
}

template <typename Number>
bool bound<Number>::operator<=(const bound<Number> &x) const {
  if (_is_infinite xor x._is_infinite) {
    if (_is_infinite) {
      return _n < 0;
    }
    return x._n > 0;
  }
  return _n <= x._n;
}

template <typename Number>
bool bound<Number>::operator>=(const bound<Number> &x) const {
  if (_is_infinite xor x._is_infinite) {
    if (_is_infinite) {
      return _n > 0;
    }
    return x._n < 0;
  }
  return _n >= x._n;
}

template <typename Number> bound<Number> bound<Number>::abs() const {
  if (operator>=(0)) {
    return *this;
  } else {
    return operator-();
  }
}

template <typename Number>
boost::optional<Number> bound<Number>::number() const {
  if (is_infinite()) {
    return boost::optional<Number>();
  } else {
    return boost::optional<Number>(_n);
  }
}

template <typename Number> void bound<Number>::write(crab::crab_os &o) const {
  if (is_plus_infinity()) {
    o << "+oo";
  } else if (is_minus_infinity()) {
    o << "-oo";
  } else {
    o << _n;
  }
}

template <typename Number> interval<Number> interval<Number>::top() {
  return interval<Number>(bound<Number>::minus_infinity(),
                          bound<Number>::plus_infinity());
}

template <typename Number> interval<Number> interval<Number>::bottom() {
  return interval<Number>();
}

template <typename Number> interval<Number>::interval() : _lb(0), _ub(-1) {}

template <typename Number> Number interval<Number>::abs(Number x) {
  return x < 0 ? -x : x;
}

template <typename Number> Number interval<Number>::max(Number x, Number y) {
  return (x <= y) ? y : x;
}

template <typename Number> Number interval<Number>::min(Number x, Number y) {
  return (x < y) ? x : y;
}

template <typename Number>
interval<Number>::interval(bound<Number> lb, bound<Number> ub)
    : _lb(lb), _ub(ub) {
  if (lb > ub) {
    _lb = 0;
    _ub = -1;
  }
}

template <typename Number>
interval<Number>::interval(bound<Number> b) : _lb(b), _ub(b) {
  if (b.is_infinite()) {
    _lb = 0;
    _ub = -1;
  }
}

template <typename Number>
interval<Number>::interval(Number n) : _lb(n), _ub(n) {}

template <typename Number>
interval<Number>::interval(std::string b) : _lb(b), _ub(b) {
  if (_lb.is_infinite()) {
    _lb = 0;
    _ub = -1;
  }
}

template <typename Number> bound<Number> interval<Number>::lb() const {
  return _lb;
}

template <typename Number> bound<Number> interval<Number>::ub() const {
  return _ub;
}

template <typename Number> bool interval<Number>::is_bottom() const {
  return (_lb > _ub);
}

template <typename Number> bool interval<Number>::is_top() const {
  return (_lb.is_infinite() && _ub.is_infinite());
}

template <typename Number>
interval<Number> interval<Number>::lower_half_line() const {
  return interval<Number>(bound<Number>::minus_infinity(), _ub);
}

template <typename Number>
interval<Number> interval<Number>::upper_half_line() const {
  return interval<Number>(_lb, bound<Number>::plus_infinity());
}

template <typename Number>
bool interval<Number>::operator==(const interval<Number> &x) const {
  if (is_bottom()) {
    return x.is_bottom();
  } else {
    return (_lb == x._lb) && (_ub == x._ub);
  }
}

template <typename Number>
bool interval<Number>::operator<=(const interval<Number> &x) const {
  if (is_bottom()) {
    return true;
  } else if (x.is_bottom()) {
    return false;
  } else {
    return (x._lb <= _lb) && (_ub <= x._ub);
  }
}

template <typename Number>
interval<Number> interval<Number>::operator|(const interval<Number> &x) const {
  if (is_bottom()) {
    return x;
  } else if (x.is_bottom()) {
    return *this;
  } else {
    return interval<Number>(bound<Number>::min(_lb, x._lb),
                            bound<Number>::max(_ub, x._ub));
  }
}

template <typename Number>
interval<Number> interval<Number>::operator&(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return interval<Number>(bound<Number>::max(_lb, x._lb),
                            bound<Number>::min(_ub, x._ub));
  }
}

template <typename Number>
interval<Number> interval<Number>::operator||(const interval<Number> &x) const {
  if (is_bottom()) {
    return x;
  } else if (x.is_bottom()) {
    return *this;
  } else {
    return interval<Number>(x._lb < _lb ? bound<Number>::minus_infinity() : _lb,
                            _ub < x._ub ? bound<Number>::plus_infinity() : _ub);
  }
}

template <typename Number>
interval<Number> interval<Number>::operator&&(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return interval<Number>(
        _lb.is_infinite() && x._lb.is_finite() ? x._lb : _lb,
        _ub.is_infinite() && x._ub.is_finite() ? x._ub : _ub);
  }
}

template <typename Number>
interval<Number> interval<Number>::operator+(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return interval<Number>(_lb + x._lb, _ub + x._ub);
  }
}

template <typename Number>
interval<Number> &interval<Number>::operator+=(const interval<Number> &x) {
  return operator=(operator+(x));
}

template <typename Number>
interval<Number> interval<Number>::operator-() const {
  if (is_bottom()) {
    return bottom();
  } else {
    return interval<Number>(-_ub, -_lb);
  }
}

template <typename Number>
interval<Number> interval<Number>::operator-(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return interval<Number>(_lb - x._ub, _ub - x._lb);
  }
}

template <typename Number>
interval<Number> &interval<Number>::operator-=(const interval<Number> &x) {
  return operator=(operator-(x));
}

template <typename Number>
interval<Number> interval<Number>::operator*(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    bound<Number> ll = _lb * x._lb;
    bound<Number> lu = _lb * x._ub;
    bound<Number> ul = _ub * x._lb;
    bound<Number> uu = _ub * x._ub;
    return interval<Number>(bound<Number>::min(ll, lu, ul, uu),
                            bound<Number>::max(ll, lu, ul, uu));
  }
}

template <typename Number>
interval<Number> &interval<Number>::operator*=(const interval<Number> &x) {
  return operator=(operator*(x));
}

template <typename Number>
interval<Number> &interval<Number>::operator/=(const interval<Number> &x) {
  return operator=(operator/(x));
}

template <typename Number>
boost::optional<Number> interval<Number>::singleton() const {
  if (!is_bottom() && _lb == _ub) {
    return _lb.number();
  } else {
    return boost::optional<Number>();
  }
}

template <typename Number> bool interval<Number>::operator[](Number n) const {
  if (is_bottom()) {
    return false;
  } else {
    bound<Number> b(n);
    return (_lb <= b) && (b <= _ub);
  }
}

template <typename Number>
void interval<Number>::write(crab::crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else {
    o << "[" << _lb << ", " << _ub << "]";
  }
}

template <typename Number>
interval<Number> interval<Number>::UDiv(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::SRem(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::URem(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::And(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::Or(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::Xor(const interval<Number> &x) const {
  return Or(x);
}

template <typename Number>
interval<Number> interval<Number>::Shl(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::LShr(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}

template <typename Number>
interval<Number> interval<Number>::AShr(const interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    return top();
  }
}
} // namespace ikos
