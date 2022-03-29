#pragma once

/*******************************************************************************
 * A simple class for representing intervals and performing interval
 * arithmetic.
 ******************************************************************************/

#include <crab/domains/linear_interval_solver.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <string>

namespace ikos {

template <typename Number> class bound {
public:
  using bound_t = bound<Number>;

private:
  bool _is_infinite;
  Number _n;

  bound();

  bound(bool is_infinite, Number n);

public:
  static bound_t min(const bound_t &x, const bound_t &y);

  static bound_t min(const bound_t &x, const bound_t &y, const bound_t &z);

  static bound_t min(const bound_t &x, const bound_t &y, const bound_t &z,
                     const bound_t &t);

  static bound_t max(const bound_t &x, const bound_t &y);

  static bound_t max(const bound_t &x, const bound_t &y, const bound_t &z);

  static bound_t max(const bound_t &x, const bound_t &y, const bound_t &z,
                     const bound_t &t);

  static bound_t plus_infinity();

  static bound_t minus_infinity();

  bound(int n);

  bound(std::string s);

  bound(Number n);

  bound(const bound_t &o) = default;

  bound_t &operator=(const bound_t &o) = default;

  bound(bound_t &&o) = default;

  bound_t &operator=(bound_t &&o) = default;

  bool is_infinite() const;

  bool is_finite() const;

  bool is_plus_infinity() const;

  bool is_minus_infinity() const;

  bound_t operator-() const;

  bound_t operator+(const bound_t &x) const;

  bound_t &operator+=(const bound_t &x);

  bound_t operator-(const bound_t &x) const;

  bound_t &operator-=(const bound_t &x);

  bound_t operator*(const bound_t &x) const;

  bound_t &operator*=(const bound_t &x);

  bound_t operator/(const bound_t &x) const;

  bound_t &operator/=(const bound_t &x);

  bool operator<(const bound_t &x) const;

  bool operator>(const bound_t &x) const;

  bool operator==(const bound_t &x) const;

  bool operator!=(const bound_t &x) const;

  bool operator<=(const bound_t &x) const;

  bool operator>=(const bound_t &x) const;

  bound_t abs() const;

  boost::optional<Number> number() const;

  void write(crab::crab_os &o) const;

  friend crab::crab_os &operator<<(crab::crab_os &o, const bound<Number> &b) {
    b.write(o);
    return o;
  }
}; // class bound

namespace bounds_impl {
// Conversion between z_bound and q_bound
template <class B1, class B2> void convert_bounds(B1 b1, B2 &b2);
} // namespace bounds_impl

template <typename Number> class interval {
public:
  using bound_t = bound<Number>;
  using interval_t = interval<Number>;

private:
  bound_t _lb;
  bound_t _ub;

  interval();

  static Number abs(Number x);

  static Number max(Number x, Number y);

  static Number min(Number x, Number y);

public:
  static interval_t top();

  static interval_t bottom();

  interval(bound_t lb, bound_t ub);

  interval(bound_t b);

  interval(Number n);

  interval(std::string b);

  interval(const interval_t &i) = default;

  interval_t &operator=(const interval_t &i) = default;

  interval(interval_t &&i) = default;

  interval_t &operator=(interval_t &&i) = default;

  bound_t lb() const;

  bound_t ub() const;

  bool is_bottom() const;

  bool is_top() const;

  interval_t lower_half_line() const;

  interval_t upper_half_line() const;

  bool operator==(const interval_t &x) const;

  bool operator!=(const interval_t &x) const { return !operator==(x); }

  bool operator<=(const interval_t &x) const;

  interval_t operator|(const interval_t &x) const;

  interval_t operator&(const interval_t &x) const;

  interval_t operator||(const interval_t &x) const;

  template <typename Thresholds>
  interval_t widening_thresholds(const interval_t &x,
                                 const Thresholds &ts) const {
    if (is_bottom()) {
      return x;
    } else if (x.is_bottom()) {
      return *this;
    } else {
      bound<Number> lb = (x._lb < _lb ? ts.get_prev(x._lb) : _lb);
      bound<Number> ub = (_ub < x._ub ? ts.get_next(x._ub) : _ub);
      return interval<Number>(lb, ub);
    }
  }

  interval_t operator&&(const interval_t &x) const;

  interval_t operator+(const interval_t &x) const;

  interval_t &operator+=(const interval_t &x);

  interval_t operator-() const;

  interval_t operator-(const interval_t &x) const;

  interval_t &operator-=(const interval_t &x);

  interval_t operator*(const interval_t &x) const;

  interval_t &operator*=(const interval_t &x);

  interval_t operator/(const interval_t &x) const;

  interval_t &operator/=(const interval_t &x);

  boost::optional<Number> singleton() const;

  bool operator[](Number n) const;

  // division and remainder operations
  interval_t UDiv(const interval_t &x) const;
  interval_t SRem(const interval_t &x) const;
  interval_t URem(const interval_t &x) const;

  // bitwise operations
  interval_t And(const interval_t &x) const;
  interval_t Or(const interval_t &x) const;
  interval_t Xor(const interval_t &x) const;
  interval_t Shl(const interval_t &x) const;
  interval_t LShr(const interval_t &x) const;
  interval_t AShr(const interval_t &x) const;

  void write(crab::crab_os &o) const;
  friend crab::crab_os &operator<<(crab::crab_os &o, const interval_t &i) {
    i.write(o);
    return o;
  }
}; //  class interval

template <typename Number>
inline interval<Number> operator+(Number c, const interval<Number> &x) {
  return interval<Number>(c) + x;
}

template <typename Number>
inline interval<Number> operator+(const interval<Number> &x, Number c) {
  return x + interval<Number>(c);
}

template <typename Number>
inline interval<Number> operator*(Number c, const interval<Number> &x) {
  return interval<Number>(c) * x;
}

template <typename Number>
inline interval<Number> operator*(const interval<Number> &x, Number c) {
  return x * interval<Number>(c);
}

template <typename Number>
inline interval<Number> operator/(Number c, const interval<Number> &x) {
  return interval<Number>(c) / x;
}

template <typename Number>
inline interval<Number> operator/(const interval<Number> &x, Number c) {
  return x / interval<Number>(c);
}

template <typename Number>
inline interval<Number> operator-(Number c, const interval<Number> &x) {
  return interval<Number>(c) - x;
}

template <typename Number>
inline interval<Number> operator-(const interval<Number> &x, Number c) {
  return x - interval<Number>(c);
}

namespace bounds_impl {
void convert_bounds(bound<z_number> b1, bound<z_number> &b2);
void convert_bounds(bound<q_number> b1, bound<q_number> &b2);
void convert_bounds(bound<z_number> b1, bound<q_number> &b2);
void convert_bounds(bound<q_number> b1, bound<z_number> &b2);
} // namespace bounds_impl

namespace linear_interval_solver_impl {
template <>
interval<z_number> trim_interval(const interval<z_number> &i,
                                 const interval<z_number> &j);
template <>
interval<q_number> trim_interval(const interval<q_number> &i,
                                 const interval<q_number> &j);
template <>
interval<z_number> lower_half_line(const interval<z_number> &i, bool is_signed);
template <>
interval<q_number> lower_half_line(const interval<q_number> &i, bool is_signed);
template <>
interval<z_number> upper_half_line(const interval<z_number> &i, bool is_signed);
template <>
interval<q_number> upper_half_line(const interval<q_number> &i, bool is_signed);
} // namespace linear_interval_solver_impl

} // namespace ikos
