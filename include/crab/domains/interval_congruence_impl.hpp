#include <boost/optional.hpp>
#include <crab/domains/interval_congruence.hpp>

namespace crab {
namespace domains {

template <typename Number>
interval_congruence<Number>::interval_congruence(bool is_bottom)
    : m_first(is_bottom ? ikos::interval<Number>::bottom()
                        : ikos::interval<Number>::top()),
      m_second(is_bottom ? ikos::congruence<Number>::bottom()
                         : ikos::congruence<Number>::top()) {}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::top() {
  return interval_congruence(false);
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::bottom() {
  return interval_congruence(true);
}

template <typename Number>
Number interval_congruence<Number>::abs(Number x) const {
  return x < 0 ? -x : x;
}

// operator % can return a negative number
// mod(a, b) always returns a positive number
template <typename Number>
Number interval_congruence<Number>::mod(Number a, Number b) const {
  Number m = a % b;
  if (m < 0)
    return m + b;
  else
    return m;
}

// R(c,a) is the least element of c greater or equal than a
template <typename Number>
Number interval_congruence<Number>::R(ikos::congruence<Number> c,
                                      Number a) const {
  Number m = c.get_modulo();
  Number p = c.get_remainder();
  return a + mod(p - a, abs(m));
}

// L(c,a) is the greatest element of c smaller or equal than a
template <typename Number>
Number interval_congruence<Number>::L(ikos::congruence<Number> c,
                                      Number a) const {
  Number m = c.get_modulo();
  Number p = c.get_remainder();
  return a - mod(a - p, abs(m));
}

template <typename Number>
interval_congruence<Number>::interval_congruence(Number n)
    : m_first(ikos::interval<Number>(n)),
      m_second(ikos::congruence<Number>(n)) {}

template <typename Number>
interval_congruence<Number>::interval_congruence(ikos::interval<Number> &&i,
                                                 ikos::congruence<Number> &&c)
    : m_first(std::move(i)), m_second(std::move(c)) {
  reduce();
}

template <typename Number>
interval_congruence<Number>::interval_congruence(ikos::interval<Number> &&i)
    : m_first(std::move(i)), m_second(ikos::congruence<Number>::top()) {
  reduce();
}

template <typename Number>
interval_congruence<Number>::interval_congruence(ikos::congruence<Number> &&c)
    : m_first(ikos::interval<Number>::top()), m_second(std::move(c)) {
  reduce();
}

template <typename Number> bool interval_congruence<Number>::is_bottom() {
  return m_first.is_bottom() || m_second.is_bottom();
}

template <typename Number> bool interval_congruence<Number>::is_top() {
  return m_first.is_top() && m_second.is_top();
}

template <typename Number>
ikos::interval<Number> &interval_congruence<Number>::first() {
  return m_first;
}
template <typename Number>
const ikos::interval<Number> &interval_congruence<Number>::first() const {
  return m_first;
}

template <typename Number>
ikos::congruence<Number> &interval_congruence<Number>::second() {
  return m_second;
}
template <typename Number>
const ikos::congruence<Number> &interval_congruence<Number>::second() const {
  return m_second;
}

/*
   Let (i,c) be a pair of interval and congruence these are the
   main rules described by Granger:

   if (c.is_bottom() || i.is_bottom()) (bottom(), bottom());
   if (c = 0Z+a and a \notin i)        (bottom(), bottom());
   if (c = 0Z+a)                       ([a,a]   , c);
   if (i=[a,b] and R(c,a) > L(c,b))    (bottom(), bottom());
   if (i=[a,b])                        ([R(c,a), L(c,b)], c);
   if (i=[a,+oo])                      ([R(c,a), +oo], c);
   if (i=[-oo,b])                      ([-oo, L(c,b)], c);
   otherwise                           (i,c)
 */
template <typename Number> void interval_congruence<Number>::reduce() {
  ikos::interval<Number> &i = first();
  ikos::congruence<Number> &c = second();

  if (i.is_bottom() || c.is_bottom()) {
    i = ikos::interval<Number>::bottom();
    c = ikos::congruence<Number>::bottom();
  }

  // congruence is top and interval is a singleton
  if (c.is_top()) {
    boost::optional<Number> n = i.singleton();
    if (n) {
      c = ikos::congruence<Number>(*n);
    }
    return;
  }

  Number modulo = c.get_modulo();
  if (modulo == 0) {
    // congruence is a singleton so we refine the interval
    ikos::interval<Number> a(c.get_remainder());
    if (!(a <= i)) {
      i = ikos::interval<Number>::bottom();
      c = ikos::congruence<Number>::bottom();
    } else {
      i = a;
    }
  } else {
    // refine lower and upper bounds of the interval using
    // congruences
    bound_t lb = i.lb();
    bound_t ub = i.ub();

    if (lb.is_finite() && ub.is_finite()) {
      Number x = R(c, *(lb.number()));
      Number y = L(c, *(ub.number()));
      if (x > y) {
        i = ikos::interval<Number>::bottom();
        c = ikos::congruence<Number>::bottom();
      } else if (x == y) {
        i = ikos::interval<Number>(x);
        c = ikos::congruence<Number>(x);
      } else {
        i = ikos::interval<Number>(bound_t(x), bound_t(y));
      }
    } else if (lb.is_finite()) {
      Number x = R(c, *(lb.number()));
      i = ikos::interval<Number>(bound_t(x), bound_t::plus_infinity());
    } else if (ub.is_finite()) {
      Number y = L(c, *(ub.number()));
      i = ikos::interval<Number>(bound_t::minus_infinity(), bound_t(y));
    } else {
      // interval is top
    }
  }
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator+(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.operator+(x.first())),
                                     std::move(m_second.operator+(x.second())));
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator-(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.operator-(x.first())),
                                     std::move(m_second.operator-(x.second())));
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator*(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.operator*(x.first())),
                                     std::move(m_second.operator*(x.second())));
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator/(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.operator/(x.first())),
                                     std::move(m_second.operator/(x.second())));
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator|(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first | x.m_first),
                                     std::move(m_second | x.m_second));
}

template <typename Number>
interval_congruence<Number> interval_congruence<Number>::operator&(
    const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first & x.m_first),
                                     std::move(m_second & x.m_second));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::SDiv(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first / x.first()),
                                     std::move(m_second.SDiv(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::UDiv(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.UDiv(x.first())),
                                     std::move(m_second.UDiv(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::SRem(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.SRem(x.first())),
                                     std::move(m_second.SRem(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::URem(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.URem(x.first())),
                                     std::move(m_second.URem(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::Trunc(unsigned width) const {
  interval_congruence<Number> res = interval_congruence<Number>::top();
  return res;
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::ZExt(unsigned width) const {
  interval_congruence<Number> res = interval_congruence<Number>::top();
  return res;
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::SExt(unsigned width) const {
  interval_congruence<Number> res = interval_congruence<Number>::top();
  return res;
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::And(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.And(x.first())),
                                     std::move(m_second.And(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::Or(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.Or(x.first())),
                                     std::move(m_second.Or(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::Xor(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.Xor(x.first())),
                                     std::move(m_second.Xor(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::Shl(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.Shl(x.first())),
                                     std::move(m_second.Shl(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::LShr(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.LShr(x.first())),
                                     std::move(m_second.LShr(x.second())));
}

template <typename Number>
interval_congruence<Number>
interval_congruence<Number>::AShr(const interval_congruence<Number> &x) const {
  return interval_congruence<Number>(std::move(m_first.AShr(x.first())),
                                     std::move(m_second.AShr(x.second())));
}

template <typename Number>
void interval_congruence<Number>::write(crab_os &o) const {
  o << "(";
  m_first.write(o);
  o << ", ";
  m_second.write(o);
  o << ")";
}

} // end namespace domains
} // end namespace crab
