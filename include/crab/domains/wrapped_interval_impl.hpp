#pragma once

#include <crab/domains/wrapped_interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

#define PRINT_WRAPINT_AS_SIGNED

namespace crab {
namespace domains {

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::default_implementation(
    const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  } else {
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>::wrapped_interval(wrapint start, wrapint end,
                                           bool is_bottom)
    : m_start(start), m_end(end), m_is_bottom(is_bottom) {
  if (start.get_bitwidth() != end.get_bitwidth()) {
    CRAB_ERROR("inconsistent bitwidths in wrapped interval");
  }
}

template <typename Number>
void wrapped_interval<Number>::signed_split(
    std::vector<wrapped_interval<Number>> &intervals) const {
  if (is_bottom())
    return;

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);
  if (is_top()) {
    intervals.push_back(wrapped_interval<Number>(wrapint::get_unsigned_min(b),
                                                 wrapint::get_signed_max(b)));
    intervals.push_back(wrapped_interval<Number>(wrapint::get_signed_min(b),
                                                 wrapint::get_unsigned_max(b)));
  } else {
    if (signed_limit(b) <= *this) {
      intervals.push_back(
          wrapped_interval<Number>(m_start, wrapint::get_signed_max(b)));
      intervals.push_back(
          wrapped_interval<Number>(wrapint::get_signed_min(b), m_end));
    } else {
      intervals.push_back(*this);
    }
  }
}

template <typename Number>
void wrapped_interval<Number>::unsigned_split(
    std::vector<wrapped_interval<Number>> &intervals) const {
  if (is_bottom())
    return;

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);
  if (is_top()) {
    intervals.push_back(wrapped_interval<Number>(wrapint::get_signed_min(b),
                                                 wrapint::get_unsigned_max(b)));
    intervals.push_back(wrapped_interval<Number>(wrapint::get_unsigned_min(b),
                                                 wrapint::get_signed_max(b)));
  } else {
    if (unsigned_limit(b) <= *this) {
      intervals.push_back(
          wrapped_interval<Number>(m_start, wrapint::get_unsigned_max(b)));
      intervals.push_back(
          wrapped_interval<Number>(wrapint::get_unsigned_min(b), m_end));
    } else {
      intervals.push_back(*this);
    }
  }
}

template <typename Number>
void wrapped_interval<Number>::signed_and_unsigned_split(
    std::vector<wrapped_interval<Number>> &out) const {
  std::vector<wrapped_interval<Number>> ssplit;
  signed_split(ssplit);
  for (unsigned i = 0, e = ssplit.size(); i < e; ++i) {
    ssplit[i].unsigned_split(out);
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::signed_mul(const wrapped_interval<Number> &x) const {
  assert(!is_bottom() && !x.is_bottom());

  bool msb_start = m_start.msb();
  bool msb_end = m_end.msb();
  bool msb_x_start = x.m_start.msb();
  bool msb_x_end = x.m_end.msb();
  wrapint::bitwidth_t b = get_bitwidth(__LINE__);

  wrapped_interval<Number> res = wrapped_interval<Number>::top();
  if (msb_start == msb_end && msb_end == msb_x_start &&
      msb_x_start == msb_x_end) {
    // the two intervals are in the same hemisphere
    if (!msb_start) {
      return unsigned_mul(x);
    } else {
      // -- check if multiplication will overflow
      if ((m_start.get_unsigned_bignum() * x.m_start.get_unsigned_bignum()) -
              (m_end.get_unsigned_bignum() * x.m_end.get_unsigned_bignum()) <
          wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
        res = wrapped_interval<Number>(m_end * x.m_end, m_start * x.m_start);
      }
      goto EXIT_SIGNED_MUL;
    }
  }

  // each interval cannot cross the limit: one interval in a
  // different hemisphere.
  if (!(msb_start != msb_end || msb_x_start != msb_x_end)) {
    if (msb_start && !msb_x_start) {
      // -- check if multiplication will overflow
      if ((m_end.get_unsigned_bignum() * x.m_start.get_unsigned_bignum()) -
              (m_start.get_unsigned_bignum() * x.m_end.get_unsigned_bignum()) <
          wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
        res = wrapped_interval<Number>(m_start * x.m_end, m_end * x.m_start);
      }
    } else if (!msb_start && msb_x_start) {
      // -- check if multiplication will overflow
      if ((m_start.get_unsigned_bignum() * x.m_end.get_unsigned_bignum()) -
              (m_end.get_unsigned_bignum() * x.m_start.get_unsigned_bignum()) <
          wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
        res = wrapped_interval<Number>(m_end * x.m_start, m_start * x.m_end);
      }
    }
  }

EXIT_SIGNED_MUL:
  CRAB_LOG("wrapped-int-mul", crab::outs() << "Signed " << *this << " * " << x
                                           << "=" << res << "\n";);
  return res;
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::unsigned_mul(
    const wrapped_interval<Number> &x) const {
  assert(!is_bottom() && !x.is_bottom());

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);
  wrapped_interval<Number> res = wrapped_interval<Number>::top();
  // -- check if multiplication will overflow
  if ((m_end.get_unsigned_bignum() * x.m_end.get_unsigned_bignum()) -
          (m_start.get_unsigned_bignum() * x.m_start.get_unsigned_bignum()) <
      wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
    res = wrapped_interval<Number>(m_start * x.m_start, m_end * x.m_end);
  }
  CRAB_LOG("wrapped-int-mul", crab::outs() << "Unsigned " << *this << " * " << x
                                           << "=" << res << "\n";);
  return res;
}

// if out is empty then the intersection is empty
template <typename Number>
void wrapped_interval<Number>::exact_meet(
    const wrapped_interval<Number> &x,
    std::vector<wrapped_interval<Number>> &out) const {
  if (is_bottom() || x.is_bottom()) {
    // bottom
  } else if (*this == x || is_top()) {
    out.push_back(x);
  } else if (x.is_top()) {
    out.push_back(*this);
  } else if (x.at(m_start) && x.at(m_end) && at(x.m_start) && at(x.m_end)) {
    out.push_back(wrapped_interval<Number>(m_start, x.m_end));
    out.push_back(wrapped_interval<Number>(x.m_start, m_end));
  } else if (x.at(m_start) && x.at(m_end)) {
    out.push_back(*this);
  } else if (at(x.m_start) && at(x.m_end)) {
    out.push_back(x);
  } else if (x.at(m_start) && at(x.m_end) && !x.at(m_end) && at(x.m_start)) {
    out.push_back(wrapped_interval<Number>(m_start, x.m_end));
  } else if (x.at(m_end) && at(x.m_start) && !x.at(m_start) && at(x.m_end)) {
    out.push_back(wrapped_interval<Number>(x.m_start, m_end));
  } else {
    // bottom
  }
}

// Perform the reduced product of signed and unsigned multiplication.
// It uses exact meet rather than abstract meet.
template <typename Number>
void wrapped_interval<Number>::reduced_signed_unsigned_mul(
    const wrapped_interval<Number> &x,
    std::vector<wrapped_interval<Number>> &out) const {
  if (is_bottom() || x.is_bottom()) {
    return;
  }

  wrapped_interval<Number> s = signed_mul(x);
  wrapped_interval<Number> u = unsigned_mul(x);
  s.exact_meet(u, out);
  CRAB_LOG(
      "wrapped-int-mul",
      crab::outs() << "Exact signed x unsigned " << s << " * " << u << "=\n";
      for (unsigned i = 0; i < out.size();
           ++i) { crab::outs() << "\t" << out[i] << "\n"; });
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::unsigned_div(
    const wrapped_interval<Number> &x) const {
  CRAB_LOG("wrapped-int-div", crab::outs() << *this << " /_u " << x << "=";);
  assert(!x.at(wrapint(0, x.get_bitwidth(__LINE__))));
  assert(!is_bottom() && !x.is_bottom());
  wrapped_interval<Number> res =
      wrapped_interval<Number>(m_start.udiv(x.m_end), m_end.udiv(x.m_start));
  CRAB_LOG("wrapped-int-div", crab::outs() << res << "\n";);
  return res;
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::signed_div(const wrapped_interval<Number> &x) const {
  CRAB_LOG("wrapped-int-div", crab::outs() << *this << " /_s " << x << "=";);
  assert(!x.at(wrapint(0, x.get_bitwidth(__LINE__))));
  assert(!is_bottom() && !x.is_bottom());

  bool msb_start = m_start.msb();
  bool msb_x_start = x.m_start.msb();

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);
  wrapint smin = wrapint::get_signed_min(b);
  wrapint minus_one(-1, b);

  wrapped_interval<Number> res = wrapped_interval<Number>::top();
  if (msb_start == msb_x_start) {
    if (msb_start) { // both negative
      // -- check if division will overflow
      if (!((m_end == smin && x.m_start == minus_one) ||
            (m_start == smin && x.m_end == minus_one))) {
        res = wrapped_interval<Number>(m_end.sdiv(x.m_start),
                                       m_start.sdiv(x.m_end));
      }
    } else { // both positive
      // -- check if division will overflow
      if (!((m_start == smin && x.m_end == minus_one) ||
            (m_end == smin && x.m_start == minus_one))) {
        res = wrapped_interval<Number>(m_start.sdiv(x.m_end),
                                       m_end.sdiv(x.m_start));
      }
    }
  } else {
    if (msb_start) {
      assert(!msb_x_start);
      // -- check if division will overflow
      if (!((m_start == smin && x.m_start == minus_one) ||
            (m_end == smin && x.m_end == minus_one))) {
        res = wrapped_interval<Number>(m_start.sdiv(x.m_start),
                                       m_end.sdiv(x.m_end));
      }
    } else {
      assert(msb_x_start);
      // -- check if division will overflow
      if (!((m_end == smin && x.m_end == minus_one) ||
            (m_start == smin && x.m_start == minus_one))) {
        res = wrapped_interval<Number>(m_end.sdiv(x.m_end),
                                       m_start.sdiv(x.m_start));
      }
    }
  }
  CRAB_LOG("wrapped-int-div", crab::outs() << res << "\n";);
  return res;
}

/// This is sound only if wrapped interval defined over z_number.
template <typename Number>
void wrapped_interval<Number>::trim_zero(
    std::vector<wrapped_interval<Number>> &out) const {
  wrapint zero(0, get_bitwidth(__LINE__));
  if (!is_bottom() && (!(*this == zero))) {
    if (start() == zero) {
      out.push_back(
          wrapped_interval<Number>(wrapint(1, get_bitwidth(__LINE__)), end()));
    } else if (end() == zero) {
      out.push_back(wrapped_interval<Number>(
          start(), wrapint(-1, get_bitwidth(__LINE__))));
    } else if (at(zero)) {
      out.push_back(wrapped_interval<Number>(
          start(), wrapint(-1, get_bitwidth(__LINE__))));
      out.push_back(
          wrapped_interval<Number>(wrapint(1, get_bitwidth(__LINE__)), end()));
    } else {
      out.push_back(*this);
    }
  }
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::Shl(uint64_t k) const {
  if (is_bottom())
    return *this;

  // XXX: we need the check is_top before calling get_bitwidth();
  if (is_top())
    return *this;

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);
  wrapped_interval<Number> y = Trunc(b - k);
  if (!y.is_top()) {
    wrapint wk(k, b);
    return wrapped_interval<Number>(start() << wk, end() << wk);
  } else {
    // TODO: return (0, 1..10..0) where #1's = b-k and #0's= k
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::LShr(uint64_t k) const {
  if (is_bottom())
    return *this;

  // XXX: we need the check is_top before cross_signed_limit calls
  // get_bitwidth();
  if (is_top())
    return *this;

  if (!cross_unsigned_limit()) {
    wrapint::bitwidth_t b = get_bitwidth(__LINE__);
    wrapint wk(k, b);
    return wrapped_interval<Number>(start().lshr(wk), end().lshr(wk));
  } else {
    // TODO: return (0, 0..01..1 where #0's = k and #1's = b-k
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::AShr(uint64_t k) const {
  if (is_bottom())
    return *this;

  // XXX: we need the check is_top before cross_signed_limit calls
  // get_bitwidth();
  if (is_top())
    return *this;

  if (!cross_signed_limit()) {
    wrapint::bitwidth_t b = get_bitwidth(__LINE__);
    wrapint wk(k, b);
    return wrapped_interval<Number>(start().ashr(wk), end().ashr(wk));
  } else {
    // TODO: return (1...10..0,        0..01..1)
    //               #1's=k #0's=b-k   #0's=k #1's=b-k
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>::wrapped_interval(wrapint n)
    : m_start(n), m_end(n), m_is_bottom(false) {}

template <typename Number>
wrapped_interval<Number>::wrapped_interval(wrapint start, wrapint end)
    : m_start(start), m_end(end), m_is_bottom(false) {
  if (start.get_bitwidth() != end.get_bitwidth()) {
    CRAB_ERROR("inconsistent bitwidths in wrapped interval");
  }
}

// To represent top, the particular bitwidth here is irrelevant. We
// just make sure that m_end - m_start == get_max()
template <typename Number>
wrapped_interval<Number>::wrapped_interval()
    : m_start(wrapint(0, 3)), m_end(7, 3), m_is_bottom(false) {}

// return top if n does not fit into a wrapint. No runtime errors.
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::mk_winterval(Number n, wrapint::bitwidth_t width) {
  if (wrapint::fits_wrapint(n, width)) {
    return wrapped_interval<Number>(wrapint(n, width));
  } else {
    CRAB_WARN(n, " does not fit into a wrapint. Returned top wrapped interval");
    return wrapped_interval<Number>::top();
  }
}

// Return top if lb or ub do not fit into a wrapint. No runtime errors.
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::mk_winterval(Number lb, Number ub,
                                       wrapint::bitwidth_t width) {
  if (!wrapint::fits_wrapint(lb, width)) {
    CRAB_WARN(lb,
              " does not fit into a wrapint. Returned top wrapped interval");
    return wrapped_interval<Number>::top();
  } else if (!wrapint::fits_wrapint(ub, width)) {
    CRAB_WARN(ub,
              " does not fit into a wrapint. Returned top wrapped interval");
    return wrapped_interval<Number>::top();
  } else {
    return wrapped_interval<Number>(wrapint(lb, width), wrapint(ub, width));
  }
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::top() {
  return wrapped_interval<Number>(wrapint(0, 3), wrapint(7, 3), false);
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::bottom() {
  // the wrapint is irrelevant.
  wrapint i(0, 1);
  return wrapped_interval<Number>(i, i, true);
}

// return interval [0111...1, 1000....0]
// In the APLAS'12 paper "signed limit" corresponds to "north pole".
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::signed_limit(wrapint::bitwidth_t b) {
  return wrapped_interval<Number>(wrapint::get_signed_max(b),
                                  wrapint::get_signed_min(b));
}

// return interval [1111...1, 0000....0]
// In the APLAS'12 paper "unsigned limit" corresponds to "south pole".
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::unsigned_limit(wrapint::bitwidth_t b) {
  return wrapped_interval<Number>(wrapint::get_unsigned_max(b),
                                  wrapint::get_unsigned_min(b));
}

template <typename Number>
bool wrapped_interval<Number>::cross_signed_limit() const {
  return (signed_limit(get_bitwidth(__LINE__)) <= *this);
}

template <typename Number>
bool wrapped_interval<Number>::cross_unsigned_limit() const {
  return (unsigned_limit(get_bitwidth(__LINE__)) <= *this);
}

template <typename Number>
wrapint::bitwidth_t wrapped_interval<Number>::get_bitwidth(int line) const {
  if (is_bottom()) {
    CRAB_ERROR("get_bitwidth() cannot be called from a bottom element at line ",
               line);
  } else if (is_top()) {
    CRAB_ERROR("get_bitwidth() cannot be called from a top element at line ",
               line);
  } else {
    assert(m_start.get_bitwidth() == m_end.get_bitwidth());
    return m_start.get_bitwidth();
  }
}

template <typename Number> wrapint wrapped_interval<Number>::start() const {
  if (is_top()) {
    CRAB_ERROR("method start() cannot be called if top");
  }
  return m_start;
}

template <typename Number> wrapint wrapped_interval<Number>::end() const {
  if (is_top()) {
    CRAB_ERROR("method end() cannot be called if top");
  }
  return m_end;
}

template <typename Number> bool wrapped_interval<Number>::is_bottom() const {
  return m_is_bottom;
}

template <typename Number> bool wrapped_interval<Number>::is_top() const {
  wrapint maxspan = wrapint::get_unsigned_max(m_start.get_bitwidth());
  return (!m_is_bottom && (m_end - m_start == maxspan));
}

template <typename Number>
ikos::interval<Number> wrapped_interval<Number>::to_interval() const {
  using interval_t = ikos::interval<Number>;
  if (is_bottom()) {
    return interval_t::bottom();
  } else if (is_top() || (cross_signed_limit())) {
    return interval_t::top();
  } else {
    return interval_t(m_start.get_signed_bignum(), m_end.get_signed_bignum());
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::lower_half_line(bool is_signed) const {
  if (is_top() || is_bottom())
    return *this;

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);

  if (is_signed) {
    if (at(wrapint::get_signed_max(b))) {
      return wrapped_interval<Number>::top();
    }
  } else {
    if (at(wrapint::get_unsigned_max(b))) {
      return wrapped_interval<Number>::top();
    }
  }

  wrapint smin = wrapint::get_signed_min(b);
  if (!is_signed) {
    smin = wrapint::get_unsigned_min(b);
  }
  wrapped_interval<Number> res = wrapped_interval<Number>(smin, m_end);
  return res;
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::upper_half_line(bool is_signed) const {
  if (is_top() || is_bottom())
    return *this;

  wrapint::bitwidth_t b = get_bitwidth(__LINE__);

  if (is_signed) {
    if (at(wrapint::get_signed_min(b))) {
      return wrapped_interval<Number>::top();
    }
  } else {
    if (at(wrapint::get_unsigned_min(b))) {
      return wrapped_interval<Number>::top();
    }
  }

  wrapint smax = wrapint::get_signed_max(b);
  if (!is_signed) {
    smax = wrapint::get_unsigned_max(b);
  }
  wrapped_interval<Number> res = wrapped_interval<Number>(m_start, smax);
  return res;
}

template <typename Number> bool wrapped_interval<Number>::is_singleton() const {
  return (!is_bottom() && !is_top() && m_start == m_end);
}

// Starting from m_start and going clock-wise x is encountered
// before than m_end.
template <typename Number> bool wrapped_interval<Number>::at(wrapint x) const {
  if (is_bottom()) {
    return false;
  } else if (is_top()) {
    return true;
  } else {
    return ((x - m_start) <= (m_end - m_start));
  }
}

template <typename Number>
bool wrapped_interval<Number>::operator<=(
    const wrapped_interval<Number> &x) const {
  if (x.is_top() || is_bottom()) {
    return true;
  } else if (x.is_bottom() || is_top()) {
    return false;
  } else if (m_start == x.m_start && m_end == x.m_end)
    return true;
  else {
    return x.at(m_start) && x.at(m_end) && (!(at(x.m_start)) || !(at(x.m_end)));
  }
}

template <typename Number>
bool wrapped_interval<Number>::operator==(
    const wrapped_interval<Number> &x) const {
  return (*this <= x && x <= *this);
}

template <typename Number>
bool wrapped_interval<Number>::operator!=(
    const wrapped_interval<Number> &x) const {
  return !(this->operator==(x));
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator|(const wrapped_interval<Number> &x) const {
  if (*this <= x) {
    return x;
  } else if (x <= *this) {
    return *this;
  } else {
    if (x.at(m_start) && x.at(m_end) && at(x.m_start) && at(x.m_end)) {
      return wrapped_interval<Number>::top();
    } else if (x.at(m_end) && at(x.m_start)) {
      return wrapped_interval<Number>(m_start, x.m_end);
    } else if (at(x.m_end) && x.at(m_start)) {
      return wrapped_interval<Number>(x.m_start, m_end);
    } else {
      wrapint span_a = x.m_start - m_end;
      wrapint span_b = m_start - x.m_end;
      if (span_a < span_b || (span_a == span_b && m_start <= x.m_start)) {
        return wrapped_interval<Number>(m_start, x.m_end);
      } else {
        return wrapped_interval<Number>(x.m_start, m_end);
      }
    }
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator&(const wrapped_interval<Number> &x) const {
  if (*this <= x) {
    return *this;
  } else if (x <= *this) {
    return x;
  } else {
    if (x.at(m_start)) {
      if (at(x.m_start)) {
        wrapint span_a = m_end - m_start;
        wrapint span_b = x.m_end - x.m_start;
        if (span_a < span_b || (span_a == span_b && m_start <= x.m_start)) {
          return *this;
        } else {
          return x; // imprecision here {(m_end, x.m_end), (x.m_start, start)}
        }
      } else {
        if (x.at(m_end)) {
          return *this;
        } else {
          return wrapped_interval<Number>(m_start, x.m_end);
        }
      }
    } else {
      if (at(x.m_start)) {
        if (at(x.m_end)) {
          return x; // imprecision here {(x.m_start, m_end), (m_start, x.m_end)}
        } else {
          return wrapped_interval<Number>(x.m_start, m_end);
        }
      } else {
        return wrapped_interval<Number>::bottom();
      }
    }
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator||(const wrapped_interval<Number> &x) const {
  if (is_bottom()) {
    return x;
  } else if (x.is_bottom()) {
    return *this;
  } else if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else if (x <= *this) {
    return *this;
  }
  assert(get_bitwidth(__LINE__) == x.get_bitwidth(__LINE__));
  wrapint::bitwidth_t w = x.get_bitwidth(__LINE__);

  /**
     TODO: make it template parameter.
     growth_rate should be power of 2. Currently only 2,4,8, or 16.
  **/
  const unsigned growth_rate = 8;
  wrapint max(0, w);
  switch (growth_rate) {
  case 4:
    if (w > 2) {
      max = wrapint(1 << (w - 2), w);
      break;
    }
    BOOST_FALLTHROUGH;     
  case 8:
    if (w > 3) {
      max = wrapint(1 << (w - 3), w);
      break;
    }
    BOOST_FALLTHROUGH;     
  case 16:
    if (w > 4) {
      max = wrapint(1 << (w - 4), w);
      break;
    }
    BOOST_FALLTHROUGH;     
  case 2:
    BOOST_FALLTHROUGH;         
  default:
    assert(w > 1);
    max = wrapint(1 << (w - 1), w);
  }
  if ((m_end - m_start) >= max) {
    return wrapped_interval<Number>::top();
  }
  wrapped_interval<Number> join = *this | x;
  if (join == wrapped_interval<Number>(m_start, x.m_end)) {
    // increase by some power of 2 the size of the old interval
    wrapint new_end(0, w);
    switch (growth_rate) {
    case 4:
      new_end =
          (m_end * wrapint(4, w)) - (m_start * wrapint(3, w)) + wrapint(3, w);
      break;
    case 8:
      new_end =
          (m_end * wrapint(8, w)) - (m_start * wrapint(7, w)) + wrapint(7, w);
      break;
    case 16:
      new_end = (m_end * wrapint(16, w)) - (m_start * wrapint(15, w)) +
                wrapint(15, w);
      break;
    case 2:
    default:
      new_end = (m_end * wrapint(2, w)) - m_start + wrapint(1, w);
    }
    return join | wrapped_interval<Number>(m_start, new_end);
  } else if (join == wrapped_interval<Number>(x.m_start, m_end)) {
    // decrease by some power of 2 the size of the old interval
    wrapint new_start(0, w);
    switch (growth_rate) {
    case 4:
      new_start =
          (m_start * wrapint(4, w)) - (m_end * wrapint(3, w)) - wrapint(3, w);
      break;
    case 8:
      new_start =
          (m_start * wrapint(8, w)) - (m_end * wrapint(7, w)) - wrapint(7, w);
      break;
    case 16:
      new_start = (m_start * wrapint(16, w)) - (m_end * wrapint(15, w)) -
                  wrapint(15, w);
      break;
    case 2:
    default:
      new_start = (m_start * wrapint(2, w)) - m_end - wrapint(1, w);
    }
    return join | wrapped_interval<Number>(new_start, m_end);
  } else if (x.at(m_start) && x.at(m_end)) {
    // in principle we should increase by some power of two the end
    // point while reducing by the same power of two the start
    // one. We just increase the end and this will eventually reach
    // the start point.
    wrapint delta(0, w);
    switch (growth_rate) {
    case 4:
      delta =
          (m_end * wrapint(4, w)) - (m_start * wrapint(4, w)) + wrapint(3, w);
      break;
    case 8:
      delta =
          (m_end * wrapint(8, w)) - (m_start * wrapint(8, w)) + wrapint(7, w);
      break;
    case 16:
      delta = (m_end * wrapint(16, w)) - (m_start * wrapint(16, w)) +
              wrapint(15, w);
      break;
    case 2:
    default:
      delta =
          (m_end * wrapint(2, w)) - (m_start * wrapint(2, w)) + wrapint(1, w);
    }
    return x | wrapped_interval<Number>(x.m_start, x.m_start + delta);
  } else {
    return wrapped_interval<Number>::top();
  }
}

// TODO: factorize code with operator||
template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::widening_thresholds(
    const wrapped_interval<Number> &x, const thresholds<Number> &ts) const {

  if (is_bottom()) {
    return x;
  } else if (x.is_bottom()) {
    return *this;
  } else if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else if (x <= *this) {
    return *this;
  }
  assert(get_bitwidth(__LINE__) == x.get_bitwidth(__LINE__));
  wrapint::bitwidth_t w = x.get_bitwidth(__LINE__);

  /**
     TODO: make it template parameter.
     growth_rate should be power of 2. Currently only 2,4,8, or 16.
  **/
  const unsigned growth_rate = 8;
  wrapint max(0, w);
  switch (growth_rate) {
  case 4:
    if (w > 2) {
      max = wrapint(1 << (w - 2), w);
      break;
    }
    BOOST_FALLTHROUGH; 
  case 8:
    if (w > 3) {
      max = wrapint(1 << (w - 3), w);
      break;
    }
    BOOST_FALLTHROUGH; 
  case 16:
    if (w > 4) {
      max = wrapint(1 << (w - 4), w);
      break;
    }
    BOOST_FALLTHROUGH; 
  case 2:
    BOOST_FALLTHROUGH; 
  default:
    assert(w > 1);
    max = wrapint(1 << (w - 1), w);
  }
  if ((m_end - m_start) >= max) {
    return wrapped_interval<Number>::top();
  }
  wrapped_interval<Number> join = *this | x;
  if (join == wrapped_interval<Number>(m_start, x.m_end)) {
    // increase by some power of 2 the size of the old interval
    wrapint new_end(0, w);
    switch (growth_rate) {
    case 4:
      new_end =
          (m_end * wrapint(4, w)) - (m_start * wrapint(3, w)) + wrapint(3, w);
      break;
    case 8:
      new_end =
          (m_end * wrapint(8, w)) - (m_start * wrapint(7, w)) + wrapint(7, w);
      break;
    case 16:
      new_end = (m_end * wrapint(16, w)) - (m_start * wrapint(15, w)) +
                wrapint(15, w);
      break;
    case 2:
    default:
      new_end = (m_end * wrapint(2, w)) - m_start + wrapint(1, w);
    }
    // Apply thresholds
    using bound_t = typename ikos::interval<Number>::bound_t;
    bound_t next_end_bound_guess =
        ts.get_next(bound_t(x.m_end.get_unsigned_bignum()));
    if (boost::optional<Number> n = next_end_bound_guess.number()) {
      if (wrapint::fits_wrapint(*n, w)) {
        wrapint new_end_guess(*n, w);
        if (new_end_guess <= new_end) {
          CRAB_LOG("wrapped-int-widening-thresholds",
                   crab::outs()
                       << "Widening with thresholds jumped to " << new_end_guess
                       << " instead of " << new_end << "\n";);
          new_end = new_end_guess;
        }
      }
    }
    return join | wrapped_interval<Number>(m_start, new_end);
  } else if (join == wrapped_interval<Number>(x.m_start, m_end)) {
    // decrease by some power of 2 the size of the old interval
    wrapint new_start(0, w);
    switch (growth_rate) {
    case 4:
      new_start =
          (m_start * wrapint(4, w)) - (m_end * wrapint(3, w)) - wrapint(3, w);
      break;
    case 8:
      new_start =
          (m_start * wrapint(8, w)) - (m_end * wrapint(7, w)) - wrapint(7, w);
      break;
    case 16:
      new_start = (m_start * wrapint(16, w)) - (m_end * wrapint(15, w)) -
                  wrapint(15, w);
      break;
    case 2:
    default:
      new_start = (m_start * wrapint(2, w)) - m_end - wrapint(1, w);
    }
    // TODO: apply thresholds
    return join | wrapped_interval<Number>(new_start, m_end);
  } else if (x.at(m_start) && x.at(m_end)) {
    // in principle we should increase by some power of two the end
    // point while reducing by the same power of two the start
    // one. We just increase the end and this will eventually reach
    // the start point.
    wrapint delta(0, w);
    switch (growth_rate) {
    case 4:
      delta =
          (m_end * wrapint(4, w)) - (m_start * wrapint(4, w)) + wrapint(3, w);
      break;
    case 8:
      delta =
          (m_end * wrapint(8, w)) - (m_start * wrapint(8, w)) + wrapint(7, w);
      break;
    case 16:
      delta = (m_end * wrapint(16, w)) - (m_start * wrapint(16, w)) +
              wrapint(15, w);
      break;
    case 2:
    default:
      delta =
          (m_end * wrapint(2, w)) - (m_start * wrapint(2, w)) + wrapint(1, w);
    }
    // TODO: apply thresholds
    return x | wrapped_interval<Number>(x.m_start, x.m_start + delta);
  } else {
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator&&(const wrapped_interval<Number> &x) const {
  // TODO: for now we call the meet operator.
  return (this->operator&(x));
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator+(const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  } else if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    // -- check if the addition will overflow
    wrapint x_sz = x.m_end - x.m_start;
    wrapint sz = m_end - m_start;
    wrapint one(1, x_sz.get_bitwidth());
    if (x_sz + sz + one <= x_sz) {
      return wrapped_interval<Number>::top();
    }

    return wrapped_interval<Number>(m_start + x.m_start, m_end + x.m_end);
  }
}

template <typename Number>
wrapped_interval<Number> &
wrapped_interval<Number>::operator+=(const wrapped_interval<Number> &x) {
  return this->operator=(this->operator+(x));
}

template <typename Number>
wrapped_interval<Number> wrapped_interval<Number>::operator-() const {
  if (is_bottom()) {
    return wrapped_interval<Number>::bottom();
  } else if (is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    return wrapped_interval<Number>(-m_end, -m_start);
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator-(const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  } else if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    // -- check if the subtraction will overflow
    wrapint x_sz = x.m_end - x.m_start;
    wrapint sz = m_end - m_start;
    wrapint one(1, x_sz.get_bitwidth());
    if (x_sz + sz + one <= x_sz) {
      return wrapped_interval<Number>::top();
    }
    return wrapped_interval<Number>(m_start - x.m_end, m_end - x.m_start);
  }
}

template <typename Number>
wrapped_interval<Number> &
wrapped_interval<Number>::operator-=(const wrapped_interval<Number> &x) {
  return this->operator=(this->operator-(x));
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::operator*(const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  }
  if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    std::vector<wrapped_interval<Number>> cuts, x_cuts;
    signed_and_unsigned_split(cuts);
    x.signed_and_unsigned_split(x_cuts);
    assert(!cuts.empty());
    assert(!x_cuts.empty());

    CRAB_LOG(
        "wrapped-int-mul", crab::outs() << "cuts for " << *this << "\n";
        for (unsigned i = 0, ie = cuts.size(); i < ie;
             ++i) { crab::outs() << "\t" << cuts[i] << "\n"; } crab::outs()
        << "cuts for " << x << "\n";
        for (unsigned i = 0, ie = x_cuts.size(); i < ie;
             ++i) { crab::outs() << "\t" << x_cuts[i] << "\n"; });

    wrapped_interval<Number> res = wrapped_interval<Number>::bottom();
    for (unsigned i = 0, ie = cuts.size(); i < ie; ++i) {
      for (unsigned j = 0, je = x_cuts.size(); j < je; ++j) {
        std::vector<wrapped_interval<Number>> exact_reduct;
        cuts[i].reduced_signed_unsigned_mul(x_cuts[j], exact_reduct);
        for (unsigned k = 0; k < exact_reduct.size(); ++k) {
          res = res | exact_reduct[k];
        }
      }
    }

    CRAB_LOG("wrapped-int-mul",
             crab::outs() << *this << " * " << x << " = " << res << "\n");
    return res;
  }
}

template <typename Number>
wrapped_interval<Number> &
wrapped_interval<Number>::operator*=(const wrapped_interval<Number> &x) {
  return this->operator=(this->operator*(x));
}

/** division and remainder operations **/
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::SDiv(const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  }
  if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    std::vector<wrapped_interval<Number>> cuts, x_cuts;
    signed_and_unsigned_split(cuts);
    x.signed_and_unsigned_split(x_cuts);
    assert(!cuts.empty());
    assert(!x_cuts.empty());
    wrapped_interval<Number> res = wrapped_interval<Number>::bottom();
    for (unsigned i = 0, ie = cuts.size(); i < ie; ++i) {
      for (unsigned j = 0, je = x_cuts.size(); j < je; ++j) {
        std::vector<wrapped_interval<Number>> trimmed_divisors;
        x_cuts[j].trim_zero(trimmed_divisors);
        for (unsigned k = 0, ke = trimmed_divisors.size(); k < ke; ++k) {
          wrapped_interval<Number> d = trimmed_divisors[k];
          res = res | cuts[i].signed_div(d);
        }
      }
    }
    return res;
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::UDiv(const wrapped_interval<Number> &x) const {
  if (is_bottom() || x.is_bottom()) {
    return wrapped_interval<Number>::bottom();
  }
  if (is_top() || x.is_top()) {
    return wrapped_interval<Number>::top();
  } else {
    std::vector<wrapped_interval<Number>> ssplits, x_ssplits;
    signed_split(ssplits);
    x.signed_split(x_ssplits);
    assert(!ssplits.empty());
    assert(!x_ssplits.empty());
    wrapped_interval<Number> res = wrapped_interval<Number>::bottom();
    for (unsigned i = 0, ie = ssplits.size(); i < ie; ++i) {
      for (unsigned j = 0, je = x_ssplits.size(); j < je; ++j) {
        std::vector<wrapped_interval<Number>> trimmed_divisors;
        x_ssplits[j].trim_zero(trimmed_divisors);
        for (unsigned k = 0, ke = trimmed_divisors.size(); k < ke; ++k) {
          wrapped_interval<Number> d = trimmed_divisors[k];
          res = res | ssplits[i].unsigned_div(d);
        }
      }
    }
    return res;
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::SRem(const wrapped_interval<Number> &x) const {
  return default_implementation(x);
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::URem(const wrapped_interval<Number> &x) const {
  return default_implementation(x);
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::ZExt(unsigned bits_to_add) const {
  std::vector<wrapped_interval<Number>> intervals;
  unsigned_split(intervals);

  wrapped_interval<Number> res = wrapped_interval<Number>::bottom();
  for (typename std::vector<wrapped_interval<Number>>::iterator
           it = intervals.begin(),
           et = intervals.end();
       it != et; ++it) {
    // this should not happen
    if ((*it).is_bottom() || (*it).is_top())
      continue;

    wrapint a = (*it).start();
    wrapint b = (*it).end();
    res = res |
          wrapped_interval<Number>(a.zext(bits_to_add), b.zext(bits_to_add));
  }
  return res;
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::SExt(unsigned bits_to_add) const {
  std::vector<wrapped_interval<Number>> intervals;
  signed_split(intervals);

  wrapped_interval<Number> res = wrapped_interval<Number>::bottom();
  for (typename std::vector<wrapped_interval<Number>>::iterator
           it = intervals.begin(),
           et = intervals.end();
       it != et; ++it) {
    // this should not happen
    if ((*it).is_bottom() || (*it).is_top())
      continue;

    wrapint a = (*it).start();
    wrapint b = (*it).end();
    res = res |
          wrapped_interval<Number>(a.sext(bits_to_add), b.sext(bits_to_add));
  }
  return res;
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::Trunc(unsigned bits_to_keep) const {
  if (is_bottom() || is_top()) {
    return *this;
  }

  wrapint::bitwidth_t w = get_bitwidth(__LINE__);
  if (m_start.ashr(wrapint(bits_to_keep, w)) ==
      m_end.ashr(wrapint(bits_to_keep, w))) {
    wrapint lower_start = m_start.keep_lower(bits_to_keep);
    wrapint lower_end = m_end.keep_lower(bits_to_keep);
    if (lower_start <= lower_end) {
      return wrapped_interval<Number>(lower_start, lower_end);
    }
  } else {
    // note that m_start is a wrapint so after +1 it can wraparound
    wrapint y(m_start.ashr(wrapint(bits_to_keep, w)));
    ++y;
    if (y == m_end.ashr(wrapint(bits_to_keep, w))) {
      wrapint lower_start = m_start.keep_lower(bits_to_keep);
      wrapint lower_end = m_end.keep_lower(bits_to_keep);
      if (!(lower_start <= lower_end)) {
        return wrapped_interval<Number>(lower_start, lower_end);
      }
    }
  }
  return wrapped_interval<Number>::top();
}

/// Shl, LShr, and AShr shifts are treated as unsigned numbers
template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::Shl(const wrapped_interval<Number> &x) const {
  if (is_bottom())
    return *this;

  // only if shift is constant
  if (x.is_singleton()) {
    return Shl(x.start().get_uint64_t());
  } else {
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::LShr(const wrapped_interval<Number> &x) const {
  if (is_bottom())
    return *this;

  // only if shift is constant
  if (x.is_singleton()) {
    return LShr(x.start().get_uint64_t());
  } else {
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::AShr(const wrapped_interval<Number> &x) const {
  if (is_bottom()) {
    return *this;
  }

  // only if shift is constant
  if (x.is_singleton()) {
    return AShr(x.start().get_uint64_t());
  } else {
    return wrapped_interval<Number>::top();
  }
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::And(const wrapped_interval<Number> &x) const {
  return default_implementation(x);
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::Or(const wrapped_interval<Number> &x) const {
  return default_implementation(x);
}

template <typename Number>
wrapped_interval<Number>
wrapped_interval<Number>::Xor(const wrapped_interval<Number> &x) const {
  return default_implementation(x);
}

template <typename Number>
void wrapped_interval<Number>::write(crab::crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else if (is_top()) {
    o << "top";
  } else {
#ifdef PRINT_WRAPINT_AS_SIGNED
    // print the wrapints as a signed number (easier to read)
    uint64_t x = m_start.get_uint64_t();
    uint64_t y = m_end.get_uint64_t();
    if (get_bitwidth(__LINE__) == 32) {
      o << "[[" << (int)x << ", " << (int)y << "]]_";
    } else if (get_bitwidth(__LINE__) == 8) {
      o << "[[" << (int)static_cast<signed char>(x) << ", "
        << (int)static_cast<signed char>(y) << "]]_";
    } else {
      o << "[[" << m_start.get_signed_bignum() << ", "
        << m_end.get_signed_bignum() << "]]_";
    }
#else
    o << "[[" << m_start << ", " << m_end << "]]_";
#endif
    o << (int)get_bitwidth(__LINE__);
  }
}

} // namespace domains
} // namespace crab
