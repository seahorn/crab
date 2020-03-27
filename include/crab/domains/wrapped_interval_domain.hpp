#pragma once

/**
 ** Machine arithmetic interval domain based on the paper
 ** "Signedness-Agnostic Program Analysis: Precise Integer Bounds for
 ** Low-Level Code" by J.A.Navas, P.Schachte, H.Sondergaard, and
 ** P.J.Stuckey published in APLAS'12.
 **/

#include <boost/optional.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/numbers/wrapint.hpp>

#define PRINT_WRAPINT_AS_SIGNED

namespace crab {
namespace domains {

template <typename Number> class wrapped_interval {

  wrapint _start;
  wrapint _stop;
  bool _is_bottom;

  typedef wrapped_interval<Number> wrapped_interval_t;

  wrapped_interval_t default_implementation(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } else {
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval(wrapint start, wrapint stop, bool is_bottom)
      : _start(start), _stop(stop), _is_bottom(is_bottom) {
    if (start.get_bitwidth() != stop.get_bitwidth()) {
      CRAB_ERROR("inconsistent bitwidths in wrapped interval");
    }
  }

  // nsplit in the APLAS'12 paper
  void signed_split(std::vector<wrapped_interval_t> &intervals) const {
    if (is_bottom())
      return;

    bitwidth_t b = get_bitwidth(__LINE__);
    if (is_top()) {
      intervals.push_back(wrapped_interval_t(wrapint::get_unsigned_min(b),
                                             wrapint::get_signed_max(b)));
      intervals.push_back(wrapped_interval_t(wrapint::get_signed_min(b),
                                             wrapint::get_unsigned_max(b)));
    } else {
      if (signed_limit(b) <= *this) {
        intervals.push_back(
            wrapped_interval_t(_start, wrapint::get_signed_max(b)));
        intervals.push_back(
            wrapped_interval_t(wrapint::get_signed_min(b), _stop));
      } else {
        intervals.push_back(*this);
      }
    }
  }

  // ssplit in the APLAS'12 paper
  void unsigned_split(std::vector<wrapped_interval_t> &intervals) const {
    if (is_bottom())
      return;

    bitwidth_t b = get_bitwidth(__LINE__);
    if (is_top()) {
      intervals.push_back(wrapped_interval_t(wrapint::get_signed_min(b),
                                             wrapint::get_unsigned_max(b)));
      intervals.push_back(wrapped_interval_t(wrapint::get_unsigned_min(b),
                                             wrapint::get_signed_max(b)));
    } else {
      if (unsigned_limit(b) <= *this) {
        intervals.push_back(
            wrapped_interval_t(_start, wrapint::get_unsigned_max(b)));
        intervals.push_back(
            wrapped_interval_t(wrapint::get_unsigned_min(b), _stop));
      } else {
        intervals.push_back(*this);
      }
    }
  }

  // cut in the APLAS'12 paper
  void signed_and_unsigned_split(std::vector<wrapped_interval_t> &out) const {
    std::vector<wrapped_interval_t> ssplit;
    signed_split(ssplit);
    for (unsigned i = 0, e = ssplit.size(); i < e; ++i) {
      ssplit[i].unsigned_split(out);
    }
  }

  wrapped_interval_t signed_mul(wrapped_interval_t x) const {
    assert(!is_bottom() && !x.is_bottom());

    bool msb_start = _start.msb();
    bool msb_stop = _stop.msb();
    bool msb_x_start = x._start.msb();
    bool msb_x_stop = x._stop.msb();
    bitwidth_t b = get_bitwidth(__LINE__);

    wrapped_interval_t res = wrapped_interval_t::top();
    if (msb_start == msb_stop && msb_stop == msb_x_start &&
        msb_x_start == msb_x_stop) {
      // the two intervals are in the same hemisphere
      if (!msb_start) {
        return unsigned_mul(x);
      } else {
        // -- check if multiplication will overflow
        if ((_start.get_unsigned_bignum() * x._start.get_unsigned_bignum()) -
                (_stop.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) <
            wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
          res = wrapped_interval_t(_stop * x._stop, _start * x._start);
        }
        goto EXIT_SIGNED_MUL;
      }
    }

    // each interval cannot cross the limit: one interval in a
    // different hemisphere.
    if (!(msb_start != msb_stop || msb_x_start != msb_x_stop)) {
      if (msb_start && !msb_x_start) {
        // -- check if multiplication will overflow
        if ((_stop.get_unsigned_bignum() * x._start.get_unsigned_bignum()) -
                (_start.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) <
            wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
          res = wrapped_interval_t(_start * x._stop, _stop * x._start);
        }
      } else if (!msb_start && msb_x_start) {
        // -- check if multiplication will overflow
        if ((_start.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) -
                (_stop.get_unsigned_bignum() * x._start.get_unsigned_bignum()) <
            wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
          res = wrapped_interval_t(_stop * x._start, _start * x._stop);
        }
      }
    }

  EXIT_SIGNED_MUL:
    CRAB_LOG("wrapped-int-mul", crab::outs() << "Signed " << *this << " * " << x
                                             << "=" << res << "\n";);
    return res;
  }

  wrapped_interval_t unsigned_mul(wrapped_interval_t x) const {
    assert(!is_bottom() && !x.is_bottom());

    bitwidth_t b = get_bitwidth(__LINE__);
    wrapped_interval_t res = wrapped_interval_t::top();
    // -- check if multiplication will overflow
    if ((_stop.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) -
            (_start.get_unsigned_bignum() * x._start.get_unsigned_bignum()) <
        wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
      res = wrapped_interval_t(_start * x._start, _stop * x._stop);
    }
    CRAB_LOG("wrapped-int-mul", crab::outs() << "Unsigned " << *this << " * "
                                             << x << "=" << res << "\n";);
    return res;
  }

  // if out is empty then the intersection is empty
  void exact_meet(wrapped_interval_t x,
                  std::vector<wrapped_interval_t> &out) const {
    if (is_bottom() || x.is_bottom()) {
      // bottom
    } else if (*this == x || is_top()) {
      out.push_back(x);
    } else if (x.is_top()) {
      out.push_back(*this);
    } else if (x[_start] &&
               x[_stop] && operator[](x._start) && operator[](x._stop)) {
      out.push_back(wrapped_interval_t(_start, x._stop));
      out.push_back(wrapped_interval_t(x._start, _stop));
    } else if (x[_start] && x[_stop]) {
      out.push_back(*this);
    } else if (operator[](x._start) && operator[](x._stop)) {
      out.push_back(x);
    } else if (x[_start] && operator[](x._stop) &&
               !x[_stop] && operator[](x._start)) {
      out.push_back(wrapped_interval_t(_start, x._stop));
    } else if (x[_stop] && operator[](x._start) &&
               !x[_start] && operator[](x._stop)) {
      out.push_back(wrapped_interval_t(x._start, _stop));
    } else {
      // bottom
    }
  }

  // Perform the reduced product of signed and unsigned multiplication.
  // It uses exact meet rather than abstract meet.
  void reduced_signed_unsigned_mul(wrapped_interval_t x,
                                   std::vector<wrapped_interval_t> &out) const {
    if (is_bottom() || x.is_bottom()) {
      return;
    }

    wrapped_interval_t s = signed_mul(x);
    wrapped_interval_t u = unsigned_mul(x);
    s.exact_meet(u, out);
    CRAB_LOG("wrapped-int-mul", crab::outs() << "Exact signed x unsigned " << s
                                             << " * " << u << "=\n";
             for (unsigned i = 0; i < out.size();
                  ++i) { crab::outs() << "\t" << out[i] << "\n"; });
  }

  wrapped_interval_t unsigned_div(wrapped_interval_t x) const {
    CRAB_LOG("wrapped-int-div", crab::outs() << *this << " /_u " << x << "=";);
    assert(!x[wrapint(0, x.get_bitwidth(__LINE__))]);
    assert(!is_bottom() && !x.is_bottom());
    wrapped_interval_t res =
        wrapped_interval_t(_start.udiv(x._stop), _stop.udiv(x._start));
    CRAB_LOG("wrapped-int-div", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_t signed_div(wrapped_interval_t x) const {
    CRAB_LOG("wrapped-int-div", crab::outs() << *this << " /_s " << x << "=";);
    assert(!x[wrapint(0, x.get_bitwidth(__LINE__))]);
    assert(!is_bottom() && !x.is_bottom());

    bool msb_start = _start.msb();
    bool msb_x_start = x._start.msb();

    bitwidth_t b = get_bitwidth(__LINE__);
    wrapint smin = wrapint::get_signed_min(b);
    wrapint minus_one(-1, b);

    wrapped_interval_t res = wrapped_interval_t::top();
    if (msb_start == msb_x_start) {
      if (msb_start) { // both negative
        // -- check if division will overflow
        if (!((_stop == smin && x._start == minus_one) ||
              (_start == smin && x._stop == minus_one))) {
          res = wrapped_interval_t(_stop.sdiv(x._start), _start.sdiv(x._stop));
        }
      } else { // both positive
        // -- check if division will overflow
        if (!((_start == smin && x._stop == minus_one) ||
              (_stop == smin && x._start == minus_one))) {
          res = wrapped_interval_t(_start.sdiv(x._stop), _stop.sdiv(x._start));
        }
      }
    } else {
      if (msb_start) {
        assert(!msb_x_start);
        // -- check if division will overflow
        if (!((_start == smin && x._start == minus_one) ||
              (_stop == smin && x._stop == minus_one))) {
          res = wrapped_interval_t(_start.sdiv(x._start), _stop.sdiv(x._stop));
        }
      } else {
        assert(msb_x_start);
        // -- check if division will overflow
        if (!((_stop == smin && x._stop == minus_one) ||
              (_start == smin && x._start == minus_one))) {
          res = wrapped_interval_t(_stop.sdiv(x._stop), _start.sdiv(x._start));
        }
      }
    }
    CRAB_LOG("wrapped-int-div", crab::outs() << res << "\n";);
    return res;
  }

  // FIXME: this is sound only if wrapped interval defined over
  //        z_number.
  void trim_zero(std::vector<wrapped_interval_t> &out) const {
    wrapint zero(0, get_bitwidth(__LINE__));
    if (!is_bottom() && (!(*this == zero))) {
      if (start() == zero) {
        out.push_back(
            wrapped_interval_t(wrapint(1, get_bitwidth(__LINE__)), stop()));
      } else if (stop() == zero) {
        out.push_back(
            wrapped_interval_t(start(), wrapint(-1, get_bitwidth(__LINE__))));
      } else if (operator[](zero)) {
        out.push_back(
            wrapped_interval_t(start(), wrapint(-1, get_bitwidth(__LINE__))));
        out.push_back(
            wrapped_interval_t(wrapint(1, get_bitwidth(__LINE__)), stop()));
      } else {
        out.push_back(*this);
      }
    }
  }

  wrapped_interval_t Shl(uint64_t k) const {
    if (is_bottom())
      return *this;

    // XXX: we need the check is_top before calling get_bitwidth();
    if (is_top())
      return *this;

    bitwidth_t b = get_bitwidth(__LINE__);
    wrapped_interval_t y = Trunc(b - k);
    if (!y.is_top()) {
      wrapint wk(k, b);
      return wrapped_interval_t(start() << wk, stop() << wk);
    } else {
      // TODO: return (0, 1..10..0) where #1's = b-k and #0's= k
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t LShr(uint64_t k) const {
    if (is_bottom())
      return *this;

    // XXX: we need the check is_top before cross_signed_limit calls
    // get_bitwidth();
    if (is_top())
      return *this;

    if (!cross_unsigned_limit()) {
      bitwidth_t b = get_bitwidth(__LINE__);
      wrapint wk(k, b);
      return wrapped_interval_t(start().lshr(wk), stop().lshr(wk));
    } else {
      // TODO: return (0, 0..01..1 where #0's = k and #1's = b-k
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t AShr(uint64_t k) const {
    if (is_bottom())
      return *this;

    // XXX: we need the check is_top before cross_signed_limit calls
    // get_bitwidth();
    if (is_top())
      return *this;

    if (!cross_signed_limit()) {
      bitwidth_t b = get_bitwidth(__LINE__);
      wrapint wk(k, b);
      return wrapped_interval_t(start().ashr(wk), stop().ashr(wk));
    } else {
      // TODO: return (1...10..0,        0..01..1)
      //               #1's=k #0's=b-k   #0's=k #1's=b-k
      return wrapped_interval_t::top();
    }
  }

  // // Can raise runtime error if n does not fit into a wrapint
  // wrapped_interval(Number n, bitwidth_t w)
  //   : _start(wrapint(n,w)), _stop(wrapint(n, w)), _is_bottom(false) {}

public:
  typedef wrapint::bitwidth_t bitwidth_t;

  wrapped_interval(wrapint n) : _start(n), _stop(n), _is_bottom(false) {}

  wrapped_interval(wrapint start, wrapint stop)
      : _start(start), _stop(stop), _is_bottom(false) {
    if (start.get_bitwidth() != stop.get_bitwidth()) {
      CRAB_ERROR("inconsistent bitwidths in wrapped interval");
    }
  }

  // To represent top, the particular bitwidth here is irrelevant. We
  // just make sure that _stop - _start == get_max()
  wrapped_interval() : _start(wrapint(0, 3)), _stop(7, 3), _is_bottom(false) {}

  // return top if n does not fit into a wrapint. No runtime errors.
  static wrapped_interval_t mk_winterval(Number n, bitwidth_t width) {
    if (wrapint::fits_wrapint(n, width)) {
      return wrapped_interval_t(wrapint(n, width));
    } else {
      CRAB_WARN(n,
                " does not fit into a wrapint. Returned top wrapped interval");
      return wrapped_interval_t::top();
    }
  }

  // Return top if lb or ub do not fit into a wrapint. No runtime errors.
  static wrapped_interval_t mk_winterval(Number lb, Number ub,
                                         bitwidth_t width) {
    if (!wrapint::fits_wrapint(lb, width)) {
      CRAB_WARN(lb,
                " does not fit into a wrapint. Returned top wrapped interval");
      return wrapped_interval_t::top();
    } else if (!wrapint::fits_wrapint(ub, width)) {
      CRAB_WARN(ub,
                " does not fit into a wrapint. Returned top wrapped interval");
      return wrapped_interval_t::top();
    } else {
      return wrapped_interval_t(wrapint(lb, width), wrapint(ub, width));
    }
  }

  static wrapped_interval_t top() {
    return wrapped_interval_t(wrapint(0, 3), wrapint(7, 3), false);
  }

  static wrapped_interval_t bottom() {
    // the wrapint is irrelevant.
    wrapint i(0, 1);
    return wrapped_interval_t(i, i, true);
  }

  // return interval [0111...1, 1000....0]
  // In the APLAS'12 paper "signed limit" corresponds to "north pole".
  static wrapped_interval_t signed_limit(bitwidth_t b) {
    return wrapped_interval_t(wrapint::get_signed_max(b),
                              wrapint::get_signed_min(b));
  }

  // return interval [1111...1, 0000....0]
  // In the APLAS'12 paper "unsigned limit" corresponds to "south pole".
  static wrapped_interval_t unsigned_limit(bitwidth_t b) {
    return wrapped_interval_t(wrapint::get_unsigned_max(b),
                              wrapint::get_unsigned_min(b));
  }

  bool cross_signed_limit() const {
    return (signed_limit(get_bitwidth(__LINE__)) <= *this);
  }

  bool cross_unsigned_limit() const {
    return (unsigned_limit(get_bitwidth(__LINE__)) <= *this);
  }

  bitwidth_t get_bitwidth(int line) const {
    if (is_bottom()) {
      CRAB_ERROR(
          "get_bitwidth() cannot be called from a bottom element at line ",
          line);
    } else if (is_top()) {
      CRAB_ERROR("get_bitwidth() cannot be called from a top element at line ",
                 line);
    } else {
      assert(_start.get_bitwidth() == _stop.get_bitwidth());
      return _start.get_bitwidth();
    }
  }

  wrapint start() const {
    if (is_top())
      CRAB_ERROR("method start() cannot be called if top");
    return _start;
  }

  wrapint stop() const {
    if (is_top())
      CRAB_ERROR("method start() cannot be called if top");
    return _stop;
  }

  bool is_bottom() const { return _is_bottom; }

  bool is_top() const {
    wrapint maxspan = wrapint::get_unsigned_max(_start.get_bitwidth());
    return (!_is_bottom && (_stop - _start == maxspan));
  }

  // Important: we make the choice here that we interpret wrapint as
  // signed mathematical integers.
  interval<Number> to_interval() const {
    typedef interval<Number> interval_t;
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top() || (cross_signed_limit())) {
      return interval_t::top();
    } else {
      return interval_t(_start.get_signed_bignum(), _stop.get_signed_bignum());
    }
  }

  wrapped_interval_t lower_half_line(bool is_signed) const {
    if (is_top() || is_bottom())
      return *this;

    bitwidth_t b = get_bitwidth(__LINE__);

    if (is_signed) {
      if (this->operator[](wrapint::get_signed_max(b))) {
        return wrapped_interval_t::top();
      }
    } else {
      if (this->operator[](wrapint::get_unsigned_max(b))) {
        return wrapped_interval_t::top();
      }
    }

    wrapint smin = wrapint::get_signed_min(b);
    if (!is_signed) {
      smin = wrapint::get_unsigned_min(b);
    }
    wrapped_interval_t res = wrapped_interval_t(smin, _stop);
    return res;
  }

  wrapped_interval_t upper_half_line(bool is_signed) const {
    if (is_top() || is_bottom())
      return *this;

    bitwidth_t b = get_bitwidth(__LINE__);

    if (is_signed) {
      if (this->operator[](wrapint::get_signed_min(b))) {
        return wrapped_interval_t::top();
      }
    } else {
      if (this->operator[](wrapint::get_unsigned_min(b))) {
        return wrapped_interval_t::top();
      }
    }

    wrapint smax = wrapint::get_signed_max(b);
    if (!is_signed) {
      smax = wrapint::get_unsigned_max(b);
    }
    wrapped_interval_t res = wrapped_interval_t(_start, smax);
    return res;
  }

  bool is_singleton() const {
    return (!is_bottom() && !is_top() && _start == _stop);
  }

  // Starting from _start and going clock-wise x is encountered
  // before than _stop.
  bool operator[](wrapint x) const {
    if (is_bottom()) {
      return false;
    } else if (is_top()) {
      return true;
    } else {
      return ((x - _start) <= (_stop - _start));
    }
  }

  bool operator<=(wrapped_interval_t x) const {
    if (x.is_top() || is_bottom()) {
      return true;
    } else if (x.is_bottom() || is_top()) {
      return false;
    } else if (_start == x._start && _stop == x._stop)
      return true;
    else {
      return x[_start] && x[_stop] &&
             (!(operator[](x._start)) || !(operator[](x._stop)));
    }
  }

  bool operator==(wrapped_interval_t x) const {
    return (*this <= x && x <= *this);
  }

  bool operator!=(wrapped_interval_t x) const { return !(this->operator==(x)); }

  wrapped_interval_t operator|(wrapped_interval_t x) const {
    if (*this <= x) {
      return x;
    } else if (x <= *this) {
      return *this;
    } else {
      if (x[_start] &&
          x[_stop] && operator[](x._start) && operator[](x._stop)) {
        return wrapped_interval_t::top();
      } else if (x[_stop] && operator[](x._start)) {
        return wrapped_interval_t(_start, x._stop);
      } else if (operator[](x._stop) && x[_start]) {
        return wrapped_interval_t(x._start, _stop);
      } else {
        wrapint span_a = x._start - _stop;
        wrapint span_b = _start - x._stop;
        if (span_a < span_b || (span_a == span_b && _start <= x._start)) {
          return wrapped_interval_t(_start, x._stop);
        } else {
          return wrapped_interval_t(x._start, _stop);
        }
      }
    }
  }

  wrapped_interval_t operator&(wrapped_interval_t x) const {
    if (*this <= x) {
      return *this;
    } else if (x <= *this) {
      return x;
    } else {
      if (x[_start]) {
        if (operator[](x._start)) {
          wrapint span_a = _stop - _start;
          wrapint span_b = x._stop - x._start;
          if (span_a < span_b || (span_a == span_b && _start <= x._start)) {
            return *this;
          } else {
            return x; // imprecision here {(_stop, x._stop), (x._start, start)}
          }
        } else {
          if (x[_stop]) {
            return *this;
          } else {
            return wrapped_interval_t(_start, x._stop);
          }
        }
      } else {
        if (operator[](x._start)) {
          if (operator[](x._stop)) {
            return x; // imprecision here {(x._start, _stop), (_start, x._stop)}
          } else {
            return wrapped_interval_t(x._start, _stop);
          }
        } else {
          return wrapped_interval_t::bottom();
        }
      }
    }
  }

  wrapped_interval_t operator||(wrapped_interval_t x) const {
    if (is_bottom()) {
      return x;
    } else if (x.is_bottom()) {
      return *this;
    } else if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else if (x <= *this) {
      return *this;
    }
    assert(get_bitwidth(__LINE__) == x.get_bitwidth(__LINE__));
    bitwidth_t w = x.get_bitwidth(__LINE__);

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
    case 8:
      if (w > 3) {
        max = wrapint(1 << (w - 3), w);
        break;
      }
    case 16:
      if (w > 4) {
        max = wrapint(1 << (w - 4), w);
        break;
      }
    case 2:
    default:
      assert(w > 1);
      max = wrapint(1 << (w - 1), w);
    }
    if ((_stop - _start) >= max) {
      return wrapped_interval_t::top();
    }
    wrapped_interval_t join = *this | x;
    if (join == wrapped_interval_t(_start, x._stop)) {
      // increase by some power of 2 the size of the old interval
      wrapint new_stop(0, w);
      switch (growth_rate) {
      case 4:
        new_stop =
            (_stop * wrapint(4, w)) - (_start * wrapint(3, w)) + wrapint(3, w);
        break;
      case 8:
        new_stop =
            (_stop * wrapint(8, w)) - (_start * wrapint(7, w)) + wrapint(7, w);
        break;
      case 16:
        new_stop = (_stop * wrapint(16, w)) - (_start * wrapint(15, w)) +
                   wrapint(15, w);
        break;
      case 2:
      default:
        new_stop = (_stop * wrapint(2, w)) - _start + wrapint(1, w);
      }
      return join | wrapped_interval_t(_start, new_stop);
    } else if (join == wrapped_interval_t(x._start, _stop)) {
      // decrease by some power of 2 the size of the old interval
      wrapint new_start(0, w);
      switch (growth_rate) {
      case 4:
        new_start =
            (_start * wrapint(4, w)) - (_stop * wrapint(3, w)) - wrapint(3, w);
        break;
      case 8:
        new_start =
            (_start * wrapint(8, w)) - (_stop * wrapint(7, w)) - wrapint(7, w);
        break;
      case 16:
        new_start = (_start * wrapint(16, w)) - (_stop * wrapint(15, w)) -
                    wrapint(15, w);
        break;
      case 2:
      default:
        new_start = (_start * wrapint(2, w)) - _stop - wrapint(1, w);
      }
      return join | wrapped_interval_t(new_start, _stop);
    } else if (x[_start] && x[_stop]) {
      // in principle we should increase by some power of two the stop
      // point while reducing by the same power of two the start
      // one. We just increase the stop and this will eventually reach
      // the start point.
      wrapint delta(0, w);
      switch (growth_rate) {
      case 4:
        delta =
            (_stop * wrapint(4, w)) - (_start * wrapint(4, w)) + wrapint(3, w);
        break;
      case 8:
        delta =
            (_stop * wrapint(8, w)) - (_start * wrapint(8, w)) + wrapint(7, w);
        break;
      case 16:
        delta = (_stop * wrapint(16, w)) - (_start * wrapint(16, w)) +
                wrapint(15, w);
        break;
      case 2:
      default:
        delta =
            (_stop * wrapint(2, w)) - (_start * wrapint(2, w)) + wrapint(1, w);
      }
      return x | wrapped_interval_t(x._start, x._start + delta);
    } else {
      return wrapped_interval_t::top();
    }
  }

  // TODO: factorize code with operator||
  template <typename Thresholds>
  wrapped_interval_t widening_thresholds(wrapped_interval_t x,
                                         const Thresholds &ts) {

    if (is_bottom()) {
      return x;
    } else if (x.is_bottom()) {
      return *this;
    } else if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else if (x <= *this) {
      return *this;
    }
    assert(get_bitwidth(__LINE__) == x.get_bitwidth(__LINE__));
    bitwidth_t w = x.get_bitwidth(__LINE__);

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
    case 8:
      if (w > 3) {
        max = wrapint(1 << (w - 3), w);
        break;
      }
    case 16:
      if (w > 4) {
        max = wrapint(1 << (w - 4), w);
        break;
      }
    case 2:
    default:
      assert(w > 1);
      max = wrapint(1 << (w - 1), w);
    }
    if ((_stop - _start) >= max) {
      return wrapped_interval_t::top();
    }
    wrapped_interval_t join = *this | x;
    if (join == wrapped_interval_t(_start, x._stop)) {
      // increase by some power of 2 the size of the old interval
      wrapint new_stop(0, w);
      switch (growth_rate) {
      case 4:
        new_stop =
            (_stop * wrapint(4, w)) - (_start * wrapint(3, w)) + wrapint(3, w);
        break;
      case 8:
        new_stop =
            (_stop * wrapint(8, w)) - (_start * wrapint(7, w)) + wrapint(7, w);
        break;
      case 16:
        new_stop = (_stop * wrapint(16, w)) - (_start * wrapint(15, w)) +
                   wrapint(15, w);
        break;
      case 2:
      default:
        new_stop = (_stop * wrapint(2, w)) - _start + wrapint(1, w);
      }
      // Apply thresholds
      typedef typename interval<Number>::bound_t bound_t;
      bound_t next_stop_bound_guess =
          ts.get_next(bound_t(x._stop.get_unsigned_bignum()));
      if (boost::optional<Number> n = next_stop_bound_guess.number()) {
        if (wrapint::fits_wrapint(*n, w)) {
          wrapint new_stop_guess(*n, w);
          if (new_stop_guess <= new_stop) {
            CRAB_LOG("wrapped-int-widening-thresholds",
                     crab::outs() << "Widening with thresholds jumped to "
                                  << new_stop_guess << " instead of "
                                  << new_stop << "\n";);
            new_stop = new_stop_guess;
          }
        }
      }
      return join | wrapped_interval_t(_start, new_stop);
    } else if (join == wrapped_interval_t(x._start, _stop)) {
      // decrease by some power of 2 the size of the old interval
      wrapint new_start(0, w);
      switch (growth_rate) {
      case 4:
        new_start =
            (_start * wrapint(4, w)) - (_stop * wrapint(3, w)) - wrapint(3, w);
        break;
      case 8:
        new_start =
            (_start * wrapint(8, w)) - (_stop * wrapint(7, w)) - wrapint(7, w);
        break;
      case 16:
        new_start = (_start * wrapint(16, w)) - (_stop * wrapint(15, w)) -
                    wrapint(15, w);
        break;
      case 2:
      default:
        new_start = (_start * wrapint(2, w)) - _stop - wrapint(1, w);
      }
      // TODO: apply thresholds
      return join | wrapped_interval_t(new_start, _stop);
    } else if (x[_start] && x[_stop]) {
      // in principle we should increase by some power of two the stop
      // point while reducing by the same power of two the start
      // one. We just increase the stop and this will eventually reach
      // the start point.
      wrapint delta(0, w);
      switch (growth_rate) {
      case 4:
        delta =
            (_stop * wrapint(4, w)) - (_start * wrapint(4, w)) + wrapint(3, w);
        break;
      case 8:
        delta =
            (_stop * wrapint(8, w)) - (_start * wrapint(8, w)) + wrapint(7, w);
        break;
      case 16:
        delta = (_stop * wrapint(16, w)) - (_start * wrapint(16, w)) +
                wrapint(15, w);
        break;
      case 2:
      default:
        delta =
            (_stop * wrapint(2, w)) - (_start * wrapint(2, w)) + wrapint(1, w);
      }
      // TODO: apply thresholds
      return x | wrapped_interval_t(x._start, x._start + delta);
    } else {
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t operator&&(wrapped_interval_t x) const {
    // TODO: for now we call the meet operator.
    return (this->operator&(x));
  }

  wrapped_interval_t operator+(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } else if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      // -- check if the addition will overflow
      wrapint x_sz = x._stop - x._start;
      wrapint sz = _stop - _start;
      wrapint one(1, x_sz.get_bitwidth());
      if (x_sz + sz + one <= x_sz) {
        return wrapped_interval_t::top();
      }

      return wrapped_interval_t(_start + x._start, _stop + x._stop);
    }
  }

  wrapped_interval_t &operator+=(wrapped_interval_t x) {
    return this->operator=(this->operator+(x));
  }

  wrapped_interval_t operator-() const {
    if (is_bottom()) {
      return wrapped_interval_t::bottom();
    } else if (is_top()) {
      return wrapped_interval_t::top();
    } else {
      return wrapped_interval_t(-_stop, -_start);
    }
  }

  wrapped_interval_t operator-(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } else if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      // -- check if the subtraction will overflow
      wrapint x_sz = x._stop - x._start;
      wrapint sz = _stop - _start;
      wrapint one(1, x_sz.get_bitwidth());
      if (x_sz + sz + one <= x_sz) {
        return wrapped_interval_t::top();
      }
      return wrapped_interval_t(_start - x._stop, _stop - x._start);
    }
  }

  wrapped_interval_t &operator-=(wrapped_interval_t x) {
    return this->operator=(this->operator-(x));
  }

  wrapped_interval_t operator*(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    }
    if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> cuts, x_cuts;
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

      wrapped_interval_t res = wrapped_interval_t::bottom();
      for (unsigned i = 0, ie = cuts.size(); i < ie; ++i) {
        for (unsigned j = 0, je = x_cuts.size(); j < je; ++j) {
          std::vector<wrapped_interval_t> exact_reduct;
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

  wrapped_interval_t &operator*=(wrapped_interval_t x) {
    return this->operator=(this->operator*(x));
  }

  // signed division
  wrapped_interval_t operator/(wrapped_interval_t x) const { return SDiv(x); }

  wrapped_interval_t &operator/=(wrapped_interval_t x) {
    return this->operator=(this->operator/(x));
  }

  /** division and remainder operations **/
  wrapped_interval_t SDiv(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    }
    if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> cuts, x_cuts;
      signed_and_unsigned_split(cuts);
      x.signed_and_unsigned_split(x_cuts);
      assert(!cuts.empty());
      assert(!x_cuts.empty());
      wrapped_interval_t res = wrapped_interval_t::bottom();
      for (unsigned i = 0, ie = cuts.size(); i < ie; ++i) {
        for (unsigned j = 0, je = x_cuts.size(); j < je; ++j) {
          std::vector<wrapped_interval_t> trimmed_divisors;
          x_cuts[j].trim_zero(trimmed_divisors);
          for (unsigned k = 0, ke = trimmed_divisors.size(); k < ke; ++k) {
            wrapped_interval_t d = trimmed_divisors[k];
            res = res | cuts[i].signed_div(d);
          }
        }
      }
      return res;
    }
  }

  wrapped_interval_t UDiv(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    }
    if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> ssplits, x_ssplits;
      signed_split(ssplits);
      x.signed_split(x_ssplits);
      assert(!ssplits.empty());
      assert(!x_ssplits.empty());
      wrapped_interval_t res = wrapped_interval_t::bottom();
      for (unsigned i = 0, ie = ssplits.size(); i < ie; ++i) {
        for (unsigned j = 0, je = x_ssplits.size(); j < je; ++j) {
          std::vector<wrapped_interval_t> trimmed_divisors;
          x_ssplits[j].trim_zero(trimmed_divisors);
          for (unsigned k = 0, ke = trimmed_divisors.size(); k < ke; ++k) {
            wrapped_interval_t d = trimmed_divisors[k];
            res = res | ssplits[i].unsigned_div(d);
          }
        }
      }
      return res;
    }
  }

  wrapped_interval_t SRem(wrapped_interval_t x) const {
    return default_implementation(x);
  }

  wrapped_interval_t URem(wrapped_interval_t x) const {
    return default_implementation(x);
  }

  /** cast operations **/
  wrapped_interval_t ZExt(unsigned bits_to_add) const {
    std::vector<wrapped_interval_t> intervals;
    unsigned_split(intervals);

    wrapped_interval_t res = wrapped_interval_t::bottom();
    for (typename std::vector<wrapped_interval_t>::iterator
             it = intervals.begin(),
             et = intervals.end();
         it != et; ++it) {
      // this should not happen
      if ((*it).is_bottom() || (*it).is_top())
        continue;

      wrapint a = (*it).start();
      wrapint b = (*it).stop();
      res = res | wrapped_interval_t(a.zext(bits_to_add), b.zext(bits_to_add));
    }
    return res;
  }

  wrapped_interval_t SExt(unsigned bits_to_add) const {
    std::vector<wrapped_interval_t> intervals;
    signed_split(intervals);

    wrapped_interval_t res = wrapped_interval_t::bottom();
    for (typename std::vector<wrapped_interval_t>::iterator
             it = intervals.begin(),
             et = intervals.end();
         it != et; ++it) {
      // this should not happen
      if ((*it).is_bottom() || (*it).is_top())
        continue;

      wrapint a = (*it).start();
      wrapint b = (*it).stop();
      res = res | wrapped_interval_t(a.sext(bits_to_add), b.sext(bits_to_add));
    }
    return res;
  }

  wrapped_interval_t Trunc(unsigned bits_to_keep) const {
    if (is_bottom() || is_top()) {
      return *this;
    }

    bitwidth_t w = get_bitwidth(__LINE__);
    if (_start.ashr(wrapint(bits_to_keep, w)) ==
        _stop.ashr(wrapint(bits_to_keep, w))) {
      wrapint lower_start = _start.keep_lower(bits_to_keep);
      wrapint lower_stop = _stop.keep_lower(bits_to_keep);
      if (lower_start <= lower_stop) {
        return wrapped_interval_t(lower_start, lower_stop);
      }
    } else {
      // note that _start is a wrapint so after +1 it can wraparound
      wrapint y(_start.ashr(wrapint(bits_to_keep, w)));
      ++y;
      if (y == _stop.ashr(wrapint(bits_to_keep, w))) {
        wrapint lower_start = _start.keep_lower(bits_to_keep);
        wrapint lower_stop = _stop.keep_lower(bits_to_keep);
        if (!(lower_start <= lower_stop)) {
          return wrapped_interval_t(lower_start, lower_stop);
        }
      }
    }
    return wrapped_interval_t::top();
  }

  /** bitwise operations **/
  // Shl, LShr, and AShr shifts are treated as unsigned numbers
  wrapped_interval_t Shl(wrapped_interval_t x) const {
    if (is_bottom())
      return *this;

    // only if shift is constant
    if (x.is_singleton()) {
      return Shl(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t LShr(wrapped_interval_t x) const {
    if (is_bottom())
      return *this;

    // only if shift is constant
    if (x.is_singleton()) {
      return LShr(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t AShr(wrapped_interval_t x) const {
    if (is_bottom())
      return *this;

    // only if shift is constant
    if (x.is_singleton()) {
      return AShr(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }

  wrapped_interval_t And(wrapped_interval_t x) const {
    return default_implementation(x);
  }

  wrapped_interval_t Or(wrapped_interval_t x) const {
    return default_implementation(x);
  }

  wrapped_interval_t Xor(wrapped_interval_t x) const {
    return default_implementation(x);
  }

  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
#ifdef PRINT_WRAPINT_AS_SIGNED
      // print the wrapints as a signed number (easier to read)
      uint64_t x = _start.get_uint64_t();
      uint64_t y = _stop.get_uint64_t();
      if (get_bitwidth(__LINE__) == 32) {
        o << "[[" << (int)x << ", " << (int)y << "]]_";
      } else if (get_bitwidth(__LINE__) == 8) {
        o << "[[" << (int)static_cast<signed char>(x) << ", "
          << (int)static_cast<signed char>(y) << "]]_";
      } else {
        o << "[[" << _start.get_signed_bignum() << ", "
          << _stop.get_signed_bignum() << "]]_";
      }
#else
      o << "[[" << _start << ", " << _stop << "]]_";
#endif
      o << (int)get_bitwidth(__LINE__);
    }
  }
};

template <typename Number>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const wrapped_interval<Number> &i) {
  i.write(o);
  return o;
}
} // namespace domains
} // namespace crab

namespace ikos {
namespace linear_interval_solver_impl {
typedef crab::domains::wrapped_interval<z_number> z_wrapped_interval_t;
typedef crab::domains::wrapped_interval<q_number> q_wrapped_interval_t;

template <>
inline z_wrapped_interval_t mk_interval(z_number c,
                                        typename crab::wrapint::bitwidth_t w) {
  return z_wrapped_interval_t::mk_winterval(c, w);
}

template <>
inline q_wrapped_interval_t mk_interval(q_number c,
                                        typename crab::wrapint::bitwidth_t w) {
  return q_wrapped_interval_t::mk_winterval(c, w);
}

template <>
inline z_wrapped_interval_t trim_interval(z_wrapped_interval_t i,
                                          z_wrapped_interval_t j) {
  if (i.is_bottom())
    return i;
  // XXX: TODO: gamma(top()) \ gamma(j) can be expressed in a
  //            wrapped interval.
  if (i.is_top())
    return i;
  if (!j.is_singleton())
    return i;

  crab::wrapint k = j.start();
  if (i.start() == k) {
    if (i.is_singleton()) {
      return z_wrapped_interval_t::bottom();
    }
    crab::wrapint k_plus(k);
    ++k_plus;
    z_wrapped_interval_t trimmed_res = z_wrapped_interval_t(k_plus, i.stop());
    return trimmed_res;
  } else if (i.stop() == k) {
    if (i.is_singleton()) {
      return z_wrapped_interval_t::bottom();
    }
    crab::wrapint k_minus(k);
    --k_minus;
    z_wrapped_interval_t trimmed_res = z_wrapped_interval_t(i.start(), k_minus);
    return trimmed_res;
  } else {
    return i;
  }
}

template <>
inline q_wrapped_interval_t trim_interval(q_wrapped_interval_t i,
                                          q_wrapped_interval_t j) {
  // No refinement possible for disequations over rational numbers
  return i;
}

template <>
inline z_wrapped_interval_t lower_half_line(z_wrapped_interval_t i,
                                            bool is_signed) {
  return i.lower_half_line(is_signed);
}

template <>
inline q_wrapped_interval_t lower_half_line(q_wrapped_interval_t i,
                                            bool is_signed) {
  return i.lower_half_line(is_signed);
}

template <>
inline z_wrapped_interval_t upper_half_line(z_wrapped_interval_t i,
                                            bool is_signed) {
  return i.upper_half_line(is_signed);
}

template <>
inline q_wrapped_interval_t upper_half_line(q_wrapped_interval_t i,
                                            bool is_signed) {
  return i.upper_half_line(is_signed);
}

} // namespace linear_interval_solver_impl
} // end namespace ikos

namespace crab {
namespace domains {

template <typename Number, typename VariableName,
          std::size_t max_reduction_cycles = 10>
class wrapped_interval_domain final
    : public abstract_domain<
          wrapped_interval_domain<Number, VariableName, max_reduction_cycles>> {
  typedef wrapped_interval_domain<Number, VariableName, max_reduction_cycles>
      wrapped_interval_domain_t;
  typedef crab::domains::abstract_domain<wrapped_interval_domain_t>
      abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef interval<number_t> interval_t;
  typedef wrapped_interval<number_t> wrapped_interval_t;
  typedef typename wrapped_interval_t::bitwidth_t bitwidth_t;

private:
  typedef separate_domain<variable_t, wrapped_interval_t> separate_domain_t;
  typedef linear_interval_solver<number_t, varname_t, separate_domain_t>
      solver_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  wrapped_interval_domain(separate_domain_t env) : _env(env) {}

  void add(linear_constraint_system_t csts,
           std::size_t threshold = max_reduction_cycles) {
    if (!this->is_bottom()) {
      solver_t solver(csts, threshold);
      solver.run(this->_env);
    }
  }

  wrapped_interval_t eval_expr(linear_expression_t expr, bitwidth_t width) {
    wrapped_interval_t r =
        wrapped_interval_t::mk_winterval(expr.constant(), width);
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      wrapped_interval_t c = wrapped_interval_t::mk_winterval(it->first, width);
      // eval_expr should be "const" but operator[] in _env is not marked as
      // "const"
      r += c * this->_env[it->second];
    }
    return r;
  }

public:
  void set_to_top() {
    wrapped_interval_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    wrapped_interval_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  wrapped_interval_domain() : _env(separate_domain_t::top()) {}

  wrapped_interval_domain(const wrapped_interval_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  wrapped_interval_domain_t &operator=(const wrapped_interval_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }

  iterator begin() { return this->_env.begin(); }

  iterator end() { return this->_env.end(); }

  bool is_bottom() { return this->_env.is_bottom(); }

  bool is_top() { return this->_env.is_top(); }

  bool operator<=(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    // CRAB_LOG("wrapped-int",
    //       crab::outs()<< *this << " <= " << e << "=";);
    bool res = (this->_env <= e._env);
    // CRAB_LOG("wrapped-int",
    //	     crab::outs() << (res ? "yes": "not") << "\n";);
    return res;
  }

  void operator|=(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    CRAB_LOG("wrapped-int", crab::outs() << *this << " U " << e << " = ");
    this->_env = this->_env | e._env;
    CRAB_LOG("wrapped-int", crab::outs() << *this << "\n";);
  }

  wrapped_interval_domain_t operator|(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    CRAB_LOG("wrapped-int", crab::outs() << *this << " U " << e << " = ");
    wrapped_interval_domain_t res(this->_env | e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t operator&(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");
    CRAB_LOG("wrapped-int", crab::outs() << *this << " n " << e << " = ");
    wrapped_interval_domain_t res(this->_env & e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t operator||(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    CRAB_LOG("wrapped-int",
             crab::outs() << "WIDENING " << *this << " and " << e << " = ");
    wrapped_interval_domain_t res(this->_env || e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t
  widening_thresholds(wrapped_interval_domain_t e,
                      const iterators::thresholds<number_t> &ts) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    CRAB_LOG("wrapped-int",
             crab::outs() << "WIDENING " << *this << " and " << e << " = ");
    wrapped_interval_domain_t res(this->_env.widening_thresholds(e._env, ts));
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t operator&&(wrapped_interval_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    return (this->_env && e._env);
  }

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
    this->_env -= v;
  }

  void set(variable_t v, wrapped_interval_t i) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, i);
    CRAB_LOG("wrapped-int", crab::outs()
                                << v << ":=" << i << "=" << _env[v] << "\n");
  }

  void set(variable_t v, interval_t i) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    if (i.lb().is_finite() && i.ub.is_finite()) {
      wrapped_interval_t rhs =
          wrapped_interval_t::mk_winterval(i.lb(), i.ub(), v.get_bitwidth());
      this->_env.set(v, rhs);
      CRAB_LOG("wrapped-int", crab::outs()
                                  << v << ":=" << i << "=" << _env[v] << "\n");
    } else {
      CRAB_WARN(
          "ignored assignment of an open interval in wrapped interval domain");
      *this -= v;
    }
  }

  void set(variable_t v, number_t n) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, wrapped_interval_t::mk_winterval(n, v.get_bitwidth()));
    CRAB_LOG("wrapped-int", crab::outs()
                                << v << ":=" << n << "=" << _env[v] << "\n");
  }

  // Return unlimited interval
  interval_t operator[](variable_t v) const {
    wrapped_interval_t w_i = this->_env[v];
    return w_i.to_interval();
  }

  // Return wrapped interval
  wrapped_interval_t get_wrapped_interval(variable_t v) const {
    return this->_env[v];
  }

  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    if (boost::optional<variable_t> v = e.get_variable()) {
      this->_env.set(x, this->_env[*v]);
    } else {
      wrapped_interval_t r = eval_expr(e, x.get_bitwidth());
      this->_env.set(x, r);
    }
    CRAB_LOG("wrapped-int", crab::outs()
                                << x << ":=" << e << "=" << _env[x] << "\n");
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi = this->_env[z];
    wrapped_interval_t xi = wrapped_interval_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
    CRAB_LOG("wrapped-int", crab::outs() << x << ":=" << y << " " << op << " "
                                         << z << "=" << _env[x] << "\n");
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi =
        wrapped_interval_t::mk_winterval(k, x.get_bitwidth());
    wrapped_interval_t xi = wrapped_interval_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
    CRAB_LOG("wrapped-int", crab::outs() << x << ":=" << y << " " << op << " "
                                         << k << "=" << _env[x] << "\n");
  }

  // cast operations

  void apply(crab::domains::int_conv_operation_t op, variable_t dst,
             variable_t src) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    wrapped_interval_t src_i = this->_env[src];
    wrapped_interval_t dst_i;

    if (src_i.is_bottom() || src_i.is_top()) {
      dst_i = src_i;
    } else {
      switch (op) {
      case crab::domains::OP_ZEXT:
      case crab::domains::OP_SEXT: {
        if (dst.get_bitwidth() < src.get_bitwidth()) {
          CRAB_ERROR("destination must be larger than source in sext/zext");
        }
        unsigned bits_to_add = dst.get_bitwidth() - src.get_bitwidth();
        dst_i = (op == crab::domains::OP_SEXT ? src_i.SExt(bits_to_add)
                                              : src_i.ZExt(bits_to_add));
      } break;
      case crab::domains::OP_TRUNC: {
        if (src.get_bitwidth() < dst.get_bitwidth()) {
          CRAB_ERROR("destination must be smaller than source in truncate");
        }
        unsigned bits_to_keep = dst.get_bitwidth();
        wrapped_interval_t dst_i;
        dst_i = src_i.Trunc(bits_to_keep);
      } break;
      default:
        CRAB_ERROR("unexpected operation: ", op);
      }
    }
    set(dst, dst_i);
  }

  // bitwise operations

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi = this->_env[z];
    wrapped_interval_t xi = wrapped_interval_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi =
        wrapped_interval_t::mk_winterval(k, x.get_bitwidth());
    wrapped_interval_t xi = wrapped_interval_t::bottom();
    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }

  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");
    this->add(csts);
    CRAB_LOG("wrapped-int", crab::outs()
                                << "Added " << csts << " = " << *this << "\n");
  }

  wrapped_interval_domain_t operator+(linear_constraint_system_t csts) {
    wrapped_interval_domain_t e(this->_env);
    e += csts;
    return e;
  }

  // backward arithmetic operations
  void backward_assign(variable_t x, linear_expression_t e,
                       wrapped_interval_domain_t inv) {
    this->operator-=(x);
    CRAB_WARN("Backward assign for wrapped intervals not implemented");
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      wrapped_interval_domain_t inv) {
    this->operator-=(x);
    CRAB_WARN("Backward apply for wrapped intervals not implemented");
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      wrapped_interval_domain_t inv) {
    this->operator-=(x);
    CRAB_WARN("Backward apply for wrapped intervals not implemented");
  }

  /*
     Begin unimplemented operations.

     wrapped_interval_domain implements only standard abstract operations of a
     numerical domain.  The implementation of boolean, array, or
     pointer operations is empty because they should never be
     called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                wrapped_interval_domain_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                wrapped_interval_domain_t invariant) {}
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  wrapped_interval_domain_t invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           wrapped_interval_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           wrapped_interval_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            wrapped_interval_domain_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            wrapped_interval_domain_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  wrapped_interval_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  wrapped_interval_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             wrapped_interval_domain_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs) {}
  void pointer_store(variable_t lhs, variable_t rhs) {}
  void pointer_assign(variable_t lhs, variable_t rhs,
                      linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    separate_domain_t projected_env;
    for (variable_t v : variables) {
      projected_env.set(v, this->_env[v]);
    }
    std::swap(this->_env, projected_env);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set(new_x, this->_env[x]);
  }

  void normalize() {}

  void minimize() {}

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    _env.rename(from, to);            
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  wrapped_interval_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void write(crab::crab_os &o) {
    crab::CrabStats::count(getDomainName() + ".count.write");
    crab::ScopedCrabStats __st__(getDomainName() + ".write");

    this->_env.write(o);
  }

  // Important: we make the choice here that we interpret wrapint as
  // signed mathematical integers.
  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (this->is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      variable_t v = it->first;
      wrapped_interval_t i = it->second;
      if (!i.is_top() && !i.cross_signed_limit()) {
        csts += linear_constraint_t(v >= i.start().get_signed_bignum());
        csts += linear_constraint_t(v <= i.stop().get_signed_bignum());
      }
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  static std::string getDomainName() { return "WrappedIntervals"; }

}; // class wrapped_interval_domain

template <typename Number, typename VariableName>
struct abstract_domain_traits<wrapped_interval_domain<Number, VariableName>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

template <typename Number, typename VariableName>
class constraint_simp_domain_traits<
    wrapped_interval_domain<Number, VariableName>> {
public:
  typedef ikos::linear_constraint<Number, VariableName> linear_constraint_t;
  typedef ikos::linear_constraint_system<Number, VariableName>
      linear_constraint_system_t;

  static void lower_equality(linear_constraint_t cst,
                             linear_constraint_system_t &csts) {
    // We cannot convert an equality into inequalities because we
    // don't know the interpretation (signed/unsigned) for those
    // inequalities.
    csts += cst;
  }
};

} // namespace domains
} // namespace crab

/*
 *
 EXPERIMENTAL CODE: USE IT ON YOUR OWN RISK!
 *
 */
namespace crab {
namespace domains {

// Simple lattice to represent which limits (if any) have been crossed
// by a wrapped interval.
class wrapped_interval_limit_value {
  /*
                    csu
                    /  \
                  cs    cu
                   \    /
                     nc
                     |
                   bottom

   nc: no cross either signed or unsigned limits.
   cs: cross signed limit.
   cu: cross unsigned limit.
   csu: cross both signed and unsigned limits.

   where signed limit   is the interval [0111...1, 1000....0]
         unsigned limit is the interval [1111...1, 0000....0]
  */

  // bottom is left outside intentionally so the join (meet) is simply
  // bitwise-or (and).
  typedef enum {
    NC = 0x0,
    CS = 0x1,
    CU = 0x2,
    CSU = 0x3 /*top*/
  } kind_t;

  kind_t _value;
  bool _is_bottom;

  wrapped_interval_limit_value(kind_t v, bool is_bottom)
      : _value(v), _is_bottom(is_bottom) {}

public:
  wrapped_interval_limit_value() : _value(CSU), _is_bottom(false) {}

  static wrapped_interval_limit_value bottom() {
    return wrapped_interval_limit_value(NC /*any value*/, true);
  }

  static wrapped_interval_limit_value top() {
    return wrapped_interval_limit_value(CSU, false);
  }

  static wrapped_interval_limit_value cross_signed_limit() {
    return wrapped_interval_limit_value(CS, false);
  }

  static wrapped_interval_limit_value cross_unsigned_limit() {
    return wrapped_interval_limit_value(CU, false);
  }

  static wrapped_interval_limit_value do_not_cross() {
    return wrapped_interval_limit_value(NC, false);
  }

  template <typename N>
  static wrapped_interval_limit_value convert(const wrapped_interval<N> &i) {
    if (i.is_bottom()) {
      return wrapped_interval_limit_value::bottom();
    } else if (i.is_top()) {
      return wrapped_interval_limit_value::top();
    } else if (i.cross_unsigned_limit()) {
      return wrapped_interval_limit_value::cross_unsigned_limit();
    } else if (i.cross_signed_limit()) {
      return wrapped_interval_limit_value::cross_signed_limit();
    } else {
      return wrapped_interval_limit_value::do_not_cross();
    }
  }

  wrapped_interval_limit_value(const wrapped_interval_limit_value &o)
      : _value(o._value), _is_bottom(o._is_bottom) {}

  wrapped_interval_limit_value &
  operator=(const wrapped_interval_limit_value &o) {
    if (this != &o) {
      _value = o._value;
      _is_bottom = o._is_bottom;
    }
    return *this;
  }

  bool is_bottom() const { return _is_bottom; }

  // the wrapped interval might have crossed both limits
  bool is_top() const { return !_is_bottom && _value == CSU; }

  // the wrapped interval might have crossed the signed limit
  bool is_crossing_signed_limit() const {
    return (!is_bottom() && (_value == CS || _value == CSU));
  }

  // the wrapped interval might have crossed the unsigned limit
  bool is_crossing_unsigned_limit() const {
    return (!is_bottom() && (_value == CU || _value == CSU));
  }

  bool is_not_crossing_limit() const {
    return (!is_bottom() && (_value == NC));
  }

  bool operator<=(const wrapped_interval_limit_value &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (o.is_bottom()) {
      return false;
    } else if (is_top()) {
      return o.is_top();
    } else if (_value == NC) {
      return true;
    } else if (_value == CS) {
      return (o._value == CS || o.is_top());
    } else if (_value == CU) {
      return (o._value == CU || o.is_top());
    }

    assert(false && "unreachable");
    return false;
  }

  bool operator==(const wrapped_interval_limit_value &o) const {
    return (_value == o._value && is_bottom() == o.is_bottom());
  }

  wrapped_interval_limit_value operator|(wrapped_interval_limit_value o) {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      return wrapped_interval_limit_value(
          static_cast<kind_t>(static_cast<int>(_value) |
                              static_cast<int>(o._value)),
          false);
    }
  }

  // the lattice satisfy ACC so join is the widening
  wrapped_interval_limit_value operator||(wrapped_interval_limit_value o) {
    return this->operator|(o);
  }

  wrapped_interval_limit_value operator&(wrapped_interval_limit_value o) {
    if (is_bottom() || o.is_bottom()) {
      return bottom();
    } else {
      return wrapped_interval_limit_value(
          static_cast<kind_t>(static_cast<int>(_value) &
                              static_cast<int>(o._value)),
          false);
    }
  }

  // the lattice satisfy DCC so meet is the narrowing
  wrapped_interval_limit_value operator&&(wrapped_interval_limit_value o) {
    return this->operator&(o);
  }

  void write(crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else {
      switch (_value) {
      case NC:
        o << "no-cross";
        break;
      case CS:
        o << "cross-signed";
        break;
      case CU:
        o << "cross-unsigned";
        break;
      default: /*top*/
        o << "top";
      }
    }
  }
};

inline crab_os &operator<<(crab_os &o, const wrapped_interval_limit_value &v) {
  v.write(o);
  return o;
}

/**
    Wrapped interval domain augmented with an abstraction of the
    execution history: it keeps track of which variable crossed which
    signed/unsigned limits.

    The only case where the history of a variable is reset is when it
    is assigned to a constant value.
**/
template <typename Number, typename VariableName,
          std::size_t max_reduction_cycles = 10>
class wrapped_interval_with_history_domain final
    : public abstract_domain<wrapped_interval_with_history_domain<
          Number, VariableName, max_reduction_cycles>> {
  typedef wrapped_interval_with_history_domain<Number, VariableName,
                                               max_reduction_cycles>
      this_type;
  typedef abstract_domain<this_type> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef interval<number_t> interval_t;
  typedef wrapped_interval<number_t> wrapped_interval_t;
  typedef wrapped_interval_domain<number_t, varname_t, max_reduction_cycles>
      wrapped_interval_domain_t;

private:
  typedef separate_domain<variable_t, wrapped_interval_limit_value>
      separate_domain_t;
  typedef discrete_domain<variable_t> discrete_domain_t;
  typedef typename linear_constraint_system_t::variable_set_t variable_set_t;

  wrapped_interval_domain_t _w_int_dom;
  // Map each variable to which limit was crossed.
  separate_domain_t _limit_env;
  // Set of may-initialized variables
  discrete_domain_t _init_set;

  wrapped_interval_with_history_domain(const wrapped_interval_domain_t &dom,
                                       const separate_domain_t &limit_env,
                                       const discrete_domain_t &init_set)
      : _w_int_dom(dom), _limit_env(limit_env), _init_set(init_set) {}

  inline bool may_be_initialized(variable_t x) {
    return (discrete_domain_t(x) <= _init_set);
  }

  // Decide whether x might cross a pole based on the intervals before
  // and after an operation occurred.
  inline void update_limits(variable_t x, wrapped_interval_t old_i,
                            wrapped_interval_t new_i) {
    if (is_bottom()) {
      return;
    }

    if (may_be_initialized(x) && (old_i.is_top() && !new_i.is_top())) {
      // This tries to capture the following pattern:
      //  bb:  x:=0; goto bb1;
      //  bb1: ...   goto bb2;
      //  bb2: y:=x+1; x:=y; goto bb1;
      //
      //  y is alive only in bb2. The first time we analyze y:=x+1, y
      //  is uninitialized. However, the second time, the interval for
      //  y before the assignment is top since at bb1 we joined bb
      //  (where y is not defined and hence top) and bb2 (where y is
      //  precisely captured). Because of this, we cannot keep track
      //  of y.
      _limit_env.set(x, wrapped_interval_limit_value::convert(new_i));
      CRAB_LOG("wrapped-int-hist", auto v = _limit_env[x];
               crab::outs()
               << x
               << " may be initialized, old val=top,  and new val != top) = "
               << v << "\n";);
      return;
    }

    wrapped_interval_limit_value old_l = wrapped_interval_limit_value::bottom();
    wrapped_interval_limit_value new_l;
    if (may_be_initialized(x)) {
      old_l = _limit_env[x];
      // XXX: it's not enough to convert only new_i. E.g., char x = 127; x++;
      //  convert([127,127])   = no-cross
      //  convert([-128,-128]) = no-cross
      new_l = wrapped_interval_limit_value::convert(old_i | new_i);
      CRAB_LOG("wrapped-int-hist", crab::outs()
                                       << x << " may be initialized. "
                                       << "old val=" << old_l
                                       << " U new val=" << new_l << "\n";);
    } else {
      CRAB_LOG("wrapped-int-hist", crab::outs()
                                       << x << " is not initialized. "
                                       << "new val=" << new_l << "\n";);
      new_l = wrapped_interval_limit_value::convert(new_i);
    }

    // -- weak update to keep past history
    _limit_env.set(x, old_l | new_l);
  }

  void update_limits(const variable_set_t &vars,
                     const std::vector<wrapped_interval_t> &old_intervals,
                     const std::vector<wrapped_interval_t> &new_intervals) {
    assert(vars.size() == old_intervals.size());
    assert(old_intervals.size() == new_intervals.size());

    unsigned i = 0;
    for (auto const &v : vars) {
      update_limits(v, old_intervals[i], new_intervals[i]);
      ++i;
    }
  }

public:
  void set_to_top() {
    this_type abs(wrapped_interval_domain_t::top(), separate_domain_t::top(),
                  discrete_domain_t::bottom() /*empty set*/);
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    this_type abs(wrapped_interval_domain_t::bottom(),
                  separate_domain_t::bottom(),
                  discrete_domain_t::bottom() /*empty set*/);
    std::swap(*this, abs);
  }

  wrapped_interval_with_history_domain()
      : _w_int_dom(), _limit_env(),
        _init_set(discrete_domain_t::bottom() /*empty set*/) {}

  wrapped_interval_with_history_domain(const this_type &o)
      : _w_int_dom(o._w_int_dom), _limit_env(o._limit_env),
        _init_set(o._init_set) {}

  wrapped_interval_with_history_domain(const this_type &&o)
      : _w_int_dom(std::move(o._w_int_dom)),
        _limit_env(std::move(o._limit_env)), _init_set(std::move(o._init_set)) {
  }

  this_type &operator=(const this_type &o) {
    if (this != &o) {
      _w_int_dom = o._w_int_dom;
      _limit_env = o._limit_env;
      _init_set = o._init_set;
    }
    return *this;
  }

  bool is_bottom() {
    // XXX: ignore _limit_env
    return _w_int_dom.is_bottom();
  }

  bool is_top() {
    // XXX: ignore _limit_env
    return _w_int_dom.is_top();
  }

  bool operator<=(this_type o) {
    return (_w_int_dom <= o._w_int_dom && _limit_env <= o._limit_env);
  }

  bool operator==(this_type o) { return (*this <= o && o <= *this); }

  void operator|=(this_type o) {
    _w_int_dom |= o._w_int_dom;
    _limit_env = _limit_env | o._limit_env;
    _init_set = _init_set | o._init_set;
  }

  this_type operator|(this_type o) {
    return this_type(_w_int_dom | o._w_int_dom, _limit_env | o._limit_env,
                     _init_set | o._init_set);
  }

  this_type operator&(this_type o) {
    return this_type(_w_int_dom & o._w_int_dom, _limit_env & o._limit_env,
                     _init_set & o._init_set);
  }

  this_type operator||(this_type o) {
    return this_type(_w_int_dom || o._w_int_dom, _limit_env || o._limit_env,
                     _init_set || o._init_set);
  }

  this_type operator&&(this_type o) {
    return this_type(_w_int_dom && o._w_int_dom, _limit_env && o._limit_env,
                     _init_set && o._init_set);
  }

  this_type widening_thresholds(this_type o,
                                const iterators::thresholds<number_t> &ts) {
    return this_type(_w_int_dom.widening_thresholds(o._w_int_dom, ts),
                     _limit_env || o._limit_env, _init_set || o._init_set);
  }

  void set(variable_t x, interval_t i) {
    _w_int_dom.set(x, i);
    if (i.singleton()) {
      // XXX: x's history is reset
      _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());
    } else {
      CRAB_WARN("TODO: set operation with unlimited interval");
      _limit_env -= x;
    }
    _init_set += x;

    CRAB_LOG("wrapped-int2", crab::outs()
                                 << x << ":=" << i << " => " << *this << "\n";);
  }

  void set(variable_t x, wrapped_interval_t i) {
    _w_int_dom.set(x, i);
    // XXX: x's history is reset
    _limit_env.set(x, wrapped_interval_limit_value::convert(i));
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs()
                                 << x << ":=" << i << " => " << *this << "\n";);
  }

  void set(variable_t x, number_t n) {
    _w_int_dom.set(x, n);
    // XXX: x's history is reset
    _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs()
                                 << x << ":=" << n << " => " << *this << "\n";);
  }

  interval_t operator[](variable_t v) const { return _w_int_dom[v]; }

  wrapped_interval_t get_wrapped_interval(variable_t v) const {
    return _w_int_dom.get_wrapped_interval(v);
  }

  wrapped_interval_limit_value get_limit_value(variable_t x) const {
    return _limit_env[x];
  }

  wrapped_interval_domain_t &get_wrapped_interval_domain() {
    return _w_int_dom;
  }

  const wrapped_interval_domain_t &get_wrapped_interval_domain() const {
    return _w_int_dom;
  }

  void operator-=(variable_t v) {
    _w_int_dom -= v;
    _limit_env -= v;
    // XXX: we never remove a variable from _init_set.
    //      This avoids, e.g. to mark as uninitialized a variable that
    //      have been havoc'ed.
    //_init_set -= v;
  }

  // filter variables that do not cross any limit
  void get_not_cross_variables(std::vector<variable_t> &out) const {
    for (typename separate_domain_t::iterator it = _limit_env.begin(),
                                              et = _limit_env.end();
         it != et; ++it) {
      wrapped_interval_limit_value val = it->second;
      if (val.is_not_crossing_limit()) {
        out.push_back(it->first);
      }
    }
  }

  // filter variables that cross only signed limit
  void get_cross_signed_variables(std::vector<variable_t> &out) const {
    for (typename separate_domain_t::iterator it = _limit_env.begin(),
                                              et = _limit_env.end();
         it != et; ++it) {
      wrapped_interval_limit_value val = it->second;
      if (val.is_crossing_signed_limit()) {
        out.push_back(it->first);
      }
    }
  }

  // filter variables that cross only unsigned limit
  void get_cross_unsigned_variables(std::vector<variable_t> &out) const {
    for (typename separate_domain_t::iterator it = _limit_env.begin(),
                                              et = _limit_env.end();
         it != et; ++it) {
      wrapped_interval_limit_value val = it->second;
      if (val.is_crossing_unsigned_limit()) {
        out.push_back(it->first);
      }
    }
  }

  // numerical_domains_api

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, z);
    wrapped_interval_t new_i = get_wrapped_interval(x);
    update_limits(x, old_i, new_i);
    // -- mark x as initialized
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs() << x << ":=" << y << " " << op << " "
                                          << z << " => " << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, k);
    wrapped_interval_t new_i = get_wrapped_interval(x);
    update_limits(x, old_i, new_i);
    // -- mark x as initialized
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs() << x << ":=" << y << " " << op << " "
                                          << k << " => " << *this << "\n";);
  }

  void assign(variable_t x, linear_expression_t e) {
    if (e.is_constant()) {
      // XXX: x's history is reset
      _w_int_dom.assign(x, e);
      _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());
    } else {

      wrapped_interval_t old_i = get_wrapped_interval(x);
      _w_int_dom.assign(x, e);
      wrapped_interval_t new_i = get_wrapped_interval(x);
      update_limits(x, old_i, new_i);
    }
    // -- mark x as initialized
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs()
                                 << x << ":=" << e << " => " << *this << "\n";);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       this_type invariant) {
    _w_int_dom.backward_assign(x, e, invariant._w_int_dom);
    // XXX: ignore _limit_env
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      this_type invariant) {
    _w_int_dom.backward_apply(op, x, y, z, invariant._w_int_dom);
    // XXX: ignore _limit_env
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      this_type invariant) {
    _w_int_dom.backward_apply(op, x, y, z, invariant._w_int_dom);
    // XXX: ignore _limit_env
  }

  void operator+=(linear_constraint_system_t csts) {
    variable_set_t variables = csts.variables();

    std::vector<wrapped_interval_t> old_intervals;
    old_intervals.reserve(variables.size());
    for (auto v : variables) {
      old_intervals.push_back(get_wrapped_interval(v));
    }

    _w_int_dom += csts;

    std::vector<wrapped_interval_t> new_intervals;
    new_intervals.reserve(variables.size());
    for (auto v : variables) {
      new_intervals.push_back(get_wrapped_interval(v));
    }

    update_limits(variables, old_intervals, new_intervals);
    CRAB_LOG("wrapped-int2",
             crab::outs() << "assume(" << csts << ") => " << *this << "\n";);
  }

  // cast operations

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    wrapped_interval_t old_i = get_wrapped_interval(dst);
    _w_int_dom.apply(op, dst, src);
    wrapped_interval_t new_i = get_wrapped_interval(dst);
    update_limits(dst, old_i, new_i);
    // -- mark x as initialized
    _init_set += dst;
    CRAB_LOG("wrapped-int2", crab::outs() << dst << ":=" << op << " " << src
                                          << " => " << *this << "\n";);
  }

  // bitwise operations

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, z);
    wrapped_interval_t new_i = get_wrapped_interval(x);
    update_limits(x, old_i, new_i);
    // -- mark x as initialized
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs() << x << ":=" << y << " " << op << " "
                                          << z << " => " << *this << "\n";);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, k);
    wrapped_interval_t new_i = get_wrapped_interval(x);
    update_limits(x, old_i, new_i);
    // -- mark x as initialized
    _init_set += x;
    CRAB_LOG("wrapped-int2", crab::outs() << x << ":=" << y << " " << op << " "
                                          << k << " => " << *this << "\n";);
  }

  /* Begin unimplemented operations */
  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                this_type invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                this_type invariant) {}
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  this_type invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           this_type invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           this_type invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, this_type invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, this_type invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, this_type invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, this_type invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             this_type invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs) {}
  void pointer_store(variable_t lhs, variable_t rhs) {}
  void pointer_assign(variable_t lhs, variable_t rhs,
                      linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void write(crab_os &o) {
    // o << "(" << _w_int_dom << "," << _limit_env << "," << _init_set << ")";
    o << "(" << _w_int_dom << "," << _limit_env << ")";
  }

  linear_constraint_system_t to_linear_constraint_system() {
    return _w_int_dom.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    return _w_int_dom.to_disjunctive_linear_constraint_system();
  }

  static std::string getDomainName() {
    return "WrappedIntervals+HistoryAbstraction";
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  this_type invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  // checker_domain_traits

  bool entail(const linear_constraint_t &cst) {
    if (is_bottom())
      return true;
    if (cst.is_tautology())
      return true;
    if (cst.is_contradiction())
      return false;

    this_type cst_inv;
    cst_inv += cst;
    // cst cannot be represented by the domain.
    if (cst_inv.is_top())
      return false;

    return get_wrapped_interval_domain() <=
           cst_inv.get_wrapped_interval_domain();
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    _w_int_dom.expand(x, new_x);
    _limit_env.set(new_x, _limit_env[x]);
    if (may_be_initialized(x)) {
      _init_set += new_x;
    }
  }

  void project(const variable_vector_t &vars) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    _w_int_dom.project(vars);

    separate_domain_t projected_env = separate_domain_t::top();
    discrete_domain_t projected_init_set = discrete_domain_t::bottom();
    for (variable_t v : vars) {
      projected_env.set(v, _limit_env[v]);
      if (may_be_initialized(v)) {
        projected_init_set += v;
      }
    }
    std::swap(_limit_env, projected_env);
    std::swap(_init_set, projected_init_set);
  }

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void normalize() {}

  void minimize() {}

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    assert(from.size() == to.size());

    if (is_top() || is_bottom())
      return;

    CRAB_WARN(getDomainName(), "::rename not implemented");
  }
};

template <typename N, typename V, std::size_t M>
struct abstract_domain_traits<wrapped_interval_with_history_domain<N, V, M>> {
  typedef N number_t;
  typedef V varname_t;
};

template <typename N, typename V, std::size_t M>
class checker_domain_traits<wrapped_interval_with_history_domain<N, V, M>> {
public:
  typedef wrapped_interval_with_history_domain<N, V, M> this_type;
  typedef typename this_type::linear_constraint_t linear_constraint_t;
  typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    CRAB_ERROR(
        "TODO: entail operation in wrapper_interval_with_history_domain");
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    CRAB_ERROR(
        "TODO: entail operation in wrapper_interval_with_history_domain");
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    return lhs.entail(rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    // default code
    if (inv.is_bottom() || cst.is_contradiction())
      return false;
    if (inv.is_top() || cst.is_tautology())
      return true;
    this_type cst_inv;
    cst_inv += cst;
    return (!(cst_inv & inv).is_bottom());
  }
};

/*
   Combine the wrapped interval domain + history abstraction with a
   (relational) numerical domain defined over (signed) mathematical
   integers. The result of this combination preserves the nice
   features of the wrapped interval domain (mainly it's sound wrt
   machine arithmetics) while gaining some extra (sound)
   (in)equalities inferred by the numerical domain.
*/
template <typename NumDom, std::size_t max_reduction_cycles = 10>
class wrapped_numerical_domain final
    : public abstract_domain<wrapped_numerical_domain<NumDom>> {
  typedef wrapped_numerical_domain<NumDom> wrapped_numerical_domain_t;
  typedef abstract_domain<wrapped_numerical_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef typename linear_constraint_system_t::variable_set_t variable_set_t;
  typedef typename NumDom::number_t number_t;
  typedef typename NumDom::varname_t varname_t;
  typedef interval<number_t> interval_t;
  typedef crab::pointer_constraint<variable_t> ptr_cst_t;
  typedef typename variable_t::bitwidth_t bitwidth_t;

private:
  typedef wrapped_interval_with_history_domain<number_t, varname_t,
                                               max_reduction_cycles>
      wrapped_interval_domain_t;
  typedef wrapped_interval<number_t> wrapped_interval_t;
  typedef enum { UNKNOWN_SIGNEDNESS, SIGNED, UNSIGNED } signedness_t;
  typedef domain_product2<number_t, varname_t, wrapped_interval_domain_t,
                          NumDom>
      domain_product2_t;

  domain_product2_t _product;

  wrapped_numerical_domain(const domain_product2_t &product)
      : _product(product) {}

  // return true if v may have overflow in the past
  // Use the history abstraction to answer the query.
  inline bool may_have_overflow(variable_t v, signedness_t signedness) {
    wrapped_interval_domain_t &wrapped_intervals = _product.first();
    wrapped_interval_limit_value val = wrapped_intervals.get_limit_value(v);
    return (val.is_top() ||
            (val.is_crossing_signed_limit() &&
             (signedness == SIGNED || signedness == UNKNOWN_SIGNEDNESS)) ||
            (val.is_crossing_unsigned_limit() &&
             (signedness == UNSIGNED || signedness == UNKNOWN_SIGNEDNESS)));
  }

  // return true if interval i fits in [min(b), max(b)]
  inline bool fit(interval_t i, bitwidth_t b, signedness_t signedness) const {
    // TODO: cache min and max
    if (signedness == SIGNED) {
      auto max = wrapint::get_signed_max(b).get_signed_bignum();
      auto min = wrapint::get_signed_min(b).get_signed_bignum();
      CRAB_LOG("wrapped-num-reduction",
               crab::outs() << "\t** checking " << i
                            << " <= " << interval_t(min, max) << "\n";);
      return (i <= interval_t(min, max));
    } else if (signedness == UNSIGNED) {
      auto max = wrapint::get_unsigned_max(b).get_unsigned_bignum();
      auto min = wrapint::get_unsigned_min(b).get_unsigned_bignum();
      CRAB_LOG("wrapped-num-reduction",
               crab::outs() << "\t** checking " << i
                            << " <= " << interval_t(min, max) << "\n";);
      return (i <= interval_t(min, max));
    } else {
      CRAB_LOG("wrapped-num-reduction",
               crab::outs()
                   << "\t** no signedness available. Assume may not fit\n";);
      return false;
    }
  }

  // E.g., if cst is c1*x1 + c2*x2 <= k then we have two residuals:
  //     x1 <= (k - c2*x2) / c1
  //     x2 <= (k - c1*x1) / c2
  // We check whether any of the computations on the rhs of each
  // residual may overflow.
  // We might not need to check for overflow after each intermediate
  // operation, specially if the constraint is interpreted over
  // unsigned integers.
  bool may_overflow_residuals(const linear_constraint_t &cst,
                              number_t coef_pivot, variable_t pivot,
                              signedness_t signedness) {
    bitwidth_t b = pivot.get_bitwidth();
    interval_t residual = cst.constant();
    if (!fit(residual, b, signedness)) {
      // If the constant is to large we bail out
      return true;
    }

    bool res = false;
    for (typename linear_constraint_t::iterator it = cst.begin();
         it != cst.end(); ++it) {
      variable_t v = it->second;
      number_t coef_v = it->first;
      if (!(v == pivot)) {
        CRAB_LOG("wrapped-num-reduction", linear_constraint_t cst_tmp(cst);
                 crab::outs() << "Checking overflow of residual " << cst_tmp
                              << " and " << v << "\n";);
        interval_t tmp = coef_v * _product.first()[v];
        // check if multiplication can overflow
        CRAB_LOG("wrapped-num-reduction",
                 crab::outs() << "\tChecking overflow of " << coef_v << " * "
                              << _product.first()[v] << "=" << tmp << "\n";);
        if (!fit(tmp, b, signedness)) {
          res = true;
          break;
        }

        CRAB_LOG("wrapped-num-reduction",
                 crab::outs() << "\tChecking overflow of " << residual << " - "
                              << tmp << "=";);
        residual = residual - tmp;
        CRAB_LOG("wrapped-num-reduction", crab::outs() << residual << "\n";);
        // check if subtraction can overflow
        if (!fit(residual, b, signedness)) {
          res = true;
          break;
        }

        CRAB_LOG("wrapped-num-reduction",
                 crab::outs() << "\tChecking overflow of " << residual << " / "
                              << coef_pivot << "=";);
        residual = residual / coef_pivot;
        CRAB_LOG("wrapped-num-reduction", crab::outs() << residual << "\n";);
        // check if division can overflow
        if (!fit(residual, b, signedness)) {
          res = true;
          break;
        }
      }
    }
    CRAB_LOG("wrapped-num-reduction",
             if (res) {
               crab::outs() << "Last intermediate computation may overflow!\n";
             } else { crab::outs() << "None of the residuals overflow\n"; });

    return res;
  }

  // return true iff any residual computation of cst may overflow.
  bool may_overflow(const linear_constraint_t &cst, signedness_t signedness) {
    if (cst.is_tautology() || cst.is_contradiction()) {
      return false;
    }

    for (typename linear_constraint_t::iterator it = cst.begin(),
                                                et = cst.end();
         it != et; ++it) {
      number_t c = it->first;
      variable_t pivot = it->second;
      if (may_overflow_residuals(cst, c, pivot, signedness)) {
        return true;
      }
    }
    return false;
  }

  // return true iff all variables in vars have the same bitwidth
  static bool has_same_bitwidth(const variable_set_t &vars) {
    bool first = true;
    typename variable_t::bitwidth_t b;
    bool same_bitwidth = true;
    for (auto const &v : vars) {
      if (first) {
        b = v.get_bitwidth();
        first = false;
      } else {
        if (b != v.get_bitwidth()) {
          same_bitwidth = false;
          break;
        }
      }
    }
    return same_bitwidth;
  }

  // If the wrapped domain + history abstraction cannot tell that v
  // does not overflow then v is forgotten from the numerical domain.
  //
  // This operation can be called before any operation. However, since
  // we are only interested in preserving reachability and for
  // accuracy, we only apply it on operations with explicit signedness
  // interpretation (branches, div, and rem).
  void rectify(const variable_t &v, signedness_t signedness) {
    CRAB_LOG(
        "wrapped-num", crab::outs() << "Correcting " << v;
        if (signedness == SIGNED) {
          crab::outs() << " signed: \n";
        } else if (signedness == UNSIGNED) {
          crab::outs() << " unsigned: \n";
        } else { crab::outs() << " unknown unsignedness: \n"; });
    if (may_have_overflow(v, signedness)) {
      _product.second() -= v;
      CRAB_LOG("wrapped-num", crab::outs() << "\t" << v << " may overflow!\n";);
    } else {
      CRAB_LOG("wrapped-num", crab::outs()
                                  << "\t" << v << " cannot overflow\n";);
    }
  }

  // Reduction from a potentially "unsound" numerical domain to the
  // wrapped interval domain. It consists of propagating linear
  // (in)equalities. We quote unsound because the domain is unsound
  // wrt to machine arithmetic although the domain is sound wrt
  // mathematical integers.
  void strengthen(const variable_set_t rel_vars) {
    if (is_bottom() || is_top()) {
      return;
    }
    // -- extract all constraints involving any variable in rel_vars
    linear_constraint_system_t csts;
    for (auto const &v : rel_vars) {
      reduced_domain_traits<NumDom>::extract(_product.second(), v, csts,
                                             /* only equalities=*/false);
    }

    // IMPORTANT: we mark as SIGNED here because all the "unsound"
    // numerical domains in Crab interpret numbers as signed.
    const signedness_t signedness = SIGNED;
    // -- filter out constraints that shouldn't be propagated to the
    // -- wrapped interval domain
    for (auto const &c : csts) {
      variable_set_t vars = c.variables();
      if (!has_same_bitwidth(vars)) {
        continue;
      }
      if (may_overflow(c, signedness)) {
        continue;
      }
      CRAB_LOG("wrapped-num", linear_constraint_t tmp(c);
               crab::outs() << "** reduction propagating " << tmp
                            << " to wrapped intervals\n";);
      CRAB_LOG("wrapped-num-reduction", crab::outs()
                                            << "\tBEFORE" << *this << "\n";);
      _product.first() += c;
      CRAB_LOG("wrapped-num-reduction", crab::outs()
                                            << "\tAFTER" << *this << "\n";);
    }

#if 0
    { // Propagate wrapped intervals to the relational domain. 
      linear_constraint_system_t csts;
      for (auto const& v: rel_vars) {
	interval_t i = _product.first()[v];
	boost::optional<number_t> lb = i.lb().number();
	boost::optional<number_t> ub = i.ub().number();
	if (lb) csts += linear_constraint_t(v >= *lb);
	if (ub) csts += linear_constraint_t(v <= *ub);
      }
      _product.second() += csts;
    }
#endif
  }

  void strengthen(const variable_t &x) {
    variable_set_t vars;
    vars += x;
    strengthen(vars);
  }

public:
  void set_to_top() {
    wrapped_numerical_domain_t abs(domain_product2_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    wrapped_numerical_domain_t abs(domain_product2_t::bottom());
    std::swap(*this, abs);
  }

  wrapped_numerical_domain() : _product() {}

  wrapped_numerical_domain(const wrapped_numerical_domain_t &other)
      : _product(other._product) {}

  wrapped_numerical_domain_t &
  operator=(const wrapped_numerical_domain_t &other) {
    if (this != &other) {
      _product = other._product;
    }
    return *this;
  }

  bool is_bottom() { return _product.is_bottom(); }

  bool is_top() { return _product.is_top(); }

  wrapped_numerical_domain_t &first() { return _product.first(); }

  NumDom &second() { return _product.second(); }

  bool operator<=(wrapped_numerical_domain_t other) {
    return _product <= other._product;
  }

  bool operator==(wrapped_numerical_domain_t other) {
    return _product == other._product;
  }

  void operator|=(wrapped_numerical_domain_t other) {
    CRAB_LOG("wrapped-num",
             crab::outs() << _product << " U " << other._product << " = ");
    _product |= other._product;
    CRAB_LOG("wrapped-num", crab::outs() << _product << "\n";);
  }

  wrapped_numerical_domain_t operator|(wrapped_numerical_domain_t other) {
    CRAB_LOG("wrapped-num",
             crab::outs() << _product << " U " << other._product << " = ");
    wrapped_numerical_domain_t res(_product | other._product);
    CRAB_LOG("wrapped-num", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_numerical_domain_t operator&(wrapped_numerical_domain_t other) {
    CRAB_LOG("wrapped-num", crab::outs() << "MEET " << _product << " "
                                         << other._product << " = ");
    wrapped_numerical_domain_t res(_product & other._product);
    CRAB_LOG("wrapped-num", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_numerical_domain_t operator||(wrapped_numerical_domain_t other) {
    CRAB_LOG("wrapped-num", crab::outs() << "WIDENING " << _product << "and "
                                         << other._product << " = ");
    wrapped_numerical_domain_t res(_product || other._product);
    CRAB_LOG("wrapped-num", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_numerical_domain_t
  widening_thresholds(wrapped_numerical_domain_t other,
                      const iterators::thresholds<number_t> &ts) {
    CRAB_LOG("wrapped-num", crab::outs() << "WIDENING " << _product << "and "
                                         << other._product << " = ");
    wrapped_numerical_domain_t res(
        _product.widening_thresholds(other._product, ts));
    CRAB_LOG("wrapped-num", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_numerical_domain_t operator&&(wrapped_numerical_domain_t other) {
    CRAB_LOG("wrapped-num", crab::outs() << "NARROWING " << _product << "and "
                                         << other._product << " = ");
    wrapped_numerical_domain_t res(_product && other._product);
    CRAB_LOG("wrapped-num", crab::outs() << res << "\n";);
    return res;
  }

  // numerical_domains_api

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    if (op == OP_SDIV || op == OP_SREM) {
      // signed division/rem
      rectify(y, SIGNED);
      rectify(z, SIGNED);
    } else if (op == OP_UDIV || op == OP_UREM) {
      // unsigned division/rem
      rectify(y, UNSIGNED);
      rectify(z, UNSIGNED);
    }
    if (op == OP_UDIV || op == OP_UREM) {
      // if unsigned division then we only apply it on wrapped intervals
      _product.first().apply(op, x, y, z);
    } else {
      _product.apply(op, x, y, z);
    }
    CRAB_LOG("wrapped-num", crab::outs() << x << ":=" << y << " " << op << " "
                                         << z << "=" << _product << "\n";);
    strengthen(x);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    if (op == OP_SDIV || op == OP_SREM) {
      // signed division/rem
      rectify(y, SIGNED);
    } else if (op == OP_UDIV || op == OP_UREM) {
      // unsigned division/rem
      rectify(y, UNSIGNED);
    }

    if (op == OP_UDIV || op == OP_UREM) {
      // if unsigned division then we only apply it on wrapped intervals
      _product.first().apply(op, x, y, k);
    } else {

      _product.apply(op, x, y, k);
    }
    CRAB_LOG("wrapped-num", crab::outs() << x << ":=" << y << " " << op << " "
                                         << k << "=" << _product << "\n";);
    strengthen(x);
  }

  void assign(variable_t x, linear_expression_t e) {
    _product.assign(x, e);
    CRAB_LOG("wrapped-num", crab::outs()
                                << x << ":=" << e << "=" << _product << "\n";);

    if (!e.is_constant()) {
      strengthen(x);
    }
  }

  void operator+=(linear_constraint_system_t csts) {
    // Add first the constraint in the wrapped interval domain
    _product.first() += csts;

    linear_constraint_system_t non_overflow_csts;
    CRAB_LOG(
        "wrapped-num", crab::outs() << "BEGIN add constraints {"; for (auto c
                                                                       : csts) {
          crab::outs() << c << ";";
        } crab::outs() << "}\n";);

    for (auto c : csts) {
      // rectify of the "unsound" numerical domain
      signedness_t signedness = UNKNOWN_SIGNEDNESS;
      if (c.is_inequality() || c.is_strict_inequality()) {
        signedness = c.is_signed() ? SIGNED : UNSIGNED;
      }
      for (auto const &v : c.variables()) {
        rectify(v, signedness);
      }
      // add constraint in the "unsound" numerical domain only if
      // c cannot overflow.
      if (!may_overflow(c, signedness)) {
        non_overflow_csts += c;
        CRAB_LOG("wrapped-num", crab::outs()
                                    << "** added constraint: " << c << "\n");
      }
    }
    _product.second() += non_overflow_csts;

    CRAB_LOG("wrapped-num", crab::outs()
                                << "END add constraints: " << _product << "\n");
  }

  void set(variable_t x, interval_t intv) {
    // domain_product2 does not define set method
    _product.first().set(x, intv);
    _product.second().set(x, intv);
  }

  interval_t operator[](variable_t v) {
    // domain_product2 does not define [] method
    return (_product.first()[v] & _product.second()[v]);
  }

  void operator-=(variable_t v) { _product -= v; }

  // backward arithmetic operations
  void backward_assign(variable_t x, linear_expression_t e,
                       wrapped_numerical_domain_t invariant) override {

    CRAB_WARN("backward assign not implemented");
    this->operator-=(x);

    //_product.backward_assign(x,e,invariant._product);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      wrapped_numerical_domain_t invariant) override {
    CRAB_WARN("backward apply not implemented");
    this->operator-=(x);

    // _product.backward_apply(op,x,y,z,invariant._product);
    // if (op == OP_SDIV) {
    //   _product.second() -= x;
    // }
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      wrapped_numerical_domain_t invariant) override {
    CRAB_WARN("backward apply not implemented");
    this->operator-=(x);

    // _product.backward_apply(op,x,y,z,invariant._product);
    // if (op == OP_SDIV) {
    //   _product.second() -= x;
    // }
  }

  // boolean_operators

  void assign_bool_cst(variable_t x, linear_constraint_t cst) override {
    // Add first the constraint in the wrapped interval domain
    _product.first().assign_bool_cst(x, cst);

    signedness_t signedness = UNKNOWN_SIGNEDNESS;
    if (cst.is_inequality() || cst.is_strict_inequality()) {
      signedness = cst.is_signed() ? SIGNED : UNSIGNED;
    }
    for (auto const &v : cst.variables()) {
      rectify(v, signedness);
    }
    if (!may_overflow(cst, signedness)) {
      _product.second().assign_bool_cst(x, cst);
    }
  }

  void assign_bool_var(variable_t x, variable_t y, bool is_not_y) override {
    _product.assign_bool_var(x, y, is_not_y);
  }

  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) override {
    _product.apply_binary_bool(op, x, y, z);
  }

  void assume_bool(variable_t x, bool is_negated) override {
    _product.assume_bool(x, is_negated);
  }

  // backward boolean operators
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                wrapped_numerical_domain_t inv) {
    CRAB_WARN("backward assign bool constraint not implemented");
    this->operator-=(lhs);
  }

  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                wrapped_numerical_domain_t inv) {
    CRAB_WARN("backward assign bool variable not implemented");
    this->operator-=(lhs);
  }

  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  wrapped_numerical_domain_t inv) {
    CRAB_WARN("backward apply binary bool not implemented");
    this->operator-=(x);
  }

  // cast operations

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    // FIXME/TODO: we might need to throw away relationships between dst and src
    _product.apply(op, dst, src);
  }

  // bitwise operations

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    _product.apply(op, x, y, z);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    _product.apply(op, x, y, k);
  }

  // array operations

  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {
    _product.array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {
    _product.array_load(lhs, a, elem_size, i);
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    _product.array_store(a, elem_size, i, val, is_strong_update);
  }

  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t v) override {
    _product.array_store_range(a, elem_size, i, j, v);
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    _product.array_assign(lhs, rhs);
  }

  // backward array operations
  virtual void
  backward_array_init(variable_t a, linear_expression_t elem_size,
                      linear_expression_t lb_idx, linear_expression_t ub_idx,
                      linear_expression_t val,
                      wrapped_numerical_domain_t invariant) override {}
  virtual void
  backward_array_load(variable_t lhs, variable_t a,
                      linear_expression_t elem_size, linear_expression_t i,
                      wrapped_numerical_domain_t invariant) override {}
  virtual void
  backward_array_store(variable_t a, linear_expression_t elem_size,
                       linear_expression_t i, linear_expression_t val,
                       bool is_strong_update,
                       wrapped_numerical_domain_t invariant) override {}
  virtual void
  backward_array_store_range(variable_t a, linear_expression_t elem_size,
                             linear_expression_t i, linear_expression_t j,
                             linear_expression_t v,
                             wrapped_numerical_domain_t invariant) override {}
  virtual void
  backward_array_assign(variable_t lhs, variable_t rhs,
                        wrapped_numerical_domain_t invariant) override {}
  // pointer operations
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    _product.pointer_load(lhs, rhs);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    _product.pointer_store(lhs, rhs);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    _product.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    _product.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    _product.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    _product.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(ptr_cst_t cst) override {
    _product.pointer_assume(cst);
  }

  virtual void pointer_assert(ptr_cst_t cst) override {
    _product.pointer_assert(cst);
  }

  void write(crab_os &o) { _product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    linear_constraint_system_t res;
    res += _product.first().to_linear_constraint_system();
    res += _product.second().to_linear_constraint_system();
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    disjunctive_linear_constraint_system_t res;
    res += _product.first().to_disjunctive_linear_constraint_system();
    res += _product.second().to_disjunctive_linear_constraint_system();
    return res;
  }

  static std::string getDomainName() {
    return domain_product2_t::getDomainName();
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    _product.rename(from, to);
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    _product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  wrapped_numerical_domain_t invariant) override {
    _product.backward_intrinsic(name, inputs, outputs, invariant._product);    
  }
  /* end intrinsics operations */
  
  void normalize() { _product.normalize(); }

  void minimize() { _product.minimize(); }

  void expand(variable_t x, variable_t new_x) { _product.expand(x, new_x); }

  void forget(const variable_vector_t &vars) { _product.forget(vars); }

  void project(const variable_vector_t &vars) { _product.project(vars); }
}; // class wrapped_numerical_domain

template <typename AbsDom>
struct abstract_domain_traits<wrapped_numerical_domain<AbsDom>> {
  typedef typename AbsDom::number_t number_t;
  typedef typename AbsDom::varname_t varname_t;
};

template <typename AbsDom>
class checker_domain_traits<wrapped_numerical_domain<AbsDom>> {
public:
  typedef wrapped_numerical_domain<AbsDom> this_type;
  typedef typename this_type::linear_constraint_t linear_constraint_t;
  typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    return checker_domain_traits<AbsDom>::entail(lhs.second(), rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    return checker_domain_traits<AbsDom>::entail(lhs, rhs.second());
  }

  static bool entail(this_type &inv, const linear_constraint_t &cst) {
    return checker_domain_traits<AbsDom>::entail(inv.second(), cst);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    return checker_domain_traits<AbsDom>::intersect(inv.second(), cst);
  }
};

} // namespace domains
} // namespace crab
