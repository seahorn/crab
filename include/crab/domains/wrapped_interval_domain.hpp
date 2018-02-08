#pragma once

/**
 ** Machine arithmetic interval domain based on the paper
 ** "Signedness-Agnostic Program Analysis: Precise Integer Bounds for
 ** Low-Level Code" by J.A.Navas, P.Schachte, H.Sondergaard, and
 ** P.J.Stuckey published in APLAS'12.
 **/

#include <crab/common/wrapint.hpp>
#include <crab/common/types.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/operators_api.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/domain_traits.hpp>
#include <boost/optional.hpp>

#define PRINT_WRAPINT_AS_SIGNED

namespace crab{
namespace domains {

template<typename Number>  
class wrapped_interval {
  
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
    if (is_bottom()) return;

    bitwidth_t b = get_bitwidth(__LINE__);
    if (is_top()) {
      intervals.push_back(wrapped_interval_t(wrapint::get_unsigned_min(b),
					     wrapint::get_signed_max(b)));
      intervals.push_back(wrapped_interval_t(wrapint::get_signed_min(b),
					     wrapint::get_unsigned_max(b)));
    } else {
      if (signed_limit(b) <= *this) {
	intervals.push_back(wrapped_interval_t(_start, wrapint::get_signed_max(b)));
	intervals.push_back(wrapped_interval_t(wrapint::get_signed_min(b), _stop));
      } else {
	intervals.push_back(*this);
      }
    }
  }

  // ssplit in the APLAS'12 paper
  void unsigned_split(std::vector<wrapped_interval_t> &intervals) const {
    if (is_bottom()) return;

    bitwidth_t b = get_bitwidth(__LINE__);
    if (is_top()) {
      intervals.push_back(wrapped_interval_t(wrapint::get_signed_min(b),
					     wrapint::get_unsigned_max(b)));
      intervals.push_back(wrapped_interval_t(wrapint::get_unsigned_min(b),
					     wrapint::get_signed_max(b)));
    } else {
      if (unsigned_limit(b) <= *this) {
	intervals.push_back(wrapped_interval_t(_start, wrapint::get_unsigned_max(b)));
	intervals.push_back(wrapped_interval_t(wrapint::get_unsigned_min(b), _stop));
      } else {
	intervals.push_back(*this);
      }
    }
  }

  // cut in the APLAS'12 paper
  void signed_and_unsigned_split(std::vector<wrapped_interval_t>& out) const {
    std::vector<wrapped_interval_t> ssplit;
    signed_split(ssplit);
    for(unsigned i=0, e=ssplit.size(); i<e; ++i) {
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
    if (msb_start == msb_stop && msb_stop == msb_x_start && msb_x_start == msb_x_stop) {
      // the two intervals are in the same hemisphere
      if (!msb_start) {
	return unsigned_mul(x);
      } else {
	if ((_start.get_unsigned_bignum() * x._start.get_unsigned_bignum()) -
	    (_stop.get_unsigned_bignum() * x._stop.get_unsigned_bignum())
	    < wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
	  res = wrapped_interval_t(_stop * x._stop, _start * x._start);
	}
	goto EXIT_SIGNED_MUL;
      }
    }

    // each interval cannot cross the limit: one interval in a
    // different hemisphere.
    if (!(msb_start != msb_stop || msb_x_start != msb_x_stop)) {
      if (msb_start && !msb_x_start) {
	if ((_stop.get_unsigned_bignum() * x._start.get_unsigned_bignum()) -
	    (_start.get_unsigned_bignum() * x._stop.get_unsigned_bignum())
	    < wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
	  res = wrapped_interval_t(_start * x._stop, _stop * x._start);
	}
      } else if (!msb_start && msb_x_start) {
	if ((_start.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) -
	    (_stop.get_unsigned_bignum() * x._start.get_unsigned_bignum())
	    < wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
	  res = wrapped_interval_t(_stop * x._start, _start * x._stop);
	}
      } 
    }
    
  EXIT_SIGNED_MUL:
    CRAB_LOG("wrapped-int-mul",
	     crab::outs() << "Signed " << *this << " * " << x << "=" << res << "\n";);    
    return res;
  }

  wrapped_interval_t unsigned_mul(wrapped_interval_t x) const {
    assert(!is_bottom() && !x.is_bottom());
 
    bitwidth_t b = get_bitwidth(__LINE__);
    wrapped_interval_t res = wrapped_interval_t::top();
    if ((_stop.get_unsigned_bignum() * x._stop.get_unsigned_bignum()) -
	(_start.get_unsigned_bignum() * x._start.get_unsigned_bignum())
	< wrapint::get_unsigned_max(b).get_unsigned_bignum()) {
      res = wrapped_interval_t(_start * x._start, _stop * x._stop);
    }
    CRAB_LOG("wrapped-int-mul",
	     crab::outs() << "Unsigned " << *this << " * " << x << "=" << res << "\n";);        
    return res;
  }

  // if out is empty then the intersection is empty
  void exact_meet(wrapped_interval_t x, std::vector<wrapped_interval_t>& out) const {
    if (is_bottom() || x.is_bottom()){
      // bottom
    } else if (*this == x || is_top()) {
      out.push_back(x);
    } else if (x.is_top()) {
      out.push_back(*this);
    } else if (x[_start] && x[_stop] && operator[](x._start) && operator[](x._stop)) {
      out.push_back(wrapped_interval_t(_start, x._stop));
      out.push_back(wrapped_interval_t(x._start, _stop));
    } else if (x[_start] && x[_stop]) {
      out.push_back(*this);
    } else if (operator[](x._start) && operator[](x._stop)) {
      out.push_back(x);
    } else if (x[_start] && operator[](x._stop) && !x[_stop] && operator[](x._start)) {
      out.push_back(wrapped_interval_t(_start, x._stop));
    } else if (x[_stop] && operator[](x._start) && !x[_start] && operator[](x._stop)) {
      out.push_back(wrapped_interval_t(x._start, _stop));
    } else {
      // bottom
    }
  }

  // Perform the reduced product of signed and unsigned multiplication.
  // It uses exact meet rather than abstract meet.
  void reduced_signed_unsigned_mul(wrapped_interval_t x,
				   std::vector<wrapped_interval_t>& out) const {
    if (is_bottom() || x.is_bottom()) {
      return;
    }
    
    wrapped_interval_t s = signed_mul(x);
    wrapped_interval_t u = unsigned_mul(x);
    s.exact_meet(u, out);
    CRAB_LOG("wrapped-int-mul",
	     crab::outs() << "Exact signed x unsigned " << s << " * " << u << "=\n";
	     for(unsigned i=0;i<out.size();++i) {
	       crab::outs () << "\t" << out[i] << "\n";
	     });
  }

  wrapped_interval_t unsigned_div(wrapped_interval_t x) const {
    CRAB_LOG("wrapped-int-div", crab::outs () << *this << " /_u " << x  << "=";);
    assert(!x[wrapint(0, x.get_bitwidth(__LINE__))]);
    assert(!is_bottom() && !x.is_bottom());
    wrapped_interval_t res = wrapped_interval_t(_start.udiv(x._stop), _stop.udiv(x._start));
    CRAB_LOG("wrapped-int-div", crab::outs () << res << "\n";);
    return res;
  }

  wrapped_interval_t signed_div(wrapped_interval_t x) const {
    CRAB_LOG("wrapped-int-div", crab::outs () << *this << " /_s " << x  << "=";);    
    assert(!x[wrapint(0, x.get_bitwidth(__LINE__))]);
    assert(!is_bottom() && !x.is_bottom());

    bool msb_start = _start.msb();
    bool msb_x_start = x._start.msb();

    wrapped_interval_t res;
    if (msb_start == msb_x_start) {
      if (msb_start) { //both negative
	res = wrapped_interval_t(_stop.sdiv(x._start), _start.sdiv(x._stop));	
      } else {  //both positive
	res = wrapped_interval_t(_start.sdiv(x._stop), _stop.sdiv(x._start));
      }
    } else {
      if (msb_start) { 
	assert(!msb_x_start);
	res = wrapped_interval_t(_start.sdiv(x._start), _stop.sdiv(x._stop));
      } else {
	assert(msb_x_start);
	res = wrapped_interval_t(_stop.sdiv(x._stop), _start.sdiv(x._start));
      }
    }
    CRAB_LOG("wrapped-int-div", crab::outs () << res << "\n";);    
    return res;
  }
  
  // FIXME: this is sound only if wrapped interval defined over
  //        z_number.
  void trim_zero(std::vector<wrapped_interval_t>& out) const {
    wrapint zero(0, get_bitwidth(__LINE__));
    if (!is_bottom() && (!(*this == zero))) {
      if (start() == zero) {
	out.push_back(wrapped_interval_t(wrapint(1, get_bitwidth(__LINE__)), stop()));
      } else if (stop() == zero) {
	out.push_back(wrapped_interval_t(start(),wrapint(-1, get_bitwidth(__LINE__))));
      } else if (operator[](zero)) {
	out.push_back(wrapped_interval_t(start(), wrapint(-1, get_bitwidth(__LINE__))));
	out.push_back(wrapped_interval_t(wrapint(1, get_bitwidth(__LINE__)), stop()));
      } else {
	out.push_back(*this);
      }
    }
  }

  wrapped_interval_t Shl(uint64_t k) const {
    if (is_bottom()) return *this;

    // XXX: we need the check is_top before calling get_bitwidth();
    if (is_top()) return *this;
    
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
    if (is_bottom()) return *this;

    // XXX: we need the check is_top before cross_signed_limit calls
    // get_bitwidth();
    if (is_top()) return *this;
    
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
    if (is_bottom()) return *this;

    // XXX: we need the check is_top before cross_signed_limit calls
    // get_bitwidth();
    if (is_top()) return *this;
    
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
  
  
public:

  typedef wrapint::bitwidth_t bitwidth_t;

  wrapped_interval(wrapint n)
    : _start(n), _stop(n), _is_bottom(false) { }
  
  wrapped_interval(wrapint start, wrapint stop)
    : _start(start), _stop(stop), _is_bottom(false) {
    if (start.get_bitwidth() != stop.get_bitwidth()) {
      CRAB_ERROR("inconsistent bitwidths in wrapped interval");
    }
  }

  wrapped_interval(Number n, bitwidth_t w)
    : _start(wrapint(n,w)), _stop(wrapint(n, w)), _is_bottom(false) {}
  
  // To represent top, the particular bitwidth here is irrelevant. We
  // just make sure that _stop - _start == get_max()
  wrapped_interval()
    : _start(wrapint(0,3)), _stop(7,3), _is_bottom(false) { }
  
  static wrapped_interval_t top () {
    return wrapped_interval_t(wrapint(0,3), wrapint(7,3), false); 
  }

  static wrapped_interval_t bottom () {
    // the wrapint is irrelevant.
    wrapint i(0,1);
    return wrapped_interval_t(i,i, true);
  }

  // return interval [0111...1, 1000....0]
  // In the APLAS'12 paper "signed limit" corresponds to "north pole".
  static wrapped_interval_t signed_limit(bitwidth_t b) {
    return wrapped_interval_t (wrapint::get_signed_max(b), wrapint::get_signed_min(b));
  }

  // return interval [1111...1, 0000....0]
  // In the APLAS'12 paper "unsigned limit" corresponds to "south pole".  
  static wrapped_interval_t unsigned_limit(bitwidth_t b) {
    return wrapped_interval_t(wrapint::get_unsigned_max(b), wrapint::get_unsigned_min(b));
  }
  
  bool cross_signed_limit() const {
    return (signed_limit(get_bitwidth(__LINE__)) <= *this);
  }

  bool cross_unsigned_limit() const {
    return (unsigned_limit(get_bitwidth(__LINE__)) <= *this);
  }
  
  bitwidth_t get_bitwidth(int line) const {
    if (is_bottom()) {
      CRAB_ERROR("get_bitwidth() cannot be called from a bottom element at line ", line);
    } else if (is_top()) {
      CRAB_ERROR("get_bitwidth() cannot be called from a top element at line ", line);
    } else {
      assert (_start.get_bitwidth() == _stop.get_bitwidth());
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
  
  bool is_bottom () const {
    return _is_bottom;
  }
  
  bool is_top () const {
    wrapint maxspan = wrapint::get_unsigned_max(_start.get_bitwidth());
    return (!_is_bottom && (_stop - _start == maxspan));
  }

  // Important: we make the choice here that we interpret wrapint as
  // signed mathematical integers.
  interval<Number> to_interval() const {
    typedef interval<Number> interval_t;
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top() || (cross_signed_limit ())) {
      return interval_t::top();
    } else {
      return interval_t(_start.get_signed_bignum(), _stop.get_signed_bignum());
    }
  }
  
  wrapped_interval_t lower_half_line(bool is_signed) const {
    if (is_top() || is_bottom()) return *this;
    wrapint smin = wrapint::get_signed_min(get_bitwidth(__LINE__));
    if (!is_signed) {
      smin = wrapint::get_unsigned_min(get_bitwidth(__LINE__));
    }
    wrapped_interval_t res = wrapped_interval_t(smin, _stop);
    return res;
  }
  
  wrapped_interval_t upper_half_line(bool is_signed) const {
    if (is_top() || is_bottom()) return *this;
    wrapint smax = wrapint::get_signed_max(get_bitwidth(__LINE__));
    if (!is_signed) {
      smax = wrapint::get_unsigned_max(get_bitwidth(__LINE__));
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
    } else  {
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
      return x[_start] && x[_stop] && (!(operator[](x._start)) || !(operator[](x._stop)));
    }
  }

  bool operator==(wrapped_interval_t x) const {
    return (*this <= x && x <= *this);
  }

  bool operator!=(wrapped_interval_t x) const {
    return !(this->operator==(x));
  }

  wrapped_interval_t operator|(wrapped_interval_t x) const {
    if (*this <= x) {
      return x;
    } else if (x <= *this) { 
      return *this;
    } else {
      if (x[_start] && x[_stop] && operator[](x._start) && operator[](x._stop)) {
	return wrapped_interval_t::top();
      } else if (x[_stop] && operator[](x._start)) {
	return wrapped_interval_t(_start, x._stop); 
      } else if (operator[](x._stop) && x[_start]) {
	return wrapped_interval_t(x._start, _stop);
      } else {
	wrapint span_a = x._start - _stop;
	wrapint span_b = _start - x._stop;
	if (span_a < span_b  || (span_a == span_b &&  _start <= x._start)) {
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
	  if (span_a < span_b  || (span_a == span_b &&  _start <= x._start)) {
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
    assert (get_bitwidth(__LINE__) == x.get_bitwidth(__LINE__));
    bitwidth_t w = x.get_bitwidth(__LINE__);
    if ((_stop - _start) >= wrapint::get_signed_min(w))
      return wrapped_interval_t::top();
    
    wrapped_interval_t join = *this | x;
    if (join == wrapped_interval_t(_start, x._stop)) {
      // increase by twice the size of the old interval
      wrapped_interval_t doubled_old(_start, (_stop * wrapint(2,w)) - _start + wrapint(1,w));
      return join | doubled_old;
    } else if (join == wrapped_interval_t(x._start, _stop)) {
      // decrease by twice the size of the old interval      
      wrapped_interval_t doubled_old((_start * wrapint(2,w)) - _stop - wrapint(1,w), _stop);
      return join | doubled_old;
    } else if (x[_start] && x[_stop]) {
      wrapped_interval_t y(x._start, x._start +  (_stop * wrapint(2,w)) - (_start * wrapint(2,w)) + wrapint(1,w));
      return x | y;
    } else {
      return wrapped_interval_t::top();
    }
  }
  
  template<typename Thresholds>
  wrapped_interval_t widening_thresholds (wrapped_interval_t x, const Thresholds &ts) {
    // TODO: for now we call the widening operator without thresholds
    return (this->operator||(x));
  }
  
  wrapped_interval_t operator&&(wrapped_interval_t x) const {
    // TODO: for now we call the meet operator.
    return (this->operator&(x));
  }

  wrapped_interval_t operator+(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } else if (is_top () || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      // -- check if the addition will overflow
      wrapint x_sz = x._stop - x._start;
      wrapint sz = _stop - _start;
      wrapint one(1, x_sz.get_bitwidth());
      if (x_sz  + sz + one <= x_sz) {
	return wrapped_interval_t::top();
      }
      
      return wrapped_interval_t(_start + x._start, _stop + x._stop);
    }
  }
  
  wrapped_interval_t& operator+=(wrapped_interval_t x) {
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
      if (x_sz  + sz + one <= x_sz) {
	return wrapped_interval_t::top();
      }
      return wrapped_interval_t(_start - x._stop, _stop - x._start);
    }
  }
  
  wrapped_interval_t& operator-=(wrapped_interval_t x) {
    return this->operator=(this->operator-(x));
  }
  
  wrapped_interval_t operator*(wrapped_interval_t x) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> cuts, x_cuts;
      signed_and_unsigned_split(cuts);
      x.signed_and_unsigned_split(x_cuts);
      assert(!cuts.empty());
      assert(!x_cuts.empty());

      CRAB_LOG("wrapped-int-mul",
	       crab::outs () << "cuts for " << *this << "\n";
	       for(unsigned i=0, ie=cuts.size(); i < ie; ++i) {
		 crab::outs() << "\t" << cuts[i] << "\n";
	       }
	       crab::outs () << "cuts for " << x << "\n";
	       for(unsigned i=0, ie=x_cuts.size(); i < ie; ++i) {
		 crab::outs() << "\t" << x_cuts[i] << "\n";
	       });
		 
      wrapped_interval_t res = wrapped_interval_t::bottom();      
      for(unsigned i=0, ie=cuts.size(); i < ie; ++i) {
	for(unsigned j=0, je=x_cuts.size(); j < je; ++j) {
	  std::vector<wrapped_interval_t> exact_reduct;
	  cuts[i].reduced_signed_unsigned_mul(x_cuts[j], exact_reduct);
	  for (unsigned k=0; k < exact_reduct.size(); ++k) {
	    res = res | exact_reduct[k];
	  }
	}
      }

      CRAB_LOG("wrapped-int-mul",
	       crab::outs () << *this << " * " << x << " = " << res << "\n");
      return res;
    }
  }
  
  wrapped_interval_t& operator*=(wrapped_interval_t x) {
    return this->operator=(this->operator*(x));
  }

  // signed division
  wrapped_interval_t operator/(wrapped_interval_t x) const {
    return SDiv(x);
  }
  
  wrapped_interval_t& operator/=(wrapped_interval_t x) {
    return this->operator=(this->operator/(x));
  }   
    
  /** division and remainder operations **/
  wrapped_interval_t SDiv(wrapped_interval_t x ) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> cuts, x_cuts;
      signed_and_unsigned_split(cuts);
      x.signed_and_unsigned_split(x_cuts);
      assert(!cuts.empty());
      assert(!x_cuts.empty());
      wrapped_interval_t res = wrapped_interval_t::bottom();      
      for(unsigned i=0, ie=cuts.size(); i < ie; ++i) {
	for(unsigned j=0, je=x_cuts.size(); j < je; ++j) {
	  std::vector<wrapped_interval_t> trimmed_divisors;
	  x_cuts[j].trim_zero(trimmed_divisors);
	  for(unsigned k=0, ke = trimmed_divisors.size(); k < ke; ++k) {
	    wrapped_interval_t d = trimmed_divisors[k];
	    res = res | cuts[i].signed_div(d);	    
	  }
	}
      }
      return res;
    }
  }
  
  wrapped_interval_t UDiv(wrapped_interval_t x ) const {
    if (is_bottom() || x.is_bottom()) {
      return wrapped_interval_t::bottom();
    } if (is_top() || x.is_top()) {
      return wrapped_interval_t::top();
    } else {
      std::vector<wrapped_interval_t> ssplits, x_ssplits;
      signed_split(ssplits);
      x.signed_split(x_ssplits);
      assert(!ssplits.empty());
      assert(!x_ssplits.empty());
      wrapped_interval_t res = wrapped_interval_t::bottom();      
      for(unsigned i=0, ie=ssplits.size(); i < ie; ++i) {
	for(unsigned j=0, je=x_ssplits.size(); j < je; ++j) {
	  std::vector<wrapped_interval_t> trimmed_divisors;
	  x_ssplits[j].trim_zero(trimmed_divisors);
	  for(unsigned k=0, ke = trimmed_divisors.size(); k < ke; ++k) {
	    wrapped_interval_t d = trimmed_divisors[k];
	    res = res | ssplits[i].unsigned_div(d);	    
	  }
	}
      }
      return res;
    }
  }
  
  wrapped_interval_t SRem(wrapped_interval_t x)  const
  { return default_implementation(x); }    
  
  wrapped_interval_t URem(wrapped_interval_t x)  const
  { return default_implementation(x); }    

  /** cast operations **/
  wrapped_interval_t ZExt(unsigned bits_to_add)  const {
    std::vector<wrapped_interval_t> intervals;
    unsigned_split(intervals);
    
    wrapped_interval_t res = wrapped_interval_t::bottom();
    for (typename std::vector<wrapped_interval_t>::iterator it = intervals.begin(),
	   et = intervals.end();
	 it!=et; ++it) {
      // this should not happen
      if ((*it).is_bottom() || (*it).is_top()) continue;
      
      wrapint a = (*it).start();
      wrapint b = (*it).stop();
      res = res | wrapped_interval_t(a.zext(bits_to_add), b.zext(bits_to_add));
    }
    return res;
  }

  wrapped_interval_t SExt(unsigned bits_to_add)  const {
    std::vector<wrapped_interval_t> intervals;
    signed_split(intervals);
    
    wrapped_interval_t res = wrapped_interval_t::bottom();
    for (typename std::vector<wrapped_interval_t>::iterator it = intervals.begin(),
	   et = intervals.end();
	 it!=et; ++it) {
      // this should not happen
      if ((*it).is_bottom() || (*it).is_top()) continue;
      
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
    if (_start.ashr(wrapint(bits_to_keep,w)) == _stop.ashr(wrapint(bits_to_keep,w))) {
      wrapint lower_start = _start.keep_lower(bits_to_keep);
      wrapint lower_stop = _stop.keep_lower(bits_to_keep);
      if (lower_start <= lower_stop) {
	return wrapped_interval_t(lower_start, lower_stop);
      }
    } else {
      // note that _start is a wrapint so after +1 it can wraparound
      wrapint y(_start.ashr(wrapint(bits_to_keep, w)));
      ++y;
      if (y == _stop.ashr(wrapint(bits_to_keep,w))) {
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
    if (is_bottom()) return *this;
    
    // only if shift is constant    
    if (x.is_singleton()) {
      return Shl(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }
  
  wrapped_interval_t LShr(wrapped_interval_t  x) const {
    if (is_bottom()) return *this;
 
    // only if shift is constant
    if (x.is_singleton()) {
      return LShr(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }
  
  wrapped_interval_t AShr(wrapped_interval_t  x) const {
    if (is_bottom()) return *this;
 
    // only if shift is constant    
    if (x.is_singleton()) {
      return AShr(x.start().get_uint64_t());
    } else {
      return wrapped_interval_t::top();
    }
  }
  
  wrapped_interval_t And(wrapped_interval_t x) const
  { return default_implementation(x); }    
    
  wrapped_interval_t Or(wrapped_interval_t x)  const
  { return default_implementation(x); }    
    
  wrapped_interval_t Xor(wrapped_interval_t x) const
  { return default_implementation(x); }    
    

  void write(crab::crab_os& o) const {
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
	o << "[[" << (int) x << ", " << (int) y << "]]_";
      } else if (get_bitwidth(__LINE__) == 8) {
	o << "[[" << (int) static_cast<signed char>(x) << ", "
	          << (int) static_cast<signed char>(y) << "]]_";	
      } else {
	o << "[[" << _start.get_signed_bignum() << ", "
	  << _stop.get_signed_bignum() << "]]_";
      }
      #else
      o << "[[" << _start << ", " << _stop << "]]_";
      #endif
      o << (int) get_bitwidth(__LINE__);
    }
  }    
  
};

template<typename Number>
inline crab::crab_os& operator<<(crab::crab_os& o, const wrapped_interval<Number> &i) {
  i.write(o);
  return o;
}
}// end namespace
}// end namespace

namespace ikos{ 
namespace linear_interval_solver_impl {
  typedef crab::domains::wrapped_interval<z_number> z_wrapped_interval_t;
  typedef crab::domains::wrapped_interval<q_number> q_wrapped_interval_t;

  template<>
  inline z_wrapped_interval_t
  mk_interval(z_number c, typename crab::wrapint::bitwidth_t w)
  { return z_wrapped_interval_t(c, w); }

  template<>
  inline q_wrapped_interval_t
  mk_interval(q_number c, typename crab::wrapint::bitwidth_t w)
  { return q_wrapped_interval_t(c, w); }
  
  template<>
  inline z_wrapped_interval_t trim_interval(z_wrapped_interval_t i, z_wrapped_interval_t j) {
    if (i.is_bottom()) return i;
    // XXX: TODO: gamma(top()) \ gamma(j) can be expressed in a
    //            wrapped interval.
    if (i.is_top()) return i;
    if (!j.is_singleton()) return i;
    
    crab::wrapint k = j.start();
    if (i.start() == k) {
      if (i.is_singleton()){
	return z_wrapped_interval_t::bottom();
      }
      crab::wrapint k_plus(k);
      ++k_plus;
      z_wrapped_interval_t trimmed_res = z_wrapped_interval_t(k_plus, i.stop());
      return trimmed_res;
    } else if (i.stop() == k) {
      if (i.is_singleton()){
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
  
  template<>
  inline q_wrapped_interval_t trim_interval(q_wrapped_interval_t i, q_wrapped_interval_t j) { 
    // No refinement possible for disequations over rational numbers
    return i;
  }

  template<>
  inline z_wrapped_interval_t lower_half_line(z_wrapped_interval_t i, bool is_signed) {
    return i.lower_half_line(is_signed);
  }
  
  template<>
  inline q_wrapped_interval_t lower_half_line(q_wrapped_interval_t i, bool is_signed) {
    return i.lower_half_line(is_signed);    
  }
  
  template<>
  inline z_wrapped_interval_t upper_half_line(z_wrapped_interval_t i, bool is_signed) {
    return i.upper_half_line(is_signed);        
  }
  
  template<>
  inline q_wrapped_interval_t upper_half_line(q_wrapped_interval_t i, bool is_signed) {
    return i.upper_half_line(is_signed);        
  }
  
} // namespace linear_interval_solver_impl
} // end namespace ikos

namespace crab{
namespace domains {
  
template<typename Number, typename VariableName, std::size_t max_reduction_cycles = 10>
class wrapped_interval_domain:
    public crab::domains::
    abstract_domain<Number, VariableName,
		    wrapped_interval_domain<Number,VariableName,max_reduction_cycles> >  {
public:
  typedef wrapped_interval_domain<Number, VariableName, max_reduction_cycles> wrapped_interval_domain_t;
  typedef crab::domains::abstract_domain<Number,VariableName,wrapped_interval_domain_t> abstract_domain_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::number_t;
  using typename abstract_domain_t::varname_t;
  typedef interval<Number> interval_t;
  typedef wrapped_interval<Number> wrapped_interval_t;
  typedef typename wrapped_interval_t::bitwidth_t bitwidth_t;
  
private:
  typedef separate_domain<variable_t, wrapped_interval_t> separate_domain_t;
  typedef linear_interval_solver<Number, VariableName, separate_domain_t> solver_t;
  
public:
  typedef typename separate_domain_t::iterator iterator;
  
private:
  separate_domain_t _env;
  
  wrapped_interval_domain(separate_domain_t env): _env(env) { }

  void add(linear_constraint_system_t csts, std::size_t threshold = max_reduction_cycles) {
    if (!this->is_bottom()) {
      solver_t solver(csts, threshold);
      solver.run(this->_env);
    }
  }

  wrapped_interval_t eval_expr(linear_expression_t expr, bitwidth_t width) {
    wrapped_interval_t r(wrapint(expr.constant(), width));
    for (typename linear_expression_t::iterator it = expr.begin(); 
  	 it != expr.end(); ++it) {
      wrapped_interval_t c(wrapint(it->first,  width));
      // eval_expr should be "const" but operator[] in _env is not marked as "const"
      r += c * this->_env[it->second];
    }
    return r;
  }
  
public:
  static wrapped_interval_domain_t top() {
    return wrapped_interval_domain(separate_domain_t::top());
  }
  
  static wrapped_interval_domain_t bottom() {
    return wrapped_interval_domain(separate_domain_t::bottom());
  }
  
public:
  wrapped_interval_domain(): _env(separate_domain_t::top()) { }

  wrapped_interval_domain(const wrapped_interval_domain_t& e): 
    _env(e._env) { 
    crab::CrabStats::count (getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }
  
  wrapped_interval_domain_t& operator=(const wrapped_interval_domain_t& o) {
    crab::CrabStats::count (getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }
  
  iterator begin() {
    return this->_env.begin();
  }
  
  iterator end() {
    return this->_env.end();
  }
  
  bool is_bottom() {
    return this->_env.is_bottom();
  }
  
  bool is_top() {
    return this->_env.is_top();
  }

  bool operator<=(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    //CRAB_LOG("wrapped-int",
    //       crab::outs()<< *this << " <= " << e << "=";);
    bool res = (this->_env <= e._env);
    //CRAB_LOG("wrapped-int",
    //	     crab::outs () << (res ? "yes": "not") << "\n";);
    return res;    
  }
  
  void operator|=(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    CRAB_LOG("wrapped-int",
	     crab::outs() << *this << " U " << e << " = ");
    this->_env = this->_env | e._env;
    CRAB_LOG("wrapped-int", crab::outs() << *this << "\n";);
  }

  wrapped_interval_domain_t operator|(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    CRAB_LOG("wrapped-int",
	     crab::outs() << *this << " U " << e << " = ");
    wrapped_interval_domain_t res(this->_env | e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t operator&(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");
    CRAB_LOG("wrapped-int",
	     crab::outs() << *this << " n " << e << " = ");    
    wrapped_interval_domain_t res(this->_env & e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }

  wrapped_interval_domain_t operator||(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    CRAB_LOG("wrapped-int",
	     crab::outs() << "WIDENING " << *this << " and " << e << " = ");    
    wrapped_interval_domain_t res(this->_env || e._env);
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }
  
  template<typename Thresholds>
  wrapped_interval_domain_t widening_thresholds (wrapped_interval_domain_t e, const Thresholds &ts) {
    crab::CrabStats::count (getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    CRAB_LOG("wrapped-int",
	     crab::outs() << "WIDENING " << *this << " and " << e << " = ");    
    wrapped_interval_domain_t res(this->_env.widening_thresholds (e._env, ts));
    CRAB_LOG("wrapped-int", crab::outs() << res << "\n";);
    return res;
  }
  
  wrapped_interval_domain_t operator&&(wrapped_interval_domain_t e) {
    crab::CrabStats::count (getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    return (this->_env && e._env);
    }
  

  void operator-=(variable_t v) {
    crab::CrabStats::count (getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
    this->_env -= v;
  }

  template<typename Iterator>
  void project(Iterator it, Iterator et) {
    separate_domain_t projected_env = separate_domain_t::top();
    for (; it!=et; ++it) {
      variable_t v = *it;
      projected_env.set(v, this->_env[v]);
    }
    std::swap(this->_env, projected_env);
  }


  void expand(variable_t x, variable_t new_x) {
    set(new_x, this->_env[x]);
  }

  void set(variable_t v, wrapped_interval_t i) {
    crab::CrabStats::count (getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, i);
    CRAB_LOG("wrapped-int",
	     crab::outs() << v << ":=" << i << "=" << _env[v] << "\n");    
  }

  void set(variable_t v, interval_t i) {
    crab::CrabStats::count (getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    if (i.lb().is_finite() && i.ub.is_finite()) {
      wrapint start(i.lb(), v.get_bitwidth());
      wrapint stop(i.ub(), v.get_bitwidth());
      this->_env.set(v, wrapped_interval_t(start, stop));
      CRAB_LOG("wrapped-int",
	       crab::outs() << v << ":=" << i << "=" << _env[v] << "\n");    
    } else {
      CRAB_WARN("ignored assignment of an open interval in wrapped interval domain");
      *this -= v;
    }
  }
  
  void set(variable_t v, Number n) {
    crab::CrabStats::count (getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    wrapint start(n, v.get_bitwidth());
    wrapint stop(n, v.get_bitwidth());
    this->_env.set(v, wrapped_interval_t(start, stop));
    CRAB_LOG("wrapped-int",
	     crab::outs() << v << ":=" << n << "=" << _env[v] << "\n");    
  }

  // Return unlimited interval
  interval_t operator[](variable_t v) {
    return this->_env[v].to_interval();
  }

  // Return wrapped interval
  wrapped_interval_t get_wrapped_interval(variable_t v) {
    return this->_env[v];
  }
  
  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count (getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    if (boost::optional<variable_t> v = e.get_variable ()) {
      this->_env.set(x, this->_env [*v]);
    } else {
      wrapped_interval_t r = eval_expr(e, x.get_bitwidth());
      this->_env.set(x, r);
    }
    CRAB_LOG("wrapped-int",
	     crab::outs() << x << ":=" << e << "=" << _env[x] << "\n");        
  }
  
  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count (getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");
    
    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi = this->_env[z];
    wrapped_interval_t xi = wrapped_interval_t::bottom();
    
    switch (op) {
    case OP_ADDITION: {
      xi = yi + zi;
      break;
    }
    case OP_SUBTRACTION: {
      xi = yi - zi;
      break;
    }
    case OP_MULTIPLICATION: {
      xi = yi * zi;
      break;
    }
    case OP_DIVISION: {
      xi = yi / zi;
      break;
    }
    }
    this->_env.set(x, xi);
    CRAB_LOG("wrapped-int",
	     crab::outs() << x << ":=" << y << " " << op << " " << z << "=" << _env[x] << "\n");        
  }
  
  void apply(operation_t op, variable_t x, variable_t y, Number k) {
    crab::CrabStats::count (getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");
    
    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi(wrapint(k, x.get_bitwidth()), wrapint(k,x.get_bitwidth()));
    wrapped_interval_t xi = wrapped_interval_t::bottom();
    
    switch (op) {
    case OP_ADDITION: {
      xi = yi + zi;
      break;
    }
    case OP_SUBTRACTION: {
      xi = yi - zi;
      break;
    }
    case OP_MULTIPLICATION: {
      xi = yi * zi;
      break;
    }
    case OP_DIVISION: { // signed division
      xi = yi / zi;
      break;
    }
    }
    this->_env.set(x, xi);
    CRAB_LOG("wrapped-int",
	     crab::outs() << x << ":=" << y << " " << op << " " << k << "="
	                  << _env[x] << "\n");        
  }
  
  void backward_assign (variable_t x, linear_expression_t e,
			wrapped_interval_domain_t inv) {
    // crab::domains::BackwardAssignOps<wrapped_interval_domain_t>::
    // 	assign (*this, x, e, inv);
  }      
  
  void backward_apply (operation_t op,
		       variable_t x, variable_t y, Number z,
		       wrapped_interval_domain_t inv) {
    // crab::domains::BackwardAssignOps<wrapped_interval_domain_t>::
    // 	apply(*this, op, x, y, z, inv);
  }      
  
  void backward_apply(operation_t op,
		      variable_t x, variable_t y, variable_t z,
		      wrapped_interval_domain_t inv) {
    // crab::domains::BackwardAssignOps<wrapped_interval_domain_t>::
    // 	apply(*this, op, x, y, z, inv);
  }
  
  // cast_operators_api
  
  void apply(crab::domains::int_conv_operation_t op, variable_t dst, variable_t src){

    wrapped_interval_t src_i = this->_env[src];
    wrapped_interval_t dst_i;

    if (src_i.is_bottom() || src_i.is_top()) {
      dst_i = src_i;
    } else {
      switch (op) {
      case crab::domains::OP_ZEXT:
      case crab::domains::OP_SEXT:
	{
	  if (dst.get_bitwidth() < src.get_bitwidth()) {
	    CRAB_ERROR("destination must be larger than source in sext/zext");
	  }
	  unsigned bits_to_add = dst.get_bitwidth() - src.get_bitwidth();
	  dst_i = (op == crab::domains::OP_SEXT ?
		   src_i.SExt(bits_to_add) : src_i.ZExt(bits_to_add));
	}
	break;
      case crab::domains::OP_TRUNC:
	{
	  if (src.get_bitwidth() < dst.get_bitwidth()) {
	    CRAB_ERROR("destination must be smaller than source in truncate");
	  }
	  unsigned bits_to_keep = dst.get_bitwidth();
	  wrapped_interval_t dst_i;
	  dst_i = src_i.Trunc(bits_to_keep);
	}
	break;
      default:
	CRAB_ERROR("unexpected operation: ", op);
      }
    }
    set(dst, dst_i);
  }
  
  // bitwise_operators_api
  
  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z){
    crab::CrabStats::count (getDomainName() + ".count.apply");
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
  
  void apply(bitwise_operation_t op, variable_t x, variable_t y, Number k){
    crab::CrabStats::count (getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");
    
    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi(wrapint(k, x.get_bitwidth()));
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
  
  // division_operators_api
  
  void apply(div_operation_t op, variable_t x, variable_t y, variable_t z){
    crab::CrabStats::count (getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");
    
    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi = this->_env[z];
    wrapped_interval_t xi = wrapped_interval_t::bottom();
    
    switch (op) {
      case OP_SDIV: {
    	xi = yi / zi;
    	break;
      }
      case OP_UDIV: {
    	xi = yi.UDiv(zi);
    	break;
      }
      case OP_SREM: {
    	xi = yi.SRem(zi);
    	break;
      }
      case OP_UREM: {
    	xi = yi.URem(zi);
    	break;
      }
      default: 
        CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }
  
  void apply(div_operation_t op, variable_t x, variable_t y, Number k){
    crab::CrabStats::count (getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");
    
    wrapped_interval_t yi = this->_env[y];
    wrapped_interval_t zi(wrapint(k, x.get_bitwidth()));
    wrapped_interval_t xi = wrapped_interval_t::bottom();
    switch (op) {
      case OP_SDIV: {
    	xi = yi / zi;
    	break;
      }
      case OP_UDIV: {
    	xi = yi.UDiv(zi);
    	break;
      }
      case OP_SREM: {
    	xi = yi.SRem(zi);
    	break;
      }
      case OP_UREM: {
    	xi = yi.URem(zi);
    	break;
      }
      default: 
        CRAB_ERROR("unreachable");
    }
    this->_env.set(x, xi);
  }

  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count (getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");
    this->add(csts);
    CRAB_LOG("wrapped-int",
	     crab::outs() << "Added " << csts << " = " << *this << "\n");        
    
  }
  
  wrapped_interval_domain_t operator+(linear_constraint_system_t csts) {
    wrapped_interval_domain_t e(this->_env);
    e += csts;
    return e;
  }
  
  void write(crab::crab_os& o) {
    this->_env.write(o);
  }

  // Important: we make the choice here that we interpret wrapint as
  // signed mathematical integers.
  linear_constraint_system_t to_linear_constraint_system () {
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
  
  static std::string getDomainName () {
    return "Wrapped Intervals";
  }
  
}; // class wrapped_interval_domain


template<typename Number, typename VariableName>
class domain_traits <wrapped_interval_domain<Number,VariableName> > {
public:

  typedef wrapped_interval_domain<Number,VariableName> wrapped_interval_domain_t;
  typedef ikos::variable<Number, VariableName> variable_t;
  
  template<class CFG>
  static void do_initialization (CFG cfg) { }

  static void expand (wrapped_interval_domain_t& inv, variable_t x, variable_t new_x) {
    inv.expand(x, new_x);
  }
  
  static void normalize (wrapped_interval_domain_t& inv) {}
  
  template <typename Iter>
  static void forget (wrapped_interval_domain_t& inv, Iter it, Iter end){
    for(;it!=end; ++it) {
      inv -= *it;
    }
  }
  
  template <typename Iter>
  static void project (wrapped_interval_domain_t& inv, Iter it, Iter end) {
    inv.project(it, end);
  }
  
};
  
}
}


namespace crab {
namespace domains {

// Simple lattice to represent which limits (if any) have been crossed
// by a wrapped interval.
class wrapped_interval_limit_value: public ikos::writeable {
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
    NC  = 0x0,
    CS  = 0x1,
    CU  = 0x2,
    CSU = 0x3 /*top*/
  } kind_t;
  
  kind_t _value;
  bool   _is_bottom;
  
  wrapped_interval_limit_value(kind_t v, bool is_bottom)
    : _value(v), _is_bottom(is_bottom) { }
  
public:
  
  wrapped_interval_limit_value(): _value(CSU), _is_bottom(false) {}
  
  static wrapped_interval_limit_value bottom() {
    return wrapped_interval_limit_value(NC/*any value*/, true);
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
  
  template<typename N>
  static wrapped_interval_limit_value convert(const wrapped_interval<N>& i) {
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
  
  wrapped_interval_limit_value(const wrapped_interval_limit_value& o)
    : _value(o._value), _is_bottom(o._is_bottom) {}
  
  wrapped_interval_limit_value& operator=(const wrapped_interval_limit_value& o) {
    if (this != &o) {
      _value = o._value;
      _is_bottom = o._is_bottom;
    }
    return *this;
  }
  
  bool is_bottom() const { return _is_bottom;}
  
  // the wrapped interval might have crossed both limits
  bool is_top() const { return !_is_bottom && _value == CSU;}
  
  // the wrapped interval might have crossed the signed limit
  bool is_crossing_signed_limit() const
  { return (!is_bottom() && (_value == CS || _value == CSU)); }
  
  // the wrapped interval might have crossed the unsigned limit  
  bool is_crossing_unsigned_limit() const
  { return (!is_bottom() && (_value == CU || _value == CSU)); }
  
  bool operator<=(const wrapped_interval_limit_value& o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (o.is_bottom()) {
      return false;
    } else if (is_top()) {
      return o.is_top();
    } else if (_value == CS) {
      return (o._value == CS || o.is_top());
    } else if (_value == CU) {
      return (o._value == CU || o.is_top());
    }
    
    assert(false && "unreachable");   
    return false;
  }
  
  bool operator==(const wrapped_interval_limit_value& o) const {
    return (_value == o._value && is_bottom() == o.is_bottom());
  }
    
  wrapped_interval_limit_value operator|(wrapped_interval_limit_value o) {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      return wrapped_interval_limit_value(static_cast<kind_t>(
	     static_cast<int>(_value) | static_cast<int>(o._value)), false);
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
      return wrapped_interval_limit_value(static_cast<kind_t>(
	     static_cast<int>(_value) & static_cast<int>(o._value)), false);
    }
  }
  
  // the lattice satisfy DCC so meet is the narrowing
  wrapped_interval_limit_value operator&&(wrapped_interval_limit_value o) { 
    return this->operator&(o); 
  }
  
  void write(crab_os& o) {
    if (is_bottom()) {
      o << "_|_";
    } else {
      switch (_value) {
      case NC:  o << "no-cross"; break;
      case CS:  o << "cross-signed"; break;
      case CU:  o << "cross-unsigned"; break;
      default:/*top*/ o << "top";  
      }
    }
  }
};
    
/** 
    Wrapped interval domain augmented with an abstraction of the
    execution history: it keeps track of which variable crossed which
    signed/unsigned limits.

    The only case where the history of a variable is reset is when it
    is assigned to a constant value.
**/
template <typename Number, typename VariableName, std::size_t max_reduction_cycles = 10>
class wrapped_interval_limits_domain:
    public abstract_domain<Number, VariableName,
	   wrapped_interval_limits_domain<Number, VariableName, max_reduction_cycles> > {
  
  typedef wrapped_interval_limits_domain<Number, VariableName, max_reduction_cycles> this_type;
  typedef abstract_domain<Number,VariableName,this_type> abstract_domain_t;
  
public:
    
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef interval<number_t> interval_t;
  typedef wrapped_interval<number_t> wrapped_interval_t;  
  
private:

  typedef wrapped_interval_domain<Number, VariableName, max_reduction_cycles>
  wrapped_interval_domain_t;

  typedef separate_domain<variable_t, wrapped_interval_limit_value> separate_domain_t; 
  typedef discrete_domain<variable_t> discrete_domain_t;
  typedef typename linear_constraint_system_t::variable_set_t variable_set_t;
  
  wrapped_interval_domain_t _w_int_dom;
  // Map each variable to which limit was crossed.
  separate_domain_t _limit_env;
  // Set of may-initialized variables
  discrete_domain_t _init_set;

  wrapped_interval_limits_domain(const wrapped_interval_domain_t& dom,
				 const separate_domain_t& limit_env,
				 const discrete_domain_t& init_set)
    : _w_int_dom(dom), _limit_env(limit_env), _init_set(init_set) {}

  inline bool may_be_initialized(variable_t x)
  { return (discrete_domain_t(x) <= _init_set); }
  
  // This hoare triple holds {x = old_i} x := ... { x = new_i}
  inline void update_limits(variable_t x, wrapped_interval_t old_i, wrapped_interval_t new_i) {
    wrapped_interval_limit_value old_l = wrapped_interval_limit_value::bottom();
    wrapped_interval_limit_value new_l;
    
    if (may_be_initialized(x)) {
      old_l = _limit_env[x];
      // XXX: it's not enough to convert only new_i. E.g., char x = 127; x++;
      //  convert([127,127])   = no-cross
      //  convert([-128,-128]) = no-cross
      new_l = wrapped_interval_limit_value::convert(old_i | new_i);
    } else {
      new_l = wrapped_interval_limit_value::convert(new_i);
    }

    // -- weak update to keep past history
    _limit_env.set(x, old_l | new_l);

    // -- mark x as initialized
    _init_set += x;
  }

  void update_limits(const variable_set_t &vars,
		     const std::vector<wrapped_interval_t>& old_intervals,
		     const std::vector<wrapped_interval_t>& new_intervals) {
    assert(vars.size() == old_intervals.size());
    assert(old_intervals.size() == new_intervals.size());
    
    unsigned i=0;
    for(auto const& v: vars) {
      update_limits(v, old_intervals[i], new_intervals[i]);
      ++i;
    }
  }
		     
public:
  
  static this_type top() {
    return this_type(wrapped_interval_domain_t::top(),
		     separate_domain_t::top(),
		     discrete_domain_t::bottom() /*empty set*/);
  }
  
  static this_type bottom() {
    return this_type(wrapped_interval_domain_t::bottom(),
		     separate_domain_t::bottom(),
		     discrete_domain_t::bottom() /*empty set*/);
  }
  
  wrapped_interval_limits_domain()
    : _w_int_dom(), _limit_env(), _init_set(discrete_domain_t::bottom() /*empty set*/) {}
    
  wrapped_interval_limits_domain(const this_type& o)
    : _w_int_dom(o._w_int_dom),
      _limit_env(o._limit_env),
      _init_set(o._init_set) {}

  wrapped_interval_limits_domain(const this_type&& o)
    : _w_int_dom(std::move(o._w_int_dom)),
      _limit_env(std::move(o._limit_env)),
      _init_set(std::move(o._init_set)) {}
  
  this_type& operator=(const this_type& o) {
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
    return (_w_int_dom <= o._w_int_dom &&
	    _limit_env <= o._limit_env);
  }
  
  bool operator==(this_type o) {
    return (*this <= o && o <= *this);
  }
  
  void operator|=(this_type o) {
    _w_int_dom |= o._w_int_dom;
    _limit_env = _limit_env | o._limit_env;
    _init_set = _init_set | o._init_set;
  }
  
  this_type operator|(this_type o) {
    return this_type(_w_int_dom | o._w_int_dom,
		     _limit_env | o._limit_env,
		     _init_set  | o._init_set);
  }

  this_type operator&(this_type o) {
    return this_type(_w_int_dom & o._w_int_dom,
		     _limit_env & o._limit_env,
		     _init_set  & o._init_set);    
  }

  this_type operator||(this_type o) {
    return this_type(_w_int_dom || o._w_int_dom,
		     _limit_env || o._limit_env,
		     _init_set  || o._init_set);    
  }

  this_type operator&&(this_type o) {
    return this_type(_w_int_dom && o._w_int_dom,
		     _limit_env && o._limit_env,
		     _init_set  && o._init_set);
  }
  
  template<typename Thresholds>
  this_type widening_thresholds(this_type o, const Thresholds& ts) {
    return this_type(_w_int_dom.widening_thresholds (o._w_int_dom, ts),
		     _limit_env || o._limit_env,
		     _init_set  || o._init_set);
  } 

  void set (variable_t x, interval_t i) {
    _w_int_dom.set(x, i);    
    if (i.singleton()) {
      // XXX: x's history is reset
      _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());            
    } else {
      CRAB_WARN("TODO: set operation with unlimited interval");
      _limit_env -= x;
    }
    _init_set += x;

    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << i << " => " << *this<<"\n";);
  }
  
  void set (variable_t x, wrapped_interval_t i) {
    _w_int_dom.set(x, i);
    // XXX: x's history is reset    
    _limit_env.set(x, wrapped_interval_limit_value::convert(i));
    _init_set += x;
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << i << " => " << *this<<"\n";);
  }

  void set (variable_t x, Number n) {
    _w_int_dom.set(x, n);
    // XXX: x's history is reset    
    _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());
    _init_set += x; 
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << n << " => " << *this<<"\n";);
  }
    
  interval_t operator[](variable_t v) {
    return _w_int_dom[v];
  }

  wrapped_interval_t get_wrapped_interval(variable_t v) {
    return _w_int_dom.get_wrapped_interval(v);
  }
  
  wrapped_interval_limit_value get_limit_value(variable_t x) const {
    return _limit_env[x];
  }
  
  void operator-=(variable_t v) {
    _w_int_dom -= v;
    _limit_env -= v;
    // XXX: we never remove a variable from _init_set.
    //      This avoids, e.g. to mark as uninitialized a variable that
    //      have been havoc'ed.
    // _init_set -= v;
  }
  
  // numerical_domains_api
  
  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, z);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << z
	                  << " => " << *this<<"\n";);
  }
  
  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, k);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << k
	                  << " => " << *this <<"\n";);    
  }
  
  void assign(variable_t x, linear_expression_t e) {
    if (e.is_constant()) {
      // XXX: x's history is reset          
      _w_int_dom.assign(x, e);
      _limit_env.set(x, wrapped_interval_limit_value::do_not_cross());
      _init_set += x; 
    } else {
      wrapped_interval_t old_i = get_wrapped_interval(x);
      _w_int_dom.assign(x, e);
      wrapped_interval_t new_i = get_wrapped_interval(x);    
      update_limits(x, old_i, new_i);
    }
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << e << " => " << *this<<"\n";);    
  }
  
  void backward_assign (variable_t x, linear_expression_t e, this_type invariant) {
    _w_int_dom.backward_assign(x, e, invariant._w_int_dom);
    // XXX: ignore _limit_env
  }
  
  void backward_apply (operation_t op, variable_t x, variable_t y, number_t z,
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
    for (auto v: variables) {
      old_intervals.push_back(get_wrapped_interval(v));
    }
    
    _w_int_dom += csts;

    std::vector<wrapped_interval_t> new_intervals;
    new_intervals.reserve(variables.size());
    for (auto v: variables) {
      new_intervals.push_back(get_wrapped_interval(v));
    }

    update_limits(variables, old_intervals, new_intervals);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << "assume(" << csts << ") => " << *this << "\n";);        
  }
  
  // cast_operators_api
  
  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    wrapped_interval_t old_i = get_wrapped_interval(dst);
    _w_int_dom.apply(op, dst, src);
    wrapped_interval_t new_i = get_wrapped_interval(dst);    
    update_limits(dst, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << dst << ":=" << op << " " << src << " => " << *this<<"\n";);    
  }
      
  // bitwise_operators_api
  
  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, z);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << z
	                  << " => " << *this<<"\n";);
  }
  
  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, k);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << k
	                  << " => " << *this<<"\n";);    
  }
      
  // division_operators_api
  
  void apply(div_operation_t op, variable_t x, variable_t y, variable_t z) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, z);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << z
	                  << " => " << *this<<"\n";);    
  }
  
  void apply(div_operation_t op, variable_t x, variable_t y, number_t k) {
    wrapped_interval_t old_i = get_wrapped_interval(x);
    _w_int_dom.apply(op, x, y, k);
    wrapped_interval_t new_i = get_wrapped_interval(x);    
    update_limits(x, old_i, new_i);
    CRAB_LOG("wrapped-int2",
	     crab::outs() << x << ":=" << y << " " << op << " " << k
	                  << " => " << *this<<"\n";);    
  }
  
  void write(crab_os& o) {
    //o << "(" << _w_int_dom << "," << _limit_env << "," << _init_set << ")";
    o << "(" << _w_int_dom << "," << _limit_env << ")";
  }
  
  linear_constraint_system_t to_linear_constraint_system() {
    return _w_int_dom.to_linear_constraint_system();
  }
      
  static std::string getDomainName() 
  { return "Wrapped Intervals with abstraction of execution history"; }
  
  // domain_traits_api
  
  void expand(variable_t x, variable_t new_x) {
    _w_int_dom.expand(x, new_x);
    _limit_env.set(new_x, _limit_env[x]);
    if (may_be_initialized(x)) {
      _init_set += new_x;
    }
  }
  
  template <typename Range>
  void project(Range vars) {
    _w_int_dom.project(vars.begin(), vars.end());

    separate_domain_t projected_env = separate_domain_t::top();
    discrete_domain_t projected_init_set = discrete_domain_t::bottom();    
    for (variable_t v: vars) {
      projected_env.set(v, _limit_env[v]);
      if (may_be_initialized(v)) {
	projected_init_set += v;
      }
    }
    std::swap(_limit_env, projected_env);
    std::swap(_init_set, projected_init_set);
  }
  
}; 
  

template<typename Number, typename VariableName>
class domain_traits <wrapped_interval_limits_domain<Number,VariableName> > {
public:

  typedef wrapped_interval_limits_domain<Number,VariableName> wrapped_interval_domain_t;
  typedef ikos::variable<Number, VariableName> variable_t;
  
  template<class CFG>
  static void do_initialization (CFG cfg) { }

  static void expand (wrapped_interval_domain_t& inv, variable_t x, variable_t new_x) {
    inv.expand(x, new_x);
  }
  
  static void normalize (wrapped_interval_domain_t& inv) {}
  
  template <typename Iter>
  static void forget (wrapped_interval_domain_t& inv, Iter it, Iter end){
    for(;it!=end; ++it) {
      inv -= *it;
    }
  }
  
  template <typename Iter>
  static void project (wrapped_interval_domain_t& inv, Iter it, Iter end) {
    inv.project(it, end);
  }
  
};

}
}
