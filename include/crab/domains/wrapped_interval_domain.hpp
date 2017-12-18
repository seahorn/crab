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

  // return true if interval [0111...1, 1000....0] is included
  bool cross_signed_limits() const {
    wrapped_interval_t i(wrapint::get_signed_max(get_bitwidth()),
			 wrapint::get_signed_min(get_bitwidth()));
    return (i <= *this);
  }

  // return true if interval [1111...1, 0000....0] is included
  bool cross_unsigned_limits() const {
    wrapped_interval_t i(wrapint::get_unsigned_max(get_bitwidth()),
			 wrapint::get_unsigned_min(get_bitwidth()));
    return (i <= *this);
  }
  
  bitwidth_t get_bitwidth() const {
    if (is_bottom()) {
      CRAB_ERROR("get_bitwidth() cannot be called from a bottom element");
    } else if (is_top()) {
      CRAB_ERROR("get_bitwidth() cannot be called from a top element");
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

  interval<Number> to_interval() const {
    typedef interval<Number> interval_t;
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top() || (cross_signed_limits () || cross_unsigned_limits())) {
      return interval_t::top();
    } else {
      return interval_t(_start.get_bignum(), _stop.get_bignum());
    }
  }
  
  /** begin needed by interval constraint solver **/
  wrapped_interval_t lower_half_line() const {
    if (is_top() || is_bottom()) return *this;
    // FIXME: we assume here signed intervals
    wrapint smin = wrapint::get_signed_min(get_bitwidth());
    wrapped_interval_t res = wrapped_interval_t(smin, _stop);
    return res;
  }
  
  wrapped_interval_t upper_half_line() const {
    if (is_top() || is_bottom()) return *this;
    // FIXME: we assume here signed intervals
    wrapint smax = wrapint::get_signed_max(get_bitwidth());
    wrapped_interval_t res = wrapped_interval_t(_start, smax);
    return res;
  }

  boost::optional<Number> singleton() const {
    if (!is_bottom() && !is_top() && _start == _stop) {
      return Number(_start.get_str());
    } else {
      return boost::optional<Number>();
    }
  }
  /** end needed by interval constraint solver **/

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
	    return x;
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
	    return x;
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
    assert (get_bitwidth() == x.get_bitwidth());
    bitwidth_t w = x.get_bitwidth();
    if ((_stop - _start) >= wrapint::get_signed_min(w))
      return wrapped_interval_t::top();
    
    wrapped_interval_t join = *this | x;
    if (join == wrapped_interval_t(_start, x._stop)) {
      wrapped_interval_t doubled_old(_start, (_stop * wrapint(2,w)) - _start + wrapint(1,w));
      return join | doubled_old;
    } else if (join == wrapped_interval_t(x._start, _stop)) {
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
      // -- check if the addition will wrap
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
      // -- check if the subtraction will wrap
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
      // temporary special cases 
      wrapped_interval_t one(1, x.get_bitwidth());
      wrapped_interval_t minus_one(-1, x.get_bitwidth());          
      if (x == one)
	return *this;
      else if (*this == one)
	return x;
      else if (x == minus_one) {
	return -(*this);
      } else if (*this == minus_one) {
	return -x;
      } else {
	// TODOX
	CRAB_WARN("Skipped ", *this, " * ", x);
	return wrapped_interval_t::top();
      }
    }
  }
  
  wrapped_interval_t& operator*=(wrapped_interval_t x) {
    return this->operator=(this->operator*(x));
  }

  // signed division
  wrapped_interval_t operator/(wrapped_interval_t x) const {
    wrapped_interval_t one(1, x.get_bitwidth());
    wrapped_interval_t minus_one(-1, x.get_bitwidth());    
    // temporary special case;
    if (x == one)
      return *this;
    else if (x == minus_one) {
      return -(*this);
    } else {
      // TODOX
      CRAB_WARN("Skipped ", *this, " / ", x);      
      return default_implementation(x);
    }
  }
  
  wrapped_interval_t& operator/=(wrapped_interval_t x) {
    return this->operator=(this->operator/(x));
  }   

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
      if (get_bitwidth() == 32) {
	o << "[[" << (int) x << ", " << (int) y << "]]_";
      } else if (get_bitwidth() == 8) {
	o << "[[" << (int) static_cast<signed char>(x) << ", "
	          << (int) static_cast<signed char>(y) << "]]_";	
      } else {
	// TODO: print the wrapint as a signed number
	o << "[[" << _start << ", " << _stop << "]]_";
      }
      #else
      o << "[[" << _start << ", " << _stop << "]]_";
      #endif
      o << (int) get_bitwidth();
    }
  }    
    
  /** division and remainder operations **/
  wrapped_interval_t UDiv(wrapped_interval_t x ) const
  { return default_implementation(x); }
  
  wrapped_interval_t SRem(wrapped_interval_t x)  const
  { return default_implementation(x); }    
  
  wrapped_interval_t URem(wrapped_interval_t x)  const
  { return default_implementation(x); }    

  /** bitwise operations **/
  
  wrapped_interval_t And(wrapped_interval_t x) const
  { return default_implementation(x); }    
    
  wrapped_interval_t Or(wrapped_interval_t x)  const
  { return default_implementation(x); }    
    
  wrapped_interval_t Xor(wrapped_interval_t x) const
  { return default_implementation(x); }    
    
  wrapped_interval_t Shl(wrapped_interval_t x) const
  { return default_implementation(x); }    
  
  wrapped_interval_t LShr(wrapped_interval_t  x) const
  { return default_implementation(x); }    

  wrapped_interval_t AShr(wrapped_interval_t  x) const
  { return default_implementation(x); }    
  
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

  template<>
  inline crab::domains::wrapped_interval<z_number>
  trim_bound(crab::domains::wrapped_interval<z_number> i, z_number c) {
    // pre: c is a singleton
    
    if (i.is_bottom()) return i;
    // XXX: not sure if we can be more precise here.
    // trim_bound is used by the linear interval solver which knows whether
    // the interval is over signed or unsigned. In that case we can refine top.
    if (i.is_top()) return i;

    crab::wrapint k(c, i.get_bitwidth());
    if (i.start() == k) {
      return crab::domains::wrapped_interval<z_number>(k++, i.stop());
    } else if (i.stop() == k) {
      return crab::domains::wrapped_interval<z_number>(i.start(), k--);
    } else {
      return i;
    }
  }

  template<>
  inline crab::domains::wrapped_interval<q_number>
  trim_bound(crab::domains::wrapped_interval<q_number> i, q_number /* c */) { 
    // No refinement possible for disequations over rational numbers
    return i;
  }

  template<>
  inline crab::domains::wrapped_interval<z_number>
  mk_interval(z_number c, typename crab::wrapint::bitwidth_t w) {
    return crab::domains::wrapped_interval<z_number>(c, w);
  }

  template<>
  inline crab::domains::wrapped_interval<q_number>
  mk_interval(q_number c, typename crab::wrapint::bitwidth_t w) {
    return crab::domains::wrapped_interval<q_number>(c, w);
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
  
  interval_t operator[](variable_t v) {
    return this->_env[v].to_interval();
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
    // TODOX
    *this -= dst;
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
  
  linear_constraint_system_t to_linear_constraint_system () {
    linear_constraint_system_t csts;
    if (this->is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }
    
    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      variable_t v = it->first;
      wrapped_interval_t i = it->second;
      if (!i.is_top()) {	
	csts += linear_constraint_t(v >= Number(i.start().get_str()));
	csts += linear_constraint_t(v <= Number(i.stop().get_str()));
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
