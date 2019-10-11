#include <crab/numbers/bignums.hpp>
#include <crab/common/debug.hpp>


#include <gmpxx.h>
#include <climits>
#include <string>

namespace ikos {

z_number::z_number(mpz_class n) : _n(n) {}

z_number z_number::from_ulong(unsigned long n) {
  mpz_class b(n);
  return z_number(b);
}

z_number z_number::from_slong(signed long n) {
  mpz_class b(n);
  return z_number(b);
}
  
z_number::operator long() const { 
  if (_n.fits_slong_p ()) {
    return _n.get_si ();
  }
    else {
      CRAB_ERROR("mpz_class ", _n.get_str(), " does not fit into a signed long integer");
    }
} 

z_number::operator int() const { 
  if (_n.fits_sint_p ()) {
    // get_si returns a signed long so we cast it to int
      return (int) _n.get_si ();
  }
  else {
    CRAB_ERROR("mpz_class ", _n.get_str(), " does not fit into a signed integer");
  }
} 

z_number::operator mpz_class() const { 
  return _n;
} 

z_number::z_number(): _n(0) {}

z_number::z_number(std::string s): _n(s) {
  /// We want to be compiled with -fno_exceptions  
  // try {
  //   _n = s;
  // } catch (std::invalid_argument& e) {
  //   CRAB_ERROR ("z_number: invalid string in constructor", s);
  // }
}

z_number::z_number(signed long long int n) : _n((signed long int) n) {
  if (n > LONG_MAX) {
    CRAB_ERROR(n, " cannot fit into a signed long int: use another mpz_class constructor");
  }
}

std::string z_number::get_str () const {
  return _n.get_str();
}

std::size_t z_number::hash() const {
  boost::hash<std::string> hasher;
  return hasher(_n.get_str());
}

bool z_number::fits_sint() const {
  return _n.fits_sint_p();
}

bool z_number::fits_slong() const {
  return _n.fits_slong_p();
}

z_number z_number::operator+(z_number x) const {
  mpz_class r = _n + x._n;
    return z_number(r);
}

z_number z_number::operator*(z_number x) const {
  mpz_class r = _n * x._n;
  return z_number(r);
}

z_number z_number::operator-(z_number x) const {
  mpz_class r = _n - x._n;
  return z_number(r);
}

z_number z_number::operator-() const {
  mpz_class r = -_n;
  return z_number(r);
}

z_number z_number::operator/(z_number x) const {
  if (x._n == 0) {
    CRAB_ERROR("z_number: division by zero [1]");
  } else {
    mpz_class r = _n / x._n;
    return z_number(r);
  }
}

z_number z_number::operator%(z_number x) const {
  if (x._n == 0) {
    CRAB_ERROR("z_number: division by zero [2]");
  } else {
    mpz_class r = _n % x._n;
    return z_number(r);
  }
}

z_number& z_number::operator+=(z_number x) {
  _n += x._n;
  return *this;
}

z_number& z_number::operator*=(z_number x) {
  _n *= x._n;
  return *this;
}

z_number& z_number::operator-=(z_number x) {
  _n -= x._n;
  return *this;
}

z_number& z_number::operator/=(z_number x) {
  if (x._n == 0) {
    CRAB_ERROR("z_number: division by zero [3]");
  } else {
    _n /= x._n;
    return *this;
  }
}

z_number& z_number::operator%=(z_number x) {
  if (x._n == 0) {
    CRAB_ERROR("z_number: division by zero [4]");
  } else {
    _n %= x._n;
      return *this;
  }
}

z_number& z_number::operator--() {
  --(_n);
  return *this;
}

z_number& z_number::operator++() {
  ++(_n);
  return *this;
}

z_number z_number::operator++(int) {
  z_number r(*this);
  ++(*this);
  return r;
}

z_number z_number::operator--(int) {
  z_number r(*this);
  --(*this);
  return r;
}

bool z_number::operator==(z_number x) const { return _n == x._n; }

bool z_number::operator!=(z_number x) const { return _n != x._n; }

bool z_number::operator<(z_number x) const { return _n < x._n; }

bool z_number::operator<=(z_number x) const { return _n <= x._n; }

bool z_number::operator>(z_number x) const { return _n > x._n; }

bool z_number::operator>=(z_number x) const { return _n >= x._n; }

z_number z_number::operator&(z_number x) const { return z_number(_n & x._n); }

z_number z_number::operator|(z_number x) const { return z_number(_n | x._n); }

z_number z_number::operator^(z_number x) const { return z_number(_n ^ x._n); }

z_number z_number::operator<<(z_number x) const {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_mul_2exp(tmp, _n.get_mpz_t(), mpz_get_ui(x._n.get_mpz_t()));
  mpz_class result(tmp);
  return z_number(result);
}

z_number z_number::operator>>(z_number x) const {
  mpz_class tmp(_n);
  return z_number(tmp.operator>>=(mpz_get_ui(x._n.get_mpz_t())));
}

z_number z_number::fill_ones() const {
  assert(_n >= 0);
  if (_n == 0) {
    return z_number(0);
  }
  
  mpz_class result;
  for (result = 1; result < _n; result = 2 * result + 1)
    ;
  return z_number(result);
}

void z_number::write(crab::crab_os& o) const {
  o << _n.get_str();
}


q_number::q_number(): _n(0) {}

q_number::q_number(mpq_class n): _n(n) {}
  
q_number::q_number(std::string s): _n(s) {
  /// We want to be compiled with -fno_exceptions
  // try {
  //   _n = s;
  //   } catch (std::invalid_argument& e) {
  //   CRAB_ERROR("q_number: invalid string in constructor ",s);
  // }
  _n.canonicalize();
}

q_number::q_number(double n)
  : _n(n) {
  _n.canonicalize();
}

q_number::q_number(z_number n)
  : _n(n._n) {
  _n.canonicalize();
}

q_number::q_number(z_number n, z_number d)
  : _n(n._n, d._n) {
  _n.canonicalize();
}

q_number::operator mpq_class() const { 
  return _n;
} 

std::string q_number::get_str () const {
  return _n.get_str();
}

std::size_t q_number::hash() const {
  boost::hash<std::string> hasher;
  return hasher(_n.get_str());
}

q_number q_number::operator+(q_number x) const {
  mpq_class r = _n + x._n;
  return q_number(r);
}

q_number q_number::operator*(q_number x) const {
  mpq_class r = _n * x._n;
  return q_number(r);
}

q_number q_number::operator-(q_number x) const {
  mpq_class r = _n - x._n;
  return q_number(r);
}

q_number q_number::operator-() const {
  mpq_class r = -_n;
  return q_number(r);
}

q_number q_number::operator/(q_number x) const {
  if (x._n == 0) {
    CRAB_ERROR("q_number: division by zero [1]");
  } else {
    mpq_class r = _n / x._n;
    return q_number(r);
  }
}

q_number& q_number::operator+=(q_number x) {
  _n += x._n;
  return *this;
}

q_number& q_number::operator*=(q_number x) {
  _n *= x._n;
  return *this;
}

q_number& q_number::operator-=(q_number x) {
  _n -= x._n;
  return *this;
}

q_number& q_number::operator/=(q_number x) {
  if (x._n == 0) {
    CRAB_ERROR("q_number: division by zero [2]");
  } else {
    _n /= x._n;
    return *this;
  }
}

q_number& q_number::operator--() {
  --(_n);
  return *this;
}

q_number& q_number::operator++() {
  ++(_n);
  return *this;
}

q_number q_number::operator--(int) {
  q_number r(*this);
  --(*this);
  return r;
}

q_number q_number::operator++(int) {
  q_number r(*this);
  ++(*this);
  return r;
}

bool q_number::operator==(q_number x) const { return _n == x._n; }

bool q_number::operator!=(q_number x) const { return _n != x._n; }

bool q_number::operator<(q_number x) const { return _n < x._n; }

bool q_number::operator<=(q_number x) const { return _n <= x._n; }

bool q_number::operator>(q_number x) const { return _n > x._n; }

bool q_number::operator>=(q_number x) const { return _n >= x._n; }

z_number q_number::numerator() const { return z_number(_n.get_num()); }

z_number q_number::denominator() const { return z_number(_n.get_den()); }

z_number q_number::round_to_upper() const {
  z_number num = numerator();
  z_number den = denominator();
  z_number q = num / den;
  z_number r = num % den;
  if (r == 0 || *this < 0) {
    return q;
  } else {
    return q + 1;
  }
}

z_number q_number::round_to_lower() const {
  z_number num = numerator();
  z_number den = denominator();
  z_number q = num / den;
  z_number r = num % den;
  if (r == 0 || *this > 0) {
    return q;
  } else {
    return q - 1;
  }
}

void q_number::write(crab::crab_os& o) const {
  o << _n.get_str();
}

} // end namespace

