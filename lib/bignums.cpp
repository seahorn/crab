#include <crab/common/debug.hpp>
#include <crab/numbers/bignums.hpp>

#include <limits>
#include <string>

namespace ikos {
namespace bignums_impl {
struct scoped_cstring {
  typedef void (*__gmp_freefunc_t)(void *, size_t);
  char *m_str;
  scoped_cstring(char *s) { m_str = s; }
  ~scoped_cstring() {
    __gmp_freefunc_t freefunc;
    mp_get_memory_functions(nullptr, nullptr, &freefunc);
    (*freefunc)(m_str, std::strlen(m_str) + 1);
  }
};
} // namespace bignums_impl

/* Wrapper for mpz */

z_number::operator int64_t() const {
  if (fits_sint()) {
    return (int64_t)mpz_get_si(_n);
  } else if (fits_int64()) {
    int64_t res = 0;
    mpz_export(&res, 0 /*NULL*/, 1, sizeof(int64_t), 0, 0, _n);
    return ((mpz_sgn(_n) < 0) ? -res : res);
  } else {
    CRAB_ERROR("z_number ", get_str(), " does not fit into int64_t");
  }
}

z_number::z_number() { mpz_init(_n); }

z_number::z_number(int64_t n) {
  if (n >= std::numeric_limits<signed long int>::min() &&
      n <= std::numeric_limits<signed long int>::max()) {
    mpz_init_set_si(_n, static_cast<signed long int>(n));
  } else {
    mpz_init(_n);
    mpz_import(_n, 1, 1, sizeof(int64_t), 0, 0, &n);
    if (n < 0) {
      mpz_neg(_n, _n);
    }
  }
}

z_number::z_number(const std::string &s, unsigned base) {
  int res = mpz_init_set_str(_n, s.c_str(), base);
  if (res == -1) {
    CRAB_ERROR("z_number: invalid string in constructor", s);
  }
}

z_number z_number::from_mpz_t(mpz_t n) {
  z_number z;
  mpz_set(z._n, n);
  return z;
}

z_number z_number::from_mpz_srcptr(mpz_srcptr n) {
  z_number z;
  mpz_set(z._n, n);
  return z;
}

z_number z_number::from_uint64(uint64_t n) {
  z_number z;
  if (n <= std::numeric_limits<unsigned long>::max()) {
    mpz_set_ui(z._n, static_cast<unsigned long>(n));
  } else {
    mpz_import(z._n, 1, 1, sizeof(uint64_t), 0, 0, &n);
  }
  return z;
}

z_number::z_number(const z_number &o) { mpz_init_set(_n, o._n); }

z_number::z_number(z_number &&o) {
  *_n = *o._n;
  mpz_init(o._n);
}

z_number &z_number::operator=(const z_number &o) {
  if (this != &o) {
    mpz_set(_n, o._n);
  }
  return *this;
}

z_number &z_number::operator=(z_number &&o) {
  if (this != &o) {
    std::swap(*_n, *o._n);
  }
  return *this;
}

z_number::~z_number() { mpz_clear(_n); }

std::string z_number::get_str(unsigned base) const {
  bignums_impl::scoped_cstring res(mpz_get_str(0, base, _n));
  return std::string(res.m_str);
}

std::size_t z_number::hash() const {
  boost::hash<std::string> hasher;
  return hasher(get_str());
}

bool z_number::fits_sint() const { return mpz_fits_sint_p(_n); }

bool z_number::fits_slong() const { return mpz_fits_slong_p(_n); }

bool z_number::fits_int64() const {
  return (*this >= std::numeric_limits<int64_t>::min() &&
          *this <= std::numeric_limits<int64_t>::max());
}

z_number z_number::operator+(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_add(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator*(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_mul(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator-(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_sub(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator-() const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_neg(mp_r, _n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator/(z_number x) const {
  if (x == 0) {
    CRAB_ERROR("z_number: division by zero [1]");
  } else {
    mpz_t mp_r;
    mpz_init(mp_r);
    mpz_tdiv_q(mp_r, _n, x._n);
    z_number res = from_mpz_t(mp_r);
    mpz_clear(mp_r);
    return res;
  }
}

z_number z_number::operator%(z_number x) const {
  if (x == 0) {
    CRAB_ERROR("z_number: division by zero [2]");
  } else {
    mpz_t mp_r;
    mpz_init(mp_r);
    mpz_tdiv_r(mp_r, _n, x._n);
    z_number res = from_mpz_t(mp_r);
    mpz_clear(mp_r);
    return res;
  }
}

z_number &z_number::operator+=(z_number x) {
  mpz_add(_n, _n, x._n);
  return *this;
}

z_number &z_number::operator*=(z_number x) {
  mpz_mul(_n, _n, x._n);
  return *this;
}

z_number &z_number::operator-=(z_number x) {
  mpz_sub(_n, _n, x._n);
  return *this;
}

z_number &z_number::operator/=(z_number x) {
  if (x == 0) {
    CRAB_ERROR("z_number: division by zero [3]");
  } else {
    mpz_tdiv_q(_n, _n, x._n);
    return *this;
  }
}

z_number &z_number::operator%=(z_number x) {
  if (x == 0) {
    CRAB_ERROR("z_number: division by zero [4]");
  } else {
    mpz_tdiv_r(_n, _n, x._n);
    return *this;
  }
}

z_number &z_number::operator--() {
  mpz_sub_ui(_n, _n, 1);
  return *this;
}

z_number &z_number::operator++() {
  mpz_add_ui(_n, _n, 1);
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

bool z_number::operator==(z_number x) const { return mpz_cmp(_n, x._n) == 0; }

bool z_number::operator!=(z_number x) const { return !(*this == x); }

bool z_number::operator<(z_number x) const { return (mpz_cmp(_n, x._n) < 0); }

bool z_number::operator<=(z_number x) const { return (mpz_cmp(_n, x._n) <= 0); }

bool z_number::operator>(z_number x) const { return (mpz_cmp(_n, x._n) > 0); }

bool z_number::operator>=(z_number x) const { return (mpz_cmp(_n, x._n) >= 0); }

z_number z_number::operator&(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_and(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator|(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_ior(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator^(z_number x) const {
  mpz_t mp_r;
  mpz_init(mp_r);
  mpz_xor(mp_r, _n, x._n);
  z_number res = from_mpz_t(mp_r);
  mpz_clear(mp_r);
  return res;
}

z_number z_number::operator<<(z_number x) const {
  if (x == 0) {
    return x;
  } else {
    mpz_t mp_r;
    mpz_init(mp_r);
    mpz_mul_2exp(mp_r, _n, mpz_get_ui(x._n));
    z_number res = from_mpz_t(mp_r);
    mpz_clear(mp_r);
    return res;
  }
}

z_number z_number::operator>>(z_number x) const {
  if (x == 0) {
    return x;
  } else {
    mpz_t mp_r;
    mpz_init(mp_r);
    mpz_fdiv_q_2exp(mp_r, _n, mpz_get_ui(x._n));
    z_number res = from_mpz_t(mp_r);
    mpz_clear(mp_r);
    return res;
  }
}

z_number z_number::fill_ones() const {
  z_number x(*this);

  assert(x >= 0);
  if (x == 0) {
    return x;
  }

  z_number result(1);
  z_number z1(1);
  z_number z2(2);
  for (; result < x; result = (z2 * result) + z1)
    ;
  return result;
}

void z_number::write(crab::crab_os &o) const { o << get_str(); }

/* Wrapper for mpq */

q_number::q_number() { mpq_init(_n); }

q_number::q_number(double d) {
  mpq_init(_n);
  mpq_set_d(_n, d);
}

q_number::q_number(const std::string &s, unsigned base) {
  mpq_init(_n);
  int res = mpq_set_str(_n, s.c_str(), base);
  if (res != 0) {
    CRAB_ERROR("q_number: invalid string in constructor", s);
  }
}

q_number::q_number(const z_number &z) {
  mpq_init(_n);
  mpq_set_z(_n, z._n);
}

q_number::q_number(const z_number &num, const z_number &den) {
  mpz_init_set(mpq_numref(_n), num._n);
  mpz_init_set(mpq_denref(_n), den._n);
}

q_number q_number::from_mpq_t(mpq_t mp) {
  q_number q;
  mpq_set(q._n, mp);
  return q;
}

q_number q_number::from_mpz_t(mpz_t mp) {
  q_number q;
  mpq_set_z(q._n, mp);
  return q;
}

q_number q_number::from_mpq_srcptr(mpq_srcptr mp) {
  q_number q;
  mpz_init_set(mpq_numref(q._n), mpq_numref(mp));
  mpz_init_set(mpq_denref(q._n), mpq_denref(mp));
  return q;
}

q_number::q_number(const q_number &q) {
  mpz_init_set(mpq_numref(_n), mpq_numref(q._n));
  mpz_init_set(mpq_denref(_n), mpq_denref(q._n));
}

q_number::q_number(q_number &&q) {
  *_n = *q._n;
  mpq_init(q._n);
}

q_number &q_number::operator=(const q_number &q) {
  if (this != &q) {
    mpq_set(_n, q._n);
  }
  return *this;
}

q_number &q_number::operator=(q_number &&q) {
  std::swap(_n, q._n);
  return *this;
}

q_number::~q_number() { mpq_clear(_n); }

double q_number::get_double() const { return mpq_get_d(_n); }

std::string q_number::get_str(unsigned base) const {
  bignums_impl::scoped_cstring res(mpq_get_str(0, base, _n));
  return std::string(res.m_str);
}

std::size_t q_number::hash() const {
  boost::hash<std::string> hasher;
  return hasher(get_str());
}

q_number q_number::operator+(q_number x) const {
  // the price to keep const the method
  mpq_t mp, mp_r;
  mpq_init(mp);
  mpq_set(mp, _n);

  mpq_canonicalize(x._n);
  mpq_canonicalize(mp);

  mpq_init(mp_r);
  mpq_add(mp_r, mp, x._n);
  q_number res = from_mpq_t(mp_r);
  mpq_clear(mp_r);
  return res;
}

q_number q_number::operator*(q_number x) const {
  // the price to keep const the method
  mpq_t mp, mp_r;
  mpq_init(mp);
  mpq_set(mp, _n);

  mpq_canonicalize(x._n);
  mpq_canonicalize(mp);

  mpq_init(mp_r);
  mpq_mul(mp_r, mp, x._n);
  q_number res = from_mpq_t(mp_r);
  mpq_clear(mp_r);
  return res;
}

q_number q_number::operator-(q_number x) const {
  // the price to keep const the method
  mpq_t mp, mp_r;
  mpq_init(mp);
  mpq_set(mp, _n);

  mpq_canonicalize(x._n);
  mpq_canonicalize(mp);

  mpq_init(mp_r);
  mpq_sub(mp_r, mp, x._n);
  q_number res = from_mpq_t(mp_r);
  mpq_clear(mp_r);
  return res;
}

q_number q_number::operator-() const {
  // the price to keep const the method
  mpq_t mp, mp_r;
  mpq_init(mp);
  mpq_set(mp, _n);

  mpq_canonicalize(mp);

  mpq_init(mp_r);
  mpq_neg(mp_r, mp);
  q_number res = from_mpq_t(mp_r);
  mpq_clear(mp_r);
  return res;
}

q_number q_number::operator/(q_number x) const {

  // the price to keep const the method
  mpq_t mp, mp_r;
  mpq_init(mp);
  mpq_set(mp, _n);

  mpq_canonicalize(x._n);
  mpq_canonicalize(mp);

  if (x == 0) {
    CRAB_ERROR("q_number: division by zero [1]");
  } else {
    mpq_init(mp_r);
    mpq_div(mp_r, mp, x._n);
    q_number res = from_mpq_t(mp_r);
    mpq_clear(mp_r);
    return res;
  }
}

q_number &q_number::operator+=(q_number x) {
  mpq_canonicalize(_n);
  mpq_canonicalize(x._n);

  mpq_add(_n, _n, x._n);
  return *this;
}

q_number &q_number::operator*=(q_number x) {
  mpq_canonicalize(_n);
  mpq_canonicalize(x._n);

  mpq_mul(_n, _n, x._n);
  return *this;
}

q_number &q_number::operator-=(q_number x) {
  mpq_canonicalize(_n);
  mpq_canonicalize(x._n);

  mpq_sub(_n, _n, x._n);
  return *this;
}

q_number &q_number::operator/=(q_number x) {
  mpq_canonicalize(_n);
  mpq_canonicalize(x._n);

  if (x == 0) {
    CRAB_ERROR("q_number: division by zero [3]");
  } else {
    mpq_div(_n, _n, x._n);
    return *this;
  }
}

q_number &q_number::operator--() {
  mpq_canonicalize(_n);

  mpz_sub(mpq_numref(_n), mpq_numref(_n), mpq_denref(_n));
  return *this;
}

q_number &q_number::operator++() {
  mpq_canonicalize(_n);

  mpz_add(mpq_numref(_n), mpq_numref(_n), mpq_denref(_n));
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

bool q_number::operator==(q_number x) const { return mpq_cmp(_n, x._n) == 0; }

bool q_number::operator!=(q_number x) const { return !(*this == x); }

bool q_number::operator<(q_number x) const { return mpq_cmp(_n, x._n) < 0; }

bool q_number::operator<=(q_number x) const { return mpq_cmp(_n, x._n) <= 0; }

bool q_number::operator>(q_number x) const { return mpq_cmp(_n, x._n) > 0; }

bool q_number::operator>=(q_number x) const { return mpq_cmp(_n, x._n) >= 0; }

z_number q_number::numerator() const {
  return z_number::from_mpz_srcptr(mpq_numref(_n));
}

z_number q_number::denominator() const {
  return z_number::from_mpz_srcptr(mpq_denref(_n));
}

z_number q_number::round_to_upper() const {
  z_number num = numerator();
  z_number den = denominator();
  z_number q = num / den;
  z_number r = num % den;
  if (r == 0 || *this < 0) {
    return q;
  } else {
    return q + z_number(1);
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
    return q - z_number(1);
  }
}

void q_number::write(crab::crab_os &o) const { o << get_str(); }

} // namespace ikos
