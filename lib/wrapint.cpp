#include <crab/common/types.hpp>
#include <crab/numbers/wrapint.hpp>
//#include <boost/functional/hash.hpp>

#include <climits>
#include <sstream>
#include <string>

namespace crab {

void wrapint::sanity_check_bitwidth() const {
  if (_width == 0) {
    CRAB_ERROR("no bitwidth found for a wrapint");
  }
  if (_width > 64) {
    CRAB_ERROR(_width, " is a too big bitwidth for a wrapint");
  }
}

void wrapint::sanity_check_bitwidths(const wrapint &other) const {
  if (_width != other._width) {
    CRAB_ERROR("two wrapint numbers with different bitwidths");
  }
}

void wrapint::compute_mod() {
  assert(_width <= 64);
  switch (_width) {
  case 8:
    _mod = mod_8;
    break;
  case 16:
    _mod = mod_16;
    break;
  case 32:
    _mod = mod_32;
    break;
  case 64:
    break;
  default:
    _mod = (uint64_t)1 << (uint64_t)_width;
  }
}

wrapint::wrapint(uint64_t n, bitwidth_t width, uint64_t mod)
    : _n(n), _width(width), _mod(mod) {}

wrapint::wrapint(uint64_t n, bitwidth_t width) : _n(n), _width(width), _mod(0) {

  sanity_check_bitwidth();
  assert(_width <= 64);
  compute_mod();

  if (_width < 64) {
    _n = _n % _mod;
  }
}

wrapint::wrapint(ikos::z_number n, bitwidth_t width)
    : _n(0), _width(width), _mod(0) {

  sanity_check_bitwidth();
  assert(_width <= 64);
  compute_mod();
  if (!n.fits_int64()) {
    CRAB_ERROR(n, " does not fit in an int64_t.");
  }
  int64_t x = static_cast<int64_t>(n);
  if (_width == 64) {
    _n = static_cast<uint64_t>(x);
  } else {
    _n = static_cast<uint64_t>(x) % _mod;
  }
}

wrapint::wrapint(ikos::q_number n, bitwidth_t width)
    : _n(0), _width(width), _mod(0) {

  sanity_check_bitwidth();
  assert(_width <= 64);
  compute_mod();
  ikos::z_number i = n.round_to_upper();
  if (!i.fits_int64()) {
    CRAB_ERROR(n, " does not fit in an int64_t.");
  }
  int64_t x = static_cast<int64_t>(i);
  if (_width == 64) {
    _n = static_cast<uint64_t>(x);
  } else {
    _n = static_cast<uint64_t>(x) % _mod;
  }
}

wrapint::wrapint(std::string s, bitwidth_t width)
    : _n(0), _width(width), _mod(0) {

  sanity_check_bitwidth();
  assert(_width <= 64);
  compute_mod();
  std::istringstream iss(s);
  iss >> _n;

  if (_width < 64) {
    _n = _n % _mod;
  }
}

wrapint::bitwidth_t wrapint::get_bitwidth() const { return _width; }

bool wrapint::fits_wrapint(ikos::z_number n, bitwidth_t width) {
  if (width > 64)
    return false;
  return n.fits_int64();
}

bool wrapint::fits_wrapint(ikos::q_number n, bitwidth_t width) {
  return fits_wrapint(n.round_to_upper(), width);
}

// return true iff most significant bit is 1.
bool wrapint::msb() const {
  return (_n & ((uint64_t)1 << (uint64_t)(_width - 1)));
}

// return 01111...1
wrapint wrapint::get_signed_max(bitwidth_t w) {
  return wrapint(((uint64_t)1 << (uint64_t)(w - 1)) - 1, w);
}

// return 1000....0
wrapint wrapint::get_signed_min(bitwidth_t w) {
  return wrapint((uint64_t)1 << (uint64_t)(w - 1), w);
}

// return 1111....1
wrapint wrapint::get_unsigned_max(bitwidth_t w) {
  switch (w) {
  case 8:
    return wrapint(mod_8 - 1, w);
  case 16:
    return wrapint(mod_16 - 1, w);
  case 32:
    return wrapint(mod_32 - 1, w);
  case 64:
    return wrapint(UINT64_MAX, w);
  default:
    return wrapint(((uint64_t)1 << (uint64_t)w) - 1, w);
  }
}

// return 0000....0
wrapint wrapint::get_unsigned_min(bitwidth_t w) { return wrapint(0, w); }

// return the wrapint as an unsigned number
std::string wrapint::get_unsigned_str() const {
  return get_unsigned_bignum().get_str();
}

// return the wrapint as a signed number
std::string wrapint::get_signed_str() const {
  return get_signed_bignum().get_str();
}

uint64_t wrapint::get_uint64_t() const { return _n; }

// return the wrapint as an unsigned big number
ikos::z_number wrapint::get_unsigned_bignum() const {
  // XXX: cannot use here ikos::z_number(_n) because it will cast _n
  // implicitly to a signed integer.
  return ikos::z_number::from_uint64(_n);
}

// return the wrapint as a signed big number
ikos::z_number wrapint::get_signed_bignum() const {
  if (msb()) {
    // get two's complement and negate
    wrapint r = *this ^ get_unsigned_max(get_bitwidth());
    return ikos::z_number(-(r.get_unsigned_bignum() + 1));
  } else {
    return get_unsigned_bignum();
  }
}

bool wrapint::is_zero() const { return _n == 0; }

wrapint wrapint::operator+(wrapint x) const {
  sanity_check_bitwidths(x);

  uint64_t r = (_width == 64 ? (_n + x._n) : (_n + x._n) % _mod);
  return wrapint(r, _width, _mod);
}

wrapint wrapint::operator*(wrapint x) const {
  sanity_check_bitwidths(x);

  uint64_t r = (_width == 64 ? (_n * x._n) : (_n * x._n) % _mod);
  return wrapint(r, _width, _mod);
}

wrapint wrapint::operator-(wrapint x) const {
  sanity_check_bitwidths(x);

  uint64_t r = (_width == 64 ? (_n - x._n) : (_n - x._n) % _mod);
  return wrapint(r, _width, _mod);
}

wrapint wrapint::operator-() const {
  uint64_t r = (_width == 64 ? -_n : -_n % _mod);
  return wrapint(r, _width, _mod);
}

// signed division
wrapint wrapint::operator/(wrapint x) const { return sdiv(x); }

// signed division: rounding towards 0
wrapint wrapint::sdiv(wrapint x) const {
  sanity_check_bitwidths(x);
  if (x.is_zero()) {
    CRAB_ERROR("wrapint: signed division by zero ", __LINE__);
  } else {
    ikos::z_number dividend = get_signed_bignum();
    ikos::z_number divisor = x.get_signed_bignum();
    ikos::z_number r = dividend / divisor;
    return wrapint(r, get_bitwidth());
  }
}

// unsigned division: rounding towards 0
wrapint wrapint::udiv(wrapint x) const {
  sanity_check_bitwidths(x);
  if (x.is_zero()) {
    CRAB_ERROR("wrapint: unsigned division by zero ", __LINE__);
  } else {
    uint64_t r = (_width == 64 ? _n / x._n : (_n / x._n) % _mod);
    return wrapint(r, _width, _mod);
  }
}

// signed remainder
wrapint wrapint::operator%(wrapint x) const { return srem(x); }

// signed rem: is the remainder of the signed division so rounding
// also towards 0.
wrapint wrapint::srem(wrapint x) const {
  sanity_check_bitwidths(x);
  if (x.is_zero()) {
    CRAB_ERROR("wrapint: signed division by zero ", __LINE__);
  } else {
    ikos::z_number dividend = get_signed_bignum();
    ikos::z_number divisor = x.get_signed_bignum();
    ikos::z_number r = dividend % divisor;
    return wrapint(r, get_bitwidth());
  }
}

// unsigned rem: is the remainder of unsigned division so rounding
// also towards 0.
wrapint wrapint::urem(wrapint x) const {
  sanity_check_bitwidths(x);
  if (x.is_zero()) {
    CRAB_ERROR("wrapint: unsigned division by zero ", __LINE__);
  } else {
    uint64_t r = (_width == 64 ? _n % x._n : (_n % x._n) % _mod);
    return wrapint(r, _width, _mod);
  }
}

wrapint &wrapint::operator+=(wrapint x) {
  sanity_check_bitwidths(x);
  _n += x._n;
  if (_width < 64)
    _n = _n % _mod;
  return *this;
}

wrapint &wrapint::operator*=(wrapint x) {
  sanity_check_bitwidths(x);
  _n *= x._n;
  if (_width < 64)
    _n = _n % _mod;
  return *this;
}

wrapint &wrapint::operator-=(wrapint x) {
  sanity_check_bitwidths(x);
  _n -= x._n;
  if (_width < 64)
    _n = _n % _mod;
  return *this;
}

wrapint &wrapint::operator--() {
  --(_n);
  if (_width < 64)
    _n = _n % _mod;
  return *this;
}

wrapint &wrapint::operator++() {
  ++(_n);
  if (_width < 64)
    _n = _n % _mod;
  return *this;
}

wrapint wrapint::operator++(int) {
  wrapint r(*this);
  ++(*this);
  return r;
}

wrapint wrapint::operator--(int) {
  wrapint r(*this);
  --(*this);
  return r;
}

bool wrapint::operator==(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n == x._n;
}

bool wrapint::operator!=(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n != x._n;
}

bool wrapint::operator<(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n < x._n;
}

bool wrapint::operator<=(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n <= x._n;
}

bool wrapint::operator>(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n > x._n;
}

bool wrapint::operator>=(wrapint x) const {
  sanity_check_bitwidths(x);
  return _n >= x._n;
}

wrapint wrapint::operator&(wrapint x) const {
  sanity_check_bitwidths(x);
  return wrapint(_n & x._n, _width, _mod);
}

wrapint wrapint::operator|(wrapint x) const {
  sanity_check_bitwidths(x);
  return wrapint(_n | x._n, _width, _mod);
}

wrapint wrapint::operator^(wrapint x) const {
  sanity_check_bitwidths(x);
  return wrapint(_n ^ x._n, _width, _mod);
}

wrapint wrapint::operator<<(wrapint x) const {
  sanity_check_bitwidths(x);

  uint64_t r = (_width == 64 ? (_n << x._n) : (_n << x._n) % _mod);
  return wrapint(r, _width, _mod);
}

// logical right shift: blanks filled by 0's
wrapint wrapint::lshr(wrapint x) const {
  sanity_check_bitwidths(x);
  return wrapint(_n >> x._n, _width, _mod);
}

// arithmetic right shift
wrapint wrapint::ashr(wrapint x) const {
  sanity_check_bitwidths(x);
  if (!msb()) {
    return wrapint(_n >> x._n, _width, _mod);
  } else {
    // fill blanks with 1's
    uint64_t all_ones =
        (_width < 64 ? ((uint64_t)1 << (uint64_t)_width) - 1 : UINT64_MAX);
    // 1110..0
    uint64_t only_upper_bits_ones = all_ones << (uint64_t)(_width - x._n);
    return wrapint(only_upper_bits_ones | (_n >> x._n), _width, _mod);
  }
}

wrapint wrapint::sext(bitwidth_t bits_to_add) const {
  bitwidth_t new_width = _width + bits_to_add;
  if (new_width > 64) {
    CRAB_ERROR("cannot signed extend: ", new_width,
               " is a too big bitwidth for a wrapint");
  }

  if (msb()) {
    // -- fill upper bits with ones
    // 111...1
    uint64_t all_ones =
        (new_width < 64 ? ((uint64_t)1 << (uint64_t)new_width) - 1
                        : UINT64_MAX);
    // 1110..0
    uint64_t only_upper_bits_ones = all_ones << (uint64_t)_width;
    return wrapint(_n | only_upper_bits_ones, new_width);
  } else {
    // -- fill upper bits with zeros
    return wrapint(_n, _width + bits_to_add);
  }
}

wrapint wrapint::zext(bitwidth_t bits_to_add) const {
  bitwidth_t new_width = _width + bits_to_add;
  if (new_width > 64) {
    CRAB_ERROR("cannot unsigned extend: ", new_width,
               " is a too big bitwidth for a wrapint");
  }
  return wrapint(_n, new_width);
}

wrapint wrapint::keep_lower(bitwidth_t bits_to_keep) const {
  if (bits_to_keep >= _width)
    return *this;
  return wrapint(_n & (((uint64_t)1 << (uint64_t)(bits_to_keep + 1)) - 1),
                 bits_to_keep);
}

void wrapint::write(crab::crab_os &o) const { o << get_unsigned_str(); }

// inline std::size_t hash_value(const wrapint& n) {
//   boost::hash<std::string> hasher;
//   return hasher(n.get_unsigned_str());
// }

} // end namespace crab
