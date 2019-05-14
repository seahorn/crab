#pragma once 

/** 
 *  A class for small, arbitrary-precision unsigned integers.
**/

#include <stdint.h>
#include <sstream>
#include <boost/functional/hash.hpp>
#include <crab/common/types.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {

class wrapint {

public:

  // bitwidth cannot be greater than 64 so even a char can represent
  // all possible bitwidths. However, using uint64_t avoids uintended
  // cast to smaller types that lead to unintended overflows.
  typedef uint64_t bitwidth_t;
  
private:
  
  uint64_t _n;        // 0 <= _n <= 2^_width - 1
  bitwidth_t _width;  // 1 <= _width <= 64
  uint64_t _mod;      // 0 if _width=64 otherwise 2^_width

  static const uint64_t mod_8  = 256;
  static const uint64_t mod_16 = 65536;
  static const uint64_t mod_32 = 4294967296;
  
  void sanity_check_bitwidth() const {
    if (_width == 0) {
      CRAB_ERROR("no bitwidth found for a wrapint");
    }
    if (_width > 64) {
      CRAB_ERROR(_width, " is a too big bitwidth for a wrapint");
    }
  }

  void sanity_check_bitwidths(const wrapint &other) const {
    if (_width != other._width) {
      CRAB_ERROR("two wrapint numbers with different bitwidths");
    }
  }
	   
  void compute_mod() {
    assert (_width <= 64);
    switch (_width){
    case 8:  _mod = mod_8;  break;
    case 16: _mod = mod_16; break;
    case 32: _mod = mod_32; break;
    case 64: break; 
    default: _mod = (uint64_t) 1 << (uint64_t) _width;
    }
  }

  wrapint(uint64_t n, bitwidth_t width, uint64_t mod)
    : _n(n), _width(width), _mod(mod) {}
  
public:

  wrapint(uint64_t n, bitwidth_t width) : _n(n), _width(width), _mod(0) {
    sanity_check_bitwidth();
    assert(_width <= 64);
    compute_mod();

    if (_width < 64) {
      _n = _n % _mod;
    }
  }

  wrapint(ikos::z_number n, bitwidth_t width)
    : _n(0), _width(width), _mod(0) {
    sanity_check_bitwidth();
    assert(_width <= 64);
    compute_mod();
    if (!n.fits_slong()) {
      // n is a signed integer over 64 bits that cannot be represented
      // by the signed long type.
      CRAB_ERROR(n, " does not fit in a signed long integer.");
    }
    if (_width == 64) {
      _n = (long) n;
    } else {
      _n = ((long) n) % _mod;
    }
  }

  wrapint(ikos::q_number n, bitwidth_t width) : _n(0), _width(width), _mod(0) {
    sanity_check_bitwidth();
    assert(_width <= 64);
    compute_mod();
    ikos::z_number i = n.round_to_upper();
    if (!i.fits_slong()) {
      // n is a signed integer over 64 bits that cannot be represented
      // by the signed long type.      
      CRAB_ERROR(n, " does not fit in a signed long integer.");
    }
    if (_width == 64) {
      _n = (long) i;
    } else {
      _n = ((long) i) % _mod;
    }
    
  }
  
  wrapint(std::string s, bitwidth_t width): _n(0), _width(width), _mod(0) {
    sanity_check_bitwidth();
    assert(_width <= 64);
    compute_mod();
    std::istringstream iss(s);
    iss >> _n;

    if (_width < 64) {
      _n = _n % _mod;
    }
  }

  bitwidth_t get_bitwidth() const { return _width;}

  // Needed because wrapint has limited precision
  static bool fits_wrapint(ikos::z_number n, bitwidth_t width) {
    if (width > 64) return false;
    return n.fits_slong();
  }

  // Needed because wrapint has limited precision
  static bool fits_wrapint(ikos::q_number n, bitwidth_t width) {
    return fits_wrapint(n.round_to_upper(), width);
  }

  // return true iff most significant bit is 1.
  bool msb() const {
    return (_n & ((uint64_t)1 << (uint64_t)(_width - 1)));
  }
  
  // return 01111...1
  static wrapint get_signed_max(bitwidth_t w)  {
    return wrapint(((uint64_t)1 << (uint64_t)(w-1)) - 1, w);
  }

  // return 1000....0
  static wrapint get_signed_min(bitwidth_t w) {
    return wrapint((uint64_t)1 << (uint64_t)(w-1), w);
  }

  // return 1111....1
  static wrapint get_unsigned_max(bitwidth_t w) {
    switch (w) {
    case 8:  return wrapint(mod_8 - 1, w);
    case 16: return wrapint(mod_16 - 1, w);
    case 32: return wrapint(mod_32 - 1, w);
    case 64: return wrapint(UINT64_MAX, w);
    default: return wrapint(((uint64_t) 1 << (uint64_t)w) - 1, w);
    }
  }

  // return 0000....0
  static wrapint get_unsigned_min(bitwidth_t w) {
    return wrapint(0, w);
  }
  
public:

  // return the wrapint as an unsigned number
  std::string get_unsigned_str () const {
    return get_unsigned_bignum().get_str();
  }

  // return the wrapint as a signed number
  std::string get_signed_str () const {
    return get_signed_bignum().get_str();
  }
  
  uint64_t get_uint64_t() const {
    return _n;
  }

  // return the wrapint as an unsigned big number
  ikos::z_number get_unsigned_bignum() const {
    // XXX: cannot use here ikos::z_number(_n) because it will cast _n
    // implicitly to a signed number. Thus, max unsigned will become
    // signed -1.
    //
    // FIXME: In some platforms uint64_t may not be the same as
    // unsigned long.
    return ikos::z_number::from_ulong(_n);
  }

  // return the wrapint as a signed big number
  ikos::z_number get_signed_bignum() const {
    if (msb()) {
      // get two's complement and negate
      wrapint r = *this ^ get_unsigned_max(get_bitwidth());
      return ikos::z_number(-(r.get_unsigned_bignum()+1));
    } else {
      return get_unsigned_bignum();
    }
  }
  
  bool is_zero() const {
    return _n == 0;
  }
  
  wrapint operator+(wrapint x) const {
    sanity_check_bitwidths(x);
    
    uint64_t r = (_width == 64 ? (_n + x._n) : (_n + x._n) % _mod);
    return wrapint(r, _width, _mod);
  }

  wrapint operator*(wrapint x) const {
    sanity_check_bitwidths(x);
    
    uint64_t r = (_width == 64 ? (_n * x._n) : (_n * x._n) % _mod);
    return wrapint(r, _width, _mod);
  }

  wrapint operator-(wrapint x) const {
    sanity_check_bitwidths(x);
    
    uint64_t r = (_width == 64 ? (_n - x._n) : (_n - x._n) % _mod);
    return wrapint(r, _width, _mod);
  }

  wrapint operator-() const {
    uint64_t r = (_width == 64 ? -_n : -_n % _mod);
    return wrapint(r, _width, _mod);
  }

  // signed division
  wrapint operator/(wrapint x) const {
    return sdiv(x);
  }

  // signed division: rounding towards 0
  wrapint sdiv(wrapint x) const {
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
  wrapint udiv(wrapint x) const {
    sanity_check_bitwidths(x);    
    if (x.is_zero()) {
      CRAB_ERROR("wrapint: unsigned division by zero ", __LINE__);
    } else {
      uint64_t r = (_width == 64 ? _n / x._n : (_n / x._n) % _mod);
      return wrapint(r, _width, _mod);
    }
  }

  // signed remainder
  wrapint operator%(wrapint x) const {
    return srem(x);
  }

  // signed rem: is the remainder of the signed division so rounding
  // also towards 0.
  wrapint srem(wrapint x) const {
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
  wrapint urem(wrapint x) const {
    sanity_check_bitwidths(x);
    if (x.is_zero()) {
      CRAB_ERROR("wrapint: unsigned division by zero ", __LINE__);
    } else {
      uint64_t r = (_width == 64 ? _n % x._n: (_n % x._n) % _mod);
      return wrapint(r, _width, _mod);
    }
  }
  
  wrapint& operator+=(wrapint x) {
    sanity_check_bitwidths(x);
    _n += x._n;
    if (_width < 64) _n = _n % _mod;
    return *this;
  }

  wrapint& operator*=(wrapint x) {
    sanity_check_bitwidths(x);
    _n *= x._n;
    if (_width < 64) _n = _n % _mod;
    return *this;
  }

  wrapint& operator-=(wrapint x) {
    sanity_check_bitwidths(x);
    _n -= x._n;
    if (_width < 64) _n = _n % _mod;
    return *this;
  }


  wrapint& operator--() {
    --(_n);
    if (_width < 64) _n = _n % _mod;
    return *this;
  }

  wrapint& operator++() {
    ++(_n);
    if (_width < 64) _n = _n % _mod;
    return *this;
  }

  wrapint operator++(int) {
    wrapint r(*this);
    ++(*this);
    return r;
  }

  wrapint operator--(int) {
    wrapint r(*this);
    --(*this);
    return r;
  }

  bool operator==(wrapint x) const {
    sanity_check_bitwidths(x);    
    return _n == x._n;
  }

  bool operator!=(wrapint x) const {
    sanity_check_bitwidths(x);        
    return _n != x._n;
  }

  bool operator<(wrapint x) const {
    sanity_check_bitwidths(x);            
    return _n < x._n;
  }

  bool operator<=(wrapint x) const {
    sanity_check_bitwidths(x);                
    return _n <= x._n;
  }

  bool operator>(wrapint x) const {
    sanity_check_bitwidths(x);                    
    return _n > x._n;
  }

  bool operator>=(wrapint x) const {
    sanity_check_bitwidths(x);                        
    return _n >= x._n;
  }

  wrapint operator&(wrapint x) const {
    sanity_check_bitwidths(x);    
    return wrapint(_n & x._n, _width, _mod);
  }

  wrapint operator|(wrapint x) const {
    sanity_check_bitwidths(x);        
    return wrapint(_n | x._n, _width, _mod);
  }

  wrapint operator^(wrapint x) const {
    sanity_check_bitwidths(x);            
    return wrapint(_n ^ x._n, _width, _mod);
  }

  wrapint operator<<(wrapint x) const {
    sanity_check_bitwidths(x);
    
    uint64_t r = (_width == 64 ? (_n << x._n) : (_n << x._n) % _mod);
    return wrapint(r, _width, _mod);
  }

  // logical right shift: blanks filled by 0's
  wrapint lshr(wrapint x) const {
    sanity_check_bitwidths(x);
    return wrapint(_n >> x._n, _width, _mod);
  }

  // arithmetic right shift
  wrapint ashr(wrapint x) const {
    sanity_check_bitwidths(x);
    if (!msb()) {
      return wrapint(_n >> x._n, _width, _mod);      
    } else {
      // fill blanks with 1's
      uint64_t all_ones = (_width < 64
			   ? ((uint64_t)1 << (uint64_t)_width) - 1
			   : UINT64_MAX);
      // 1110..0
      uint64_t only_upper_bits_ones = all_ones << (uint64_t)(_width - x._n);
      return wrapint(only_upper_bits_ones | (_n >> x._n), _width, _mod);
    }
  }
  
  wrapint sext(bitwidth_t bits_to_add) {
    bitwidth_t new_width = _width + bits_to_add;
    if (new_width > 64) {
      CRAB_ERROR("cannot signed extend: ", new_width, " is a too big bitwidth for a wrapint");
    }

    if (msb()) {
      // -- fill upper bits with ones
      // 111...1
      uint64_t all_ones = (new_width < 64
			   ? ((uint64_t) 1 << (uint64_t) new_width) - 1
			   : UINT64_MAX);
      // 1110..0
      uint64_t only_upper_bits_ones = all_ones << (uint64_t) _width;
      return wrapint(_n | only_upper_bits_ones, new_width);
    } else {
      // -- fill upper bits with zeros
      return wrapint(_n, _width + bits_to_add);       
    }
  }

  wrapint zext(bitwidth_t bits_to_add) const {
    bitwidth_t new_width = _width + bits_to_add;
    if (new_width > 64) {
      CRAB_ERROR("cannot unsigned extend: ", new_width, " is a too big bitwidth for a wrapint");
    }
    return wrapint(_n, new_width); 
  }

  wrapint keep_lower(bitwidth_t bits_to_keep) const {
    if (bits_to_keep >= _width) return *this;
    return wrapint(_n & (((uint64_t) 1 << (uint64_t) (bits_to_keep+1)) -1),
		   bits_to_keep);
  }
  
  void write(crab::crab_os& o) const {
    o << get_unsigned_str();
  }

}; // class wrapint


inline crab::crab_os& operator<<(crab::crab_os& o, const wrapint& z) {
  z.write(o);
  return o;
}

inline std::size_t hash_value(const wrapint& n) {
  boost::hash<std::string> hasher;
  return hasher(n.get_unsigned_str());
}

} // end namespace crab
