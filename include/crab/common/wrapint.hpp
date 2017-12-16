#pragma once 

/** 
 *  A class for small, arbitrary-precision unsigned integers.
**/

#include <stdint.h>
#include <sstream>
#include <boost/functional/hash.hpp>
#include <crab/common/types.hpp>
#include <crab/common/bignums.hpp>

namespace crab {

class wrapint {

public:
  
  typedef unsigned bitwidth_t;
  
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
    default: _mod = 1 << _width;
    }
  }

  wrapint(uint64_t n, bitwidth_t width, uint64_t mod)
    : _n(n), _width(width), _mod(mod) {}
  
public:

  wrapint(uint64_t n, bitwidth_t width) : _n(n), _width(width), _mod(0) {
    sanity_check_bitwidth();
    compute_mod();
    if (_width < 64) _n = _n % _mod;
  }

  wrapint(ikos::z_number n, bitwidth_t width) : _n((long) n), _width(width), _mod(0) {
    sanity_check_bitwidth();
    compute_mod();
    if (_width < 64) _n = _n % _mod;
  }

  wrapint(ikos::q_number n, bitwidth_t width)
    : _n((long) n.round_to_upper()), _width(width), _mod(0) {
    sanity_check_bitwidth();
    compute_mod();
    if (_width < 64) _n = _n % _mod;
  }
  
  wrapint(std::string s, bitwidth_t width): _n(0), _width(width), _mod(0) {
    sanity_check_bitwidth();
    std::istringstream iss(s);
    iss >> _n;
    compute_mod();
    if (_width < 64) _n = _n % _mod;
  }

  bitwidth_t get_bitwidth() const { return _width;}

  // return 01111...1
  static wrapint get_signed_max(bitwidth_t w)  {
    return wrapint((1 << (w-1)) - 1, w);
  }

  // return 1000....0
  static wrapint get_signed_min(bitwidth_t w) {
    return wrapint(1 << (w-1), w);
  }

  // return 1111....1
  static wrapint get_unsigned_max(bitwidth_t w) {
    switch (w) {
    case 8:  return wrapint(mod_8 - 1, w);
    case 16: return wrapint(mod_16 - 1, w);
    case 32: return wrapint(mod_32 - 1, w);
    case 64: return wrapint(UINT64_MAX, w);
    default: return wrapint((1 << w) - 1, w);
    }
  }

  // return 0000....0
  static wrapint get_unsigned_min(bitwidth_t w) {
    return wrapint(0, w);
  }
  
public:

  std::string get_str () const {
    return std::to_string(_n);
  }

  uint64_t get_uint64_t() const {
    return _n;
  }
  
  ikos::z_number get_bignum() const {
    return ikos::z_number(_n);
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

  wrapint operator/(wrapint x) const {
    sanity_check_bitwidths(x);    
    if (x._n == 0) {
      CRAB_ERROR("wrapint: division by zero ", __LINE__);
    } else {
      
      uint64_t r = (_width == 64 ? _n / x._n : (_n / x._n) % _mod);
      return wrapint(r, _width, _mod);
    }
  }

  wrapint operator%(wrapint x) const {
    sanity_check_bitwidths(x);
    if (x._n == 0) {
      CRAB_ERROR("wrapint: division by zero ", __LINE__);
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

  wrapint& operator/=(wrapint x) {
    sanity_check_bitwidths(x);
    if (x._n == 0) {
      CRAB_ERROR("wraping: division by zero ", __LINE__);
    } else {
      _n /= x._n;
      if (_width < 64) _n = _n % _mod;
      return *this;
    }
  }

  wrapint& operator%=(wrapint x) {
    sanity_check_bitwidths(x);
    if (x._n == 0) {
      CRAB_ERROR("wrapint: division by zero ", __LINE__);
    } else {
      _n %= x._n;
      if (_width < 64) _n = _n % _mod;
      return *this;
    }
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

  wrapint operator>>(wrapint x) const {
    sanity_check_bitwidths(x);
    return wrapint(_n >> x._n, _width, _mod);
  }

  void write(crab::crab_os& o) {
    o << get_str();    
  }

}; // class wrapint


inline crab::crab_os& operator<<(crab::crab_os& o, wrapint z) {
  z.write(o);
  return o;
}

inline std::size_t hash_value(const wrapint& n) {
  boost::hash<std::string> hasher;
  return hasher(n.get_str());
}

} // end namespace crab
