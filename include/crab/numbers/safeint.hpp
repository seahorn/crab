#pragma once

/**
 *  Safe signed integers.
 **/

#include <crab/numbers/bignums.hpp>
#include <limits>
#include <cstdint>

#ifndef __GNUC__
#include <boost/multiprecision/cpp_int.hpp>
#endif

namespace crab {

class safe_i64 {

  // Current implementation is based on
  // https://blog.regehr.org/archives/1139 using wider integers.

#ifdef __GNUC__  
  // TODO/FIXME: the current code compiles assuming the type __int128
  // exists. Both clang and gcc supports __int128 if the targeted
  // architecture is x86/64, but it wont' work with 32 bits.
  using wideint_t = __int128;
#else
    using wideint_t = boost::multiprecision::int128_t;
#endif

  static int64_t get_max();
  static int64_t get_min();

  static int checked_add(int64_t a, int64_t b, int64_t *rp);
  static int checked_sub(int64_t a, int64_t b, int64_t *rp);
  static int checked_mul(int64_t a, int64_t b, int64_t *rp);
  static int checked_div(int64_t a, int64_t b, int64_t *rp);

public:
  safe_i64();

  safe_i64(int64_t num);

  safe_i64(ikos::z_number n);

  operator int64_t() const;

  // TODO: output parameters whether operation overflows
  safe_i64 operator+(safe_i64 x) const;

  // TODO: output parameters whether operation overflows
  safe_i64 operator-(safe_i64 x) const;

  // TODO: output parameters whether operation overflows
  safe_i64 operator*(safe_i64 x) const;

  // TODO: output parameters whether operation overflows
  safe_i64 operator/(safe_i64 x) const;

  // TODO: output parameters whether operation overflows
  safe_i64 operator-() const;

  // TODO: output parameters whether operation overflows
  safe_i64 &operator+=(safe_i64 x);

  // TODO: output parameters whether operation overflows
  safe_i64 &operator-=(safe_i64 x);

  std::size_t hash() const;
  
  bool operator==(safe_i64 x) const;

  bool operator!=(safe_i64 x) const;

  bool operator<(safe_i64 x) const;

  bool operator<=(safe_i64 x) const;

  bool operator>(safe_i64 x) const;

  bool operator>=(safe_i64 x) const;

  void write(crab::crab_os &os) const;

  friend crab_os &operator<<(crab_os &o, const safe_i64 &n) {
    n.write(o);
    return o;
  }

private:
  int64_t m_num;
};

} // end namespace crab
