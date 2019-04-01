#pragma once 

/** 
 *  Safe signed integers.
**/

#include <limits>

namespace crab {

class safe_i64 {
  // Current implementation is based on
  // https://blog.regehr.org/archives/1139 using wider integers.

  // Both clang and gcc supports __int128 if the targeted architecture
  // is x86/64, but it wont' work with 32 bits.
  typedef __int128 wideint_t;
  
  inline int64_t MAX() const {
    return std::numeric_limits<int64_t>::max();
  }
  
  inline int64_t MIN() const {
    return std::numeric_limits<int64_t>::min();
  }
  
  int checked_add(int64_t a, int64_t b, int64_t *rp) const {
    #if 1
    wideint_t lr = (wideint_t)a + (wideint_t)b;
    *rp = lr;
    return lr > MAX() || lr < MIN();
    #else
    // without wider integers
    if (b > 0 && a > MAX() - b) {
      return 1;
    }
    if (b < 0 && a < MIN() - b) {
      return 1;
    }
    int64_t lr = a + b;
    *rp = lr;
    return 0;
    #endif 
  }
  
  int checked_sub(int64_t a, int64_t b, int64_t *rp) const {
    wideint_t lr = (wideint_t)a - (wideint_t)b;
    *rp = lr;
    return lr > MAX() || lr < MIN();
  }
  
  int checked_mul(int64_t a, int64_t b, int64_t *rp) const {
    wideint_t lr = (wideint_t)a * (wideint_t)b;
    *rp = lr;
    return lr > MAX() || lr < MIN();
  }
  
  int checked_div(int64_t a, int64_t b, int64_t *rp) const {
    wideint_t lr = (wideint_t)a / (wideint_t)b;
    *rp = lr;
    return lr > MAX() || lr < MIN();
  }
  
public:
  
  safe_i64(): m_num(0) {}
  
  safe_i64(int64_t num): m_num(num) {}
  
  safe_i64(ikos::z_number n): m_num((long) n) {}
  
  operator long() const {
    return (long) m_num;
  }

  // FIXME: operation should not raise an error.
  safe_i64 operator+(safe_i64 x) const {
    int64_t z;
    int err = checked_add(m_num, x.m_num, &z);
    if (err) {
      CRAB_ERROR("Integer overflow during addition");
    }
    return safe_i64(z);
  }

  // FIXME: operation should not raise an error.  
  safe_i64 operator-(safe_i64 x) const {
    int64_t z;
    int err = checked_sub(m_num, x.m_num, &z);
    if (err) {
      CRAB_ERROR("Integer overflow during subtraction");
    }
    return safe_i64(z);
  }

  // FIXME: operation should not raise an error.  
  safe_i64 operator*(safe_i64 x) const {
    int64_t z;
    int err = checked_mul(m_num, x.m_num, &z);
    if (err) {
      CRAB_ERROR("Integer overflow during multiplication");
    }
    return safe_i64(z);
  }

  // FIXME: operation should not raise an error.  
  safe_i64 operator/(safe_i64 x) const {
    int64_t z;
    int err = checked_div(m_num, x.m_num, &z);
    if (err) {
      CRAB_ERROR("Integer overflow during multiplication");
    }
    return safe_i64(z);
  }

  // FIXME: operation should not raise an error.  
  safe_i64 operator-() const {
    return safe_i64(0) - *this;
  }

  // FIXME: operation should not raise an error.  
  safe_i64& operator+=(safe_i64 x) {
    int err = checked_add(m_num, x.m_num, &m_num);
    if (err) {
      CRAB_ERROR("Integer overflow during addition");
    }
    return *this;
  }

  // FIXME: operation should not raise an error.  
  safe_i64& operator-=(safe_i64 x) {
    int err = checked_sub(m_num, x.m_num, &m_num);
    if (err) {
      CRAB_ERROR("Integer overflow during subtraction");
    }
    return *this;
  }
  
  bool operator==(safe_i64 x) const {
    return m_num == x.m_num;
  }
  
  bool operator!=(safe_i64 x) const {
    return m_num != x.m_num;
  }
  
  bool operator<(safe_i64 x) const {
    return m_num < x.m_num;
  }
  
  bool operator<=(safe_i64 x) const {
    return m_num <= x.m_num;
  }

  bool operator>(safe_i64 x) const {
    return m_num > x.m_num;
  }
  
  bool operator>=(safe_i64 x) const {
    return m_num >= x.m_num;
  }
  
  void write(crab::crab_os& os) const {
    os << m_num;
  }
  
  friend crab_os& operator<<(crab_os& o, const safe_i64& n) {
    n.write(o);
    return o;
  }
  
private:
  
  int64_t m_num;
};

} // end namespace crab
