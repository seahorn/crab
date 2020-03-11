#pragma once

#include <boost/functional/hash.hpp>
#include <crab/common/os.hpp>
#include <cstdint>
#include <gmp.h>

// TODO: replace ikos with crab namespace. This class has nothing to
// do with the ikos one. Kept for now for compatibility issues with
// some clients.
namespace ikos {

class z_number {
  friend class q_number;

private:
  mpz_t _n;

  bool fits_sint() const;
  bool fits_slong() const;

public:
  // overloaded typecast operators
  explicit operator int64_t() const;

  z_number();
  z_number(int64_t n);
  z_number(const std::string &s, unsigned base = 10);

  static z_number from_uint64(uint64_t n);
  static z_number from_mpz_t(mpz_t n);
  static z_number from_mpz_srcptr(mpz_srcptr n);

  z_number(const z_number &o);
  z_number(z_number &&o);
  z_number &operator=(const z_number &o);
  z_number &operator=(z_number &&o);

  ~z_number();

  mpz_srcptr get_mpz_t() const { return _n; }

  mpz_ptr get_mpz_t() { return _n; }

  std::string get_str(unsigned base = 10) const;

  std::size_t hash() const;

  bool fits_int64() const;

  z_number operator+(z_number x) const;

  z_number operator*(z_number x) const;

  z_number operator-(z_number x) const;

  z_number operator-() const;

  z_number operator/(z_number x) const;

  z_number operator%(z_number x) const;

  z_number &operator+=(z_number x);

  z_number &operator*=(z_number x);

  z_number &operator-=(z_number x);

  z_number &operator/=(z_number x);

  z_number &operator%=(z_number x);

  z_number &operator--();

  z_number &operator++();

  z_number operator++(int);

  z_number operator--(int);

  bool operator==(z_number x) const;

  bool operator!=(z_number x) const;

  bool operator<(z_number x) const;

  bool operator<=(z_number x) const;

  bool operator>(z_number x) const;

  bool operator>=(z_number x) const;

  z_number operator&(z_number x) const;

  z_number operator|(z_number x) const;

  z_number operator^(z_number x) const;

  z_number operator<<(z_number x) const;

  z_number operator>>(z_number x) const;

  z_number fill_ones() const;

  void write(crab::crab_os &o) const;

}; // class z_number

class q_number {

private:
  mpq_t _n;

public:
  q_number();
  q_number(double n);

  q_number(const std::string &s, unsigned base = 10);
  q_number(const z_number &n);
  q_number(const z_number &n, const z_number &d);

  static q_number from_mpq_t(mpq_t n);
  static q_number from_mpz_t(mpz_t n);
  static q_number from_mpq_srcptr(mpq_srcptr q);

  q_number(const q_number &o);
  q_number(q_number &&o);
  q_number &operator=(const q_number &o);
  q_number &operator=(q_number &&o);

  ~q_number();

  mpq_srcptr get_mpq_t() const { return _n; }

  mpq_ptr get_mpq_t() { return _n; }

  double get_double() const;

  std::string get_str(unsigned base = 10) const;

  std::size_t hash() const;

  q_number operator+(q_number x) const;

  q_number operator*(q_number x) const;

  q_number operator-(q_number x) const;

  q_number operator-() const;

  q_number operator/(q_number x) const;

  q_number &operator+=(q_number x);

  q_number &operator*=(q_number x);

  q_number &operator-=(q_number x);

  q_number &operator/=(q_number x);

  q_number &operator--();

  q_number &operator++();

  q_number operator--(int);

  q_number operator++(int);

  bool operator==(q_number x) const;

  bool operator!=(q_number x) const;

  bool operator<(q_number x) const;

  bool operator<=(q_number x) const;

  bool operator>(q_number x) const;

  bool operator>=(q_number x) const;

  z_number numerator() const;

  z_number denominator() const;

  z_number round_to_upper() const;

  z_number round_to_lower() const;

  void write(crab::crab_os &o) const;

}; // class q_number

inline crab::crab_os &operator<<(crab::crab_os &o, const z_number &z) {
  z.write(o);
  return o;
}

inline crab::crab_os &operator<<(crab::crab_os &o, const q_number &q) {
  q.write(o);
  return o;
}

/** for boost::hash_combine **/
inline std::size_t hash_value(const z_number &z) { return z.hash(); }

inline std::size_t hash_value(const q_number &q) { return q.hash(); }
} // namespace ikos

/** for specializations of std::hash **/
namespace std {
template <> struct hash<ikos::z_number> {
  size_t operator()(const ikos::z_number &z) const { return z.hash(); }
};

template <> struct hash<ikos::q_number> {
  size_t operator()(const ikos::q_number &q) const { return q.hash(); }
};
} // namespace std
