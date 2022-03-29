#include <crab/domains/congruence.hpp>
#include <crab/domains/interval.hpp>

namespace crab {
namespace domains {

/*
 *  The reduce operator based on "Static Analysis of Arithmetical
 *  Congruences" by P. Granger published in International Journal of
 *  Computer Mathematics, 1989.
 */
template <typename Number> class interval_congruence {
  using interval_congruence_t = interval_congruence<Number>;
  using interval_t = ikos::interval<Number>;
  using congruence_t = ikos::congruence<Number>;
  using bound_t = ikos::bound<Number>;

  interval_t m_first;
  congruence_t m_second;

  interval_congruence(bool is_bottom);

  Number abs(Number x) const;

  Number mod(Number a, Number b) const;

  // R(c,a) is the least element of c greater or equal than a
  Number R(congruence_t c, Number a) const;

  // L(c,a) is the greatest element of c smaller or equal than a
  Number L(congruence_t c, Number a) const;

public:
  interval_congruence(Number n);

  interval_congruence(interval_t &&i, congruence_t &&c);

  interval_congruence(interval_t &&i);

  interval_congruence(congruence_t &&c);

  interval_congruence(const interval_congruence &other) = default;
  interval_congruence(interval_congruence &&other) = default;
  interval_congruence_t &
  operator=(const interval_congruence_t &other) = default;
  interval_congruence_t &operator=(interval_congruence_t &&other) = default;

  static interval_congruence_t top();
  static interval_congruence_t bottom();

  bool is_bottom();
  bool is_top();

  interval_t &first();
  const interval_t &first() const;

  congruence_t &second();
  const congruence_t &second() const;

  /*
     Let (i,c) be a pair of interval and congruence these are the
     main rules described by Granger:

     if (c.is_bottom() || i.is_bottom()) (bottom(), bottom());
     if (c = 0Z+a and a \notin i)        (bottom(), bottom());
     if (c = 0Z+a)                       ([a,a]   , c);
     if (i=[a,b] and R(c,a) > L(c,b))    (bottom(), bottom());
     if (i=[a,b])                        ([R(c,a), L(c,b)], c);
     if (i=[a,+oo])                      ([R(c,a), +oo], c);
     if (i=[-oo,b])                      ([-oo, L(c,b)], c);
     otherwise                           (i,c)
   */

  void reduce();

  // join
  interval_congruence_t operator|(const interval_congruence_t &x) const;
  // meet
  interval_congruence_t operator&(const interval_congruence_t &x) const;
  // arithmetic/bitwise operations
  interval_congruence_t operator+(const interval_congruence_t &x) const;
  interval_congruence_t operator-(const interval_congruence_t &x) const;
  interval_congruence_t operator*(const interval_congruence_t &x) const;
  interval_congruence_t operator/(const interval_congruence_t &x) const;
  interval_congruence_t SDiv(const interval_congruence_t &x) const;
  interval_congruence_t UDiv(const interval_congruence_t &x) const;
  interval_congruence_t SRem(const interval_congruence_t &x) const;
  interval_congruence_t URem(const interval_congruence_t &x) const;
  interval_congruence_t Trunc(unsigned width) const;
  interval_congruence_t ZExt(unsigned width) const;
  interval_congruence_t SExt(unsigned width) const;
  interval_congruence_t And(const interval_congruence_t &x) const;
  interval_congruence_t Or(const interval_congruence_t &x) const;
  interval_congruence_t Xor(const interval_congruence_t &x) const;
  interval_congruence_t Shl(const interval_congruence_t &x) const;
  interval_congruence_t LShr(const interval_congruence_t &x) const;
  interval_congruence_t AShr(const interval_congruence_t &x) const;

  void write(crab_os &o) const;

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const interval_congruence<Number> &v) {
    v.write(o);
    return o;
  }
};

} // end namespace domains
} // end namespace crab
