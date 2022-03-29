/*******************************************************************************
 *
 * Standard domain of numerical congruences extended with bitwise
 * operations.
 *
 * Author: Alexandre C. D. Wimmers (alexandre.c.wimmers@nasa.gov)
 *
 * Contributors: Jorge A. Navas (jorge.a.navaslaserna@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#pragma once

#include <crab/support/os.hpp>

#include <boost/optional.hpp>

namespace ikos {

template <typename Number> class congruence {
  using congruence_t = congruence<Number>;

  bool m_is_bottom;

  /// A congruence is denoted by aZ + b, where b \in Z and a \in N.
  /// The abstract state aZ + b represents all numbers that are
  /// congruent to b modulo a.
  Number m_a; // modulo
  Number m_b; // remainder

  void normalize(void);

  congruence(bool b);

  congruence(int n);

  congruence(Number a, Number b);

  Number abs(Number x) const;

  Number max(Number x, Number y) const;

  Number min(Number x, Number y) const;

  Number gcd(Number x, Number y, Number z) const;

  Number gcd_helper(Number x, Number y) const;

  Number gcd(Number x, Number y) const;

  Number lcm(Number x, Number y) const;

  bool is_zero() const;

  bool all_ones() const;

public:
  congruence();
  congruence(Number n);
  congruence(const congruence_t &o) = default;
  congruence_t &operator=(const congruence_t &o) = default;
  congruence(congruence_t &&o) = default;
  congruence_t &operator=(congruence_t &&o) = default;

  static congruence_t top();
  static congruence_t bottom();

  bool is_bottom() const;
  bool is_top() const;

  boost::optional<Number> singleton() const;

  Number get_modulo() const;
  Number get_remainder() const;

  bool operator==(const congruence_t &o) const;
  bool operator!=(const congruence_t &x) const;

  /** Lattice Operations **/
  bool operator<=(const congruence_t &o) const;
  congruence_t operator|(const congruence_t &o) const;
  congruence_t operator&(const congruence_t &o) const;
  congruence_t operator||(const congruence_t &o) const;
  congruence_t operator&&(const congruence_t &o) const;

  /** Arithmetic Operators **/
  congruence_t operator+(const congruence_t &o) const;
  congruence_t operator-(const congruence_t &o) const;
  congruence_t operator-() const;
  congruence_t operator*(const congruence_t &o) const;
  // signed division
  congruence_t operator/(const congruence_t &o) const;
  // signed remainder operator
  congruence_t operator%(const congruence_t &o) const;
  congruence_t SDiv(const congruence_t &x) const;
  congruence_t UDiv(const congruence_t &x) const;
  congruence_t SRem(const congruence_t &x) const;
  congruence_t URem(const congruence_t &x) const;

  /**
   * Bitwise operators.
   * They are very imprecise because we ignore bitwidth.
   *
   * Bitwise operation can be implemented more precisely based on
   * Stefan Bygde's paper: Static WCET analysis based on abstract
   * interpretation and counting of elements, Vasteras : School of
   * Innovation, Design and Engineering, Malardalen University (2010).
   **/

  congruence_t And(const congruence_t &o) const;
  congruence_t Or(const congruence_t &o) const;
  congruence_t Xor(const congruence_t &o) const;
  congruence_t Shl(const congruence_t &o) const;
  congruence_t AShr(const congruence_t &o) const;
  congruence_t LShr(const congruence_t &o) const;

  void write(crab::crab_os &o) const;
  friend crab::crab_os &operator<<(crab::crab_os &o, const congruence_t &c) {
    c.write(o);
    return o;
  }

}; // end class congruence

template <typename Number>
inline congruence<Number> operator+(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) + x;
}

template <typename Number>
inline congruence<Number> operator+(const congruence<Number> &x, Number c) {
  return x + congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator*(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) * x;
}

template <typename Number>
inline congruence<Number> operator*(const congruence<Number> &x, Number c) {
  return x * congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator/(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) / x;
}

template <typename Number>
inline congruence<Number> operator/(const congruence<Number> &x, Number c) {
  return x / congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator-(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) - x;
}

template <typename Number>
inline congruence<Number> operator-(const congruence<Number> &x, Number c) {
  return x - congruence<Number>(c);
}

} // namespace ikos
