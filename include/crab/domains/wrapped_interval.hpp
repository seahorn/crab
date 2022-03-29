#pragma once

#include <crab/domains/interval.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/fixpoint/thresholds.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/numbers/wrapint.hpp>
#include <crab/support/os.hpp>

/**
 * Machine arithmetic intervals based on the paper
 * "Signedness-Agnostic Program Analysis: Precise Integer Bounds for
 * Low-Level Code" by J.A.Navas, P.Schachte, H.Sondergaard, and
 * P.J.Stuckey published in APLAS'12.
 **/

namespace crab {
namespace domains {

template <typename Number> class wrapped_interval {

  wrapint m_start;
  wrapint m_end;
  bool m_is_bottom;

  using wrapped_interval_t = wrapped_interval<Number>;
  wrapped_interval_t default_implementation(const wrapped_interval_t &x) const;
  wrapped_interval(wrapint start, wrapint stop, bool is_bottom);

  // nsplit in the APLAS'12 paper
  void signed_split(std::vector<wrapped_interval_t> &intervals) const;
  // ssplit in the APLAS'12 paper
  void unsigned_split(std::vector<wrapped_interval_t> &intervals) const;
  // cut in the APLAS'12 paper
  void signed_and_unsigned_split(std::vector<wrapped_interval_t> &out) const;
  wrapped_interval_t signed_mul(const wrapped_interval_t &x) const;
  wrapped_interval_t unsigned_mul(const wrapped_interval_t &x) const;
  // if out is empty then the intersection is empty
  void exact_meet(const wrapped_interval_t &x,
                  std::vector<wrapped_interval_t> &out) const;
  // Perform the reduced product of signed and unsigned multiplication.
  // It uses exact meet rather than abstract meet.
  void reduced_signed_unsigned_mul(const wrapped_interval_t &x,
                                   std::vector<wrapped_interval_t> &out) const;
  wrapped_interval_t unsigned_div(const wrapped_interval_t &x) const;
  wrapped_interval_t signed_div(const wrapped_interval_t &x) const;
  // This is sound only if wrapped interval defined over z_number.
  void trim_zero(std::vector<wrapped_interval_t> &out) const;
  wrapped_interval_t Shl(uint64_t k) const;
  wrapped_interval_t LShr(uint64_t k) const;
  wrapped_interval_t AShr(uint64_t k) const;

public:
  using bitwidth_t = wrapint::bitwidth_t;

  wrapped_interval();
  wrapped_interval(wrapint n);
  wrapped_interval(wrapint start, wrapint stop);

  // return top if n does not fit into a wrapint. No runtime errors.
  static wrapped_interval_t mk_winterval(Number n, bitwidth_t width);
  // Return top if lb or ub do not fit into a wrapint. No runtime errors.
  static wrapped_interval_t mk_winterval(Number lb, Number ub,
                                         bitwidth_t width);

  static wrapped_interval_t top();
  static wrapped_interval_t bottom();

  // return interval [0111...1, 1000....0]
  // In the APLAS'12 paper "signed limit" corresponds to "north pole".
  static wrapped_interval_t signed_limit(bitwidth_t b);
  // return interval [1111...1, 0000....0]
  // In the APLAS'12 paper "unsigned limit" corresponds to "south pole".
  static wrapped_interval_t unsigned_limit(bitwidth_t b);

  bool cross_signed_limit() const;
  bool cross_unsigned_limit() const;

  bitwidth_t get_bitwidth(int line) const;

  wrapint start() const;

  wrapint end() const;

  bool is_bottom() const;

  bool is_top() const;

  // Interpret wrapint as signed mathematical integers.
  ikos::interval<Number> to_interval() const;

  wrapped_interval_t lower_half_line(bool is_signed) const;

  wrapped_interval_t upper_half_line(bool is_signed) const;

  bool is_singleton() const;

  // Starting from m_start and going clock-wise x is encountered
  // before than m_stop.
  bool at(wrapint x) const;

  bool operator<=(const wrapped_interval_t &x) const;
  bool operator==(const wrapped_interval_t &x) const;
  bool operator!=(const wrapped_interval_t &x) const;
  wrapped_interval_t operator|(const wrapped_interval_t &x) const;
  wrapped_interval_t operator&(const wrapped_interval_t &x) const;
  wrapped_interval_t operator||(const wrapped_interval_t &x) const;
  wrapped_interval_t widening_thresholds(const wrapped_interval_t &x,
                                         const thresholds<Number> &ts) const;
  wrapped_interval_t operator&&(const wrapped_interval_t &x) const;

  wrapped_interval_t operator+(const wrapped_interval_t &x) const;
  wrapped_interval_t &operator+=(const wrapped_interval_t &x);
  wrapped_interval_t operator-() const;
  wrapped_interval_t operator-(const wrapped_interval_t &x) const;
  wrapped_interval_t &operator-=(const wrapped_interval_t &x);
  wrapped_interval_t operator*(const wrapped_interval_t &x) const;
  wrapped_interval_t &operator*=(const wrapped_interval_t &x);
  wrapped_interval_t operator/(const wrapped_interval_t &x) const {
    return SDiv(x);
  }
  wrapped_interval_t &operator/=(const wrapped_interval_t &x) {
    return this->operator=(this->operator/(x));
  }

  wrapped_interval_t SDiv(const wrapped_interval_t &x) const;
  wrapped_interval_t UDiv(const wrapped_interval_t &x) const;
  wrapped_interval_t SRem(const wrapped_interval_t &x) const;
  wrapped_interval_t URem(const wrapped_interval_t &x) const;

  wrapped_interval_t ZExt(unsigned bits_to_add) const;
  wrapped_interval_t SExt(unsigned bits_to_add) const;
  wrapped_interval_t Trunc(unsigned bits_to_keep) const;
  wrapped_interval_t Shl(const wrapped_interval_t &x) const;
  wrapped_interval_t LShr(const wrapped_interval_t &x) const;
  wrapped_interval_t AShr(const wrapped_interval_t &x) const;
  wrapped_interval_t And(const wrapped_interval_t &x) const;
  wrapped_interval_t Or(const wrapped_interval_t &x) const;
  wrapped_interval_t Xor(const wrapped_interval_t &x) const;

  void write(crab::crab_os &o) const;
  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const wrapped_interval<Number> &i) {
    i.write(o);
    return o;
  }
};
} // namespace domains
} // namespace crab

namespace ikos {
namespace linear_interval_solver_impl {
template <>
crab::domains::wrapped_interval<ikos::z_number>
mk_interval(ikos::z_number c, typename crab::wrapint::bitwidth_t w);

template <>
crab::domains::wrapped_interval<ikos::q_number>
mk_interval(ikos::q_number c, typename crab::wrapint::bitwidth_t w);

template <>
crab::domains::wrapped_interval<ikos::z_number>
trim_interval(const crab::domains::wrapped_interval<ikos::z_number> &i,
              const crab::domains::wrapped_interval<ikos::z_number> &j);

template <>
crab::domains::wrapped_interval<ikos::q_number>
trim_interval(const crab::domains::wrapped_interval<ikos::q_number> &i,
              const crab::domains::wrapped_interval<ikos::q_number> &j);

template <>
crab::domains::wrapped_interval<ikos::z_number>
lower_half_line(const crab::domains::wrapped_interval<ikos::z_number> &i,
                bool is_signed);

template <>
crab::domains::wrapped_interval<ikos::q_number>
lower_half_line(const crab::domains::wrapped_interval<ikos::q_number> &i,
                bool is_signed);

template <>
crab::domains::wrapped_interval<ikos::z_number>
upper_half_line(const crab::domains::wrapped_interval<ikos::z_number> &i,
                bool is_signed);

template <>
crab::domains::wrapped_interval<ikos::q_number>
upper_half_line(const crab::domains::wrapped_interval<ikos::q_number> &i,
                bool is_signed);
} // namespace linear_interval_solver_impl
} // end namespace ikos
