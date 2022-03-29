#pragma once

#include <crab/domains/interval.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {
namespace domains {

enum class sign_interval {
  BOT, // empty interval
  LTZ, // [-oo,0)
  GTZ, // (0, +oo]
  EQZ, // [0,0]
  NEZ, // [-oo,0) U (0, +oo]
  GEZ, // [0, +oo]
  LEZ, // [-oo,0]
  TOP, // [-oo,+oo]
       /*
                TOP
               / | \
              /  |  \
            LEZ NEZ  GEZ
             |\ /  \/|
             |/ EQZ \|
            LTZ  |  GTZ
              \  |  /
               \ | /
                BOT
        */
};

template <typename Number> class sign {
  static_assert(std::is_same<Number, ikos::z_number>::value,
                "Class sign only defined over ikos::z_number");

  using sign_t = sign<Number>;
  using interval_t = ikos::interval<Number>;
  using bound_t = ikos::bound<Number>;
  sign_interval m_sign;

  constexpr Number get_zero() const { return Number(0); }
  constexpr Number get_plus_one() const { return Number(1); }
  constexpr Number get_minus_one() const { return Number(-1); }

  sign_t shiftOp(const sign_t &o) const;

  sign_t defaultOp(const sign_t &o) const;

  explicit sign(sign_interval s);

public:
  explicit sign(bool is_bottom);
  explicit sign(Number c);

  sign_t from_interval(const interval_t &i) const;
  interval_t to_interval() const;

  static sign_t bottom();
  static sign_t top();
  static sign_t mk_equal_zero();
  static sign_t mk_less_than_zero();
  static sign_t mk_greater_than_zero();
  static sign_t mk_less_or_equal_than_zero();
  static sign_t mk_greater_or_equal_than_zero();
  static sign_t mk_not_equal_zero();

  bool is_bottom() const;
  bool is_top() const;
  bool equal_zero() const;
  bool less_than_zero() const;
  bool greater_than_zero() const;
  bool less_or_equal_than_zero() const;
  bool greater_or_equal_than_zero() const;
  bool not_equal_zero() const;

  // lattice operations
  bool operator<=(const sign_t &o) const;
  bool operator==(const sign_t &o) const;
  sign_t operator|(const sign_t &o) const;
  sign_t operator&(const sign_t &o) const;

  // addition
  sign_t operator+(const sign_t &o) const;
  // subtraction
  sign_t operator-(const sign_t &o) const;
  // multiplication
  sign_t operator*(const sign_t &o) const;
  // signed division
  sign_t operator/(const sign_t &o) const;
  // division and remainder operations
  sign_t UDiv(const sign_t &o) const;
  sign_t SRem(const sign_t &o) const;
  sign_t URem(const sign_t &o) const;
  // bitwise operations
  sign_t And(const sign_t &o) const;
  sign_t Or(const sign_t &o) const;
  sign_t Xor(const sign_t &o) const;
  sign_t Shl(const sign_t &o) const;
  sign_t LShr(const sign_t &o) const;
  sign_t AShr(const sign_t &o) const;

  void write(crab::crab_os &o) const;

  friend inline crab_os &operator<<(crab_os &o, const sign_t &c) {
    c.write(o);
    return o;
  }
};
} // namespace domains
} // namespace crab
