#pragma once

#include <crab/fixpoint/thresholds.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

template <typename Number> class constant {
  using constant_t = constant<Number>;

  boost::optional<Number> m_constant;
  bool m_is_bottom;

  constant(bool is_bottom);

public:
  constant(Number c);

  static constant_t bottom();
  static constant_t top();
  static constant_t zero();

  bool is_bottom() const;
  bool is_top() const;

  bool is_constant() const;
  Number get_constant() const;

  /** lattice operations **/
  bool operator<=(const constant_t &o) const;
  bool operator==(const constant_t &o) const;
  constant_t operator|(const constant_t &o) const;
  constant_t operator||(const constant_t &o) const;
  template <typename Thresholds>
  constant_t widening_thresholds(const constant_t &o,
                                 const Thresholds &ts /*unused*/) const {
    return *this | o;
  }
  constant_t operator&(const constant_t &o) const;
  constant_t operator&&(const constant_t &o) const;

  /** arithmetic operations **/
  constant_t Add(const constant_t &o) const;
  constant_t Sub(const constant_t &o) const;
  constant_t Mul(const constant_t &o) const;
  constant_t SDiv(const constant_t &o) const;
  constant_t SRem(const constant_t &o) const;
  constant_t UDiv(const constant_t &o) const;
  constant_t URem(const constant_t &o) const;

  /** bitwise operations **/
  constant_t BitwiseAnd(const constant_t &o) const;
  constant_t BitwiseOr(const constant_t &o) const;
  constant_t BitwiseXor(const constant_t &o) const;
  constant_t BitwiseShl(const constant_t &o) const;
  constant_t BitwiseLShr(const constant_t &o) const;
  constant_t BitwiseAShr(const constant_t &o) const;

  void write(crab::crab_os &o) const;
  friend crab_os &operator<<(crab_os &o, const constant_t &c) {
    c.write(o);
    return o;
  }
};

} // namespace domains
} // namespace crab
