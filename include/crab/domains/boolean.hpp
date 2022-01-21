#pragma once

/* A 3-valued boolean value */

#include <crab/support/debug.hpp>
#include <crab/domains/lattice_domain.hpp>

namespace crab {
namespace domains {

class boolean_value: public lattice_domain_api<boolean_value> {
  /*
            Top
            / \
           /   \
        True  False
           \   /
            \ /
           Bottom
  */
  using kind_t = enum { False = 0x0, True = 0x1, Bottom = 0x2, Top = 0x3 };

  kind_t _value;

  boolean_value(kind_t v);

public:
  boolean_value();

  static boolean_value bottom();

  static boolean_value top();

  boolean_value make_top() const override;

  boolean_value make_bottom() const override;

  void set_to_top() override;

  void set_to_bottom() override;
  
  static boolean_value get_true();

  static boolean_value get_false();

  boolean_value(const boolean_value &other);

  boolean_value &operator=(const boolean_value &other);

  bool is_bottom() const override;

  bool is_top() const override;

  bool is_true() const;

  bool is_false() const;

  bool operator<=(const boolean_value &other) const override;

  bool operator==(const boolean_value &other) const;

  void operator|=(const boolean_value &other) override;
  
  boolean_value operator|(const boolean_value &other) const override;

  boolean_value operator&(const boolean_value &other) const override;

  // the lattice satisfy ACC so join is the widening
  boolean_value operator||(const boolean_value &other) const override;

  // the lattice satisfy DCC so meet is the narrowing
  boolean_value operator&&(const boolean_value &other) const override;

  boolean_value And(boolean_value other) const;

  boolean_value Or(boolean_value other) const;

  boolean_value Xor(boolean_value other) const;

  boolean_value Negate() const;

  void write(crab_os &o) const override;

  std::string domain_name() const override;
}; // end class boolean_value

inline crab_os &operator<<(crab_os &o, const boolean_value &v) {
  v.write(o);
  return o;
}
  
} //end namespace domains
} //end namespace crab
