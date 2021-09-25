#pragma once

/* Type lattice */

#include <crab/support/os.hpp>
#include <crab/types/variable.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
/**
 *  Lattice for variable types which is almost flat except regions
 *  that can be unknown. An unknown region unifies with any other
 *  region.
 *
 *  Num    = bool | int | real
 *  Scalar = Num | reference
 * 
 *                 -------- top
 *                /         /   \
 *               /         /    region(unknown)-------- 
 *              /         /           |                |
 *            Scalar array(Num)    region(Scalar) region(array(Num))
 *               |          \      /                   |
 *               ---------- bottom----------------------
 *
 **/
class type_value {
  
  boost::optional<variable_type> m_type;
  bool m_is_bottom;
  
  type_value(bool is_bottom);
  
public:

  type_value(variable_type ty);

  static type_value bottom();

  static type_value top();

  bool is_bottom() const;

  bool is_top() const;

  bool operator<=(const type_value &o) const;

  bool operator==(const type_value &o) const;
  
  type_value operator|(const type_value &o) const;

  type_value operator||(const type_value &o) const;

  template <typename Thresholds>
  type_value widening_thresholds(const type_value &o,
				 const Thresholds &ts /*unused*/) const;

  type_value operator&(const type_value &o) const;

  type_value operator&&(const type_value &o) const;

  variable_type get() const;

  void set(variable_type ty);
    
  void write(crab::crab_os &o) const;

  friend inline crab_os &operator<<(crab_os &o, const type_value &c) {
    c.write(o);
    return o;
  }
  
};
} // namespace domains
} // namespace crab
