#pragma once

/* Type lattice */

#include <crab/domains/lattice_domain.hpp>
#include <crab/support/os.hpp>
#include <crab/types/variable.hpp>

#include <boost/optional.hpp>
#include <string>

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
class type_value: public lattice_domain_api<type_value> {
  
  boost::optional<variable_type> m_type;
  bool m_is_bottom;
  
  type_value(bool is_bottom);
  
public:

  type_value();
  
  type_value(variable_type ty);

  static type_value bottom();

  static type_value top();

  type_value make_top() const override;

  type_value make_bottom() const override;

  void set_to_top() override;

  void set_to_bottom() override;
  
  bool is_bottom() const override;

  bool is_top() const override;

  bool operator<=(const type_value &o) const override;

  bool operator==(const type_value &o) const;

  void operator|=(const type_value &o) override;
  
  type_value operator|(const type_value &o) const override;

  type_value operator||(const type_value &o) const override;

  template <typename Thresholds>
  type_value widening_thresholds(const type_value &o,
				 const Thresholds &ts /*unused*/) const;

  type_value operator&(const type_value &o) const override;

  type_value operator&&(const type_value &o) const override;

  variable_type get() const;

  void set(variable_type ty);
    
  void write(crab::crab_os &o) const override;

  friend inline crab_os &operator<<(crab_os &o, const type_value &c) {
    c.write(o);
    return o;
  }

  std::string domain_name() const override;
  
};
} // namespace domains
} // namespace crab
