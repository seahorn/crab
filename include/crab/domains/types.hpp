#pragma once

/* Type lattice */

#include <crab/support/os.hpp>
#include <crab/types/variable.hpp>

#include <boost/optional.hpp>
#include <vector>

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
  
  type_value(bool is_bottom)
    : m_type(boost::none), m_is_bottom(is_bottom) {}
  
public:

  type_value(variable_type ty) : m_type(ty), m_is_bottom(false) {}

  static type_value bottom() { return type_value(true); }

  static type_value top() { return type_value(false); }

  bool is_bottom() const { return m_is_bottom; }

  bool is_top() const { return (!is_bottom() && !m_type); }

  bool operator<=(const type_value &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (o.is_bottom() || is_top()) {
      return false;
    } else {
      return ((get().is_region() && o.get().is_unknown_region()) ||
	      get() == o.get());
    }
  }

  bool operator==(const type_value &o) const {
    return *this <= o && o <= *this;
  }
  
  type_value operator|(const type_value &o) const {
    if (is_bottom() || o.is_top())
      return o;
    else if (is_top() || o.is_bottom())
      return *this;
    else {
      if (get().is_region() && o.get().is_unknown_region()) {
	return o;
      } else if (get().is_unknown_region() && o.get().is_region()) {
	return *this;
      } else if (get() == o.get()) {
	return *this;
      } else {
	return type_value::top();
      }
    }
  }

  type_value operator||(const type_value &o) const {
    return *this | o;
  }

  template <typename Thresholds>
  type_value widening_thresholds(const type_value &o,
				 const Thresholds &ts /*unused*/) const {
    return *this | o;
  }

  type_value operator&(const type_value &o) const {
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      if (get().is_region() && o.get().is_unknown_region()) {
	return *this;
      } else if (get().is_unknown_region() && o.get().is_region()) {
	return o;
      } else if (get() == o.get()) {
	return *this;
      } else {
	return type_value::bottom();
      }
    }
  }

  type_value operator&&(const type_value &o) const {
    return *this & o;
  }

  variable_type get() const {
    if (is_bottom()) {
      CRAB_ERROR("type_value::get() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("type_value::get() on top");
    }
    return *m_type;
  }

  void set(variable_type ty) {
    m_type = ty;
  }
    
  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      variable_type ty = get();
      if (ty.is_bool()) {
	o << "bool";
      } else if (ty.is_integer()) {
	o << "int" << ty.get_integer_bitwidth();
      } else if (ty.is_real()) {
	o << "real";
      } else if (ty.is_reference()) {
	o << "ref";
      } else if (ty.is_bool_array()) {
	o << "arr(bool)";
      } else if (ty.is_integer_array()) {
	o << "arr(int)";
      } else if (ty.is_real_array()) {
	o << "arr(real)";
      } else if (ty.is_unknown_region()) {
	o << "region(unknown)";
      } else if (ty.is_bool_region()) {
	o << "region(bool)";
      } else if (ty.is_integer_region()) {
	o << "region(int" << ty.get_integer_region_bitwidth() << ")";
      } else if (ty.is_real_region()) {
	o << "region(real)";
      } else if (ty.is_reference_region()) {
	o << "region(ref)";
      } else if (ty.is_bool_array_region()) {
	o << "region(arr(bool))";
      } else if (ty.is_int_array_region()) {
	o << "region(arr(int))";
      } else if (ty.is_real_array_region()) {
	o << "region(arr(real))";
      } else {
	// this shouldn't happen
	o << "top";
      }
    }
  }

  friend inline crab_os &operator<<(crab_os &o, const type_value &c) {
    c.write(o);
    return o;
  }
  
};
} // namespace domains
} // namespace crab
