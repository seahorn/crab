#include <crab/domains/types.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace domains {

type_value::type_value()
  : m_type(boost::none), m_is_bottom(false) {}
  
type_value::type_value(bool is_bottom)
  : m_type(boost::none), m_is_bottom(is_bottom) {}
  
type_value::type_value(variable_type ty) : m_type(ty), m_is_bottom(false) {}

type_value type_value::bottom() { return type_value(true); }

type_value type_value::top() { return type_value(false); }

bool type_value::is_bottom() const { return m_is_bottom; }

bool type_value::is_top() const { return (!is_bottom() && !m_type); }

type_value type_value::make_bottom() const {
  type_value res(true);
  return res;
}
type_value type_value::make_top() const {
  type_value res(false);
  return res;
}
void type_value::set_to_top()  {
  *this = make_top();
}
void type_value::set_to_bottom()  {
  *this = make_bottom();
}
  
bool type_value::operator<=(const type_value &o) const {
  if (is_bottom() || o.is_top()) {
    return true;
  } else if (o.is_bottom() || is_top()) {
    return false;
  } else {
    return ((get().is_region() && o.get().is_unknown_region()) ||
	    get() == o.get());
  }
}

bool type_value::operator==(const type_value &o) const {
  return *this <= o && o <= *this;
}

void type_value::operator|=(const type_value &o)  {
  *this = *this | o;   
}
  
type_value type_value::operator|(const type_value &o) const {
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

type_value type_value::operator||(const type_value &o) const {
  return *this | o;
}

template <typename Thresholds>
type_value type_value::widening_thresholds(const type_value &o,
					   const Thresholds &ts /*unused*/) const {
  return *this | o;
}

type_value type_value::operator&(const type_value &o) const {
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

type_value type_value::operator&&(const type_value &o) const {
  return *this & o;
}

variable_type type_value::get() const {
  if (is_bottom()) {
    CRAB_ERROR("type_value::get() on bottom");
  }
  if (is_top()) {
      CRAB_ERROR("type_value::get() on top");
  }
  return *m_type;
}

void type_value::set(variable_type ty) {
  m_type = ty;
}
    
void type_value::write(crab::crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else if (is_top()) {
    o << "top";
  } else {
    // variable_type::write is the single definition of how a type is printed.
    get().write(o);
  }  
}
  
std::string  type_value::domain_name() const {
  return "TypeDomain";
}  
} // end namespace domains
} // end namespace crab
