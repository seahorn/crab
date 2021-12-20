#include <crab/domains/boolean.hpp>

namespace crab {
namespace domains {

boolean_value::boolean_value(kind_t v) : _value(v){};

boolean_value::boolean_value() : _value(Top) {}

boolean_value boolean_value::bottom() { return boolean_value(Bottom); }

boolean_value boolean_value::top() { return boolean_value(Top); }

boolean_value boolean_value::get_true() { return boolean_value(True); }

boolean_value boolean_value::get_false() { return boolean_value(False); }

boolean_value::boolean_value(const boolean_value &other)
    : _value(other._value) {}

boolean_value &boolean_value::operator=(const boolean_value &other) {
  if (this != &other)
    _value = other._value;
  return *this;
}

boolean_value boolean_value::make_bottom() const {
  boolean_value res(Bottom);
  return res;
}
boolean_value boolean_value::make_top() const {
  boolean_value res(Top);
  return res;
}
void boolean_value::set_to_top()  {
  _value = Top;
}
  void boolean_value::set_to_bottom()  {
  _value = Bottom;
}
  
bool boolean_value::is_bottom() const { return (_value == Bottom); }

bool boolean_value::is_top() const { return (_value == Top); }

bool boolean_value::is_true() const { return (_value == True); }

bool boolean_value::is_false() const { return (_value == False); }

bool boolean_value::operator<=(const boolean_value &other) const {

  if (_value == Bottom || other._value == Top)
    return true;
  else if (_value == Top)
    return (other._value == Top);
  else if (_value == True)
    return ((other._value == True) || (other._value == Top));
  else if (_value == False)
    return ((other._value == False) || (other._value == Top));
  // this should be unreachable
  return false;
}

bool boolean_value::operator==(const boolean_value &other) const {
  return (_value == other._value);
}
  
void boolean_value::operator|=(const boolean_value &other) {
  *this = *this | other;
}
  
boolean_value boolean_value::operator|(const boolean_value &other) const {
  if (is_bottom())
    return other;
  if (other.is_bottom())
    return *this;
  if (is_top() || other.is_top())
    return top();
  if (_value == other._value)
    return *this;

  // othewise true | false or false | true ==> top
  return top();
}

boolean_value boolean_value::operator&(const boolean_value &other) const {
  if (is_bottom())
    return *this;
  if (other.is_bottom())
    return other;
  if (is_top())
    return other;
  if (other.is_top())
    return *this;
  if (_value == other._value)
    return *this;

  // othewise true & false or false & true ==> bottom
  return bottom();
}

// the lattice satisfy ACC so join is the widening
boolean_value boolean_value::operator||(const boolean_value &other) const {
  return this->operator|(other);
}

// the lattice satisfy DCC so meet is the narrowing
boolean_value boolean_value::operator&&(const boolean_value &other) const {
  return this->operator&(other);
}

// Boolean operations
/*
          And  Or  X0r
     0 0   0   0    0
     0 1   0   1    1
     1 0   0   1    1
     1 1   1   1    0
     0 *   0   *    *
     * 0   0   *    *
     1 *   *   1    *
     * 1   *   1    *
     * *   *   *    *
*/

boolean_value boolean_value::And(boolean_value other) const {
  if (is_bottom() || other.is_bottom())
    return bottom();

  if (!is_top() && !other.is_top())
    return boolean_value(static_cast<kind_t>(static_cast<int>(_value) &
                                             static_cast<int>(other._value)));

  int x = static_cast<int>(_value);
  int y = static_cast<int>(other._value);
  if (x == 0 || y == 0)
    return get_false();
  else
    return top();
}

boolean_value boolean_value::Or(boolean_value other) const {
  if (is_bottom() || other.is_bottom())
    return bottom();

  if (!is_top() && !other.is_top())
    return boolean_value(static_cast<kind_t>(static_cast<int>(_value) |
                                             static_cast<int>(other._value)));

  int x = static_cast<int>(_value);
  int y = static_cast<int>(other._value);
  if (x == 1 || y == 1)
    return get_true();
  else
    return top();
}

boolean_value boolean_value::Xor(boolean_value other) const {
  if (is_bottom() || other.is_bottom())
    return bottom();

  if (!is_top() && !other.is_top())
    return boolean_value(static_cast<kind_t>(static_cast<int>(_value) ^
                                             static_cast<int>(other._value)));
  else
    return top();
}

boolean_value boolean_value::Negate() const {
  if (is_bottom())
    return bottom();
  if (_value == True)
    return get_false();
  if (_value == False)
    return get_true();
  return top();
}

void boolean_value::write(crab_os &o) const {
  switch (_value) {
  case Bottom:
    o << "_|_";
    break;
  case Top:
    o << "*";
    break;
  case True:
    o << "true";
    break;
  default: /*False*/
    o << "false";
  }
}

std::string boolean_value::domain_name() const {
  return "Bool";
}
  
} // end namespace domains
} // end namespace crab
