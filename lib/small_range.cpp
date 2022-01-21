#include <crab/domains/small_range.hpp>
#include <crab/support/debug.hpp>
#include <assert.h>

namespace crab {
namespace domains {

small_range::small_range(kind_t v) : m_value(v){};

/*
   [0,0] | [0,0] = [0,0]
   [0,0] | [1,1] = [0,1]
   [0,0] | [0,1] = [0,1]
   [0,0] | [0,+oo] = [0,+oo]
   [0,0] | [1,+oo] = [0,+oo]
*/
small_range small_range::join_zero_with(const small_range &other) const {
  if (other.m_value == ExactlyZero) {
    return other;
  } else if (other.m_value == ExactlyOne) {
    return small_range(ZeroOrOne);
    } else if (other.m_value == ZeroOrOne) {
    return other;
  } else {
    return small_range(ZeroOrMore);
  }
}

/*
   [1,1] | [0,0] = [0,1]
   [1,1] | [1,1] = [1,+oo]
   [1,1] | [0,1] = [0,+oo]
   [1,1] | [0,+oo] = [0,+oo]
   [1,1] | [1,+oo] = [1,+oo]
*/
small_range small_range::join_one_with(const small_range &other) const {
  if (other.m_value == ExactlyZero) {
    return small_range(ZeroOrOne);
  } else if (other.m_value == ExactlyOne) {
      return small_range(OneOrMore);
  } else if (other.m_value == OneOrMore) {
    return other;
  } else {
    return small_range(ZeroOrMore);
  }
}

/*
   [0,1] | [0,0] = [0,1]
   [0,1] | [1,1] = [0,+oo]
   [0,1] | [0,1] = [0,+oo]
   [0,1] | [0,+oo] = [0,+oo]
   [0,1] | [1,+oo] = [0,+oo]
*/
small_range small_range::join_zero_or_one_with(const small_range &other) const {
  if (other.m_value == ExactlyZero) {
    return small_range(ZeroOrOne);
  } else {
    return small_range(ZeroOrMore);
  }
}

/*
   [1,+oo] | [0,0] = [0,+oo]
   [1,+oo] | [1,1] = [1,+oo]
   [1,+oo] | [0,1] = [0,+oo]
   [1,+oo] | [0,+oo] = [0,+oo]
   [1,+oo] | [1,+oo] = [1,+oo]
*/
small_range small_range::join_one_or_more_with(const small_range &other) const {
  if (other.m_value == ExactlyOne || other.m_value == OneOrMore) {
    return small_range(OneOrMore);
  } else {
    return small_range(ZeroOrMore);
    }
}

/*
   [0,0] & [0,0] = [0,0]
   [0,0] & [1,1] = _|_
   [0,0] & [0,1] = [0,0]
   [0,0] & [0,+oo] = [0,0]
   [0,0] & [1,+oo] = _|_
*/
small_range small_range::meet_zero_with(const small_range &other) const {
  assert(is_zero());
  if (other.m_value == ExactlyOne || other.m_value == OneOrMore) {
    return small_range::bottom();
  } else {
    return *this;
  }
}

/*
   [1,1] & [0,0] = _|_
   [1,1] & [1,1] = [1,1]
   [1,1] & [0,1] = [1,1]
   [1,1] & [0,+oo] = [1,1]
   [1,1] & [1,+oo] = [1,1]
*/
small_range small_range::meet_one_with(const small_range &other) const {
  assert(is_one());
  if (other.m_value == ExactlyZero) {
    return small_range::bottom();
    } else {
    return *this;
  }
}

/*
   [0,1] & [0,0] = [0,0]
   [0,1] & [1,1] = [1,1]
   [0,1] & [0,1] = [0,1]
   [0,1] & [0,+oo] = [0,1]
   [0,1] & [1,+oo] = [1,1]
*/
small_range small_range::meet_zero_or_one_with(const small_range &other) const {
  assert(m_value == ZeroOrOne);
  if (other.is_zero() || other.is_one()) {
    return other;
  } else if (other.m_value == OneOrMore) {
    return one();
  } else {
    return *this;
  }
}

/*
   [1,+oo] & [0,0] = _|_
   [1,+oo] & [1,1] = [1,1]
   [1,+oo] & [0,1] = [1,1]
   [1,+oo] & [0,+oo] = [1,+oo]
   [1,+oo] & [1,+oo] = [1,+oo]
*/
small_range small_range::meet_one_or_more_with(const small_range &other) const {
  if (other.is_zero()) {
    return small_range::bottom();
  } else if (other.is_one()) {
    return other;
  } else if (other.m_value == ZeroOrOne) {
      return one();
  } else {
    assert(other.is_top() || other.m_value == OneOrMore);
    return *this;
  }
}

small_range::small_range() : m_value(ZeroOrMore) {}

small_range small_range::bottom() { return small_range(Bottom); }
small_range small_range::top() { return small_range(ZeroOrMore); }
small_range small_range::zero() { return small_range(ExactlyZero); }
small_range small_range::one() { return small_range(ExactlyOne); }
small_range small_range::oneOrMore() { return small_range(OneOrMore); }
small_range small_range::make_bottom() const {
  small_range res(Bottom);
  return res;
}
small_range small_range::make_top() const {
  small_range res(ZeroOrMore);
  return res;
}
void small_range::set_to_top()  {
  m_value = ZeroOrMore;
}
void small_range::set_to_bottom()  {
  m_value = Bottom;
}  
    
bool small_range::is_bottom() const { return (m_value == Bottom); }

bool small_range::is_top() const { return (m_value == ZeroOrMore); }

bool small_range::is_zero() const { return (m_value == ExactlyZero); }

bool small_range::is_one() const { return (m_value == ExactlyOne); }

/*
  [0,+oo]
    |   \
    |    \
   [0,1] [1,+oo]
    / \ /
   0   1
*/
bool small_range::operator<=(const small_range &other) const {
  if (m_value == other.m_value) {
    return true;
  } else if (m_value == Bottom || other.m_value == ZeroOrMore) {
    return true;
  } else if (m_value == ExactlyZero) {
    return other.m_value != ExactlyOne && other.m_value != OneOrMore;
  } else if (m_value == ExactlyOne) {
    return other.m_value != ExactlyZero;
  } else if (m_value == ZeroOrOne || m_value == OneOrMore) {
    return other.m_value == ZeroOrMore;
  } else if (m_value == ZeroOrMore) {
    assert(other.m_value != ZeroOrMore);
    return false;
  }
  // should be unreachable
  return false;
}

bool small_range::operator==(const small_range &other) const { return m_value == other.m_value; }

void small_range::operator|=(const small_range &other) {
  *this = *this | other;     
}
  
small_range small_range::operator|(const small_range &other) const {
  if (is_bottom()) {
    return other;
  } else if (other.is_bottom()) {
    return *this;
  } else if (is_zero()) {
    return join_zero_with(other);
  } else if (other.is_zero()) {
    return join_zero_with(*this);
  } else if (is_one()) {
    return join_one_with(other);
  } else if (other.is_one()) {
    return join_one_with(*this);
  } else if (m_value == ZeroOrOne) {
    return join_zero_or_one_with(other);
  } else if (other.m_value == ZeroOrOne) {
    return join_zero_or_one_with(*this);
  } else if (m_value == OneOrMore) {
    return join_one_or_more_with(other);
  } else if (other.m_value == OneOrMore) {
    return join_one_or_more_with(*this);
  } else {
    return small_range(ZeroOrMore);
  }
}

small_range small_range::operator||(const small_range &other) const { return *this | other; }

small_range small_range::operator&(const small_range &other) const {
  if (is_bottom() || other.is_top()) {
    return *this;
  } else if (other.is_bottom() || is_top()) {
    return other;
  }

  if (is_zero()) {
    return meet_zero_with(other);
  } else if (other.is_zero()) {
    return other.meet_zero_with(*this);
  } else if (is_one()) {
    return meet_one_with(other);
  } else if (other.is_one()) {
    return other.meet_one_with(*this);
  } else if (m_value == ZeroOrOne) {
    return meet_zero_or_one_with(other);
  } else if (other.m_value == ZeroOrOne) {
    return other.meet_zero_or_one_with(*this);
  } else if (m_value == OneOrMore) {
    return meet_one_or_more_with(other);
  } else if (other.m_value == OneOrMore) {
    return other.meet_one_or_more_with(*this);
  } else { // unreachable because top cases handled above
    CRAB_ERROR("unexpected small_range::meet operands");
  }
}

small_range small_range::operator&&(const small_range &other) const { return *this & other; }

small_range small_range::increment(void) {
  if (!is_bottom()) {
    if (m_value == ExactlyZero) {
      m_value = ExactlyOne;
    } else if (m_value == ExactlyOne || m_value == ZeroOrMore ||
	       m_value == ZeroOrOne || m_value == OneOrMore) {
      m_value = OneOrMore;
    } else {
      CRAB_ERROR("small_range::increment unreachable");
    }
  }
  return *this;
}

void small_range::write(crab_os &o) const {
  switch (m_value) {
  case Bottom:
    o << "_|_";
    break;
  case ExactlyZero:
    o << "[0,0]";
    break;
  case ExactlyOne:
    o << "[1,1]";
    break;
  case ZeroOrOne:
    o << "[0,1]";
    break;
  case ZeroOrMore:
    o << "[0,+oo]";
    break;
  case OneOrMore:
    o << "[1,+oo]";
    break;
  default:
    CRAB_ERROR("unexpected small_range value");
  }
}

std::string small_range::domain_name() const {
  return "SmallRange";
}
} // end namespace domains
} // end namespace crab
