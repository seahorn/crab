#include <crab/domains/small_range.hpp>
#include <assert.h>

namespace crab {
namespace domains {

#define UNREACHABLE_BOTTOM \
  case Bottom: CRAB_ERROR("cannot be bottom");
  
small_range::small_range(kind_t k)
  : m_kind(k), m_value(boost::none) {
  if (k == kind_t::ExactlyOne || k == kind_t::ZeroOrOne) {
    CRAB_ERROR("cannot call small_range constructor without variable");
  }
}

small_range::small_range(kind_t k, ikos::index_t v)
  : m_kind(k), m_value(v) {
  if (k != kind_t::ExactlyOne && k != kind_t::ZeroOrOne) {
    CRAB_ERROR("call small_range constructor only with either ZeroOrOne or ExactlyOne");
  }
}  

small_range::small_range(): m_kind(ZeroOrMore), m_value(boost::none) {}
  
/*
   0 | 0       = 0
   0 | 1       = [0,1]
   0 | [0,1]   = [0,1]
   0 | [0,+oo] = [0,+oo]
   0 | [1,+oo] = [0,+oo]
*/
small_range small_range::join_zero_with(const small_range &other) const {
  assert(m_kind == ExactlyZero);
  assert(other.m_kind != Bottom);

  switch(other.m_kind) {
  case ExactlyZero:
    return other;
  case ExactlyOne:
    return small_range(ZeroOrOne, other.m_value.value());
  case ZeroOrOne:
    return other;
  case ZeroOrMore:
  case OneOrMore:  
    return zeroOrMore();
  UNREACHABLE_BOTTOM    
  }
}

/*
   1     | 0        = [0,1]
   1(V1) | 1(V2)    = 1       if V1==V2
                      [1,+oo] otherwise
   1(V1) | [0,1](V2)= [0,1]   if V1==V2
                      [0,+oo] otherwise
   1     | [0,+oo]  = [0,+oo]
   1     | [1,+oo]  = [1,+oo]
*/
small_range small_range::join_one_with(const small_range &other) const {
  assert(m_kind == ExactlyOne);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyZero:
    return small_range(ZeroOrOne, m_value.value());
  case ExactlyOne:
    if (m_value.value() == other.m_value.value()) {
      return *this;
    } else {
      return oneOrMore();
    }
  case ZeroOrOne: 
    if (m_value.value() == other.m_value.value()) {
      return other;
    } else {
      return zeroOrMore();
    }
  case ZeroOrMore:
  case OneOrMore:  
    return other;
  UNREACHABLE_BOTTOM    
  }
}

/*
   [0,1]     | 0         = [0,1]
   [0,1](V1) | 1(V2)     = [0,1]   if V1==V2
                         = [0,+oo] otherwise
   [0,1](V1) | [0,1](V2) = [0,1]   if V1==V2
                         = [0,+oo] otherwise
   [0,1]     | [0,+oo]   = [0,+oo]
   [0,1]     | [1,+oo]   = [0,+oo]
*/
small_range small_range::join_zero_or_one_with(const small_range &other) const {
  assert(m_kind == ZeroOrOne);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyZero:
    return *this;
  case ExactlyOne:
    BOOST_FALLTHROUGH;     
  case ZeroOrOne:
    if (m_value.value() == other.m_value.value()) {
      return *this;
    }
    BOOST_FALLTHROUGH;             
  case ZeroOrMore:
    BOOST_FALLTHROUGH;         
  case OneOrMore:  
    return zeroOrMore();
  UNREACHABLE_BOTTOM  
  }
}

/*
   [1,+oo] | 0 = [0,+oo]
   [1,+oo] | 1 = [1,+oo]
   [1,+oo] | [0,1] = [0,+oo]
   [1,+oo] | [0,+oo] = [0,+oo]
   [1,+oo] | [1,+oo] = [1,+oo]
*/
small_range small_range::join_one_or_more_with(const small_range &other) const {
  assert(m_kind == OneOrMore);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyOne:
    BOOST_FALLTHROUGH;                 
  case OneOrMore:
    return oneOrMore();
  case ExactlyZero:
    BOOST_FALLTHROUGH;                     
  case ZeroOrOne:
    BOOST_FALLTHROUGH;                     
  case ZeroOrMore:
    return zeroOrMore();
  UNREACHABLE_BOTTOM      
  }
}

/*
   0 & 0 = 0
   0 & 1 = _|_
   0 & [0,1] = 0
   0 & [0,+oo] = 0
   0 & [1,+oo] = _|_
*/
small_range small_range::meet_zero_with(const small_range &other) const {
  assert(m_kind == ExactlyZero);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyOne:
    BOOST_FALLTHROUGH;     
  case OneOrMore:
    return bottom();
  case ExactlyZero:
    BOOST_FALLTHROUGH;         
  case ZeroOrOne:
    BOOST_FALLTHROUGH;         
  case ZeroOrMore:
    return *this;
  UNREACHABLE_BOTTOM      
  }
}

/*
   1     & 0         = _|_
   1(V1) & 1(V2)     = 1 if V1==V2
                       _|_   otherwise
   1(V1) & [0,1](V2) = 1 if V1==V2
                       _|_   otherwise
   1     & [0,+oo]   = 1
   1     & [1,+oo]   = 1
*/
small_range small_range::meet_one_with(const small_range &other) const {
  assert(m_kind == ExactlyOne);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyZero:
    return bottom();
  case ExactlyOne:
    BOOST_FALLTHROUGH;             
  case ZeroOrOne: 
    if (m_value.value() == other.m_value.value()) {
      return *this;
    } else {
      return bottom();
    }
  case ZeroOrMore:
    BOOST_FALLTHROUGH;                 
  case OneOrMore:
    return *this;
  UNREACHABLE_BOTTOM      
  }
}

/*
   [0,1]     & 0         = 0
   [0,1](V1) & 1(V2)     = 1 if V1==V2
                         = _|_   otherwise
   [0,1](V1) & [0,1](V2) = [0,1] if V1==V2
                         = 0 otherwise 
   [0,1]     & [0,+oo]   = [0,1]
   [0,1]     & [1,+oo]   = 1 
*/
small_range small_range::meet_zero_or_one_with(const small_range &other) const {
  assert(m_kind == ZeroOrOne);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyZero:
    return other;
  case ExactlyOne: 
    if (m_value.value() == other.m_value.value()) {
      return other;
    } else {
      return bottom();
    }         
  case ZeroOrOne: 
    if (m_value.value() == other.m_value.value()) {
      return other;
    } else {
      return zero();
    }     
  case OneOrMore:
    return small_range(ExactlyOne, m_value.value());
  case ZeroOrMore:
    return *this;
  UNREACHABLE_BOTTOM      
  }
}

/*
   [1,+oo] & 0 = _|_
   [1,+oo] & 1 = 1
   [1,+oo] & [0,1] = 1
   [1,+oo] & [0,+oo] = [1,+oo]
   [1,+oo] & [1,+oo] = [1,+oo]
*/
small_range small_range::meet_one_or_more_with(const small_range &other) const {
  assert(m_kind == OneOrMore);
  assert(other.m_kind != Bottom);
  
  switch(other.m_kind) {
  case ExactlyZero:
    return bottom();
  case ExactlyOne:
    return other;
  case ZeroOrOne:
    return small_range(ExactlyOne, other.m_value.value());
  case ZeroOrMore:
    BOOST_FALLTHROUGH;  
  case OneOrMore:    
    return *this;
  UNREACHABLE_BOTTOM      
  } 
}


small_range small_range::bottom() { return small_range(Bottom); }
small_range small_range::top() { return small_range(ZeroOrMore); }
small_range small_range::zero() { return small_range(ExactlyZero); }
small_range small_range::oneOrMore() { return small_range(OneOrMore); }
small_range small_range::zeroOrMore() { return small_range(ZeroOrMore); }  
small_range small_range::make_bottom() const {
  small_range res(Bottom);
  return res;
}
small_range small_range::make_top() const { return zeroOrMore();}
  
void small_range::set_to_top()  {
  m_kind = ZeroOrMore;
  m_value = boost::none;
}
void small_range::set_to_bottom()  {
  m_kind = Bottom;
  m_value = boost::none;  
}  
    
bool small_range::is_bottom() const { return (m_kind == Bottom); }

bool small_range::is_top() const { return (m_kind == ZeroOrMore); }

bool small_range::is_zero() const { return (m_kind == ExactlyZero); }

bool small_range::is_one() const { return (m_kind == ExactlyOne); }

  
/*
      [0,+oo]
       /    \
      /      \
   [0,1](V2) [1,+oo]
    /  \    /   /
   0   1(V1)   /
    \   |     /
     \__|____/
        |
      Bottom
*/
bool small_range::operator<=(const small_range &other) const {
  if (*this == other) {
    return true;
  }

  if (is_bottom() || other.is_top()) {
    return true;
  }

  switch(m_kind) {
  case ExactlyZero:
    return other.m_kind != ExactlyOne && other.m_kind != OneOrMore;
  case ExactlyOne: {
    switch(other.m_kind) {
    case ExactlyZero:
      return false;
    case ExactlyOne:
      assert(m_value.value() != other.m_value.value());
      return false;
    case ZeroOrOne:
      return m_value.value() == other.m_value.value();
    case ZeroOrMore:
    case OneOrMore:
      return true;
    UNREACHABLE_BOTTOM        
    }
  }
  case ZeroOrOne:
    BOOST_FALLTHROUGH;    
  case OneOrMore:
    return other.m_kind == ZeroOrMore;
  case ZeroOrMore:
    assert(other.m_kind != ZeroOrMore);
    return false;
  UNREACHABLE_BOTTOM      
  }
  
}

bool small_range::operator==(const small_range &other) const {
  return (m_kind == other.m_kind && m_value == other.m_value);
}

void small_range::operator|=(const small_range &other) {
  *this = *this | other;     
}
  
small_range small_range::operator|(const small_range &other) const {
  if (is_bottom() || other.is_top()) {
    return other;
  } else if (other.is_bottom() || is_top()) {
    return *this;
  }

  if (is_zero()) {
    return join_zero_with(other);
  } else if (other.is_zero()) {
    return other.join_zero_with(*this);
  } else if (is_one()) {
    return join_one_with(other);
  } else if (other.is_one()) {
    return other.join_one_with(*this);
  } else if (m_kind == ZeroOrOne) {
    return join_zero_or_one_with(other);
  } else if (other.m_kind == ZeroOrOne) {
    return other.join_zero_or_one_with(*this);
  } else if (m_kind == OneOrMore) {
    return join_one_or_more_with(other);
  } else if (other.m_kind == OneOrMore) {
    return other.join_one_or_more_with(*this);
  } else {
    return small_range(ZeroOrMore);
  }
}

// the lattice has finite height so widening is the join
small_range small_range::operator||(const small_range &other) const {
  return *this | other;
}

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
  } else if (m_kind == ZeroOrOne) {
    return meet_zero_or_one_with(other);
  } else if (other.m_kind == ZeroOrOne) {
    return other.meet_zero_or_one_with(*this);
  } else if (m_kind == OneOrMore) {
    return meet_one_or_more_with(other);
  } else if (other.m_kind == OneOrMore) {
    return other.meet_one_or_more_with(*this);
  } else { // unreachable because top cases handled above
    CRAB_ERROR("unexpected small_range::meet operands");
  }
}

// the lattice has finite height so narrowing is the meet
small_range small_range::operator&&(const small_range &other) const {
  return *this & other;
}

void small_range::write(crab_os &o) const {
  switch (m_kind) {
  case Bottom:
    o << "_|_";
    break;
  case ExactlyZero:
    o << "[0,0]";
    break;
  case ExactlyOne:
    o << "[1,1](" << m_value.value() << ")";
    break;
  case ZeroOrOne:
    o << "[0,1](" << m_value.value() << ")";
    break;
  case ZeroOrMore:
    o << "[0,+oo]";
    break;
  case OneOrMore:
    o << "[1,+oo]";
    break;
  }
}

std::string small_range::domain_name() const {
  return "SmallRange";
}
} // end namespace domains
} // end namespace crab
