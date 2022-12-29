#pragma once

#include <crab/domains/lattice_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp> // for ikos::index_t

#include <boost/optional.hpp>

namespace crab {
namespace domains {
/*
 * This class models an abstract counter for variables. The typical
 * application is to count how many variables satisfy a property.
 * 
 * Initially, a counter can be set to 0. This class only
 * provides one operation to increment abstractly the counter.
 *
 * The possible abstract values for a counter are:
 * 
 *  {0, 1(V), [0,1](V), [1,+oo], [0,+oo]}
 *
 * where:
 *
 * 0       : the counter value is zero.
 * 1(V)    : the counter value is one and the variable is known to be V
 * [0,1](V): zero or one. If it is one then the variable is V
 * [1,+oo] : one or more
 * [0,+oo] : zero or more  (this is the top element)
 */
class small_range: public lattice_domain_api<small_range> {
  using kind_t = enum {
    Bottom,
    ExactlyZero,
    ExactlyOne,
    ZeroOrOne,
    ZeroOrMore,
    OneOrMore
  };

  // --- The counter value
  kind_t m_kind;
  // --- It represents V if [0,1] or 1
  boost::optional<ikos::index_t> m_value;
  
  small_range(kind_t k);
  small_range(kind_t k, ikos::index_t val);
  static small_range zeroOrMore();

  /** helpers for join and meet **/
  small_range join_zero_with(const small_range &other) const;
  small_range join_one_with(const small_range &other) const;
  small_range join_zero_or_one_with(const small_range &other) const;
  small_range join_one_or_more_with(const small_range &other) const;
  small_range meet_zero_with(const small_range &other) const;
  small_range meet_one_with(const small_range &other) const;
  small_range meet_zero_or_one_with(const small_range &other) const;
  small_range meet_one_or_more_with(const small_range &other) const;

public:
  small_range(); // counter initialized to top
  static small_range bottom();
  static small_range top();
  static small_range zero();
  static small_range oneOrMore();

  small_range(const small_range &other) = default;
  small_range(small_range &&other) = default;
  small_range &operator=(const small_range &other) = default;
  small_range &operator=(small_range &&other) = default;

  bool is_zero() const;
  bool is_one() const;

  /** lattice operations **/
  small_range make_top() const override;
  small_range make_bottom() const override;
  void set_to_top() override;
  void set_to_bottom() override;
  bool is_bottom() const override;
  bool is_top() const override;  
  bool operator<=(const small_range &other) const override;
  bool operator==(const small_range &other) const;
  void operator|=(const small_range &other) override;
  small_range operator|(const small_range &other) const override;
  small_range operator||(const small_range &other) const override;
  small_range operator&(const small_range &other) const override;
  small_range operator&&(const small_range &other) const override;

  /** Increment a counter **/
  template<typename VariableName>
  small_range increment(const VariableName &v) {
    if (!is_bottom()) {
      if (m_kind == ExactlyZero) {
	m_kind = ExactlyOne;
	m_value = v.index();
      } else if (m_kind == ExactlyOne || m_kind == ZeroOrMore ||
		 m_kind == ZeroOrOne || m_kind == OneOrMore) {
	if (!(m_kind == ExactlyOne && m_value.value() == v.index())) {
	  m_kind = OneOrMore;
	  m_value = boost::none;
	}
      } else {
	CRAB_ERROR("small_range::increment unreachable");
      }
    }
    return *this;    
  }

  std::string domain_name() const override;

  void write(crab_os &o) const override;  
  friend crab_os &operator<<(crab_os &o, const small_range &v) {
    v.write(o);
    return o;
  }
}; // end class small_range
} // end namespace domains
} // end namespace crab 
