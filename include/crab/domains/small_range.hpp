#pragma once

#include <crab/support/os.hpp>
#include <crab/domains/lattice_domain.hpp>

namespace crab {
namespace domains {
/*
 * Class that represents the range of intervals
 * {[0,0], [1,1], [0,1], [0,+oo], [1,+oo]}
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

  kind_t m_value;

  small_range(kind_t v);

  /*
      [0,0] | [0,0] = [0,0]
      [0,0] | [1,1] = [0,1]
      [0,0] | [0,1] = [0,1]
      [0,0] | [0,+oo] = [0,+oo]
      [0,0] | [1,+oo] = [0,+oo]
  */
  small_range join_zero_with(const small_range &other) const;

  /*
      [1,1] | [0,0] = [0,1]
      [1,1] | [1,1] = [1,+oo]
      [1,1] | [0,1] = [0,+oo]
      [1,1] | [0,+oo] = [0,+oo]
      [1,1] | [1,+oo] = [1,+oo]
  */
  small_range join_one_with(const small_range &other) const;

  /*
      [0,1] | [0,0] = [0,1]
      [0,1] | [1,1] = [0,+oo]
      [0,1] | [0,1] = [0,+oo]
      [0,1] | [0,+oo] = [0,+oo]
      [0,1] | [1,+oo] = [0,+oo]
  */
  small_range join_zero_or_one_with(const small_range &other) const;

  /*
      [1,+oo] | [0,0] = [0,+oo]
      [1,+oo] | [1,1] = [1,+oo]
      [1,+oo] | [0,1] = [0,+oo]
      [1,+oo] | [0,+oo] = [0,+oo]
      [1,+oo] | [1,+oo] = [1,+oo]
  */
  small_range join_one_or_more_with(const small_range &other) const;

  /*
      [0,0] & [0,0] = [0,0]
      [0,0] & [1,1] = _|_
      [0,0] & [0,1] = [0,0]
      [0,0] & [0,+oo] = [0,0]
      [0,0] & [1,+oo] = _|_
  */
  small_range meet_zero_with(const small_range &other) const;

  /*
      [1,1] & [0,0] = _|_
      [1,1] & [1,1] = [1,1]
      [1,1] & [0,1] = [1,1]
      [1,1] & [0,+oo] = [1,1]
      [1,1] & [1,+oo] = [1,1]
  */
  small_range meet_one_with(const small_range &other) const;

  /*
      [0,1] & [0,0] = [0,0]
      [0,1] & [1,1] = [1,1]
      [0,1] & [0,1] = [0,1]
      [0,1] & [0,+oo] = [0,1]
      [0,1] & [1,+oo] = [1,1]
  */
  small_range meet_zero_or_one_with(const small_range &other) const;

  /*
      [1,+oo] & [0,0] = _|_
      [1,+oo] & [1,1] = [1,1]
      [1,+oo] & [0,1] = [1,1]
      [1,+oo] & [0,+oo] = [1,+oo]
      [1,+oo] & [1,+oo] = [1,+oo]
  */
  small_range meet_one_or_more_with(const small_range &other) const;

public:
  small_range();

  static small_range bottom();
  static small_range top();
  static small_range zero();
  static small_range one();
  static small_range oneOrMore();

  small_range(const small_range &other) = default;
  small_range(small_range &&other) = default;
  small_range &operator=(const small_range &other) = default;
  small_range &operator=(small_range &&other) = default;

  small_range make_top() const override;

  small_range make_bottom() const override;

  void set_to_top() override;

  void set_to_bottom() override;
  
  bool is_bottom() const override;

  bool is_top() const override;

  bool is_zero() const;

  bool is_one() const;

  /*
     [0,+oo]
       |   \
       |    \
      [0,1] [1,+oo]
       / \ /
      0   1
   */
  bool operator<=(const small_range &other) const override;

  bool operator==(const small_range &other) const;

  void operator|=(const small_range &other) override;
  
  small_range operator|(const small_range &other) const override;

  small_range operator||(const small_range &other) const override;

  small_range operator&(const small_range &other) const override;

  small_range operator&&(const small_range &other) const override;

  small_range increment(void);

  void write(crab_os &o) const override;

  std::string domain_name() const override;
  
  friend crab_os &operator<<(crab_os &o, const small_range &v) {
    v.write(o);
    return o;
  }
}; // end class small_range
} // end namespace domains
} // end namespace crab 
