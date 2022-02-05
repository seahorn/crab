#pragma once

#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/types.hpp>

namespace crab {
namespace domains {
namespace region_domain_impl {

/* This class contains any non-relational information about region
   variables */
class region_info {
  // Simple reduced product of (SmallRange x Boolean x Type)
  using product_t =
      basic_domain_product2<small_range,
                            basic_domain_product2<boolean_value, type_value>>;

  product_t m_product;
  region_info(product_t &&product);

public:
  region_info();
  region_info(small_range count, boolean_value init, type_value type);
  region_info(const region_info &other) = default;
  region_info(region_info &&other) = default;
  region_info &operator=(const region_info &other) = default;
  region_info &operator=(region_info &&other) = default;

  static region_info bottom();
  static region_info top();
  bool is_bottom() const;
  bool is_top() const;

  // Number of references that may point to the region.
  small_range &refcount_val();
  // Whether the region has been initialized.
  boolean_value &init_val();
  // Type of the region. The type of a region is always known
  // statically except if its type is unknown.
  type_value &type_val();
  const small_range &refcount_val() const;
  const boolean_value &init_val() const;
  const type_value &type_val() const;

  bool operator<=(const region_info &other) const;
  bool operator==(const region_info &other) const;
  void operator|=(const region_info &other);
  region_info operator|(const region_info &other) const;
  region_info operator||(const region_info &other) const;
  region_info operator&(const region_info &other) const;
  region_info operator&&(const region_info &other) const;
  void write(crab::crab_os &o) const;
  friend crab::crab_os &operator<<(crab::crab_os &o, region_info &dom) {
    dom.write(o);
    return o;
  }
};
} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
