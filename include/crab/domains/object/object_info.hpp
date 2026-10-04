#pragma once

#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/small_range.hpp>

namespace crab {
namespace domains {
namespace object_domain_impl {

/* This class contains domain reduction between obj's cache and base domains */
class reduction_info {
  // Simple pair of (bool x bool)
  using product_t = std::pair<bool, bool>;

  product_t m_product;

public:
  reduction_info() { m_product = {false, false}; }
  reduction_info(const bool &is_loaded, const bool &is_stored) {
    cache_reg_loaded_val() = is_loaded;
    cache_reg_stored_val() = is_stored;
  }
  reduction_info(const reduction_info &other) = default;
  reduction_info(reduction_info &&other) = default;
  reduction_info &operator=(const reduction_info &other) = default;
  reduction_info &operator=(reduction_info &&other) = default;

  bool &cache_reg_loaded_val() { return m_product.first; }
  bool &cache_reg_stored_val() { return m_product.second; }
  const bool &cache_reg_loaded_val() const { return m_product.first; }
  const bool &cache_reg_stored_val() const { return m_product.second; }

  void write(crab::crab_os &o) const {
    auto print_bool = [&o](const bool &b) {
      o << (b ? "true" : "false");
    };
    o << "CachLoadedBYReg=";
    print_bool(cache_reg_loaded_val());
    o << ",CacheStoredFromReg=";
    print_bool(cache_reg_stored_val());
  }
};

/* This class contains any non-relational information about abstract objects */
class object_info {
  // Product (no reduction):
  //   <refcount, <obj-init, sum-present>> x <cache-used, cache-dirty>
  using cache_info_t = basic_domain_product2<boolean_value, boolean_value>;
  using status_info_t = basic_domain_product2<boolean_value, boolean_value>;
  using object_info_t = basic_domain_product2<small_range, status_info_t>;
  using product_t = basic_domain_product2<object_info_t, cache_info_t>;

  using reduction_flag_t = reduction_info;

  product_t m_product;
  reduction_flag_t m_flags;

  // NOTE: drops m_flags; only reachable from the lattice operators below,
  // which reset the flags to false anyway (a merged cache is never
  // reg-loaded/stored).
  object_info(product_t &&product) : m_product(std::move(product)) {}

public:
  object_info() { m_product.set_to_top(); }
  object_info(const small_range &count, const boolean_value &obj_init,
              const boolean_value &no_sum, const boolean_value &cache_used,
              const boolean_value &cache_dirty, const bool &is_loaded,
              const bool &is_stored) {
    refcount_val() = count;
    objinit_val() = obj_init;
    sumpresence_val() = no_sum;
    cacheused_val() = cache_used;
    cachedirty_val() = cache_dirty;
    m_flags = reduction_flag_t(is_loaded, is_stored);
  }
  object_info(const object_info &other) = default;
  object_info(object_info &&other) = default;
  object_info &operator=(const object_info &other) = default;
  object_info &operator=(object_info &&other) = default;

  void set_to_top() { m_product.set_to_top(); }
  void set_to_bottom() { m_product.set_to_bottom(); }
  bool is_bottom() const { return m_product.is_bottom(); }
  bool is_top() const { return m_product.is_top(); }

  // Number of references that may point to an object.
  small_range &refcount_val() { return m_product.first().first(); }
  // Whether the object is inited through the store_ref
  boolean_value &objinit_val() { return m_product.first().second().first(); }
  // Object summary is presence.
  boolean_value &sumpresence_val() {
    return m_product.first().second().second();
  }
  // Whether the cache domain is used.
  boolean_value &cacheused_val() { return m_product.second().first(); }
  // Whether the cache domain has been updated.
  boolean_value &cachedirty_val() {
    return m_product.second().second();
  }
  bool &cache_reg_loaded_val() {
    return m_flags.cache_reg_loaded_val();
  }
  bool &cache_reg_stored_val() {
    return m_flags.cache_reg_stored_val();
  }
  const small_range &refcount_val() const { return m_product.first().first(); }
  const boolean_value &objinit_val() const {
    return m_product.first().second().first();
  }
  const boolean_value &sumpresence_val() const {
    return m_product.first().second().second();
  }
  const boolean_value &cacheused_val() const {
    return m_product.second().first();
  }
  const boolean_value &cachedirty_val() const {
    return m_product.second().second();
  }
  const bool &cache_reg_loaded_val() const {
    return m_flags.cache_reg_loaded_val();
  }
  const bool &cache_reg_stored_val() const {
    return m_flags.cache_reg_stored_val();
  }

  /// @brief exact equality of the whole info tuple, including the two
  /// reduction flags (unlike the old operator==, which compared only the
  /// refcount). Used for O(1) leaf-equality in the odi map.
  bool equals(const object_info &o) const {
    return refcount_val() == o.refcount_val() &&
           objinit_val() == o.objinit_val() &&
           sumpresence_val() == o.sumpresence_val() &&
           cacheused_val() == o.cacheused_val() &&
           cachedirty_val() == o.cachedirty_val() &&
           cache_reg_loaded_val() == o.cache_reg_loaded_val() &&
           cache_reg_stored_val() == o.cache_reg_stored_val();
  }

  // The lattice interface below is required to instantiate
  // basic_domain_product2<object_info, ...>: the product's lattice operators
  // are virtual overrides, so they are instantiated with the class whether or
  // not the object domain ever calls them (the odi map combines infos
  // explicitly in its join/meet operators instead).
  bool operator<=(const object_info &other) const {
    return refcount_val() <= other.refcount_val();
  }
  void operator|=(const object_info &other) { m_product |= other.m_product; }
  object_info operator|(const object_info &other) const {
    return object_info(m_product | other.m_product);
  }
  object_info operator||(const object_info &other) const {
    return object_info(m_product || other.m_product);
  }
  object_info operator&(const object_info &other) const {
    return object_info(m_product & other.m_product);
  }
  object_info operator&&(const object_info &other) const {
    return object_info(m_product && other.m_product);
  }
  std::string domain_name(void) const { return "Object Info"; }

  void write(crab::crab_os &o) const {

    o << "RefCount=" << refcount_val() << ",";
    CRAB_LOG("object-print", o << "ObjInit=" << objinit_val() << ","
                               << "SumPresence=" << sumpresence_val() << ",");
    o << "CacheUsed=" << cacheused_val() << ","
      << "CacheDirty=" << cachedirty_val() << ",";
    m_flags.write(o);
  }
  friend crab::crab_os &operator<<(crab::crab_os &o, const object_info &dom) {
    dom.write(o);
    return o;
  }
};
} // namespace object_domain_impl
} // end namespace domains
} // end namespace crab
