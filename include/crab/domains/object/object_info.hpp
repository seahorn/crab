#pragma once

#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/term/term_operators.hpp>

namespace crab {
namespace domains {
namespace object_domain_impl {

/* This class contains domain reduction between obj's cache and base domains */
class reduction_info {
  // Simple pair of (bool x bool)
  using product_t = std::pair<bool, bool>;

  product_t m_product;
  reduction_info(product_t &&product) : m_product(std::move(product)) {}

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
    o << ","<< "CacheStoredFromReg=";
    print_bool(cache_reg_stored_val());
  }
  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const reduction_info &dom) {
    dom.write(o);
    return o;
  }
  std::string domain_name(void) const { return "Reduction Info"; }
};

/* This class contains any non-relational information about abstract objects */
class object_info {
  // Simple reduced product of (SmallRange x Boolean x Boolean)
  using product_t = basic_domain_product2<
      small_range, basic_domain_product2<boolean_value, boolean_value>>;

  using reduction_flag_t = reduction_info;

  product_t m_product;
  reduction_flag_t m_flags;
  object_info(product_t &&product) : m_product(std::move(product)) {}

public:
  object_info() { m_product.set_to_top(); }
  object_info(const small_range &count, const boolean_value &cache_used,
              const boolean_value &cache_dirty, const bool &is_loaded,
              const bool &is_stored) {
    refcount_val() = count;
    cacheused_val() = cache_used;
    cachedirty_val() = cache_dirty;
    m_flags = reduction_flag_t(is_loaded, is_stored);
  }
  object_info(const object_info &other) = default;
  object_info(object_info &&other) = default;
  object_info &operator=(const object_info &other) = default;
  object_info &operator=(object_info &&other) = default;

  static object_info bottom() {
    object_info res;
    res.m_product.set_to_bottom();
    return res;
  }
  static object_info top() {
    object_info res;
    return res;
  }
  void set_to_top() { m_product.set_to_top(); }
  void set_to_bottom() { m_product.set_to_bottom(); }
  bool is_bottom() const { return m_product.is_bottom(); }
  bool is_top() const { return m_product.is_top(); }

  // Number of references that may point to an object.
  small_range &refcount_val() { return m_product.first(); }
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
  const small_range &refcount_val() const { return m_product.first(); }
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
    return m_flags.cache_reg_loaded_val();
  }

  bool operator<=(const object_info &other) const {
    return refcount_val() <= other.refcount_val();
  }
  bool operator==(const object_info &other) const {
    return (*this <= other && other <= *this);
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
  void write(crab::crab_os &o) const {
    o << "RefCount=" << refcount_val() << ","
      << "CacheUsed=" << cacheused_val() << ","
      << "CacheDirty=" << cachedirty_val() << ",";
    m_flags.write(o);
  }
  friend crab::crab_os &operator<<(crab::crab_os &o, const object_info &dom) {
    dom.write(o);
    return o;
  }
  std::string domain_name(void) const { return "Object Info"; }
};

class id_val_generator_t {
private:
  static uint32_t m_next_val;

public:
  static uint32_t get_next_val() { return m_next_val++; }
};
} // namespace object_domain_impl
} // end namespace domains
} // end namespace crab
