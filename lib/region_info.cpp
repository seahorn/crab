#include <crab/domains/region/region_info.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace domains {
namespace region_domain_impl {

region_info::region_info() { m_product.set_to_top(); }

region_info::region_info(region_info::product_t &&product)
    : m_product(std::move(product)) {}

region_info::region_info(small_range count, boolean_value init,
                         type_value type) {
  refcount_val() = count;
  init_val() = init;
  type_val() = type;
}

region_info region_info::bottom() {
  region_info res;
  res.m_product.set_to_bottom();
  return res;
}

region_info region_info::top() {
  region_info res;
  return res;
}

bool region_info::is_bottom() const { return m_product.is_bottom(); }

bool region_info::is_top() const { return m_product.is_top(); }

small_range &region_info::refcount_val() { return m_product.first(); }

boolean_value &region_info::init_val() { return m_product.second().first(); }

type_value &region_info::type_val() { return m_product.second().second(); }

const small_range &region_info::refcount_val() const { return m_product.first(); }

const boolean_value &region_info::init_val() const {
  return m_product.second().first();
}

const type_value &region_info::type_val() const {
  return m_product.second().second();
}

bool region_info::operator<=(const region_info &other) const {
  return m_product <= other.m_product;
}

bool region_info::operator==(const region_info &other) const {
  return (operator<=(other) && other.operator<=(*this));
}

void region_info::operator|=(const region_info &other) {
  m_product |= other.m_product;
}

region_info region_info::operator|(const region_info &other) const {
  return region_info(m_product | other.m_product);
}

region_info region_info::operator||(const region_info &other) const {
  return region_info(m_product || other.m_product);
}

region_info region_info::operator&(const region_info &other) const {
  return region_info(m_product & other.m_product);
}

region_info region_info::operator&&(const region_info &other) const {
  return region_info(m_product && other.m_product);
}

void region_info::write(crab::crab_os &o) const {
  o << "RefCount=" << refcount_val() << ","
    << "Init=" << init_val() << ","
    << "DynType=" << type_val();
}
} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
