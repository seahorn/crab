#include <crab/support/os.hpp>
#include <crab/types/allocation_sites.hpp>

namespace crab {

allocation_site::allocation_site(size_t id) : m_id(id) {}

bool allocation_site::operator<(const allocation_site &as) const {
  return m_id < as.m_id;
}

bool allocation_site::operator==(const allocation_site &as) const {
  return m_id == as.m_id;
}

ikos::index_t allocation_site::index() const { return m_id; }

void allocation_site::write(crab_os &o) const { o << "as_" << m_id; }

allocation_site_man::allocation_site_man() : m_id(0) {}

allocation_site allocation_site_man::mk_allocation_site() {
  allocation_site as(m_id++);
  return as;
}

size_t allocation_site_man::size() const { return m_id; }

} // end namespace crab
