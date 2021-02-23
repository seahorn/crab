#include <crab/support/os.hpp>
#include <crab/types/tag.hpp>

namespace crab {

tag::tag(size_t id) : m_id(id) {}

bool tag::operator<(const tag &as) const {
  return m_id < as.m_id;
}

bool tag::operator==(const tag &as) const {
  return m_id == as.m_id;
}

ikos::index_t tag::index() const { return m_id; }

void tag::write(crab_os &o) const { o << "as_" << m_id; }

tag_manager::tag_manager() : m_id(0) {}

tag tag_manager::mk_tag() {
  tag as(m_id++);
  return as;
}

size_t tag_manager::size() const { return m_id; }

} // end namespace crab
