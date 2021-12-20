#pragma once

#include <crab/domains/types.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace domains {
namespace region_domain_impl {
// A simple class for tags (i.e., numerical identifiers). We don't
// use crab::tag because we want to have the flexibility of creating
// tags without a tag manager.
template <typename Number> class tag : public indexable {
  ikos::index_t m_id;

public:
  tag(Number n) : m_id(0) {
    if (n < 0) {
      CRAB_ERROR("Cannot use negative numbers for tags");
    }
    if (!n.fits_int64()) {
      CRAB_ERROR("Too large value for a tag");
    }
    m_id = (int64_t)n;
  }
  bool operator<(const tag &as) const { return m_id < as.m_id; }
  bool operator==(const tag &as) const { return m_id == as.m_id; }
  virtual ikos::index_t index() const override { return m_id; }
  void write(crab_os &o) const override { o << "TAG_" << m_id; }
  friend crab_os &operator<<(crab_os &o, const tag &as) {
    as.write(o);
    return o;
  }
}; /* end class tag */
} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
