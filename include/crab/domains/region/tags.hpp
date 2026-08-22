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

// A class with more information about the source of a tag.
template <typename Number, typename Location>
class tag_with_info : public indexable {
  ikos::index_t m_id; // Tag id
  Location m_src_loc; // location that generated the tag
  using this_class_t = tag_with_info<Number, Location>;
  void check_id(Number n) {
    if (n < 0) {
      CRAB_ERROR("Cannot use negative numbers for tags");
    }
    if (!n.fits_int64()) {
      CRAB_ERROR("Too large value for a tag");
    }
  }

public:
  tag_with_info(Number n) : m_id(0), m_src_loc(Location()) {
    check_id(n);
    m_id = (int64_t)n;
  }
  tag_with_info(Number n, Location loc) : m_id(0), m_src_loc(loc) {
    check_id(n);
    m_id = (int64_t)n;
  }
  // Lexicographic comparison over (m_id, is_valid, location id) so
  // that operator< is a strict weak order. An invalid location
  // carries no meaningful id, so it is ordered before any valid one
  // and compared as id-only otherwise.
  bool operator<(const this_class_t &as) const {
    if (m_id != as.m_id) {
      return m_id < as.m_id;
    }
    const bool this_valid = m_src_loc.is_valid();
    const bool other_valid = as.m_src_loc.is_valid();
    if (this_valid != other_valid) {
      return !this_valid;
    }
    if (!this_valid) {
      // both locations are invalid: compare by tag id only
      return false;
    }
    return m_src_loc.get_id() < as.m_src_loc.get_id();
  }
  bool operator==(const this_class_t &as) const {
    if (m_id != as.m_id) {
      return false;
    }
    const bool this_valid = m_src_loc.is_valid();
    const bool other_valid = as.m_src_loc.is_valid();
    if (this_valid != other_valid) {
      return false;
    }
    // both invalid locations compare equal (id-only comparison)
    return !this_valid || (m_src_loc.get_id() == as.m_src_loc.get_id());
  }
  virtual ikos::index_t index() const override { return m_id; }
  void write(crab_os &o) const override {
    o << "TAG_" << m_id;
    if (m_src_loc.is_valid()) {
      o << " [src: " << m_src_loc << "]";
    } else {
      o << " [src: <unknown>]";
    }
  }
  friend crab_os &operator<<(crab_os &o, const this_class_t &as) {
    as.write(o);
    return o;
  }
}; /* end class tag_with_info */
} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
