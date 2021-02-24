#pragma once

#include <crab/types/indexable.hpp>
#include <cstddef>

/* A simple class to represent tags (i.e., labels) */

namespace crab {
class crab_os;
} // end namespace crab

namespace crab {
class tag_manager;

class tag : public indexable {
  std::size_t m_id;
  explicit tag(std::size_t id);
  friend class tag_manager;

public:
  tag(const tag &as) = default;
  tag(tag &&as) = default;
  tag &operator=(const tag &as) = default;
  tag &operator=(tag &&as) = default;
  bool operator<(const tag &as) const;
  bool operator==(const tag &as) const;
  virtual ikos::index_t index() const override;
  void write(crab_os &o) const override;
  friend crab_os &operator<<(crab_os &o, const tag &as) {
    as.write(o);
    return o;
  }
};

class tag_manager {
  std::size_t m_id;

public:
  tag_manager();
  // create a tag
  tag mk_tag();
  // return the number of tags
  std::size_t size() const;
};

} // namespace crab
