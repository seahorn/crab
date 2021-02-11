#pragma once

#include <crab/types/indexable.hpp>
#include <cstddef>

/* A simple class to represent an allocation site */

namespace crab {
class crab_os;
} // end namespace crab

namespace crab {
class allocation_site_man;

class allocation_site : public indexable {
  std::size_t m_id;
  explicit allocation_site(std::size_t id);
  friend class allocation_site_man;

public:
  allocation_site(const allocation_site &as) = default;
  allocation_site(allocation_site &&as) = default;
  allocation_site &operator=(const allocation_site &as) = default;
  allocation_site &operator=(allocation_site &&as) = default;
  bool operator<(const allocation_site &as) const;
  bool operator==(const allocation_site &as) const;
  virtual ikos::index_t index() const override;
  void write(crab_os &o) const override;
  friend crab_os &operator<<(crab_os &o, const allocation_site &as) {
    as.write(o);
    return o;
  }
};

class allocation_site_man {
  std::size_t m_id;

public:
  allocation_site_man();
  // create an allocation site
  allocation_site mk_allocation_site();
  // return the number of allocation sites
  std::size_t size() const;
};

} // namespace crab
