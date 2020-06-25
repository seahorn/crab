#pragma once

#include <cstdint>

namespace ikos {
// Numerical type for indexed objects used by patricia trees
using index_t = uint64_t;
} // end namespace ikos

namespace crab {
class indexable {
public:
  virtual ~indexable() = default;
  virtual ikos::index_t index() const = 0;
};
} // end namespace crab
