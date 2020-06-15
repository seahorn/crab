#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <functional>
#include <iosfwd>
#include <memory>

/* Basic type definitions */

namespace ikos {
// Numerical type for indexed objects
typedef uint64_t index_t;
} // end namespace ikos

namespace crab {
namespace domains {
using namespace ikos;
}
} // namespace crab


