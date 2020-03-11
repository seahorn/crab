#pragma once

#include <string>

namespace crab {
namespace transforms {

/**
 * API for Crab Control-Flow Graph Transformations.
 **/
template <typename CFG> class transform {
public:
  virtual bool run(CFG &cfg) = 0;
  virtual std::string get_name() const = 0;
};

} // namespace transforms
} // namespace crab
