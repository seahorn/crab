#pragma once

#include <string>

namespace crab {

template <typename BasicBlock> class basic_block_traits {
public:
  using basic_block_label_t = typename BasicBlock::basic_block_label_t;
  static std::string to_string(const basic_block_label_t &bb);
};

} // namespace crab
