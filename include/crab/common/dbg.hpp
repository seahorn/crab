#ifndef __DBG_HPP__
#define __DBG_HPP__

#include <crab/common/types.hpp>

#undef CRAB_DEBUG

#define CRAB_DEBUG(...)              \
    do {                             \
      ___print___(__VA_ARGS__);      \
      std::cerr << "\n";             \
    } while (0)

#endif 
