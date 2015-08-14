#ifndef __DBG_HPP__
#define __DBG_HPP__
#include <ikos/common/types.hpp>

namespace ikos {
#ifdef _IKOS_DEBUG_
#define IKOS_DEBUG(...)              \
    do {                             \
      ___print___(__VA_ARGS__);      \
      std::cerr << "\n";             \
    } while (0)
#else
#define IKOS_DEBUG(...)
#endif
} // end namespace

#endif 
