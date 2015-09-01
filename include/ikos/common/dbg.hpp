#ifndef __DBG_HPP__
#define __DBG_HPP__
#include <ikos/common/types.hpp>

namespace ikos {

#undef IKOS_DEBUG

#define IKOS_DEBUG(...)              \
    do {                             \
      ___print___(__VA_ARGS__);      \
      std::cerr << "\n";             \
    } while (0)

} // end namespace

#endif 
