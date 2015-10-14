#include <crab/config.h>

#ifdef HAVE_LDD
#include <crab/domains/boxes.hpp>
using namespace crab::domains::ldd;
namespace crab { 
   namespace domains {
     LddManager* LddManagerWrapper::m_ldd_man = nullptr;
   } 
}
#endif 
