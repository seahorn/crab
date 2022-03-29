#include <crab/domains/interval_congruence_impl.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {
namespace domains {  
// Default instantiations
template class interval_congruence<ikos::z_number>;
} // end namespace domains
} // end namespace crab
