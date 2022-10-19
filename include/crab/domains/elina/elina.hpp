#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "elina_abstract0.h"
#include "elina_coeff.h"
#include "elina_dimension.h"
#include "elina_interval.h"
#include "elina_lincons0.h"
#include "elina_linexpr0.h"
#include "elina_manager.h"
#include "elina_scalar.h"
#include "elina_tcons0.h"
#include "elina_texpr0.h"
#include "opt_oct.h"
#include "opt_pk.h"
#include "opt_zones.h"
#pragma GCC diagnostic pop

#include <boost/shared_ptr.hpp>

namespace crab {
namespace domains {
namespace elina {

using elina_state = elina_abstract0_t *;
using elina_state_ptr = boost::shared_ptr<elina_abstract0_t>;

/** Custom deleter for elina_state_ptr */
class elina_state_deleter {
private:
  elina_manager_t *man;

public:
  elina_state_deleter(elina_manager_t *_m) : man(_m) {}

  void operator()(elina_state s) {
    if (s)
      elina_abstract0_free(man, s);
  }
};

/** Create a elina_state_ptr from elina_state */
inline elina_state_ptr elinaPtr(elina_manager_t *_m, elina_state _s) {
  if (_s) {
    elina_state_ptr p(_s, elina_state_deleter(_m));
    return p;
  }
  CRAB_ERROR("elinaPtr is taking a null pointer!");
}

} // namespace elina
} // namespace domains
} // namespace crab
