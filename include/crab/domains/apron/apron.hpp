#pragma once

#include "ap_global0.h"
#include "box.h"
#include "oct.h"
#include "pk.h"
#include "pkeq.h"

#ifdef HAVE_PPLITE
#include "ap_pplite.h"
#endif

#include <memory>

namespace crab {

namespace domains {

namespace apron {

using ap_state = ap_abstract0_t *;
using ap_state_ptr = std::shared_ptr<ap_abstract0_t>;

/** Custom deleter for ap_state_ptr */
class ap_state_deleter {
private:
  ap_manager_t *man;

public:
  ap_state_deleter(ap_manager_t *_m) : man(_m) {}

  void operator()(ap_state s) {
    if (s)
      ap_abstract0_free(man, s);
  }

  // ap_manager_t** getManager () const { return man; }
};

/** Create a ap_state_ptr from ap_state */
inline ap_state_ptr apPtr(ap_manager_t *_m, ap_state _s) {
  if (_s) {
    ap_state_ptr p(_s, ap_state_deleter(_m));
    return p;
  }
  CRAB_ERROR("apPtr is taking a null pointer!");
  // return ap_state_ptr();
}
} // namespace apron
} // namespace domains
} // namespace crab
