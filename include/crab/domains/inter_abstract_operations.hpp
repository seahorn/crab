#pragma once

#include <crab/domains/inter_abstract_operations_callsite_info.hpp>
#include <string>

namespace crab {
namespace domains {

template <class Domain, bool UseDefaultImpl> class inter_abstract_operations {
public:
  using variable_t = typename Domain::variable_t;

  /**
   * Upon completion, @callee is the abstract state at the callee for
   * the entry block.
   *
   *  - @caller: abstract state at the caller before the call.
   *  - @callee: initial state at the callee (usually top)
   **/
  static void callee_entry(const callsite_info<variable_t> &callsite,
                           const Domain &caller, Domain &callee);

  /**
   * Upon completion, @caller is the abstract state at the
   * caller after the execution of the call finished.
   *
   * - @caller: abstract state at the caller before the call
   * - @callee: abstract state at the exit block of the callee
   *   but already projected onto formal parameters of the function.
   *
   * This code should work even if input and output parameters at the
   * callsite are not disjoint.
   **/
  static void caller_continuation(const callsite_info<variable_t> &callsite,
                                  const Domain &callee, Domain &caller);
};

} // end namespace domains
} // end namespace crab

#include "inter_abstract_operations_impl.hpp"
