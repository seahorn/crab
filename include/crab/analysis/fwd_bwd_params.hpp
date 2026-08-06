#pragma once

namespace crab {
namespace analyzer {
  
/**
 *  Class to group together all the parameters of the iterative
 *  forward-backward analyzer
 **/
class fwd_bwd_parameters {
  // If false, then it behaves like a standard forward analyzer
  bool enabled_backward;
  // Maximum number of refinement iterations (if enabled_backward=true)
  //
  // IMPORTANT: this bound is what guarantees termination of the
  // forward-backward refinement loop. The loop combines the old
  // assumptions with the new backward preconditions using meet (not
  // narrowing), so the sequence of assumptions is strictly decreasing
  // but need not stabilize: abstract domains such as intervals have
  // infinite strictly-decreasing chains. Keep this number small. Apart
  // from termination, each extra iteration costs a full forward plus a
  // full backward fixpoint, and in practice refinement rarely makes
  // progress after a handful of iterations.
  unsigned max_refine_iterations;
  // If disabled the iterative forward-backward analyzer produces the
  // same invariants than a classical forward analyzer. The difference
  // is that it might prove more assertions.
  bool use_refined_invariants;
  
public:
  
  fwd_bwd_parameters():
    enabled_backward(false),
    max_refine_iterations(5),
    use_refined_invariants(false) {}

  bool is_enabled_backward() const { return enabled_backward; }
  bool& enable_backward() { return enabled_backward; }  

  unsigned get_max_refine_iterations() const { return max_refine_iterations; }
  unsigned& get_max_refine_iterations() { return max_refine_iterations; }  

  bool get_use_refined_invariants() const { return use_refined_invariants; }
  bool& get_use_refined_invariants() { return use_refined_invariants; }  
};

} // end namespace analyzer  
} // end namespace crab 
