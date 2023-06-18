#pragma once

namespace crab {

/** Class to group together all the fixpoint parameters **/
class fixpoint_parameters {
  // number of iterations until widening operator is called.
  unsigned widening_delay;
  // number of descending iterations. This is needed in case the
  // abstract domain does not implement a narrowing operator.
  unsigned descending_iterations;
  // Used to implement the technique of widening with thresholds.  The
  // bigger is this number the slower can be the fixpoint convergence.
  // 
  // The set of thresholds is currently computed by a fixed strategy
  // that collects at most max_thresholds per wto cycle. This
  // means that each wto cycle has its own set of thresholds.
  unsigned max_thresholds;

public:
  
  fixpoint_parameters():
    widening_delay(2),
    descending_iterations(1),
    // Set to 0 to disable widening with thresholds
    max_thresholds(0) {}

  unsigned get_widening_delay() const { return widening_delay; }
  unsigned& get_widening_delay() { return widening_delay; }  

  unsigned get_descending_iterations() const { return descending_iterations; }
  unsigned& get_descending_iterations() { return descending_iterations; }  

  unsigned get_max_thresholds() const { return max_thresholds; }
  unsigned& get_max_thresholds() { return max_thresholds; }  
};
  
} // end namespace crab 
    
