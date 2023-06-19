#pragma once

namespace crab {

/** Class to group together all the fixpoint parameters **/
class fixpoint_parameters {
  // number of iterations until widening operator is called.
  unsigned widening_delay;
  // number of descending iterations. This is needed in case the
  // abstract domain does not implement a narrowing operator.
  unsigned descending_iterations;
  // Used to implement the technique of widening with thresholds.
  // If 0 then widening with thresholds is disabled
  // 
  // The set of thresholds is currently computed by a fixed strategy
  // that collects at most max_thresholds per wto cycle. This means
  // that each wto cycle has its own set of thresholds.  Note that the
  // bigger is this number the slower can be the fixpoint convergence.
  unsigned max_thresholds;
  // Enable decoupling of ascending/descending fixpoint phases
  bool enabled_decoupling;
  
public:
  
  fixpoint_parameters():
    widening_delay(2),
    descending_iterations(1),
    max_thresholds(0),
    enabled_decoupling(false) {}

  unsigned get_widening_delay() const { return widening_delay; }
  unsigned& get_widening_delay() { return widening_delay; }  

  unsigned get_descending_iterations() const { return descending_iterations; }
  unsigned& get_descending_iterations() { return descending_iterations; }  

  unsigned get_max_thresholds() const { return max_thresholds; }
  unsigned& get_max_thresholds() { return max_thresholds; }

  bool is_enabled_decoupling() const { return enabled_decoupling; }
  bool& enable_decoupling() { return enabled_decoupling; }
  
};
  
} // end namespace crab 
    
