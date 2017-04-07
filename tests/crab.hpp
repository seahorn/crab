#ifndef __CRAB_INTRA_ANALYZERS__
#define __CRAB_INTRA_ANALYZERS__

#include "./crab_lang.hpp"

// To run abstract domains defined over integers
template<typename Dom>
extern void run (crab::cfg_impl::z_cfg_t* cfg, 
		 bool run_liveness,
		 unsigned widening, 
		 unsigned narrowing, 
		 unsigned jump_set_size,
		 bool enable_stats);

// To run abstract domains defined over rationals
template<typename Dom>
extern void run (crab::cfg_impl::q_cfg_t* cfg, 
		 bool run_liveness,
		 unsigned widening, 
		 unsigned narrowing, 
		 unsigned jump_set_size,
		 bool enable_stats);
  
   
#endif  /*__CRAB_INTRA_ANALYZERS__*/

