#ifndef __CRAB_BACKWARD_ANALYZERS__
#define __CRAB_BACKWARD_ANALYZERS__

#include "./crab_lang.hpp"

// To run abstract domains defined over integers
template<typename Dom>
extern void backward_run (crab::cfg_impl::z_cfg_t* cfg,
			  crab::cfg_impl::basic_block_label_t entry,
			  Dom initial_states,
			  Dom final_states,
			  unsigned widening, 
			  unsigned narrowing, 
			  unsigned jump_set_size,
			  bool enable_stats);

#endif  

