#pragma once
#include "./crab_lang.hpp"

// To run abstract domains defined over integers
template<typename BUDom, typename TDDom>
extern void bu_inter_run(crab::cg_impl::z_cg_t* cg,
			 BUDom bu_top, TDDom td_top,
			 bool run_liveness,
			 unsigned widening, 
			 unsigned narrowing, 
			 unsigned jump_set_size,
			 bool enable_stats);
   

