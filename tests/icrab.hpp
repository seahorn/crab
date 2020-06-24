#pragma once
#include "./crab_lang.hpp"

// For top-down inter-procedural analysis parameters
#include <crab/analysis/inter/top_down_inter_analyzer.hpp>


// To run abstract domains defined over integers
template<typename BUDom, typename TDDom>
extern void bu_inter_run(crab::cg_impl::z_cg_t* cg,
			 BUDom bu_top, TDDom td_top,
			 bool run_liveness,
			 unsigned widening, 
			 unsigned narrowing, 
			 unsigned jump_set_size,
			 bool enable_stats);

// To run abstract domains defined over integers
typedef crab::analyzer::top_down_inter_analyzer_parameters
<crab::cg_impl::z_cg_ref_t> td_inter_params_t;

template<typename Dom>
extern void td_inter_run(crab::cg_impl::z_cg_t* cg,
			 Dom init,
			 td_inter_params_t params,
			 bool print_checks,
			 bool print_invariants,
			 bool enable_stats);


