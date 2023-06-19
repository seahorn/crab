#pragma once

#include "./crab_lang.hpp"
#include <crab/fixpoint/fixpoint_params.hpp>

template <typename Dom>
extern void z_intra_run(crab::cfg_impl::z_cfg_t *cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, const crab::fixpoint_parameters &params,
                        bool enable_stats, bool enable_checker, bool print_invariants);

template <typename Dom>
extern void q_intra_run(crab::cfg_impl::q_cfg_t *cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, const crab::fixpoint_parameters &params,
                        bool enable_stats, bool enable_checker, bool print_invariants);

#ifdef USE_GENERIC_WRAPPER
template<typename Dom>
crab::domain_impl::z_abs_domain_t make_z_absval(Dom val) {
  return crab::domain_impl::z_abs_domain_t(val);
}
template<typename Dom>
crab::domain_impl::q_abs_domain_t make_q_absval(Dom val) {
  return crab::domain_impl::q_abs_domain_t(val);
}
#else
template<typename Dom>
Dom make_z_absval(Dom val) {
  return val;
}
template<typename Dom>
Dom make_q_absval(Dom val) {
  return val;
}
#endif


// To run abstract domains defined over integers
template <typename Dom>
void run(crab::cfg_impl::z_cfg_t *cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening;
  fixpo_params.get_descending_iterations() = narrowing;
  fixpo_params.get_max_thresholds() = jump_set_size; 
  z_intra_run(cfg, entry, make_z_absval(init), run_liveness, fixpo_params, enable_stats, false, true);
}

template <typename Dom>
void run(crab::cfg_impl::z_cfg_t *cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         const crab::fixpoint_parameters &fixpo_params, bool enable_stats) {
  z_intra_run(cfg, entry, make_z_absval(init), run_liveness, fixpo_params, enable_stats, false, true);
}


template <typename Dom>
void run_and_check(crab::cfg_impl::z_cfg_t *cfg,
                   crab::cfg_impl::basic_block_label_t entry, Dom init,
                   bool run_liveness, unsigned widening, unsigned narrowing,
                   unsigned jump_set_size, bool enable_stats) {
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening;
  fixpo_params.get_descending_iterations() = narrowing;
  fixpo_params.get_max_thresholds() = jump_set_size;  
  z_intra_run(cfg, entry, make_z_absval(init), run_liveness, fixpo_params, enable_stats, true, false);
}


template <typename Dom>
void run_and_check(crab::cfg_impl::z_cfg_t *cfg,
                   crab::cfg_impl::basic_block_label_t entry, Dom init,
                   bool run_liveness, const crab::fixpoint_parameters &fixpo_params, bool enable_stats) {
  z_intra_run(cfg, entry, make_z_absval(init), run_liveness, fixpo_params, enable_stats, true, false);
}

// To run abstract domains defined over rationals
template <typename Dom>
void run(crab::cfg_impl::q_cfg_t *cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening;
  fixpo_params.get_descending_iterations() = narrowing;
  fixpo_params.get_max_thresholds() = jump_set_size;  
  q_intra_run(cfg, entry, make_q_absval(init), run_liveness, fixpo_params, enable_stats, false, true);
}
