#pragma once

#include "./crab_lang.hpp"

template <typename Dom>
extern void z_intra_run(crab::cfg_impl::z_cfg_t *cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, unsigned widening,
                        unsigned narrowing, unsigned jump_set_size,
                        bool enable_stats, bool enable_checker, bool print_invariants);

template <typename Dom>
extern void q_intra_run(crab::cfg_impl::q_cfg_t *cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, unsigned widening,
                        unsigned narrowing, unsigned jump_set_size,
                        bool enable_stats, bool enable_checker, bool print_invariants);

// To run abstract domains defined over integers
template <typename Dom>
void run(crab::cfg_impl::z_cfg_t *cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  z_abs_domain_t init_wrapper(init);
  z_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#else
  z_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#endif
}

template <typename Dom>
void run_and_check(crab::cfg_impl::z_cfg_t *cfg,
                   crab::cfg_impl::basic_block_label_t entry, Dom init,
                   bool run_liveness, unsigned widening, unsigned narrowing,
                   unsigned jump_set_size, bool enable_stats) {

#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  z_abs_domain_t init_wrapper(init);
  z_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, true, false);

#else
  z_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, true, false);

#endif
}

// To run abstract domains defined over rationals
template <typename Dom>
void run(crab::cfg_impl::q_cfg_t *cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  q_abs_domain_t init_wrapper(init);
  q_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#else
  q_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#endif
}
