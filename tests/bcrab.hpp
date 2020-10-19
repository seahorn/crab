#ifndef __CRAB_BACKWARD_ANALYZERS__
#define __CRAB_BACKWARD_ANALYZERS__

#include "./crab_lang.hpp"
#include "./crab_dom.hpp"

template <typename Dom>
extern void z_backward_run(crab::cfg_impl::z_cfg_t *cfg,
                           crab::cfg_impl::basic_block_label_t entry,
                           Dom initial_states, unsigned widening,
                           unsigned narrowing, unsigned jump_set_size,
                           bool enable_stats);

template <typename Dom>
void backward_run(crab::cfg_impl::z_cfg_t *cfg,
                  crab::cfg_impl::basic_block_label_t entry, Dom initial_states,
                  unsigned widening, unsigned narrowing, unsigned jump_set_size,
                  bool enable_stats) {
#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  z_abs_domain_t init(initial_states);
  z_backward_run(cfg, entry, init, widening, narrowing, jump_set_size,
                 enable_stats);
#else
  z_backward_run(cfg, entry, initial_states, widening, narrowing, jump_set_size,
                 enable_stats);
#endif
}

#endif
