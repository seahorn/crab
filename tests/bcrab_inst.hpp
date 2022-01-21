#pragma once
// Here all explicit instantiations

#define Z_BWD_RUNNER(DOM)                                                      \
  template void z_backward_run<DOM>(crab::cfg_impl::z_cfg_t *,                 \
                                    crab::cfg_impl::basic_block_label_t, DOM,  \
                                    unsigned, unsigned, unsigned, bool);

#ifdef USE_GENERIC_WRAPPER
Z_BWD_RUNNER(crab::domain_impl::z_abs_domain_t)
#else
Z_BWD_RUNNER(crab::domain_impl::z_interval_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_boxes_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_box_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_pk_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_oct_elina_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_aa_int_t)
Z_BWD_RUNNER(crab::domain_impl::z_aa_box_apron_t)
Z_BWD_RUNNER(crab::domain_impl::z_aa_zones_elina_t)
#endif
