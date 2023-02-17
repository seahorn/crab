#pragma once

#include "./crab.hpp"

// Here all explicit instantiations

#define Z_RUNNER(DOM)                                                          \
  template void z_intra_run<DOM>(                                              \
      crab::cfg_impl::z_cfg_t *, crab::cfg_impl::basic_block_label_t, DOM,     \
      bool, unsigned, unsigned, unsigned, bool, bool, bool);

#define Q_RUNNER(DOM)                                                          \
  template void q_intra_run<DOM>(                                              \
      crab::cfg_impl::q_cfg_t *, crab::cfg_impl::basic_block_label_t, DOM,     \
      bool, unsigned, unsigned, unsigned, bool, bool, bool);

#ifdef USE_GENERIC_WRAPPER
Z_RUNNER(crab::domain_impl::z_abs_domain_t)
Q_RUNNER(crab::domain_impl::q_abs_domain_t)
#else
Z_RUNNER(crab::domain_impl::z_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_constant_domain_t)
Z_RUNNER(crab::domain_impl::z_ric_domain_t)
Z_RUNNER(crab::domain_impl::z_dbm_domain_t)
Z_RUNNER(crab::domain_impl::z_sdbm_domain_t)
Z_RUNNER(crab::domain_impl::z_soct_domain_t)
Z_RUNNER(crab::domain_impl::z_soct_domain_lw_t)
Z_RUNNER(crab::domain_impl::z_fixed_tvpi_domain_t)
Z_RUNNER(crab::domain_impl::z_term_domain_t)
Z_RUNNER(crab::domain_impl::z_term_dis_int_t)
Z_RUNNER(crab::domain_impl::z_num_domain_t)
Z_RUNNER(crab::domain_impl::z_bool_num_domain_t)
Z_RUNNER(crab::domain_impl::z_bool_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_boxes_domain_t)
Z_RUNNER(crab::domain_impl::z_dis_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_box_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_pk_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_poly_pplite_domain_t)
Z_RUNNER(crab::domain_impl::z_fpoly_pplite_domain_t)
Z_RUNNER(crab::domain_impl::z_pset_pplite_domain_t)
Z_RUNNER(crab::domain_impl::z_dummy_decoupled_box_domain_t)
Z_RUNNER(crab::domain_impl::z_decoupled_box_poly_domain_t)
Z_RUNNER(crab::domain_impl::z_decoupled_box_pset_domain_t)
Z_RUNNER(crab::domain_impl::z_zones_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_oct_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_pk_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_wrapped_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_as_dis_int_t)
Z_RUNNER(crab::domain_impl::z_as_sdbm_t)
Z_RUNNER(crab::domain_impl::z_as_bool_num_t)
Z_RUNNER(crab::domain_impl::z_aa_term_int_t)
Z_RUNNER(crab::domain_impl::z_aa_bool_int_t)
Z_RUNNER(crab::domain_impl::z_pow_aa_int_t)
Z_RUNNER(crab::domain_impl::z_rgn_aa_int_t)
Z_RUNNER(crab::domain_impl::z_rgn_bool_int_t)
Z_RUNNER(crab::domain_impl::z_rgn_constant_t)
Z_RUNNER(crab::domain_impl::z_rgn_int_t)
Z_RUNNER(crab::domain_impl::z_rgn_sign_t)
Z_RUNNER(crab::domain_impl::z_rgn_sign_constant_t)
Z_RUNNER(crab::domain_impl::z_rgn_sdbm_t)

Q_RUNNER(crab::domain_impl::q_interval_domain_t)
Q_RUNNER(crab::domain_impl::q_oct_apron_domain_t)
Q_RUNNER(crab::domain_impl::q_pk_apron_domain_t)
Q_RUNNER(crab::domain_impl::q_poly_pplite_domain_t)
Q_RUNNER(crab::domain_impl::q_fpoly_pplite_domain_t)
Q_RUNNER(crab::domain_impl::q_pset_pplite_domain_t)
Q_RUNNER(crab::domain_impl::q_oct_elina_domain_t)
Q_RUNNER(crab::domain_impl::q_pk_elina_domain_t)
Q_RUNNER(crab::domain_impl::q_boxes_domain_t)
#endif
