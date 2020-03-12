#pragma once

// Here all explicit instantiations

#define Z_RUNNER(DOM) \
  template void run<DOM > \
  (crab::cfg_impl::z_cfg_t*, crab::cfg_impl::basic_block_label_t, bool, unsigned, unsigned, unsigned, bool);

#define Z_CHECKER(DOM) \
  template void run_and_check<DOM > \
  (crab::cfg_impl::z_cfg_t*, crab::cfg_impl::basic_block_label_t, bool, unsigned, unsigned, unsigned, bool);

#define Q_RUNNER(DOM) \
template void run<DOM > \
(crab::cfg_impl::q_cfg_t*, crab::cfg_impl::basic_block_label_t, bool, unsigned, unsigned, unsigned, bool);

#define Q_CHECKER(DOM) \
template void run_and_check<DOM > \
(crab::cfg_impl::q_cfg_t*, crab::cfg_impl::basic_block_label_t, bool, unsigned, unsigned, unsigned, bool);

Z_RUNNER(crab::domain_impl::z_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_ric_domain_t)
Z_RUNNER(crab::domain_impl::z_dbm_domain_t)
Z_RUNNER(crab::domain_impl::z_sdbm_domain_t)
Z_RUNNER(crab::domain_impl::z_term_domain_t)
Z_RUNNER(crab::domain_impl::z_term_dis_int_t)
Z_RUNNER(crab::domain_impl::z_num_domain_t)
Z_RUNNER(crab::domain_impl::z_bool_num_domain_t)
Z_RUNNER(crab::domain_impl::z_bool_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_boxes_domain_t)
Z_RUNNER(crab::domain_impl::z_dis_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_nullity_domain_t)
Z_RUNNER(crab::domain_impl::z_box_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_pk_apron_domain_t)
Z_RUNNER(crab::domain_impl::z_zones_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_oct_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_pk_elina_domain_t)
Z_RUNNER(crab::domain_impl::z_wrapped_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_as_dis_int_t)
Z_RUNNER(crab::domain_impl::z_as_sdbm_t)
Z_RUNNER(crab::domain_impl::z_as_num_null_t) 
Z_RUNNER(crab::domain_impl::z_as_bool_num_t)
Z_RUNNER(crab::domain_impl::z_ag_sdbm_intv_t)
Z_RUNNER(crab::domain_impl::z_ag_num_null_t)
Z_RUNNER(crab::domain_impl::z_ae_term_int_t)
Z_RUNNER(crab::domain_impl::z_aa_term_int_t)
Z_RUNNER(crab::domain_impl::z_aa_bool_int_t)

Z_CHECKER(crab::domain_impl::z_interval_domain_t)
Z_CHECKER(crab::domain_impl::z_sdbm_domain_t)
Z_CHECKER(crab::domain_impl::z_ae_term_int_t)
Z_CHECKER(crab::domain_impl::z_aa_term_int_t)

Q_RUNNER(crab::domain_impl::q_interval_domain_t)
Q_RUNNER(crab::domain_impl::q_oct_apron_domain_t)
Q_RUNNER(crab::domain_impl::q_pk_apron_domain_t)
Q_RUNNER(crab::domain_impl::q_oct_elina_domain_t)
Q_RUNNER(crab::domain_impl::q_pk_elina_domain_t)
Q_RUNNER(crab::domain_impl::q_boxes_domain_t)  
