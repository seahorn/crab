#pragma once
// Here all explicit instantiations

#define Z_BWD_RUNNER(DOM) \
template void backward_run<DOM> \
 (crab::cfg_impl::z_cfg_t*, crab::cfg_impl::basic_block_label_t, DOM, unsigned, unsigned, unsigned, bool);

Z_BWD_RUNNER(crab::domain_impl::z_interval_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_sdbm_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_term_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_term_dis_int_t)
//Z_BWD_RUNNER(crab::domain_impl::z_num_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_bool_num_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_boxes_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_dis_interval_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_nullity_domain_t)
//Z_BWD_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_box_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_pk_apron_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_oct_elina_domain_t)
Z_BWD_RUNNER(crab::domain_impl::z_ae_int_t)
Z_BWD_RUNNER(crab::domain_impl::z_ae_box_apron_t)
Z_BWD_RUNNER(crab::domain_impl::z_ae_zones_elina_t)
