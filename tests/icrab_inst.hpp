#pragma once

#include "./icrab.hpp"

// Here all explicit instantiations

#define Z_BU_RUNNER(BUDOM,TDDOM)	    \
template void bu_inter_run<BUDOM,TDDOM>	    \
(crab::cg_impl::z_cg_t*, BUDOM, TDDOM, bool, unsigned, unsigned, unsigned, bool);

#define Z_TD_RUNNER(DOM)	    \
template void td_inter_run<DOM>	    \
(crab::cg_impl::z_cg_t*, DOM, td_inter_params_t, bool, bool, bool);


Z_BU_RUNNER(crab::domain_impl::z_dbm_domain_t, crab::domain_impl::z_interval_domain_t)
#ifdef HAVE_APRON
Z_BU_RUNNER(crab::domain_impl::z_oct_apron_domain_t, crab::domain_impl::z_interval_domain_t)
#endif
Z_BU_RUNNER(crab::domain_impl::z_term_domain_t, crab::domain_impl::z_interval_domain_t)
Z_BU_RUNNER(crab::domain_impl::z_num_domain_t, crab::domain_impl::z_num_domain_t)


Z_TD_RUNNER(crab::domain_impl::z_dbm_domain_t)
Z_TD_RUNNER(crab::domain_impl::z_sdbm_domain_t)
#ifdef HAVE_APRON
Z_TD_RUNNER(crab::domain_impl::z_oct_apron_domain_t)
#endif
#ifdef HAVE_ELINA
Z_TD_RUNNER(crab::domain_impl::z_oct_elina_domain_t)
#endif



