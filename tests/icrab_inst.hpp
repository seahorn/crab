#pragma once

// Here all explicit instantiations

#define Z_RUNNER(BUDOM,TDDOM)		    \
template void inter_run<BUDOM,TDDOM>	    \
(crab::cg_impl::z_cg_t*, bool, unsigned, unsigned, unsigned, bool);

Z_RUNNER(crab::domain_impl::z_dbm_domain_t, crab::domain_impl::z_interval_domain_t)
#ifdef HAVE_APRON
Z_RUNNER(crab::domain_impl::z_oct_apron_domain_t, crab::domain_impl::z_interval_domain_t)
#endif
Z_RUNNER(crab::domain_impl::z_term_domain_t, crab::domain_impl::z_interval_domain_t)
Z_RUNNER(crab::domain_impl::z_num_domain_t, crab::domain_impl::z_num_domain_t)

