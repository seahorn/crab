#include "../common.hpp"

namespace crab
{
   namespace cfg_impl
   {
     typedef Cfg< basic_block_label_t, varname_t, q_number> q_cfg_t;
     typedef q_cfg_t::basic_block_t                         q_basic_block_t;
     
     typedef cfg_ref<q_cfg_t>                               q_cfg_ref_t;
     typedef cfg_rev<q_cfg_ref_t>                           q_cfg_rev_t;

   }
   namespace domain_impl
   {
    
     using namespace crab::cfg_impl;
     using namespace crab::domains; 
     using namespace ikos;
     
     typedef variable< q_number, varname_t> q_var;
     typedef linear_expression<q_number, varname_t> q_lin_t;
     typedef linear_constraint<q_number, varname_t> q_lin_cst_t;
     typedef linear_constraint_system<q_number, varname_t> q_lin_cst_sys_t;
     typedef interval<q_number> q_interval_t;
     typedef bound<q_number> q_bound_t;

     // Numerical domains
     typedef interval_domain<q_number, varname_t> q_interval_domain_t;     
     typedef apron_domain<q_number, varname_t, apron_domain_id_t::APRON_PK> q_pk_apron_domain_t;
     typedef boxes_domain< q_number, varname_t > q_boxes_domain_t;     
   }
}

// To run abstract domains defined over reals
template<typename Dom>
void run (crab::cfg_impl::q_cfg_t* cfg, 
	  crab::cfg_impl::variable_factory_t& vfac, 
	  bool run_liveness,
	  unsigned widening, 
	  unsigned narrowing, 
	  unsigned jump_set_size)
{
  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::q_cfg_ref_t,Dom> intra_fwd_analyzer_t;
  intra_run_impl<crab::cfg_impl::q_cfg_t, Dom, intra_fwd_analyzer_t>
    (cfg, vfac, run_liveness, widening, narrowing, jump_set_size);
}

using namespace std;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
q_cfg_t* prog (variable_factory_t &vfac)  {

  // Definining program variables
  q_var x (vfac ["x"]);
  q_var y (vfac ["y"]);  
  // entry and exit block
  q_cfg_t* cfg = new q_cfg_t("entry","bb1");
  // adding blocks
  q_basic_block_t& entry = cfg->insert ("entry");
  q_basic_block_t& bb1   = cfg->insert ("bb1");
  // adding control flow 
  entry >> bb1;
  bb1 >> bb1;
  // adding statements
  entry.assign (x, q_number(0.0));
  entry.assign (y, q_number(0.5));
  bb1.add(x, x, y);
  bb1.div(y, y, q_number(2));  

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  q_cfg_t* cfg = prog (vfac);
  crab::outs() << *cfg << "\n";

  run<q_interval_domain_t>(cfg,vfac,false,1,2,20);    
  run<q_pk_apron_domain_t>(cfg,vfac,false,1,2,20);  
  #ifdef HAVE_LDD
  run<q_boxes_domain_t>(cfg,vfac,false,1,2,20);
  #endif   
  return 0;
}
