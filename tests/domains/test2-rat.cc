#include "../program_options.hpp"
#include "../common.hpp"


// To run abstract domains defined over reals

using namespace crab::cfg;
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

  run<q_interval_domain_t>(cfg,false,1,2,20,stats_enabled);
  #ifdef HAVE_APRON
  run<q_pk_apron_domain_t>(cfg,false,1,2,20,stats_enabled);
  #endif
  #ifdef HAVE_LDD
  run<q_boxes_domain_t>(cfg,false,1,2,20,stats_enabled);
  #endif   
  return 0;
}
