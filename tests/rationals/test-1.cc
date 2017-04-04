#include "../common.hpp"

using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
q_cfg_t* prog (variable_factory_t &vfac)  {

  // Definining program variables
  q_var i (vfac ["i"]);
  // entry and exit block
  q_cfg_t* cfg = new q_cfg_t("entry","ret");
  // adding blocks
  q_basic_block_t& entry = cfg->insert ("entry");
  q_basic_block_t& bb1   = cfg->insert ("bb1");
  q_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  q_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  q_basic_block_t& bb2   = cfg->insert ("bb2");
  q_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow 
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, q_number(0.0));
  bb1_t.assume (i <= q_number(9.9));
  bb1_f.assume (i >= q_number(10));
  bb2.add(i, i, q_number(1));

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  q_cfg_t* cfg = prog (vfac);
  cfg->simplify (); // this is optional
  crab::outs() << *cfg << "\n";

  run<q_interval_domain_t>(cfg,vfac,false,1,2,20);
  run<q_pk_apron_domain_t>(cfg,vfac,false,1,2,20);  

  return 0;
}
