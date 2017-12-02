#include "../../program_options.hpp"
#include "../../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog (variable_factory_t &vfac)  {

  // Defining program variables
  z_var x (vfac["x"], crab::INT_TYPE, 8);
  z_var y (vfac["y"], crab::INT_TYPE, 8);
  z_var nd (vfac["nd"], crab::INT_TYPE, 32);
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb_nd = cfg->insert ("bb_nd");
  z_basic_block_t& bb_nd_tt = cfg->insert ("bb_nd_tt");
  z_basic_block_t& bb_nd_ff = cfg->insert ("bb_nd_ff");      
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb_nd;
  bb_nd >> bb_nd_tt;
  bb_nd >> bb_nd_ff;
  bb_nd_tt >> bb1;
  bb_nd_ff >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign(y, z_number(-10));
  entry.havoc(nd);
  bb_nd_tt.assume(nd >= 1);
  bb_nd_tt.assign(x, z_number(0));
  bb_nd_ff.assume(nd <= 0);
  bb_nd_ff.assign(x, z_number(100));		  
  bb1_t.assume (x >= y);
  bb1_f.assume (x <= y-1);
  bb2.sub(x, x, y);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  // unsound result
  run<z_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
  // sound result
  run<z_wrapped_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);

  // free the CFG
  delete cfg;

  return 0;
}
