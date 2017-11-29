#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog (variable_factory_t &vfac)  {

  // Defining program variables
  z_var i(vfac ["i"]);
  z_var b1(vfac ["b1"]);
  z_var n(vfac ["n"]);
  z_var n1(vfac ["n1"]);
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, z_number(0));
  entry.havoc(n);
  entry.truncate(n, 64, n1, 32);
  entry.bool_assign(b1, n1 == z_number(9));
  entry.bool_assume(b1);
  bb1_t.assume (i <= n);
  bb2.add(i,i,1);
  bb1_f.assume (i >= n + 1);
  ret.assertion(i == 10);

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  run<z_bool_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
  run<z_bool_num_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);  
  
  // free the CFG
  delete cfg;

  return 0;
}
