#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog1(variable_factory_t &vfac)  {

  /*
    x = 0;
    y = 1;
    if (nd()) {
      x += 1;
      y += 2;
    }
    if (nd()) {
      x += 2;
      y += 3;
    }
    if (nd()) {
      y += 1;
    }
    assert(y > x);
   */
  
  // Defining program variables
  z_var x (vfac ["x"], crab::INT_TYPE, 32);
  z_var y (vfac ["y"], crab::INT_TYPE, 32);  
  z_var nd1 (vfac ["nd1"], crab::INT_TYPE, 32);
  z_var nd2 (vfac ["nd2"], crab::INT_TYPE, 32);
  z_var nd3 (vfac ["nd3"], crab::INT_TYPE, 32);
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","exit");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& b1 = cfg->insert ("b1");
  z_basic_block_t& b1_tt = cfg->insert ("b1_tt");
  z_basic_block_t& b1_ff = cfg->insert ("b1_ff");
  z_basic_block_t& b2 = cfg->insert ("b2");
  z_basic_block_t& b2_tt = cfg->insert ("b2_tt");
  z_basic_block_t& b2_ff = cfg->insert ("b2_ff");
  z_basic_block_t& b3 = cfg->insert ("b3");
  z_basic_block_t& b3_tt = cfg->insert ("b3_tt");
  z_basic_block_t& b3_ff = cfg->insert ("b3_ff");
  z_basic_block_t& exit   = cfg->insert ("exit");
  
  // adding control flow
  entry >> b1;
  b1 >> b1_tt;
  b1 >> b1_ff;
  b1_tt >> b2;
  b1_ff >> b2;  

  b2 >> b2_tt;
  b2 >> b2_ff;
  b2_tt >> b3;
  b2_ff >> b3;  

  b3 >> b3_tt;
  b3 >> b3_ff;
  b3_tt >> exit;
  b3_ff >> exit;  
  
  // adding statements
  entry.assign(x,0);
  entry.assign(y,1);

  b1_tt.assume(nd1 >= 1);
  b1_tt.add(x, x, 1);
  b1_tt.add(y, y, 2);
  b1_ff.assume(nd1 <= 0);

  b2_tt.assume(nd2 >= 1);
  b2_tt.add(x, x, 2);
  b2_tt.add(y, y, 3);
  b2_ff.assume(nd2 <= 0);

  b3_tt.assume(nd3 >= 1);
  b3_tt.add(y, y, 1);
  b3_ff.assume(nd3 <= 0);

  exit.assertion(y >= x + 1);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  variable_factory_t vfac;
  z_cfg_t* cfg = prog1(vfac);
  crab::outs() << *cfg << "\n";

  run_and_check<z_pow_aa_int_t>(cfg,cfg->entry(),false,1,2,20,stats_enabled);  

  // free the CFG
  delete cfg;

  return 0;
}
