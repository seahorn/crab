#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  /*
    x := 0; 
    y := 0;
    while (true) {
      if (x<=50) {
        y++;
      } else {
        y--;
      }
      if (y < 0) break;
      x++;
    }
   */
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  z_basic_block_t &bb5 = cfg->insert("bb5");
  z_basic_block_t &bb6 = cfg->insert("bb6");
  z_basic_block_t &bb7 = cfg->insert("bb7");    
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> bb1;
  bb1 >> bb2;
  bb1 >> bb3;
  bb2 >> bb4;
  bb3 >> bb4;
  bb4 >> bb5;
  bb4 >> bb7;
  bb5 >> bb6;
  bb6 >> bb1;
  bb7 >> exit;
  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 0);
  bb2.assume(x <= 50);
  bb2.add(y, y, 1);
  bb3.assume(x >= 51);
  bb3.sub(y, y, 1);
  bb5.assume(y >= 0);
  bb6.add(x, x, 1);
  //bb6.assertion(x+y <= 102);
  bb7.assume(y <= -1);

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  {
    z_soct_domain_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 0, stats_enabled);
  }

  {
    z_soct_domain_lw_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 0, stats_enabled);
  }
  
  // free the CFG
  delete cfg;

  return 0;
}
