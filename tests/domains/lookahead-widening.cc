#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

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
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb2);
  BB(cfg, bb3);
  BB(cfg, bb4);
  BB(cfg, bb5);
  BB(cfg, bb6);
  BB(cfg, bb7);    
  BB(cfg, exit);
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
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    {
      z_soct_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_widening(2).with_narrowing(1).with_jump_set_size(0));
    }

    {
      z_soct_domain_lw_t init;
      run(cfg, init, stats_enabled, run_config().with_widening(2).with_narrowing(1).with_jump_set_size(0));
    }
  
    // free the CFG

    return 0;
  });
}
