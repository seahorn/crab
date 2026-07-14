#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog1(variable_factory_t &vfac) {

  /* APLAS'12 example
     x and y are int8
     
     y = -10;
     assume(x >= 0 && x <= 100);
     while (x >= y) {
        x = x - y;
     }
       
     The expected result at the end is x=[-128,-119]
  */
  
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 8);
  z_var y(vfac["y"], crab::INT_TYPE, 8);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb_nd);
  BB(cfg, bb_nd_tt);
  BB(cfg, bb_nd_ff);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
  // adding control flow
  entry >> bb_nd;
  bb_nd >> bb_nd_tt;
  bb_nd >> bb_nd_ff;
  bb_nd_tt >> bb1;
  bb_nd_ff >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(y, z_number(-10));
  entry.havoc(nd);
  bb_nd_tt.assume(nd >= 1);
  bb_nd_tt.assign(x, z_number(0));
  bb_nd_ff.assume(nd <= 0);
  bb_nd_ff.assign(x, z_number(100));
  bb1_t.assume(x >= y);
  bb1_f.assume(x <= y - 1);
  bb2.sub(x, x, y);
  return cfg;
}

std::unique_ptr<z_cfg_t> prog2(variable_factory_t &vfac) {
  /*
     x= 127;
     x= x+1;
     if (x <= 1) {
       x = 10;
     } else {
       x = -10;
     }

   */
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 8);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  z_basic_block_t &bb_if = cfg->insert("if");
  z_basic_block_t &bb_then = cfg->insert("then");
  BB(cfg, ret);
  // adding control flow
  entry >> bb_if;
  entry >> bb_then;
  bb_if >> ret;
  bb_then >> ret;
  // adding statements
  entry.assign(x, z_number(127));
  entry.add(x, x, 1);
  z_lin_cst_t c1(x <= z_number(1));
  bb_if.assume(c1);
  bb_if.assign(x, z_number(10));
  z_lin_cst_t c2(x >= z_number(2));
  bb_then.assume(c2);
  bb_then.assign(x, z_number(-10));
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    {
      variable_factory_t vfac;
      auto cfg = prog1(vfac);
      crab::outs() << *cfg << "\n";
      {
        // unsound result
        z_interval_domain_t init;
        run(cfg, init, stats_enabled);
      }
      {
        // sound result
        z_wrapped_interval_domain_t init;
        run(cfg, init, stats_enabled);
      }
    }
    {
      variable_factory_t vfac;
      auto cfg = prog2(vfac); 
      crab::outs() << *cfg << "\n";
      z_wrapped_interval_domain_t init;
      run(cfg, init, stats_enabled);
    }
    return 0;
  });
}
