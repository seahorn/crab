#include "../common.hpp"
#include "../program_options.hpp"

// Modified version from Rival'05 paper

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "bb3");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  // adding control flow
  entry >> bb1;
  entry >> bb2;
  bb1 >> bb3;
  bb2 >> bb3;

  // adding statements
  bb1.assume(x >= 1);
  bb1.assign(y, x);

  bb2.assume(x <= 0);
  bb2.assign(tmp, 0);
  bb2.sub(y, tmp, x);

  bb3.assume(y >= 6);
  bb3.assertion(x <= 5);
  bb3.assertion(x >= -5);
  return cfg;
}

int main(int argc, char **argv) {
#ifdef HAVE_LDD  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  // A forward+backward analysis cannot prove the two assertion holds,
  // even used disjunctive boxes.
  // - for classical intervals, we lose all information going
  //   backwards after the assert(x<=5)
  // - for disjunctive intervals, we can keep that error can happen
  //   either x <= -6 or x >= 6.  However, when going backwards
  //   through blocks bb1 and bb2 one of the two constraints is true
  //   so we never infer bottom.

  z_boxes_domain_t initial_states;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  backward_run<z_boxes_domain_t>(cfg, cfg->entry(), initial_states, 1, 2, 20,
				 stats_enabled);

  // free the CFG
  delete cfg;
#endif

  return 0;
}
