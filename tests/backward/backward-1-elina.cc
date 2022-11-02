#include "../common.hpp"
#include "../program_options.hpp"

// Example from Monniaux's slides.

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
  bb1.assume(x >= 0);
  bb1.assign(y, x);

  bb2.assume(x <= -1);
  bb2.assign(tmp, 0);
  bb2.sub(y, tmp, x);

  bb3.assume(y >= 1);
  bb3.assertion(x != 0);

  return cfg;
}

int main(int argc, char **argv) {
#ifdef HAVE_ELINA  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

// A forward+backward analysis should prove the assertion holds.
  z_oct_elina_domain_t initial_states;
  backward_run<z_oct_elina_domain_t>(cfg, cfg->entry(), initial_states, 1, 2,
                                     20, stats_enabled);

  // free the CFG
  delete cfg;
#endif

  return 0;
}
