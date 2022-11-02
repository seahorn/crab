#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("bb1", "bb4");
  // adding blocks
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  // adding control flow
  bb1 >> bb2;
  bb2 >> bb3;
  bb3 >> bb2;
  bb2 >> bb4;

  // adding statements
  bb1.assign(x, 0);
  bb3.assume(x <= 99);
  bb3.add(x, x, 1);
  bb4.assume(x >= 100);
  bb4.assertion(x <= 100);
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

  z_oct_elina_domain_t initial_states;
  // no thresholds, no narrowing
  backward_run<z_oct_elina_domain_t>(cfg, cfg->entry(), initial_states, 1, 0,
				     0, stats_enabled);

  // free the CFG
  delete cfg;
#endif

  return 0;
}
