#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("bb1", "bb4");
  // adding blocks
  BB(cfg, bb1);
  BB(cfg, bb2);
  BB(cfg, bb3);
  BB(cfg, bb4);
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
  auto cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  z_oct_elina_domain_t initial_states;
  // no thresholds, no narrowing
  backward_run<z_oct_elina_domain_t>(cfg, cfg->entry(), initial_states, 1, 0,
				     0, stats_enabled);

  // free the CFG
#endif

  return 0;
}
