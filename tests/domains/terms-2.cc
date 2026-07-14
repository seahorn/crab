#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {
  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var z0(vfac["z0"], crab::INT_TYPE, 32);
  z_var y0(vfac["y0"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("p0", "ret");
  // adding blocks
  BB(cfg, p0);
  BB(cfg, p_neg);
  BB(cfg, p_pos);
  BB(cfg, exit);
  BB(cfg, ret);
  // adding control flow
  p0 >> p_pos;
  p0 >> p_neg;
  p_neg >> exit;
  p_pos >> exit;
  exit >> ret;

  // adding statements
  p0.assign(x, 50);
  p0.havoc(y);
  p0.assume(y >= -1);
  p0.assume(y <= 1);
  p0.mul(z, x, y);

  p_neg.assume(y <= -1);
  p_neg.mul(z, z, -1);

  p_pos.assume(y >= 0);

  exit.assign(z0, z);
  exit.assign(y0, y);

  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    {
      z_interval_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_term_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }

    return 0;
  });
}
