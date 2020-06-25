#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {
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
  z_cfg_t *cfg = new z_cfg_t("p0", "ret");
  // adding blocks
  z_basic_block_t &p0 = cfg->insert("p0");
  z_basic_block_t &p_neg = cfg->insert("p_neg");
  z_basic_block_t &p_pos = cfg->insert("p_pos");
  z_basic_block_t &exit = cfg->insert("exit");
  z_basic_block_t &ret = cfg->insert("ret");
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
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  cfg->simplify();
  crab::outs() << *cfg << "\n";

  {
    z_interval_domain_t init;
    run(cfg, cfg->entry(), init, true, 1, 2, 20, stats_enabled);
  }
  {
    z_term_domain_t init;
    run(cfg, cfg->entry(), init, true, 1, 2, 20, stats_enabled);
  }

  return 0;
}
