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
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("x0", "ret");
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &x1 = cfg->insert("x1");
  z_basic_block_t &x2 = cfg->insert("x2");
  z_basic_block_t &x3 = cfg->insert("x3");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> x1;
  x1 >> x2;
  x2 >> x3;
  x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  x0.assign(k, 50);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);

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
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_term_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  return 0;
}
