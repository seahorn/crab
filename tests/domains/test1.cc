#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  /*
    k := 2147483648;
    i := 0;
    while (i <= 99) {
      i += (* ? 1: 2);
    }
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("x0", "ret");
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb1;
  bb1_f >> ret;
  // adding statements
  x0.assign(k, 2147483648);
  x0.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb1_t.havoc(nd);
  bb1_t.select(inc, nd, 1, 2);
  bb1_t.add(i, i, inc);

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
    z_interval_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_dbm_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_sdbm_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_ric_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_term_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_dis_interval_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  // free the CFG
  delete cfg;

  return 0;
}
