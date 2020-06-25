#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {

  /*
    i := 0;
    k := 30;
    while (i <= 9) {
      i := i + 1;
    }
    j := 0;
    while (j <= 9) {
      j := i + 1;
    }
 */

  z_cfg_t *cfg = new z_cfg_t("loop1_entry", "ret");
  // z_cfg_t cfg ("loop1_entry");
  z_basic_block_t &loop1_entry = cfg->insert("loop1_entry");
  z_basic_block_t &loop1_bb1 = cfg->insert("loop1_bb1");
  z_basic_block_t &loop1_bb1_t = cfg->insert("loop1_bb1_t");
  z_basic_block_t &loop1_bb1_f = cfg->insert("loop1_bb1_f");
  z_basic_block_t &loop1_bb2 = cfg->insert("loop1_bb2");
  z_basic_block_t &loop2_entry = cfg->insert("loop2_entry");
  z_basic_block_t &loop2_bb1 = cfg->insert("loop2_bb1");
  z_basic_block_t &loop2_bb1_t = cfg->insert("loop2_bb1_t");
  z_basic_block_t &loop2_bb1_f = cfg->insert("loop2_bb1_f");
  z_basic_block_t &loop2_bb2 = cfg->insert("loop2_bb2");
  z_basic_block_t &ret = cfg->insert("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t;
  loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2;
  loop1_bb2 >> loop1_bb1;
  loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t;
  loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2;
  loop2_bb2 >> loop2_bb1;
  loop2_bb1_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  loop1_entry.assign(i, 0);
  loop1_entry.assign(k, 30);
  loop1_bb1_t.assume(i <= 9);
  loop1_bb1_f.assume(i >= 10);
  loop1_bb2.add(i, i, 1);

  loop2_entry.assign(j, 0);
  loop2_bb1_t.assume(j <= 9);
  loop2_bb1_f.assume(j >= 10);
  loop2_bb2.add(j, j, 1);
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

  delete cfg;
  return 0;
}
