#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {

  /*
     i := 0;
     while (i<=3) {
       j:=0;
       while (j <= 3) {
         assert (i <= j + 3);
         i++;
         j++;
       }
       i := i - j + 1;
     }
   */
  z_cfg_t *cfg = new z_cfg_t("entry");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &l1 = cfg->insert("l1");
  z_basic_block_t &l1_entry = cfg->insert("l1_entry");
  z_basic_block_t &l1_cont = cfg->insert("l1_cont");
  z_basic_block_t &l2 = cfg->insert("l2");
  z_basic_block_t &l2_body = cfg->insert("l2_body");
  z_basic_block_t &l2_exit = cfg->insert("l2_exit");

  entry >> l1;
  // outer loop
  l1 >> l1_entry;
  l1_entry >> l2;
  // inner loop
  l2 >> l2_body;
  l2_body >> l2;
  l2 >> l2_exit;
  // outer loop again
  l2_exit >> l1_cont;
  l1_cont >> l1;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  l1.assume(i <= 3);
  l1_entry.assign(j, 0);

  l2_body.assume(j <= 3);
  l2_body.assertion(i <= j + 3);
  l2_body.add(i, i, 1);
  l2_body.add(j, j, 1);

  l2_exit.assume(j >= 4);
  l1_cont.assign(i, i - j + 1);

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
    run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_sdbm_domain_t init;
    run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  delete cfg;
  return 0;
}
