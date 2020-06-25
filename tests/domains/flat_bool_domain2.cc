#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/*
 * Crab distinguishes between integer and booleans.
 * Example of program with booleans
 */

z_cfg_t *prog(variable_factory_t &vfac) {

  /*
     i := 0;
     n := *;
     b := (n == 10);
     assume(b);
     while (i <= n) {
       i++;
     }
     assert(i == 10);
   */

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::BOOL_TYPE, 1);
  z_var n(vfac["n"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(i, z_number(0));
  entry.havoc(n);
  entry.bool_assign(b, n == z_number(10));
  // entry.assign(n, z_number(1));
  entry.bool_assume(b);
  bb1_t.assume(i <= n);
  bb2.add(i, i, 1);
  bb1_f.assume(i >= n + 1);
  ret.assertion(i == z_number(10));

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  // precise
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";
  z_bool_interval_domain_t init;
  run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  delete cfg;

  return 0;
}
