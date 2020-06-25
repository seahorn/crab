#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {

  /*
    i, n: 64 integer
    n1  : 32 integer

    i  := 0;
    n  := *;
    n1 := truncate(n);
    b1 := (n1 == 9);
    assume(b1);
    while (i <= n) {
       i++;
    }
    assert(i == 10);
   */

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 64);
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var n(vfac["n"], crab::INT_TYPE, 64);
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);

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
  entry.truncate(n, n1);
  entry.bool_assign(b1, n1 == z_number(9));
  entry.bool_assume(b1);
  bb1_t.assume(i <= n);
  bb2.add(i, i, 1);
  bb1_f.assume(i >= n + 1);
  ret.assertion(i == 10);

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
    z_bool_interval_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_bool_num_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  // free the CFG
  delete cfg;

  return 0;
}
