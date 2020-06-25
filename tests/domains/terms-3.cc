#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var s(vfac["s"], crab::INT_TYPE, 32);
  z_var t(vfac["t"], crab::INT_TYPE, 32);
  z_var nd1(vfac["nd1"], crab::INT_TYPE, 32);
  z_var nd2(vfac["nd2"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
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
  entry.assign(i, 0);
  entry.assign(x, 5);
  entry.assign(y, 5);
  entry.assign(z, 3);
  entry.assign(w, 3);
  entry.assign(s, 0);
  entry.assign(t, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.havoc(nd1);
  bb2.havoc(nd2);
  bb2.bitwise_and(x, x, nd1);
  bb2.bitwise_and(y, y, nd1);
  bb2.bitwise_or(z, z, nd1);
  bb2.bitwise_or(w, w, nd1);
  bb2.bitwise_xor(s, nd1, nd2);
  bb2.bitwise_xor(t, nd1, nd2);
  bb2.add(i, i, 1);

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
  cfg->simplify(); // this is optional
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
    z_term_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  return 0;
}
