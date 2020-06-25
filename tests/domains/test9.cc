#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *cfg1(variable_factory_t &vfac) {

  /*
   i := 0;
   x := 1;
   y := 0;

    while(i <= 99) {
     x = x + y;
     y = y + 1;
     i = i + 1;
    }
    z := x;
    w := y;
    assert(z>= w);
   */

  // === Define program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  // === Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // === Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &ret = cfg->insert("ret");
  // === Adding CFG edges
  entry.add_succ(bb1);
  bb1.add_succ(bb1_t);
  bb1.add_succ(bb1_f);
  bb1_t.add_succ(bb2);
  bb2.add_succ(bb1);
  bb1_f.add_succ(bb3);
  bb3.add_succ(ret);
  // === Adding statements
  entry.assign(i, 0);
  entry.assign(x, 1);
  entry.assign(y, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(x, x, y);
  bb2.add(y, y, 1);
  bb2.add(i, i, 1);
  ret.assign(z, x);
  ret.assign(w, y);
  ret.assertion(z >= w);
  return cfg;
}

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_cfg_t *p1 = cfg1(vfac);
  crab::outs() << *p1 << "\n";
  z_sdbm_domain_t init;
  run_and_check(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
  delete p1;

  return 0;
}
