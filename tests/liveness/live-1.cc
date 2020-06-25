#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/dataflow/liveness.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {
  /*
     i := 0;
     x := 1;
     y := 0;
     z := 3;
     w := 3;
     while (i < 100) {
       x  := x + y;
       y  := y + 1;
       nd := *;
       z  := z xor nd;
       w  := w xor nd;
       i  := i + 1;
     }
   */
  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
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
  z_basic_block_t &exit = cfg->insert("exit");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> exit;
  exit >> ret;
  // adding statements
  entry.assign(i, 0);
  entry.assign(x, 1);
  entry.assign(y, 0);
  entry.assign(z, 3);
  entry.assign(w, 3);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.havoc(nd1);
  bb2.havoc(nd2);
  bb2.add(x, x, y);
  bb2.add(y, y, 1);
  bb2.bitwise_xor(z, z, nd1);
  bb2.bitwise_xor(w, w, nd1);
  bb2.add(i, i, 1);
  exit.assume(x <= y);

  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  using liveness_t = crab::analyzer::live_and_dead_analysis<z_cfg_ref_t>;
  liveness_t live(*cfg);
  live.exec();

  // use --log=Liveness --log=liveness

  delete cfg;
  return 0;
}
