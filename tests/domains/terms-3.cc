#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

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
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, ret);
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb1;
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
  bb1_t.havoc(nd1);
  bb1_t.havoc(nd2);
  bb1_t.bitwise_and(x, x, nd1);
  bb1_t.bitwise_and(y, y, nd1);
  bb1_t.bitwise_or(z, z, nd1);
  bb1_t.bitwise_or(w, w, nd1);
  bb1_t.bitwise_xor(s, nd1, nd2);
  bb1_t.bitwise_xor(t, nd1, nd2);
  bb1_t.add(i, i, 1);

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    run_all<z_interval_domain_t, z_dbm_domain_t, z_term_domain_t>(cfg, stats_enabled);

    return 0;
  });
}
