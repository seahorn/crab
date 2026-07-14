#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

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
  // entry and exit block (#3: unique_ptr, no manual delete)
  auto cfg = std::make_unique<z_cfg_t>("x0", "ret");
  // adding blocks (#3: BB writes the label once)
  BB(cfg, x0);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, ret);
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
  // (#4) test_main handles the parse/stats/early-exit preamble.
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    // (#2) sweep the CFG through several domains; (#1) fixpoint knobs default
    // to widening=1, narrowing=2, jump_set_size=20, no liveness.
    run_all<z_interval_domain_t, z_dbm_domain_t, z_sdbm_domain_t, z_ric_domain_t,
            z_term_domain_t, z_dis_interval_domain_t>(cfg, stats_enabled);
    return 0;
  });
}
