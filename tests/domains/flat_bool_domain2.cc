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

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

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
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
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
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    // precise
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";
    z_bool_interval_domain_t init;
    run(cfg, init, stats_enabled);

    return 0;
  });
}
