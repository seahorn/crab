#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

  /*
     i := 0;

     while (true) {
       i := i + 1;
       j := 0;

       while (j <= 9) {
          assert(i >= 0);
          assert(i <= 10);
          j := j + 1;
       }
       if (i >= 10) {
          i := 0;
       }
     }
   */
  auto cfg = std::make_unique<z_cfg_t>("entry");
  BB(cfg, entry);
  BB(cfg, l1);
  BB(cfg, l1_entry);
  BB(cfg, l1_reset);
  BB(cfg, l1_dont_reset);
  BB(cfg, l2);
  BB(cfg, l2_body);
  BB(cfg, l2_exit);

  entry >> l1;
  // outer loop
  l1 >> l1_entry;
  l1_entry >> l2;
  // inner loop
  l2 >> l2_body;
  l2_body >> l2;
  l2 >> l2_exit;
  // outer loop again
  l2_exit >> l1_reset;
  l2_exit >> l1_dont_reset;
  l1_reset >> l1;
  l1_dont_reset >> l1;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  l1_entry.add(i, i, 1);
  l1_entry.assign(j, 0);

  l2_body.assume(j <= 9);
  l2_body.assertion(i >= 0);
  l2_body.assertion(i <= 10);
  l2_body.add(j, j, 1);

  l2_exit.assume(j >= 10);

  l1_reset.assume(i >= 10);
  l1_reset.assign(i, 0);
  l1_dont_reset.assume(i <= 9);

  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    cfg->simplify();
    crab::outs() << *cfg << "\n";
    z_interval_domain_t init;
    run_and_check(cfg, init, stats_enabled);

    return 0;
  });
}
