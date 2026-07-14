#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

  /*
    i := 0;
    k := 30;
    while (i <= 9) {
      i := i + 1;
    }
    j := 0;
    while (j <= 9) {
      j := i + 1;
    }
 */

  auto cfg = std::make_unique<z_cfg_t>("loop1_entry", "ret");
  // z_cfg_t cfg ("loop1_entry");
  BB(cfg, loop1_entry);
  BB(cfg, loop1_bb1);
  BB(cfg, loop1_bb1_t);
  BB(cfg, loop1_bb1_f);
  BB(cfg, loop2_bb1);
  BB(cfg, loop2_bb1_t);
  BB(cfg, loop2_bb1_f);
  BB(cfg, ret);

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t;
  loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >>  loop1_bb1;
  loop1_bb1_f >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t;
  loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb1;
  loop2_bb1_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  loop1_entry.assign(i, 0);
  loop1_entry.assign(k, 30);
  loop1_bb1_t.assume(i <= 9);
  loop1_bb1_f.assume(i >= 10);
  loop1_bb1_f.assign(j, 0);
  
  loop1_bb1_t .add(i, i, 1);

  loop2_bb1_t.assume(j <= 9);
  loop2_bb1_t.add(j, j, 1);  
  loop2_bb1_f.assume(j >= 10);
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    run_all<z_interval_domain_t, z_dbm_domain_t, z_sdbm_domain_t, z_ric_domain_t, z_term_domain_t, z_dis_interval_domain_t>(cfg, stats_enabled);

    return 0;
  });
}
