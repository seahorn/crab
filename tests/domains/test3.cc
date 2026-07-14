#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {
  /*
      i := 0;
      while (i <= 10) {
        i:= i+1;
        if (i >= 9) {
          i:=0;
        }
      }
      while (i <= 100) {
        i:=i-1;
      }
   */
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  BB(cfg, entry);
  BB(cfg, loop1_head);
  BB(cfg, loop1_t);
  BB(cfg, loop1_f);
  BB(cfg, loop1_body);

  BB(cfg, loop1_body_t);
  BB(cfg, loop1_body_f);
  BB(cfg, loop1_body_x);

  BB(cfg, cont);
  BB(cfg, loop2_head);
  BB(cfg, loop2_t);
  BB(cfg, loop2_f);
  BB(cfg, loop2_body);
  BB(cfg, ret);

  entry >> loop1_head;
  loop1_head >> loop1_t;
  loop1_head >> loop1_f;
  loop1_t >> loop1_body;

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t;
  loop2_head >> loop2_f;
  loop2_t >> loop2_body;
  loop2_body >> loop2_head;
  loop2_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  loop1_t.assume(i <= 10);
  loop1_f.assume(i >= 11);
  loop1_body.add(i, i, 1);

  loop1_body_t.assume(i >= 9);
  loop1_body_t.assign(i, 0);
  loop1_body_f.assume(i <= 8);

  loop2_t.assume(i <= 100);
  loop2_f.assume(i >= 101);
  loop2_body.sub(i, i, 1);
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    // cfg->simplify ();
    crab::outs() << *cfg << "\n";

    {
      z_interval_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_dbm_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_sdbm_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_ric_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_term_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_dis_interval_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }

    return 0;
  });
}
