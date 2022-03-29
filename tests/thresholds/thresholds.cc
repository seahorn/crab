#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog1(variable_factory_t &vfac) {
  // Definining program variables
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &l_1 = cfg->insert("l_1");
  z_basic_block_t &l_2 = cfg->insert("l_2");
  z_basic_block_t &l_3 = cfg->insert("l_3");
  z_basic_block_t &l_4 = cfg->insert("l_4");
  z_basic_block_t &l_5 = cfg->insert("l_5");
  z_basic_block_t &l_6 = cfg->insert("l_6");

  // adding control flow
  entry >> l_1;
  l_1 >> l_2;
  l_2 >> l_1;
  l_2 >> l_3;
  l_3 >> l_4;
  l_3 >> l_5;
  l_4 >> l_6;
  l_5 >> l_6;
  l_6 >> l_1;

  // adding statements
  entry.assign(n, 0);
  l_4.assume(n <= 59);
  l_4.add(n, n, 1);
  l_5.assume(n >= 60);
  l_5.assign(n, 0);

  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {
  // Definining program variables
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop1_header = cfg->insert("loop1_header");
  z_basic_block_t &loop1_bb0 = cfg->insert("loop1_bb0");
  z_basic_block_t &loop1_bb1 = cfg->insert("loop1_bb1");
  z_basic_block_t &loop1_bb1_t = cfg->insert("loop1_bb1_t");
  z_basic_block_t &loop1_bb1_f = cfg->insert("loop1_bb1_f");
  z_basic_block_t &loop1_bb2 = cfg->insert("loop1_bb2");
  z_basic_block_t &loop1_bb3 = cfg->insert("loop1_bb3");
  z_basic_block_t &loop1_bb4 = cfg->insert("loop1_bb4");
  z_basic_block_t &cont = cfg->insert("cont");
  z_basic_block_t &loop2_header = cfg->insert("loop2_header");
  z_basic_block_t &loop2_bb0 = cfg->insert("loop2_bb0");
  z_basic_block_t &loop2_bb1 = cfg->insert("loop2_bb1");
  z_basic_block_t &loop2_bb1_t = cfg->insert("loop2_bb1_t");
  z_basic_block_t &loop2_bb1_f = cfg->insert("loop2_bb1_f");
  z_basic_block_t &loop2_bb2 = cfg->insert("loop2_bb2");
  z_basic_block_t &loop2_bb3 = cfg->insert("loop2_bb3");
  z_basic_block_t &loop2_bb4 = cfg->insert("loop2_bb4");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> loop1_header;

  loop1_header >> loop1_bb0;
  loop1_header >> loop1_bb1;
  loop1_bb0 >> loop1_header;
  loop1_bb1 >> loop1_bb1_t;
  loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2;
  loop1_bb1_f >> loop1_bb3;
  loop1_bb2 >> loop1_bb4;
  loop1_bb3 >> loop1_bb4;
  loop1_bb4 >> loop1_header;
  loop1_header >> cont;

  cont >> loop2_header;
  loop2_header >> loop2_bb0;
  loop2_header >> loop2_bb1;
  loop2_bb0 >> loop2_header;
  loop2_bb1 >> loop2_bb1_t;
  loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2;
  loop2_bb1_f >> loop2_bb3;
  loop2_bb2 >> loop2_bb4;
  loop2_bb3 >> loop2_bb4;
  loop2_bb4 >> loop2_header;
  loop2_header >> ret;

  // adding statements
  entry.assign(n, 0);
  loop1_bb0.havoc(x);
  loop1_bb0.assume(x >= 0);
  loop1_bb1_t.assume(n <= 59);
  loop1_bb1_f.assume(n >= 60);
  loop1_bb2.add(n, n, 1);
  loop1_bb3.assign(n, 0);

  cont.assign(n, 0);
  loop2_bb0.havoc(x);
  loop2_bb0.assume(x >= 10);
  loop2_bb1_t.assume(n <= 160);
  loop2_bb1_f.assume(n >= 161);
  loop2_bb2.add(n, n, 2);
  loop2_bb3.assign(n, 0);

  return cfg;
}

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  // Simple operations on thresholds
  crab::thresholds<ikos::z_number> thresholds;
  using z_bound_t = ikos::bound<ikos::z_number>;
  thresholds.add(z_bound_t(-89));
  thresholds.add(z_bound_t(1000));
  thresholds.add(z_bound_t(100));
  thresholds.add(z_bound_t(-999));
  thresholds.add(z_bound_t(5));
  thresholds.add(z_bound_t(10));
  thresholds.add(z_bound_t(-10));

  crab::outs() << "Thresholds= " << thresholds << "\n";

  auto t1 = thresholds.get_next(z_bound_t(3));
  crab::outs() << "next threshold for 3:" << t1 << "\n";
  auto t2 = thresholds.get_next(z_bound_t(8));
  crab::outs() << "next threshold for 8: " << t2 << "\n";
  auto t3 = thresholds.get_next(z_bound_t(100));
  crab::outs() << "next threshold for 100: " << t3 << "\n";
  auto t4 = thresholds.get_next(z_bound_t(500));
  crab::outs() << "next threshold for 500: " << t4 << "\n";
  auto t5 = thresholds.get_next(z_bound_t(10000));
  crab::outs() << "next threshold for 10000: " << t5 << "\n";
  auto t6 = thresholds.get_next(z_bound_t(-4));
  crab::outs() << "next threshold for -4: " << t6 << "\n";
  auto t7 = thresholds.get_prev(z_bound_t(-4));
  crab::outs() << "prev threshold for -4:" << t7 << "\n";
  auto t8 = thresholds.get_prev(z_bound_t(-78));
  crab::outs() << "prev threshold for -78:" << t8 << "\n";
  auto t9 = thresholds.get_prev(z_bound_t(-10000));
  crab::outs() << "prev threshold for -10000:" << t9 << "\n";

  { // Example where widening with thresholds is better than narrowing
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    crab::outs() << *cfg << "\n";

    {
      ////////////////////////////////////////////////////////////
      // -- w/o thresholds
      ////////////////////////////////////////////////////////////
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 0, stats_enabled);
    }

    {
      ////////////////////////////////////////////////////////////
      // -- w/thresholds
      ////////////////////////////////////////////////////////////
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }

    delete cfg;
  }

  { // Print thresholds per loop
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac);
    crab::outs() << *cfg << "\n";

    // Print thresholds
    using z_cfg_ref_t = crab::cfg::cfg_ref<z_cfg_t>;
    using wto_t = ikos::wto<z_cfg_ref_t>;
    using wto_thresholds_t = crab::wto_thresholds<z_cfg_ref_t>;

    z_cfg_ref_t cfg_ref(*cfg);
    wto_t wto(cfg_ref);
    wto_thresholds_t thresholds(cfg_ref, 50);
    wto.accept(&thresholds);
    crab::outs() << "Thresholds per loop \n" << thresholds << "\n";
  }
}
