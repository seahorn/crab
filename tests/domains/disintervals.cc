#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
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
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    // using z_dis_interval_t = dis_interval <z_number>;
    // using interval_t = interval <z_number>;

    // {
    //   z_dis_interval_t x (interval_t (0,5));
    //   crab::outs() << x << "\n";
    //   z_dis_interval_t y (interval_t (10,20));
    //   crab::outs() << y << "\n";
    //   z_dis_interval_t z (interval_t (7,8));
    //   crab::outs() << z << "\n";

    //   z_dis_interval_t r = (x | y) | z;
    //   crab::outs() << r << "\n";

    //   z_dis_interval_t s (interval_t (2,3));
    //   crab::outs() << "after adding " << s << "=" ;
    //   r = r + s;
    //   crab::outs() << r << "\n";
    // }

    // {
    //   interval_t zero (5,5);
    //   z_dis_interval_t x (zero.lower_half_line ());
    //   crab::outs() << x << "\n";
    //   z_dis_interval_t y (zero.upper_half_line ());
    //   crab::outs() << y << "\n";
    //   z_dis_interval_t r = (x | y);
    //   crab::outs() << r << "\n";
    // }

    // {
    //   z_dis_interval_t x (interval_t (0,5));
    //   crab::outs() << x << "\n";
    //   z_dis_interval_t y (interval_t (10,20));
    //   crab::outs() << y << "\n";
    //   z_dis_interval_t z (interval_t (7,8));
    //   crab::outs() << z << "\n";

    //   z_dis_interval_t r = (x & y) & z;
    //   crab::outs() << r << "\n";
    // }

    // {
    //   z_dis_interval_t x (interval_t (0,5));
    //   crab::outs() << x << "\n";
    //   z_dis_interval_t y (interval_t (2,9));
    //   crab::outs() << y << "\n";
    //   z_dis_interval_t z (interval_t (2,8));
    //   crab::outs() << z << "\n";

    //   z_dis_interval_t r = (x & y) & z;
    //   crab::outs() << r << "\n";
    // }

    {
      variable_factory_t vfac;
      auto cfg = prog(vfac);
      crab::outs() << *cfg << "\n";
      z_dis_interval_domain_t init;
      run(cfg, init, stats_enabled);
    }
    return 0;
  });
}
