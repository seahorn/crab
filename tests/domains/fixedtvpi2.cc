#include "../program_options.hpp"
#include "../common.hpp"


using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

z_cfg_t *prog1(variable_factory_t &vfac) {

  /*
    int x = nd_int();
    int y = 2 * x;
    int z = nd_int();
    int z2 = 2 * z;

    __CRAB_assert(z - x >= 3);   // EXPECTED OK
    __CRAB_assert(y == 2 * x);   // EXPECTED OK
    __CRAB_assert(z2 - y >= 6);  // EXPECTED OK
   */
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var z2(vfac["z2"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header_1 = cfg->insert("header_1");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> header_1;
  header_1 >> exit;
  // adding statements
  entry.havoc(x);
  entry.mul(y, x, 2);
  entry.havoc(z);
  entry.mul(z2, z, 2);
  header_1.assume(z - x >= 3);
  exit.assertion(y == 2 * x);
  exit.assertion(z2 - y >= 6);
  exit.assertion(2*z - y >= 6);

  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {

  /*
    int x = nd_int();
    __CRAB_assume(x >= 0);
    __CRAB_assume(x <= 4);
    int y = nd_int();
    __CRAB_assume(y = 2);

    int z = nd_int();
    z = x * y;

    int m = nd_int();
    __CRAB_assume(m >= z);

    int j = x - 1;
    int k = j * y;

    __CRAB_assert(k + y <= m); // EXPECTED OK
   */
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var m(vfac["m"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb_1 = cfg->insert("bb_1");
  z_basic_block_t &bb_2 = cfg->insert("bb_2");
  z_basic_block_t &bb_3 = cfg->insert("bb_3");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> bb_1;
  bb_1 >> bb_2;
  bb_2 >> bb_3;
  bb_3 >> exit;
  // adding statements
  entry.havoc(x);
  entry.assume(x >= 0);
  entry.assume(x <= 4);
  entry.havoc(y);
  entry.assume(y == 2);

  bb_1.havoc(z);
  bb_1.mul(z, x, y);

  bb_2.havoc(m);
  bb_2.assume(m >= z);

  bb_3.sub(j, x, 1);
  bb_3.mul(k, j, y);
  exit.assertion(k + y <= m);

  return cfg;
}

z_cfg_t *prog3(variable_factory_t &vfac) {

  /*
    int n = nd_int();
    int x = nd_int();
    int z = nd_int();
    int z' = nd_int();
    int m = nd_int();
    int y = 2 * n;

    x = 1 * y;
    z = x;
    m = y / 1;
    z' = m;

    __CRAB_assert(z == 2 * n);   // EXPECTED OK
    __CRAB_assert(z' == 2 * n);   // EXPECTED OK
   */
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var z_p(vfac["z'"], crab::INT_TYPE, 32);
  z_var m(vfac["m"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header_1 = cfg->insert("header_1");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> header_1;
  header_1 >> exit;
  // adding statements
  entry.havoc(n);
  entry.havoc(x);
  entry.havoc(z);
  entry.havoc(m);
  entry.mul(y, n, 2);
  header_1.mul(x, y, 1);
  header_1.assign(z, x);
  header_1.div(m, y, 1);
  header_1.assign(z_p, m);
  exit.assertion(z == 2 * n);
  exit.assertion(z_p == 2 * n);
  return cfg;
}

int main (int argc, char** argv) {
  bool stats_enabled = false;  
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  {
    crab_domain_params_man::get().coefficients().push_back(2);
    crab_domain_params_man::get().coefficients().push_back(3); 
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    z_fixed_tvpi_domain_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    crab_domain_params_man::get().coefficients().clear();
    delete cfg;
  }

  {
    crab_domain_params_man::get().coefficients().push_back(2);
    crab_domain_params_man::get().coefficients().push_back(3); 
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac);
    z_fixed_tvpi_domain_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    crab_domain_params_man::get().coefficients().clear();
    delete cfg;
  }

  {
    crab_domain_params_man::get().coefficients().push_back(2);
    variable_factory_t vfac;
    z_cfg_t *cfg = prog3(vfac);
    z_fixed_tvpi_domain_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    crab_domain_params_man::get().coefficients().clear();
    delete cfg;
  }

  return 0;
}
