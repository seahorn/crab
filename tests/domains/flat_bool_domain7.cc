#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
/*
bb:
  assume(x >= 0)
  b1 = (y = 1);
  b2 = (x >= 0);
  b3 = (x <= 10);
  b4 = true;
  b5 = ite(b2,b3,b4);
  b6 = b1 & b5;
  assume(b6);
exit
  assert(y == 1);
  assert(x >= 0);
  assert(x <= 10);
*/

  variable_factory_t vfac;  
  // entry and exit block
  z_cfg_t cfg("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg.insert("entry");
  z_basic_block_t &exit = cfg.insert("exit");
  // adding control flow
  entry >> exit;
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
  z_var b5(vfac["b5"], crab::BOOL_TYPE, 1);
  z_var b6(vfac["b6"], crab::BOOL_TYPE, 1);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  
  // adding statements
  entry.assume(z_lin_exp_t(x) >= 0);
  entry.bool_assign(b1, z_lin_exp_t(y) == 1);
  entry.bool_assign(b2, z_lin_exp_t(x) >= 0);
  entry.bool_assign(b3, z_lin_exp_t(x) <= 10);
  entry.bool_assign(b4, z_lin_cst_t::get_true());
  entry.bool_select(b5, b2, b3, b4);
  entry.bool_and(b6, b1, b5);
  entry.bool_assume(b6);

  // All assertions are OK
  exit.assertion(z_lin_exp_t(y) >= 1);  
  exit.assertion(z_lin_exp_t(x) >= 0);
  exit.assertion(z_lin_exp_t(x) <= 10);
  z_bool_num_domain_t init;
  run_and_check(&cfg, cfg.entry(), init, false, 1, 2, 20, stats_enabled);
  
  return 0;
  
}
