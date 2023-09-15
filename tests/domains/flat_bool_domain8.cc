#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog1(variable_factory_t &vfac) {
/*
entry:
  havoc(x)
  havoc(y)
  b1 = x <= y
  assume(b1)
  assume(y >= 0)
  goto bb1, bb2
bb1:
  assume x >= 0
  b2 = item(y <= -1, true, b1)
  goto exit
bb2
  assume x <= -1
  b2 = item(y <= -1, b1, false)
  goto exit
exit:
  assert(b2) // EXPECTED FALSE (e.g., x=-1, y=1)
*/

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  // adding control flow
  entry >> bb1;
  entry >> bb2;
  bb1 >> exit;
  bb2 >> exit;
  
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var bTrue(vfac["bTrue"], crab::BOOL_TYPE, 1);
  z_var bFalse(vfac["bFalse"], crab::BOOL_TYPE, 1);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  // adding statements
  entry.havoc(x);
  entry.havoc(y);
  entry.bool_assign(b1, z_lin_exp_t(x) <=  z_lin_exp_t(y));
  entry.bool_assume(b1);
  entry.assume(z_lin_exp_t(y) >= 0);
  entry.bool_assign(b3, z_lin_exp_t(y) <= -1);
  entry.bool_assign(bTrue, z_lin_cst_t::get_true());
  entry.bool_assign(bFalse, z_lin_cst_t::get_false());
  bb1.assume(z_lin_exp_t(x) >= 0);
  bb1.bool_select(b2, b3, bTrue, b1);
  bb2.assume(z_lin_exp_t(x) <= -1);
  bb2.bool_select(b2, b3, b1, bFalse);
  // The assertion is NOT provable
  exit.bool_assert(b2);

  return cfg;
}


int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  {
    z_cfg_t *cfg = prog1(vfac);
    //crab::outs() << *cfg << "\n";
    z_bool_num_domain_t init;
    run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  
  return 0;
  
}
