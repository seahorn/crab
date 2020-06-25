#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog1(variable_factory_t &vfac) {

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 8);
  z_var y(vfac["y"], crab::INT_TYPE, 8);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb_nd = cfg->insert("bb_nd");
  z_basic_block_t &bb_nd_tt = cfg->insert("bb_nd_tt");
  z_basic_block_t &bb_nd_ff = cfg->insert("bb_nd_ff");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb_nd;
  bb_nd >> bb_nd_tt;
  bb_nd >> bb_nd_ff;
  bb_nd_tt >> bb1;
  bb_nd_ff >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(y, z_number(-10));
  entry.havoc(nd);
  bb_nd_tt.assume(nd >= 1);
  bb_nd_tt.assign(x, z_number(0));
  bb_nd_ff.assume(nd <= 0);
  bb_nd_ff.assign(x, z_number(100));
  bb1_t.assume(x >= y);
  bb1_f.assume(x <= y - 1);
  bb2.sub(x, x, y);
  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac, bool is_signed) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 8);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb_if = cfg->insert("if");
  z_basic_block_t &bb_then = cfg->insert("then");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb_if;
  entry >> bb_then;
  bb_if >> ret;
  bb_then >> ret;
  // adding statements
  entry.assign(x, z_number(127));
  entry.add(x, x, 1);
  z_lin_cst_t c1(x <= z_number(1));
  if (is_signed) {
    c1.set_signed();
  } else {
    c1.set_unsigned();
  }
  bb_if.assume(c1);
  bb_if.assign(x, z_number(10));
  z_lin_cst_t c2(x >= z_number(2));
  if (is_signed) {
    c2.set_signed();
  } else {
    c2.set_unsigned();
  }
  bb_then.assume(c2);
  bb_then.assign(x, z_number(-10));
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    crab::outs() << *cfg << "\n";
    {
      // unsound result
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      // sound result
      z_wrapped_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }
  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac, true);
    crab::outs() << *cfg << "\n";
    z_wrapped_interval_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac, false);
    crab::outs() << *cfg << "\n";
    z_wrapped_interval_domain_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  return 0;
}
