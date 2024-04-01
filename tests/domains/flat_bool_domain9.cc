#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/*
 * Example that shows how the flat boolean domain deals with sext/zext
 */

z_cfg_t *prog(variable_factory_t &vfac) {

  /* 
        x  := 2;
        b1 := true;
	b2 := false;
	b3 := havoc();
	x1 := sext b1 to i64
	x2 := sext b2 to i64
	x3 := sext b3 to i64
	x4 := zext b1 to i64
	x5 := zext b2 to i64
	x6 := zext b3 to i64
	
	y1 := x + x1
	y2 := x + x2
	y3 := x + x3
	y4 := x + x4
	y5 := x + x5
	y6 := x + x6

	assert(y1 == 1); // OK
	assert(y2 == 2); // OK 
	assert(y3 == 2); // FAIL
	assert(y4 == 3); // OK
	assert(y5 == 2); // OK 
	assert(y6 == 2); // FAIL
  */ 
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var x1(vfac["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac["x2"], crab::INT_TYPE, 32);
  z_var x3(vfac["x3"], crab::INT_TYPE, 32);
  z_var x4(vfac["x4"], crab::INT_TYPE, 32);
  z_var x5(vfac["x5"], crab::INT_TYPE, 32);
  z_var x6(vfac["x6"], crab::INT_TYPE, 32);
  z_var y1(vfac["y1"], crab::INT_TYPE, 32);
  z_var y2(vfac["y2"], crab::INT_TYPE, 32);
  z_var y3(vfac["y3"], crab::INT_TYPE, 32);
  z_var y4(vfac["y4"], crab::INT_TYPE, 32);
  z_var y5(vfac["y5"], crab::INT_TYPE, 32);
  z_var y6(vfac["y6"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, z_number(2));
  entry.bool_assign(b1, z_lin_cst_t::get_true());
  entry.bool_assign(b2, z_lin_cst_t::get_false());
  entry.havoc(b3);
  entry.sext(b1, x1);
  entry.sext(b2, x2);
  entry.sext(b3, x3);
  entry.zext(b1, x4);
  entry.zext(b2, x5);
  entry.zext(b3, x6);

  entry.add(y1, x, x1);
  entry.add(y2, x, x2);
  entry.add(y3, x, x3);
  entry.add(y4, x, x4);
  entry.add(y5, x, x5);
  entry.add(y6, x, x6);
  
  entry.assertion(y1 == z_number(1)); // OK
  entry.assertion(y2 == z_number(2)); // OK
  entry.assertion(y3 == z_number(2)); // FAIL

  entry.assertion(y4 == z_number(3)); // OK
  entry.assertion(y5 == z_number(2)); // OK
  entry.assertion(y6 == z_number(2)); // FAIL    
  
  return cfg;
}


/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  z_bool_interval_domain_t init;
  run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);

  delete cfg;
  return 0;
}
