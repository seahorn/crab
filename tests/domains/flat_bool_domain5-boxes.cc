#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* test boolean-to-non-boolean propagation through select */

// test propagation through select
z_cfg_t *prog2(variable_factory_t &vfac) {

  /*
    havoc(x) 
    b2 = (-x <= 0)
    b3 = (x <= 10)
    b4 = false
    b5 = ite(b2,b3,b4)
    assume(b5)

    assert(x >= 0);
    assert(x <= 10);

    EXPECTED RESULT: SAFE

    This program is safe because b5 is true and b4 is false so b2 must
    be true. However, our propagation is not that strong. Note that
    the boxes domain can prove this program.
   */
  
  // Defining program variables
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
  z_var b5(vfac["b5"], crab::BOOL_TYPE, 1);  
  z_var x(vfac["x"], crab::INT_TYPE, 64);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry.add_succ(exit);
  // adding statements
  entry.havoc(x);
  entry.bool_assign(b2, x >= z_number(0));
  entry.bool_assign(b3, x <= z_number(10));
  entry.bool_assign(b4, z_lin_cst_t::get_false());
  entry.bool_select(b5,b2,b3,b4);
  entry.bool_assume(b5);
  exit.assertion(x >= z_number(0));  
  exit.assertion(x <= z_number(10));
  return cfg;
}

// test propagation through select
z_cfg_t *prog3(variable_factory_t &vfac) {

  /*
    havoc(x) 
    b2 = (-x <= 0)
    b3 = (x <= 10)
    b4 = false
    b5 = ite(b2,b3,b4)
    assume_not(b5)

    assert(x <= 10);

    EXPECTED RESULT: UNSAFE

    The program is not safe. Because b5 is false then either
    1) if b2 is false (x < 0) then b4=false
    2) if b2 is true then b3=false that implies x > 10
   */
  
  // Defining program variables
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
  z_var b5(vfac["b5"], crab::BOOL_TYPE, 1);  
  z_var x(vfac["x"], crab::INT_TYPE, 64);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry.add_succ(exit);
  // adding statements
  entry.havoc(x);
  entry.bool_assign(b2, x >= z_number(0));
  entry.bool_assign(b3, x <= z_number(10));
  entry.bool_assign(b4, z_lin_cst_t::get_false());
  entry.bool_select(b5,b2,b3,b4);
  entry.bool_not_assume(b5);
  exit.assertion(x <= z_number(10));
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
#ifdef HAVE_LDD  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  {
    z_cfg_t *cfg = prog2(vfac);
    crab::outs() << *cfg << "\n";
    // Boxes can prove this program while bool_interval_domain cannot
    z_boxes_domain_t boxes_init;
    run_and_check(cfg, cfg->entry(), boxes_init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }

  {
    z_cfg_t *cfg = prog3(vfac);
    crab::outs() << *cfg << "\n";
    z_boxes_domain_t boxes_init;
    run_and_check(cfg, cfg->entry(), boxes_init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
#endif     
  return 0;
}
