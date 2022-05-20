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

bb1:
    havoc(b1) 
    havoc(b2) 
    b3 = b1&_b2;
    assume(b3);
    b4 = b1&_b3;
    goto bb2, bb4;
bb2:
    assume(b4);
    goto bb3;
bb4:
    assume(not(b4));
    assume(false);
    goto bb3;
bb3:
    assert(b2)   
    assert(b4)  
*/

  variable_factory_t vfac;  
  // entry and exit block
  z_cfg_t cfg("bb1", "bb3");
  // adding blocks
  z_basic_block_t &bb1 = cfg.insert("bb1");
  z_basic_block_t &bb2 = cfg.insert("bb2");
  z_basic_block_t &bb3 = cfg.insert("bb3");
  z_basic_block_t &bb4 = cfg.insert("bb4");
  // adding control flow
  bb1 >> bb2;
  bb1 >> bb4;
  bb2 >> bb3;
  bb4 >> bb3;
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
  // adding statements
  bb1.havoc(b1);
  bb1.havoc(b2);
  bb1.bool_and(b3, b1, b2);
  bb1.bool_assume(b3);
  bb1.bool_and(b4, b1, b3);
  /////
  bb2.bool_assume(b4);
  /////
  bb4.bool_not_assume(b4);
  bb4.assume(z_lin_cst_t::get_false());
  /////
  bb3.bool_assert(b2);
  bb3.bool_assert(b4);

  z_bool_num_domain_t init;
  run_and_check(&cfg, cfg.entry(), init, false, 1, 2, 20, stats_enabled);  
  
  return 0;
  
}
