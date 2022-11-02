#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/*
 * Crab distinguishes between integer and booleans.
 * Example of program with booleans
 */

z_cfg_t *prog(variable_factory_t &vfac) {

  /*
        i := 0;
        b_true  := true;
        b_false := false;
        b1 := (i <= 99);
        while (b1) {
          assert(b1);
          i  += ( * ? 1: 2);
          b3 := (i >= 1);
          b4 := b3
          // b4 won't change after these four statements
          b4 := b4 or b_false
          b4 := b4 and b_true
          b4 := b4 xor b_true;
          b4 := b4 xor b_true;
          b1 := (i <= 99);
        }
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
  z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
  z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
  z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
  z_var bfalse(vfac["bf"], crab::BOOL_TYPE, 1);
  z_var btrue(vfac["bt"], crab::BOOL_TYPE, 1);

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(i, z_number(0));
  entry.bool_assign(bfalse, z_lin_cst_t::get_false());
  entry.bool_assign(btrue, z_lin_cst_t::get_true());
  bb1_t.bool_assign(b1, i <= z_number(99));
  bb1_t.bool_assume(b1);
  bb1_t.bool_assert(b1); // trivial
  bb1_f.bool_assign(b2, i >= z_number(100));
  bb1_f.bool_assume(b2);
  bb2.havoc(nd);
  bb2.select(inc, nd, 1, 2);
  bb2.add(i, i, inc);
  bb2.bool_assign(b3, i >= z_number(1));
  bb2.bool_assign(b4, b3);
  bb2.bool_or(b4, b4, bfalse); // tautology
  bb2.bool_and(b4, b4, btrue); // tautology
  bb2.bool_xor(b4, b4, btrue); // complement
  bb2.bool_xor(b4, b4, btrue); // complement

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
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  z_boxes_domain_t init;
  run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  delete cfg;
#endif  
  return 0;
}
