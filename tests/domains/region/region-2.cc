#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *cfg1(variable_factory_t &vfac) {

  /*
   int v1,v2,v3;
   int *i = &v1;
   int *x = &v2;
   int *y = &v3;

   assume(i != NULL);
   assume(x != NULL);
   assume(y != NULL);

   *i := 0;
   *x := 1;
   *y := 0;

    while(*i <= 99) {
     *x+=*y;
     *y++;
     *i++;
    }

    assert(x != NULL);
    assert(y != NULL);
    assert(*x >= *y);
   */

  // Define program variables
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var y(vfac["y"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var deref_x(vfac["*x"], crab::INT_TYPE, 32);
  z_var deref_y(vfac["*y"], crab::INT_TYPE, 32);
  // Define memory regions
  z_var mem1(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var mem2(vfac["region_1"], crab::REG_INT_TYPE, 32);
  z_var mem3(vfac["region_2"], crab::REG_INT_TYPE, 32);
  // === Create allocation sites
  crab::tag_manager as_man;  
  // Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // Adding CFG edges
  entry.add_succ(bb1);
  bb1.add_succ(bb1_t);
  bb1.add_succ(bb1_f);
  bb1_t.add_succ(bb2);
  bb2.add_succ(bb1);
  bb1_f.add_succ(ret);

  // === adding statements

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
 
  // Intialization of memory regions
  entry.region_init(mem1);
  entry.region_init(mem2);
  entry.region_init(mem3);

  //// Create references
  entry.make_ref(i, mem1, size4, as_man.mk_tag());
  entry.make_ref(x, mem2, size4, as_man.mk_tag());
  entry.make_ref(y, mem3, size4, as_man.mk_tag());
  #if 0
  // Disequalities too hard
  entry.assume_ref(z_ref_cst_t::mk_not_null(i));
  #else
  entry.assume_ref(z_ref_cst_t::mk_gt_null(i));
  #endif
  #if 0
  // Disequalities too hard
  entry.assume_ref(z_ref_cst_t::mk_not_null(x));
  #else
  entry.assume_ref(z_ref_cst_t::mk_gt_null(x));
  #endif 
  #if 0
  // Disequalities too hard
  entry.assume_ref(z_ref_cst_t::mk_not_null(y));
  #else
  entry.assume_ref(z_ref_cst_t::mk_gt_null(y));
  #endif 
  //// *i := 0;
  entry.store_to_ref(i, mem1, zero32);
  //// *x := 1;
  entry.store_to_ref(x, mem2, one32);
  //// *y := 0;
  entry.store_to_ref(y, mem3, zero32);
  //// assume(*i <= 99);
  bb1_t.load_from_ref(deref_i, i, mem1);
  bb1_t.assume(deref_i <= 99);
  //// assume(*i >= 100);
  bb1_f.load_from_ref(deref_i, i, mem1);
  bb1_f.assume(deref_i >= 100);
  //// *x = *x + *y
  bb2.load_from_ref(deref_x, x, mem2);
  bb2.load_from_ref(deref_y, y, mem3);
  bb2.add(deref_x, deref_x, deref_y);
  bb2.store_to_ref(x, mem2, deref_x);
  //// *y = *y + 1
  bb2.load_from_ref(deref_y, y, mem3);
  bb2.add(deref_y, deref_y, 1);
  bb2.store_to_ref(y, mem3, deref_y);
  //// *i = *i + 1
  bb2.load_from_ref(deref_i, i, mem1);
  bb2.add(deref_i, deref_i, 1);
  bb2.store_to_ref(i, mem1, deref_i);
  //// assume(*x <= *y)
  ret.assert_ref(z_ref_cst_t::mk_not_null(x));
  ret.assert_ref(z_ref_cst_t::mk_not_null(y));
  ret.load_from_ref(deref_x, x, mem2);
  ret.load_from_ref(deref_y, y, mem3);
  ret.assertion(deref_x >= deref_y);

  return cfg;
}

int main(int argc, char **argv) {
  region_domain_params p(true/*allocation_sites*/,
			 true/*deallocation*/,
			 true/*tag_analysis*/,
			 false/*is_dereferenceable*/,
			 true/*skip_unknown_regions*/);
  crab_domain_params_man::get().update_params(p);

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_cfg_t *p1 = cfg1(vfac);
  crab::outs() << *p1 << "\n";
  z_rgn_sdbm_t init;
  run_and_check(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
  delete p1;

  return 0;
}
