#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *cfg1(variable_factory_t &vfac, bool add_free) {

/*
==== C program ====
typedef struct node{
  int f;
  int *s;
  struct node *n;
} *List;
#define N 10000
int main() {
  List l = 0;
  int i;
  for (i=0; i<N;i++) {
    List tmp = (List) xmalloc(sizeof(struct node));
    tmp->f = i;
    tmp->s = (int*) xmalloc(sizeof(int)); 
    tmp->n = l ;
    l = tmp;
  }
  List aux = l;
  while (aux) {
    assert(is_not_dangling(aux->s))
    // if add_free is true
    // 
    // NOTE: the region domain cannot prove that aux->s is not
    // dangling if we free inside the loop because the region domain
    // smashes all the memory objects within the region.
    free(aux->s); 
    aux = aux->n;
  }
}
*/

  // Define program variables
  z_var b1(vfac["b1"], crab::BOOL_TYPE);
  z_var b2(vfac["b2"], crab::BOOL_TYPE);  
  z_var i0(vfac["i0"], crab::INT_TYPE, 32);
  z_var x1(vfac["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac["x2"], crab::INT_TYPE, 32);    
  z_var ref0(vfac["ref0"], crab::REF_TYPE, 32);
  z_var ref1(vfac["ref1"], crab::REF_TYPE, 32);
  z_var ref2(vfac["ref2"], crab::REF_TYPE, 32);
  z_var ref3(vfac["ref3"], crab::REF_TYPE, 32);
  z_var ref4(vfac["ref4"], crab::REF_TYPE, 32);
  z_var ref5(vfac["ref5"], crab::REF_TYPE, 32);
  z_var ref6(vfac["ref6"], crab::REF_TYPE, 32);
  z_var ref7(vfac["ref7"], crab::REF_TYPE, 32);
  z_var ref8(vfac["ref8"], crab::REF_TYPE, 32);
  z_var ref9(vfac["ref9"], crab::REF_TYPE, 32);
  z_var ref10(vfac["ref10"], crab::REF_TYPE, 32);  
  
  // Define memory regions
  z_var mem_field_f(vfac["region_field_f"], crab::REG_INT_TYPE, 32);
  z_var mem_field_s(vfac["region_field_s"], crab::REG_REF_TYPE, 32);
  z_var mem_field_deref_s(vfac["region_field_deref_s"], crab::REG_INT_TYPE, 32);    
  z_var mem_field_next(vfac["region_field_next"], crab::REG_REF_TYPE, 32);

  // Create allocation sites
  crab::tag_manager as_man;
  
  // Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry", "exit");
  // Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb2_loop = cfg->insert("bb2_loop");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  z_basic_block_t &bb5 = cfg->insert("bb5");
  z_basic_block_t &bb6_loop = cfg->insert("bb6_loop");
  z_basic_block_t &bb7 = cfg->insert("bb7");
  z_basic_block_t &bb8 = cfg->insert("bb8");
  z_basic_block_t &bb9 = cfg->insert("bb9");
  z_basic_block_t &bb9_assert = cfg->insert("bb9_assert");        
  z_basic_block_t &exit = cfg->insert("exit");
  // Adding CFG edges
  entry.add_succ(bb2_loop);
  bb2_loop.add_succ(bb3);
  bb2_loop.add_succ(bb4);
  bb3.add_succ(bb5);
  bb5.add_succ(bb2_loop);  
  bb4.add_succ(bb6_loop);
  bb6_loop.add_succ(bb7);
  bb6_loop.add_succ(bb8);
  bb7.add_succ(exit);
  bb8.add_succ(bb9);
  bb9.add_succ(bb9_assert);
  bb9_assert.add_succ(bb6_loop);  

  z_var_or_cst_t n30000_32(z_number(30000), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size16(z_number(16), crab::variable_type(crab::INT_TYPE, 32));
  
  entry.region_init(mem_field_next);
  entry.region_init(mem_field_s);
  entry.region_init(mem_field_deref_s);  
  entry.region_init(mem_field_f);
  // i=0;
  entry.assign(i0, 0);
  // List l = 0;
  entry.assume_ref(z_ref_cst_t::mk_null(ref0));
  // i < N
  bb3.assume(i0 <= 9999);
  // i >= N
  bb4.assume(i0 >= 10000);
  // List aux = l
  bb4.gep_ref(ref1, mem_field_f, ref0, mem_field_f);
  // tmp = malloc(...)
  bb5.make_ref(ref2, mem_field_f, size16, as_man.mk_tag());
  bb5.assume_ref(z_ref_cst_t::mk_gt_null(ref2));
  /// tmp->f = i  
  bb5.gep_ref(ref3, mem_field_f, ref2, mem_field_f);
  bb5.store_to_ref(ref3, mem_field_f, i0);
  // tmp->s = xmalloc(...)
  bb5.gep_ref(ref5, mem_field_s, ref2, mem_field_f, 4);
  bb5.make_ref(ref10, mem_field_deref_s, size4, as_man.mk_tag()); 
  bb5.assume_ref(z_ref_cst_t::mk_gt_null(ref10));
  bb5.store_to_ref(ref5, mem_field_s, ref10);
  // tmp->n = l
  bb5.gep_ref(ref4, mem_field_next, ref2, mem_field_f, 8);
  bb5.store_to_ref(ref4, mem_field_next, ref0);
  // i=i+1
  bb5.add(i0, i0, 1);
  // l = tmp;
  bb5.gep_ref(ref0, mem_field_f, ref2, mem_field_f);
  bb7.assume_ref(z_ref_cst_t::mk_null(ref1));
  // while(aux)
  bb8.assume_ref(z_ref_cst_t::mk_not_null(ref1));
  // __CRAB_assert(is_not_dangling(aux->s))
  bb9.gep_ref(ref7, mem_field_s, ref1, mem_field_f, 4);
  bb9.load_from_ref(ref6, ref7, mem_field_s);
  bb9.intrinsic("is_unfreed_or_null",{b2},{mem_field_deref_s, ref6});
  if (add_free) 
    bb9.remove_ref(mem_field_deref_s, ref6);
  bb9_assert.bool_assert(b2);
  // aux = aux->n;
  bb9_assert.gep_ref(ref8, mem_field_next, ref1, mem_field_f, 8);
  bb9_assert.load_from_ref(ref9, ref8, mem_field_next);
  bb9_assert.gep_ref(ref1, mem_field_f, ref9, mem_field_f);
  
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
  {
    z_cfg_t *p1 = cfg1(vfac, false);
    crab::outs() << *p1 << "\n";
    z_rgn_bool_int_t init;
    run_and_check(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
    //run(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
    // EXPECTED TO PROVE ASSERTION 
    delete p1;
  }

  {
    z_cfg_t *p1 = cfg1(vfac, true);
    crab::outs() << *p1 << "\n";
    z_rgn_bool_int_t init;
    run_and_check(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
    //run(p1, p1->entry(), init, false, 2, 2, 20, stats_enabled);
    // NOT EXPECTED TO PROVE ASSERTION
    delete p1;
  }
  

  return 0;
}
