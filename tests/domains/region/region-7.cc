#include "../../common.hpp"
#include "../../program_options.hpp"

#include <crab/analysis/fwd_analyzer.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *cfg1(variable_factory_t &vfac) {
  /*
    void (*fun_ptr)(int) = fun3;

    int main() {

      int x = nd();
      assume (x >= 0);
      if (x > 0)
        fun_ptr = fun1;
      else
        fun_ptr = fun2;

      if (x < 0) fun_ptr = fun3;

      sassert(fun_ptr != fun3);
      (*fun_ptr)(x);
      return 0;
    }
  */
  // Define program variables
  z_var f1(vfac["fun1"], crab::REF_TYPE);
  z_var f2(vfac["fun2"], crab::REF_TYPE);
  z_var f3(vfac["fun3"], crab::REF_TYPE);
  z_var f_ptr(vfac["fun_ptr"], crab::REF_TYPE);
  z_var res(vfac["res"], crab::REF_TYPE);
  z_var orphan_ptr(vfac["orphan_ptr"], crab::REF_TYPE);  
  z_var x(vfac["x"], crab::INT_TYPE, 32);    
  z_var m1(vfac["rgn_0"], crab::REG_REF_TYPE, 32);
  z_var m2(vfac["rgn_1"], crab::REG_REF_TYPE, 32);
  z_var m3(vfac["rgn_2"], crab::REG_REF_TYPE, 32);
  z_var m4(vfac["rgn_3"], crab::REG_REF_TYPE, 32);  
  // === Create allocation sites
  crab::tag_manager as_man;  
  // Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  z_basic_block_t &bb5 = cfg->insert("bb5");
  z_basic_block_t &bb6 = cfg->insert("bb6");    
  z_basic_block_t &ret = cfg->insert("ret");
  // Adding CFG edges
  entry.add_succ(bb1);
  bb1.add_succ(bb2);
  bb1.add_succ(bb3);
  bb2.add_succ(bb4);
  bb3.add_succ(bb4);
  bb4.add_succ(bb5);
  bb4.add_succ(bb6);
  bb5.add_succ(bb6);
  bb6.add_succ(ret);

  // === adding statements

  // Intialization of memory regions
  entry.region_init(m1);
  entry.region_init(m2);
  entry.region_init(m3);
  entry.region_init(m4);  

  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  
  // Create references
  entry.make_ref(f1, m1, size8, as_man.mk_tag());
  entry.make_ref(f2, m2, size8, as_man.mk_tag());
  entry.make_ref(f3, m3, size8, as_man.mk_tag());
  entry.make_ref(f_ptr, m4, size8, as_man.mk_tag());
  entry.havoc(orphan_ptr);
  entry.havoc(x);
  bb1.assume(x >= 0);
  bb2.assume(x > 0);
  bb2.store_to_ref(f_ptr, m4, f1); 
  bb3.assume(x <= 0);
  bb3.store_to_ref(f_ptr, m4, f2);  
  bb5.assume(x < 0);
  bb5.store_to_ref(f_ptr, m4, f3);
  bb6.load_from_ref(res, f_ptr, m4);
  // true but cannot be proven by simply comparing addresses because
  // the region domain is not strong enough for that. However, the
  // region domain can tell that &f3 is not an allocation site of res.
  bb6.assert_ref(z_ref_cst_t::mk_not_eq(res, f3)); 
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
  z_rgn_sdbm_t absval_fac, init;
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = 2;
  fixpo_params.get_descending_iterations() = 2;
  fixpo_params.get_max_thresholds() = 20;
  
  using intra_fwd_analyzer_t =
    crab::analyzer::intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, z_rgn_sdbm_t>;
  intra_fwd_analyzer_t a(*p1, absval_fac, nullptr, fixpo_params);
  a.run(p1->entry(), init, typename intra_fwd_analyzer_t::assumption_map_t());
  
  /* 
   * Example of how to ask an abstract domain about the allocation
   * sites associated to a program reference
   *  
   * We ask about the allocation sites of *fun_ptr at the exit of the
   * function
   */
  
  z_var f_ptr(vfac["fun_ptr"], crab::REF_TYPE);
  z_var m(vfac["rgn_3"], crab::REG_REF_TYPE, 32);  // region of f_ptr
  z_var deref_f_ptr(vfac["*fun_ptr"], crab::REF_TYPE);
  z_var orphan_ptr(vfac["orphan_ptr"], crab::REF_TYPE);    
  z_rgn_sdbm_t exit_inv = a["ret"];

  auto print_alloc_sites = [](const z_var &ref, z_rgn_sdbm_t &inv) {
			     crab::outs() << "Allocation sites for " << ref << " at the exit: ";    
			     std::vector<crab::allocation_site> alloc_sites;
			     if (inv.get_allocation_sites(ref, alloc_sites)) {
			       crab::outs() << "{";
			       for (unsigned i=0,sz=alloc_sites.size();i<sz;) {
				 crab::outs() << alloc_sites[i];
				 ++i;
				 if (i < sz) {
				   crab::outs() << ",";
				 }
			       }
			       crab::outs() << "}\n";
			     } else {
			       crab::outs() << "unknown\n";
			     }
			   };

  exit_inv.ref_load(f_ptr, m, deref_f_ptr);  
  print_alloc_sites(deref_f_ptr, exit_inv);
  print_alloc_sites(orphan_ptr, exit_inv);
  delete p1;
  return 0;
}
