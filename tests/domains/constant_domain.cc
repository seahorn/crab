#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog1(variable_factory_t &vfac, crab::tag_manager &as_man) {

  /*
    i := 0;
    while (*) {
      i += (* ? 1: 1);
      if (i > 0)
        i=0;
    }
    assert(i==0);
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");  
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb3;
  bb2 >> bb4;
  bb3 >> bb1;
  bb4 >> bb1;  
  bb1_f >> ret;
  // adding statements
  entry.assign(i, 0);
  bb2.havoc(nd);
  bb2.select(inc, nd, 1, 1);
  bb2.add(i, i, inc);
  bb3.assume(i >= 1);
  bb3.assign(i, 0);
  bb4.assume(i <= 0);
  ret.assertion(i == 0);
  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac, crab::tag_manager &as_man) {

  /*
   *i := 0;
    while (*) {
      *i += (* ? 1: 1);
      if (*i > 0) {
         *i=0;
      }
    }
    assert(*i==0);
   */
  // Defining program variables
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var mem1(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));  

  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");  
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb3;
  bb2 >> bb4;
  bb3 >> bb1;
  bb4 >> bb1;  
  bb1_f >> ret;

  // adding statements  
  entry.region_init(mem1);
  entry.make_ref(i, mem1, size4, as_man.mk_tag());
  entry.store_to_ref(i, mem1, zero32);
  bb2.havoc(nd);
  bb2.select(inc, nd, 1, 1);
  bb2.load_from_ref(deref_i, i, mem1);  
  bb2.add(deref_i, deref_i, inc);
  bb2.store_to_ref(i, mem1, deref_i);
  bb3.load_from_ref(deref_i, i, mem1);  
  bb3.assume(deref_i >= 1);
  bb3.store_to_ref(i, mem1, zero32);
  bb4.load_from_ref(deref_i, i, mem1);  
  bb4.assume(deref_i <= 0);
  ret.load_from_ref(deref_i, i, mem1);    
  ret.assertion(deref_i == 0);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
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

  {
    variable_factory_t vfac;
    crab::tag_manager as_man;
    z_cfg_t *cfg = prog1(vfac, as_man);
    crab::outs() << *cfg << "\n";
  
    z_constant_domain_t init;
    run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }

  {
    variable_factory_t vfac;
    crab::tag_manager as_man;    
    z_cfg_t *cfg = prog2(vfac, as_man);
    crab::outs() << *cfg << "\n";
    z_rgn_constant_t init;  
    run_and_check(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }

#if 0  
  {
    using constant_t = crab::domains::constant<z_number>;
    constant_t c1(z_number(5));
    constant_t c2(z_number(10));
    constant_t c3(z_number(-10));
    constant_t one(z_number(1));    
    constant_t two(z_number(2));
    
    crab::outs() << c1 << " + "  << c2 << "=" << c1.Add(c2) << "\n";
    crab::outs() << c1 << " - "  << c2 << "=" << c1.Sub(c2) << "\n";
    crab::outs() << c1 << " * "  << c2 << "=" << c1.Mul(c2) << "\n";
    crab::outs() << c2 << " / "  << c1 << "=" << c2.SDiv(c1) << "\n";
    crab::outs() << c1 << " rem " << c2 << "=" << c1.SRem(c2) << "\n";
    crab::outs() << c1 << " | "  << c2 << "=" << c1.BitwiseOr(c2) << " (expected 15) \n";
    crab::outs() << c1 << " & "  << c2 << "=" << c1.BitwiseAnd(c2) << " (expected 0)\n";
    crab::outs() << c2 << " << "    << two << "=" << c2.BitwiseShl(two) << " (expected 40)\n";
    crab::outs() << c3 << " >>_l "  << one << "=" << c3.BitwiseLShr(one) << " (expected top)\n";
    crab::outs() << c3 << " >>_a "  << one << "=" << c3.BitwiseAShr(one) << " (expected -5)\n";        
    
  }
#endif   
  return 0;
}
