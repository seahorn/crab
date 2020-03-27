#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog (variable_factory_t &vfac)  {

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("x0","ret");
  // adding blocks
  z_basic_block_t& x0 = cfg->insert ("x0");
  z_basic_block_t& x1 = cfg->insert ("x1");
  z_basic_block_t& x2 = cfg->insert ("x2");
  z_basic_block_t& x3 = cfg->insert ("x3");
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  x0 >> x1; x1 >> x2; x2 >> x3; x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  x0.assign (k, 2147483648);
  x0.intrinsic("foo1", {i}, {i,k});
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.havoc(nd);
  bb2.select(inc,nd,1,2);
  bb2.add(i, i, inc);
  bb2.intrinsic("foo2", {i}, {i,k});

  return cfg;
}

int main (int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  /* Show that a CFG can have intrinsics */
  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  cfg->simplify (); 
  crab::outs() << *cfg << "\n";
  delete cfg;

  /* For now all domains ignore intrinsics so nothing is printed */  
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var o1(vfac["o1"], crab::INT_TYPE, 32);
  z_var o2(vfac["o2"], crab::INT_TYPE, 32);

  { 
    z_interval_domain_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";
  }

  { 
    z_term_domain_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }

  { 
    z_bool_num_domain_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }

  { 
    z_aa_bool_int_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }

  { 
    domain_product2<z_number, varname_t,
		    z_sdbm_domain_t, z_dis_interval_domain_t> inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }
#ifdef HAVE_APRON
  { 
    z_pk_apron_domain_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }
#endif
#ifdef HAVE_ELINA
  { 
    z_pk_elina_domain_t inv;
    inv.intrinsic("foo1", {o1,o2}, {i1,i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";    
  }
#endif   
  
  return 0;
}
