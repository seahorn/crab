#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/fwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog1(variable_factory_t &vfac)  {
  // Definining program variables
  z_var i(vfac ["i"], crab::INT_TYPE, 32);
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var x1(vfac ["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac ["x2"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  //  entry.assign(x1, 1);
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);
  //bb2.add(x2, x1, 1);
  bb2.add(k, k, 1);
  return cfg;
}

z_cfg_t* prog2(variable_factory_t &vfac) {
  z_cfg_t* cfg = new z_cfg_t("loop1_entry","ret");
  z_basic_block_t& loop1_entry = cfg->insert("loop1_entry");
  z_basic_block_t& loop1_bb1   = cfg->insert("loop1_bb1");
  z_basic_block_t& loop1_bb1_t = cfg->insert("loop1_bb1_t");
  z_basic_block_t& loop1_bb1_f = cfg->insert("loop1_bb1_f");
  z_basic_block_t& loop1_bb2   = cfg->insert("loop1_bb2");
  z_basic_block_t& loop2_entry = cfg->insert("loop2_entry");
  z_basic_block_t& loop2_bb1   = cfg->insert("loop2_bb1");
  z_basic_block_t& loop2_bb1_t = cfg->insert("loop2_bb1_t");
  z_basic_block_t& loop2_bb1_f = cfg->insert("loop2_bb1_f");
  z_basic_block_t& loop2_bb2   = cfg->insert("loop2_bb2");
  z_basic_block_t& ret         = cfg->insert("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  loop1_entry.assign(i, 0);
  loop1_entry.assign(k, 30);
  loop1_bb1_t.assume(i <= 9);
  loop1_bb1_f.assume(i >= 10);
  loop1_bb2.add(i, i, 1);

  loop2_entry.assign(j, 0);
  loop2_bb1_t.assume(j <= 9);
  loop2_bb1_f.assume(j >= 10);
  loop2_bb2.add(j, j, 1);
  return cfg;
}

template<class Dom1, class Dom2>
void run_and_compare(crab::cfg_impl::z_cfg_t* cfg,
		     crab::cfg_impl::basic_block_label_t entry,
		     unsigned widening, 
		     unsigned narrowing, 
		     unsigned jump_set_size,
		     bool enable_stats){

  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom1> analyzer_1_t;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom2> analyzer_2_t;  
  typename analyzer_1_t::assumption_map_t assumptions1;
  typename analyzer_2_t::assumption_map_t assumptions2;

  typedef crab::cfg::cfg_ref<crab::cfg_impl::z_cfg_t> cfg_ref_t;

  
  // Run with Dom1
  Dom1 inv1;
  analyzer_1_t a1(*cfg, inv1,
		  nullptr, nullptr,
		  widening, narrowing, jump_set_size);
  a1.run(entry, assumptions1);
  
  // Run with Dom2
  Dom2 inv2;        
  analyzer_2_t a2(*cfg, inv2,
		  nullptr, nullptr,
		  widening, narrowing, jump_set_size);
  a2.run(entry, assumptions2);

  for (auto &bb: *cfg) {
    Dom1 inv1 = a1.get_pre(bb.label());
    Dom2 inv2 = a2.get_pre(bb.label());
    #ifdef HAVE_APRON
    // use apron as baseline
    z_oct_apron_domain_t baseline1, baseline2;
    #else
    // run_and_compare is called if HAVE_APRON or HAVE_ELINA
    // use elina as baseline    
    z_oct_elina_domain_t baseline1, baseline2;
    #endif 
    
    baseline1 += inv1.to_linear_constraint_system();
    baseline2 += inv2.to_linear_constraint_system();
    bool error = false;
    if (!(baseline1 <= baseline2)) {
      crab::outs() << "FAIL. Different invariant at block " << bb << ". ";
      crab::outs() << inv1.domain_name() << " is not included in "
		   << inv2.domain_name() << ":\n";
      crab::outs() << inv1.domain_name() << ":" << inv1 << "\n";
      crab::outs() << inv2.domain_name() << ":" << inv2 << "\n";
      error = true;
    }
    if (!(baseline2 <= baseline1)) {
      crab::outs() << "FAIL. Different invariant at block " << bb << ". ";
      crab::outs() << inv2.domain_name() << " is not included in "
		   << inv1.domain_name() << ":\n";
      crab::outs() << inv1.domain_name() << ":" << inv1 << "\n";
      crab::outs() << inv2.domain_name() << ":" << inv2 << "\n";
      error = true;      
    }
    if (error) {
      crab::outs() << "ERROR: Domains " << inv1.domain_name() << " and "
		   << inv2.domain_name() << " disagree agree on block " 
		   << bb.label () << "\n";
    }
    
  }
}

z_cfg_t* prog3(variable_factory_t &vfac) {
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  z_basic_block_t& entry       = cfg->insert("entry");
  z_basic_block_t& loop1_head  = cfg->insert("loop1_head");
  z_basic_block_t& loop1_t     = cfg->insert("loop1_t");
  z_basic_block_t& loop1_f     = cfg->insert("loop1_f");
  z_basic_block_t& loop1_body  = cfg->insert("loop1_body");

  z_basic_block_t& loop1_body_t  = cfg->insert("loop1_body_t");
  z_basic_block_t& loop1_body_f  = cfg->insert("loop1_body_f");
  z_basic_block_t& loop1_body_x  = cfg->insert("loop1_body_x");

  z_basic_block_t& cont        = cfg->insert("cont");
  z_basic_block_t& loop2_head  = cfg->insert("loop2_head");
  z_basic_block_t& loop2_t     = cfg->insert("loop2_t");
  z_basic_block_t& loop2_f     = cfg->insert("loop2_f");
  z_basic_block_t& loop2_body  = cfg->insert("loop2_body");
  z_basic_block_t& ret         = cfg->insert("ret");

  entry >> loop1_head;
  loop1_head >> loop1_t; 
  loop1_head >> loop1_f; 
  loop1_t >>    loop1_body; 

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t; 
  loop2_head >> loop2_f; 
  loop2_t >>    loop2_body; 
  loop2_body >> loop2_head;
  loop2_f >> ret;
  
  z_var i(vfac["i"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  loop1_t.assume(i <= 10);
  loop1_f.assume(i >= 11);
  loop1_body.add(i, i, 1);

  loop1_body_t.assume(i >= 9);
  loop1_body_t.assign(i , 0);
  loop1_body_f.assume(i <= 8);

  loop2_t.assume(i <= 100);
  loop2_f.assume(i >= 101);
  loop2_body.sub(i, i, 1);
  return cfg;
}

z_cfg_t* prog4(variable_factory_t &vfac) {
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  z_basic_block_t& entry      = cfg->insert("entry");
  z_basic_block_t& loop_head  = cfg->insert("loop_head");
  z_basic_block_t& loop_t     = cfg->insert("loop_t");
  z_basic_block_t& loop_f     = cfg->insert("loop_f");
  z_basic_block_t& loop_body  = cfg->insert("loop_body");
  z_basic_block_t& ret        = cfg->insert("ret");

  entry >> loop_head;
  loop_head >> loop_t; 
  loop_head >> loop_f; 
  loop_t >> loop_body; 
  loop_body >> loop_head;
  loop_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  entry.assign(p, 0);

  loop_t.assume(i <= 9);
  loop_f.assume(i >= 10);
  loop_body.add(i, i, 1);
  loop_body.add(p, p, 4);

  return cfg;
}

/* Example of how to build a CFG */
z_cfg_t* prog5(variable_factory_t &vfac)  {
  // Definining program variables
  z_var i(vfac ["i"], crab::INT_TYPE, 32);
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var nd(vfac ["nd"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i != 9);
  bb1_f.assume(i == 9);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

// from TOPLAS paper Fig 1.b: version 1 
// (computing expression y:= 200-2*x in two instructions:
//  tmp:= 2*x; y:= 200-tmp;)
z_cfg_t* prog6(variable_factory_t &vfac)  {
  // Definining program variables
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var n(vfac ["n"], crab::INT_TYPE, 32);
  z_var x(vfac ["x"], crab::INT_TYPE, 32);
  z_var y(vfac ["y"], crab::INT_TYPE, 32);
  z_var t(vfac ["t"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& loop = cfg->insert("loop");
  z_basic_block_t& loop_body_1 = cfg->insert("loop_body_1");
  z_basic_block_t& loop_body_2 = cfg->insert("loop_body_2");
  z_basic_block_t& loop_body_3 = cfg->insert("loop_body_3");
  z_basic_block_t& loop_body_4 = cfg->insert("loop_body_4");  
  z_basic_block_t& ret = cfg->insert("ret");
  // adding control flow
  entry >> loop;
  loop >> loop_body_1;
  loop_body_1 >> loop_body_2;
  loop_body_2 >> loop_body_3;
  loop_body_3 >> loop_body_4;    
  loop_body_4 >> loop;
  loop >> ret;
  // adding statements
  //  entry.assign(x1, 1);
  entry.assign(k, 200);
  entry.assign(n, 100);
  entry.assign(x, 0);
  entry.assign(y, k);
  loop_body_1.assume(x <= n - 1);
  loop_body_2.add(x, x, 1);
  loop_body_3.assign(t, 2*x);
  //loop_body_4.sub(y, k , t);
  loop_body_4.assign(y, 200 - t);
  ret.assume(x >= n);
  ret.assertion(x + y <= k);
  return cfg;
}

// from TOPLAS paper Fig 1.b: version 2
// (computing expression y:= 200-2*x in one instruction)
z_cfg_t* prog7(variable_factory_t &vfac)  {
  // Definining program variables
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var n(vfac ["n"], crab::INT_TYPE, 32);
  z_var x(vfac ["x"], crab::INT_TYPE, 32);
  z_var y(vfac ["y"], crab::INT_TYPE, 32);
  z_var t(vfac ["t"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& loop = cfg->insert("loop");
  z_basic_block_t& loop_body_1 = cfg->insert("loop_body_1");
  z_basic_block_t& loop_body_2 = cfg->insert("loop_body_2");
  z_basic_block_t& loop_body_3 = cfg->insert("loop_body_3");
  z_basic_block_t& loop_body_4 = cfg->insert("loop_body_4");
  z_basic_block_t& ret = cfg->insert("ret");
  // adding control flow
  entry >> loop;
  loop >> loop_body_1;
  loop_body_1 >> loop_body_2;
  loop_body_2 >> loop_body_3;
  loop_body_3 >> loop_body_4;
  loop_body_4 >> loop;
  loop >> ret;
  // adding statements
  entry.assign(k, 200);
  entry.assign(n, 100);
  entry.assign(x, 0);
  entry.assign(y, k);
  loop_body_1.assume(x <= n - 1);
  loop_body_2.add(x, x, 1);
  loop_body_4.assign(y, k -2*x);
  ret.assume(x >= n);
  ret.assertion(x + y <= k);
  return cfg;
}

// from TOPLAS paper Fig 1.b: version 3
// similar to prog7 but simulating SSA
z_cfg_t* prog8(variable_factory_t &vfac)  {
  // Definining program variables
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var n(vfac ["n"], crab::INT_TYPE, 32);
  z_var x(vfac ["x"], crab::INT_TYPE, 32);
  z_var y(vfac ["y"], crab::INT_TYPE, 32);
  z_var x1(vfac ["x1"], crab::INT_TYPE, 32);
  z_var y1(vfac ["y1"], crab::INT_TYPE, 32);
  
  z_var t(vfac ["t"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& loop = cfg->insert("loop");
  z_basic_block_t& loop_body_1 = cfg->insert("loop_body_1");
  z_basic_block_t& loop_body_2 = cfg->insert("loop_body_2");
  z_basic_block_t& loop_body_3 = cfg->insert("loop_body_3");
  z_basic_block_t& loop_body_4 = cfg->insert("loop_body_4");  
  z_basic_block_t& ret = cfg->insert("ret");
  // adding control flow
  entry >> loop;
  loop >> loop_body_1;
  loop_body_1 >> loop_body_2;
  loop_body_2 >> loop_body_3;
  loop_body_3 >> loop_body_4;
  loop_body_4 >> loop;
  loop >> ret;
  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 200);
  loop_body_1.assume(x <= 99);
  #if 0
  loop_body_2.add(x, x, 1);
  loop_body_3.assign(y, 200 -2*x);
  #else
  // loop_body_2.add(x1, x, 1);
  // loop_body_3.assign(y, 200 -2*x1);
  // loop_body_4.assign(x, x1);
  loop_body_2.assign(x1, x);
  loop_body_3.add(x1, x1, 1);  
  loop_body_4.assign(y, 200 -2*x1);
  loop_body_4.assign(x, x1);
  #endif 
  ret.assume(x >= 100);
  ret.assertion(x + y <= 200);
  return cfg;
}

// from TOPLAS paper Fig 1.a: 
// (computing expression y:= 200+2*x in two instructions:
//  tmp:= 2*x; y:= 200+tmp;)
z_cfg_t* prog9(variable_factory_t &vfac)  {
  // Definining program variables
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var n(vfac ["n"], crab::INT_TYPE, 32);
  z_var x(vfac ["x"], crab::INT_TYPE, 32);
  z_var y(vfac ["y"], crab::INT_TYPE, 32);
  z_var t(vfac ["t"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& loop = cfg->insert("loop");
  z_basic_block_t& loop_body_1 = cfg->insert("loop_body_1");
  z_basic_block_t& loop_body_2 = cfg->insert("loop_body_2");
  z_basic_block_t& loop_body_3 = cfg->insert("loop_body_3");
  z_basic_block_t& loop_body_4 = cfg->insert("loop_body_4");  
  z_basic_block_t& ret = cfg->insert("ret");
  // adding control flow
  entry >> loop;
  loop >> loop_body_1;
  loop_body_1 >> loop_body_2;
  loop_body_2 >> loop_body_3;
  loop_body_3 >> loop_body_4;    
  loop_body_4 >> loop;
  loop >> ret;
  // adding statements
  //  entry.assign(x1, 1);
  entry.assign(k, 200);
  entry.assign(n, 100);
  entry.assign(x, 0);
  entry.assign(y, k);
  loop_body_1.assume(x <= n - 1);
  loop_body_2.add(x, x, 1);
  loop_body_3.assign(t, 2*x);
  loop_body_4.add(y, k , t);
  ret.assume(x >= n);
  ret.assertion(y - x >= k);
  return cfg;
}

// from TOPLAS paper Fig 1.b as translated by crab-llvm.
z_cfg_t* prog10(variable_factory_t &vfac)  {
  // Definining program variables
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var n(vfac ["n"], crab::INT_TYPE, 32);
  z_var x(vfac ["x"], crab::INT_TYPE, 32);
  z_var y(vfac ["y"], crab::INT_TYPE, 32);
  z_var x1(vfac ["x1"], crab::INT_TYPE, 32);
  z_var y1(vfac ["y1"], crab::INT_TYPE, 32);
  z_var t(vfac ["t"], crab::INT_TYPE, 32);
  
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& loop = cfg->insert("loop");
  z_basic_block_t& loop_body = cfg->insert("loop_body");
  z_basic_block_t& ret = cfg->insert("ret");
  // adding control flow
  entry >> loop;
  loop >> loop_body;
  loop_body >> loop;
  loop >> ret;
  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 200);
  loop_body.assume(x <= 99);
  loop_body.add(x1, x, 1);
  loop_body.assign(y1, 200 - (2 * x1));
  loop_body.assign(x, x1);
  loop_body.assign(y, y1);
  ret.assume(x >= 100);
  ret.assertion(x + y <= 200);
  return cfg;
}


int main(int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  #if 1
  // From zones test
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog1(vfac);
    // crab::outs() << *cfg << "\n";
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_compare<z_oct_apron_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    // {
    //   z_soct_domain_t inv;
    //   run<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    // }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog2(vfac);
    // crab::outs() << *cfg << "\n";
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_compare<z_oct_apron_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    // {
    //   z_soct_domain_t inv;
    //   run<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    // }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog3(vfac);
    // crab::outs() << *cfg << "\n";
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_compare<z_oct_apron_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif    
    // {
    //   z_soct_domain_t inv;      
    //   run<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    // }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog4(vfac);
    // crab::outs() << *cfg << "\n";
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_compare<z_oct_apron_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif    
    // {
    //   z_soct_domain_t inv;      
    //   run<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    // }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog5(vfac);
    // crab::outs() << *cfg << "\n";
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_compare<z_oct_apron_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg, cfg->entry(), 1, 2, 20, stats_enabled);      
    }
    #endif
    // {
    //   z_soct_domain_t inv;      
    //   run<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    // }
    delete cfg;
  }
  #endif
  
  /*** TOPLAS example: several variants ***/
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog6(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    {
      z_soct_domain_t inv;
      run_and_check<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_check<z_oct_elina_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif    
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  // {
  //   variable_factory_t vfac;
  //   z_cfg_t* cfg = prog7(vfac);
  //   crab::outs() << *cfg << "\n";
  //   // EXPECTED: SAFE
  //   #ifdef HAVE_APRON
  //   run_and_compare<z_soct_domain_t, z_oct_apron_domain_t>
  //     (cfg, cfg->entry(), 1, 2, 20, stats_enabled);
  //   #endif 
  //   delete cfg;
  // }
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog7(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    {
      z_soct_domain_t inv;
      run_and_check<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_check<z_oct_elina_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif     
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog8(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    {
      z_soct_domain_t inv;
      run_and_check<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_check<z_oct_elina_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif         
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog9(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    {
      z_soct_domain_t inv;
      run_and_check<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_check<z_oct_elina_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif        
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  // FIXME: soct is imprecise if DefaultParams.special_assign = 1   
  { 
    variable_factory_t vfac;
    z_cfg_t* cfg = prog10(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    {
      z_soct_domain_t inv;
      run_and_check<z_soct_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    } 
    #ifdef HAVE_APRON
    {
      z_oct_apron_domain_t inv;
      run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif
    #ifdef HAVE_ELINA
    {
      z_oct_elina_domain_t inv;
      run_and_check<z_oct_elina_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    }
    #endif        
    delete cfg;
  }
  
  return 0;
}
