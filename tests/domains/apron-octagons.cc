#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/fwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

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
#ifdef HAVE_APRON  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  /*** TOPLAS example: several variants ***/
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog6(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    z_oct_apron_domain_t inv;
    run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog7(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    z_oct_apron_domain_t inv;
    run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog8(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    z_oct_apron_domain_t inv;
    run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog9(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    z_oct_apron_domain_t inv;
    run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
  // crab::outs() << "##============================================##\n";
  { 
    variable_factory_t vfac;
    z_cfg_t* cfg = prog10(vfac);
    crab::outs() << *cfg << "\n";
    // EXPECTED: SAFE
    z_oct_apron_domain_t inv;
    run_and_check<z_oct_apron_domain_t>(cfg, cfg->entry(), inv, false, 1, 2, 20, stats_enabled);
    delete cfg;
  }
#endif
  return 0;
}
