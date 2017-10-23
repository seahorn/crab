#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog1 (variable_factory_t &vfac)  {

  // Defining program variables
  z_var i (vfac ["i"]);
  z_var nd (vfac ["nd"]);
  z_var inc (vfac ["inc"]);
  auto b1 = vfac ["b1"];
  auto b2 = vfac ["b2"];
  auto b3 = vfac ["b3"];
  auto b4 = vfac ["b4"];    
  auto bfalse = vfac ["bf"];
  auto btrue = vfac ["bt"];      
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, 0);
  entry.bool_assign (bfalse, z_lin_cst_t::get_false ());
  entry.bool_assign (btrue, z_lin_cst_t::get_true ());  
  bb1_t.bool_assign (b1, i <= 99);  
  bb1_t.bool_assume (b1);
  bb1_t.bool_assert (b1);        // trivial
  bb1_f.bool_assign (b2, i >= 100);  
  bb1_f.bool_assume (b2);
  bb2.havoc(nd.name());
  bb2.select(inc,nd,1,2);
  bb2.add(i, i, inc);
  bb2.bool_assign (b3, i >= 1);
  bb2.bool_assign (b4, b3);  
  bb2.bool_or (b4, b4, bfalse);  // tautology
  bb2.bool_and (b4, b4, btrue);  // tautology
  bb2.bool_xor (b4, b4, btrue);  // complement
  bb2.bool_xor (b4, b4, btrue);  // complement  

  return cfg;
}


z_cfg_t* prog2 (variable_factory_t &vfac)  {

  // Defining program variables
  z_var i (vfac ["i"]);
  auto b = vfac ["b"];
  z_var n (vfac ["n"]);
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, 0);
  entry.havoc(n.name()); 
  entry.bool_assign(b,  n == z_number(10));
  //entry.assign(n, z_number(1));
  entry.bool_assume(b);
  bb1_t.assume (i <= n);
  bb2.add(i,i,1);
  bb1_f.assume (i >= n + 1);
  ret.assertion(i == 10);

  return cfg;
}

z_cfg_t* prog3 (variable_factory_t &vfac)  {

  // Defining program variables
  z_var i (vfac ["i"]);
  auto b = vfac ["b"];
  z_var n (vfac ["n"]);
  
  // entry and exit block
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb0 = cfg->insert ("bb0");
  z_basic_block_t& entry_cnt = cfg->insert ("entry_cnt");    
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb0;
  bb0 >> entry_cnt;
  entry >>  entry_cnt;
  entry_cnt >>	bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, 0);
  entry.havoc(n.name()); 
  entry.bool_assign(b,  n == z_number(10));
  bb0.assign(n, z_number(1));
  entry_cnt.bool_assume(b);
  bb1_t.assume (i <= n);
  bb2.add(i,i,1);
  bb1_f.assume (i >= n + 1);
  ret.assertion(i == 10);

  return cfg;
}


/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  
  {
    z_cfg_t* cfg = prog1(vfac);
    crab::outs() << *cfg << "\n";
    
    run<z_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
    run<z_bool_num_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
    #ifdef HAVE_LDD
    run<z_boxes_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
    #endif
    delete cfg;    
  }
  { // precise
    z_cfg_t* cfg = prog2(vfac);
    crab::outs() << *cfg << "\n";
    run<z_bool_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
    delete cfg;    
  }
  {
    // imprecise
    z_cfg_t* cfg = prog3(vfac);
    crab::outs() << *cfg << "\n";
    run<z_bool_interval_domain_t>(cfg,  false, 1, 2, 20, stats_enabled);
    delete cfg;    
  }
  
  return 0;
}
