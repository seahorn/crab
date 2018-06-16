#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* cfg1 (variable_factory_t &vfac) 
{

  ////
  // Building the CFG
  ////

  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("b0","b3",PTR);
  // adding blocks
  z_basic_block_t& b0 = cfg->insert ("b0");
  z_basic_block_t& b1 = cfg->insert ("b1");
  z_basic_block_t& b2 = cfg->insert ("b2");
  z_basic_block_t& b3 = cfg->insert ("b3");
  // adding control flow
  b0 >> b1; b0 >> b2; b1 >> b3; b2 >> b3;

  // definining program variables
  z_var p(vfac ["p"], crab::PTR_TYPE);
  z_var q(vfac ["q"], crab::PTR_TYPE);
  z_var r(vfac ["r"], crab::PTR_TYPE);
  z_var s(vfac ["s"], crab::PTR_TYPE);
  // adding statements
  b0.ptr_new_object (p, 1);
  b0.havoc (r);
  b0.havoc (s);
  b1.ptr_assume (z_ptr_cst_t::mk_eq_null (r));
  b1.ptr_assume (z_ptr_cst_t::mk_eq (r,s));
  b2.ptr_assume (z_ptr_cst_t::mk_diseq (r,s));
  b2.ptr_new_object (q, 2);
  b2.ptr_assign (p, q, z_number(4));
  return cfg;
}

z_cfg_t* cfg2 (variable_factory_t &vfac) 
{

  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("b0","b3",PTR);
  // adding blocks
  z_basic_block_t& b0 = cfg->insert ("b0");
  z_basic_block_t& b1 = cfg->insert ("b1");
  z_basic_block_t& b2 = cfg->insert ("b2");
  z_basic_block_t& b3 = cfg->insert ("b3");
  // adding control flow
  b0 >> b1; b0 >> b2; b1 >> b3; b2 >> b3;
  

  // definining program variables
  z_var p(vfac ["p"], crab::PTR_TYPE);
  z_var q1(vfac ["q1"], crab::PTR_TYPE);
  z_var q2(vfac ["q2"], crab::PTR_TYPE);
  z_var r(vfac ["r"], crab::PTR_TYPE);
  z_var nd(vfac ["nd"], crab::INT_TYPE, 32);
  // adding statements
  b0.ptr_new_object (p , 1);  // p = malloc (...);
  b0.ptr_new_object (q1, 2);  // q1 = malloc (...);
  b0.ptr_new_object (q2, 3);  // q2 = malloc (...);
  b0.havoc (nd);
  b1.assume (nd >= 1);
  b2.assume (nd <= 0);
  b1.ptr_store (p, q1);  // *p = q1
  b2.ptr_store (p, q2);  // *p = q2
  b3.ptr_load (r, p);    // r = *p
  return cfg;
}


z_cfg_t* cfg3 (variable_factory_t &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"], crab::INT_TYPE, 32);
  z_var x (vfac ["x"], crab::INT_TYPE, 32);
  z_var y (vfac ["y"], crab::INT_TYPE, 32);
  z_var p (vfac ["p"], crab::PTR_TYPE);
  z_var q (vfac ["q"], crab::PTR_TYPE);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret",PTR);
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
  entry.assign (x, 1);
  entry.assign (y, 0);
  entry.ptr_null (p);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.add(x,x,y);
  bb2.add(y,y,1);
  bb2.ptr_new_object (q, 1);
  bb2.ptr_assign (p, q, z_number(4));
  bb2.add(i, i, 1);
  ret.assume (x <= y);

  return cfg;
}


int main (int argc, char** argv) {
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg_1 = cfg1 (vfac);
  crab::outs () << *cfg_1 << "\n";
  run<z_nullity_domain_t>(cfg_1,cfg_1->entry(),false,1,2,20,stats_enabled);
  z_cfg_t* cfg_2 = cfg2 (vfac);
  crab::outs () << *cfg_2 << "\n";
  run<z_nullity_domain_t>(cfg_2,cfg_2->entry(),false,1,2,20,stats_enabled);
  z_cfg_t* cfg_3 = cfg3 (vfac);
  crab::outs () << *cfg_3 << "\n";
  run<z_nullity_domain_t>(cfg_3,cfg_3->entry(),false,1,2,20,stats_enabled);
  
  delete cfg_1;
  delete cfg_2;
  delete cfg_3;

  return 0;
}
