#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t cfg1 (VariableFactory &vfac) 
{

  ////
  // Building the CFG
  ////

  // entry and exit block
  cfg_t cfg ("b0","b3",PTR);
  // adding blocks
  basic_block_t& b0 = cfg.insert ("b0");
  basic_block_t& b1 = cfg.insert ("b1");
  basic_block_t& b2 = cfg.insert ("b2");
  basic_block_t& b3 = cfg.insert ("b3");
  // adding control flow
  b0 >> b1; b0 >> b2; b1 >> b3; b2 >> b3;

  // definining program variables
  varname_t p = vfac ["p"];
  varname_t q = vfac ["q"];
  varname_t r = vfac ["r"];
  varname_t s = vfac ["s"];
  // adding statements
  b0.new_object (p, 1);
  b0.havoc (r);
  b0.havoc (s);
  b1.ptr_assume (ptr_cst_t::mk_eq_null (r));
  b1.ptr_assume (ptr_cst_t::mk_eq (r,s));
  b2.ptr_assume (ptr_cst_t::mk_diseq (r,s));
  b2.new_object (q, 2);
  b2.ptr_assign (p, q, z_number(4));
  return cfg;
}

cfg_t cfg2 (VariableFactory &vfac) 
{

  // entry and exit block
  cfg_t cfg ("b0","b3",PTR);
  // adding blocks
  basic_block_t& b0 = cfg.insert ("b0");
  basic_block_t& b1 = cfg.insert ("b1");
  basic_block_t& b2 = cfg.insert ("b2");
  basic_block_t& b3 = cfg.insert ("b3");
  // adding control flow
  b0 >> b1; b0 >> b2; b1 >> b3; b2 >> b3;
  

  // definining program variables
  varname_t p = vfac ["p"];
  varname_t q1 = vfac ["q1"];
  varname_t q2 = vfac ["q2"];
  varname_t r = vfac ["r"];
  z_var nd (vfac ["nd"]);
  // adding statements
  b0.new_object (p , 1);  // p = malloc (...);
  b0.new_object (q1, 2);  // q1 = malloc (...);
  b0.new_object (q2, 3);  // q2 = malloc (...);
  b0.havoc (nd.name ());
  b1.assume (nd >= 1);
  b2.assume (nd <= 0);
  b1.ptr_store (p, q1, z_interval (0,3));  // *p = q1
  b2.ptr_store (p, q2, z_interval (0,3));  // *p = q2
  b3.ptr_load (r, p, z_interval (0,3));    // r = *p
  return cfg;
}



cfg_t cfg3 (VariableFactory &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  varname_t p (vfac ["p"]);
  varname_t q (vfac ["q"]);
  // entry and exit block
  cfg_t cfg ("entry","ret",PTR);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& bb1   = cfg.insert ("bb1");
  basic_block_t& bb1_t = cfg.insert ("bb1_t");
  basic_block_t& bb1_f = cfg.insert ("bb1_f");
  basic_block_t& bb2   = cfg.insert ("bb2");
  basic_block_t& ret   = cfg.insert ("ret");
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
  bb2.new_object (q, 1);
  bb2.ptr_assign (p, q, z_number(4));
  bb2.add(i, i, 1);
  ret.assume (x <= y);

  return cfg;
}


void run (cfg_t &cfg, VariableFactory& vfac) {
  typedef NullityAnalyzer<cfg_t, VariableFactory>::analyzer_t nullity_analyzer_t;
  typedef NullityAnalyzer<cfg_t, VariableFactory>::nullity_domain_t nullity_domain_t;

  Liveness<cfg_t> live (cfg);
  //live.exec ();

  nullity_analyzer_t a (cfg, vfac, &live);
  cout << cfg << "\n";

  a.Run (nullity_domain_t::top ());

  std::cout << "Nullity analysis \n";
  for (auto &b : cfg)  {
    auto pre = a.get_pre (b.label ());
    auto post = a.get_post (b.label ());
    std::cout << get_label_str (b.label ()) << "=" 
              << pre 
              << " ==> "
              << post << "\n";
  }

}

int main () {
  VariableFactory vfac;
  cfg_t cfg_1 = cfg1 (vfac);
  run (cfg_1, vfac);
  cfg_t cfg_2 = cfg2 (vfac);
  run (cfg_2, vfac);
  cfg_t cfg_3 = cfg3 (vfac);
  run (cfg_3, vfac);
  
  return 0;
}
