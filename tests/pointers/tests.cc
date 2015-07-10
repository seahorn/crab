#include <ikos/tests/Cfg_impl.hpp>
#include <ikos/cfg/VarFactory.hpp>

#include <ikos/common/types.hpp>
#include <ikos/analysis/Pointer.hpp>

using namespace std;

using namespace cfg_impl;
using namespace analyzer;

typedef interval <z_number> z_interval;

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
  z_var nd (vfac ["nd"]);
  // adding statements
  b0.new_object (p, 1);
  b0.havoc (nd.name ());
  b1.assume (nd >= 1);
  b2.assume (nd <= 0);
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


cfg_t cfg3 (VariableFactory &vfac) 
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
  b2.ptr_store (p, q2, z_interval (4,7));  // *p = q2
  b3.ptr_load (r, p, z_interval (0,3));    // r = *p
  return cfg;
}

cfg_t foo (VariableFactory &vfac) 
{
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["p"], PTR_TYPE));
  params.push_back (make_pair (vfac["q"], PTR_TYPE));
  FunctionDecl<varname_t> decl (PTR_TYPE, vfac["foo"], params);
  // entry and exit block
  cfg_t cfg ("b0","b6", decl, PTR);
  // adding blocks
  basic_block_t& b0 = cfg.insert ("b0");
  basic_block_t& b1 = cfg.insert ("b1");
  basic_block_t& b2 = cfg.insert ("b2");
  basic_block_t& b3 = cfg.insert ("b3");
  basic_block_t& b4 = cfg.insert ("b4");
  basic_block_t& b5 = cfg.insert ("b5");
  basic_block_t& b6 = cfg.insert ("b6");
  // adding control flow
  b0 >> b1; b0 >> b2; b1 >> b3; b2 >> b3; b3>>b4; b3>>b5; b4>>b6; b5>>b6;
  
  // definining program variables
  varname_t p = vfac ["p"];
  varname_t q = vfac ["q"];
  varname_t r = vfac ["r"];
  z_var nd (vfac ["nd"]);
  // adding statements
  b0.havoc (nd.name ());
  b1.assume (nd >= 1);
  b2.assume (nd <= 0);
  b1.new_object (p , 1);  // p = &(1);
  b2.new_object (p , 2);  // p = &(2);
  b3.havoc (nd.name ());
  b4.assume (nd >= 1);
  b5.assume (nd <= 0);
  b4.new_object (q , 3);  // p = &(3);
  b5.new_object (q , 4);  // p = &(4);

  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (p, PTR_TYPE));
  args.push_back (make_pair (q, PTR_TYPE));
  b6.callsite (make_pair (r, PTR_TYPE), vfac ["bar"], args);
  b6.ret (r, PTR_TYPE);
  return cfg;
}


cfg_t bar (VariableFactory &vfac) 
{
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["x1"], PTR_TYPE));
  params.push_back (make_pair (vfac["x2"], PTR_TYPE));
  FunctionDecl<varname_t> decl (PTR_TYPE, vfac["bar"], params);

  // entry and exit block
  cfg_t cfg ("b0","b1",decl,PTR);
  // adding blocks
  basic_block_t& b0 = cfg.insert ("b0");
  basic_block_t& b1 = cfg.insert ("b1");
  // adding control flow
  b0 >> b1; 
  
  // definining program variables
  varname_t x1 = vfac ["x1"];
  varname_t x2 = vfac ["x2"];
  varname_t x3 = vfac ["x3"];
  // adding statements
  b0.ptr_assign (x3, x1, z_number(0));
  b0.ptr_assign (x3, x2, z_number(0));
  b1.ret (x3, PTR_TYPE);
  return cfg;
}

void run (vector<cfg_t> cfgs, VariableFactory& vfac)
{

  Pointer<cfg_t, VariableFactory> pta (vfac);
  for (auto cfg: cfgs)
  {
    cfg.simplify ();
    cout << cfg << endl;
    pta.gen_constraints  (cfg);
  }

  pta.solve ();
  cout << "Pointer information\n";
  cout << pta << endl; 
}


int main ()
{
  {
    // intra-procedural example
    VariableFactory vfac;
    vector<cfg_t> cfgs;
    cfgs.push_back (cfg1 (vfac));
    run (cfgs, vfac);
  }

  {
    // intra-procedural example
    VariableFactory vfac;
    vector<cfg_t> cfgs;
    cfgs.push_back (cfg2 (vfac));
    run (cfgs, vfac);
  }

  {
    // intra-procedural example
    VariableFactory vfac;
    vector<cfg_t> cfgs;
    cfgs.push_back (cfg3 (vfac));
    run (cfgs, vfac);
  }

  {
    // inter-procedural example
    VariableFactory vfac;
    vector<cfg_t> cfgs;
    cfgs.push_back (foo (vfac));
    cfgs.push_back (bar (vfac));
    run (cfgs, vfac);
  }
  
  return 0;
}
