#include "../common.hpp"

#include <crab/cg/Cg.hpp>
#include <crab/analysis/InterFwdAnalyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

cfg_t foo (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["x"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["foo"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add (y, x, 1);
  exit.add (z , y , 2);
  exit.ret (vfac ["z"], INT_TYPE);
  return cfg;
}

cfg_t rec1 (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["s"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["rec1"], params);
  // Defining program variables
  z_var r (vfac ["r"]);
  z_var s (vfac ["s"]);
  z_var t (vfac ["t"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.sub (r, s, 1);
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["r"], INT_TYPE));
  exit.callsite (make_pair (vfac["t"], INT_TYPE), vfac ["rec2"], args);
  exit.ret (vfac["t"], INT_TYPE);
  return cfg;
}

cfg_t rec2 (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["s1"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["rec2"], params);
  // Defining program variables
  z_var r (vfac ["r1"]);
  z_var s (vfac ["s1"]);
  z_var t (vfac ["t1"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.sub (r, s, 1);
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["r1"], INT_TYPE));
  exit.callsite (make_pair (vfac["t1"], INT_TYPE), vfac ["rec1"], args);
  args.clear ();
  //args.push_back (make_pair (vfac["t1"], INT_TYPE));
  //exit.callsite (make_pair (vfac["t1"], INT_TYPE), vfac ["foo"], args);
  exit.ret (vfac["t1"], INT_TYPE);
  return cfg;
}


cfg_t bar (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["a"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["bar"], params);
  // Defining program variables
  z_var a (vfac ["a"]);
  z_var x (vfac ["x1"]);
  z_var y (vfac ["y1"]);
  z_var w (vfac ["w1"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["x1"], INT_TYPE));
  exit.callsite (make_pair (vfac["y1"], INT_TYPE), vfac ["foo"], args);
  entry.assign (x, a);
  entry.assign (w, 5);
  exit.ret (vfac["y1"], INT_TYPE);
  return cfg;
}

cfg_t m (VariableFactory &vfac)  {
  vector<pair<varname_t,VariableType> > params;
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["main"], params);
  // Defining program variables
  z_var x (vfac ["x2"]);
  z_var y (vfac ["y2"]);
  z_var z (vfac ["z2"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["x2"], INT_TYPE));
  entry.callsite (make_pair (vfac["y2"], INT_TYPE), vfac ["bar"], args);
  args.clear ();
  /////
  args.push_back (make_pair (vfac["y2"], INT_TYPE));
  entry.callsite (make_pair (vfac["z3"], INT_TYPE), vfac ["rec1"], args);
  args.clear ();
  /////
  exit.add (z, y, 2);
  args.push_back (make_pair (vfac["z2"], INT_TYPE));
  exit.callsite (make_pair (vfac["w2"], INT_TYPE), vfac ["foo"], args);
  exit.ret (vfac["w2"], INT_TYPE);
  return cfg;
}

int main (int argc, char** argv ) {

  VariableFactory vfac;
  cfg_t t1 = foo (vfac);
  cfg_t t2 = bar (vfac);
  cfg_t t3 = rec1 (vfac);
  cfg_t t4 = rec2 (vfac);
  cfg_t t5 = m (vfac);

  cout << t1 << endl;
  cout << t2 << endl;
  cout << t3 << endl;
  cout << t4 << endl;
  cout << t5 << endl;

  vector<cfg_t> cfgs;
  cfgs.push_back(t1);
  cfgs.push_back(t2);
  cfgs.push_back(t3);
  cfgs.push_back(t4);
  cfgs.push_back(t5);

  CallGraph<cfg_t> cg (cfgs);

  const bool run_live = false;

  InterFwdAnalyzer<CallGraph<cfg_t>, VariableFactory,
                   dbm_domain_t, interval_domain_t>
      a (cg, vfac, run_live);

  a.Run ();
  
  // Print invariants
  for (auto &cfg : cfgs) {
    auto fdecl_opt = cfg.get_func_decl ();
    assert (fdecl_opt);
    cout << *fdecl_opt << "\n"; 
    for (auto &b : cfg) {
      auto inv = a.get_post (cfg, b.label ());
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
    cout << "=================================\n";
  }

  // Print summaries
  for (auto &cfg : cfgs) {
    if (a.has_summary (cfg)) {
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      cout << "Summary for " << *fdecl_opt << ": "; 
      auto sum = a.get_summary (cfg);
      cout << sum << "\n";
    }
  }

  return 0;
}
