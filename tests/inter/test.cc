#include "../common.hpp"

#include <crab/cg/cg_bgl.hpp>
#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter_fwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

cfg_t* foo (variable_factory_t &vfac) {
  vector<pair<varname_t,variable_type> > params;
  params.push_back (make_pair (vfac["x"], INT_TYPE));
  function_decl<varname_t> decl (INT_TYPE, vfac["foo"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add (y, x, 1);
  exit.add (z , y , 2);
  exit.ret (vfac ["z"], INT_TYPE);
  return cfg;
}

cfg_t* rec1 (variable_factory_t &vfac) {
  vector<pair<varname_t,variable_type> > params;
  params.push_back (make_pair (vfac["s"], INT_TYPE));
  function_decl<varname_t> decl (INT_TYPE, vfac["rec1"], params);
  // Defining program variables
  z_var r (vfac ["r"]);
  z_var s (vfac ["s"]);
  z_var t (vfac ["t"]);
  // entry and exit block
  cfg_t* cfg  = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.sub (r, s, 1);
  vector<pair<varname_t,variable_type> > args;
  args.push_back (make_pair (vfac["r"], INT_TYPE));
  exit.callsite (make_pair (vfac["t"], INT_TYPE), vfac ["rec2"], args);
  exit.ret (vfac["t"], INT_TYPE);
  return cfg;
}

cfg_t* rec2 (variable_factory_t &vfac) {
  vector<pair<varname_t,variable_type> > params;
  params.push_back (make_pair (vfac["s1"], INT_TYPE));
  function_decl<varname_t> decl (INT_TYPE, vfac["rec2"], params);
  // Defining program variables
  z_var r (vfac ["r1"]);
  z_var s (vfac ["s1"]);
  z_var t (vfac ["t1"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.sub (r, s, 1);
  vector<pair<varname_t,variable_type> > args;
  args.push_back (make_pair (vfac["r1"], INT_TYPE));
  exit.callsite (make_pair (vfac["t1"], INT_TYPE), vfac ["rec1"], args);
  args.clear ();
  //args.push_back (make_pair (vfac["t1"], INT_TYPE));
  //exit.callsite (make_pair (vfac["t1"], INT_TYPE), vfac ["foo"], args);
  exit.ret (vfac["t1"], INT_TYPE);
  return cfg;
}


cfg_t* bar (variable_factory_t &vfac) {
  vector<pair<varname_t,variable_type> > params;
  params.push_back (make_pair (vfac["a"], INT_TYPE));
  function_decl<varname_t> decl (INT_TYPE, vfac["bar"], params);
  // Defining program variables
  z_var a (vfac ["a"]);
  z_var x (vfac ["x1"]);
  z_var y (vfac ["y1"]);
  z_var w (vfac ["w1"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  vector<pair<varname_t,variable_type> > args;
  args.push_back (make_pair (vfac["x1"], INT_TYPE));
  exit.callsite (make_pair (vfac["y1"], INT_TYPE), vfac ["foo"], args);
  entry.assign (x, a);
  entry.assign (w, 5);
  exit.ret (vfac["y1"], INT_TYPE);
  return cfg;
}

cfg_t* m (variable_factory_t &vfac)  {
  vector<pair<varname_t,variable_type> > params;
  function_decl<varname_t> decl (INT_TYPE, vfac["main"], params);
  // Defining program variables
  z_var x (vfac ["x2"]);
  z_var y (vfac ["y2"]);
  z_var z (vfac ["z2"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  vector<pair<varname_t,variable_type> > args;
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

  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  cfg_t* t1 = foo (vfac);
  cfg_t* t2 = bar (vfac);
  cfg_t* t3 = rec1 (vfac);
  cfg_t* t4 = rec2 (vfac);
  cfg_t* t5 = m (vfac);

  crab::outs() << *t1 << "\n";
  crab::outs() << *t2 << "\n";
  crab::outs() << *t3 << "\n";
  crab::outs() << *t4 << "\n";
  crab::outs() << *t5 << "\n";

  vector<cfg_ref_t> cfgs;
  cfgs.push_back(*t1);
  cfgs.push_back(*t2);
  cfgs.push_back(*t3);
  cfgs.push_back(*t4);
  cfgs.push_back(*t5);

  typedef call_graph<cfg_ref_t> callgraph_t;
  typedef call_graph_ref<callgraph_t> callgraph_ref_t;

  boost::scoped_ptr<callgraph_t> cg(new callgraph_t(cfgs));
  {
    inter_fwd_analyzer<callgraph_ref_t, variable_factory_t,
                     dbm_domain_t, interval_domain_t> a (*cg, vfac, nullptr); 
    crab::outs() << "Running" 
         << " summary domain=" << dbm_domain_t::getDomainName () 
         << " and forward domain=" << interval_domain_t::getDomainName () << "\n";

    a.Run ();
    
    if (stats_enabled) {
      crab::CrabStats::Print(crab::outs());
      crab::CrabStats::reset();
    }  

    // Print invariants
    for (auto cfg : cfgs) {
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      crab::outs() << *fdecl_opt << "\n"; 
      for (auto &b : cfg) {
        auto inv = a.get_post (cfg, b.label ());
        crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
      }
      crab::outs() << "=================================\n";
    }
    
    // Print summaries
    for (auto cfg : cfgs) {
      if (a.has_summary (cfg)) {
        auto fdecl_opt = cfg.get_func_decl ();
        assert (fdecl_opt);
        crab::outs() << "Summary for " << *fdecl_opt << ": "; 
        auto sum = a.get_summary (cfg);
        crab::outs() << sum << "\n";
      }
    }
  }


#ifdef HAVE_APRON
  {
    inter_fwd_analyzer<callgraph_ref_t, variable_factory_t,
                     opt_oct_apron_domain_t, interval_domain_t> a (*cg, vfac, nullptr); 

    crab::outs() << "Running" 
         << " summary domain=" << opt_oct_apron_domain_t::getDomainName () 
         << " and forward domain=" << interval_domain_t::getDomainName () << "\n";

    a.Run ();
    if (stats_enabled) {
      crab::CrabStats::Print(crab::outs());
      crab::CrabStats::reset();
    }  

    
    // Print invariants
    for (auto cfg : cfgs) {
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      crab::outs() << *fdecl_opt << "\n"; 
      for (auto &b : cfg) {
        auto inv = a.get_post (cfg, b.label ());
        crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
      }
      crab::outs() << "=================================\n";
    }
    
    // Print summaries
    for (auto cfg : cfgs) {
      if (a.has_summary (cfg)) {
        auto fdecl_opt = cfg.get_func_decl ();
        assert (fdecl_opt);
        crab::outs() << "Summary for " << *fdecl_opt << ": "; 
        auto sum = a.get_summary (cfg);
        crab::outs() << sum << "\n";
      }
    }
  }
#endif 

  {
    inter_fwd_analyzer<callgraph_ref_t,  variable_factory_t,
                     term_domain_t, interval_domain_t> a (*cg, vfac, nullptr); 

    crab::outs() << "Running" 
         << " summary domain=" << term_domain_t::getDomainName () 
         << " and forward domain=" << interval_domain_t::getDomainName () << "\n";

    a.Run ();
    if (stats_enabled) {
      crab::CrabStats::Print(crab::outs());
      crab::CrabStats::reset();
    }  

    
    // Print invariants
    for (auto cfg : cfgs) {
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      crab::outs() << *fdecl_opt << "\n"; 
      for (auto &b : cfg) {
        auto inv = a.get_post (cfg, b.label ());
        crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
      }
      crab::outs() << "=================================\n";
    }
    
    // Print summaries
    for (auto cfg : cfgs) {
      if (a.has_summary (cfg)) {
        auto fdecl_opt = cfg.get_func_decl ();
        assert (fdecl_opt);
        crab::outs() << "Summary for " << *fdecl_opt << ": "; 
        auto sum = a.get_summary (cfg);
        crab::outs() << sum << "\n";
      }
    }
  }

  {
    inter_fwd_analyzer<callgraph_ref_t,  variable_factory_t,
                       num_domain_t, num_domain_t> a (*cg, vfac, nullptr); 

    crab::outs() << "Running" 
         << " summary domain=" << num_domain_t::getDomainName () 
         << " and forward domain=" << num_domain_t::getDomainName () << "\n";

    a.Run ();

    if (stats_enabled) {
      crab::CrabStats::Print(crab::outs());
      crab::CrabStats::reset();
    }  
    
    // Print invariants
    for (auto cfg : cfgs) {
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      crab::outs() << *fdecl_opt << "\n"; 
      for (auto &b : cfg) {
        auto inv = a.get_post (cfg, b.label ());
        crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
      }
      crab::outs() << "=================================\n";
    }
    
    // Print summaries
    for (auto cfg : cfgs) {
      if (a.has_summary (cfg)) {
        auto fdecl_opt = cfg.get_func_decl ();
        assert (fdecl_opt);
        crab::outs() << "Summary for " << *fdecl_opt << ": "; 
        auto sum = a.get_summary (cfg);
        crab::outs() << sum << "\n";
      }
    }
  }

  delete t1;
  delete t2;
  delete t3;
  delete t4;
  delete t5;

  return 0;
}
