#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

// /* Example of how to build a CFG */
// cfg_t prog (VariableFactory &vfac)  {

//   // Definining program variables
//   z_var i (vfac ["i"]);
//   z_var k (vfac ["k"]);
//   z_var nd (vfac ["nd"]);
//   z_var inc (vfac ["inc"]);
//   // entry and exit block
//   cfg_t cfg ("x0","ret");
//   // adding blocks
//   basic_block_t& x0 = cfg.insert ("x0");
//   basic_block_t& x1 = cfg.insert ("x1");
//   basic_block_t& x2 = cfg.insert ("x2");
//   basic_block_t& x3 = cfg.insert ("x3");
//   basic_block_t& entry = cfg.insert ("entry");
//   basic_block_t& bb1   = cfg.insert ("bb1");
//   basic_block_t& bb1_t = cfg.insert ("bb1_t");
//   basic_block_t& bb1_f = cfg.insert ("bb1_f");
//   basic_block_t& bb2   = cfg.insert ("bb2");
//   basic_block_t& ret   = cfg.insert ("ret");
//   // adding control flow
//   x0 >> x1; x1 >> x2; x2 >> x3; x3 >> entry;
//   entry >> bb1;
//   bb1 >> bb1_t; bb1 >> bb1_f;
//   bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
//   // adding statements
//   x0.assign (k, 50);
//   entry.assign (i, 0);
//   bb1_t.assume (i <= 99);
//   bb1_f.assume (i >= 100);
//   bb2.havoc(nd.name());
//   bb2.select(inc,nd,1,2);
//   bb2.add(i, i, inc);

//   return cfg;
// }

/* Example of how to build a CFG */
cfg_t prog (VariableFactory &vfac)  {

  // Definining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var nd (vfac ["nd"]);
  z_var tmp (vfac ["tmp"]);
  // entry and exit block
  cfg_t cfg ("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& bb1   = cfg.insert ("bb1");
  basic_block_t& bb1_t = cfg.insert ("bb1_t");
  basic_block_t& bb1_f = cfg.insert ("bb1_f");
  basic_block_t& bb2   = cfg.insert ("bb2");
  basic_block_t& bb2_t = cfg.insert ("bb2_t");
  basic_block_t& bb2_f = cfg.insert ("bb2_f");

  basic_block_t& ret   = cfg.insert ("ret");
  // adding control flow

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb1_f >> bb2;
  bb2 >> bb2_t;   bb2 >> bb2_f; 
  bb2_t >> ret; bb2_f >> ret;
  
  // adding statements
  entry.havoc (nd.name());
  bb1_t.assume (nd >= 1);
  bb1_t.sub (x, 0, 10);
  bb1_f.assume (nd <= 0);
  bb1_f.assign (x, 20);
  bb2.assign (y,x);
  bb2_t.assume (y <= 0);
  bb2_t.mul(tmp, x, -1);
  bb2_t.assign(y, tmp);
  bb2_f.assume (y >= 1);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  VariableFactory vfac;
  cfg_t cfg = prog (vfac);
  cfg.simplify (); // this is optional
  cout << cfg << endl;
  const bool run_live = false;

  {
    NumFwdAnalyzer <cfg_t, dbm_domain_t,VariableFactory>::type a (cfg,vfac,run_live);
    // Run fixpoint 
    dbm_domain_t inv = dbm_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_t, interval_domain_t,VariableFactory>::type a (cfg,vfac,run_live);
    // Run fixpoint 
    interval_domain_t inv = interval_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_t,sym_cst_intv_domain_t,VariableFactory>::type a (cfg,vfac,run_live);
    // Run fixpoint 
    sym_cst_intv_domain_t inv = sym_cst_intv_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }


#if 0
  {
    NumFwdAnalyzer <cfg_t, dbm_domain_t,VariableFactory>::type a (cfg,vfac,run_live);
    // Run fixpoint 
    dbm_domain_t inv = dbm_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }
  {
    NumFwdAnalyzer <cfg_t, term_domain_t,VariableFactory>::type a (cfg,vfac,run_live);
    // Run fixpoint 
    term_domain_t inv = term_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }
#endif 

  return 0;
}
