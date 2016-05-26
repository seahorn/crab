#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
cfg_t* prog (VariableFactory &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  z_var w (vfac ["w"]);
  z_var s (vfac ["s"]);
  z_var t (vfac ["t"]);
  z_var nd1 (vfac ["nd1"]);
  z_var nd2 (vfac ["nd2"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow 
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (i, 0);
  entry.assign (x, 5);
  entry.assign (y, 5);
  entry.assign (z, 3);
  entry.assign (w, 3);
  entry.assign (s, 0);
  entry.assign (t, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.havoc(nd1.name());
  bb2.havoc(nd2.name());
  bb2.bitwise_and(x,x,nd1);
  bb2.bitwise_and(y,y,nd1);
  bb2.bitwise_or(z,z,nd1);
  bb2.bitwise_or(w,w,nd1);
  bb2.bitwise_xor(s,nd1,nd2);
  bb2.bitwise_xor(t,nd1,nd2);
  bb2.add(i, i, 1);

  return cfg;
}

template <typename Domain, typename Live>
void run(cfg_ref_t cfg, VariableFactory &vfac, Live live)
{
  const unsigned int w = 1;
  const unsigned int n = 2;

  typename NumFwdAnalyzer <cfg_ref_t, Domain, VariableFactory>::type 
      It (cfg, vfac, live, w, n, 20);
  Domain inv = Domain::top ();
  It.Run (inv);
  crab::outs() << "Invariants using " << Domain::getDomainName () << ":\n";

  for (auto &b : cfg)
  {
    // invariants at the entry of the block
    auto inv = It [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  crab::outs() << "\n";
  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  VariableFactory vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify (); // this is optional
  crab::outs() << *cfg << "\n";

  run<interval_domain_t>(*cfg,vfac,nullptr);
  run<dbm_domain_t>(*cfg,vfac,nullptr);
  run<term_domain_t>(*cfg,vfac,nullptr);

  return 0;
}
