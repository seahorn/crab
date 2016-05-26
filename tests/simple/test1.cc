#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
cfg_t* prog (VariableFactory &vfac)  {

  // Defining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  z_var nd (vfac ["nd"]);
  z_var inc (vfac ["inc"]);
  // entry and exit block
  auto cfg = new cfg_t("x0","ret");
  // adding blocks
  basic_block_t& x0 = cfg->insert ("x0");
  basic_block_t& x1 = cfg->insert ("x1");
  basic_block_t& x2 = cfg->insert ("x2");
  basic_block_t& x3 = cfg->insert ("x3");
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  x0 >> x1; x1 >> x2; x2 >> x3; x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  x0.assign (k, 2147483648);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.havoc(nd.name());
  bb2.select(inc,nd,1,2);
  bb2.add(i, i, inc);

  return cfg;
}

template <typename Domain>
void run(cfg_ref_t cfg, VariableFactory &vfac)
{
  typename NumFwdAnalyzer <cfg_ref_t, Domain, VariableFactory>::type 
      It (cfg, vfac, nullptr, 1, 2, 20);
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
  cfg_t* cfg = prog(vfac);
  cfg->simplify (); // this is optional
  crab::outs() << *cfg << "\n";

  run<interval_domain_t>(*cfg, vfac);
  run<dbm_domain_t>(*cfg, vfac);
  run<sdbm_domain_t>(*cfg, vfac);
  run<ric_domain_t>(*cfg, vfac);
  run<term_domain_t>(*cfg, vfac);
  run<dis_interval_domain_t>(*cfg, vfac);

  // free the CFG
  delete cfg;

  return 0;
}
