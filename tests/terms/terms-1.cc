#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t* prog (variable_factory_t &vfac) 
{
  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  // entry and exit block
  cfg_t* cfg  = new cfg_t("x0","ret");
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
  x0.assign (k, 50);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.add(i, i, 1);

  return cfg;
}


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";

  num_fwd_analyzer <cfg_ref_t, interval_domain_t, variable_factory_t>::type itv_a (*cfg, vfac, nullptr);
  itv_a.Run (interval_domain_t::top ());
  crab::outs() << "Results with intervals:\n";
  for (auto &b : *cfg)
  {
    interval_domain_t inv = itv_a [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }

  num_fwd_analyzer <cfg_ref_t, term_domain_t, variable_factory_t>::type term_a (*cfg, vfac, nullptr);
  term_a.Run (term_domain_t::top ());
  crab::outs() << "Results with term<interval> domain:\n";
  for (auto &b : *cfg)
  {
    term_domain_t inv = term_a [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }

  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }

  crab::outs() << "\n  as linear constraints:\n" << "\n";
  for (auto & b : *cfg)
  {
    term_domain_t inv = term_a [b.label ()];
    term_domain_t::linear_constraint_system_t cst(inv.to_linear_constraint_system());
    crab::outs() << "  " << get_label_str (b.label ()) << "=" << cst << "\n";
  }

  /*
  NumFwdAnalyzer <cfg_ref_t, dbm_domain_t, variable_factory_t>::type dom_a (*cfg, vfac, nullptr);

  dbm_a.Run (dbm_domain_t::top ());
  crab::outs() << "Results with DBMs:\n";
  for (auto &b : *cfg)
  {
    dbm_domain_t inv = dbm_a [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  */

  return 0;
}
