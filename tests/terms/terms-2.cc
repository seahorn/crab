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
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  z_var z0 (vfac ["z0"]);
  z_var y0 (vfac ["y0"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("p0","ret");
  // adding blocks
  basic_block_t& p0 = cfg->insert ("p0");
  basic_block_t& p_neg = cfg->insert ("p_neg");
  basic_block_t& p_pos = cfg->insert ("p_pos");
  basic_block_t& exit = cfg->insert ("exit");
  basic_block_t& ret = cfg->insert ("ret");
  // adding control flow
  p0 >> p_pos;
  p0 >> p_neg;
  p_neg >> exit;
  p_pos >> exit;
  exit >> ret;

  // adding statements
  p0.assign (x, 50);
  p0.havoc (y.name());
  p0.assume (y >= -1);
  p0.assume (y <= 1);
  p0.mul (z, x, y);

  p_neg.assume (y <= -1);
  p_neg.mul (z, z, -1);

  p_pos.assume(y >= 0);

  exit.assign(z0, z);
  exit.assign(y0, y);

  return cfg;
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";

  run<interval_domain_t>(cfg, vfac, true, 1, 2, 20);
  run<term_domain_t>(cfg, vfac, true, 1, 2, 20);

  /*
  term_dbm_t t_dbm_dom = term_dbm_t::top ();
  num_fwd_analyzer <cfg_ref_t, term_dbm_t, variable_factory_t>::type
      term_dbm_a (*cfg, vfac, &live);
  term_dbm_a.Run (t_dbm_dom);
  crab::outs() << "Results with term<dbm> domain:\n";
  for (auto &b : *cfg)
  {
    term_dbm_t inv = term_dbm_a [b.label ()];
    crab::outs() << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
//    auto cst = inv.to_linear_constraint_system();
//    crab::outs() << "  := " << cst << "\n";
  }


  dbm_domain_t dbm = dbm_domain_t::top ();
  num_fwd_analyzer <cfg_ref_t, dbm_domain_t, variable_factory_t>::type
      dbm_a (*cfg, vfac, &live);
  dbm_a.Run (dbm);
  crab::outs() << "Results with DBMs:\n";
  for (auto &b : *cfg)
  {
    dbm_domain_t inv = dbm_a [b.label ()];
    crab::outs() << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
  }
  */

  return 0;
}
