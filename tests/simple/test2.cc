#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t* prog (VariableFactory &vfac) 
{

  cfg_t* cfg = new cfg_t("loop1_entry","ret");
  //cfg_t cfg ("loop1_entry");
  basic_block_t& loop1_entry = cfg->insert ("loop1_entry");
  basic_block_t& loop1_bb1   = cfg->insert ("loop1_bb1");
  basic_block_t& loop1_bb1_t = cfg->insert ("loop1_bb1_t");
  basic_block_t& loop1_bb1_f = cfg->insert ("loop1_bb1_f");
  basic_block_t& loop1_bb2   = cfg->insert ("loop1_bb2");
  basic_block_t& loop2_entry = cfg->insert ("loop2_entry");
  basic_block_t& loop2_bb1   = cfg->insert ("loop2_bb1");
  basic_block_t& loop2_bb1_t = cfg->insert ("loop2_bb1_t");
  basic_block_t& loop2_bb1_f = cfg->insert ("loop2_bb1_f");
  basic_block_t& loop2_bb2   = cfg->insert ("loop2_bb2");
  basic_block_t& ret         = cfg->insert ("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var i(vfac["i"]);
  z_var j(vfac["j"]);
  z_var k(vfac["k"]);

  loop1_entry.assign (i, 0);
  loop1_entry.assign (k, 30);
  loop1_bb1_t.assume (i <= 9);
  loop1_bb1_f.assume (i >= 10);
  loop1_bb2.add (i, i, 1);

  loop2_entry.assign (j, 0);
  loop2_bb1_t.assume (j <= 9);
  loop2_bb1_f.assume (j >= 10);
  loop2_bb2.add (j, j, 1);
  return cfg;
}

template <typename Domain, typename Live>
void run(cfg_ref_t cfg, VariableFactory &vfac, Live *live)
{
  typename NumFwdAnalyzer <cfg_ref_t, Domain, VariableFactory>::type 
      It (cfg, vfac, live, 1, 2, 20);
  Domain inv = Domain::top ();
  It.Run (inv);
  crab::outs() << "Invariants using " << Domain::getDomainName () << ":\n";

  for (auto &b : cfg)
  {
    // invariants at the entry of the block
    auto inv = It [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  crab::outs() << endl;
  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  VariableFactory vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << endl;

  Liveness<cfg_ref_t> live (*cfg);
  live.exec ();

  run<interval_domain_t>(*cfg,vfac,&live);
  run<dbm_domain_t>(*cfg,vfac,&live);
  run<sdbm_domain_t>(*cfg,vfac,&live);
  run<ric_domain_t>(*cfg,vfac,&live);
  run<term_domain_t>(*cfg,vfac,&live);
  run<dis_interval_domain_t>(*cfg,vfac,&live);

  delete cfg;
  return 0;
}
