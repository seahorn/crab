#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t* prog (variable_factory_t &vfac) 
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


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";

  run<interval_domain_t>(cfg,vfac,true,1,2,20);
  run<dbm_domain_t>(cfg,vfac,true,1,2,20);
  run<sdbm_domain_t>(cfg,vfac,true,1,2,20);
  run<ric_domain_t>(cfg,vfac,true,1,2,20);
  run<term_domain_t>(cfg,vfac,true,1,2,20);
  run<dis_interval_domain_t>(cfg,vfac,true,1,2,20);

  delete cfg;
  return 0;
}
