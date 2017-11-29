#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog (variable_factory_t &vfac) 
{

  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  z_basic_block_t& entry       = cfg->insert ("entry");
  z_basic_block_t& loop1_head  = cfg->insert ("loop1_head");
  z_basic_block_t& loop1_t     = cfg->insert ("loop1_t");
  z_basic_block_t& loop1_f     = cfg->insert ("loop1_f");
  z_basic_block_t& loop1_body  = cfg->insert ("loop1_body");

  z_basic_block_t& loop1_body_t  = cfg->insert ("loop1_body_t");
  z_basic_block_t& loop1_body_f  = cfg->insert ("loop1_body_f");
  z_basic_block_t& loop1_body_x  = cfg->insert ("loop1_body_x");

  z_basic_block_t& cont        = cfg->insert ("cont");
  z_basic_block_t& loop2_head  = cfg->insert ("loop2_head");
  z_basic_block_t& loop2_t     = cfg->insert ("loop2_t");
  z_basic_block_t& loop2_f     = cfg->insert ("loop2_f");
  z_basic_block_t& loop2_body  = cfg->insert ("loop2_body");
  z_basic_block_t& ret         = cfg->insert ("ret");

  entry >> loop1_head;
  loop1_head >> loop1_t; 
  loop1_head >> loop1_f; 
  loop1_t >>    loop1_body; 

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t; 
  loop2_head >> loop2_f; 
  loop2_t >>    loop2_body; 
  loop2_body >> loop2_head;
  loop2_f >> ret;
  
  z_var i(vfac["i"], crab::INT_TYPE);

  entry.assign (i, 0);
  loop1_t.assume (i <= 10);
  loop1_f.assume (i >= 11);
  loop1_body.add (i, i, 1);

  loop1_body_t.assume (i >= 9);
  loop1_body_t.assign (i , 0);
  loop1_body_f.assume (i <= 8);

  loop2_t.assume (i <= 100);
  loop2_f.assume (i >= 101);
  loop2_body.sub (i, i, 1);
  return cfg;
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog (vfac);
  //cfg->simplify ();
  crab::outs() << *cfg << "\n";

  run<z_interval_domain_t>(cfg,true,1,2,20,stats_enabled);
  run<z_dbm_domain_t>(cfg,true,1,2,20,stats_enabled);
  run<z_sdbm_domain_t>(cfg,true,1,2,20,stats_enabled);
  run<z_ric_domain_t>(cfg,true,1,2,20,stats_enabled);
  run<z_term_domain_t>(cfg,true,1,2,20,stats_enabled);
  run<z_dis_interval_domain_t>(cfg,true,1,2,20,stats_enabled);
  delete cfg;
  return 0;
}
