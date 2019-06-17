#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog (variable_factory_t &vfac) 
{

  z_cfg_t* cfg = new z_cfg_t("entry");
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& l1 = cfg->insert ("l1");
  z_basic_block_t& l1_entry = cfg->insert ("l1_entry");
  z_basic_block_t& l1_reset = cfg->insert ("l1_reset");
  z_basic_block_t& l1_dont_reset = cfg->insert ("l1_dont_reset");
  z_basic_block_t& l2 = cfg->insert ("l2");
  z_basic_block_t& l2_body = cfg->insert ("l2_body");
  z_basic_block_t& l2_exit   = cfg->insert ("l2_exit");

  
  entry >> l1;
  // outer loop
  l1 >> l1_entry;
  l1_entry >> l2;
  // inner loop
  l2 >> l2_body;
  l2_body >> l2;
  l2 >> l2_exit;
  // outer loop again
  l2_exit >> l1_reset;
  l2_exit >> l1_dont_reset;
  l1_reset >> l1;
  l1_dont_reset >> l1;

  
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  l1_entry.add(i, i, 1);
  l1_entry.assign(j, 0);

  l2_body.assume(j <= 9);
  l2_body.assertion(i >= 0);
  l2_body.assertion(i <= 10);
  l2_body.add(j, j, 1);

  l2_exit.assume(j >= 10);
  
  l1_reset.assume(i >= 10);
  l1_reset.assign(i, 0);
  l1_dont_reset.assume(i <= 9);
  
  return cfg;
}


int main (int argc, char** argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  variable_factory_t vfac;
  z_cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";

  run_and_check<z_interval_domain_t>(cfg,cfg->entry(),false,1,2,20,stats_enabled);

  delete cfg;
  return 0;
}
