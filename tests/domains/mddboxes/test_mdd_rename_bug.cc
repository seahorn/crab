#include "crab/config.h"

#ifndef HAVE_MDD
int main(int argc, char**argv) {
  return 0;
}
#else
#include "../../program_options.hpp"
#include "../../common.hpp"
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog1 (variable_factory_t &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"], crab::INT_TYPE, 32);
  z_var k (vfac ["k"], crab::INT_TYPE, 32);
  z_var nd (vfac ["nd"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign (k, 0);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

int main (int argc, char** argv ) { 
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  z_cfg_t* cfg1 = prog1 (vfac);      
  crab::outs() << *cfg1 << "\n";
  run<z_mdd_boxes_domain_t>(cfg1,cfg1->entry(),false,10,2,20,stats_enabled);
  delete cfg1;
  return 0;
}

#endif

