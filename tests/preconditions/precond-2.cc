#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t* prog (variable_factory_t &vfac)  {

  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var tmp (vfac ["tmp"]);    
  // entry and exit block
  auto cfg = new z_cfg_t("entry","bb3");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& bb3   = cfg->insert ("bb3");
  // adding control flow
  entry >> bb1;
  entry >> bb2;
  bb1 >> bb3;
  bb2 >> bb3;
  
  // adding statements
  bb1.assume (x >= 1);
  bb1.assign (y, x);
  
  bb2.assume (x <= 0);
  bb2.assign (tmp, 0);
  bb2.sub (y, tmp, x);
  
  bb3.assume (y >= 6);

  return cfg;
}


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  {
    z_box_apron_domain_t initial_states, final_states;
    z_var x (vfac ["x"]);  
    final_states += (x <= 5);
    final_states += (x >= -5);    
    backward_run<z_box_apron_domain_t>
      (cfg, initial_states, final_states, 1, 2, 20, stats_enabled);
  }

  {
    z_boxes_domain_t initial_states, final_states;
    z_var x (vfac ["x"]);  
    final_states += (x <= 5);
    final_states += (x >= -5);    
    backward_run<z_boxes_domain_t>
      (cfg, initial_states, final_states, 1, 2, 20, stats_enabled);
  }
  
  // free the CFG
  delete cfg;

  return 0;
}
