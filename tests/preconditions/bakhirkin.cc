/* 
   Fig.5 from SAS'17 paper "Combining Forward and Backward Abstract
   Interpretation of Horn Clauses" by Bakhirkin and Monniaux.
*/

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
  auto cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry   = cfg->insert ("entry");
  z_basic_block_t& header1 = cfg->insert ("header1");
  z_basic_block_t& body1   = cfg->insert ("body1");
  z_basic_block_t& exit1    = cfg->insert ("exit1");
  z_basic_block_t& ifxpos = cfg->insert ("ifxpos");  
  z_basic_block_t& header2 = cfg->insert ("header2");
  z_basic_block_t& body2   = cfg->insert ("body2");
  z_basic_block_t& exit2   = cfg->insert ("exit2");
  z_basic_block_t& ret     = cfg->insert ("ret");
  
  // adding control flow
  entry >> header1;
  header1 >> body1;
  body1 >> header1;
  header1 >> exit1;
  exit1 >> ifxpos;
  exit1 >> ret;
  ifxpos >> header2;
  header2 >> body2;
  body2 >> header2;
  header2 >> exit2;
  exit2 >> ret;

  
  /*
    x=0;y=*;
    while(*) 
      x += y;
    if (x>0) {
      while(*)
        y += x;
      assert(y>=0);
    }
  */
  entry.assign(x, 0);
  body1.add (x, x, y);
  ifxpos.assume(x >= 1);
  body2.add(y,y,x);
  exit2.assertion(y>=0);
  return cfg;
}


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  if (false) {
    z_box_apron_domain_t initial_states, final_states;
    final_states = z_box_apron_domain_t::bottom();    
    backward_run<z_box_apron_domain_t>
      (cfg, initial_states, final_states, 1, 2, 20, stats_enabled);
  }
  
  if (true) {
    z_pk_apron_domain_t initial_states, final_states;
    final_states = z_pk_apron_domain_t::bottom();
    backward_run<z_pk_apron_domain_t>
      (cfg, initial_states, final_states, 1, 2, 20, stats_enabled);
  }
  
  // free the CFG
  delete cfg;

  return 0;
}
