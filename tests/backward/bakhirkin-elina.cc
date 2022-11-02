/*
   Fig.5 from SAS'17 paper "Combining Forward and Backward Abstract
   Interpretation of Horn Clauses" by Bakhirkin and Monniaux.
*/

#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header1 = cfg->insert("header1");
  z_basic_block_t &body1 = cfg->insert("body1");
  z_basic_block_t &exit1 = cfg->insert("exit1");
  z_basic_block_t &ifxpos = cfg->insert("ifxpos");
  z_basic_block_t &header2 = cfg->insert("header2");
  z_basic_block_t &body2 = cfg->insert("body2");
  z_basic_block_t &exit2 = cfg->insert("exit2");
  z_basic_block_t &ret = cfg->insert("ret");

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
  body1.add(x, x, y);
  ifxpos.assume(x >= 1);
  body2.add(y, y, x);
  exit2.assertion(y >= 0);
  return cfg;
}

int main(int argc, char **argv) {
#ifdef HAVE_ELINA
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  z_oct_elina_domain_t initial_states;
  backward_run<z_oct_elina_domain_t>(cfg, cfg->entry(), initial_states, 1, 2,
				     20, stats_enabled);
  // free the CFG
  delete cfg;
#endif

  return 0;
}
