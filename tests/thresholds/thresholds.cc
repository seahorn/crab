#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog (variable_factory_t &vfac) 
{

  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var n (vfac ["n"]);
  z_var x (vfac ["x"]);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& loop_header = cfg->insert ("loop_header");
  z_basic_block_t& loop_bb0 = cfg->insert ("loop_bb0");
  z_basic_block_t& loop_bb1 = cfg->insert ("loop_bb1");
  z_basic_block_t& loop_bb1_t = cfg->insert ("loop_bb1_t");
  z_basic_block_t& loop_bb1_f = cfg->insert ("loop_bb1_f");
  z_basic_block_t& loop_bb2 = cfg->insert ("loop_bb2");
  z_basic_block_t& loop_bb3 = cfg->insert ("loop_bb3");
  z_basic_block_t& loop_bb4 = cfg->insert ("loop_bb4");
  z_basic_block_t& ret = cfg->insert ("ret");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_bb1; loop_header >> ret;
  loop_header >> loop_bb0;
  loop_bb0 >> loop_header;
  loop_bb1 >> loop_bb1_t;  loop_bb1 >> loop_bb1_f; 
  loop_bb1_t >> loop_bb2;  loop_bb1_f >> loop_bb3; 
  loop_bb2 >> loop_bb4;
  loop_bb3 >> loop_bb4;
  loop_bb4 >> loop_header;

  // adding statements
  entry.assign (n , 0);
  loop_bb0.havoc (x.name ());
  loop_bb0.assume (x >= 0);
  loop_bb1_t.assume (n <= 60);
  loop_bb1_f.assume (n >= 61);
  loop_bb2.add(n, n, 1);
  loop_bb3.assign(n, 0);

  return cfg;
}

int main (int argc, char**argv){

  SET_TEST_OPTIONS(argc,argv)
  
  // Simple operations on thresholds
  typename z_cfg_t::thresholds_t thresholds;
  thresholds.add (z_bound(-89));
  thresholds.add (z_bound(1000));
  thresholds.add (z_bound(100));
  thresholds.add (z_bound(-999));
  thresholds.add (z_bound(5));
  thresholds.add (z_bound(10));
  thresholds.add (z_bound(-10));

  crab::outs() << "Thresholds= " << thresholds << "\n";

  auto t1 = thresholds.get_next (z_bound(3));
  crab::outs() << "next threshold for 3:" << t1 << "\n";
  auto t2 = thresholds.get_next (z_bound(8));
  crab::outs() << "next threshold for 8: " << t2 << "\n";
  auto t3 = thresholds.get_next (z_bound(100));
  crab::outs() << "next threshold for 100: " << t3 << "\n";
  auto t4 = thresholds.get_next (z_bound(500));
  crab::outs() << "next threshold for 500: " << t4 << "\n";
  auto t5 = thresholds.get_next (z_bound(10000));
  crab::outs() << "next threshold for 10000: " << t5 << "\n";
  auto t6 = thresholds.get_next (z_bound(-4));
  crab::outs() << "next threshold for -4: " << t6 << "\n";
  auto t7 = thresholds.get_prev (z_bound(-4));
  crab::outs() << "prev threshold for -4:" << t7 << "\n";
  auto t8 = thresholds.get_prev (z_bound(-78));
  crab::outs() << "prev threshold for -78:" << t8  << "\n";
  auto t9 = thresholds.get_prev (z_bound(-10000));
  crab::outs() << "prev threshold for -10000:" << t9 << "\n";

  {
    variable_factory_t vfac;
    z_cfg_t* cfg = prog (vfac);
    crab::outs() << *cfg << "\n";

    // w/o thresholds
    run<z_term_domain_t>(cfg, false, 1, 2, 0, stats_enabled);      

    // w/thresholds
    typename z_cfg_t::thresholds_t ts = cfg->initialize_thresholds_for_widening (50);
    crab::outs() << "Thresholds=" << ts << "\n";
    run<z_term_domain_t>(cfg, false, 1, 2, 20, stats_enabled);


    delete cfg;
  }
}
