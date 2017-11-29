#include "../program_options.hpp"
#include "../common.hpp"
#include "crab/iterators/fwd_fixpoint_iterators.hpp"

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
  z_var n (vfac ["n"], crab::INT_TYPE);
  z_var x (vfac ["x"], crab::INT_TYPE);
  z_var y (vfac ["y"], crab::INT_TYPE);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry","ret");
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& loop1_header = cfg->insert ("loop1_header");
  z_basic_block_t& loop1_bb0 = cfg->insert ("loop1_bb0");
  z_basic_block_t& loop1_bb1 = cfg->insert ("loop1_bb1");
  z_basic_block_t& loop1_bb1_t = cfg->insert ("loop1_bb1_t");
  z_basic_block_t& loop1_bb1_f = cfg->insert ("loop1_bb1_f");
  z_basic_block_t& loop1_bb2 = cfg->insert ("loop1_bb2");
  z_basic_block_t& loop1_bb3 = cfg->insert ("loop1_bb3");
  z_basic_block_t& loop1_bb4 = cfg->insert ("loop1_bb4");
  z_basic_block_t& cont = cfg->insert ("cont");  
  z_basic_block_t& loop2_header = cfg->insert ("loop2_header");
  z_basic_block_t& loop2_bb0 = cfg->insert ("loop2_bb0");
  z_basic_block_t& loop2_bb1 = cfg->insert ("loop2_bb1");
  z_basic_block_t& loop2_bb1_t = cfg->insert ("loop2_bb1_t");
  z_basic_block_t& loop2_bb1_f = cfg->insert ("loop2_bb1_f");
  z_basic_block_t& loop2_bb2 = cfg->insert ("loop2_bb2");
  z_basic_block_t& loop2_bb3 = cfg->insert ("loop2_bb3");
  z_basic_block_t& loop2_bb4 = cfg->insert ("loop2_bb4");
  z_basic_block_t& ret = cfg->insert ("ret");
  // adding control flow
  entry >> loop1_header;

  loop1_header >> loop1_bb0; 
  loop1_header >> loop1_bb1;
  loop1_bb0 >> loop1_header;
  loop1_bb1 >> loop1_bb1_t;  loop1_bb1 >> loop1_bb1_f; 
  loop1_bb1_t >> loop1_bb2;  loop1_bb1_f >> loop1_bb3; 
  loop1_bb2 >> loop1_bb4;
  loop1_bb3 >> loop1_bb4;
  loop1_bb4 >> loop1_header;
  loop1_header >> cont;

  cont >> loop2_header;
  loop2_header >> loop2_bb0;
  loop2_header >> loop2_bb1;  
  loop2_bb0 >> loop2_header;
  loop2_bb1 >> loop2_bb1_t;  loop2_bb1 >> loop2_bb1_f; 
  loop2_bb1_t >> loop2_bb2;  loop2_bb1_f >> loop2_bb3; 
  loop2_bb2 >> loop2_bb4;
  loop2_bb3 >> loop2_bb4;
  loop2_bb4 >> loop2_header;
  loop2_header >> ret;
  
  // adding statements
  entry.assign (n , 0);
  entry.assume(y <= 878);
  
  loop1_bb0.havoc (x);
  loop1_bb0.assume (x >= 0);
  loop1_bb1_t.assume (n <= 60);
  loop1_bb1_f.assume (n >= 61);
  loop1_bb2.add(n, n, 1);
  loop1_bb3.assign(n, 0);

  cont.assign (n,0);
  loop2_bb0.havoc (x);
  loop2_bb0.assume (x >= 10);
  loop2_bb1_t.assume (n <= 160);
  loop2_bb1_f.assume (n >= 161);
  loop2_bb2.add(n, n, 2);
  loop2_bb3.assign(n, 0);
  

  return cfg;
}


int main (int argc, char**argv){

  SET_TEST_OPTIONS(argc,argv)
  
  // Simple operations on thresholds
  crab::iterators::thresholds<ikos::z_number> thresholds;
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

    // Print thresholds 
    typedef crab::cfg::cfg_ref<z_cfg_t> z_cfg_ref_t;
    typedef ikos::wto<typename z_cfg_t::basic_block_label_t, z_cfg_ref_t> wto_t;
    typedef crab::iterators::wto_thresholds<typename z_cfg_t::basic_block_label_t, z_cfg_ref_t> wto_thresholds_t;
    
    z_cfg_ref_t cfg_ref(*cfg);
    wto_t wto(cfg_ref);
    wto_thresholds_t thresholds(cfg_ref, 50);
    wto.accept(&thresholds);
    crab::outs() << "Thresholds per loop \n" << thresholds << "\n";

    
    ////////////////////////////////////////////////////////////
    // -- w/o thresholds
    ////////////////////////////////////////////////////////////    
    run<z_term_domain_t>(cfg, false, 1, 2, 0, stats_enabled);

    ////////////////////////////////////////////////////////////
    // -- w/thresholds
    ////////////////////////////////////////////////////////////        
    run<z_term_domain_t>(cfg, false, 1, 2, 20, stats_enabled);


    delete cfg;
  }
}
