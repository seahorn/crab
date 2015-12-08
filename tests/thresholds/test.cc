#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t prog (VariableFactory &vfac) 
{

  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var n (vfac ["n"]);
  z_var x (vfac ["x"]);
  // entry and exit block
  cfg_t cfg ("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& loop_header   = cfg.insert ("loop_header");
  basic_block_t& loop_bb0 = cfg.insert ("loop_bb0");
  basic_block_t& loop_bb1 = cfg.insert ("loop_bb1");
  basic_block_t& loop_bb1_t = cfg.insert ("loop_bb1_t");
  basic_block_t& loop_bb1_f = cfg.insert ("loop_bb1_f");
  basic_block_t& loop_bb2 = cfg.insert ("loop_bb2");
  basic_block_t& loop_bb3 = cfg.insert ("loop_bb3");
  basic_block_t& loop_bb4 = cfg.insert ("loop_bb4");
  basic_block_t& ret   = cfg.insert ("ret");
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

int main (){
  
  // Simple operations on thresholds
  typename cfg_t::thresholds_t thresholds;
  thresholds.add (-89);
  thresholds.add (1000);
  thresholds.add (100);
  thresholds.add (-999);
  thresholds.add (5);
  thresholds.add (10);
  thresholds.add (-10);

  cout << "Thresholds= " << thresholds << "\n";

  auto t1 = thresholds.get_next (3);
  std::cout << "next threshold for 3:" << t1 << "\n";
  auto t2 = thresholds.get_next (8);
  std::cout << "next threshold for 8: " << t2 << "\n";
  auto t3 = thresholds.get_next (100);
  std::cout << "next threshold for 100: " << t3 << "\n";
  auto t4 = thresholds.get_next (500);
  std::cout << "next threshold for 500: " << t4 << "\n";
  auto t5 = thresholds.get_next (10000);
  std::cout << "next threshold for 10000: " << t5 << "\n";
  auto t6 = thresholds.get_next (-4);
  std::cout << "next threshold for -4: " << t6 << "\n";
  auto t7 = thresholds.get_prev (-4);
  std::cout << "prev threshold for -4:" << t7 << "\n";
  auto t8 = thresholds.get_prev (-78);
  std::cout << "prev threshold for -78:" << t8  << "\n";
  auto t9 = thresholds.get_prev (-10000);
  std::cout << "prev threshold for -10000:" << t9 << "\n";

  {
    VariableFactory vfac;
    cfg_t cfg = prog (vfac);
    cout << cfg << "\n";

    typename cfg_t::thresholds_t ts = cfg.initialize_thresholds_for_widening (50);
    cout << "Thresholds=" << ts << "\n";

    NumFwdAnalyzer <cfg_t, term_domain_t,VariableFactory>::type 
        a (cfg, vfac, nullptr, 1, 2, 20);
    // Run fixpoint 
    term_domain_t inv = term_domain_t::top ();
    a.Run (inv);
    // Print invariants
    cout << "Invariants using " << term_domain_t::getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }
}
