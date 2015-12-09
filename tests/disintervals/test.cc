#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t prog (VariableFactory &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  // entry and exit block
  cfg_t cfg ("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& bb1   = cfg.insert ("bb1");
  basic_block_t& bb1_t = cfg.insert ("bb1_t");
  basic_block_t& bb1_f = cfg.insert ("bb1_f");
  basic_block_t& bb2   = cfg.insert ("bb2");
  basic_block_t& ret   = cfg.insert ("ret");
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

int main (){

  // typedef dis_interval <z_number> dis_interval_t;
  // typedef interval <z_number> interval_t;

  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   std::cout << x << "\n";
  //   dis_interval_t y (interval_t (10,20));
  //   std::cout << y << "\n";
  //   dis_interval_t z (interval_t (7,8));
  //   std::cout << z << "\n";
    
  //   dis_interval_t r = (x | y) | z;
  //   std::cout << r << "\n";

  //   dis_interval_t s (interval_t (2,3));
  //   std::cout << "after adding " << s << "=" ;
  //   r = r + s;
  //   std::cout << r << "\n";
  // }

  // {  
  //   interval_t zero (5,5);
  //   dis_interval_t x (zero.lower_half_line ());
  //   std::cout << x << "\n";
  //   dis_interval_t y (zero.upper_half_line ());
  //   std::cout << y << "\n";
  //   dis_interval_t r = (x | y);
  //   std::cout << r << "\n";
  // }


  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   std::cout << x << "\n";
  //   dis_interval_t y (interval_t (10,20));
  //   std::cout << y << "\n";
  //   dis_interval_t z (interval_t (7,8));
  //   std::cout << z << "\n";
    
  //   dis_interval_t r = (x & y) & z;
  //   std::cout << r << "\n";
  // }

  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   std::cout << x << "\n";
  //   dis_interval_t y (interval_t (2,9));
  //   std::cout << y << "\n";
  //   dis_interval_t z (interval_t (2,8));
  //   std::cout << z << "\n";
    
  //   dis_interval_t r = (x & y) & z;
  //   std::cout << r << "\n";
  // }

  {
    VariableFactory vfac;
    cfg_t cfg = prog (vfac);
    cout << cfg << "\n";

    NumFwdAnalyzer <cfg_t, dis_interval_domain_t,VariableFactory>::type 
        a (cfg, vfac, nullptr, 1, 2, 20);
    // Run fixpoint 
    a.Run (dis_interval_domain_t::top ());
    // Print invariants
    cout << "Invariants using " << dis_interval_domain_t::getDomainName () << "\n";
    for (auto &b : cfg) {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }
}
