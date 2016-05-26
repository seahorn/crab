#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t* prog (VariableFactory &vfac)  {

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
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

template <typename Domain, typename Live>
void run(cfg_ref_t cfg, VariableFactory &vfac, Live live)
{
  const unsigned int w = 1;
  const unsigned int n = 2;

  typename NumFwdAnalyzer <cfg_ref_t, Domain, VariableFactory>::type 
      It (cfg, vfac, live, w, n, 20);
  Domain inv = Domain::top ();
  It.Run (inv);
  crab::outs() << "Invariants using " << Domain::getDomainName () << ":\n";

  for (auto &b : cfg)
  {
    // invariants at the entry of the block
    auto inv = It [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  crab::outs() << "\n";
  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}

int main (int argc, char**argv){

  SET_TEST_OPTIONS(argc,argv)

  // typedef dis_interval <z_number> dis_interval_t;
  // typedef interval <z_number> interval_t;

  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   crab::outs() << x << "\n";
  //   dis_interval_t y (interval_t (10,20));
  //   crab::outs() << y << "\n";
  //   dis_interval_t z (interval_t (7,8));
  //   crab::outs() << z << "\n";
    
  //   dis_interval_t r = (x | y) | z;
  //   crab::outs() << r << "\n";

  //   dis_interval_t s (interval_t (2,3));
  //   crab::outs() << "after adding " << s << "=" ;
  //   r = r + s;
  //   crab::outs() << r << "\n";
  // }

  // {  
  //   interval_t zero (5,5);
  //   dis_interval_t x (zero.lower_half_line ());
  //   crab::outs() << x << "\n";
  //   dis_interval_t y (zero.upper_half_line ());
  //   crab::outs() << y << "\n";
  //   dis_interval_t r = (x | y);
  //   crab::outs() << r << "\n";
  // }


  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   crab::outs() << x << "\n";
  //   dis_interval_t y (interval_t (10,20));
  //   crab::outs() << y << "\n";
  //   dis_interval_t z (interval_t (7,8));
  //   crab::outs() << z << "\n";
    
  //   dis_interval_t r = (x & y) & z;
  //   crab::outs() << r << "\n";
  // }

  // {  
  //   dis_interval_t x (interval_t (0,5));
  //   crab::outs() << x << "\n";
  //   dis_interval_t y (interval_t (2,9));
  //   crab::outs() << y << "\n";
  //   dis_interval_t z (interval_t (2,8));
  //   crab::outs() << z << "\n";
    
  //   dis_interval_t r = (x & y) & z;
  //   crab::outs() << r << "\n";
  // }

  {
    VariableFactory vfac;
    cfg_t* cfg = prog (vfac);
    crab::outs() << *cfg << "\n";
    run<dis_interval_domain_t>(*cfg, vfac, nullptr); 
    delete cfg;
  }
}
