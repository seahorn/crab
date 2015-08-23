#include <ikos/tests/Cfg_impl.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/intervals_congruences.hpp>                      
#include <ikos/domains/octagons.hpp>                      
#include <ikos/domains/dbm.hpp>                      

using namespace std;

namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains
  typedef interval_domain< z_number, varname_t > interval_domain_t;
  typedef interval_congruence_domain< z_number, varname_t > ric_t;
  typedef DBM<z_number, varname_t> dbm_domain_t;
  typedef octagon< z_number, varname_t > octagon_domain_t;
} // end namespace

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;

cfg_t prog (VariableFactory &vfac) 
{

  cfg_t cfg ("loop1_entry","ret");
  //cfg_t cfg ("loop1_entry");
  basic_block_t& loop1_entry = cfg.insert ("loop1_entry");
  basic_block_t& loop1_bb1   = cfg.insert ("loop1_bb1");
  basic_block_t& loop1_bb1_t = cfg.insert ("loop1_bb1_t");
  basic_block_t& loop1_bb1_f = cfg.insert ("loop1_bb1_f");
  basic_block_t& loop1_bb2   = cfg.insert ("loop1_bb2");
  basic_block_t& loop2_entry = cfg.insert ("loop2_entry");
  basic_block_t& loop2_bb1   = cfg.insert ("loop2_bb1");
  basic_block_t& loop2_bb1_t = cfg.insert ("loop2_bb1_t");
  basic_block_t& loop2_bb1_f = cfg.insert ("loop2_bb1_f");
  basic_block_t& loop2_bb2   = cfg.insert ("loop2_bb2");
  basic_block_t& ret         = cfg.insert ("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var i(vfac["i"]);
  z_var j(vfac["j"]);
  z_var k(vfac["k"]);

  loop1_entry.assign (i, 0);
  loop1_entry.assign (k, 30);
  loop1_bb1_t.assume (i <= 9);
  loop1_bb1_f.assume (i >= 10);
  loop1_bb2.add (i, i, 1);

  loop2_entry.assign (j, 0);
  loop2_bb1_t.assume (j <= 9);
  loop2_bb1_f.assume (j >= 10);
  loop2_bb2.add (j, j, 1);
  return cfg;
}


int main (int argc, char** argv )
{


  VariableFactory vfac;
  cfg_t cfg = prog (vfac);
  cfg.simplify ();
  cout << cfg << endl;

  const bool run_live = true;

  NumFwdAnalyzer <cfg_t, interval_domain_t, VariableFactory>::type itv_a (cfg,vfac,run_live);
  itv_a.Run (interval_domain_t::top ());
  cout << "Results with intervals:\n";
  for (auto &b : cfg)
  {
    interval_domain_t inv = itv_a [b.label ()];
    std::cout << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
  }

  NumFwdAnalyzer <cfg_t, dbm_domain_t, VariableFactory>::type dbm_a (cfg,vfac,run_live);
  dbm_a.Run (dbm_domain_t::top ());
  cout << "Results with DBMs:\n";
  for (auto &b : cfg)
  {
    dbm_domain_t inv = dbm_a [b.label ()];
    std::cout << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
  }

  return 0;
}
