#include <iostream>
#include <ikos/muaz_analyzer.hpp>

using namespace ikos;
using namespace ikos_impl;

cfg_t prog (VariableFactory &vfac) 
{

  cfg_t cfg("loop1_entry","ret");
  basic_block_t& loop1_entry = cfg.insert_basic_block("loop1_entry");
  basic_block_t& loop1_bb1   = cfg.insert_basic_block("loop1_bb1");
  basic_block_t& loop1_bb1_t = cfg.insert_basic_block("loop1_bb1_t");
  basic_block_t& loop1_bb1_f = cfg.insert_basic_block("loop1_bb1_f");
  basic_block_t& loop1_bb2   = cfg.insert_basic_block("loop1_bb2");
  basic_block_t& loop2_entry = cfg.insert_basic_block("loop2_entry");
  basic_block_t& loop2_bb1   = cfg.insert_basic_block("loop2_bb1");
  basic_block_t& loop2_bb1_t = cfg.insert_basic_block("loop2_bb1_t");
  basic_block_t& loop2_bb1_f = cfg.insert_basic_block("loop2_bb1_f");
  basic_block_t& loop2_bb2   = cfg.insert_basic_block("loop2_bb2");
  basic_block_t& ret   = cfg.insert_basic_block("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var j(vfac["j"]);

  loop1_entry.assign(n1, 1);
  loop1_entry.assign(i, 0);
  loop1_bb1_t.assertion(i <= 9);
  loop1_bb1_f.assertion(i >= 10);
  loop1_bb2.add(i, i, n1);
  loop1_bb2.check("end_of_loop");

  loop2_entry.assign(j, 0);
  loop2_bb1_t.assertion(j <= 9);
  loop2_bb1_f.assertion(j >= 10);
  loop2_bb2.add(j, j, n1);
  loop2_bb2.check("end_of_loop");
  ret.check("end_of_program");
  return cfg;
}


int main (int argc, char** argv )
{


  VariableFactory vfac;
  cfg_t CFG = prog (vfac);
  cout << "CFG: \n"  << CFG << endl;

  interval_domain_t intervals = interval_domain_t::top ();
  cout << "Running " << intervals.getDomainName () << endl;
  FwdAnalyzer <interval_domain_t>  itv_a (CFG, vfac);
  itv_a.run (intervals);

  dbm_domain_t dbm = dbm_domain_t::top ();
  cout << "Running " << dbm.getDomainName () << endl;
  FwdAnalyzer <dbm_domain_t>  dbm_a (CFG, vfac);
  dbm_a.run (dbm);

  return 0;
}
