#include <iostream>
#include <ikos/muaz_analyzer.hpp>

using namespace ikos;
using namespace ikos_impl;

cfg_t prog (VariableFactory &vfac) 
{

  // Definining program variables
  z_var n1 (vfac["n1"]);
  z_var i (vfac["i"]);
  ////
  // Building the CFG
  ////

  // entry and exit block
  cfg_t cfg("entry","ret");
  // adding blocks
  basic_block_t& entry = cfg.insert_basic_block("entry");
  basic_block_t& bb1   = cfg.insert_basic_block("bb1");
  basic_block_t& bb1_t = cfg.insert_basic_block("bb1_t");
  basic_block_t& bb1_f = cfg.insert_basic_block("bb1_f");
  basic_block_t& bb2   = cfg.insert_basic_block("bb2");
  basic_block_t& ret   = cfg.insert_basic_block("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign(n1, 1);
  entry.assign(i, 0);
  bb1_t.assertion(i <= 99);
  bb1_f.assertion(i >= 100);
  bb2.add(i, i, n1);
  bb2.check("end_of_loop");
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
