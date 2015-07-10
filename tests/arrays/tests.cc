#include <ikos/tests/Cfg_impl.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/dbm.hpp>                      
#include <ikos/domains/array_graph.hpp>                      
#include <ikos/domains/array_smashing.hpp>                      

using namespace std;

namespace domain_impl
{
  using namespace cfg_impl;
  // Abstract domains
  typedef interval_domain< z_number, varname_t > interval_domain_t;
  typedef DBM< z_number, varname_t > dbm_domain_t;
  typedef array_graph_domain<dbm_domain_t,
                             z_number, varname_t,
                             interval_domain_t, false> array_graph_domain_t;
  typedef array_smashing <dbm_domain_t, z_number, varname_t> array_smashing_t;

} // end namespace

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;

cfg_t prog1 (VariableFactory &vfac) 
{
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A0"]);
  z_var a1(vfac["A1"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp5(vfac["tmp5"]);
  z_var tmp6(vfac["tmp6"]);

  cfg_t cfg ("entry","ret", MEM);
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& bb1   = cfg.insert ("bb1");
  basic_block_t& bb1_t = cfg.insert ("bb1_t");
  basic_block_t& bb1_f = cfg.insert ("bb1_f");
  basic_block_t& bb2   = cfg.insert ("bb2");
  basic_block_t& ret   = cfg.insert ("ret");

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 123456);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3); // initialized
  ret.array_load(tmp6, a, i);    // top
  return cfg;
}


cfg_t prog2(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp4(vfac["tmp4"]);
  z_var tmp5(vfac["tmp5"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store (a, i, 123456);
  bb2.assign(tmp1, i);
  bb2.add(tmp2, tmp1, n1);
  bb2.assign(i, tmp2);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp4, a, tmp3); // initialized
  ret.array_load(tmp5, a, i);    // top
  return cfg;
}

cfg_t prog3(VariableFactory &vfac) 
{
  cfg_t cfg("loop1_entry","ret",MEM);
  basic_block_t& loop1_entry = cfg.insert("loop1_entry");
  basic_block_t& loop1_bb1   = cfg.insert("loop1_bb1");
  basic_block_t& loop1_bb1_t = cfg.insert("loop1_bb1_t");
  basic_block_t& loop1_bb1_f = cfg.insert("loop1_bb1_f");
  basic_block_t& loop1_bb2   = cfg.insert("loop1_bb2");
  ///
  basic_block_t& loop2_entry = cfg.insert("loop2_entry");
  basic_block_t& loop2_bb1   = cfg.insert("loop2_bb1");
  basic_block_t& loop2_bb1_t = cfg.insert("loop2_bb1_t");
  basic_block_t& loop2_bb1_f = cfg.insert("loop2_bb1_f");
  basic_block_t& loop2_bb2   = cfg.insert("loop2_bb2");
  /// 
  basic_block_t& ret   = cfg.insert("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;
  /////

  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var j(vfac["j"]);
  z_var a(vfac["A"]);
  z_var b(vfac["B"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp4(vfac["tmp4"]);

  loop1_entry.array_init (a.name ());
  loop1_entry.array_init (b.name ());
  loop1_entry.assign(n1, 1);
  loop1_entry.assign(i, 0);
  loop1_bb1_t.assume(i <= 9);
  loop1_bb1_f.assume(i >= 10);
  loop1_bb2.array_store(a, i, 123456);
  loop1_bb2.add(i, i, n1);

  loop2_entry.assign(j, 0);
  loop2_bb1_t.assume(j <= 9);
  loop2_bb1_f.assume(j >= 10);
  loop2_bb2.array_load(tmp1, a, j);    
  loop2_bb2.array_store(b, j, tmp1);
  loop2_bb2.add(j, j, n1);

  ret.sub(tmp2, j, n1);
  ret.array_load(tmp3, b, tmp2); // initialized
  ret.array_load(tmp4, b, j);    // top
  return cfg;
}


cfg_t prog4(VariableFactory &vfac) 
{

  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var b(vfac["B"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp5(vfac["tmp5"]);
  z_var tmp6(vfac["tmp6"]);
  z_var x(vfac["x"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.array_init (b.name ());
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 8);
  bb2.array_store(b, i, 5);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3); 
  ret.array_load(tmp6, b, tmp3); 
  return cfg;
}

cfg_t prog5(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  bb2.array_store(a, i, 123456);
  bb2.add(i, i, n1);
  ret.sub(tmp1, i, n1);
  ret.array_load(tmp2, a, tmp1); // initialized
  return cfg;
}

cfg_t prog6(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var c1(vfac["#1"]);
  z_var c3(vfac["#3"]);
  z_var c5(vfac["#5"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var tmp(vfac["tmp"]);
  z_var tmp1_offset_2(vfac["tmp1_offset_2"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp4(vfac["tmp4"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assign(i, 0);
  entry.assign(c1, 1);
  entry.assign(c5, 5);
  entry.assign(c3, 3);
  ///////
  bb1_t.assume(i <= 4);
  bb1_f.assume(5 <= i);
  bb2.assign(tmp, i);
  bb2.mul(tmp1_offset_2, tmp, c1); // 32
  bb2.array_store(a, tmp1_offset_2, 123456);
  bb2.add(tmp2, i, c1);
  bb2.assign(i, tmp2);
  ret.array_load(tmp4, a, c3);     // 160
  return cfg;
}

cfg_t prog7(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp4(vfac["tmp4"]);
  z_var x(vfac["x"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assume(n >= 2);
  entry.assign(n1, 1);
  entry.assign(i , 0);
  entry.array_store(a, i, 89);
  entry.assign(i , 1);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  ///////
  bb2.sub(tmp1, i, n1);
  bb2.array_load(tmp2, a, tmp1); 
  bb2.array_store(a, i, tmp2);
  bb2.add(i, i, n1);
  ///////
  ret.sub(tmp3, n, n1);
  ret.array_load(tmp4, a, tmp3); 
  return cfg;
}

// Initialize only even positions
// TODO: need a reduced of octagons and congruences
cfg_t prog8(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret", MEM);
  basic_block_t& entry = cfg.insert("entry");
  basic_block_t& bb1   = cfg.insert("bb1");
  basic_block_t& bb1_t = cfg.insert("bb1_t");
  basic_block_t& bb1_f = cfg.insert("bb1_f");
  basic_block_t& bb2   = cfg.insert("bb2");
  basic_block_t& ret   = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var n2(vfac["n2"]);
  z_var i(vfac["i"]);
  z_var i1(vfac["i1"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp3(vfac["tmp3"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(n2, 2);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 123456);
  // If we comment these two lines then we do only initialization of
  // even positions.
  bb2.add(i1, i, n1);
  bb2.array_store(a, i1, 123);
  bb2.add(i, i, n2);
  ret.assign(tmp1, 6);
  ret.array_load(tmp2, a, tmp1); // initialized
  return cfg;

}


// this is the program init_rand from Gange et.al paper.
cfg_t prog9(VariableFactory &vfac) 
{
  cfg_t cfg("entry","ret",MEM);
  basic_block_t& entry   = cfg.insert("entry");
  basic_block_t& bb1     = cfg.insert("bb1");
  basic_block_t& bb1_t   = cfg.insert("bb1_t");
  basic_block_t& bb1_f1  = cfg.insert("bb1_f1");
  basic_block_t& bb1_f2  = cfg.insert("bb1_f2");
  basic_block_t& bb1_f   = cfg.insert("bb1_f");
  basic_block_t& bb2_a   = cfg.insert("bb2a");
  basic_block_t& bb2_b   = cfg.insert("bb2b");
  basic_block_t& ret     = cfg.insert("ret");
  z_var n1(vfac["n1"]);
  z_var i1(vfac["i1"]);
  z_var i2(vfac["i2"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f1; bb1 >> bb1_f2;
  bb1_f1 >> bb1_f;   bb1_f2 >> bb1_f; 
  bb1_t >> bb2_a; bb1_t >> bb2_b; bb2_a >> bb1; bb2_b >> bb1; bb1_f >> ret;
  ////////
  entry.array_init (a.name ());
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i1, 0);
  entry.assign(i2, 0);
  ///////
  // while (i1 < n && i2 < n){
  bb1_t.assume(i1 <= n -1);
  bb1_t.assume(i2 <= n -1);
  // if (*)
  bb2_a.array_store(a, i1, 123456);
  bb2_a.add(i1, i1, n1);
  // else
  bb2_b.array_store(a, i2, 9);
  bb2_b.add(i2, i2, n1);
  // } end while
  bb1_f1.assume(i1 >= n);
  bb1_f2.assume(i2 >= n);
  ret.sub(tmp1, n, n1);
  ret.array_load(tmp2, a, tmp1); // initialized
  return cfg;
}


template <typename ArrayDomain>
void run(cfg_t cfg, string name, VariableFactory &vfac)
{
  cout << "--- " << name  << endl;
  cfg.simplify ();
  cout << cfg << endl;
  
  const bool run_live = false;
  typedef NumAbsTransformer <varname_t, ArrayDomain> arr_transformer_t;
  FwdAnalyzer <cfg_t, arr_transformer_t, VariableFactory> It (cfg,vfac,run_live);
  ArrayDomain inv = ArrayDomain::top ();
  It.Run (inv);
  cout << "Results with " << inv.getDomainName () << ":\n";

  for (auto &b : cfg)
  {
    // invariants at the entry of the block
    auto inv = It [b.label ()];
    cout << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
  }
  cout << endl;
}

void test1(){
  VariableFactory vfac;
  cfg_t cfg = prog1(vfac);
  run<array_graph_domain_t> (cfg, "Program 1", vfac);
}

void test2(){
  VariableFactory vfac;
  cfg_t cfg = prog2(vfac);
  run<array_graph_domain_t>(cfg, "Program 2", vfac);
}

void test3(){
  VariableFactory vfac;
  cfg_t cfg = prog3(vfac);
  run<array_graph_domain_t>(cfg, "Program 3", vfac);
}

void test4(){
  VariableFactory vfac;
  cfg_t cfg = prog4(vfac);
  run<array_graph_domain_t>(cfg, "Program 4", vfac);
}

void test5(){
  VariableFactory vfac;
  cfg_t cfg = prog5(vfac);
  run<array_graph_domain_t>(cfg, "Program 5", vfac);
}

void test6(){
  VariableFactory vfac;
  cfg_t cfg = prog6(vfac);
  run<array_graph_domain_t>(cfg, "Program 6", vfac);
}

void test7(){
  VariableFactory vfac;
  cfg_t cfg = prog7(vfac);
  run<array_graph_domain_t>(cfg, "Program 7", vfac);
}

void test8(){
  VariableFactory vfac;
  cfg_t cfg = prog8(vfac);
  run<array_graph_domain_t>(cfg, "Program 8", vfac);
}


void test9(){
  VariableFactory vfac;
  cfg_t cfg = prog9(vfac);
  run<array_graph_domain_t>(cfg, "Program 9", vfac);
}

void test10(){
  VariableFactory vfac;
  cfg_t cfg = prog1(vfac);
  run<array_smashing_t>(cfg, "Program 1", vfac);
}

void test11(){
  VariableFactory vfac;
  cfg_t cfg = prog2(vfac);
  run<array_smashing_t>(cfg, "Program 2", vfac);
}

void test12(){
  VariableFactory vfac;
  cfg_t cfg = prog3(vfac);
  run<array_smashing_t>(cfg, "Program 3", vfac);
}

void test13(){
  VariableFactory vfac;
  cfg_t cfg = prog4(vfac);
  run<array_smashing_t>(cfg, "Program 4", vfac);
}

void test14(){
  VariableFactory vfac;
  cfg_t cfg = prog5(vfac);
  run<array_smashing_t>(cfg, "Program 5", vfac);
}

void test15(){
  VariableFactory vfac;
  cfg_t cfg = prog6(vfac);
  run<array_smashing_t>(cfg, "Program 6", vfac);
}

void test16(){
  VariableFactory vfac;
  cfg_t cfg = prog7(vfac);
  run<array_smashing_t>(cfg, "Program 7", vfac);
}

void test17(){
  VariableFactory vfac;
  cfg_t cfg = prog8(vfac);
  run<array_smashing_t>(cfg, "Program 8", vfac);
}


void test18(){
  VariableFactory vfac;
  cfg_t cfg = prog9(vfac);
  run<array_smashing_t>(cfg, "Program 9", vfac);
}




int main(int, char **) 
{
  // array smashing
  test10 ();
  test11 ();
  test12 ();
  test13 ();
  test14 ();
  test15 ();
  test16 ();
  test17 ();
  test18 ();
  // array graph
  test1 ();
  test2 ();
  test3 ();
  test4 ();
  test5 ();
  test6 ();
  test7 ();
  test8 ();
  test9 ();
  
  return 42;
}


