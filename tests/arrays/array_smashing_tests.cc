#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

cfg_t* prog1 (variable_factory_t &vfac) 
{
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A0"]);
  z_var a1(vfac["A1"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp5(vfac["tmp5"]);
  z_var tmp6(vfac["tmp6"]);

  cfg_t* cfg = new cfg_t("entry","ret", ARR);
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");

  // assume array element of 1 byte

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.assume_array (a.name (), 0);
  /////////
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 123456, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3, 1); // initialized
  ret.array_load(tmp6, a, i, 1);    // top
  return cfg;
}


cfg_t* prog2(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
  z_var n0(vfac["n0"]);
  z_var n1(vfac["n1"]);
  z_var n9(vfac["n9"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp4(vfac["tmp4"]);
  z_var tmp5(vfac["tmp5"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);
  /////////
  entry.assign(n0, 0); // we need it to be considered as graph node
  entry.assign(n1, 1); 
  entry.assign(n9, 9); // we need it to be considered as graph node
  entry.assign(i, n9);
  ///////
  bb1_t.assume(i >= 0);
  bb1_f.assume(i <= -1);
  bb2.array_store (a, i, 123456, 1);
  bb2.sub(i, i, n1);
  ret.assign(tmp3, 5);
  ret.array_load(tmp4, a, tmp3, 1); // initialized
  ret.array_load(tmp5, a, i, 1);    // top
  return cfg;
}

cfg_t* prog3(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;

  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var b(vfac["B"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp3(vfac["tmp3"]);
  z_var tmp4(vfac["tmp4"]);

  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);
  entry.assume_array (b.name (), 0);

  entry.assign(n1, 1);
  entry.assign(i, 0);
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 123456, 1);
  bb2.array_load(tmp1, a, i, 1);    
  bb2.array_store(b, i, tmp1, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp2, i, n1);
  ret.array_load(tmp3, b, tmp2, 1); // initialized
  ret.array_load(tmp4, b, i, 1);    // top
  return cfg;
}


cfg_t* prog4(variable_factory_t &vfac) 
{

  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
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

  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);
  entry.assume_array (b.name (), 0);

  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 8, 1);
  bb2.array_store(b, i, 5, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3, 1); 
  ret.array_load(tmp6, b, tmp3, 1); 
  return cfg;
}

cfg_t* prog5(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;

  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);

  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  bb2.array_store(a, i, 123456, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp1, i, n1);
  ret.array_load(tmp2, a, tmp1, 1); // initialized
  return cfg;
}

cfg_t* prog6(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
  z_var i(vfac["i"]);
  z_var a(vfac["A"]);
  z_var tmp(vfac["tmp"]);
  z_var offset(vfac["o"]);
  z_var tmp2(vfac["tmp2"]);
  z_var tmp4(vfac["tmp4"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // assume array element of 4 bytes
  entry.assume_array (a.name (), 0);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(10 <= i);
  bb2.assign(tmp, i);
  bb2.mul(offset, tmp, 4); 
  bb2.array_store(a, offset, 123456, 4);
  bb2.add(i, i, 1);
  ret.array_load(tmp4, a, 8, 4);    
  return cfg;
}

cfg_t* prog7(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
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
  // assume array element of 1 byte

  // assume (forall i. a[i] =0);
  entry.assume_array (a.name (), 0);
  //////
  entry.assume(n >= 2);
  entry.assign(n1, 1);
  entry.assign(i , 0);
  entry.array_store(a, i, 89, 1);
  entry.assign(i , 1);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  ///////
  bb2.sub(tmp1, i, n1);
  bb2.array_load(tmp2, a, tmp1, 1); 
  bb2.array_store(a, i, tmp2, 1);
  bb2.add(i, i, n1);
  ///////
  ret.sub(tmp3, n, n1);
  ret.array_load(tmp4, a, tmp3, 1); 
  return cfg;
}

// Initialize only even positions
cfg_t* prog8(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret", ARR);
  basic_block_t& entry = cfg->insert("entry");
  basic_block_t& bb1   = cfg->insert("bb1");
  basic_block_t& bb1_t = cfg->insert("bb1_t");
  basic_block_t& bb1_f = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& ret   = cfg->insert("ret");
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
  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(n2, 2);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, 123456, 1);
  // If we comment these two lines then we do only initialization of
  // even positions.
  //bb2.add(i1, i, n1);
  //bb2.array_store(a, i1, 123, 1);
  bb2.add(i, i, n2);
  ret.assign(tmp1, 6);
  ret.array_load(tmp2, a, tmp1, 1); // initialized
  return cfg;

}


// this is the program init_rand from Gange et.al paper.
cfg_t* prog9(variable_factory_t &vfac) 
{
  cfg_t* cfg = new cfg_t("entry","ret",ARR);
  basic_block_t& entry   = cfg->insert("entry");
  basic_block_t& bb1     = cfg->insert("bb1");
  basic_block_t& bb1_t   = cfg->insert("bb1_t");
  basic_block_t& bb1_f1  = cfg->insert("bb1_f1");
  basic_block_t& bb1_f2  = cfg->insert("bb1_f2");
  basic_block_t& bb1_f   = cfg->insert("bb1_f");
  basic_block_t& bb2   = cfg->insert("bb2");
  basic_block_t& bb2_a   = cfg->insert("bb2a");
  basic_block_t& bb2_b   = cfg->insert("bb2b");
  basic_block_t& ret     = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var i1(vfac["i1"]);
  z_var i2(vfac["i2"]);
  z_var n(vfac["n"]);
  z_var a(vfac["A"]);
  z_var tmp1(vfac["tmp1"]);
  z_var tmp2(vfac["tmp2"]);
  z_var nd(vfac["nd"]);
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f1; bb1 >> bb1_f2;
  bb1_f1 >> bb1_f;   bb1_f2 >> bb1_f; 
  bb1_t >> bb2; bb2 >> bb2_a; bb2 >> bb2_b; bb2_a >> bb1; bb2_b >> bb1; bb1_f >> ret;
  // assume array element of 1 byte
  entry.assume_array (a.name (), 0);
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i1, 0);
  entry.assign(i2, 0);  
  ///////
  // while (i1 < n && i2 < n){
  bb1_t.assume(i1 <= n -1);
  bb1_t.assume(i2 <= n -1);
  bb1_t.havoc (nd.name ());

  // if (*)
  bb2_a.assume (nd >= 1);
  bb2_a.array_store(a, i1, 1, 1);
  bb2_a.add(i1, i1, n1);
  // else
  bb2_b.assume (nd <= 0);
  bb2_b.array_store(a, i2, 2, 1);
  bb2_b.add(i2, i2, n1);
  // } end while
  bb1_f1.assume(i1 >= n);
  bb1_f2.assume(i2 >= n);
  ret.sub(tmp1, n, n1);
  ret.array_load(tmp2, a, tmp1, 1); // initialized
  return cfg;
}



void test1(){
  variable_factory_t vfac;
  cfg_t* cfg = prog1(vfac);
  crab::outs () << "Program 1: forall 0<= i< 10. a[i] = 123456\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test2(){
  variable_factory_t vfac;
  cfg_t* cfg = prog3(vfac);
  crab::outs () << "Program 2: forall 0<= i< 10. a[i] = b[i] = x and x = 123456\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test3(){
  variable_factory_t vfac;
  cfg_t* cfg = prog4(vfac);
  crab::outs () << "Program 3: forall 0<= i< 10. a[i] = 8 and b[i] = 5\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test4(){
  variable_factory_t vfac;
  cfg_t* cfg = prog5(vfac);
  crab::outs () << "Program 4: forall 0<= i < n. a[i] = 123456 (unbounded loop)\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test5(){
  variable_factory_t vfac;
  cfg_t* cfg = prog6(vfac);
  crab::outs () << "Program 5: for all 0<= i< 10. a[i] = 123456 (assume elem size of 4 bytes)\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test6(){
  variable_factory_t vfac;
  cfg_t* cfg = prog7(vfac);
  crab::outs () << "Program 6: a[0] = 89 and for all 1<= i < n. a[i] = a[i-1]\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test7(){
  variable_factory_t vfac;
  cfg_t* cfg = prog8(vfac);
  crab::outs () << "Program 7: forall 0<= i< 10 and i % 2 = 0. a[i] = 123456\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}


void test8(){
  variable_factory_t vfac;
  cfg_t* cfg = prog9(vfac);
  crab::outs () << "Program 8: forall 0<= i < n. 1 <= a[i] <= 2\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

void test9(){
  variable_factory_t vfac;
  cfg_t* cfg = prog2(vfac);
  crab::outs () << "Program 9: forall 0<= i < n. a[i] == 123456 (decrementing loop)\n";
  run<array_smashing_t> (cfg, vfac, false, 1, 2, 20);
  run<array_smashing<sdbm_domain_t> > (cfg, vfac, false, 1, 2, 20);
  delete cfg;
}

int main(int argc, char **argv) 
{
  SET_TEST_OPTIONS(argc,argv)

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


