#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t* prog1 (variable_factory_t &vfac, bool temp_add) 
{
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var i1(vfac["i1"]);
  varname_t a = vfac["A0"];
  z_var tmp3(vfac["tmp3"]);
  z_var val(vfac["val"]);
  varname_t tmp5 = vfac["tmp5"];
  varname_t tmp6 = vfac["tmp6"];

  z_cfg_t* cfg = new z_cfg_t("entry","ret", ARR);
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb1_t = cfg->insert ("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert ("bb1_f");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& ret   = cfg->insert ("ret");

  // assume array element of 1 byte

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  if (temp_add) 
  {
    bb2.add(i1, i, n1);
    bb2.assign(i, i1);
  }
  else
    bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, crab::ARR_INT_TYPE, tmp3, 1); // initialized
  ret.array_load(tmp6, a, crab::ARR_INT_TYPE, i, 1);    // top
  return cfg;
}


z_cfg_t* prog2(variable_factory_t &vfac, bool temp_sub) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var n0(vfac["n0"]);
  z_var n9(vfac["n9"]);
  z_var i(vfac["i"]);
  z_var i1(vfac["i1"]);
  z_var i2(vfac["i2"]);
  varname_t  a = vfac["A"];
  z_var tmp3(vfac["tmp3"]);
  varname_t tmp4 = vfac["tmp4"];
  varname_t tmp5 = vfac["tmp5"];
  z_var val(vfac["val"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array elements of 1 byte
  ///////
  entry.assign(i,9);
  bb1_t.assume(i >= 0);
  bb1_f.assume(i <= -1);
  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  if (temp_sub) 
  {
    bb2.sub(i1, i, 1);
    bb2.assign(i2, i1);
    bb2.assign(i, i2);

  } 
  else
    bb2.sub(i, i, 1);
  ret.assign(tmp3, 5);
  ret.array_load(tmp4, a, crab::ARR_INT_TYPE,tmp3, 1); // initialized
  ret.array_load(tmp5, a, crab::ARR_INT_TYPE,i, 1);    // top
  return cfg;
}

z_cfg_t* prog3(variable_factory_t &vfac) 
{
  z_cfg_t* cfg= new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;

  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  varname_t a = vfac["A"];
  varname_t b = vfac["B"];
  varname_t tmp1 = vfac["tmp1"];
  z_var tmp2(vfac["tmp2"]);
  z_var val(vfac["val"]);
  varname_t tmp3 = vfac["tmp3"];
  varname_t tmp4 = vfac["tmp4"];

  // assume array element of 1 byte
  entry.assign(n1, 1);
  entry.assign(i, 0);
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  bb2.array_load(tmp1, a, crab::ARR_INT_TYPE, i, 1);    
  bb2.array_store(b, crab::ARR_INT_TYPE, i, z_var(tmp1), 1);
  bb2.add(i, i, n1);
  ret.sub(tmp2, i, n1);
  ret.array_load(tmp3, b, crab::ARR_INT_TYPE, tmp2, 1); // initialized
  ret.array_load(tmp4, b, crab::ARR_INT_TYPE, i, 1);    // top
  return cfg;
}


z_cfg_t* prog4(variable_factory_t &vfac) 
{

  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  varname_t a = vfac["A"];
  varname_t b = vfac["B"];
  z_var tmp3(vfac["tmp3"]);
  varname_t tmp5 = vfac["tmp5"];
  varname_t tmp6 = vfac["tmp6"];
  z_var val1(vfac["val1"]);
  z_var val2(vfac["val2"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);

  bb2.assign(val1, 8); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val1, 1);
  bb2.assign(val2, 5); 
  bb2.array_store(b, crab::ARR_INT_TYPE, i, val2, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, crab::ARR_INT_TYPE, tmp3, 1); 
  ret.array_load(tmp6, b, crab::ARR_INT_TYPE, tmp3, 1); 
  return cfg;
}

z_cfg_t* prog5(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  varname_t a= vfac["A"];
  z_var tmp1(vfac["tmp1"]);
  z_var val(vfac["val"]);
  varname_t tmp2 = vfac["tmp2"];

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  bb2.add(i, i, n1);
  ret.sub(tmp1, i, n1);
  ret.array_load(tmp2, a, crab::ARR_INT_TYPE, tmp1, 1); // initialized
  return cfg;
}

z_cfg_t* prog6(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var i(vfac["i"]);
  varname_t a = vfac["A"];
  z_var tmp(vfac["tmp"]);
  z_var offset(vfac["o"]);
  z_var tmp2(vfac["tmp2"]);
  z_var val(vfac["val"]);
  varname_t tmp4 = vfac["tmp4"];
  z_var tmp5(vfac["tmp5"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 4 bytes
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(10 <= i);
  bb2.assign(tmp, i);
  bb2.mul(offset, tmp, 4); 

  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, offset, val, 4);
  bb2.add(i, i, 1);
  ret.assign (tmp5, 8);
  ret.array_load(tmp4, a, crab::ARR_INT_TYPE, tmp5, 4);    
  return cfg;
}

z_cfg_t* prog7(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  varname_t a = vfac["A"];
  z_var tmp1(vfac["tmp1"]);
  varname_t tmp2 = vfac["tmp2"];
  z_var tmp3(vfac["tmp3"]);
  z_var val(vfac["val"]);
  varname_t tmp4 = vfac["tmp4"];
  z_var x(vfac["x"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assume(n >= 2);
  entry.assign(n1, 1);
  entry.assign(i , 0);

  entry.assign(val, 89); 
  entry.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  entry.assign(i , 1);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  ///////
  bb2.sub(tmp1, i, n1);
  bb2.array_load(tmp2, a, crab::ARR_INT_TYPE, tmp1, 1); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, z_var(tmp2), 1);
  bb2.add(i, i, n1);
  ///////
  ret.sub(tmp3, n, n1);
  ret.array_load(tmp4, a, crab::ARR_INT_TYPE, tmp3, 1); 
  return cfg;
}

// Initialize only even positions
z_cfg_t* prog8(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret", ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var n1(vfac["n1"]);
  z_var n2(vfac["n2"]);
  z_var i(vfac["i"]);
  z_var i1(vfac["i1"]);
  z_var n(vfac["n"]);
  varname_t a = vfac["A"];
  z_var tmp1(vfac["tmp1"]);
  z_var val(vfac["val"]);
  varname_t tmp2 = vfac["tmp2"];
  z_var tmp3(vfac["tmp3"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(n2, 2);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);

  bb2.assign(val, 123456); 
  bb2.array_store(a, crab::ARR_INT_TYPE, i, val, 1);
  // If we comment these two lines then we do only initialization of
  // even positions.
  //bb2.add(i1, i, n1);
  // bb2.assign(val, 123); 
  // bb2.array_store(a, crab::ARR_INT_TYPE, i1, val, 1);
  bb2.add(i, i, n2);
  ret.assign(tmp1, 6);
  ret.array_load(tmp2, a, crab::ARR_INT_TYPE, tmp1, 1); // initialized
  return cfg;

}


// this is the program init_rand from Gange et.al paper.
z_cfg_t* prog9(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry   = cfg->insert("entry");
  z_basic_block_t& bb1     = cfg->insert("bb1");
  z_basic_block_t& bb1_t   = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f1  = cfg->insert("bb1_f1");
  z_basic_block_t& bb1_f2  = cfg->insert("bb1_f2");
  z_basic_block_t& bb1_f   = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& bb2_a   = cfg->insert("bb2a");
  z_basic_block_t& bb2_b   = cfg->insert("bb2b");
  z_basic_block_t& bb3   = cfg->insert("bb3");
  z_basic_block_t& ret     = cfg->insert("ret");
  z_var i1(vfac["i1"]);
  z_var i2(vfac["i2"]);
  z_var n(vfac["n"]);
  varname_t a = vfac["A"];
  z_var tmp1(vfac["tmp1"]);
  z_var val(vfac["val"]);
  varname_t tmp2 = vfac["tmp2"];
  z_var nd(vfac["nd"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f1; bb1 >> bb1_f2;
  bb1_f1 >> bb1_f;   bb1_f2 >> bb1_f; 
  bb1_t >> bb2; bb2 >> bb2_a; bb2 >> bb2_b; bb2_a >> bb3;  bb2_b >> bb3; bb3 >> bb1; bb1_f >> ret;
  ////////
  // assume array element of 1 byte
  entry.assume(n >= 1);
  entry.assign(i1, 0);
  entry.assign(i2, 0);  
  ///////
  // while (i1 < n && i2 < n){
  bb1_t.assume(i1 <= n -1);
  bb1_t.assume(i2 <= n -1);
  bb1_t.havoc (nd.name ());

  // if (*)
  bb2_a.assume (nd >= 1);
  bb2_a.assign(val, 1); 
  bb2_a.array_store(a, crab::ARR_INT_TYPE, i1, val, 1);
  bb2_a.add(i1, i1, 1);
  // else
  bb2_b.assume (nd <= 0);
  bb2_b.assign(val, 2); 
  bb2_b.array_store(a, crab::ARR_INT_TYPE, i2, val, 1);
  bb2_b.add(i2, i2, 1);
  // } end while
  bb1_f1.assume(i1 >= n);
  bb1_f2.assume(i2 >= n);
  ret.sub(tmp1, n, 1);
  ret.array_load(tmp2, a, crab::ARR_INT_TYPE, tmp1, 1); // initialized
  return cfg;
}

z_cfg_t* prog10(variable_factory_t &vfac) 
{
  z_cfg_t* cfg = new z_cfg_t("entry","ret",ARR);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& bb1   = cfg->insert("bb1");
  z_basic_block_t& bb1_t = cfg->insert("bb1_t");
  z_basic_block_t& bb1_f = cfg->insert("bb1_f");
  z_basic_block_t& bb2   = cfg->insert("bb2");
  z_basic_block_t& ret   = cfg->insert("ret");
  z_var i(vfac["i"]);
  z_var n(vfac["n"]);
  z_var max(vfac["max"]);
  varname_t a = vfac["A"];
  varname_t b = vfac["B"];
  varname_t obj1 = vfac["obj1"];
  varname_t tmp1 = vfac["tmp1"];
  varname_t tmp2 = vfac["tmp2"];
  z_var tmp3(vfac["tmp3"]);

  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  
  uint64_t elem_size = 1;
  entry.assume(n >= 1);
  // forall i :: is_not_null(a[i])
  // forall i :: is_not_null(b[i])
  entry.ptr_new_object (obj1, 0);
  entry.sub (max, n, 1); // max = n-1 
  entry.array_assume (a, crab::ARR_PTR_TYPE, elem_size, 0, max, z_var(obj1));
  ///
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  // b[i] := a[i]
  bb2.array_load(tmp1, a, crab::ARR_PTR_TYPE, i, elem_size);
  bb2.array_store(b, crab::ARR_PTR_TYPE, i, z_var(tmp1), elem_size);
  bb2.add(i, i, 1);
  ret.sub(tmp3, i, 1);
  // read b[i-1] 
  ret.array_load(tmp2, b, crab::ARR_PTR_TYPE, tmp3, elem_size); 
  return cfg;
}


// template <typename ArrayDomain>
// void run(z_cfg_ref_t cfg, string name)
// {
//   crab::outs() << "--- " << name  << "\n";
//   cfg.simplify ();
//   crab::outs() << cfg << "\n";
  
//   ArrayDomain inv = ArrayDomain::top ();
//   using namespace crab::analyzer;
//   intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t,ArrayDomain> It (cfg, inv, nullptr, 1, 2, 20);
//   It.run ();
//   crab::outs() << "Results with " << ArrayDomain::getDomainName () << ":\n";
//   crab::outs() << "(Invariants hold at the block's entries)\n";
//   for (auto &b : cfg)
//   {
//     // invariants at the entry of the block
//     auto inv = It.get_pre(b.label ());
//     crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
//   }
//   crab::outs() << "\n";
//   if (stats_enabled) {
//     crab::CrabStats::Print(crab::outs());
//     crab::CrabStats::reset();
//   }  
// }

typedef array_sparse_graph_domain<z_sdbm_domain_t, z_interval_domain_t> array_sgraph_domain_t;

void test1(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog1(vfac, false);
  crab::outs () << "Program 1: forall 0<= i< 10. a[i] = 123456";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t> (*cfg, "Program 1: forall 0<= i< 10. a[i] = 123456");
  delete cfg;
}


void test2(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog3(vfac);
  crab::outs () << "Program 2: forall 0<= i< 10. a[i] = b[i] = x and x = 123456";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 2: forall 0<= i< 10. a[i] = b[i] = x and x = 123456");
  delete cfg;
}

void test3(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog4(vfac);
  crab::outs () << "Program 3: forall 0<= i< 10. a[i] = 8 and b[i] = 5";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //  run<array_sgraph_domain_t>(*cfg, "Program 3: forall 0<= i< 10. a[i] = 8 and b[i] = 5");
  delete cfg;
}

void test4(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog5(vfac);
  crab::outs () << "Program 4: forall 0<= i < n. a[i] = 123456 (unbounded loop)";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 4: forall 0<= i < n. a[i] = 123456 (unbounded loop)");
  delete cfg;
}

void test5(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog6(vfac);
  crab::outs () << "Program 5: for all 0<= i< 10. a[i] = 123456 (assume elem size of 4 bytes)";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 5: for all 0<= i< 10. a[i] = 123456 (assume elem size of 4 bytes)");
  delete cfg;
}

void test6(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog7(vfac);
  crab::outs () << "Program 6: a[0] = 89 and for all 1<= i < n. a[i] = a[i-1]";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 6: a[0] = 89 and for all 1<= i < n. a[i] = a[i-1]");
  delete cfg;
}

void test7(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog8(vfac);
  crab::outs () << "Program 7: forall 0<= i< 10 and i % 2 = 0. a[i] = 123456";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 7: forall 0<= i< 10 and i % 2 = 0. a[i] = 123456");
  delete cfg;
}


void test8(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog9(vfac);
  crab::outs () << "Program 8: forall 0<= i < n. 1 <= a[i] <= 2";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 8: forall 0<= i < n. 1 <= a[i] <= 2");
  delete cfg;
}

void test9(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog2(vfac, false);
  crab::outs () << "Program 9: forall 0<= i < n. a[i] == 123456 (decrementing loop)";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t>(*cfg, "Program 9: forall 0<= i < n. a[i] == 123456 (decrementing loop)");
  delete cfg;
}

void test10(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog2(vfac, true);
  crab::outs () << "Program 10: forall 0<= i < n. a[i] == 123456 (decrementing loop w/ temp vars)";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);
  //run<array_sgraph_domain_t>(*cfg, "Program 10: forall 0<= i < n. a[i] == 123456 (decrementing loop w/ temp vars)");
  delete cfg;
}

void test11(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog1(vfac, true);
  crab::outs () << "Program 11: forall 0<= i< 10. a[i] = 123456 (w/ temp vars)";
  run<array_sgraph_domain_t> (cfg, false, 1, 2, 20, stats_enabled);  
  //run<array_sgraph_domain_t> (*cfg, "Program 11: forall 0<= i< 10. a[i] = 123456 (w/ temp vars)");
  delete cfg;
}

void test12(){
  variable_factory_t vfac;
  z_cfg_t* cfg = prog10(vfac); 
  crab::outs () << "Program 12: forall 0<= i < n. is_not_null(a[i]) &&  is_not_null(b[i]) \n";
  run<array_sparse_graph_domain<z_num_null_domain_t, z_nullity_domain_t> >
    (cfg, false, 1, 2, 20, stats_enabled);
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
  test10 ();
  test11 ();
  test12 ();
  return 0;
}


