#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog1(variable_factory_t &vfac) {
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A0"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A0_prop"], crab::INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp5(vfac["tmp5"], crab::INT_TYPE, 32);
  z_var tmp6(vfac["tmp6"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");

  uint64_t elem_size = 1;

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  ////////
  entry.assign(a_p, 0);
  entry.array_init(a, 0, 9, a_p, elem_size);
  /////////
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val, 123456);
  bb2.array_store(a, i, val, elem_size);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3, elem_size); // initialized
  ret.array_load(tmp6, a, i, elem_size);    // top
  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n0(vfac["n0"], crab::INT_TYPE, 32);
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var n9(vfac["n9"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp4(vfac["tmp4"], crab::INT_TYPE, 32);
  z_var tmp5(vfac["tmp5"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  uint64_t elem_size = 1;

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  ////////
  entry.assign(a_p, 0);
  entry.array_init(a, 0, 9, a_p, elem_size);
  /////////
  entry.assign(n0, 0); // we need it to be considered as graph node
  entry.assign(n1, 1);
  entry.assign(n9, 9); // we need it to be considered as graph node
  entry.assign(i, n9);
  ///////
  bb1_t.assume(i >= 0);
  bb1_f.assume(i <= -1);
  bb2.assign(val, 123456);
  bb2.array_store(a, i, val, elem_size);
  bb2.sub(i, i, n1);
  ret.assign(tmp3, 5);
  ret.array_load(tmp4, a, tmp3, elem_size); // initialized
  ret.array_load(tmp5, a, i, elem_size);    // top
  return cfg;
}

z_cfg_t *prog3(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var b(vfac["B"], crab::ARR_INT_TYPE, 32);
  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp4(vfac["tmp4"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  uint64_t elem_size = 1;

  entry.assign(a_p, 0);
  entry.array_init(a, 0, 9, a_p, elem_size);
  entry.array_init(b, 0, 9, a_p, elem_size);

  entry.assign(n1, 1);
  entry.assign(i, 0);
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val, 123456);
  bb2.array_store(a, i, val, elem_size);
  bb2.array_load(tmp1, a, i, elem_size);
  bb2.array_store(b, i, tmp1, elem_size);
  bb2.add(i, i, n1);
  ret.sub(tmp2, i, n1);
  ret.array_load(tmp3, b, tmp2, elem_size); // initialized
  ret.array_load(tmp4, b, i, elem_size);    // top
  return cfg;
}

z_cfg_t *prog4(variable_factory_t &vfac) {

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var b(vfac["B"], crab::ARR_INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp5(vfac["tmp5"], crab::INT_TYPE, 32);
  z_var tmp6(vfac["tmp6"], crab::INT_TYPE, 32);
  z_var val1(vfac["val1"], crab::INT_TYPE, 32);
  z_var val2(vfac["val2"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;

  entry.assign(a_p, 0);
  entry.array_init(a, 0, 9, a_p, elem_size);
  entry.array_init(b, 0, 9, a_p, elem_size);

  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val1, 8);
  bb2.assign(val2, 5);
  bb2.array_store(a, i, val1, elem_size);
  bb2.array_store(b, i, val2, elem_size);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3, elem_size);
  ret.array_load(tmp6, b, tmp3, elem_size);
  return cfg;
}

z_cfg_t *prog4b(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_BOOL_TYPE, 1);
  z_var b(vfac["B"], crab::ARR_BOOL_TYPE, 1);
  z_var tt(vfac["TRUE"], crab::BOOL_TYPE, 1);
  z_var ff(vfac["FALSE"], crab::BOOL_TYPE, 1);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp5(vfac["tmp5"], crab::BOOL_TYPE, 1);
  z_var tmp6(vfac["tmp6"], crab::BOOL_TYPE, 1);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;
  entry.bool_assign(tt, z_lin_cst_t::get_true());
  entry.bool_assign(ff, z_lin_cst_t::get_false());
  entry.array_init(a, 0, 9, tt, elem_size);
  entry.array_init(b, 0, 9, ff, elem_size);

  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.array_store(a, i, tt, elem_size);
  bb2.array_store(b, i, ff, elem_size);
  bb2.add(i, i, n1);
  ret.sub(tmp3, i, n1);
  ret.array_load(tmp5, a, tmp3, elem_size);
  ret.array_load(tmp6, b, tmp3, elem_size);
  return cfg;
}

z_cfg_t *prog5(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;
  entry.assign(a_p, 0);
  entry.array_init(a, 0, n, a_p, elem_size);

  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  bb2.assign(val, 123456);
  bb2.array_store(a, i, val, elem_size);
  bb2.add(i, i, n1);
  ret.sub(tmp1, i, n1);
  ret.array_load(tmp2, a, tmp1, elem_size); // initialized
  return cfg;
}

z_cfg_t *prog6(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  z_var offset(vfac["o"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var tmp4(vfac["tmp4"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 4;
  entry.assign(a_p, 0);
  entry.array_init(a, 0, 9, a_p, elem_size);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(10 <= i);
  bb2.assign(tmp, i);
  bb2.mul(offset, tmp, 4);
  bb2.assign(val, 123456);
  bb2.array_store(a, offset, val, elem_size);
  bb2.add(i, i, 1);
  ret.array_load(tmp4, a, 8, elem_size);
  return cfg;
}

z_cfg_t *prog7(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var tmp4(vfac["tmp4"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;
  // assume (forall i. a[i] =0);
  entry.assign(a_p, 0);
  entry.array_init(a, 0, n, a_p, elem_size);
  //////
  entry.assume(n >= 2);
  entry.assign(n1, 1);
  entry.assign(i, 0);
  entry.assign(val, 89);
  entry.array_store(a, i, val, elem_size);
  entry.assign(i, 1);
  ///////
  bb1_t.assume(i <= n - 1);
  bb1_f.assume(i >= n);
  ///////
  bb2.sub(tmp1, i, n1);
  bb2.array_load(tmp2, a, tmp1, elem_size);
  bb2.array_store(a, i, tmp2, elem_size);
  bb2.add(i, i, n1);
  ///////
  ret.sub(tmp3, n, n1);
  ret.array_load(tmp4, a, tmp3, elem_size);
  return cfg;
}

// Initialize only even positions
z_cfg_t *prog8(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var n2(vfac["n2"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var tmp3(vfac["tmp3"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;
  entry.assign(a_p, 0);
  entry.array_init(a, 0, 10, a_p, elem_size);
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(n2, 2);
  entry.assign(i, 0);
  ///////
  bb1_t.assume(i <= 9);
  bb1_f.assume(i >= 10);
  bb2.assign(val, 123456);
  bb2.array_store(a, i, val, elem_size);
  // If we comment these two lines then we do only initialization of
  // even positions.
  // bb2.add(i1, i, n1);
  // bb2.assign(val, 123);
  // bb2.array_store(a,  i1, val, elem_size);
  bb2.add(i, i, n2);
  ret.assign(tmp1, 6);
  ret.array_load(tmp2, a, tmp1, elem_size); // initialized
  return cfg;
}

// this is the program init_rand from Gange et.al paper.
z_cfg_t *prog9(variable_factory_t &vfac) {
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f1 = cfg->insert("bb1_f1");
  z_basic_block_t &bb1_f2 = cfg->insert("bb1_f2");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb2_a = cfg->insert("bb2a");
  z_basic_block_t &bb2_b = cfg->insert("bb2b");
  z_basic_block_t &ret = cfg->insert("ret");
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);
  z_var a_p(vfac["A_p"], crab::INT_TYPE, 32);
  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var val(vfac["val"], crab::INT_TYPE, 32);

  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f1;
  bb1 >> bb1_f2;
  bb1_f1 >> bb1_f;
  bb1_f2 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb2_a;
  bb2 >> bb2_b;
  bb2_a >> bb1;
  bb2_b >> bb1;
  bb1_f >> ret;

  uint64_t elem_size = 1;
  entry.assign(a_p, 0);
  entry.array_init(a, 0, n, a_p, elem_size);
  entry.assume(n >= 1);
  entry.assign(n1, 1);
  entry.assign(i1, 0);
  entry.assign(i2, 0);
  ///////
  // while (i1 < n && i2 < n){
  bb1_t.assume(i1 <= n - 1);
  bb1_t.assume(i2 <= n - 1);
  bb1_t.havoc(nd);

  // if (*)
  bb2_a.assume(nd >= 1);
  bb2_a.assign(val, 1);
  bb2_a.array_store(a, i1, val, elem_size);
  bb2_a.add(i1, i1, n1);
  // else
  bb2_b.assume(nd <= 0);
  bb2_b.assign(val, 2);
  bb2_b.array_store(a, i2, val, elem_size);
  bb2_b.add(i2, i2, n1);
  // } end while
  bb1_f1.assume(i1 >= n);
  bb1_f2.assume(i2 >= n);
  ret.sub(tmp1, n, n1);
  ret.array_load(tmp2, a, tmp1, elem_size); // initialized
  return cfg;
}

void test1(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog1(vfac);
  crab::outs() << "Program 1: forall 0<= i< 10. a[i] = 123456\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test2(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog3(vfac);
  crab::outs()
      << "Program 2: forall 0<= i< 10. a[i] = b[i] = x and x = 123456\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test3(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog4(vfac);
  crab::outs() << "Program 3: forall 0<= i< 10. a[i] = 8 and b[i] = 5\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test4(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog5(vfac);
  crab::outs()
      << "Program 4: forall 0<= i < n. a[i] = 123456 (unbounded loop)\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test5(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog6(vfac);
  crab::outs() << "Program 5: for all 0<= i< 10. a[i] = 123456 (assume elem "
                  "size of 4 bytes)\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test6(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog7(vfac);
  crab::outs() << "Program 6: a[0] = 89 and for all 1<= i < n. a[i] = a[i-1]\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test7(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog8(vfac);
  crab::outs() << "Program 7: forall 0<= i< 10 and i % 2 = 0. a[i] = 123456\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test8(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog9(vfac);
  crab::outs() << "Program 8: forall 0<= i < n. 1 <= a[i] <= 2\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test9(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog2(vfac);
  crab::outs()
      << "Program 9: forall 0<= i < n. a[i] == 123456 (decrementing loop)\n";
  {
    array_smashing<z_dis_interval_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    array_smashing<z_sdbm_domain_t> init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  delete cfg;
}

void test10(bool stats_enabled) {
  variable_factory_t vfac;
  z_cfg_t *cfg = prog4b(vfac);
  crab::outs()
      << "Program 10: forall 0<= i< 10. a[i] = true and b[i] = false\n";
  array_smashing<z_bool_num_domain_t> init;
  run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  delete cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  test1(stats_enabled);
  test2(stats_enabled);
  test3(stats_enabled);
  test4(stats_enabled);
  test5(stats_enabled);
  test6(stats_enabled);
  test7(stats_enabled);
  test8(stats_enabled);
  test9(stats_enabled);
  test10(stats_enabled);

  return 0;
}
