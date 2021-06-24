#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/transforms/dce.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog1(variable_factory_t &vfac, crab::tag_manager &as_man) {
  /*
    k := 2147483648
    o := 4;
    p := malloc(...);
    ..
    i := 0;
    while (i <= 99) {
       nd := *;
       q := p + o;
       k := 5;
       inc := (nd ? 1: 2);
       i := i + inc;
    }
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);
  z_var mem(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var_or_cst_t size400(z_number(400), crab::variable_type(crab::INT_TYPE, 32));
  
  // entry and exit block
  auto cfg = new z_cfg_t("x0", "ret");
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &x1 = cfg->insert("x1");
  z_basic_block_t &x2 = cfg->insert("x2");
  z_basic_block_t &x3 = cfg->insert("x3");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> x1;
  x1 >> x2;
  x2 >> x3;
  x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  // adding statements
  x0.assign(k, 2147483648);
  x0.assign(o, 4);
  x0.make_ref(p, mem, size400, as_man.mk_tag());
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.havoc(nd);
  bb2.gep_ref(q, mem, p, mem, o);
  bb2.assign(k, 5);
  bb2.select(inc, nd, 1, 2);
  bb2.add(i, i, inc);
  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac, crab::tag_manager &as_man) {
  /*
    k := 2147483648
    o := 4;
    p := malloc(...);
    ..
    i := 0;
    while (i <= 99) {
       nd := *;
       q := p + o;
       k := 5;
       inc := (nd ? 1: 2);
       i := i + inc;
    }
    k := *q;
    return k;
   */

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);
  z_var mem(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var_or_cst_t size400(z_number(400), crab::variable_type(crab::INT_TYPE, 32));  
  // entry and exit block
  typename z_cfg_t::fdecl_t fdecl("main", {}, {k});  
  auto cfg = new z_cfg_t("x0", "ret", fdecl);
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &x1 = cfg->insert("x1");
  z_basic_block_t &x2 = cfg->insert("x2");
  z_basic_block_t &x3 = cfg->insert("x3");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> x1;
  x1 >> x2;
  x2 >> x3;
  x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  // adding statements
  x0.assign(k, 2147483648);
  x0.assign(o, 4);
  x0.make_ref(p, mem, size400, as_man.mk_tag());
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.havoc(nd);
  bb2.gep_ref(q, mem, p, mem, o);
  bb2.assign(k, 5);
  bb2.select(inc, nd, 1, 2);
  bb2.add(i, i, inc);
  ////////////////
  ret.load_from_ref(k, q, mem);
  ////////////////
  return cfg;
}

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  crab::tag_manager as_man;
  using z_cfg_ref_t = crab::cfg::cfg_ref<z_cfg_t>;
  using dce_t = crab::transforms::dead_code_elimination<z_cfg_ref_t>;

  {
    z_cfg_t *cfg = prog1(vfac, as_man);
    crab::outs() << "CFG\n" << *cfg << "\n";
    z_cfg_t *cloned_cfg = cfg->clone();
    dce_t dce;
    z_cfg_ref_t cfg_ref(*cfg);
    dce.run(cfg_ref);
    crab::outs() << "After " << dce.get_name() << "\n" << *cfg << "\n";
    crab::outs() << "Cloned CFG\n" << *cloned_cfg << "\n";
    delete cfg;
    delete cloned_cfg;
  }

  {
    z_cfg_t *cfg = prog2(vfac, as_man);
    crab::outs() << "CFG\n" << *cfg << "\n";
    z_cfg_t *cloned_cfg = cfg->clone();
    dce_t dce;
    z_cfg_ref_t cfg_ref(*cfg);
    dce.run(cfg_ref);
    crab::outs() << "After " << dce.get_name() << "\n" << *cfg << "\n";
    crab::outs() << "Cloned CFG\n" << *cloned_cfg << "\n";
    delete cfg;
    delete cloned_cfg;
  }

  return 0;
}
