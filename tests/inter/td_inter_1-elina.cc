#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter/top_down_inter_params.hpp>
#include <crab/cg/cg_bgl.hpp>


using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;
/*
foo(x) {
  y = x+1;
  z = y+2;
  return z;
}

rec1(s) {
  r = s - 1;
  t = rec2(r)
  return t;
}

rec2(s) {
  a = 10;
  r = s - 1;
  t = rec1(r)
  assert(a >= 5);
  return t;
}

bar(a) {
   x = a;
   w = 5;
   y = foo(x);
   assert(y>=6);
   assert(y<=17);
   return y;
}

main(){
   x=3;
   x3=4;
   x4=5;
   x5=6;
   y = bar(x);
   assert(y == 6);
   u = rec1(y);
   z = y +2;
   y3 = bar(x);
   assert(y3 == 6);
   z3 = y3 +z;
   w = foo(z3)
   assert(w == 17);
   y4 = bar(x3);
   assert(y4 == 7);
   y5 = bar(x4);
   assert(y5 == 8);
   y6 = bar(x5);
   assert(y6 == 9);
   res = w + y4 + y5 + y6 ;
   assert(res ==41);
}
 */

z_cfg_t *foo(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("foo", {x}, {z});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add(y, x, 1);
  exit.add(z, y, 2);
  return cfg;
}

z_cfg_t *rec1(variable_factory_t &vfac) {
  // Defining program variables
  z_var r(vfac["r"], crab::INT_TYPE, 32);
  z_var s(vfac["s"], crab::INT_TYPE, 32);
  z_var t(vfac["t"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("rec1", {s}, {t});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.sub(r, s, 1);
  exit.callsite("rec2", {t}, {r});
  return cfg;
}

z_cfg_t *rec2(variable_factory_t &vfac) {
  // Defining program variables
  z_var r(vfac["r1"], crab::INT_TYPE, 32);
  z_var s(vfac["s1"], crab::INT_TYPE, 32);
  z_var t(vfac["t1"], crab::INT_TYPE, 32);
  z_var a(vfac["a"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("rec2", {s}, {t});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(a, 10);
  entry.sub(r, s, 1);
  exit.callsite("rec1", {t}, {r});
  // exit.callsite("foo", {t}, {t});
  exit.assertion(a >= 5);
  return cfg;
}

z_cfg_t *bar(variable_factory_t &vfac) {
  // Defining program variables
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var x(vfac["x1"], crab::INT_TYPE, 32);
  z_var y(vfac["y1"], crab::INT_TYPE, 32);
  z_var w(vfac["w1"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("bar", {a}, {y});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  exit.callsite("foo", {y}, {x});
  exit.assertion(y >= 6);
  exit.assertion(y <= 17);
  entry.assign(x, a);
  entry.assign(w, 5);
  return cfg;
}

z_cfg_t *m(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x2"], crab::INT_TYPE, 32);
  z_var x3(vfac["x3"], crab::INT_TYPE, 32);
  z_var x4(vfac["x4"], crab::INT_TYPE, 32);
  z_var x5(vfac["x5"], crab::INT_TYPE, 32);
  z_var y(vfac["y2"], crab::INT_TYPE, 32);
  z_var y3(vfac["y3"], crab::INT_TYPE, 32);
  z_var y4(vfac["y4"], crab::INT_TYPE, 32);
  z_var y5(vfac["y5"], crab::INT_TYPE, 32);
  z_var y6(vfac["y6"], crab::INT_TYPE, 32);
  z_var z(vfac["z2"], crab::INT_TYPE, 32);
  z_var z3(vfac["z3"], crab::INT_TYPE, 32);
  z_var u(vfac["__"], crab::INT_TYPE, 32);
  z_var w(vfac["w2"], crab::INT_TYPE, 32);
  z_var res(vfac["res"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("main", {}, {res});

  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  entry.assign(x3, 4);
  entry.assign(x4, 5);
  entry.assign(x5, 6);
  entry.callsite("bar", {y}, {x});
  entry.assertion(y == 6);
  /////
  entry.callsite("rec1", {u}, {y});
  /////
  exit.add(z, y, 2);
  exit.callsite("bar", {y3}, {x});
  exit.assertion(y3 == 6);
  exit.add(z3, y3, z);
  exit.callsite("foo", {w}, {z3});
  exit.assertion(w == 17); // provable even if we don't join calling contexts
  /// At this callsite, z3 = 14
  /// We have the summary (after joining calling contexts): 
  ///	I={x -> [3, 14]}
  ///   O={x -> [3, 14], z -> [6, 17], z-x<=3, x-z<=-3}
  /// Without the difference constraints w should be [6,17]. However,
  /// with the difference constraints, w is 17 (i.e., no lose of
  /// precision).
  /// 
  exit.callsite("bar", {y4}, {x3});
  exit.assertion(y4 == 7);
  exit.callsite("bar", {y5}, {x4});
  exit.assertion(y5 == 8);
  exit.callsite("bar", {y6}, {x5});
  exit.assertion(y6 == 9);
  exit.add(res, w, y4);
  exit.add(res, res, y5);
  exit.add(res, res, y6);
  exit.assertion(res == 41);
  return cfg;
}

using callgraph_t = call_graph<z_cfg_ref_t>;
using inter_params_t = top_down_inter_analyzer_parameters<callgraph_t>;

int main(int argc, char **argv) {
#ifdef HAVE_ELINA
  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *t1 = foo(vfac);
  z_cfg_t *t2 = bar(vfac);
  z_cfg_t *t3 = rec1(vfac);
  z_cfg_t *t4 = rec2(vfac);
  z_cfg_t *t5 = m(vfac);

  crab::outs() << *t1 << "\n"
               << *t2 << "\n"
               << *t3 << "\n"
               << *t4 << "\n"
               << *t5 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2, *t3, *t4, *t5});
  callgraph_t cg(cfgs);
  z_oct_elina_domain_t init;
  crab::outs() << "Running top-down inter-procedural analysis with "
	       << init.domain_name() << "\n";
  
  /////////////////////////////////////////
  // it should prove all assertions
  /////////////////////////////////////////
  inter_params_t params1;
  td_inter_run(cg, init, params1, true, false, false);
  /////////////////////////////////////////
  // it should prove all assertions (see above comments)    
  /////////////////////////////////////////
  inter_params_t params2;
  params2.max_call_contexts = 3;
  params2.checker_verbosity = 1;
  td_inter_run(cg, init, params2, true, false, false);

  delete t1;
  delete t2;
  delete t3;
  delete t4;
  delete t5;
#endif

  return 0;
}
