#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter/inter_params.hpp>
#include <crab/cg/cg_bgl.hpp>


using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;
/*
  int foo(n) {
    if (n <= 0) return n;
    else return 1 + foo(n-1);
  }

  void main() {
    res = foo(10);
    assert(res == 10);
  }
 */
z_cfg_t *foo(variable_factory_t &vfac) {
  // Defining program variables
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var n1(vfac["n1"], crab::INT_TYPE, 32);  
  z_var ret_val(vfac["res"], crab::INT_TYPE, 32);
  z_var foo_ret(vfac["foo_ret"], crab::INT_TYPE, 32);
  
  function_decl<z_number, varname_t> decl("foo", {n}, {ret_val});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &base = cfg->insert("base");
  z_basic_block_t &rec = cfg->insert("rec");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> base;
  entry >> rec;
  base >> exit;
  rec >> exit;
  // adding statements
  base.assume(n <= 0);
  base.assign(ret_val, n);
  rec.assume(n >= 1);
  rec.sub(n1, n, 1);
  rec.callsite("foo", {foo_ret}, {n1});
  rec.add(ret_val, foo_ret, 1); 
  return cfg;
}


z_cfg_t *_main(variable_factory_t &vfac) {
  // Defining program variables
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var res(vfac["res"], crab::INT_TYPE, 32);
  
  function_decl<z_number, varname_t> decl("main", {}, {});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(n, 10);
  entry.callsite("foo", {res}, {n});
  exit.assertion(res == 10);
  return cfg;
}

using callgraph_t = call_graph<z_cfg_ref_t>;
using inter_params_t = inter_analyzer_parameters<callgraph_t>;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  
  variable_factory_t vfac;
  z_cfg_t *t1 = foo(vfac);  
  crab::outs() << *t1 << "\n";
  z_cfg_t *t2 = _main(vfac);  
  crab::outs() << *t2 << "\n";
  
  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  callgraph_t cg(cfgs);
  crab::outs() << "CallGraph=" << cg << "\n";
  {
    z_sdbm_domain_t init;
    crab::outs() << "Running top-down inter-procedural analysis with "
                 << init.domain_name() << "\n";
    inter_params_t params;
    params.analyze_recursive_functions = true;    
    td_inter_run(cg, init, params, true, true, false);
  }
  
  
  delete t1;

  return 0;
}
