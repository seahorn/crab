#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/inter/top_down_inter_params.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

/*
 * This test exercises the case where a function *re-assigns* its input
 * parameter in the body. This used to require a CFG pre-pass that
 * renamed the input parameter so that it was used exactly once (as the
 * rhs of an assignment at the entry block). The top-down analyzer now
 * snapshots the entry value of each input parameter internally, so the
 * summary keeps relating the input's *entry* value with the outputs.
 *
 * g(p) {
 *   p = p + 1;   // input parameter re-assigned
 *   o = p;
 *   return o;
 * }
 *
 * main() {
 *   a = 10;
 *   r = g(a);
 *   assert(r == 11);   // output uses the mutated value
 *   assert(a == 10);   // the caller's actual argument is NOT clobbered
 * }
 */

z_cfg_t *g(variable_factory_t &vfac) {
  z_var p(vfac["p"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("g", {p}, {o});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  // p := p + 1  (re-assign the input parameter)
  entry.add(p, p, 1);
  exit.assign(o, p);
  return cfg;
}

z_cfg_t *m(variable_factory_t &vfac) {
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var r(vfac["r"], crab::INT_TYPE, 32);
  z_var res(vfac["res"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("main", {}, {res});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  entry.assign(a, 10);
  entry.callsite("g", {r}, {a});
  exit.assertion(r == 11); // output reflects the mutation inside g
  exit.assertion(a == 10); // caller's actual argument stays unchanged
  exit.assign(res, r);
  return cfg;
}

using callgraph_t = call_graph<z_cfg_ref_t>;
using inter_params_t = top_down_inter_analyzer_parameters<callgraph_t>;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *t1 = g(vfac);
  z_cfg_t *t2 = m(vfac);

  crab::outs() << *t1 << "\n" << *t2 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  callgraph_t cg(cfgs);
  {
    z_dbm_domain_t init;
    crab::outs() << "Running top-down inter-procedural analysis with "
                 << init.domain_name() << "\n";
    // it should prove all assertions
    inter_params_t params;
    td_inter_run(cg, init, params, true, false, false);
  }

  delete t1;
  delete t2;
  return 0;
}
