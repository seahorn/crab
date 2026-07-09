#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/inter/inter_params.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

/*
 * Bottom-up counterpart of td_inter_mutate_input: the callee re-assigns
 * its input parameter. The bottom-up phase now snapshots the entry value
 * of each input parameter, so the computed summary relates the input's
 * *entry* value with the outputs.
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
  using callgraph_t = call_graph<z_cfg_ref_t>;
  using inter_params_t = inter_analyzer_parameters<callgraph_t>;
  auto cg = std::make_unique<callgraph_t>(cfgs);
  inter_params_t params;
  params.widening_delay = 2;
  params.descending_iters = 2;
  params.thresholds_size = 20;
  // Relational domain for the bottom-up phase so the summary can capture
  // o == p + 1.
  z_dbm_domain_t bu_top;
  z_interval_domain_t td_top;
  bu_inter_run<z_dbm_domain_t, z_interval_domain_t>(*cg, bu_top, td_top, false,
                                                    params, stats_enabled);

  delete t1;
  delete t2;
  return 0;
}
