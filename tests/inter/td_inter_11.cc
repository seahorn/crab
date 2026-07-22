#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/inter/inter_params.hpp>
#include <memory>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;
/*
  int inc(a) {
    b = a + 1;
    return b;
  }

  void main() {
    x = 5;
    y = inc(x);
    assert(y == 6);
  }
 */
std::unique_ptr<z_cfg_t> inc(variable_factory_t &vfac) {
  // Defining program variables
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("inc", {a}, {b});
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);
  // adding blocks
  BB(cfg, entry);
  BB(cfg, exit);
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add(b, a, 1);
  return cfg;
}

std::unique_ptr<z_cfg_t> _main(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("main", {}, {});
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);
  // adding blocks
  BB(cfg, entry);
  BB(cfg, exit);
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 5);
  entry.callsite("inc", {y}, {x});
  exit.assertion(y == 6);
  return cfg;
}

using callgraph_t = call_graph<z_cfg_ref_t>;
using inter_params_t = inter_analyzer_parameters<callgraph_t>;

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto t1 = inc(vfac);
    auto t2 = _main(vfac);

    crab::outs() << *t1 << "\n";
    crab::outs() << *t2 << "\n";

    vector<z_cfg_ref_t> cfgs({*t1, *t2});
    callgraph_t cg(cfgs);
    crab::outs() << "CallGraph=" << cg << "\n";
    z_sdbm_domain_t init;
    crab::outs() << "Running top-down inter-procedural analysis with "
                 << init.domain_name() << "\n";
    inter_params_t params;
    td_inter_run(cg, init, params, true, true, false);

    return 0;
  });
}
