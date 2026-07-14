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
using namespace crab::cg_impl;

/*
loop(x,y) {
  if (x <= 100) {
     (x',y') = loop(x+1,y+1)
     return (x',y');
  } else {
     return (x,y);
  }
}

main() {
  (x,y) := loop(0,0);
}
 */

std::unique_ptr<z_cfg_t> f_loop(variable_factory_t &vfac) {
  // Defining program variables
  z_var x_in(vfac["x_in"], crab::INT_TYPE, 32);
  z_var y_in(vfac["y_in"], crab::INT_TYPE, 32);
  z_var x_out(vfac["x_out"], crab::INT_TYPE, 32);
  z_var y_out(vfac["y_out"], crab::INT_TYPE, 32);
  z_var x_next(vfac["x_next"], crab::INT_TYPE, 32);
  z_var y_next(vfac["y_next"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("loop", {x_in, y_in}, {x_out, y_out});
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);  
  // adding blocks
  BB(cfg, entry);
  BB(cfg, rec_case);
  BB(cfg, base_case);
  BB(cfg, exit);
  // adding control flow
  entry >> rec_case;
  entry >> base_case;
  rec_case >> exit;
  base_case >> exit;
  // adding statements
  rec_case.assume(x_in <= 100);
  rec_case.add(x_next, x_in, 1);
  rec_case.add(y_next, y_in, 1);
  rec_case.callsite("loop", {x_out, y_out}, {x_next, y_next});
  base_case.assume(x_in >= 101);
  base_case.assign(x_out, x_in);
  base_case.assign(y_out, y_in);
  return cfg;
}

std::unique_ptr<z_cfg_t> m(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var x_res(vfac["x_res"], crab::INT_TYPE, 32);
  z_var y_res(vfac["y_res"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("main", {}, {x_res, y_res});

  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);  
  // adding blocks
  BB(cfg, entry);
  BB(cfg, exit);
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 0);

  entry.callsite("loop", {x_res, y_res}, {x, y});
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    using inter_params_t = inter_analyzer_parameters<z_cg_t>;

    variable_factory_t vfac;
    auto t1 = f_loop(vfac);
    auto t2 = m(vfac);

    crab::outs() << *t1 << "\n" << *t2 << "\n";

    vector<z_cfg_ref_t> cfgs({*t1, *t2});

    inter_params_t params;
    params.max_call_contexts = 5;
    //params.analyze_recursive_functions = true;  
    z_dbm_domain_t init;
    z_cg_t cg(cfgs);

    td_inter_run(cg, init, params, false, true, false);

    return 0;
  });
}
