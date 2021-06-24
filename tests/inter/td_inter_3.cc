#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter/inter_params.hpp>
#include <crab/cg/cg_bgl.hpp>

/*
  Example where the callsite and the declaration of the callee share
  variables.

  All the assertions should be proven.

  int foo(int x, int y) {
    int z = x;
    return z;
  }

  void  main() {
    int x = 5;
    int y = 10;
    int z = 8;
    z = call foo(x, y);
    assert(x = 5);
    assert(y = 10);
    assert(z = 5);
  }
 */
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;
using namespace crab::cg_impl;

z_cfg_t *foo(z_var x, z_var y, z_var z) {
  function_decl<z_number, varname_t> decl("foo", {x, y}, {z});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  exit.assign(z, x);
  return cfg;
}

z_cfg_t *m(z_var x, z_var y, z_var z) {
  function_decl<z_number, varname_t> decl("main", {}, {});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  entry.assign(x, 5);
  entry.assign(y, 10);
  entry.assign(z, 8);
  exit.callsite("foo", {z}, {x, y});
  exit.assertion(x == 5);
  exit.assertion(y == 10);
  exit.assertion(z == 5);
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  using inter_params_t = inter_analyzer_parameters<z_cg_t>;
  variable_factory_t vfac;
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  z_cfg_t *t1 = foo(x, y, z);
  z_cfg_t *t2 = m(x, y, z);

  crab::outs() << *t1 << "\n" << *t2 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  z_sdbm_domain_t init;
  crab::outs() << "Running top-down inter-procedural analysis with "
               << init.domain_name() << "\n";
  z_cg_t cg(cfgs);
  inter_params_t params;
  td_inter_run(cg, init, params, true, false, false);

  delete t1;
  delete t2;

  return 0;
}
