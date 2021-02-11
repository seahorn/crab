#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/graphs/sccg_bgl.hpp>
#include <crab/analysis/inter/top_down_inter_params.hpp>
#include <crab/cg/cg_bgl.hpp>

/*
  Inter-procedural analysis with region domain.
  foo has different variable names at the callsite and as formal parameters.

  All the assertions should be proven.

  int* foo(int a) {
    int *b= malloc(...);
    tmp = nd_int();
    assume(tmp > 0);
    assume(tmp <= a);
    *b = tmp;
    return b;
  }

  void  main() {
    int x = nd_int();
    int *z = call foo(x);
    assert(*z <= x);
    assert(*z > 0);
  }
 */
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;
using namespace crab::cg_impl;

z_cfg_t *foo(variable_factory_t &vfac, crab::allocation_site_man &as_man) {
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::REG_INT_TYPE, 32);  
  z_var ref(vfac["ref"], crab::REF_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);    
  
  function_decl<z_number, varname_t> decl("foo", {a}, {b});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  entry.region_init(b);
  entry.make_ref(ref, b, as_man.mk_allocation_site());
  entry.assume(tmp <= a);
  entry.assume(tmp > 0);
  exit.store_to_ref(ref, b, tmp);
  exit.ret(b);
  return cfg;
}

z_cfg_t *__main(variable_factory_t &vfac, crab::allocation_site_man &as_man) {
  function_decl<z_number, varname_t> decl("main", {}, {});
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::REG_INT_TYPE, 32);
  
  entry >> exit;
  entry.havoc(x);
  entry.assume(x > 0);
  exit.callsite("foo", {z}, {x});
  z_var ref(vfac["ref"], crab::REF_TYPE, 32);
  z_var lhs(vfac["lhs"], crab::INT_TYPE, 32);    
  exit.load_from_ref(lhs, ref, z);
  exit.assertion(z_lin_t(x) >= z_lin_t(lhs));
  exit.assertion(lhs > 0);  
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  using inter_params_t = top_down_inter_analyzer_parameters<z_cg_ref_t>;
  variable_factory_t vfac;
  crab::allocation_site_man as_man;
  // Defining program variables
  
  z_cfg_t *t1 = foo(vfac, as_man);
  z_cfg_t *t2 = __main(vfac, as_man);

  crab::outs() << *t1 << "\n" << *t2 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  z_rgn_sdbm_t init;
  crab::outs() << "Running top-down inter-procedural analysis with "
               << init.domain_name() << "\n";
  z_cg_t cg(cfgs);
  inter_params_t params;
  td_inter_run(&cg, init, params, true, false/*print invariants*/, false);

  delete t1;
  delete t2;

  return 0;
}
