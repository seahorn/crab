#include "../program_options.hpp"

#include "../crab_lang.hpp"
#include "../crab_dom.hpp"

#include <crab/analysis/inter/top_down_inter_analyzer.hpp>

/*
  Example where the callsite and the declaration of the callee share
  variables.

  All the assertions should be proven.
 */
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

z_cfg_t* foo(z_var x, z_var y, z_var z) {
  function_decl<z_number, varname_t> decl("foo", {x,y}, {z});
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& exit   = cfg->insert("exit");
  entry >> exit;
  exit.assign(z, x);
  exit.ret(z);
  return cfg;
}

z_cfg_t* m(z_var x, z_var y, z_var z)  {
  function_decl<z_number, varname_t> decl("main", {}, {});				 
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  z_basic_block_t& entry = cfg->insert("entry");
  z_basic_block_t& exit   = cfg->insert("exit");
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

typedef call_graph<z_cfg_ref_t> callgraph_t;
typedef call_graph_ref<callgraph_t> callgraph_ref_t;
typedef top_down_inter_analyzer_parameters<callgraph_ref_t> inter_analyzer_params_t;

template<typename InterAnalyzer>
void run(callgraph_t& cg) {
  InterAnalyzer default_analyzer(cg);
  default_analyzer.run();
  // Print invariants
  #if 0
  for(auto &v: boost::make_iterator_range(cg.nodes())) {
    auto cfg = v.get_cfg();
    auto fdecl = cfg.get_func_decl();
    crab::outs() << fdecl << "\n";      
    for(auto &b : cfg) {
      auto pre_inv = default_analyzer.get_pre(cfg, b.label());
      auto post_inv = default_analyzer.get_post(cfg, b.label());
      crab::outs() <<  crab::cfg_impl::get_label_str(b.label()) << ":\n"
		   << "\t" << pre_inv << " ==>\n\t" << post_inv << "\n";
    }
    crab::outs() << "=================================\n";
  }
  #endif 
  default_analyzer.print_checks(crab::outs());
}

int main(int argc, char** argv) {
  bool stats_enabled = false;
  if(!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  
  z_cfg_t* t1 = foo(x,y,z);
  z_cfg_t* t2 = m(x,y,z);

  crab::outs() << *t1 << "\n"
	       << *t2 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  callgraph_t cg(cfgs);
  {
    typedef top_down_inter_analyzer<callgraph_ref_t, z_sdbm_domain_t> inter_analyzer_t;
    crab::outs() << "Running top-down inter-procedural analysis with "
		 << z_sdbm_domain_t::getDomainName() << "\n";
    run<inter_analyzer_t>(cg);
  }
  
  
  delete t1;
  delete t2;

  return 0;
}
