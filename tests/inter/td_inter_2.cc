#include "../program_options.hpp"

#include "../crab_lang.hpp"
#include "../crab_dom.hpp"

#include <crab/analysis/inter/top_down_inter_analyzer.hpp>


using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cg;

z_cfg_t* f_loop (variable_factory_t &vfac) {
  // Defining program variables
  z_var x_in(vfac["x_in"], crab::INT_TYPE, 32);
  z_var y_in(vfac["y_in"], crab::INT_TYPE, 32);
  z_var x_out(vfac["x_out"], crab::INT_TYPE, 32);
  z_var y_out(vfac["y_out"], crab::INT_TYPE, 32);
  z_var x_next(vfac["x_next"], crab::INT_TYPE, 32);
  z_var y_next(vfac["y_next"], crab::INT_TYPE, 32);
  
  function_decl<z_number, varname_t> decl ("loop", {x_in,y_in}, {x_out, y_out});
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& rec_case = cfg->insert ("rec_case");
  z_basic_block_t& base_case = cfg->insert ("base_case");    
  z_basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> rec_case;
  entry >> base_case;
  rec_case >> exit;
  base_case >> exit;
  // adding statements
  rec_case.assume(x_in <= 100);
  rec_case.add(x_next, x_in, 1);
  rec_case.add(y_next, y_in, 1);
  rec_case.callsite("loop",{x_out, y_out}, {x_next, y_next});
  base_case.assume(x_in >= 101);
  base_case.assign(x_out, x_in);
  base_case.assign(y_out, y_in);
  exit.ret({x_out, y_out});
  return cfg;
}

z_cfg_t* m(variable_factory_t &vfac)  {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var x_res(vfac["x_res"], crab::INT_TYPE, 32);
  z_var y_res(vfac["y_res"], crab::INT_TYPE, 32);
  
  function_decl<z_number, varname_t> decl("main", {}, {x_res,y_res});
				 
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 0);
  
  entry.callsite ("loop", {x_res,y_res}, {x,y});
  exit.ret({x_res,y_res});
  return cfg;
}

int main (int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  z_cfg_t* t1 = f_loop(vfac);
  z_cfg_t* t2 = m(vfac);

  crab::outs() << *t1 << "\n"
	       << *t2 << "\n";

  vector<z_cfg_ref_t> cfgs({*t1, *t2});
  typedef call_graph<z_cfg_ref_t> callgraph_t;
  typedef call_graph_ref<callgraph_t> callgraph_ref_t;
  typedef top_down_inter_analyzer<callgraph_ref_t, z_dbm_domain_t> inter_analyzer_t;
  typedef top_down_inter_analyzer_parameters<callgraph_ref_t> inter_params_t;        

  inter_params_t params;
  params.max_call_contexts = 5;
  callgraph_t cg(cfgs);  
  inter_analyzer_t analyzer(cg, params);
  analyzer.run();

  // Print invariants
  for (auto &v: boost::make_iterator_range(cg.nodes())) {
    auto cfg = v.get_cfg();
    auto fdecl = cfg.get_func_decl();
    crab::outs() << fdecl << "\n";      
    for (auto &b : cfg) {
      auto inv = analyzer.get_pre(cfg, b.label());
      crab::outs() <<  crab::cfg_impl::get_label_str(b.label()) << "=" << inv << "\n";
    }
      crab::outs() << "=================================\n";
  }
  
  delete t1;
  delete t2;

  return 0;
}
