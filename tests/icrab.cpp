#include "./crab_lang.hpp"
#include "./crab_dom.hpp"

#include "./icrab.hpp"
#include <crab/analysis/inter/bottom_up_inter_analyzer.hpp>
#include <crab/analysis/inter/top_down_inter_analyzer.hpp>

// Helper
template<typename CG, typename BUDom, typename TDDom, typename InterFwdAnalyzer>
void inter_run_impl (CG* cg,
		     BUDom bu_top, TDDom td_top,
		     bool /*run_liveness*/,
		     unsigned widening, 
		     unsigned narrowing, 
		     unsigned jump_set_size,
		     bool enable_stats) {
  
  typedef crab::cg::call_graph_ref<CG> cg_ref_t;
  cg_ref_t cg_ref (*cg);
  
  crab::outs() << "Running " 
	       << "summary domain=" << bu_top.domain_name()
	       << " and forward domain=" << td_top.domain_name() << "\n";
  
  InterFwdAnalyzer a (cg_ref, td_top, bu_top,
		      nullptr /*live*/, widening, narrowing, jump_set_size);
  a.run (td_top);
  
  // Print invariants
  for (auto &v: boost::make_iterator_range (cg_ref.nodes())) {
    auto cfg = v.get_cfg ();
    auto fdecl = cfg.get_func_decl ();
    crab::outs() << fdecl << "\n";      
    for (auto &b : cfg) {
      auto inv = a.get_post (cfg, b.label ());
        crab::outs() <<  crab::cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
    }
      crab::outs() << "=================================\n";
  }
  
  // Print summaries
  for (auto &v: boost::make_iterator_range (cg_ref.nodes())) {
    auto cfg = v.get_cfg ();
    if (a.has_summary (cfg)) {
      auto sum = a.get_summary (cfg);
      crab::outs() << "Summary " << sum << "\n";
    }
  }
  
  if (enable_stats) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}

template<typename Dom, typename InterAnalyzer>
void td_inter_run_impl(crab::cg_impl::z_cg_t* cg,
		       Dom init,
		       td_inter_params_t params,
		       bool print_checks,
		       bool print_invariants,
		       bool enable_stats) {

  InterAnalyzer analyzer(*cg, init, params);
  analyzer.run(init);

  if (print_checks)
    analyzer.print_checks(crab::outs());
  if (print_invariants) {
    // TODO: fix order of cg traversal
    for (auto &v: boost::make_iterator_range(cg->nodes())) {
      auto cfg = v.get_cfg();
      auto fdecl = cfg.get_func_decl();
      crab::outs() << fdecl << "\n";
      
      // Print invariants in DFS to enforce a fixed order
      std::set<crab::cfg_impl::basic_block_label_t> visited;
      std::vector<crab::cfg_impl::basic_block_label_t> worklist;
      worklist.push_back(cfg.entry());
      visited.insert(cfg.entry());
      while (!worklist.empty()) {
	auto cur_label = worklist.back();
	worklist.pop_back();
	auto inv = analyzer.get_pre(cfg, cur_label);
	crab::outs() << crab::cfg_impl::get_label_str (cur_label) << "=" << inv << "\n";
	auto const &cur_node = cfg.get_node (cur_label);
	for (auto const kid_label : boost::make_iterator_range (cur_node.next_blocks ())) {
	  if (visited.insert(kid_label).second) {
	    worklist.push_back(kid_label);
	  }
	}
      }
      crab::outs() << "=================================\n";
    }
  }
  
  if (enable_stats) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}
// To run abstract domains defined over integers
template<typename BUDom, typename TDDom>
void bu_inter_run(crab::cg_impl::z_cg_t* cg,
		  BUDom bu_top, TDDom td_top,
		  bool run_liveness,
		  unsigned widening, 
		  unsigned narrowing, 
		  unsigned jump_set_size,
		  bool enable_stats) {
  using namespace crab::analyzer;
  typedef bottom_up_inter_analyzer<crab::cg_impl::z_cg_ref_t, BUDom, TDDom> inter_analyzer_t;
  inter_run_impl<crab::cg_impl::z_cg_t, BUDom, TDDom, inter_analyzer_t>
    (cg, bu_top, td_top, run_liveness, widening, narrowing, jump_set_size, enable_stats);
}


template<typename Dom>
void td_inter_run(crab::cg_impl::z_cg_t* cg,
		  Dom init,
		  td_inter_params_t params,
		  bool print_checks,
		  bool print_invariants,
		  bool enable_stats) {

  using namespace crab::analyzer;
  using inter_analyzer_t = top_down_inter_analyzer<crab::cg_impl::z_cg_ref_t, Dom>;
  td_inter_run_impl<Dom, inter_analyzer_t>(cg, init, params, print_checks, print_invariants,
					   enable_stats);
}


///////
//// Explicit instantiations here
///////
#include "./icrab_inst.hpp"

