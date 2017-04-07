#include "./crab_lang.hpp"
#include "./crab_dom.hpp"

#include <crab/analysis/liveness.hpp>
#include <crab/analysis/inter_fwd_analyzer.hpp>

// Helper
template<typename CG, typename BUDom, typename TDDom, typename InterFwdAnalyzer>
void inter_run_impl (CG* cg,
		     bool /*run_liveness*/,
		     unsigned widening, 
		     unsigned narrowing, 
		     unsigned jump_set_size,
		     bool enable_stats) {
  
  typedef crab::cg::call_graph_ref<CG> cg_ref_t;
  cg_ref_t cg_ref (*cg);
  
  crab::outs() << "Running " 
	       << "summary domain=" << BUDom::getDomainName () 
	       << "and forward domain=" << TDDom::getDomainName () << "\n";
  
  InterFwdAnalyzer a (cg_ref, nullptr /*live*/, widening, narrowing, jump_set_size);
  a.Run ();
  
  // Print invariants
  for (auto &v: boost::make_iterator_range (cg_ref.nodes())) {
    auto cfg = v.get_cfg ();
    auto fdecl_opt = cfg.get_func_decl ();
    assert (fdecl_opt);
    crab::outs() << *fdecl_opt << "\n";      
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
      auto fdecl_opt = cfg.get_func_decl ();
      assert (fdecl_opt);
      crab::outs() << "Summary for " << *fdecl_opt << ": "; 
      auto sum = a.get_summary (cfg);
      crab::outs() << sum << "\n";
    }
  }
  
  if (enable_stats) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}

// To run abstract domains defined over integers
template<typename BUDom, typename TDDom>
void inter_run (crab::cg_impl::z_cg_t* cg, 
		bool run_liveness,
		unsigned widening, 
		unsigned narrowing, 
		unsigned jump_set_size,
		bool enable_stats) {
  using namespace crab::analyzer;
  typedef inter_fwd_analyzer<crab::cg_impl::z_cg_ref_t, BUDom, TDDom> inter_fwd_analyzer_t;
  inter_run_impl<crab::cg_impl::z_cg_t, BUDom, TDDom, inter_fwd_analyzer_t>
    (cg, run_liveness, widening, narrowing, jump_set_size, enable_stats);
}

///////
//// Explicit instantiations here
///////
#include "./icrab_inst.hpp"

