#include "./crab_lang.hpp"
#include "./crab_dom.hpp"

#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>

#include <crab/checkers/base_property.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/checker.hpp>

// Helper
template<typename CFG, typename Dom, typename IntraFwdAnalyzer>
void intra_run_impl (CFG* cfg, 
		     bool run_liveness,
		     unsigned widening, 
		     unsigned narrowing, 
		     unsigned jump_set_size,
		     bool enable_stats,
		     bool enable_checker){
  typedef crab::cfg::cfg_ref<CFG> cfg_ref_t;  
  crab::analyzer::liveness<cfg_ref_t> *live = nullptr;
  if (run_liveness) {
    crab::analyzer::liveness<cfg_ref_t> live_(*cfg);
    live_.exec ();
    live=&live_;
  }
  // Run fixpoint
  Dom inv = Dom::top ();        
  crab::outs() << "Invariants using " << inv.getDomainName () << "\n";
  IntraFwdAnalyzer a (*cfg, inv, live, widening, narrowing, jump_set_size);
  
  a.run ();
  
  // Print invariants
  for (auto &b : *cfg) {
    auto inv = a[b.label ()];
    crab::outs() << crab::cfg_impl::get_label_str (b.label ()) << "="
		 << inv << "\n";
  }
  // Print abstract trace.
  // We consider as an abtract trace the WTO annotated with how many
  // times a cycle is visited by the fixpoint.
  auto &wto = a.get_wto ();
  crab::outs () << "Abstract trace: " << wto << "\n";

  if (enable_checker) {
    typedef crab::checker::intra_checker<IntraFwdAnalyzer> checker_t;
    typedef crab::checker::assert_property_checker<IntraFwdAnalyzer> assert_checker_t;  
    const int verbose = 3;
    typename checker_t::prop_checker_ptr prop (new assert_checker_t (verbose));
    checker_t checker (a, {prop});
    checker.run ();
    checker.show (crab::outs());
  }
  
  crab::outs() << "\n";
  if (enable_stats) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }

}

// To run abstract domains defined over integers
template<typename Dom>
void run (crab::cfg_impl::z_cfg_t* cfg, 
	  bool run_liveness,
	  unsigned widening, 
	  unsigned narrowing, 
	  unsigned jump_set_size,
	  bool enable_stats) {
  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom> intra_fwd_analyzer_t;
  intra_run_impl<crab::cfg_impl::z_cfg_t, Dom, intra_fwd_analyzer_t>
    (cfg, run_liveness, widening, narrowing, jump_set_size, enable_stats, false);
}

template<typename Dom>
void run_and_check (crab::cfg_impl::z_cfg_t* cfg, 
		    bool run_liveness,
		    unsigned widening, 
		    unsigned narrowing, 
		    unsigned jump_set_size,
		    bool enable_stats) {
  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom> intra_fwd_analyzer_t;
  intra_run_impl<crab::cfg_impl::z_cfg_t, Dom, intra_fwd_analyzer_t>
    (cfg, run_liveness, widening, narrowing, jump_set_size, enable_stats, true);
}

// To run abstract domains defined over rationals
template<typename Dom>
void run (crab::cfg_impl::q_cfg_t* cfg, 
	  bool run_liveness,
	  unsigned widening, 
	  unsigned narrowing, 
	  unsigned jump_set_size,
	  bool enable_stats) {
  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::q_cfg_ref_t,Dom> intra_fwd_analyzer_t;
  intra_run_impl<crab::cfg_impl::q_cfg_t, Dom, intra_fwd_analyzer_t>
    (cfg, run_liveness, widening, narrowing, jump_set_size, enable_stats, false);
}

////////
//// Explicit instantiations here
////////

#include "./crab_inst.hpp"
