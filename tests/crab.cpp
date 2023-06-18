#include "./crab_dom.hpp"
#include "./crab_lang.hpp"

#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>

#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>

// Helper
template <typename CFG, typename Dom, typename IntraFwdAnalyzer>
void intra_run_impl(CFG *cfg, crab::cfg_impl::basic_block_label_t entry,
                    Dom init, bool run_liveness, unsigned widening,
                    unsigned narrowing, unsigned jump_set_size,
                    bool enable_stats, bool enable_checker, bool print_invariants) {
  using cfg_ref_t = crab::cfg::cfg_ref<CFG>;
  using basic_block_t = typename CFG::basic_block_t;
  using assumption_map_t = typename IntraFwdAnalyzer::assumption_map_t;
  
  crab::analyzer::live_and_dead_analysis<cfg_ref_t> live(*cfg);
  if (run_liveness) {
    live.exec();
  }
  // Run fixpoint
  auto absval_fac = init.make_top();  
  crab::outs() << "Invariants using " << absval_fac.domain_name() << "\n";
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening;
  fixpo_params.get_descending_iterations() = narrowing;
  fixpo_params.get_max_thresholds() = jump_set_size;
  IntraFwdAnalyzer a(*cfg, absval_fac, (run_liveness) ? &live : nullptr, 
                     fixpo_params);

  assumption_map_t assumptions;
  a.run(entry, init, assumptions);

  CRAB_LOG("crab-tests-print-invariants",
	   print_invariants = true
	   );
  
  if (print_invariants) {
    // Print invariants in DFS to enforce a fixed order
    std::set<crab::cfg_impl::basic_block_label_t> visited;
    std::vector<crab::cfg_impl::basic_block_label_t> worklist;
    worklist.push_back(cfg->entry());
    visited.insert(cfg->entry());
    while (!worklist.empty()) {
      auto cur_label = worklist.back();
      worklist.pop_back();
      auto inv = a[cur_label];
      crab::outs() << crab::basic_block_traits<basic_block_t>::to_string(cur_label)
		   << "=" << inv << "\n";
      auto const &cur_node = cfg->get_node(cur_label);
      for (auto const& kid_label :
         boost::make_iterator_range(cur_node.next_blocks())) {
	if (visited.insert(kid_label).second) {
	  worklist.push_back(kid_label);
	}
      }
    }
  }

  // for (auto &b : *cfg) {
  //   auto inv = a[b.label ()];
  //   crab::outs() << crab::basic_block_traits<basic_block_t>::to_string(b.label()) << "="
  // 		 << inv << "\n";
  // }

  // Print abstract trace.
  // We consider as an abtract trace the WTO annotated with how many
  // times a cycle is visited by the fixpoint.
  auto &wto = a.get_wto();
  crab::outs() << "Abstract trace: " << wto << "\n";

  if (enable_checker) {
    using checker_t = crab::checker::intra_checker<IntraFwdAnalyzer>;
    using assert_checker_t =
        crab::checker::assert_property_checker<IntraFwdAnalyzer>;
    const int verbose = 0;
    typename checker_t::prop_checker_ptr prop(new assert_checker_t(verbose));
    checker_t checker(a, {prop});
    checker.run();
    checker.show(crab::outs());
  }

  crab::outs() << "\n";
  if (enable_stats) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }
}

template <typename Dom>
void z_intra_run(crab::cfg_impl::z_cfg_t *cfg,
                 crab::cfg_impl::basic_block_label_t entry, Dom init,
                 bool run_liveness, unsigned widening, unsigned narrowing,
                 unsigned jump_set_size, bool enable_stats,
                 bool enable_checker, bool print_invariants) {
  using namespace crab::analyzer;
  using intra_fwd_analyzer_t =
      intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom>;
  intra_run_impl<crab::cfg_impl::z_cfg_t, Dom, intra_fwd_analyzer_t>(
      cfg, entry, init, run_liveness, widening, narrowing, jump_set_size,
      enable_stats, enable_checker, print_invariants);
}

template <typename Dom>
void q_intra_run(crab::cfg_impl::q_cfg_t *cfg,
                 crab::cfg_impl::basic_block_label_t entry, Dom init,
                 bool run_liveness, unsigned widening, unsigned narrowing,
                 unsigned jump_set_size, bool enable_stats,
                 bool enable_checker, bool print_invariants) {
  using namespace crab::analyzer;
  using intra_fwd_analyzer_t =
      intra_fwd_analyzer<crab::cfg_impl::q_cfg_ref_t, Dom>;
  intra_run_impl<crab::cfg_impl::q_cfg_t, Dom, intra_fwd_analyzer_t>(
      cfg, entry, init, run_liveness, widening, narrowing, jump_set_size,
      enable_stats, enable_checker, print_invariants);
}

////////
//// Explicit instantiations here
////////

#include "./crab_inst.hpp"
