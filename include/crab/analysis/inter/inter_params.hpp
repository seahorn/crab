#pragma once

#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/fixpoint/wto.hpp>
#include <unordered_map>

namespace crab {
namespace analyzer {

/* Parameters for the inter-procedural analyses */
template <typename CallGraph>
struct inter_analyzer_parameters {
  using cg_node_t = typename CallGraph::node_t;
  using cfg_t = typename cg_node_t::cfg_t;
  using liveness_t = live_and_dead_analysis<cfg_t>;
  using liveness_map_t = std::unordered_map<cfg_t, const liveness_t *>;

  inter_analyzer_parameters()
      : only_main_as_entry(false), widening_delay(2),
        descending_iters(2 /*UINT_MAX*/), thresholds_size(20 /*0*/),
        live_map(nullptr), 
	run_checker(true), checker_verbosity(0), keep_cc_invariants(false),
        keep_invariants(true), max_call_contexts(UINT_MAX),
        analyze_recursive_functions(false), exact_summary_reuse(true) {}

  // Start the analysis from main
  bool only_main_as_entry;
  // Fixpoint parameters
  unsigned int widening_delay;
  unsigned int descending_iters;
  unsigned int thresholds_size;
  // Use liveness symbols if already available
  const liveness_map_t *live_map;
  
  // -- Begin parameters for top-down analysis -- //
  // whether interleave analysis with checking
  bool run_checker;
  unsigned checker_verbosity;
  // whether keep context-sensitive invariants
  bool keep_cc_invariants;
  // whether keep context-insensitive invariants
  bool keep_invariants;
  // maximum number of calling contexts per callsite
  unsigned int max_call_contexts;
  // analyze precisely recursive functions
  bool analyze_recursive_functions;
  // reuse summaries without losing precision
  bool exact_summary_reuse;
  // -- End parameters for top-down analysis -- //
  
};

} // namespace analyzer
} // namespace crab
