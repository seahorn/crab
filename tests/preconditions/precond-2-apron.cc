#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/bwd_analyzer.hpp>

// Simplified Cousots, Logozzo, and Fahndrich's example (VMCAI'13)
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("bb1", "bb5");
  // adding blocks
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  z_basic_block_t &bb5 = cfg->insert("bb5");
  // adding control flow
  bb1 >> bb2;
  bb2 >> bb3;
  bb3 >> bb4;
  bb4 >> bb2;
  bb4 >> bb5;
  bb2 >> bb5;

  // adding statements
  bb1.assign(i, 0);
  bb3.assume(i <= n);
  bb3.assertion(i >= 0);
  bb3.assertion(i <= n - 1);
  bb3.add(i, i, 1);
  bb5.assume(i >= 100);
  return cfg;
}

int main(int argc, char **argv) {
#ifdef HAVE_APRON

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";
  
  {
    using analysis_t =
        crab::analyzer::necessary_preconditions_fixpoint_iterator<
            z_cfg_ref_t, z_pk_apron_domain_t>;
    z_pk_apron_domain_t absval_fac;
    crab::fixpoint_parameters fixpo_params;
    analysis_t analyzer(*cfg, absval_fac, false /*error states*/, fixpo_params);
    analyzer.run_backward(absval_fac.make_bottom());
    crab::outs()
        << "Necessary preconditions from error states using Polyhedra:\n";
    // Print preconditions in DFS to enforce a fixed order
    std::set<crab::cfg_impl::basic_block_label_t> visited;
    std::vector<crab::cfg_impl::basic_block_label_t> worklist;
    worklist.push_back(cfg->entry());
    visited.insert(cfg->entry());
    while (!worklist.empty()) {
      auto cur_label = worklist.back();
      worklist.pop_back();
      auto inv = analyzer[cur_label];
      crab::outs() << crab::basic_block_traits<z_basic_block_t>::to_string(cur_label)
		   << "=" << inv << "\n";
      auto const &cur_node = cfg->get_node(cur_label);
      for (auto const &kid_label :
           boost::make_iterator_range(cur_node.next_blocks())) {
        if (visited.insert(kid_label).second) {
          worklist.push_back(kid_label);
        }
      }
    }
  }
  {
    using analysis_t =
        crab::analyzer::necessary_preconditions_fixpoint_iterator<
            z_cfg_ref_t, z_pk_apron_domain_t>;
    z_pk_apron_domain_t absval_fac;
    crab::fixpoint_parameters fixpo_params;    
    analysis_t analyzer(*cfg, absval_fac, true /*good states*/, fixpo_params);
    analyzer.run_backward(absval_fac.make_top());
    crab::outs()
        << "Necessary preconditions from safe states using Polyhedra:\n";
    // Print preconditions in DFS to enforce a fixed order
    std::set<crab::cfg_impl::basic_block_label_t> visited;
    std::vector<crab::cfg_impl::basic_block_label_t> worklist;
    worklist.push_back(cfg->entry());
    visited.insert(cfg->entry());
    while (!worklist.empty()) {
      auto cur_label = worklist.back();
      worklist.pop_back();
      auto inv = analyzer[cur_label];
      crab::outs() << crab::basic_block_traits<z_basic_block_t>::to_string(cur_label)
		   << "=" << inv << "\n";
      auto const &cur_node = cfg->get_node(cur_label);
      for (auto const &kid_label :
           boost::make_iterator_range(cur_node.next_blocks())) {
        if (visited.insert(kid_label).second) {
          worklist.push_back(kid_label);
        }
      }
    }
  }
  // free the CFG
  delete cfg;
#endif
  return 0;
}
