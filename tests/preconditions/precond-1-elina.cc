#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/bwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* CheckedC example: args.c */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var argc(vfac["argc"], crab::INT_TYPE, 32);
  z_var len(vfac["len"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = new z_cfg_t("bb0", "bb6");
  /*
    if (argc > 2)
       NumNodes = atoi(argv[2]);
    else
       NumNodes = 4;

    if (argc > 1)
      nbody = atoi(argv[1]);
    else
      nbody = 32;
  */

  // adding blocks
  z_basic_block_t &bb0 = cfg->insert("bb0");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &bb4 = cfg->insert("bb4");
  z_basic_block_t &bb5 = cfg->insert("bb5");
  z_basic_block_t &bb6 = cfg->insert("bb6");

  bb0 >> bb1;
  bb0 >> bb2;
  bb1 >> bb3;
  bb2 >> bb3;

  bb3 >> bb4;
  bb3 >> bb5;
  bb4 >> bb6;
  bb5 >> bb6;

  bb1.assume(argc >= 3);
  bb1.assertion(2 <= len - 1);

  bb4.assume(argc >= 2);
  bb4.assertion(1 <= len - 1);

  return cfg;
}

int main(int argc, char **argv) {
#ifdef HAVE_ELINA

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
            z_cfg_ref_t, z_pk_elina_domain_t>;
    z_pk_elina_domain_t absval_fac;
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
            z_cfg_ref_t, z_pk_elina_domain_t>;
    z_pk_elina_domain_t absval_fac;
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
