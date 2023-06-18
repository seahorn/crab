#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/dataflow/assertion_crawler.hpp>
#include <crab/analysis/dataflow/assumptions.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>
#include <crab/checkers/div_zero.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::checker;

z_cfg_t *get_cfg(variable_factory_t &vfac) {
  /*
    i := 0;
    x := 1;
    y := 0;
    while(i <= 99) {
       x+=y;
       y++;
       i++;
    }
    assume(x<=y);
    assert(i>=200);
    assert(x>=y);
    assert(x>=200);
   */

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;

  // adding statements
  entry.assign(i, 0);
  entry.assign(x, 1);
  entry.assign(y, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(x, x, y);
  bb2.add(y, y, 1);
  bb2.add(i, i, 1);
  ret.assume(x <= y);
  ret.assertion(i == 100);
  ret.assertion(i >= 200);
  ret.assertion(x >= y);
  ret.assertion(x >= 200);

  return cfg;
}

// Print invariants by traversing the cfg in dfs.
template <typename Analyzer>
static void print_invariants(z_cfg_ref_t cfg, Analyzer &analyser) {
  std::set<crab::cfg_impl::basic_block_label_t> visited;
  std::vector<crab::cfg_impl::basic_block_label_t> worklist;
  worklist.push_back(cfg.entry());
  visited.insert(cfg.entry());
  while (!worklist.empty()) {
    auto cur_label = worklist.back();
    worklist.pop_back();

    auto pre = analyser.get_pre(cur_label);
    auto post = analyser.get_post(cur_label);
    
    crab::outs() << crab::basic_block_traits<z_basic_block_t>::to_string(cur_label)
		 << "=" << pre << " ==> " << post
                 << "\n";

    auto const &cur_node = cfg.get_node(cur_label);
    for (auto const &kid_label :
         boost::make_iterator_range(cur_node.next_blocks())) {
      if (visited.insert(kid_label).second) {
        worklist.push_back(kid_label);
      }
    }
  }
}

void check(z_cfg_ref_t cfg, variable_factory_t &vfac) {

  // Each checker is associated to one analyzer
  using num_analyzer_t = intra_fwd_analyzer<z_cfg_ref_t, z_sdbm_domain_t>;
  using num_checker_t = intra_checker<num_analyzer_t>;

  // We can have multiple properties per analyzer
  using div_zero_prop_num_checker_t = div_zero_property_checker<num_analyzer_t>;
  using assert_prop_num_checker_t = assert_property_checker<num_analyzer_t>;

  crab::outs() << cfg << "\n";

  // Run analyses
  z_sdbm_domain_t absval_fac, init;
  crab::fixpoint_parameters fixpo_params;    
  num_analyzer_t num_a(cfg, absval_fac, nullptr, fixpo_params);
  num_a.run(init);
  crab::outs() << "Analysis using " << init.domain_name() << "\n";
  print_invariants(cfg, num_a);

  // for(auto &b : cfg)  {
  //   auto pre = num_a.get_pre(b.label());
  //   auto post = num_a.get_post(b.label());
  //   crab::outs() << basic_block_traits<z_basic_block_t>::to_string(b.label()) << "="
  //             << pre
  //             << " ==> "
  //             << post << "\n";
  // }

  // Run the checkers with several properties
  // A checker can take any property checker associated to same
  // analyzer.
  const int verbose = 3;

  typename num_checker_t::prop_checker_ptr prop1(
      new div_zero_prop_num_checker_t(verbose));
  typename num_checker_t::prop_checker_ptr prop2(
      new assert_prop_num_checker_t(verbose));
  num_checker_t checker(num_a, {prop1, prop2});

  checker.run();
  checker.show(crab::outs());
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  z_cfg_t *p = get_cfg(vfac);

  // To test the assertion crawler analysis
  using assertion_crawler_t = crab::analyzer::assertion_crawler<z_cfg_ref_t>;
  typename assertion_crawler_t::assert_map_t assert_map;
  typename assertion_crawler_t::summary_map_t summaries;
  
  assertion_crawler_t assert_crawler(*p, assert_map, summaries);
  assert_crawler.exec();
  crab::outs() << "\n";
  assert_crawler.write(crab::outs());
  crab::outs() << "\n";

  // To test the assumption analyses
  crab::analyzer::assumption_naive_analysis<z_cfg_ref_t>
      assumption_naive_analyzer(*p);
  assumption_naive_analyzer.exec();
  crab::outs() << "\n" << assumption_naive_analyzer << "\n";

  crab::analyzer::assumption_dataflow_analysis<z_cfg_ref_t>
      assumption_dataflow_analyzer(*p);
  assumption_dataflow_analyzer.exec();
  crab::outs() << "\n" << assumption_dataflow_analyzer << "\n";

  check(*p, vfac);

  delete p;

  return 0;
}
