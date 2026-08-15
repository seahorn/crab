#pragma once

#include <crab/analysis/graphs/dominance.hpp>
#include <crab/support/print.hpp>
//#include <crab/cfg/basic_block_traits.hpp>
/*

  Node y is control-dependent on x if y does NOT post-dominate x but
  there exists a path from x to y such that all nodes in the path
  (different from x and y) are post-dominated by y.

 */
namespace crab {
namespace analyzer {
namespace graph_algo {

// OUT: cdeps is a map from nodes to a set of nodes that
// control-dependent on it.
template <typename G, typename VectorMap>
void control_dep_graph(G g, VectorMap &cdg) {
  VectorMap pdf;
  using basic_block_t = typename G::basic_block_t;
  crab::analyzer::graph_algo::post_dominance(g, pdf);

  for (auto &kv : pdf) {
    for (auto v : kv.second) {
      auto &cdeps = cdg[v];
      if (std::find(cdeps.begin(), cdeps.end(), kv.first) == cdeps.end()) {
        cdeps.push_back(kv.first);
      }
    }
  }

  CRAB_LOG(
      "cdg", crab::outs() << "Control-dependence graph \n"; for (auto &kv
                                                                 : cdg) {
        print::print_range_with(
            crab::outs(), kv.second,
            [](crab::crab_os &o, const auto &v) {
              o << crab::basic_block_traits<basic_block_t>::to_string(v);
            },
            print::fmt_debug());
        crab::outs() << "  control-dependent on "
                     << crab::basic_block_traits<basic_block_t>::to_string(
                            kv.first)
                     << "\n";
      });
}

} // namespace graph_algo
} // namespace analyzer
} // namespace crab
