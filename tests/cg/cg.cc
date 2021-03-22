/** Test that traverses a call graph using different algorithms: dfs,
    scc graph, wto, ... **/

#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>
#include <crab/cg/cg.hpp>
#include <crab/cg/cg_bgl.hpp>

#include <unordered_map>

using namespace crab::analyzer;
using namespace crab::analyzer::graph_algo;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::cg;
using namespace crab::domain_impl;

using call_graph_t = call_graph<z_cfg_ref_t>;
using call_graph_ref_t = call_graph_ref<call_graph_t>;

z_cfg_t *foo(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("foo", {x}, {w});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add(y, x, 1);
  exit.add(z, y, 2);
  exit.callsite("barz", {w}, {z});
  return cfg;
}

z_cfg_t *bar(variable_factory_t &vfac) {
  // Defining program variables
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("bar", {a}, {y});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  exit.callsite("foo", {y}, {x});
  entry.assign(x, a);
  return cfg;
}

z_cfg_t *barz(variable_factory_t &vfac) {
  // Defining program variables
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  function_decl<z_number, varname_t> decl("barz", {a}, {y});
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  exit.callsite("bar", {y}, {x});
  entry.assign(x, a);
  return cfg;
}

z_cfg_t *m(variable_factory_t &vfac) {
  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  std::vector<z_var> inputs, outputs;
  function_decl<z_number, varname_t> decl("main", inputs, outputs);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  entry.callsite("bar", {y}, {x});
  exit.add(z, y, 2);
  return cfg;
}

struct print_visitor : public boost::default_dfs_visitor {
  void discover_vertex(boost::graph_traits<call_graph_t>::vertex_descriptor v,
                       const call_graph_t &g) {
    crab::outs() << v.name() << "\n";
  }
};

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  z_cfg_t *t1 = foo(vfac);
  z_cfg_t *t2 = bar(vfac);
  z_cfg_t *t3 = m(vfac);
  z_cfg_t *t4 = barz(vfac);

  crab::outs() << *t1 << "\n";
  crab::outs() << *t2 << "\n";
  crab::outs() << *t3 << "\n";

  std::vector<z_cfg_ref_t> cfgs_1({*t1, *t2, *t3});
  call_graph_t cg_1(cfgs_1);
  crab::outs() << cg_1 << "\n";

  /// -- Basic BGL api
  for (auto v : boost::make_iterator_range(vertices(cg_1))) {
    if (out_degree(v, cg_1) > 0) {
      crab::outs() << "number of successors " << v.name() << "="
                   << out_degree(v, cg_1) << "\n";
      boost::graph_traits<call_graph_t>::out_edge_iterator ei, ei_end;
      boost::tie(ei, ei_end) = out_edges(v, cg_1);
      for (; ei != ei_end; ++ei) {
        auto s = source(*ei, cg_1);
        auto t = target(*ei, cg_1);
        crab::outs() << s.name() << "-->" << t.name() << "\n";
      }
    }

    if (in_degree(v, cg_1) > 0) {
      crab::outs() << "number of predecessors " << v.name() << "="
                   << in_degree(v, cg_1) << "\n";
      for (auto e : boost::make_iterator_range(in_edges(v, cg_1))) {
        auto s = source(e, cg_1);
        auto t = target(e, cg_1);
        crab::outs() << s.name() << "-->" << t.name() << "\n";
      }
    }
  }

  /// --- DFS
  using color_map_t =
      std::unordered_map<boost::graph_traits<call_graph_t>::vertex_descriptor,
                         boost::default_color_type>;
  color_map_t color;
  for (auto v : boost::make_iterator_range(vertices(cg_1))) {
    color[v] = boost::default_color_type();
  }
  boost::associative_property_map<color_map_t> cm(color);

  // find root
  boost::graph_traits<call_graph_t>::vertex_descriptor root;
  for (auto v : boost::make_iterator_range(vertices(cg_1))) {
    if (in_degree(v, cg_1) == 0) {
      root = v;
      break;
    }
  }
  crab::outs() << "Found root " << root.name() << "\n";

  crab::outs() << "Printing in preorder ...\n";
  print_visitor vis;
  boost::detail::depth_first_visit_impl(cg_1, root, vis, cm,
                                        boost::detail::nontruth2());

  /// --- SccGraph
  scc_graph<call_graph_ref_t> scc_g(cg_1);
  // scc_g.write(crab::outs());

  std::vector<call_graph_ref_t::node_t> order;
  rev_topo_sort(scc_g, order);

  crab::outs() << "reverse topological sort: ";
  for (auto n : order) {
    crab::outs() << n.name() << "--";
  }
  crab::outs() << "\n";
  crab::outs() << "topological sort: ";
  for (auto n : boost::make_iterator_range(order.rbegin(), order.rend())) {
    crab::outs() << n.name() << "--";
  }
  crab::outs() << "\n";

  /// --- WTO

  std::vector<z_cfg_ref_t> cfgs_2({*t1, *t2, *t3, *t4});
  call_graph_t cg_2(cfgs_2);
  using wto_t = wto<call_graph_ref_t>;
  wto_t wto(cg_2);
  crab::outs() << "Callgraph=\n" << cg_2 << "\n";
  crab::outs() << "Weak topological ordering=" << wto << "\n";

  // delete all stuff

  delete t1;
  delete t2;
  delete t3;
  delete t4;

  return 0;
}
