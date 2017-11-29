#include "../program_options.hpp"
#include "../common.hpp"

#include <crab/cg/cg.hpp>
#include <crab/cg/cg_bgl.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>


using namespace std;
using namespace crab::analyzer;
using namespace crab::analyzer::graph_algo;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::cg;
using namespace crab::domain_impl;


using namespace boost;

typedef call_graph<z_cfg_ref_t> call_graph_t;
typedef call_graph_ref<call_graph_t> call_graph_ref_t;
typedef ikos::variable<z_number, varname_t> variable_t;

z_cfg_t* foo (variable_factory_t &vfac) {
  // Defining program variables
  z_var x (vfac ["x"], crab::INT_TYPE);
  z_var y (vfac ["y"], crab::INT_TYPE);
  z_var z (vfac ["z"], crab::INT_TYPE);
  
  function_decl<z_number, varname_t> decl (vfac["foo"], {x}, {z});
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add (y, x, 1);
  exit.add (z , y , 2);
  exit.ret (z);
  return cfg;
}


z_cfg_t* bar (variable_factory_t &vfac) {
  // Defining program variables
  z_var a (vfac ["a"], crab::INT_TYPE);
  z_var x (vfac ["x"], crab::INT_TYPE);
  z_var y (vfac ["y"], crab::INT_TYPE);
  
  function_decl<z_number, varname_t> decl (vfac["bar"], {a}, {y});
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  exit.callsite (vfac ["foo"], {y}, {x});
  entry.assign (x, a);
  exit.ret (y);
  return cfg;
}

z_cfg_t* m (variable_factory_t &vfac)  {
  // Defining program variables
  z_var x (vfac ["x"], crab::INT_TYPE);
  z_var y (vfac ["y"], crab::INT_TYPE);
  z_var z (vfac ["z"], crab::INT_TYPE);
  
  vector<variable_t> inputs, outputs;
  function_decl<z_number, varname_t> decl (vfac["main"], inputs, outputs);
  // entry and exit block
  z_cfg_t* cfg = new z_cfg_t("entry", "exit", decl);
  // adding blocks
  z_basic_block_t& entry = cfg->insert ("entry");
  z_basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  entry.callsite (vfac ["bar"], {y}, {x});
  exit.add (z, y, 2);
  return cfg;
}

struct print_visitor: public boost::default_dfs_visitor {
  void discover_vertex(graph_traits<call_graph_t>::vertex_descriptor v, 
                       const call_graph_t& g)  {
    crab::outs() << v.name () << "\n";
  }
};


int main (int argc, char** argv ) {
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  z_cfg_t* t1 = foo (vfac);
  z_cfg_t* t2 = bar (vfac);
  z_cfg_t* t3 = m (vfac);

  crab::outs() << *t1 << "\n";
  crab::outs() << *t2 << "\n";
  crab::outs() << *t3 << "\n";

  vector<z_cfg_ref_t> cfgs;
  cfgs.push_back(*t1);
  cfgs.push_back(*t2);
  cfgs.push_back(*t3);


  call_graph_t cg (cfgs);
  crab::outs() << cg << "\n";

  /// -- Basic BGL api
  for (auto v : boost::make_iterator_range (vertices (cg))) {
    if (out_degree (v, cg) > 0) {
      crab::outs() << "number of successors " << v.name () << "=" << out_degree (v,cg) << "\n";
      boost::graph_traits<call_graph_t>::out_edge_iterator ei, ei_end;
      boost::tie (ei, ei_end) = out_edges (v, cg);
      for (; ei != ei_end; ++ei) {
        auto s = source (*ei, cg);
        auto t = target (*ei, cg);
        crab::outs() << s.name () << "-->" << t.name ()  << "\n";
      }
    }

    if (in_degree (v, cg) > 0) {
      crab::outs() << "number of predecessors " << v.name () << "=" << in_degree (v,cg) << "\n";
      for (auto e: boost::make_iterator_range (in_edges (v, cg))) {
        auto s = source (e, cg);
        auto t = target (e, cg);
        crab::outs() << s.name () << "-->" << t.name () << "\n";
      }
    }
  }

#if 1
  /// --- DFS
  typedef boost::unordered_map< boost::graph_traits<call_graph_t>::vertex_descriptor, 
                                default_color_type > color_map_t;
  color_map_t color;
  for (auto v : boost::make_iterator_range (vertices (cg))) {
    color[v] = default_color_type();
  }
  boost::associative_property_map< color_map_t > cm (color);

  // find root 
  boost::graph_traits<call_graph_t>::vertex_descriptor root;
  for (auto v: boost::make_iterator_range (vertices (cg))) {
    if (in_degree (v, cg) == 0) {
      root = v;
      break;
    }
  }  
  crab::outs() << "Found root " << root.name () << "\n";
  
  //print all cg nodes in depth-first search order
  crab::outs() << "Printing in preorder ...\n";
  print_visitor vis;
  boost::detail::depth_first_visit_impl (cg, root, vis, cm, detail::nontruth2());
#endif 

#if 1
  // --- SccGraph 
  scc_graph<call_graph_ref_t> scc_g(cg);
  scc_g.write (crab::outs());
  
  std::vector<call_graph_ref_t::node_t> order;
  rev_topo_sort (scc_g, order);
 
  crab::outs() << "reverse topological sort: ";
  for(auto n: order)
    crab::outs() << n.name () << "--";
  crab::outs() << "\n";
  crab::outs() << "topological sort: ";
  for(auto n: boost::make_iterator_range (order.rbegin(), order.rend ()))
    crab::outs() << n.name () << "--";
  crab::outs() << "\n";
  
#endif 

  delete t1;
  delete t2;
  delete t3;

  return 0;
}
