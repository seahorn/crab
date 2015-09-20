#include "../common.hpp"

#include <crab/cg/CgBgl.hpp>
#include <crab/cfg/Cfg.hpp>
#include <crab/cg/Cg.hpp>
#include <crab/cg/Sccg.hpp>
#include <crab/cg/TopoOrder.hpp>


// #include <boost/graph/depth_first_search.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::cfg;
using namespace crab::cg;

using namespace boost;

cfg_t foo (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["x"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["foo"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add (y, x, 1);
  exit.add (z , y , 2);
  exit.ret (vfac ["z"], INT_TYPE);
  return cfg;
}


cfg_t bar (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["a"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["bar"], params);
  // Defining program variables
  z_var a (vfac ["a"]);
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["x"], INT_TYPE));
  exit.callsite (make_pair (vfac["y"], INT_TYPE), vfac ["foo"], args);
  entry.assign (x, a);
  exit.ret (vfac["y"], INT_TYPE);
  return cfg;
}

cfg_t m (VariableFactory &vfac)  {
  vector<pair<varname_t,VariableType> > params;
  FunctionDecl<varname_t> decl (UNK_TYPE, vfac["main"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t cfg ("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& exit   = cfg.insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.assign(x, 3);
  vector<pair<varname_t,VariableType> > args;
  args.push_back (make_pair (vfac["x"], INT_TYPE));
  entry.callsite (make_pair (vfac["y"], INT_TYPE), vfac ["bar"], args);
  exit.add (z, y, 2);
  return cfg;
}

struct print_visitor: public boost::default_dfs_visitor {
  void discover_vertex(graph_traits< CallGraph<cfg_t> >::vertex_descriptor v, 
                       CallGraph<cfg_t> g)  {
    std::cout << v.name () << endl;
  }
};


int main (int argc, char** argv ) {
  VariableFactory vfac;

  cfg_t t1 = foo (vfac);
  cfg_t t2 = bar (vfac);
  cfg_t t3 = m (vfac);

  cout << t1 << endl;
  cout << t2 << endl;
  cout << t3 << endl;

  vector<cfg_t> cfgs;
  cfgs.push_back(t1);
  cfgs.push_back(t2);
  cfgs.push_back(t3);

  CallGraph<cfg_t> cg (cfgs);
  cout << cg << endl;

  /// -- Basic BGL api
  for (auto v : boost::make_iterator_range (vertices (cg))) {
    if (out_degree (v, cg) > 0) {
      cout << "number of successors " << v.name () << "=" << out_degree (v,cg) << endl;
      boost::graph_traits< CallGraph<cfg_t> >::out_edge_iterator ei, ei_end;
      boost::tie (ei, ei_end) = out_edges (v, cg);
      for (; ei != ei_end; ++ei) {
        auto s = source (*ei, cg);
        auto t = target (*ei, cg);
        cout << s.name () << "-->" << t.name ()  << endl;
      }
    }

    if (in_degree (v, cg) > 0) {
      cout << "number of predecessors " << v.name () << "=" << in_degree (v,cg) << endl;
      for (auto e: boost::make_iterator_range (in_edges (v, cg))) {
        auto s = source (e, cg);
        auto t = target (e, cg);
        cout << s.name () << "-->" << t.name () << endl;
      }
    }
  }

#if 1
  /// --- DFS
  typedef boost::unordered_map< boost::graph_traits< CallGraph<cfg_t> >::vertex_descriptor, 
                                default_color_type > color_map_t;
  color_map_t color;
  for (auto v : boost::make_iterator_range (vertices (cg))) {
    color[v] = default_color_type();
  }
  boost::associative_property_map< color_map_t > cm (color);

  // find root 
  boost::graph_traits< CallGraph<cfg_t> >::vertex_descriptor root;
  for (auto v: boost::make_iterator_range (vertices (cg))) {
    if (in_degree (v, cg) == 0) {
      root = v;
      break;
    }
  }  
  cout << "Found root " << root.name () << endl;
  
  //print all cg nodes in depth-first search order
  cout << "Printing in preorder ...\n";
  print_visitor vis;
  boost::detail::depth_first_visit_impl (cg, root, vis, cm, detail::nontruth2());
#endif 

#if 1
  // --- SccGraph 
  SccGraph<CallGraph<cfg_t> > scc_graph(cg);
  scc_graph.write (std::cout);
  
  std::vector<CallGraph<cfg_t>::node_t> order;
  rev_topo_sort (scc_graph, order);
 
  cout << "reverse topological sort: ";
  for(auto n: order)
    cout << n.name () << "--";
  cout << endl;
  cout << "topological sort: ";
  for(auto n: boost::make_iterator_range (order.rbegin(), order.rend ()))
    cout << n.name () << "--";
  cout << endl;
  
#endif 

  return 0;
}
