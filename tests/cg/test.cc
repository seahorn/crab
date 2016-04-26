#include "../common.hpp"

#include <crab/cg/Cg.hpp>
#include <crab/cg/CgBgl.hpp>
#include <crab/analysis/graphs/Sccg.hpp>
#include <crab/analysis/graphs/TopoOrder.hpp>


using namespace std;
using namespace crab::analyzer;
using namespace crab::analyzer::graph_algo;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::cg;
using namespace crab::domain_impl;


using namespace boost;

typedef CallGraph<cfg_ref_t> call_graph_t;
typedef CallGraph_Ref<call_graph_t> call_graph_ref_t;

cfg_t* foo (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["x"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["foo"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
  // adding control flow
  entry >> exit;
  // adding statements
  entry.add (y, x, 1);
  exit.add (z , y , 2);
  exit.ret (vfac ["z"], INT_TYPE);
  return cfg;
}


cfg_t* bar (VariableFactory &vfac) {
  vector<pair<varname_t,VariableType> > params;
  params.push_back (make_pair (vfac["a"], INT_TYPE));
  FunctionDecl<varname_t> decl (INT_TYPE, vfac["bar"], params);
  // Defining program variables
  z_var a (vfac ["a"]);
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
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

cfg_t* m (VariableFactory &vfac)  {
  vector<pair<varname_t,VariableType> > params;
  FunctionDecl<varname_t> decl (UNK_TYPE, vfac["main"], params);
  // Defining program variables
  z_var x (vfac ["x"]);
  z_var y (vfac ["y"]);
  z_var z (vfac ["z"]);
  // entry and exit block
  cfg_t* cfg = new cfg_t("entry", "exit", decl);
  // adding blocks
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& exit   = cfg->insert ("exit");
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
  void discover_vertex(graph_traits<call_graph_t>::vertex_descriptor v, 
                       const call_graph_t& g)  {
    crab::outs() << v.name () << endl;
  }
};


int main (int argc, char** argv ) {
  SET_LOGGER(argc,argv)

  VariableFactory vfac;

  cfg_t* t1 = foo (vfac);
  cfg_t* t2 = bar (vfac);
  cfg_t* t3 = m (vfac);

  crab::outs() << *t1 << endl;
  crab::outs() << *t2 << endl;
  crab::outs() << *t3 << endl;

  vector<cfg_ref_t> cfgs;
  cfgs.push_back(*t1);
  cfgs.push_back(*t2);
  cfgs.push_back(*t3);


  call_graph_t cg (cfgs);
  crab::outs() << cg << endl;

  /// -- Basic BGL api
  for (auto v : boost::make_iterator_range (vertices (cg))) {
    if (out_degree (v, cg) > 0) {
      crab::outs() << "number of successors " << v.name () << "=" << out_degree (v,cg) << endl;
      boost::graph_traits<call_graph_t>::out_edge_iterator ei, ei_end;
      boost::tie (ei, ei_end) = out_edges (v, cg);
      for (; ei != ei_end; ++ei) {
        auto s = source (*ei, cg);
        auto t = target (*ei, cg);
        crab::outs() << s.name () << "-->" << t.name ()  << endl;
      }
    }

    if (in_degree (v, cg) > 0) {
      crab::outs() << "number of predecessors " << v.name () << "=" << in_degree (v,cg) << endl;
      for (auto e: boost::make_iterator_range (in_edges (v, cg))) {
        auto s = source (e, cg);
        auto t = target (e, cg);
        crab::outs() << s.name () << "-->" << t.name () << endl;
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
  crab::outs() << "Found root " << root.name () << endl;
  
  //print all cg nodes in depth-first search order
  crab::outs() << "Printing in preorder ...\n";
  print_visitor vis;
  boost::detail::depth_first_visit_impl (cg, root, vis, cm, detail::nontruth2());
#endif 

#if 1
  // --- SccGraph 
  SccGraph<call_graph_ref_t> scc_graph(cg);
  scc_graph.write (crab::outs());
  
  std::vector<call_graph_ref_t::node_t> order;
  rev_topo_sort (scc_graph, order);
 
  crab::outs() << "reverse topological sort: ";
  for(auto n: order)
    crab::outs() << n.name () << "--";
  crab::outs() << endl;
  crab::outs() << "topological sort: ";
  for(auto n: boost::make_iterator_range (order.rbegin(), order.rend ()))
    crab::outs() << n.name () << "--";
  crab::outs() << endl;
  
#endif 

  delete t1;
  delete t2;
  delete t3;

  return 0;
}
