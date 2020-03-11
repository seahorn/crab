#pragma once

/* Convert a call_graph into BGL graph */

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <crab/cg/cg.hpp>

namespace boost {

template <class CFG> struct graph_traits<crab::cg::call_graph<CFG>> {

  typedef crab::cg::call_graph<CFG> cg_t;

  typedef typename cg_t::node_t vertex_descriptor;
  typedef typename cg_t::edge_t edge_descriptor;
  typedef typename cg_t::node_iterator vertex_iterator;
  typedef typename cg_t::pred_iterator in_edge_iterator;
  typedef typename cg_t::succ_iterator out_edge_iterator;

  typedef disallow_parallel_edge_tag edge_parallel_category;
  typedef bidirectional_tag directed_category;
  struct this_graph_tag : virtual bidirectional_graph_tag,
                          virtual vertex_list_graph_tag {};
  typedef this_graph_tag traversal_category;

  typedef size_t vertices_size_type;
  typedef size_t edges_size_type;
  typedef size_t degree_size_type;

  static vertex_descriptor null_vertex() {
    if (std::is_pointer<vertex_descriptor>::value)
      return nullptr;
    else {
      // XXX: if vertex_descriptor is a basic type then
      // null_vertex will return an undefined value, otherwise it
      // will return the result of calling the default
      // constructor.
      vertex_descriptor n;
      return n;
    }
  }
}; // end class graph_traits

template <class CG> struct graph_traits<crab::cg::call_graph_ref<CG>> {
  typedef crab::cg::call_graph_ref<CG> cg_ref_t;
  typedef typename cg_ref_t::node_t vertex_descriptor;
  typedef typename cg_ref_t::edge_t edge_descriptor;
  typedef typename cg_ref_t::node_iterator vertex_iterator;
  typedef typename cg_ref_t::pred_iterator in_edge_iterator;
  typedef typename cg_ref_t::succ_iterator out_edge_iterator;

  typedef disallow_parallel_edge_tag edge_parallel_category;
  typedef bidirectional_tag directed_category;
  struct this_graph_tag : virtual bidirectional_graph_tag,
                          virtual vertex_list_graph_tag {};
  typedef this_graph_tag traversal_category;

  typedef size_t vertices_size_type;
  typedef size_t edges_size_type;
  typedef size_t degree_size_type;

  static vertex_descriptor null_vertex() {
    vertex_descriptor n;
    return n;
  }
}; // end class graph_traits

} // end namespace boost

// XXX: do not put it in the boost namespace because it won't compile
namespace crab {
namespace cg {

// --- Functions for call_graph

// this is not part of BGL but needed by Crab's algorithms
template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor
entry(const call_graph<CFG> &g) {
  return g.entry();
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor
source(typename boost::graph_traits<call_graph<CFG>>::edge_descriptor e,
       const call_graph<CFG> &g) {
  return e.src();
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor
target(typename boost::graph_traits<call_graph<CFG>>::edge_descriptor e,
       const call_graph<CFG> &g) {
  return e.dest();
}

template <class CFG>
std::pair<typename boost::graph_traits<call_graph<CFG>>::in_edge_iterator,
          typename boost::graph_traits<call_graph<CFG>>::in_edge_iterator>
in_edges(typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor v,
         const call_graph<CFG> &g) {
  return g.preds(v);
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::degree_size_type
in_degree(typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor v,
          const call_graph<CFG> &g) {
  return g.num_preds(v);
}

template <class CFG>
std::pair<typename boost::graph_traits<call_graph<CFG>>::out_edge_iterator,
          typename boost::graph_traits<call_graph<CFG>>::out_edge_iterator>
out_edges(typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor v,
          const call_graph<CFG> &g) {
  return g.succs(v);
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::degree_size_type
out_degree(typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor v,
           const call_graph<CFG> &g) {
  return g.num_succs(v);
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::degree_size_type
degree(typename boost::graph_traits<call_graph<CFG>>::vertex_descriptor v,
       const call_graph<CFG> &g) {
  return g.num_preds(v) + g.num_succs(v);
}

template <class CFG>
std::pair<typename boost::graph_traits<call_graph<CFG>>::vertex_iterator,
          typename boost::graph_traits<call_graph<CFG>>::vertex_iterator>
vertices(const call_graph<CFG> &g) {
  return g.nodes();
}

template <class CFG>
typename boost::graph_traits<call_graph<CFG>>::vertices_size_type
num_vertices(const call_graph<CFG> &g) {
  return g.num_nodes();
}

// --- Functions for call_graph_ref

// this is not part of BGL but needed by Crab's algorithms
template <class CFG>
typename boost::graph_traits<call_graph_ref<CFG>>::vertex_descriptor
entry(const call_graph_ref<CFG> &g) {
  return g.entry();
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor
source(typename boost::graph_traits<call_graph_ref<CG>>::edge_descriptor e,
       const call_graph_ref<CG> &g) {
  return e.src();
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor
target(typename boost::graph_traits<call_graph_ref<CG>>::edge_descriptor e,
       const call_graph_ref<CG> &g) {
  return e.dest();
}

template <class CG>
std::pair<typename boost::graph_traits<call_graph_ref<CG>>::in_edge_iterator,
          typename boost::graph_traits<call_graph_ref<CG>>::in_edge_iterator>
in_edges(typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor v,
         const call_graph_ref<CG> &g) {
  return g.preds(v);
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::degree_size_type
in_degree(typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor v,
          const call_graph_ref<CG> &g) {
  return g.num_preds(v);
}

template <class CG>
std::pair<typename boost::graph_traits<call_graph_ref<CG>>::out_edge_iterator,
          typename boost::graph_traits<call_graph_ref<CG>>::out_edge_iterator>
out_edges(typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor v,
          const call_graph_ref<CG> &g) {
  return g.succs(v);
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::degree_size_type out_degree(
    typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor v,
    const call_graph_ref<CG> &g) {
  return g.num_succs(v);
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::degree_size_type
degree(typename boost::graph_traits<call_graph_ref<CG>>::vertex_descriptor v,
       const call_graph_ref<CG> &g) {
  return g.num_preds(v) + g.num_succs(v);
}

template <class CG>
std::pair<typename boost::graph_traits<call_graph_ref<CG>>::vertex_iterator,
          typename boost::graph_traits<call_graph_ref<CG>>::vertex_iterator>
vertices(const call_graph_ref<CG> &g) {
  return g.nodes();
}

template <class CG>
typename boost::graph_traits<call_graph_ref<CG>>::vertices_size_type
num_vertices(const call_graph_ref<CG> &g) {
  return g.num_nodes();
}

} // namespace cg
} // namespace crab
