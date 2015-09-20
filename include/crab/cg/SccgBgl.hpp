#ifndef SCC_GRAPH_BOOST_GRAPH_TRAITS_HPP__
#define SCC_GRAPH_BOOST_GRAPH_TRAITS_HPP__

/* Convert a Strongly connected component graph into BGL graph */

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <crab/cg/Sccg.hpp>

namespace boost {

  template<class CG>
  struct graph_traits < crab::cg::SccGraph<CG> >  {

    typedef crab::cg::SccGraph<CG> sccg_t;

    typedef typename sccg_t::node_t vertex_descriptor;
    typedef typename sccg_t::edge_t edge_descriptor;
    typedef typename sccg_t::node_iterator vertex_iterator;
    typedef typename sccg_t::pred_iterator in_edge_iterator;
    typedef typename sccg_t::succ_iterator out_edge_iterator;
    
    typedef disallow_parallel_edge_tag edge_parallel_category;
    typedef bidirectional_tag directed_category;
    struct  this_graph_tag : virtual bidirectional_graph_tag, 
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

}

namespace crab { namespace cg {

   // --- Functions for crab::cg::SccGraph<CG>

   template<class CG> 
   typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor
   source (typename graph_traits< crab::cg::SccGraph<CG> >::edge_descriptor e, 
             const crab::cg::SccGraph<CG> &g) {
     return e.Src(); 
  } 

  template<class CG>
  typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor
  target (typename graph_traits< crab::cg::SccGraph<CG> >::edge_descriptor e, 
          const crab::cg::SccGraph<CG> &g) {
    return e.Dest();
  }

  template<class CG>
  std::pair< typename graph_traits< crab::cg::SccGraph<CG> >::in_edge_iterator, 
             typename graph_traits< crab::cg::SccGraph<CG> >::in_edge_iterator >
  in_edges (typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor v, 
            const crab::cg::SccGraph<CG> &g) {
    return g.preds (v);
  }

  template<class CG>
  typename graph_traits< crab::cg::SccGraph<CG> >::degree_size_type
  in_degree (typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor v, 
             const crab::cg::SccGraph<CG> &g) {
    return g.num_preds (v);
  }

  template<class CG>
  std::pair< typename graph_traits< crab::cg::SccGraph<CG> >::out_edge_iterator, 
             typename graph_traits< crab::cg::SccGraph<CG> >::out_edge_iterator >
  out_edges (typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor v, 
             const crab::cg::SccGraph<CG> &g) {
    return g.succs (v);
  }

  template<class CG>
  typename graph_traits< crab::cg::SccGraph<CG> >::degree_size_type
  out_degree (typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor v, 
              const crab::cg::SccGraph<CG> &g) { 
    return g.num_succs (v);
  }

  template<class CG>
  typename graph_traits< crab::cg::SccGraph<CG> >::degree_size_type
  degree (typename graph_traits< crab::cg::SccGraph<CG> >::vertex_descriptor v, 
          const crab::cg::SccGraph<CG> &g) {
    return g.num_preds (v) + g.num_succs (v);
  }

  template<class CG>
  std::pair<typename graph_traits< crab::cg::SccGraph<CG> >::vertex_iterator, 
            typename graph_traits< crab::cg::SccGraph<CG> >::vertex_iterator > 
  vertices (const crab::cg::SccGraph<CG> &g) {
    return g.nodes ();
  }
  
  template<class CG>
  typename graph_traits< crab::cg::SccGraph<CG> >::vertices_size_type
  num_vertices (const crab::cg::SccGraph<CG> &g) {
    return g.num_nodes ();
  }


} } // end namespaces

#endif 
