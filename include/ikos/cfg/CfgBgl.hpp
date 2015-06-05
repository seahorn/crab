#ifndef CFG_BOOST_GRAPH_TRAITS_HPP
#define CFG_BOOST_GRAPH_TRAITS_HPP

/*
 * Convert a CFG into a BGL (Boost Graph Library) graph
 */
#include <ikos/cfg/Cfg.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/transform_iterator.hpp>

using namespace std;
using namespace boost;

namespace cfg
{  
   namespace graph 
   {

     template <typename G>
     struct MkInEdge: 
      public std::unary_function<typename boost::graph_traits<G>::vertex_descriptor,
                                 typename boost::graph_traits<G>::edge_descriptor>
     {
       typedef typename boost::graph_traits<G>::vertex_descriptor Node;
       typedef typename boost::graph_traits<G>::edge_descriptor   Edge;

       Node _dst;
       MkInEdge () {}
       MkInEdge (const Node &dst) : _dst (dst) {}     
       Edge operator() (const Node &src) const { return Edge (src, _dst); }
     };
   
     template <typename G> 
     struct MkOutEdge: 
         public std::unary_function<typename boost::graph_traits<G>::vertex_descriptor,
                                    typename boost::graph_traits<G>::edge_descriptor>
     {
       typedef typename boost::graph_traits<G>::vertex_descriptor Node;
       typedef typename boost::graph_traits<G>::edge_descriptor   Edge;

       Node _src;
       MkOutEdge () {}
       MkOutEdge (const Node &src) : _src (src) { }
       Edge operator() (const Node &dst) const { return Edge (_src, dst); }
     };


   }  // end namespace graph
} // end namespace cfg

namespace boost 
{
   using namespace std;

   template<typename BasicBlockLabel, typename VariableName>
   struct graph_traits<cfg::Cfg <BasicBlockLabel, VariableName> > 
   {     
     typedef cfg::Cfg<BasicBlockLabel, VariableName> graph_t;
     typedef BasicBlockLabel vertex_descriptor;
     typedef pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
     typedef pair<const vertex_descriptor, 
                  const vertex_descriptor> const_edge_descriptor;
    
     typedef disallow_parallel_edge_tag edge_parallel_category;
     typedef bidirectional_tag directed_category;
     struct  this_graph_tag : virtual bidirectional_graph_tag, 
                              virtual vertex_list_graph_tag {};
     typedef this_graph_tag traversal_category;
     
     typedef size_t vertices_size_type;
     typedef size_t edges_size_type;
     typedef size_t degree_size_type;

     // iterator of BasicBlockLabel's
     typedef typename graph_t::label_iterator vertex_iterator;
     // iterator of BasicBlockLabel's
     typedef boost::transform_iterator<cfg::graph::MkInEdge<graph_t>, 
                                       typename graph_t::pred_iterator> in_edge_iterator;
     // iterator of BasicBlockLabel's
     typedef boost::transform_iterator<cfg::graph::MkOutEdge<graph_t>, 
                                       typename graph_t::succ_iterator> out_edge_iterator; 
   }; // end class graph_traits

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor
   source (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::edge_descriptor e, 
           cfg::Cfg<BasicBlockLabel, VariableName> /*g*/)
   {
     return e.first;
   }

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor
   target (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::edge_descriptor e, 
           cfg::Cfg<BasicBlockLabel,VariableName> /*g*/)
   {
     return e.second;
   }


   template<typename BasicBlockLabel, typename VariableName>
   inline pair<typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_iterator, 
               typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_iterator > 
   vertices (cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     return std::make_pair (g.label_begin (), g.label_end ());
   }
             
   template<typename BasicBlockLabel, typename VariableName>
   inline pair< typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::out_edge_iterator, 
                typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::out_edge_iterator >
   out_edges (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor v, 
              cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     typedef cfg::Cfg<BasicBlockLabel,VariableName> G;

     auto& node = g.get_node (v);
     auto p = node.next_blocks ();
     return std::make_pair (make_transform_iterator (p.first, 
                                                     cfg::graph::MkOutEdge<G> (v)),
                            make_transform_iterator (p.second, 
                                                     cfg::graph::MkOutEdge<G> (v)));     
   }

   template<typename BasicBlockLabel, typename VariableName>
   inline pair< typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::in_edge_iterator, 
                typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::in_edge_iterator >
   in_edges (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor v, 
             cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     typedef cfg::Cfg<BasicBlockLabel,VariableName> G;

     auto& node = g.get_node (v);
     auto p = node.prev_blocks ();
     return std::make_pair (make_transform_iterator (p.first, 
                                                     cfg::graph::MkInEdge<G> (v)),
                            make_transform_iterator (p.second, 
                                                     cfg::graph::MkInEdge<G> (v)));     
   }

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertices_size_type
   num_vertices (cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     return std::distance (g.label_begin (), g.label_end ());
   }

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::degree_size_type
   in_degree (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor v, 
              cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     auto preds = g.prev_nodes (v);
     return std::distance (preds.begin (), preds.end ());
   }

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::degree_size_type
   out_degree (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor v, 
               cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     auto succs = g.next_nodes (v);
     return std::distance (succs.begin (), succs.end ());
   }

   template<typename BasicBlockLabel, typename VariableName>
   typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::degree_size_type
   degree (typename graph_traits<cfg::Cfg<BasicBlockLabel,VariableName> >::vertex_descriptor v, 
           cfg::Cfg<BasicBlockLabel,VariableName> g)
   {
     return out_degree (v, g) + in_degree (v, g);
   }

} // namespace boost

#endif 
