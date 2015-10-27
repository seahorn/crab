#ifndef STRONGLY_CONNECTED_COMPONENT_GRAPH_HPP__
#define STRONGLY_CONNECTED_COMPONENT_GRAPH_HPP__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/strong_components.hpp>

#include <crab/cg/Cg.hpp> // for graph_algo_impl namespace

///Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

/* 
   Strongly connected component graph
 */

using namespace boost;

namespace crab {
   namespace analyzer {
     namespace graph_algo {

      namespace graph_algo_impl {
        // -- convert a scc graph node to a string label. 
        template<typename Node>
        string str (const Node& n) {
          return crab::cfg_impl::get_label_str (n);
        }
        template<typename G> 
        string str (const typename crab::cg::CallGraph<G>::CgNode& n) { 
          return n.name ();
        }
        template<> string str (const string& n) { return n;}
      } // end namespace

      template< typename G>
      class SccGraph {

       public:

        typedef typename G::node_t node_t;

       private:

        /// --- begin internal representation of the SccGraph
        struct  vertex_t { node_t m_node; };
        typedef adjacency_list<vecS, vecS, bidirectionalS, 
                               property<vertex_color_t, 
                                        default_color_type, 
                                        vertex_t> > scc_graph_t;     
        typedef boost::shared_ptr <scc_graph_t> scc_graph_ptr;
        typedef typename boost::graph_traits<scc_graph_t>::vertex_descriptor vertex_descriptor_t;
        typedef typename boost::graph_traits<scc_graph_t>::edge_descriptor edge_descriptor_t;
        typedef typename boost::graph_traits<scc_graph_t>::vertex_iterator vertex_iterator;
        typedef typename boost::graph_traits<scc_graph_t>::out_edge_iterator out_edge_iterator;
        typedef typename boost::graph_traits<scc_graph_t>::in_edge_iterator in_edge_iterator;
        /// --- end internal representation of the SccGraph

        typedef boost::unordered_map <node_t, vertex_descriptor_t> node_to_vertex_map_t;
        typedef boost::shared_ptr <node_to_vertex_map_t> node_to_vertex_map_ptr;
        typedef boost::unordered_map< node_t, node_t > root_map_t;
        typedef boost::unordered_map <node_t, std::vector<node_t> > m_comp_members_map_t;
        typedef boost::shared_ptr <m_comp_members_map_t> m_comp_members_map_ptr;

        // Wrapper for edges
        // XXX: BGL complains if we use std::pair<edge_t,edge_t>
        template <typename T>
        struct Edge  {
          T m_s;
          T m_d;
          Edge (){ }
          Edge (T s, T d): m_s (s), m_d (d) { }
          T Src () const { return m_s; }
          T Dest () const { return m_d; }
          bool operator== (const Edge<T> &o) const {
            return (m_s == o.Src () && m_d == o.Dest ());
          }
          bool operator!= (const Edge<T> &o) const {
            return !(*this == o);
          }
        };

        G m_g;
        scc_graph_ptr m_sccg;
        node_to_vertex_map_ptr m_node_to_vertex_map;
        root_map_t m_root_map;
        m_comp_members_map_ptr m_comp_members_map;

       public:

        typedef Edge<node_t> edge_t;

       private:

        struct MkNode : 
            public std::unary_function < vertex_descriptor_t, node_t > {

          scc_graph_ptr _g;
          
          MkNode () { }
          MkNode (scc_graph_ptr g): _g (g) { }
          node_t& operator()(const vertex_descriptor_t& v) const { 
            return (*_g)[v].m_node; 
          }
        };
        
        struct MkEdge :  
            public std::unary_function < edge_descriptor_t, Edge<node_t> > {

          scc_graph_ptr _g;
          
          MkEdge() {}
          MkEdge(scc_graph_ptr g): _g (g) { }
          Edge<node_t> operator()(const edge_descriptor_t& e) const { 
            node_t& s = (*_g) [boost::source (e, *_g)].m_node;
            node_t& t = (*_g) [boost::target (e, *_g)].m_node;
            return Edge<node_t> (s,t); 
          }
        };
        
       public:
        typedef boost::transform_iterator<MkNode, vertex_iterator> node_iterator; 
        typedef boost::transform_iterator<MkEdge, in_edge_iterator> pred_iterator; 
        typedef boost::transform_iterator<MkEdge, out_edge_iterator> succ_iterator; 

       public:

        SccGraph (G g): 
            m_g (g), m_sccg (new scc_graph_t ()), 
            m_node_to_vertex_map (new node_to_vertex_map_t ()),
            m_comp_members_map (new m_comp_members_map_t ()) {

          typedef boost::unordered_map< node_t, std::size_t > component_map_t;
          typedef boost::unordered_map< node_t, default_color_type > color_map_t;
          typedef boost::associative_property_map< component_map_t > property_component_map_t;
          typedef boost::associative_property_map< root_map_t > property_root_map_t;
          typedef boost::associative_property_map< color_map_t > property_color_map_t;
          
          component_map_t comp_map, discover_time;
          color_map_t color_map;

          for (auto const &v : boost::make_iterator_range (vertices (m_g))) {
            comp_map [v] = 0;
            color_map [v] = default_color_type();
            discover_time [v] = 0;
            m_root_map [v] = node_t ();
          }
          
          boost::strong_components(m_g,
                                   property_component_map_t(comp_map),
                                   root_map(property_root_map_t(m_root_map))
                                   .color_map(property_color_map_t(color_map))
                                   .discover_time_map(property_component_map_t(discover_time)));
          

          // build scc graph
          for (auto p : m_root_map) {
            auto it = m_node_to_vertex_map->find (p.second );
            if (it != m_node_to_vertex_map->end ())
              continue;
            vertex_descriptor_t v = add_vertex (*m_sccg);
            (*m_sccg) [v].m_node = p.second;
            m_node_to_vertex_map->insert (std::make_pair (p.second, v));
            CRAB_DEBUG("Added scc graph node ", 
                       graph_algo_impl::str (p.second), 
                       "--- id=", v);
          }
          
          // Build the DAG
          for (const node_t &u : boost::make_iterator_range (vertices (m_g))) {
            for (auto e : boost::make_iterator_range (out_edges (u, m_g))) {
              const node_t &d = target (e, m_g);              

              if (m_root_map [u] == m_root_map [d])
                continue;

              add_edge ((*m_node_to_vertex_map) [m_root_map [u]], 
                        (*m_node_to_vertex_map) [m_root_map [d]], 
                        *m_sccg);

              CRAB_DEBUG("Added scc graph edge ", 
                         graph_algo_impl::str (m_root_map [u]), 
                         " --> ",
                         graph_algo_impl::str (m_root_map [d]));

            }
          }

          // Build a map from scc roots to their members
          for (auto p: m_root_map) {
            node_t r = p.second;
            auto it  = m_comp_members_map->find (r);
            if (it != m_comp_members_map->end ())
              it->second.push_back (p.first);
            else {
              std::vector<node_t> comp_mems;
              comp_mems.push_back (p.first);
              m_comp_members_map->insert (std::make_pair (r, comp_mems));
            }
          }

        }

        // return the component root of n
        const node_t& getComponentRoot (const node_t& n) const{
          auto const it = m_root_map.find (n);
          assert (it != m_root_map.end ());
          return it->second;
        }

        // return the members of the scc that contains n
        std::vector<node_t>& getComponentMembers (const node_t& n) {
          const node_t &r = getComponentRoot (n);
          return (*m_comp_members_map) [r];
        }

        std::pair<node_iterator, node_iterator> 
        nodes () const {
          auto p = boost::vertices (*m_sccg);
          return std::make_pair (make_transform_iterator (p.first, MkNode (m_sccg)),
                                 make_transform_iterator (p.second, MkNode (m_sccg)));
        }


        std::pair<succ_iterator, succ_iterator> 
        succs (const node_t &n) const {
          vertex_descriptor_t v = (*m_node_to_vertex_map) [n]; 
          auto p = boost::out_edges (v, *m_sccg);
          return std::make_pair (make_transform_iterator (p.first, MkEdge (m_sccg)),
                                 make_transform_iterator (p.second, MkEdge (m_sccg)));
        }

        std::pair<pred_iterator, pred_iterator> 
        preds (const node_t &n) const {
          vertex_descriptor_t v = (*m_node_to_vertex_map) [n]; 
          auto p = boost::in_edges (v, *m_sccg);
          return std::make_pair (make_transform_iterator (p.first, MkEdge (m_sccg)),
                                 make_transform_iterator (p.second, MkEdge (m_sccg)));
        }

        std::size_t num_nodes () const {
           return boost::num_vertices (*m_sccg);
        }

        std::size_t num_succs (const node_t &n) const {
          vertex_descriptor_t v = (*m_node_to_vertex_map) [n]; 
          return boost::out_degree (v, *m_sccg);
        }

        std::size_t num_preds (const node_t &n) const {
          vertex_descriptor_t v = (*m_node_to_vertex_map) [n]; 
          return boost::in_degree (v, *m_sccg);
        }

        void write(std::ostream& o) const {
          o << "SCCG=\n";
          for (auto f: boost::make_iterator_range (nodes ())){
            if (num_succs (f) > 0) {
              for (auto e: boost::make_iterator_range (succs (f)))  {
                o << graph_algo_impl::str (e.Src ()) 
                  << "--> " 
                  << graph_algo_impl::str (e.Dest ()) << std::endl;
              }
            }
          }
          // debugging
          o << "root map: \n";
          for (auto p: m_root_map) {
            o <<"\t" 
              << graph_algo_impl::str (p.first) 
              << " --> " 
              << graph_algo_impl::str (p.second) << std::endl;
          }
        }
        
      };
     } // end namespace graph_algo
   } // end namespace analyzer
} // end namespace crab

#endif 
