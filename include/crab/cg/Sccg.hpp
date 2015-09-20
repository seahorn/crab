#ifndef STRONGLY_CONNECTED_COMPONENT_GRAPH_HPP__
#define STRONGLY_CONNECTED_COMPONENT_GRAPH_HPP__

#include <boost/graph/strong_components.hpp>

#include <crab/cg/CgBgl.hpp>

///Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

/* 
   Strongly connected component graph
 */

using namespace boost;

namespace crab {
   namespace cg {


      template< typename CG>
      class SccGraph {

        typedef typename CG::node_t cg_node_t;

        /// --- begin internal representation of the SccGraph
        struct  vertex_t { cg_node_t func; };
        typedef adjacency_list<vecS, vecS, bidirectionalS, 
                               property<vertex_color_t, 
                                        default_color_type, 
                                        vertex_t> > scc_graph_t;     
        typedef boost::shared_ptr <scc_graph_t> scc_graph_ptr;
        typedef typename graph_traits<scc_graph_t>::vertex_descriptor vertex_descriptor_t;
        typedef typename graph_traits<scc_graph_t>::edge_descriptor edge_descriptor_t;
        typedef typename graph_traits<scc_graph_t>::vertex_iterator vertex_iterator;
        typedef typename graph_traits<scc_graph_t>::out_edge_iterator out_edge_iterator;
        typedef typename graph_traits<scc_graph_t>::in_edge_iterator in_edge_iterator;
        /// --- end internal representation of the SccGraph

        typedef boost::unordered_map <cg_node_t, vertex_descriptor_t> node_to_vertex_map_t;
        typedef boost::shared_ptr <node_to_vertex_map_t> node_to_vertex_map_ptr;
        typedef boost::unordered_map< cg_node_t, cg_node_t > root_map_t;
        typedef boost::unordered_map <cg_node_t, std::vector<cg_node_t> > m_comp_members_map_t;
        typedef boost::shared_ptr <m_comp_members_map_t> m_comp_members_map_ptr;

        // Wrapper for edges
        // BGL complains if we use std::pair<edge_t,edge_t>
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

        CG m_cg;
        scc_graph_ptr m_sccg;
        node_to_vertex_map_ptr m_node_to_vertex_map;
        root_map_t m_root_map;
        m_comp_members_map_ptr m_comp_members_map;

       public:

        typedef cg_node_t node_t;
        typedef Edge<node_t> edge_t;

       private:

        struct MkNode : 
            public std::unary_function < vertex_descriptor_t, node_t > {

          scc_graph_ptr _g;
          
          MkNode () { }
          MkNode (scc_graph_ptr g): _g (g) { }
          node_t& operator()(const vertex_descriptor_t& v) const { 
            return (*_g)[v].func; 
          }
        };
        
        struct MkEdge :  
            public std::unary_function < edge_descriptor_t, Edge<node_t> > {

          scc_graph_ptr _g;
          
          MkEdge() {}
          MkEdge(scc_graph_ptr g): _g (g) { }
          Edge<node_t> operator()(const edge_descriptor_t& e) const { 
            node_t& s = (*_g) [boost::source (e, *_g)].func;
            node_t& t = (*_g) [boost::target (e, *_g)].func;
            return Edge<node_t> (s,t); 
          }
        };
        
       public:
        typedef boost::transform_iterator<MkNode, vertex_iterator> node_iterator; 
        typedef boost::transform_iterator<MkEdge, in_edge_iterator> pred_iterator; 
        typedef boost::transform_iterator<MkEdge, out_edge_iterator> succ_iterator; 

       public:

        SccGraph (CG cg): 
            m_cg (cg), m_sccg (new scc_graph_t ()), 
            m_node_to_vertex_map (new node_to_vertex_map_t ()),
            m_comp_members_map (new m_comp_members_map_t ()) {

          typedef boost::unordered_map< cg_node_t, std::size_t > component_map_t;
          typedef boost::unordered_map< cg_node_t, default_color_type > color_map_t;
          typedef boost::associative_property_map< component_map_t > property_component_map_t;
          typedef boost::associative_property_map< root_map_t > property_root_map_t;
          typedef boost::associative_property_map< color_map_t > property_color_map_t;
          
          component_map_t comp_map, discover_time;
          color_map_t color_map;

          for (auto const &v : boost::make_iterator_range (vertices (m_cg))) {
            comp_map [v] = 0;
            color_map [v] = default_color_type();
            discover_time [v] = 0;
            m_root_map [v] = cg_node_t ();
          }
          
          boost::strong_components(m_cg,
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
            (*m_sccg) [v].func = p.second;
            m_node_to_vertex_map->insert (make_pair (p.second, v));
            CRAB_DEBUG("Added scc graph node ", p.second.name (), "--- id=", v);
          }

          for (auto p : m_root_map) {
            cg_node_t s = p.second;
            for (auto e:  boost::make_iterator_range (out_edges (s, m_cg))) {
              cg_node_t d = m_root_map [target (e, m_cg)];
              
              if (s == d) {
                CRAB_DEBUG("scc ", s.name (), " is recursive");
                continue; 
              }

              add_edge ((*m_node_to_vertex_map) [s], 
                        (*m_node_to_vertex_map) [d], 
                        *m_sccg);

              CRAB_DEBUG("Added scc graph edge ", 
                         s.name (), 
                         " --> ",
                         d.name ());
            }
          }

          // Build a map from scc roots to their members
          for (auto p: m_root_map) {
            cg_node_t r = p.second;
            auto it  = m_comp_members_map->find (r);
            if (it != m_comp_members_map->end ())
              it->second.push_back (p.first);
            else {
              std::vector<cg_node_t> comp_mems;
              comp_mems.push_back (p.first);
              m_comp_members_map->insert (make_pair (r, comp_mems));
            }
          }

        }

        // return the component root of n
        const cg_node_t& getComponentRoot (const cg_node_t& n) const{
          auto const it = m_root_map.find (n);
          assert (it != m_root_map.end ());
          return it->second;
        }

        // return the members of the scc that contains n
        std::vector<cg_node_t>& getComponentMembers (const cg_node_t& n) {
          const cg_node_t &r = getComponentRoot (n);
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
                o << e.Src ().name () << "--> " << e.Dest ().name () << endl;
              }
            }
          }
          // debugging
          o << "root map: \n";
          for (auto p: m_root_map) {
            o <<"\t" << p.first.name () << " --> " << p.second.name () << endl;
          }
        }
        
      };

   } // end namespace cg
} // end namespace crab

#endif 
