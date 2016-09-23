/*******************************************************************************
 * This domain is a simplified version of the paper: 
 * 
 * "A Partial-Order Approach to Array Content Analysis" by 
 *  Gange, Navas, Schachte, Sondergaard, and Stuckey
 *  available here http://arxiv.org/pdf/1408.1754v1.pdf.
 *
 * It reasons about array contents and the idea is to compute all the
 * feasible partial orderings between array indexes. It keeps a single
 * graph where vertices are the array indexes (potentially all scalar
 * variables) and edges are labelled with weights (that includes
 * scalar and array variables).  An edge (i,j) with weight w denotes
 * that the property w holds for the all elements in the array between
 * [i,j).
 ******************************************************************************/

/* Limitations:

  - The implementation is just a proof-of-concept so it is horribly
    inefficient. I have not tried to make it more efficient yet or
    ran even with real programs.
  - Assume all array accesses are aligned wrt to the size of the array
    element (e.g., if the size of the array element is 4 bytes then
    all array accesses must be multiple of 4).
  - Assume that the size of the array element is always
    1. Therefore, if the array indexes are incremented or decremented
    by 2,4,... we will lose all the precision.

  FIXMEs:

  - Use AdaptGraph instead of boost::graph
  - Have a flag is_normalized and normalize only if the flag is false
  - Perform incremental Floyd-Warshall by keeping track of changed edges.
  - array_graph::widening is normalizing both operands. The first
    operand cannot be normalized.
  - perform common renaming before binary operations (join/widening/narrowing/meet)
  - Do no use shared_ptr because we are doing deep copies Use magic
    move semantics and also shared references with copy-on-write to
    avoid unnecessary copies.
  - etc.

 */

#ifndef ARRAY_GRAPH_HPP
#define ARRAY_GRAPH_HPP

#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_traits.hpp> 
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/lexical_cast.hpp>

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/domains/operators_api.hpp>
#include <crab/domains/domain_traits.hpp>

using namespace std;
using namespace boost;
using namespace ikos;

namespace crab {

  namespace domains {

     /*
       A weighted array graph is a graph (V,E,L) where V is the set of
       vertices, E the edges and L is a label function: E -> W such that:
       - If there is an edge e from vertex i to j whose weight is not
       bottom (i.e., L(e) != bot) then it means that i < j.
       - It is possible one non-bottom edge from i to j and another
       non-bottom one from j to i. This means that both i < j and j < i
       are possible.
       - If the both edges from i to j and from j to i are bottom then it
       must be that (i>=j) && (j>=i) = (i==j)
     */
     template< typename VertexName, typename Weight, 
               typename ScalarNumDomain, bool IsDistWeight >
     class array_graph: public writeable{
       
       typedef VertexName VertexNameKey;
       
       template < typename Any1, typename Any2, bool Any5> 
       friend class array_graph_domain;
       
       typedef index_t key_t;
       typedef boost::shared_ptr<VertexName> VertexNamePtr;
       typedef boost::shared_ptr<Weight>     WeightPtr;
       struct  graph_vertex_t { VertexNamePtr name; };
       struct  graph_edge_t   { WeightPtr weight; }; 
       
      public:
       typedef array_graph<VertexName,Weight,ScalarNumDomain,IsDistWeight> array_graph_t;
       typedef boost::tuple<VertexName, VertexName, Weight> edge_t;
       
      private:
       typedef adjacency_list<listS,listS,bidirectionalS,graph_vertex_t,graph_edge_t> graph_t;
       
       typedef typename graph_traits<graph_t>::edge_iterator     edge_iterator;
       typedef typename graph_traits<graph_t>::vertex_iterator   vertex_iterator;
       typedef typename graph_traits<graph_t>::edge_descriptor   edge_descriptor_t;
       typedef typename graph_traits<graph_t>::vertex_descriptor vertex_descriptor_t;
       typedef typename graph_traits<graph_t>::out_edge_iterator out_edge_iterator;  
       typedef typename graph_traits<graph_t>::in_edge_iterator  in_edge_iterator;  
       
       typedef typename ScalarNumDomain::linear_constraint_t     linear_constraint_t;
       typedef typename ScalarNumDomain::variable_t              variable_t;
       
       typedef boost::unordered_map<key_t, vertex_descriptor_t >  vertex_map_t;
       
       typedef std::set<VertexName> vertex_names_set_t;
       
       typedef boost::shared_ptr <graph_t> graph_ptr;
       typedef boost::shared_ptr< vertex_map_t > vertex_map_ptr;
       typedef boost::shared_ptr< vertex_names_set_t > vertex_names_set_ptr;
       
       bool _is_bottom;
       graph_ptr  _graph; 
       vertex_map_ptr _vertex_map;   //! map a VertexName to a graph vertex
       vertex_names_set_ptr _vertices_set;
       
       bool find_vertex_map (VertexName v) {
         return (_vertex_map->find(v.index()) != _vertex_map->end()); 
       }
       
       void insert_vertex_map (VertexName key, vertex_descriptor_t value)
       {
         if (find_vertex_map(key))
           CRAB_ERROR (key," already in the vertex map");
         
         _vertex_map->insert (make_pair(key.index(), value));
         _vertices_set->insert (key);
       }
       
       void remove_vertex_map (VertexName key)
       {
         _vertex_map->erase (key.index());
         _vertices_set->erase (key);
       }
    
       vertex_descriptor_t lookup_vertex_map (VertexName key) const
       {
         auto it = _vertex_map->find(key.index());
         if (it != _vertex_map->end())
           return it->second;
      
         CRAB_ERROR ("No vertex with name ",key," found in the graph");
       }
    
       // All methods that add new vertices should call this one.
       void add (vector<VertexName> vertices, vector<edge_t> edges) 
       {
         for(auto v: vertices)
         {
           vertex_descriptor_t u = add_vertex(*_graph);
           (*_graph)[u].name = VertexNamePtr (new VertexName(v));
           insert_vertex_map(v, u);
         }
      
         for(auto e: edges)
         {
           vertex_descriptor_t u =  lookup_vertex_map(e.template get<0>());
           vertex_descriptor_t v =  lookup_vertex_map(e.template get<1>());
           edge_descriptor_t k; bool b;
           boost::tie(k,b) = add_edge(u, v, *_graph);
           if (!b)
             CRAB_ERROR ("edge is already in the graph");
        
           (*_graph)[k].weight = WeightPtr(new Weight(e.template get<2>()));
         }
      
         canonical();
       }
    
       // All methods that remove vertices should call this one.
       void remove (VertexName v)
       {
         if (!find_vertex_map (v)) return ;
      
         canonical(); 
         vertex_descriptor_t u = lookup_vertex_map(v);
      
         // remove all in and out edges to/from u
         clear_vertex(u, *_graph);
      
         // remove the vertex
         remove_vertex(u, *_graph);
         remove_vertex_map(v);
       }
    
       ///////////////////////////////////////////////////////////////////////
       // For our canonical form, we would like to compute the greatest
       // fixed point to the set of inequalities:
       //    \forall i,j,k. G[i,j] \subseteq G[i,k] \cup G[k,j]
       // If the weight domain is distributive we can solve this set of
       // inequations by solving:
       //    \forall i,j,k. G[i,j] = G[i,j] \cap G[i,k] \cup G[k,j]
       // The Floyd-Warshall algorithm does exactly that.
       // Otherwise, we iterate the Floyd-Warshall algorithm until no change.
       ///////////////////////////////////////////////////////////////////////
       bool oneStep()
       {
         binary_join_op join;
         binary_meet_op meet;
         vertex_iterator ki, ke, ii, ie, ji, je;
         bool change = false;
         for (tie(ki, ke) = vertices(*_graph); ki != ke; ++ki)
         {
           for (tie(ii, ie) = vertices(*_graph); ii != ie; ++ii)
           {
             if (edge(*ii, *ki, *_graph).second)
             {
               for (tie(ji, je) = vertices(*_graph); ji != je; ++ji)
               {
                 if (edge(*ii, *ji, (*_graph)).second && 
                     edge(*ki, *ji, (*_graph)).second) 
                 {
                   auto e_ij = edge(*ii, *ji, (*_graph)).first;
                   auto e_ik = edge(*ii, *ki, (*_graph)).first;
                   auto e_kj = edge(*ki, *ji, (*_graph)).first;
                   auto Old = (*_graph)[e_ij].weight;
                   auto New = meet((*_graph)[e_ij].weight,
                                   (join((*_graph)[e_ik].weight,
                                         (*_graph)[e_kj].weight)));
                   change |= (!(*Old <= *New && *New <= *Old));
                   (*_graph)[e_ij].weight = New;
                 }
               }
             }
           }
         }
         return change;
       }
    
       void canonical()
       {
         if (IsDistWeight)
           oneStep ();
         else
         {
           bool change = true;
           while (change)
             change = oneStep ();
         }
       }
    
       void insert_vertex (VertexName u, Weight val = Weight::top())
       {
         if (!is_bottom () && !find_vertex_map(u))
         {
           vector<VertexName> new_vertices;
           new_vertices.push_back(u);
           vector<edge_t> new_edges;
           vertex_iterator i, e;
           for (tie(i, e) = vertices(*_graph); i != e; ++i)
           {
             VertexNamePtr v = (*_graph)[*i].name;
             // add two edges in both directions
             new_edges.push_back(edge_t(u, *v, val));
             new_edges.push_back(edge_t(*v, u, val));
           }
           add(new_vertices, new_edges);
         }
       }
    
       template<typename Iterator>
       void insert_vertices(array_graph_t &g, Iterator begin, Iterator end)
       {
         for(;begin!=end;++begin)
           g.insert_vertex(begin->name());
       }  
    
       // pre: caller must ensure the graph is in canonical form
       void set_incoming (const VertexName &v, const Weight &weight)
       {
         if (!is_bottom ())
         {      
           vertex_descriptor_t u = lookup_vertex_map(v);
           in_edge_iterator in_it, in_et;
           for (tie(in_it, in_et) = in_edges(u, *_graph); in_it != in_et; ++in_it)
             (*_graph)[*in_it].weight = WeightPtr( new Weight(weight));
         }
       }
    
       // pre: caller must ensure the graph is in canonical form
       void set_outgoing (const VertexName &v, const Weight &weight)
       {
         if (!is_bottom ())
         {
           vertex_descriptor_t u = lookup_vertex_map(v);
           out_edge_iterator out_it, out_et;
           for (tie(out_it, out_et) = out_edges(u, *_graph); 
                out_it != out_et; ++out_it)
             (*_graph)[*out_it].weight = WeightPtr( new Weight(weight));
         }
       }
    
       struct binary_join_op{
         WeightPtr operator()(WeightPtr w1, WeightPtr w2) {
           return WeightPtr(new Weight(*w1 | *w2));
         }
       };
    
       struct binary_meet_op{
         WeightPtr operator()(WeightPtr w1, WeightPtr w2) {
           return WeightPtr(new Weight(*w1 & *w2));
         }
       };
    
       struct binary_widening_op{
         WeightPtr operator()(WeightPtr w1, WeightPtr w2) {
           return WeightPtr(new Weight(*w1 || *w2));
         }
       };
    
       struct binary_narrowing_op{
         WeightPtr operator()(WeightPtr w1, WeightPtr w2) {
           return WeightPtr(new Weight(*w1 && *w2));
         }
       };
    
    
       template<typename Op>
       void pointwise_binop_helper (array_graph_t &g1, 
                                    const array_graph_t &g2)
       {
         // pre: g1 and g2 have the same adjacency structure
         edge_iterator it_1, et_1;
         for(tie(it_1,et_1) = edges(*g1._graph); it_1 != et_1; ++it_1)
         {
           edge_descriptor_t e_1   = *it_1;
           vertex_descriptor_t u_1 = source(e_1, *g1._graph);
           vertex_descriptor_t v_1 = target(e_1, *g1._graph);
           VertexNamePtr u_name_1  = (*g1._graph)[u_1].name;
           VertexNamePtr v_name_1  = (*g1._graph)[v_1].name;
           WeightPtr     weight_1  = (*g1._graph)[e_1].weight;
           vertex_descriptor_t u_2 = g2.lookup_vertex_map(*u_name_1);
           vertex_descriptor_t v_2 = g2.lookup_vertex_map(*v_name_1);
           if (edge(u_2, v_2, *g2._graph).second)
           {
             Op op;
             edge_descriptor_t e_2 = edge(u_2, v_2, *g2._graph).first;
             (*g1._graph)[e_1].weight = 
                 op((*g1._graph)[e_1].weight, (*g2._graph)[e_2].weight);
           }
           else
             CRAB_ERROR("unreachable");
         } 
       }
    
       template<typename Op>
       array_graph_t pointwise_binop (array_graph_t g1, 
                                      array_graph_t g2)
       {
         g1.canonical();
         g2.canonical();
      
         // if (*(g1._vertices_set) != *(g2._vertices_set)){
         //   set<VertexName> all_vs;
         //   set_union(g1._vertices_set->begin(), g1._vertices_set->end(), 
         //             g2._vertices_set->begin(), g2._vertices_set->end(), inserter(all_vs, all_vs.end()));
         //   vector<VertexName> new_g1, new_g2;
         //   set_difference(all_vs.begin(), all_vs.end(),
         //                  g1._vertices_set->begin(), g1._vertices_set->end(), inserter(new_g1, new_g1.end()));
         //   set_difference(all_vs.begin(), all_vs.end(),
         //                  g2._vertices_set->begin(), g2._vertices_set->end(), inserter(new_g2, new_g2.end()));
         //   insert_vertices<typename vector<VertexName>::iterator>(g1, new_g1.begin(), new_g1.end());
         //   insert_vertices<typename vector<VertexName>::iterator>(g2, new_g2.begin(), new_g2.end());
         // }
      
         // pre: g1 and g2 have the same set of vertices and edges at this
         // point
         pointwise_binop_helper<Op>(g1,g2);
         return g1;
       }
    
       array_graph(bool is_bot): 
           _is_bottom(is_bot), 
           _graph(new graph_t(0)), 
           _vertex_map(new vertex_map_t()), 
           _vertices_set(new vertex_names_set_t()) 
       {  }
    
      public:
    
       static array_graph_t bottom()  { return array_graph(true); }
    
       static array_graph_t top() { return array_graph(false); }
    
       // Deep copy of the array graph
       array_graph(const array_graph_t &other): 
           writeable(), 
           _is_bottom(other._is_bottom), 
           //_graph(new graph_t(*other._graph)),
           //_vertex_map(new vertex_map_t(*other._vertex_map)),
           //_vertices_set( new vertex_names_set_t(*other._vertices_set))
           _graph (new graph_t(0)),
           _vertex_map (new vertex_map_t()), 
           _vertices_set (new vertex_names_set_t())
       {
         crab::CrabStats::count ("Domain.count.copy");
         crab::ScopedCrabStats __st__("Domain.copy");

         if (!is_bottom ())
         {
           // copy vertices, _vertex_map and _vertices_set
           vertex_iterator i, e;
           for (tie(i, e) = vertices(*other._graph); i != e; ++i)
           {
             vertex_descriptor_t u  = add_vertex(*_graph);
             VertexName u_name = *((*other._graph)[*i].name);
             (*_graph)[u].name = VertexNamePtr(new VertexName(u_name));
             insert_vertex_map(u_name, u);
           }
        
           // copy edges
           edge_iterator ie, ee;
           for(tie(ie,ee) = edges(*other._graph); ie != ee; ++ie)
           {
             VertexName u = *((*other._graph)[source(*ie, *other._graph)].name);
             VertexName v = *((*other._graph)[target(*ie, *other._graph)].name);
             Weight     w = *((*other._graph)[*ie].weight);
             vertex_descriptor_t _u = lookup_vertex_map(u);
             vertex_descriptor_t _v = lookup_vertex_map(v);
             edge_descriptor_t _e; bool b;
             boost::tie(_e,b) = add_edge(_u, _v, *_graph);
             (*_graph)[_e].weight = WeightPtr(new Weight(w));      
           }
         }
       }
    
       array_graph_t& operator=(const array_graph_t &other)
       {
         crab::CrabStats::count ("Domain.count.copy");
         crab::ScopedCrabStats __st__("Domain.copy");
         if (this != &other)
         {
           _is_bottom      = other._is_bottom;
           _graph          = other._graph;
           _vertex_map     = other._vertex_map;
           _vertices_set   = other._vertices_set;
         }
         return *this;
       }
    
       bool is_bottom() { return this->_is_bottom; }
    
       bool is_top() 
       {
         if (this->is_bottom())
           return false;
         else
         {
           // FIXME: speedup this operation
           canonical();
           edge_iterator it, et;
           for(tie(it,et) = edges(*_graph); it != et; ++it)
           {
             edge_descriptor_t e = *it;
             if (!(*(*_graph)[e].weight).is_top ()) 
               return false;
           }
           return true;
         }
       }
    
       void reduce(ScalarNumDomain scalar)
       {
         if (is_bottom ()) return;
      
         canonical();
      
         edge_iterator it, et;
         for(tie(it,et) = edges(*_graph); it != et; ++it)
         {
           edge_descriptor_t e   = *it;
           VertexNamePtr u = (*_graph)[source(e, *_graph)].name;
           VertexNamePtr v = (*_graph)[target(e, *_graph)].name;
           ScalarNumDomain tmp(scalar);
           linear_constraint_t cst ( variable_t(*u) <= variable_t(*v) - 1);
           tmp += cst;
           if (tmp.is_bottom())
             (*_graph)[e].weight = WeightPtr(new Weight(Weight::bottom()));         
         }
      
         canonical();
       }
    
       // Point-wise application of <= in the weight domain
       bool operator <=(array_graph_t other)
       {
         if (is_bottom())  
           return true;
         else if (other.is_bottom()) 
           return false;
         else
         {
           canonical();
           other.canonical();
           edge_iterator it_1, et_1;
           edge_iterator it_2, et_2;
           for(tie(it_1,et_1) = edges(*_graph); it_1 != et_1; ++it_1)
           {
             edge_descriptor_t e_1 = *it_1;
             vertex_descriptor_t u_1 = source(e_1, *_graph);
             vertex_descriptor_t v_1 = target(e_1, *_graph);
             VertexNamePtr u_name_1 = (*_graph)[u_1].name;
             VertexNamePtr v_name_1 = (*_graph)[v_1].name;
             WeightPtr     weight_1 = (*_graph)[e_1].weight;
             vertex_descriptor_t u_2 = other.lookup_vertex_map(*u_name_1);
             vertex_descriptor_t v_2 = other.lookup_vertex_map(*v_name_1);
             if (edge(u_2,v_2, *other._graph).second)
             {
               edge_descriptor_t e_2 = edge(u_2,v_2, *other._graph).first;
               WeightPtr weight_2 = (*other._graph)[e_2].weight;
               if (!(*weight_1 <= *weight_2))
                 return false;
             }
             else
               CRAB_ERROR ("operator<= with graphs with different adjacency structure");
           }
           return true;
         }
       }
    
       bool operator==(array_graph_t other)
       {
         if (is_bottom()) return other.is_bottom();
         else
           return (*(_vertices_set) == *(other._vertices_set) && 
                   ( *this <= other && other <= *this));
       }
    
       void operator-=(VertexName v) 
       {
         if (!is_bottom ()) 
           remove (v);
       }

       void operator|=(array_graph_t other) {
         *this = *this | other;
       }
    
       // Point-wise join in the weight domain
       array_graph_t operator|(array_graph_t other)
       {
         if (is_bottom())
           return other;
         else if (other.is_bottom())
           return *this;
         else 
           return pointwise_binop<binary_join_op>(*this, other);
       }
    
       // Point-wise widening in the weight domain
       array_graph_t operator||(array_graph_t other)
       {
         if (is_bottom())
           return other;
         else if (other.is_bottom())
           return *this;
         else 
           return pointwise_binop<binary_widening_op>(*this, other);
       }
    
       // Point-wise meet in the weight domain
       array_graph_t operator&(array_graph_t other)
       {
         if (this->is_bottom())
           return *this;
         else if (other.is_bottom())
           return other;
         else {
           return pointwise_binop<binary_meet_op>(*this, other);
         }
       }
    
       // Point-wise narrowing in the weight domain
       array_graph_t operator&&(array_graph_t other)
       {
         if (this->is_bottom())
           return *this;
         else if (other.is_bottom())
           return other;
         else {
           return pointwise_binop<binary_narrowing_op>(*this, other);
         }
       }
    
       void meet_weight (const VertexName &src, const VertexName &dest, 
                         Weight weight)
       {
         if (find_vertex_map(src) && find_vertex_map(dest))
         {
           vertex_descriptor_t u = lookup_vertex_map(src);
           vertex_descriptor_t v = lookup_vertex_map(dest);
           if (edge(u,v,*_graph).second) {
             edge_descriptor_t e = edge(u,v,*_graph).first;
             Weight meet = weight & (*(*_graph)[e].weight);
             (*_graph)[e].weight = WeightPtr(new Weight(meet));
           }
           else {
             vector<VertexName> vertices;
             vector<edge_t>     edges;
             edges.push_back(edge_t(src,dest,weight));
             add(vertices,edges);
           }
         }
       }
    
       void set_weight (const VertexName &src, const VertexName &dest, 
                        Weight weight)
       {
         if (find_vertex_map(src) && find_vertex_map(dest))
         {
           vertex_descriptor_t u = lookup_vertex_map(src);
           vertex_descriptor_t v = lookup_vertex_map(dest);
           if (edge(u,v,*_graph).second) {
             edge_descriptor_t e = edge(u,v,*_graph).first;
             (*_graph)[e].weight = WeightPtr(new Weight(weight));
           }
           else {
             vector<VertexName> vertices;
             vector<edge_t>     edges;
             edges.push_back(edge_t(src,dest,weight));
             add(vertices,edges);
           }
         }
       }
    
       Weight& get_weight (const VertexName &src, const VertexName &dest) 
       {
         if (find_vertex_map(src) && find_vertex_map(dest))
         {
           vertex_descriptor_t u = lookup_vertex_map(src);
           vertex_descriptor_t v = lookup_vertex_map(dest);
           if (edge(u,v,*_graph).second)
           {
             edge_descriptor_t e = edge(u,v,*_graph).first;
             return *((*_graph)[e].weight);
           }
         }
         CRAB_ERROR ("No edge found with given vertices");
       }
    
       void write(crab_os& o) 
       {
         if (is_bottom())
           o << "_|_";
         else
         {
           {
             vertex_iterator it, et;
             o << "(V={";
             for (tie(it, et) = vertices(*_graph); it != et; ++it){
               vertex_descriptor_t u = *it;
               VertexNamePtr u_name  = (*_graph)[u].name;          
               o << *u_name << " ";
             }
             o << "},";
           }
           {
             edge_iterator it, et;
             o << "E={";
             for(tie(it,et) = edges(*_graph); it!= et; ++it)
             {
               edge_descriptor_t e   = *it;
               vertex_descriptor_t u = source(e, *_graph);
               vertex_descriptor_t v = target(e, *_graph);
               VertexNamePtr u_name = (*_graph)[u].name;
               VertexNamePtr v_name = (*_graph)[v].name;
               WeightPtr     weight = (*_graph)[e].weight;
               if (!weight->is_bottom())
                 o << "(" << *u_name << "," << *v_name << "," << *weight << ") ";
             }
             o << "})";
           }
         }
       }
     }; // end class array_graph

 
    /*
      Reduced product of a scalar numerical domain with a weighted array
      graph.
    */
    template<typename ScalarNumDomain, typename WeightDomain, bool IsDistWeight = false>
    class array_graph_domain: 
        public writeable, 
        public numerical_domain<typename ScalarNumDomain::number_t,
                                typename ScalarNumDomain::varname_t>,
        public bitwise_operators<typename ScalarNumDomain::number_t, 
                                 typename ScalarNumDomain::varname_t>, 
        public division_operators<typename ScalarNumDomain::number_t,
                                  typename ScalarNumDomain::varname_t>,
        public array_operators<typename ScalarNumDomain::number_t,
                               typename ScalarNumDomain::varname_t >,
        public pointer_operators<typename ScalarNumDomain::number_t,
                                 typename ScalarNumDomain::varname_t > {

      template<typename Key, typename Value>
      class merge_op_check_equal: public patricia_tree< Key, Value >::binary_op_t {
        boost::optional< Value > apply(Value x, Value y) {
          if (x == y) return x;
          else
            CRAB_ERROR("merging a key with two different values");
        };
        bool default_is_absorbing() { return false; }
      }; 
   
      template<typename Key, typename Value>
      class merge_op_first: public patricia_tree< Key, Value >::binary_op_t {
        boost::optional< Value > apply(Value x, Value y)  {
          return x;
        };
        bool default_is_absorbing() { return false; } 
      }; 
   
      template<typename Key, typename Value>
      class merge_op_second: public patricia_tree< Key, Value >::binary_op_t {  
        boost::optional< Value > apply(Value x, Value y)  {
          return y;
        };
        bool default_is_absorbing() { return false; }
      }; 
   
     template < typename Key, typename Value, 
                typename MergeOp = merge_op_check_equal <Key, Value> >
     class mergeable_map: public writeable {
       
      private:
       typedef patricia_tree< Key, Value > patricia_tree_t;
       typedef typename patricia_tree_t::binary_op_t binary_op_t;
  
      public:
       typedef mergeable_map< Key, Value > mergeable_map_t;
       typedef typename patricia_tree_t::iterator iterator;
       
      private:
       patricia_tree_t _tree;
       
       static patricia_tree_t do_union(patricia_tree_t t1, patricia_tree_t t2) {
         MergeOp o;
         t1.merge_with(t2, o);
         return t1;
       }
       
       mergeable_map(patricia_tree_t t): _tree(t) { }
       
      public:
       
       mergeable_map(): _tree(patricia_tree_t()) { }
  
       mergeable_map(const mergeable_map_t& e): writeable(), _tree(e._tree) { }
       
       mergeable_map_t& operator=(mergeable_map_t e) {
         _tree = e._tree;
         return *this;
       }
       
       iterator begin() { return _tree.begin(); }

       iterator end() { return _tree.end(); }

       std::size_t size(){ return _tree.size(); }

       void set(Key k, Value v) { _tree.insert(k, v); }

       mergeable_map_t& operator-=(Key k) {
         _tree.remove(k);
         return *this;
       }       
       mergeable_map_t operator|(mergeable_map_t e) {
         mergeable_map_t u(do_union(_tree, e._tree));
         return u;
       }
              
       boost::optional<Value> operator[](Key k) { return _tree.lookup(k); }
       void clear() { _tree = patricia_tree_t(); }
       
       void write(crab_os& o) {
         o << "{";
         for (auto it = _tree.begin(); it != _tree.end(); ) {
           Key k = it->first;
           k.write(o);
           o << " -> ";
           Value v = it->second;
           o << v;
           ++it;
           if (it != _tree.end()) {
             o << "; ";
           }
         }
         o << "}";
       }    
     }; // class mergeable_map
      
     public:
      typedef typename ScalarNumDomain::number_t Number;
      typedef typename ScalarNumDomain::varname_t VariableName;
      
      // WARNING: assumes ScalarNumDomain::number_t = WeightDomain::number_t and
      //                  ScalarNumDomain::varname_t = WeightDomain::varname_t
      using typename numerical_domain< Number, VariableName>::linear_expression_t;
      using typename numerical_domain< Number, VariableName>::linear_constraint_t;
      using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
      using typename numerical_domain< Number, VariableName>::variable_t;
      using typename numerical_domain< Number, VariableName>::number_t;
      using typename numerical_domain< Number, VariableName>::varname_t;
      typedef WeightDomain content_domain_t;      
      typedef ScalarNumDomain index_domain_t;      

     private:
      typedef array_graph< VariableName,WeightDomain,ScalarNumDomain,IsDistWeight> array_graph_t;
      typedef array_graph_domain<ScalarNumDomain,WeightDomain,IsDistWeight> array_graph_domain_t;
      
      typedef mergeable_map<VariableName,VariableName> succ_index_map_t;
      typedef boost::shared_ptr< succ_index_map_t > succ_index_map_ptr;
      
      bool  _is_bottom;
      ScalarNumDomain _scalar;        
      array_graph_t _g;        
      // for each array index i we keep track of a special index that
      // represent i+1
      succ_index_map_ptr _succ_idx_map;
      
      void abstract (VariableName v) 
      {
        if (_g.find_vertex_map(v))
        {
          _g.set_incoming(v , WeightDomain::top());
          _g.set_outgoing(v , WeightDomain::top());
          optional<VariableName> succ_v = get_succ_idx(v);
          if (succ_v)
          {
            _g.set_incoming(*succ_v, WeightDomain::top());
            _g.set_outgoing(*succ_v, WeightDomain::top());
          }
        }
      }
      
    
      optional<VariableName> get_succ_idx (VariableName v) const
      {
        return (*_succ_idx_map)[v];
      }

      template <typename VariableFactory>
      VariableName add_variable (Number n, VariableFactory &vfac)
      {
        // FIXME: really big assumption that the variable factory
        // understands strings. For instance, this is not true if the
        // factory is created by Crab-llvm.
        VariableName var_n = vfac["#C" + n.get_str()];

        if (n >= 0)
        {
          _g.insert_vertex(var_n);
          _scalar.assign(var_n, n);
        }
        return var_n;
      }

      void add_variable (VariableName v)
      {
        if (is_array_index(v))
        {
          /*assign to v_succ a fresh var must be always the same*/ 
          VariableName v_succ = v.get_var_factory().get (v.index()); 

          _g.insert_vertex (v);
          _g.insert_vertex (v_succ);
          _succ_idx_map->set (v, v_succ);
      
          /// FIXME: assume that the array element size is 1.

          // --- Enforce: i+ == i+1
          _scalar += linear_constraint_t( variable_t (v_succ) == variable_t(v) + 1);
          // needed if scalar domain is non-relational:
          _g.set_weight (v_succ,v,WeightDomain::bottom ());
        }
      }

      void meet_weight (VariableName i, VariableName j, WeightDomain w)
      {
        add_variable (i);
        add_variable (j);
        _g.meet_weight (i,j,w);
        reduce();
      }

      template <typename VariableFactory>
      void meet_weight (Number i, Number j, WeightDomain w, VariableFactory &vfac)
      {
        _g.meet_weight (add_variable (i, vfac),
                        add_variable (j, vfac) ,
                        w);
        reduce();
      }

      void meet_weight (Number i, VariableName j, WeightDomain w)
      {
        add_variable (j);
        _g.meet_weight (add_variable(i, j.get_var_factory ()),j,w);
        reduce();
      }

      void meet_weight (VariableName i, Number j, WeightDomain w)
      {
        add_variable(i);
        _g.meet_weight(i,add_variable(j, i.get_var_factory ()),w);
        reduce();
      }

      bool IsDefiniteOne (Number x)  { 
        return x == Number (1); 
      }

      bool IsDefiniteOne (VariableName x)  { 
        auto n = _scalar [x].singleton ();
        if (n) return *n == Number (1); 
        else return false;
      }
  
      // x := x op k 
      // Most of the magic happens here.
      template<typename VarNum>
      void apply_helper (operation_t op, VariableName x, VarNum k) 
      {
        if (is_bottom()) return;

        /// step 1: add x_old in the graph
        VariableName x_old = x.get_var_factory ().get (); /*fresh var*/ 
        VariableName x_old_succ = x.get_var_factory ().get (); /*fresh var*/
        _g.insert_vertex(x_old);
        _g.insert_vertex(x_old_succ);
        _succ_idx_map->set(x_old, x_old_succ);

        /// --- Enforce the following relationships:
        ///     { x_old = x, x_old+ = x+, x_old+ = x_old + 1} 

        /// x_old = x
        _scalar.assign(x_old, linear_expression_t(x)); 
        // needed if scalar domain is non-relational:
        // enforcing x_old = x
        _g.set_weight (x_old,x,WeightDomain::bottom ());
        _g.set_weight (x,x_old,WeightDomain::bottom ());

        /// x_old+ = x_old +1
        _scalar += linear_constraint_t(variable_t(x_old_succ) == variable_t(x_old) + 1);      
        // needed if scalar domain is non-relational:
        // enforcing x_old+ < x_old is false
        _g.set_weight (x_old_succ,x_old,WeightDomain::bottom ());
        // enforcing x_old+ < x is false
        _g.set_weight (x_old_succ,x,WeightDomain::bottom ());

        /// x_old+ = x+
        optional<VariableName> x_succ = get_succ_idx(x);
        if (x_succ) {
          _scalar += linear_constraint_t( variable_t(x_old_succ) == variable_t(*x_succ));      
          // needed if scalar domain is non-relational:
          // enforcing x_old + = x+
          _g.set_weight (x_old_succ,*x_succ,WeightDomain::bottom ());
          _g.set_weight (*x_succ,x_old_succ,WeightDomain::bottom ());
          // enforcing x+ < x_old is false
          _g.set_weight (*x_succ,x_old,WeightDomain::bottom ());
        }
        // propagate the scalar constraints to the graph
        reduce();

        /// step 2: abstract all incoming/outgoing edges of x
        abstract(x);

        /// step 3: update the graph with the scalar domain after applying
        ///         x = x op k.
        _scalar.apply(op, x, x, k); 

#if 1
        // This is not needed at all if the scalar domain is relational.
        // Otherwise, we would like to keep the relationship between
        // x_old+ and x. We do it in a completely adhoc way but at least
        // we cover common cases when the array is traversed forward or
        // backwards one element by one.
        if (op == OP_ADDITION && IsDefiniteOne (k)) {
          _g.set_weight (x,x_old_succ,WeightDomain::bottom ());
          _g.set_weight (x_old_succ,x,WeightDomain::bottom ());  
        }
        else if (op == OP_SUBTRACTION && IsDefiniteOne (k)) {
          _g.set_weight (x_old,*x_succ,WeightDomain::bottom ());
          _g.set_weight (*x_succ,x_old,WeightDomain::bottom ());  
        }
#endif 

        if (x_succ){
          _scalar -= *x_succ;
          /// --- Enforce x+ == x + 1 
          _scalar += linear_constraint_t(variable_t(*x_succ) == variable_t(x) + 1); 
          // needed if scalar domain is non-relational:
          _g.set_weight (*x_succ,x,WeightDomain::bottom ());
        }

        /* { x = x op k, x+ = x+1} */
        reduce();

        /// step 4: delete x_old
        _g -= x_old;
        _g -= x_old_succ;
        (*_succ_idx_map) -= x_old;
        _scalar -= x_old;
        _scalar -= x_old_succ;

        //this->reduce();
      }

      //! In case we can statically determine which variables should be
      //  considered array indexes. Note that any subset is sound but it
      //  might be imprecise. By default we consider all.
      bool is_array_index(VariableName v) const {  
        return true; 
      }

      // model array reads: return the weight from the edge i to i+
      WeightDomain array_read (VariableName i) 
      {
        if (is_bottom()) 
          return WeightDomain::bottom();
  
        if (!is_array_index(i)) 
          return WeightDomain::top();
    
        //this->reduce();
        optional<VariableName> i_succ = get_succ_idx(i);
        if (i_succ) 
          return _g.get_weight(i, *i_succ);

        CRAB_ERROR ("There is no successor index associated with ",i);
      }
 
      // model array writes
      void array_write (VariableName arr, VariableName i, WeightDomain w)
      {
        if (is_bottom()) return;
        //this->reduce();

        //--- strong update
        optional<VariableName> i_succ = get_succ_idx(i);
        if (!i_succ) 
          CRAB_ERROR ("There is no successor index associated with ",i);

        WeightDomain& old_w = _g.get_weight(i, *i_succ);
        old_w -= arr;
        _g.meet_weight(i, *i_succ, w);
        WeightDomain new_w = _g.get_weight(i, *i_succ);
    
        //--- weak update: 
        // An edge (p,q) must be weakened if p <= i <= q and p < q
        typename array_graph_t::edge_iterator it, et;
        for(tie(it,et) = edges(*_g._graph); it!= et; ++it)
        {
          typename array_graph_t::edge_descriptor_t e = *it;
          typename array_graph_t::VertexNamePtr  p = (*_g._graph)[source(e, *_g._graph)].name;
          typename array_graph_t::VertexNamePtr  q = (*_g._graph)[target(e, *_g._graph)].name;
          typename array_graph_t::WeightPtr weight = (*_g._graph)[e].weight;
          if ( ((*p == i) &&  (*q == *i_succ)) || weight->is_bottom())
            continue;
          // p < q 
          ScalarNumDomain tmp(_scalar);
          tmp += linear_constraint_t( variable_t(*p) <= variable_t(i));      
          tmp += linear_constraint_t( variable_t(*i_succ)  <= variable_t(*q));     
          if (tmp.is_bottom())
            continue;
          // p <= i <= q and p < q
          typename array_graph_t::binary_join_op join;
          (*_g._graph)[e].weight = join (weight, 
                                         typename array_graph_t::WeightPtr (new WeightDomain(new_w)));
        }
        _g.canonical();

      }
  
      void set_to_bottom()
      {
        _is_bottom = true;
        _scalar = ScalarNumDomain::bottom();
        _g = array_graph_t::bottom();
        _succ_idx_map->clear();
      }

      array_graph_domain(ScalarNumDomain scalar, 
                         array_graph_t g, 
                         succ_index_map_ptr map): 
          writeable(), 
          _is_bottom(false), 
          _scalar(scalar), 
          _g(g), 
          _succ_idx_map(new succ_index_map_t(*map)) 
      { 
        if (_scalar.is_bottom() || _g.is_bottom())
          set_to_bottom();
        else
          reduce();
      }

     public:

      array_graph_domain(): 
          writeable(), 
          _is_bottom(false), 
          _scalar(ScalarNumDomain::top()), 
          _g(array_graph_t::top()), 
          _succ_idx_map(new succ_index_map_t()) 
      { }

      static array_graph_domain_t top() 
      {
        return array_graph_domain(ScalarNumDomain::top(), 
                                  array_graph_t::top(), 
                                  succ_index_map_ptr(new succ_index_map_t()));
      }
  
      static array_graph_domain_t bottom() 
      {
        return array_graph_domain(ScalarNumDomain::bottom(), 
                                  array_graph_t::bottom(),
                                  succ_index_map_ptr(new succ_index_map_t()));
      }
      
      array_graph_domain(const array_graph_domain_t& other): 
          writeable(), 
          _is_bottom(other._is_bottom),
          _scalar(other._scalar) , 
          _g(other._g), 
          _succ_idx_map(new succ_index_map_t(*other._succ_idx_map))
      {
        crab::CrabStats::count (getDomainName() + ".count.copy");
        crab::ScopedCrabStats __st__(getDomainName() + ".copy");
      }
  
      array_graph_domain_t& operator=(array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.copy");
        crab::ScopedCrabStats __st__(getDomainName() + ".copy");

        if (this != &other)
        {
          this->_is_bottom = other._is_bottom;
          this->_scalar = other._scalar;
          this->_g = other._g;
          this->_succ_idx_map = other._succ_idx_map;
        }
        return *this;
      }
  
      bool is_bottom() { return _is_bottom; }

      // --- top operation in the graph is expensive because we need to
      //     traverse the whole graph and check each edge.
      // bool is_top() { return (_scalar.is_top() && _g.is_top()); }
    
      bool is_top() { return (_scalar.is_top()); }
  
      void reduce ()
      {
        crab::CrabStats::count (getDomainName() + ".count.reduce");
        crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

        if (is_bottom ()) return; 
      
        domain_traits<ScalarNumDomain>::normalize(_scalar);

        if (_scalar.is_bottom() || _g.is_bottom())
          set_to_bottom();
        else
          _g.reduce(_scalar);
      }

      bool operator<=(array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.leq");
        crab::ScopedCrabStats __st__(getDomainName() + ".leq");

        if (is_bottom ()) {
          return true;
        } else if (other.is_bottom ()) {
          return false;
        } else {
          return (_scalar <= other._scalar && _g <= other._g);
        }
      }

      void operator|=(array_graph_domain_t other)  {
        *this = *this | other;
      }
  
      array_graph_domain_t operator|(array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.join");
        crab::ScopedCrabStats __st__(getDomainName() + ".join");

        if (is_bottom ()) {
          return other;
        } else if (other.is_bottom ()) {
          return *this;
        } else {
          succ_index_map_ptr map(new succ_index_map_t(*(_succ_idx_map) | 
                                                      *(other._succ_idx_map)));
          return array_graph_domain_t(_scalar | other._scalar, 
                                      _g | other._g, map); 
        }
      }
  
      array_graph_domain_t operator&(array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.meet");
        crab::ScopedCrabStats __st__(getDomainName() + ".meet");

        if (is_bottom () || other.is_bottom ()) {
          return bottom();
        } else {
          succ_index_map_ptr map(new succ_index_map_t(*(_succ_idx_map) | 
                                                      *(other._succ_idx_map)));
          return array_graph_domain_t(_scalar & other._scalar, 
                                      _g & other._g, map);
        }
      }
  
      array_graph_domain_t operator||(array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.widening");
        crab::ScopedCrabStats __st__(getDomainName() + ".widening");

        if (is_bottom ())  return other;
        else if (other.is_bottom ())  return *this;
        else 
        {
          succ_index_map_ptr map (new succ_index_map_t(*(_succ_idx_map) | 
                                                       *(other._succ_idx_map)));
          array_graph_domain_t widen (_scalar || other._scalar, 
                                      _g || other._g, map);
          CRAB_LOG("array-graph" , crab::outs() << "Widening: " << *this<<"\n";);
          return widen;
        }
      }

      template<typename Thresholds>
      array_graph_domain_t widening_thresholds (array_graph_domain_t other, 
                                                const Thresholds & /*ts*/) {
        return (*this || other);
      }
        
      array_graph_domain_t operator&& (array_graph_domain_t other) 
      {
        crab::CrabStats::count (getDomainName() + ".count.narrowing");
        crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

        if (is_bottom ()|| other.is_bottom ()) 
          return bottom();
        else 
        {
          succ_index_map_ptr map (new succ_index_map_t(*(_succ_idx_map) | 
                                                       *(other._succ_idx_map)));
          return array_graph_domain_t (_scalar && other._scalar, 
                                       _g && other._g, map);
        }
      }
  
      void operator-=(VariableName var)
      {
        crab::CrabStats::count (getDomainName() + ".count.forget");
        crab::ScopedCrabStats __st__(getDomainName() + ".forget");

        if (is_bottom ()) return;

        // scalar domain
        _scalar -= var;
        _g -= var;
        optional<VariableName> var_succ = get_succ_idx(var);
        if (var_succ) {
          _scalar -= *var_succ;
          _g -= *var_succ;        
          (*_succ_idx_map) -= var;
        }

        // graph domain
        typename array_graph_t::edge_iterator it, et;
        for(tie(it,et) = edges(*(_g._graph)); it!= et; ++it) {
          auto e  = *it;
          auto weight = (*(_g._graph))[e].weight;
          (*weight) -= var;
        }      
        // this->reduce();
      }
  
      /////
      // Transfer functions
      /////

      void operator += (linear_constraint_system_t csts) 
      {
        crab::CrabStats::count (getDomainName() + ".count.add_constraints");
        crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

        if (is_bottom()) return;
    
        // graph domain: make sure that all the relevant variables
        // (included special "0") are inserted in the graph
        for (auto cst : csts) {
          // TODO
          //Number n = cst.expression().constant();
          //if (n == 0) add_variable(n, vfac);
          for (auto v : cst.variables())
            add_variable (v.name());
        }

        _scalar += csts;
        reduce();

        CRAB_LOG("array-graph", 
                 crab::outs() << "Assume("<< csts<< ") --- "<< *this<<"\n";);
      }

      void assign (VariableName x, linear_expression_t e) 
      {
        crab::CrabStats::count (getDomainName() + ".count.assign");
        crab::ScopedCrabStats __st__(getDomainName() + ".assign");

        if (is_bottom()) return;

        if (optional<variable_t> y = e.get_variable())
        {
          if ((*y).name() == x) return;
        }

        // scalar domain
        _scalar.assign(x, e);
   
        // graph domain
        if (e.is_constant() && (e.constant() == 0))
          add_variable (e.constant(), x.get_var_factory ());

        if (_g.find_vertex_map(x))
        {
          abstract(x);
          // wrong results if we do not restore the relationship between x
          // and x+ in the scalar domain
          optional<VariableName> x_succ = get_succ_idx(x);      
          if (x_succ){
            _scalar -= *x_succ;
            /// --- Enforce x+ == x+1
            _scalar += linear_constraint_t( variable_t(*x_succ) == variable_t(x) + 1);        
            // needed if scalar domain is non-relational:
            _g.set_weight (*x_succ,x,WeightDomain::bottom ());
          }
        }
        else
          add_variable(x);
    
        reduce();

        CRAB_LOG("array-graph", 
                 crab::outs() << "Assign "<<x<<" := "<<e<<" ==> "<<*this<<"\n";);
      }

      void apply (operation_t op, VariableName x, VariableName y, Number z) 
      {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        assign (x, linear_expression_t(y));
        apply_helper<Number> (op, x, z);

        CRAB_LOG("array-graph",
                 crab::outs() << "Apply "<<x<<" := "<<y<<" "<<op<<" "<<z<<" ==> "<<*this<<"\n";); 
      }

      void apply(operation_t op, VariableName x, VariableName y, VariableName z) 
      {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        assign (x, linear_expression_t(y));
        apply_helper<VariableName> (op, x, z);

        CRAB_LOG("array-graph", 
                 crab::outs() << "Apply "<<x<<" := "<<y<<" "<<op<<" "<<z<<" ==> "<<*this<<"\n";);
      }

      void apply(operation_t op, VariableName x, Number k) 
      {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        apply_helper<Number> (op, x, k);

        CRAB_LOG("array-graph",
                 crab::outs() << "Apply "<<x<<" := "<<x<<" "<<op<<" "<<k<<" ==> "<<*this<<"\n";);
      }


      // bitwise_operators_api
      void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
        assign (x, variable_t (y));
      }
      
      void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
        assign (x, k);
      }
      
      void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
        CRAB_WARN ("bitwise operations not implemented in array_graph");
      }
      
      void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
        CRAB_WARN ("bitwise operations not implemented in array_graph");
      }
      
      // division_operators_api
      void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
        CRAB_WARN ("division operations not implemented in array_graph");
      }
      
      void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
        CRAB_WARN ("division operations not implemented in array_graph");
      }

      // array_operators_api
        
      virtual void array_init (VariableName a, const vector<ikos::z_number> &values)  override {
        CRAB_WARN ("array_graph_domain init not implemented");
      }

      virtual void array_assume (VariableName a, 
                                 boost::optional<Number> lb, boost::optional<Number> ub) override {
        CRAB_WARN ("array_graph_domain assume not implemented");
      }

      virtual void array_load (VariableName lhs, VariableName arr, VariableName idx, 
                               z_number /*bytes*/) override {
        crab::CrabStats::count (getDomainName() + ".count.load");
        crab::ScopedCrabStats __st__(getDomainName() + ".load");

        WeightDomain w = array_read (idx);
        // --- Simplification wrt Gange et.al.:
        //     Only non-relational invariants are passed from the graph
        //     domain to the scalar domain.
        //     We use operator[] as the conversion function
        _scalar.set (lhs, w [arr]);

        CRAB_LOG("array-graph",
                 crab::outs() << "Array read "<<lhs<<" := "<< arr<<"["<<idx<<"] ==> "
                           << *this <<"\n";);    
      }

      virtual void array_store (VariableName arr, VariableName idx, linear_expression_t val,
                                 z_number /*n_bytes*/, bool /*is_singleton*/) override {

        crab::CrabStats::count (getDomainName() + ".count.store");
        crab::ScopedCrabStats __st__(getDomainName() + ".store");

        // --- Simplification wrt Gange et.al.:
        //     Only non-relational invariants are passed from the scalar
        //     domain to the graph domain.
        //     We use operator[] as the conversion function
        WeightDomain w = WeightDomain::top ();
        if (val.is_constant ())
          w.assign (arr, val);      
        else if (auto v = val.get_variable ()){
          w.set (arr, _scalar[(*v).name()]);      
        }
        else {
          // If you see this warning you can switch to intervals.
          crab::outs() << "Warning: scalar domain does not support assignments with arbitrary rhs.\n";
        }
        array_write (arr, idx, w);

        CRAB_LOG("array-graph",
                 crab::outs() << "Array write "<<arr<<"["<<idx<<"] := "<<val<< " ==> "<< *this <<"\n";);
      }
    
      void write(crab_os& o) 
      {
        o << "(" ;
#if 1
        // less verbose: remove the special variables i+ from the scalar
        // domain
        ScalarNumDomain inv (_scalar);
        for(auto p : *_succ_idx_map) {
          inv -= p.second;
        }
        o << inv;
#else
        o << _scalar;
#endif 
        o << "," << _g;
        o << ")";
      }

      linear_constraint_system_t to_linear_constraint_system (){
        CRAB_ERROR ("array_graph: to_linear_constraint_system not implemented");
      }

      static string getDomainName () {
        string name ("ArrayGraph(" + 
                     ScalarNumDomain::getDomainName () +  "," + 
                     WeightDomain::getDomainName () +
                     ")");
        return name;
      }

    }; // end array_graph_domain


    }// namespace domain_traits

} // namespace crab

#endif 
