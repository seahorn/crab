#ifndef TOPOLOGICAL_ORDER_HPP__
#define TOPOLOGICAL_ORDER_HPP__

#include <boost/graph/topological_sort.hpp>

///Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/cg/Sccg.hpp>
#include <crab/cg/SccgBgl.hpp>


/* 
   Topological order of a call graph
 */

using namespace boost;

namespace crab {

   namespace cg {

     // res contains the reverse topological order of the SCC graph g
     template <typename CG>
     void rev_topo_sort (const SccGraph<CG> &g, std::vector<typename CG::node_t>& res) {

       typedef boost::unordered_map< typename CG::node_t, default_color_type > color_map_t;
       typedef boost::associative_property_map< color_map_t > property_color_map_t;

       color_map_t colormap;

       for (auto const &v: boost::make_iterator_range (vertices (g))) 
         colormap [v] = default_color_type();

       res.reserve(num_vertices (g));
       
       boost::topological_sort(g,
                               std::back_inserter(res),
                               color_map(property_color_map_t(colormap)));
     }
   
   } // end namespace cg
} // end namespace crab

#endif 
