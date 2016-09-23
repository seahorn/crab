/*******************************************************************************
 * Extend abstract domains with non-standard operations.
 * Some of them might be moved into the domains later, others might
 * stay here.
 ******************************************************************************/

#ifndef DOMAINS_TRAITS_HPP
#define DOMAINS_TRAITS_HPP

#include <crab/common/bignums.hpp>
#include <crab/domains/intervals.hpp>

#include <boost/range/iterator_range.hpp>
#include <vector>

namespace crab {

 namespace domains {

   template<typename Domain>
   class domain_traits {

    public:

     // Initialization of static data 
     // XXX: it doesn't take inv as argument because this method
     // should only access to static data.
     template<class CFG>
     static void do_initialization(CFG cfg) { }

     // Normalize the abstract domain if such notion exists.
     static void normalize (Domain& inv) { }

     // Remove all variables [begin, end)
     template<typename Iter>
     static void forget (Domain& inv, Iter begin, Iter end) {
       // -- inefficient if after each forget the domain requires
       //    normalization
       for (auto v : boost::make_iterator_range (begin,end)){
         inv -= v; 
       }
     }

     // Forget all variables except [begin, end)
     template <typename Iter>
     static void project(Domain& inv, Iter begin, Iter end){
       // -- lose precision if relational or disjunctive domain
       Domain res = Domain::top ();
       for (auto v : boost::make_iterator_range (begin, end)){
         res.set (v, inv[v]); 
       }
       std::swap (inv, res);
     }
         
     // Make a new copy of x without relating x with new_x
     template <typename VariableName>
     static void expand (Domain& inv, VariableName x, VariableName new_x) {
       // -- lose precision if relational or disjunctive domain
       inv.set (new_x , inv [x]);
     }

   };

   template<typename Domain1, typename Domain2>
   class product_domain_traits {
    public:

     // To perform reduction between domains
     template <typename VariableName>
     static void push (const VariableName& x, Domain1 from, Domain2& to){ }
     
   };

   // Special operations needed by the array_sparse_graph domain
   template<typename Domain>
   class array_sgraph_domain_traits {

    public:

     static bool is_unsat(Domain &inv, typename Domain::linear_constraint_t cst) { 
       Domain copy(inv);
       copy += cst;
       return copy.is_bottom();
     }

     static std::vector<typename Domain::varname_t> active_variables(Domain &inv) {
       CRAB_ERROR("operation active_variables not implemented");
     }
     
   };


 } // end namespace domains   
}// end namespace crab

#endif 
