/*******************************************************************************
 * Extend abstract domains with non-standard operations.
 * Some of them might be moved into the domains later, others might
 * stay here.
 ******************************************************************************/

#pragma once

#include <crab/common/bignums.hpp>
#include <crab/domains/intervals.hpp>

#include <boost/range/iterator_range.hpp>
#include <vector>

namespace crab {

 namespace domains {

   template<typename Domain>
   class domain_traits {
    public:
     typedef typename Domain::variable_t variable_t;

     // Initialization of static data 
     // XXX: it doesn't take inv as argument because this method
     // should only access to static data.
     template<class CFG>
     static void do_initialization(CFG cfg) { }

     // Normalize the abstract domain if such notion exists.
     static void normalize (Domain& inv) { }

     // Remove all variables [begin, end)
     template<typename VarIter>
     static void forget (Domain& inv, VarIter begin, VarIter end) {
       // -- inefficient if after each forget the domain requires
       //    normalization
       for (auto v : boost::make_iterator_range (begin,end)){
         inv -= v; 
       }
     }

     // Forget all variables except [begin, end)
     template <typename VarIter>
     static void project(Domain& inv, VarIter begin, VarIter end){
       // -- lose precision if relational or disjunctive domain
       Domain res = Domain::top ();
       for (auto v : boost::make_iterator_range (begin, end)){
         res.set (v, inv[v]); 
       }
       std::swap (inv, res);
     }
         
     // Make a new copy of x without relating x with new_x
     static void expand (Domain& inv, variable_t x, variable_t new_x) {
       // -- lose precision if relational or disjunctive domain
       inv.set (new_x , inv [x]);
     }
   };

   // For domains used by the checker
   template<typename Domain>
   class checker_domain_traits{
   public:
     typedef typename Domain::linear_constraint_t linear_constraint_t;
     
     // Return true if inv entails cst.
     // if cst cannot be represented by Domain then return false.    
     static bool entail(Domain& inv, const linear_constraint_t& cst) {
       if (inv.is_bottom()) return true;
       if (cst.is_tautology ()) return true;
       if (cst.is_contradiction ()) return false;
       
       Domain cst_inv;
       cst_inv += cst;
       // cst cannot be represented by the domain.
       if (cst_inv.is_top ()) return false;
      
       return inv <= cst_inv; 
     }
     
     // Return true if inv intersects with cst.
     static bool intersect(Domain& inv, const linear_constraint_t& cst) {
       if (inv.is_bottom () || cst.is_contradiction ()) return false;
       if (inv.is_top () || cst.is_tautology ()) return true;
       
       Domain cst_inv;
       cst_inv += cst;
       return (!(cst_inv & inv).is_bottom ());
     }
   };
   
   template<typename Domain1, typename Domain2>
   class product_domain_traits {
    public:
     typedef typename Domain1::variable_t variable_t;
     
     // Perform customized reduction between domains.
     static void push (const variable_t& x, Domain1 from, Domain2& to){ }
   };

   // Special operations needed by the array_sparse_graph domain
   // TO BE REMOVED
   template<typename Domain>
   class array_sgraph_domain_traits {
    public:
     typedef typename Domain::linear_constraint_t linear_constraint_t;
     typedef typename Domain::variable_t variable_t;

     // This does similar thing to checker_domain_traits<Domain>::entail
     static bool is_unsat(Domain &inv, linear_constraint_t cst) { 
       Domain copy(inv);
       copy += cst;
       return copy.is_bottom();
     }

     static std::vector<variable_t> active_variables(Domain &inv)
     { CRAB_ERROR("operation active_variables not implemented"); }
   };


 } // end namespace domains   
}// end namespace crab

