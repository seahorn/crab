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

   // Special operations needed by the checker
   template<typename Domain>
   class checker_domain_traits{
   public:
     typedef typename Domain::linear_constraint_t linear_constraint_t;
     
     // Return true if inv entails cst.
     static bool entail(Domain& inv, const linear_constraint_t& cst) {
       if (inv.is_bottom()) return true;
       if (cst.is_tautology ()) return true;
       if (cst.is_contradiction ()) return false;

       CRAB_LOG("checker-entailment",
		linear_constraint_t tmp(cst);
		crab::outs() << "Checking whether\n" << inv << "\nentails " << tmp << "\n";);
       
       linear_constraint_t neg_cst = cst.negate();
       Domain cst_inv = inv;
       cst_inv += neg_cst;
       bool res = cst_inv.is_bottom();

       CRAB_LOG("checker-entailment",
		if (res) {
		  crab::outs() << "\t**entailment holds.\n";
		} else  {
		  crab::outs() << "\t**entailment does not hold.\n";
		});

       // Note: we cannot convert cst into Domain and then use the <=
       //       operator. The problem is that we cannot know for sure
       //       whether Domain can represent precisely cst. It is not
       //       enough to do something like
       // 
       //       Dom dom = cst;
       //       if (dom.is_top()) { ... }
       
       return res; 
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

   // Special operations for applying reduction between domains.
   template<typename Domain>
   class reduced_domain_traits {
    public:
     typedef typename Domain::variable_t variable_t;     
     typedef typename Domain::linear_constraint_t linear_constraint_t;     
     typedef typename Domain::linear_constraint_system_t linear_constraint_system_t;

     // extract linear constraints from dom involving x and store in
     // ctsts
     static void extract(Domain& dom, const variable_t& x, 
			 linear_constraint_system_t& csts,
			 bool only_equalities){ 
       auto all_csts = dom.to_linear_constraint_system();
       for(auto cst: all_csts) {
       	 if (only_equalities && (!cst.is_equality())) {
       	   continue;
       	 }
       	 auto vars = cst.variables();
       	 if (vars[x]) {
       	   csts += cst;
       	 }
       }
     }
   };

   // Experimental:
   // 
   // Special operations needed by the array_sparse_graph domain.
   // TO BE REMOVED
   template<typename Domain>
   class array_sgraph_domain_traits {
    public:
     typedef typename Domain::linear_constraint_t linear_constraint_t;
     typedef typename Domain::variable_t variable_t;

     // FIXME: this does similar thing to checker_domain_traits<Domain>::entail
     static bool is_unsat(Domain& inv, linear_constraint_t cst) { 
       Domain copy(inv);
       copy += cst;
       return copy.is_bottom();
     }

     static void active_variables(Domain& inv, std::vector<variable_t>& out)
     { CRAB_ERROR("operation active_variables not implemented"); }
   };


 } // end namespace domains   
}// end namespace crab

