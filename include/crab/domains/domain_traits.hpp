/*******************************************************************************
 * Extend abstract domains with non-standard operations or types
 ******************************************************************************/

#ifndef DOMAINS_TRAITS_HPP
#define DOMAINS_TRAITS_HPP

// If a domain provides a different implementation from the default
// one (available in domain_traits_impl.hpp) then its header file
// should be included here
#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/dis_intervals.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/dbm.hpp>
#include <crab/domains/var_packing_naive_dbm.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/domains/naive_dbm.hpp>
#include <crab/domains/term_equiv.hpp>
#ifdef HAVE_APRON
#include <crab/domains/apron_domains.hpp>
#endif 
#ifdef HAVE_LDD
#include <crab/domains/boxes.hpp>
#endif 
#include <crab/domains/array_graph.hpp>
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/combined_domains.hpp>

namespace crab {

   namespace domain_traits {

      // Normalize the abstract domain if such notion exists.
      template <typename NumDomain>
      void normalize(NumDomain& inv); 
    
      // Remove all variables [it, end)
      template <typename NumDomain, typename Iterator>
      void forget(NumDomain& inv, Iterator begin, Iterator end); 

      // Forget all variables except [begin, end)
      template <typename NumDomain, typename Iterator>
      void project(NumDomain& inv, Iterator begin, Iterator end); 
    
      // Make a new copy of x without relating x with new_x
      template <typename Domain, typename VariableName>
      void expand (Domain& inv, VariableName x, VariableName new_x);

      // To perform reduction between domains
      template <typename VariableName, typename NumDomain1, typename NumDomain2>
      void push (const VariableName& x, NumDomain1 from, NumDomain2& to);
   
      ////// 
      /// Special operations for arrays
      //////

      template <typename Domain, typename VariableName>
      void array_init (Domain& inv, VariableName a, 
                       const vector<ikos::z_number>& vals); 
    
      template <typename Domain, typename VariableName, typename Number>
      void assume_array (Domain& inv, VariableName a, Number val); 
    
      template <typename Domain, typename VariableName, typename Number>
      void assume_array (Domain& inv, VariableName a, interval<Number> val);
      
      template <typename Domain, typename VariableName>
      void array_load (Domain& inv, VariableName lhs, 
                       VariableName a, VariableName i,
                       ikos::z_number n_bytes);

      template <typename Domain, typename VariableName>
      void array_store (Domain& inv, VariableName a, 
                        VariableName i, typename Domain::linear_expression_t v,
                        ikos::z_number n_bytes, bool is_singleton); 

       // temporary for profiling domains
      template <typename Domain>
      void print_stats (Domain inv); 

   } // end namespace domain_traits
}// end namespace crab

#endif 
