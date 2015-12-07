/*******************************************************************************
 * Extend abstract domains with non-standard operations or types
 ******************************************************************************/

#ifndef DOMAINS_TRAITS_HPP
#define DOMAINS_TRAITS_HPP

#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/intervals_congruences.hpp>
#include <crab/domains/dbm.hpp>
#include <crab/domains/term_equiv.hpp>
#ifdef HAVE_LDD
#include <crab/domains/boxes.hpp>
#endif 
#include <crab/domains/array_graph.hpp>
#include <crab/domains/array_smashing.hpp>

namespace crab {

   namespace domain_traits {

      // Normalize the abstract domain if such notion exists.
      template <typename AbsNumDomain>
      void normalize(AbsNumDomain& inv); 
    
      // Remove all variables [it, end)
      template <typename AbsNumDomain, typename Iterator>
      void forget(AbsNumDomain& inv, Iterator begin, Iterator end); 

      // Forget all variables except [begin, end)
      template <typename AbsNumDomain, typename Iterator>
      void project(AbsNumDomain& inv, Iterator begin, Iterator end); 
    
      // Make a new copy of x without relating x with new_x
      template <typename AbsDomain, typename VariableName>
      void expand (AbsDomain& inv, VariableName x, VariableName new_x);
   
      ////// 
      /// Special operations for arrays
      //////

      template <typename AbsDomain, typename VariableName>
      void array_init (AbsDomain& inv, VariableName a, 
                       const vector<ikos::z_number>& vals); 
    
      template <typename AbsDomain, typename VariableName, typename Number>
      void assume_array (AbsDomain& inv, VariableName a, Number val); 
    
      template <typename AbsDomain, typename VariableName, typename Number>
      void assume_array (AbsDomain& inv, VariableName a, interval<Number> val);
      
      template <typename AbsDomain, typename VariableName>
      void array_load (AbsDomain& inv, VariableName lhs, 
                       VariableName a, VariableName i,
                       ikos::z_number n_bytes);

      template <typename AbsDomain, typename VariableName>
      void array_store (AbsDomain& inv, VariableName a, 
                        VariableName i, typename AbsDomain::linear_expression_t v,
                        ikos::z_number n_bytes, bool is_singleton); 

   } // end namespace domain_traits
}// end namespace crab

#endif 
