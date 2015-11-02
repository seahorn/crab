/*******************************************************************************
 * Extend abstract domains with non-standard operations or types
 ******************************************************************************/

#ifndef DOMAINS_TRAITS_HPP
#define DOMAINS_TRAITS_HPP

#include <crab/common/types.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/intervals_congruences.hpp>
#include <crab/domains/dbm.hpp>
#include <crab/domains/term_equiv.hpp>
#include <crab/domains/array_graph.hpp>
#include <crab/domains/array_smashing.hpp>

namespace crab {

   namespace domain_traits {

      template <typename AbsNumDomain>
      void normalize(AbsNumDomain& inv); 
    
      template <typename AbsNumDomain, typename Iterator>
      void forget(AbsNumDomain& inv, Iterator begin, Iterator end); 

      template <typename AbsNumDomain, typename Iterator>
      void project(AbsNumDomain& inv, Iterator begin, Iterator end); 
    
      template <typename AbsDomain, typename VariableName>
      void expand (AbsDomain& inv, VariableName x, VariableName new_x);
    
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

      // Stateless class to convert the internal representation of an
      // abstract domain to a formula. By default, no conversion is
      // done. See domain_traits_impl.hpp for a conversion to e.g.,
      // conjunctive linear arithmetic.
      template<typename AbsDomain, typename Formula>
      class absdom_to_formula {

        static void check_same_type () {
          if (!std::is_same<AbsDomain,Formula>::value)
            CRAB_ERROR ("absdom_to_formula method expects two parameters with same types");
        }

       public:
        typedef AbsDomain type_from;
        typedef Formula type_to;
                
        //! Convert the abstract domain to a formula
        static Formula marshall (AbsDomain& abs) {
          //  FIXME: should be `const AbsDomain&` but we need to fix first
          //  constness of the abstract domain methods.
          check_same_type ();
          return abs;
        }

        //! Convert from a formula to the abstract domain.
        static AbsDomain unmarshall (Formula inv) {
          check_same_type ();
          return inv;
        }

        //! Return true 
        static Formula mkTrue () {
          check_same_type ();
          return AbsDomain::top ();
        }

        //! Return false
        static Formula mkFalse () {
          check_same_type ();
          return AbsDomain::bottom ();
        }
      };
   } // end namespace domain_traits
}// end namespace crab

#endif 
