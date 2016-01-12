/*******************************************************************************
 * Default implementations for non-standard abstract operations
 *
 * We separate domain_traits.hpp from domain_traits_impl.hpp so that
 * we do not create cyclic dependencies when including headers. In
 * this way, abstract domains can also call domain_traits operations.
 *
 ******************************************************************************/


#ifndef DOMAINS_TRAITS_IMPL_HPP
#define DOMAINS_TRAITS_IMPL_HPP

#include <crab/domains/domain_traits.hpp>

namespace crab {

   namespace domain_traits {

        template <typename NumDomain >
        void normalize(NumDomain& inv) {
        }
     
        template <typename NumDomain, typename Iterator >
        void forget(NumDomain& inv, Iterator it, Iterator end){
          // inefficient if after each forget the domain requires
          // normalization
          for (auto v : boost::make_iterator_range (it,end)){
            inv -= v; 
          }
        }

        template <typename NumDomain, typename Iterator>
        void project(NumDomain& inv, Iterator begin, Iterator end) {
          // CRAB_WARN(NumDomain::getDomainName (), 
          //           " project may lose if relation or disjunctive domain");
          // lose precision if relational or disjunctive domain
          NumDomain res = NumDomain::top ();
          for (auto v : boost::make_iterator_range (begin, end)){
            res.set (v, inv[v]); 
          }
          std::swap (inv, res);
        }
     
        template <typename Domain, typename VariableName >
        void expand (Domain& inv, VariableName x, VariableName new_x) {
          // CRAB_WARN(Domain::getDomainName (),
          //           " expand may lose if relation or disjunctive domain");
          // lose precision if relational or disjunctive domain
          inv.set (new_x , inv [x]);
        }

        template <typename VariableName, typename NumDomain1, typename NumDomain2>
        void push (const VariableName& x, NumDomain1 from, NumDomain2& to) {
          CRAB_WARN(" push operator from ", NumDomain1::getDomainName (), " to ",
                    NumDomain2::getDomainName (), " has no effect ");
        }
   
        template <typename Domain, typename VariableName>
        void array_init (Domain& inv, VariableName a,
                         const vector<ikos::z_number>& vals) {
        }
     
        template <typename Domain, typename VariableName, typename Number>
        void assume_array (Domain& inv, VariableName a, Number val) {
        }

        template <typename Domain, typename VariableName, typename Number>
        void assume_array (Domain& inv, VariableName a, interval<Number> val) {
        }
     
        template <typename Domain, typename VariableName >
        void array_load (Domain& inv, VariableName lhs, 
                         VariableName a, VariableName i,
                         ikos::z_number n_bytes) {
        }

        template <typename Domain, typename VariableName >
        void array_store (Domain& inv, VariableName a, 
                          VariableName i, typename Domain::linear_expression_t v,
                          ikos::z_number n_bytes, bool is_singleton) {
        }

       // temporary for profiling domains
       template <typename Domain>
       void print_stats (Domain inv) { }

     
   } //end namespace domain_traits
}// end crab
#endif 
