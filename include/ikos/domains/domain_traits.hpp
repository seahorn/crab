/*******************************************************************************
 * Extend abstract domains with non-standard operations or types
 ******************************************************************************/

#ifndef IKOS_DOMAINS_TRAITS_HPP
#define IKOS_DOMAINS_TRAITS_HPP

#include <ikos/common/types.hpp>
#include <ikos/domains/intervals.hpp>
#include <ikos/domains/congruences.hpp>
#include <ikos/domains/intervals_congruences.hpp>
#include <ikos/domains/octagons.hpp>
#include <ikos/domains/dbm.hpp>
#include <ikos/domains/array_graph.hpp>
#include <ikos/domains/array_smashing.hpp>

namespace ikos {

namespace domain_traits {

/// / Default implementation
// template <typename AbsNumDomain >
// void normalize(AbsNumDomain& inv) {
// }

// // Specialized version
// template <typename Number, typename VariableName>
// void normalize(octagon<Number, VariableName>& inv) {
//    inv.normalize();
// }

// // Specialized version
// template <typename Number, typename VariableName>
// void normalize(DBM<Number, VariableName>& inv) {
//    inv.normalize();
// }

// Default implementation
template <typename AbsDomain, typename VariableName >
void array_init (AbsDomain& inv, VariableName arr) {
}

////
// (Partial) specialized versions
////
template <typename BaseDomain, typename VariableName>
void array_init (array_smashing<BaseDomain,z_number,VariableName>& inv, 
                 VariableName arr) {
  inv.array_init (arr);
}

// Default implementation
template <typename AbsDomain, typename VariableName >
void array_load (AbsDomain& inv, VariableName lhs, 
                 VariableName arr, VariableName idx) {
}

/////
// (Partial) Specialized versions
/////
template <typename ScalarDomain, typename WeightDomain, typename VariableName>
void array_load (array_graph_domain<ScalarDomain, z_number, VariableName, 
                                    WeightDomain, false>& inv, 
                 VariableName lhs, VariableName arr, VariableName idx) {
   inv.load (lhs, arr, idx);
}

template <typename BaseDomain, typename VariableName>
void array_load (array_smashing<BaseDomain, z_number, VariableName>& inv, 
                 VariableName lhs, VariableName arr, VariableName idx) {
   inv.load (lhs, arr, idx);
}


// Default implementation
template <typename AbsDomain, typename VariableName >
void array_store (AbsDomain& inv, VariableName arr_out, 
                  VariableName arr_in, VariableName idx,
                  typename AbsDomain::linear_expression_t val,
                  bool is_singleton) {
}

/////
// (Partial) Specialized versions
////
template <typename ScalarDomain, typename WeightDomain, typename VariableName>
void array_store (array_graph_domain<ScalarDomain, z_number, VariableName, 
                                     WeightDomain, false>& inv, 
                  VariableName arr_out, VariableName arr_in, VariableName idx,
                  typename ScalarDomain::linear_expression_t val,
                  bool /*is_singleton*/) {
   inv.store (arr_out, arr_in, idx, val);
}

template <typename BaseDomain, typename VariableName>
void array_store (array_smashing<BaseDomain, z_number, VariableName>& inv, 
                  VariableName arr_out, VariableName arr_in, VariableName idx,
                  typename BaseDomain::linear_expression_t val,
                  bool is_singleton) {
   inv.store (arr_out, arr_in, idx, val, is_singleton);
}

}
} // end namespace ikos

#endif 
