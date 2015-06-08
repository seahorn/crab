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
void array_load (AbsDomain& inv, VariableName lhs, VariableName arr, VariableName idx) {
}

// (Partial) Specialized version
template <typename VariableName>
void array_load (array_graph_domain<DBM <z_number, VariableName>, 
                                    z_number, 
                                    VariableName, 
                                    interval_domain <z_number, VariableName>, false>& inv, 
                 VariableName lhs, VariableName arr, VariableName idx) {
   inv.load (lhs, arr, idx);
}


// Default implementation
template <typename AbsDomain, typename VariableName >
void array_store (AbsDomain& inv, VariableName arr_out, VariableName arr_in, VariableName idx,
                  typename AbsDomain::linear_expression_t val) {
}

// (Partial) Specialized version
template <typename VariableName>
void array_store (array_graph_domain<DBM <z_number, VariableName>, 
                                     z_number, 
                                     VariableName, 
                                     interval_domain <z_number, VariableName>, false>& inv, 
                  VariableName arr_out, VariableName arr_in, VariableName idx,
                  typename interval_domain <z_number, VariableName>::linear_expression_t val) {
   inv.store (arr_out, arr_in, idx, val);
}

}
} // end namespace ikos

#endif 
