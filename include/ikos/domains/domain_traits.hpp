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

template <typename AbsNumDomain>
void normalize(AbsNumDomain& inv); 

template <typename AbsNumDomain, typename Iterator>
void forget(AbsNumDomain& inv, Iterator begin, Iterator end); 

template <typename AbsDomain, typename VariableName>
void expand (AbsDomain& inv, VariableName x, VariableName new_x);

template <typename AbsDomain, typename VariableName>
void array_init (AbsDomain& inv, VariableName arr); 

template <typename AbsDomain, typename VariableName>
void array_load (AbsDomain& inv, VariableName lhs, 
                 VariableName arr, VariableName idx);

template <typename AbsDomain, typename VariableName>
void array_store (AbsDomain& inv, VariableName arr, 
                  VariableName idx, typename AbsDomain::linear_expression_t val,
                  bool is_singleton); 
}
} // end namespace ikos

#endif 
