/**
 * Extend abstract domains with very specialized operations
 **/

#pragma once
#include <crab/domains/abstract_domain_operators.hpp>

#include <boost/optional.hpp>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>

namespace crab {
namespace domains {

// Only named by the static_asserts below, which reject the generic
// (type-erased) domains. Declaring them is enough for std::is_same, so
// this header does not need to include generic_abstract_domain.hpp --
// which every domain would then pull in just for those asserts.
template <typename Variable> class abstract_domain;
template <typename Variable> class abstract_domain_ref;


// Perform constraint simplifications depending on the abstract domain
template <typename Domain> class constraint_simp_domain_traits {
public:
  static_assert(
      !std::is_same<Domain, abstract_domain<typename Domain::variable_t>>::value,
      "constraint_simp_domain_traits not supported for generic domain");
  static_assert(
      !std::is_same<Domain,
                   abstract_domain_ref<typename Domain::variable_t>>::value,
      "constraint_simp_domain_traits not supported for generic domain");
  
  using number_t = typename Domain::number_t;
  using variable_t = typename Domain::variable_t;  
  using linear_expression_t = typename Domain::linear_expression_t;  
  using linear_constraint_t = typename Domain::linear_constraint_t;
  using linear_constraint_system_t = typename Domain::linear_constraint_system_t; 

  // Convert a disequality into a strict inequality:
  // - if cst is x!=y and abs_val |= x <= y then add x < y
  // - if cst is x!=y and abs_val |= x >= y then add x > y
  static void lower_disequality(const Domain &abs_val,
				const linear_constraint_t &cst,
				linear_constraint_system_t &out_csts) {

    // TODO: we could use abs_val to infer more disequalities from cst.
    auto get_binary_operands = [](const linear_constraint_t &c) -> 
      boost::optional<std::pair<variable_t, variable_t>> {
	if (c.is_disequation()) {
	  if (c.size() == 2 && c.constant() == 0) {
	    auto it = c.begin();
	    auto nx = it->first;
	    auto vx = it->second;
	    ++it;
	    assert(it != c.end());
	    auto ny = it->first;
	    auto vy = it->second;
	    if (nx == (ny * -1)) {
	      return std::make_pair(vx, vy);
	    }
	  }
	}
      return boost::none;
    };
    
    if (auto pair = get_binary_operands(cst)) {
      variable_t x = (*pair).first;
      variable_t y = (*pair).second;
      linear_constraint_t x_le_y(x <= y);
      linear_constraint_t x_lt_y(linear_expression_t(x) < linear_expression_t(y));
      if (abs_val.entails(x_le_y)) {
	out_csts += x_lt_y;
      } else {
	linear_constraint_t x_ge_y(x >= y);
	linear_constraint_t x_gt_y(linear_expression_t(x) > linear_expression_t(y));	    
	if (abs_val.entails(x_ge_y)) {
	  out_csts += x_gt_y;	      
	} 
      }
    }
  }  
};

// Default implementation of integer cast instructions:
// signed-extension, zero-extension and truncation.
// 
// This implementation assumes that the abstract domain Domain models
// integers as mathematical integers and hence, bit-widths are mostly
// ignored. Moreover, it assumes that Domain does not model Booleans.
template<typename Domain> class int_cast_domain_traits {
public:
  using number_t = typename Domain::number_t;
  using variable_t = typename Domain::variable_t;

  static_assert(
      !std::is_same<Domain, abstract_domain<typename Domain::variable_t>>::value,
      "int_cast_domain_traits not supported for generic domain");
  static_assert(
      !std::is_same<Domain,
                   abstract_domain_ref<typename Domain::variable_t>>::value,
      "int_cast_domain_traits not supported for generic domain");
  
  static void apply(Domain &dom, int_conv_operation_t op,
		    const variable_t &dst, const variable_t &src) {

    if (!(dst.get_type().is_bool() || src.get_type().is_bool())) {
      dom.assign(dst, src);
    } else {
      dom -= dst;
    }
    
    /// Refine dst based on src's type
    if (op == crab::domains::OP_ZEXT) {
      if (src.get_type().is_bool()) {
	dom += (dst >= 0);
	dom += (dst <= 1);
      } else if (src.get_type().is_integer()) {
	unsigned bitwidth = src.get_type().get_integer_bitwidth();
	number_t upper_bound = (number_t(1) << number_t(bitwidth)) - number_t(1);
	dom += (dst <= number_t(upper_bound));
      }
    }
  }
};

// DEPRECATED: the boxes domain has global state that should be reset
// by the client if multiple crab instances will be run. This is just
// a temporary hack. The proper solution is to avoid global state.
template<typename Domain> class special_domain_traits {
public:
  static void clear_global_state(void) {}
};

} // end namespace domains
} // end namespace crab
