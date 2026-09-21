/**
 * Extend abstract domains with very specialized operations
 **/

#pragma once
#include <crab/domains/abstract_domain_operators.hpp>

#include <type_traits>

namespace crab {
namespace domains {

// Only named by the static_asserts below, which reject the generic
// (type-erased) domains. Declaring them is enough for std::is_same, so
// this header does not need to include generic_abstract_domain.hpp --
// which every domain would then pull in just for those asserts.
template <typename Variable> class abstract_domain;
template <typename Variable> class abstract_domain_ref;


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
      } else if (src.get_type().is_fixed_width_integer()) {
	      // This is the "mostly" in the class comment above: the only place where
	      // a bitwidth is not ignored. The bound is sound only if src holds a value
	      // representable in its declared bitwidth, which Crab does not enforce --
	      // no domain reaching this code constrains an iN variable to N bits. When
	      // that does not hold, the bound contradicts the assign above and the
	      // state silently becomes bottom, i.e. the analysis concludes that the
	      // code is unreachable:
	      //
	      //   x:i32 == 5000000000; y:i64 := zext x
	      //   assign gives y == 5000000000, then y <= 2^32-1 makes it bottom.
	      //
	      // Clients that derive iN from a real machine type never build such a
	      // state, and for them the bound is useful precision.
	      //
	      // No lower bound is added here, unlike the Boolean case above: that case
	      // havocs dst first, whereas here dst has just been assigned src, so
	      // dst >= 0 would contradict a negative src.
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
