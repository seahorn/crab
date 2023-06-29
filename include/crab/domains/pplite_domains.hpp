#pragma once

#include <crab/config.h>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {
namespace domains {

/*
  The values of this enumeration need to be kept in sync with
    enum pplite::dynamic::Abs_Poly::Kind
*/
enum pplite_domain_id_t : unsigned int;

template <typename Number> class PPLiteDefaultParams {
public:
  // use integers with truncation rounding
  enum { use_integers = 1 };
  enum { implement_inter_transformers = 0 };
};

template <> class PPLiteDefaultParams<ikos::q_number> {
public:
  // use reals
  enum { use_integers = 0 };
  enum { implement_inter_transformers = 0 };
};

} // namespace domains
} // namespace crab

#ifndef HAVE_PPLITE

/*
 * Dummy implementation if PPLite not found
 */

#include <crab/domains/dummy_abstract_domain.hpp>

namespace crab {
namespace domains {

template <typename N, typename V, pplite_domain_id_t K,
          class P = PPLiteDefaultParams<N>>
class pplite_domain final
  : public dummy_abstract_domain<pplite_domain<N,V,K,P>> {
public:
  std::string not_implemented_msg() const override {
    return "No PPLite. Run cmake with -DCRAB_USE_PPLITE=ON";
  }
};

} // namespace domains
} // namespace crab

#else // defined(HAVE_PPLITE)

#include <crab/domains/pplite/pplite_native_wrapper.hpp>

#endif // HAVE_PPLITE

namespace crab {
namespace domains {
template <typename Number, typename VariableName, pplite_domain_id_t K,
          class Params>
struct abstract_domain_traits<pplite_domain<Number, VariableName, K, Params>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // namespace domains
} // namespace crab
