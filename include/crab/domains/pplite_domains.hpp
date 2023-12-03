#pragma once

#include <crab/config.h>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {
namespace domains {

namespace pplite_domains {
/*
  The values of this enumeration need to be kept in sync with
    enum pplite::dynamic::Abs_Poly::Kind
*/
//enum pplite_domain_id_t : unsigned int;

/*
  This enumeration needs to be kept in sync with
    enum pplite::dynamic::Abs_Poly::Kind
  The listed enum values are up-to-date with PPLite 0.11.
  Note: even values are for domains; odd values for _STATS variants.
*/
enum pplite_domain_id_t : unsigned int {
  POLY, POLY_STATS,
  B_POLY, B_POLY_STATS,
  F_POLY, F_POLY_STATS,
  U_POLY, U_POLY_STATS,
  UF_POLY, UF_POLY_STATS,
  P_SET, P_SET_STATS,
  FP_SET, FP_SET_STATS
};
} // pplite_domains

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

#if not defined(HAVE_PPLITE) || not defined(HAVE_PPLITE_DOMAINS)

/*
 * Dummy implementation if PPLite not found
 */

#include <crab/domains/dummy_abstract_domain.hpp>

namespace crab {
namespace domains {

template <typename N, typename V, pplite_domains::pplite_domain_id_t K,
          class P = PPLiteDefaultParams<N>>
class pplite_domain final
  : public dummy_abstract_domain<pplite_domain<N,V,K,P>> {
public:
  std::string not_implemented_msg() const override {
#ifndef HAVE_PPLITE
    return "No PPLite. Run cmake with -DCRAB_USE_PPLITE=ON";
#else
    return "No PPLite. Make sure that pplite_domains can be found";
#endif     
  }
};

} // namespace domains
} // namespace crab

#else // defined(HAVE_PPLITE)

#include <crab/domains/pplite/pplite_native_wrapper.hpp>

#endif // HAVE_PPLITE

namespace crab {
namespace domains {
template <typename Number, typename VariableName, pplite_domains::pplite_domain_id_t K,
          class Params>
struct abstract_domain_traits<pplite_domain<Number, VariableName, K, Params>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // namespace domains
} // namespace crab
