#include <crab/domains/wrapped_interval_impl.hpp>
#include <crab/numbers/bignums.hpp>

namespace ikos {
namespace linear_interval_solver_impl {
using z_wrapped_interval_t = crab::domains::wrapped_interval<z_number>;
using q_wrapped_interval_t = crab::domains::wrapped_interval<q_number>;

template <>
z_wrapped_interval_t mk_interval(z_number c,
				 typename crab::wrapint::bitwidth_t w) {
  return z_wrapped_interval_t::mk_winterval(c, w);
}
  
template <>
q_wrapped_interval_t mk_interval(q_number c,
				 typename crab::wrapint::bitwidth_t w) {
  return q_wrapped_interval_t::mk_winterval(c, w);
}

template <>
z_wrapped_interval_t trim_interval(const z_wrapped_interval_t &i,
				   const z_wrapped_interval_t &j) {
  if (i.is_bottom())
    return i;
  // XXX: TODO: gamma(top()) \ gamma(j) can be expressed in a
  //            wrapped interval.
  if (i.is_top())
    return i;
  if (!j.is_singleton())
    return i;

  crab::wrapint k = j.start();
  if (i.start() == k) {
    if (i.is_singleton()) {
      return z_wrapped_interval_t::bottom();
    }
    crab::wrapint k_plus(k);
    ++k_plus;
    z_wrapped_interval_t trimmed_res = z_wrapped_interval_t(k_plus, i.end());
    return trimmed_res;
  } else if (i.end() == k) {
    if (i.is_singleton()) {
      return z_wrapped_interval_t::bottom();
    }
    crab::wrapint k_minus(k);
    --k_minus;
    z_wrapped_interval_t trimmed_res = z_wrapped_interval_t(i.start(), k_minus);
    return trimmed_res;
  } else {
    return i;
  }
}

template <>
q_wrapped_interval_t trim_interval(const q_wrapped_interval_t &i,
				   const q_wrapped_interval_t &j) {
  // No refinement possible for disequations over rational numbers
  return i;
}

template <>
z_wrapped_interval_t lower_half_line(const z_wrapped_interval_t &i,
				     bool is_signed) {
  return i.lower_half_line(is_signed);
}

template <>
q_wrapped_interval_t lower_half_line(const q_wrapped_interval_t &i,
				     bool is_signed) {
  return i.lower_half_line(is_signed);
}

template <>
z_wrapped_interval_t upper_half_line(const z_wrapped_interval_t &i,
				     bool is_signed) {
  return i.upper_half_line(is_signed);
}
  
template <>
q_wrapped_interval_t upper_half_line(const q_wrapped_interval_t &i,
				     bool is_signed) {
  return i.upper_half_line(is_signed);
}

} // namespace linear_interval_solver_impl
} // end namespace ikos

namespace crab {
namespace domains {
// Default instantiations
template class wrapped_interval<ikos::z_number>;
} // end namespace domains
} // end namespace crab
