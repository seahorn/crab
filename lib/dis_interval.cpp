#include <crab/domains/dis_interval_impl.hpp>
#include <crab/numbers/bignums.hpp>

namespace ikos {
namespace linear_interval_solver_impl {

using z_dis_interval_t = crab::domains::dis_interval<z_number>;
using q_dis_interval_t = crab::domains::dis_interval<q_number>;

template <>
z_dis_interval_t trim_interval(const z_dis_interval_t &x,
                               const z_dis_interval_t &y) {
  using z_interval_t = ikos::interval<z_number>;

  if (x.is_bottom())
    return x;

  boost::optional<z_number> s = y.singleton();
  if (!s)
    return x;
  z_number c = *s;

  z_dis_interval_t res = z_dis_interval_t::bottom();
  if (x.is_top()) {
    res = res | z_dis_interval_t(z_interval_t(c - 1).lower_half_line());
    res = res | z_dis_interval_t(z_interval_t(c + 1).upper_half_line());
  } else {
    for (auto i : boost::make_iterator_range(x.begin(), x.end())) {
      if (!(z_interval_t(c) <= i)) {
        res = res | i;
        continue;
      }

      if (i.lb() == c) {
        res = res | z_dis_interval_t(z_interval_t(c + 1, i.ub()));
      } else if (i.ub() == c) {
        res = res | z_dis_interval_t(z_interval_t(i.lb(), c - 1));
      } else {
        res = res | z_dis_interval_t(z_interval_t(i.lb(), c - 1));
        res = res | z_dis_interval_t(z_interval_t(c + 1, i.ub()));
      }
    }
  }

  return res;
}

template <>
q_dis_interval_t trim_interval(const q_dis_interval_t &i,
                               const q_dis_interval_t & /* j */) {
  // No refinement possible for disequations over rational numbers
  return i;
}

template <>
z_dis_interval_t lower_half_line(const z_dis_interval_t &i,
                                 bool /*is_signed*/) {
  return i.lower_half_line();
}

template <>
q_dis_interval_t lower_half_line(const q_dis_interval_t &i,
                                 bool /*is_signed*/) {
  return i.lower_half_line();
}

template <>
z_dis_interval_t upper_half_line(const z_dis_interval_t &i,
                                 bool /*is_signed*/) {
  return i.upper_half_line();
}

template <>
q_dis_interval_t upper_half_line(const q_dis_interval_t &i,
                                 bool /*is_signed*/) {
  return i.upper_half_line();
}

} // namespace linear_interval_solver_impl
} // end namespace ikos

namespace crab {
namespace domains {
// Default instantiations of dis_interval class
template class dis_interval<ikos::z_number>;
template class dis_interval<ikos::q_number>;
} // end namespace domains
} // end namespace crab
