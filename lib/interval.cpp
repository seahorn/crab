#include <crab/domains/interval_impl.hpp>
#include <crab/numbers/bignums.hpp>

namespace ikos {

using z_bound_t = bound<z_number>;
using z_interval_t = interval<z_number>;
using q_bound_t = bound<q_number>;
using q_interval_t = interval<q_number>;

template <> q_interval_t q_interval_t::operator/(const q_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    boost::optional<q_number> d = x.singleton();
    if (d && *d == 0) {
      // [_, _] / 0 = _|_
      return bottom();
    } else if (x[0]) {
      boost::optional<q_number> n = singleton();
      if (n && *n == 0) {
        // 0 / [_, _] = 0
        return interval_t(q_number(0));
      } else {
        return top();
      }
    } else {
      bound_t ll = _lb / x._lb;
      bound_t lu = _lb / x._ub;
      bound_t ul = _ub / x._lb;
      bound_t uu = _ub / x._ub;
      return interval_t(bound_t::min(ll, lu, ul, uu),
                        bound_t::max(ll, lu, ul, uu));
    }
  }
}

template <> z_interval_t z_interval_t::operator/(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    // Divisor is a singleton:
    //   the linear interval solver can perform many divisions where
    //   the divisor is a singleton interval. We optimize for this case.
    if (boost::optional<z_number> n = x.singleton()) {
      z_number c = *n;
      if (c == 1) {
        return *this;
      } else if (c > 0) {
        return interval_t(_lb / c, _ub / c);
      } else if (c < 0) {
        return interval_t(_ub / c, _lb / c);
      } else {
      }
    }
    // Divisor is not a singleton
    if (x[0]) {
      z_interval_t l(x._lb, z_bound_t(-1));
      z_interval_t u(z_bound_t(1), x._ub);
      return (operator/(l) | operator/(u));
    } else if (operator[](0)) {
      z_interval_t l(_lb, z_bound_t(-1));
      z_interval_t u(z_bound_t(1), _ub);
      return ((l / x) | (u / x) | z_interval_t(z_number(0)));
    } else {
      // Neither the dividend nor the divisor contains 0
      z_interval_t a =
          (_ub < 0) ? (*this + ((x._ub < 0) ? (x + z_interval_t(z_number(1)))
                                            : (z_interval_t(z_number(1)) - x)))
                    : *this;
      bound_t ll = a._lb / x._lb;
      bound_t lu = a._lb / x._ub;
      bound_t ul = a._ub / x._lb;
      bound_t uu = a._ub / x._ub;
      return interval_t(bound_t::min(ll, lu, ul, uu),
                        bound_t::max(ll, lu, ul, uu));
    }
  }
}

template <> z_interval_t z_interval_t::SRem(const z_interval_t &x) const {
  // note that the sign of the divisor does not matter

  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else if (singleton() && x.singleton()) {
    z_number dividend = *singleton();
    z_number divisor = *x.singleton();

    if (divisor == 0) {
      return bottom();
    }

    return interval_t(dividend % divisor);
  } else if (x.ub().is_finite() && x.lb().is_finite()) {
    z_number max_divisor = max(abs(*x.lb().number()), abs(*x.ub().number()));

    if (max_divisor == 0) {
      return bottom();
    }

    if (lb() < 0) {
      if (ub() > 0) {
        return interval_t(-(max_divisor - 1), max_divisor - 1);
      } else {
        return interval_t(-(max_divisor - 1), 0);
      }
    } else {
      return interval_t(0, max_divisor - 1);
    }
  } else {
    return top();
  }
}

template <> z_interval_t z_interval_t::URem(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else if (singleton() && x.singleton()) {
    z_number dividend = *singleton();
    z_number divisor = *x.singleton();

    if (divisor < 0) {
      return top();
    } else if (divisor == 0) {
      return bottom();
    } else if (dividend < 0) {
      // dividend is treated as an unsigned integer.
      // we would need the size to be more precise
      return interval_t(0, divisor - 1);
    } else {
      return interval_t(dividend % divisor);
    }
  } else if (x.ub().is_finite() && x.lb().is_finite()) {
    z_number max_divisor = *x.ub().number();

    if (x.lb() < 0 || x.ub() < 0) {
      return top();
    } else if (max_divisor == 0) {
      return bottom();
    }

    return interval_t(0, max_divisor - 1);
  } else {
    return top();
  }
}

template <> z_interval_t z_interval_t::And(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    boost::optional<z_number> left_op = singleton();
    boost::optional<z_number> right_op = x.singleton();

    if (left_op && right_op) {
      return interval_t((*left_op) & (*right_op));
    } else if (lb() >= 0 && x.lb() >= 0) {
      return interval_t(0, bound_t::min(ub(), x.ub()));
    } else {
      return top();
    }
  }
}

template <> z_interval_t z_interval_t::Or(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    boost::optional<z_number> left_op = singleton();
    boost::optional<z_number> right_op = x.singleton();

    if (left_op && right_op) {
      return interval_t((*left_op) | (*right_op));
    } else if (lb() >= 0 && x.lb() >= 0) {
      boost::optional<z_number> left_ub = ub().number();
      boost::optional<z_number> right_ub = x.ub().number();

      if (left_ub && right_ub) {
        z_number m = (*left_ub > *right_ub ? *left_ub : *right_ub);
        return interval_t(0, m.fill_ones());
      } else {
        return interval_t(0, bound_t::plus_infinity());
      }
    } else {
      return top();
    }
  }
}

template <> z_interval_t z_interval_t::Xor(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    boost::optional<z_number> left_op = singleton();
    boost::optional<z_number> right_op = x.singleton();

    if (left_op && right_op) {
      return interval_t((*left_op) ^ (*right_op));
    } else {
      return Or(x);
    }
  }
}

template <> z_interval_t z_interval_t::Shl(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    if (boost::optional<z_number> shift = x.singleton()) {
      z_number k = *shift;
      if (k < 0) {
        // CRAB_ERROR("lshr shift operand cannot be negative");
        return top();
      }
      // Some crazy linux drivers generate shl instructions with
      // huge shifts.  We limit the number of times the loop is run
      // to avoid wasting too much time on it.
      if (k <= 128) {
        z_number factor = 1;
        for (int i = 0; k > i; i++) {
          factor *= 2;
        }
        return (*this) * factor;
      }
    }
    return top();
  }
}

template <> z_interval_t z_interval_t::AShr(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    if (boost::optional<z_number> shift = x.singleton()) {
      z_number k = *shift;
      if (k < 0) {
        return top();
      }
      // Some crazy linux drivers generate ashr instructions with
      // huge shifts.  We limit the number of times the loop is run
      // to avoid wasting too much time on it.
      if (k <= 128) {
        z_number factor = 1;
        for (int i = 0; k > i; i++) {
          factor *= 2;
        }
        return (*this) / factor;
      }
    }
    return top();
  }
}

template <> z_interval_t z_interval_t::LShr(const z_interval_t &x) const {
  if (is_bottom() || x.is_bottom()) {
    return bottom();
  } else {
    if (boost::optional<z_number> shift = x.singleton()) {
      z_number k = *shift;
      if (k < 0) {
        return top();
      }
      if (lb() >= 0 && ub().is_finite()) {
        z_number lb = *this->lb().number();
        z_number ub = *this->ub().number();
        return z_interval_t(lb >> k, ub >> k);
      }
    }
    return this->top();
  }
}

namespace bounds_impl {
// Conversion between z_bound_t and q_bound_t
void convert_bounds(z_bound_t b1, z_bound_t &b2) { std::swap(b1, b2); }
void convert_bounds(q_bound_t b1, q_bound_t &b2) { std::swap(b1, b2); }
void convert_bounds(z_bound_t b1, q_bound_t &b2) {
  if (b1.is_plus_infinity())
    b2 = q_bound_t::plus_infinity();
  else if (b1.is_minus_infinity())
    b2 = q_bound_t::minus_infinity();
  else
    b2 = q_bound_t(q_number(*b1.number()));
}
void convert_bounds(q_bound_t b1, z_bound_t &b2) {
  if (b1.is_plus_infinity())
    b2 = z_bound_t::plus_infinity();
  else if (b1.is_minus_infinity())
    b2 = z_bound_t::minus_infinity();
  else
    b2 = z_bound_t((*(b1.number())).round_to_lower());
}
} // namespace bounds_impl

namespace linear_interval_solver_impl {
template <>
z_interval_t trim_interval(const z_interval_t &i, const z_interval_t &j) {
  if (boost::optional<z_number> c = j.singleton()) {
    if (i.lb() == *c) {
      return z_interval_t(*c + 1, i.ub());
    } else if (i.ub() == *c) {
      return z_interval_t(i.lb(), *c - 1);
    } else {
    }
  }
  return i;
}

template <>
q_interval_t trim_interval(const q_interval_t &i,
                           const q_interval_t & /* j */) {
  // No refinement possible for disequations over rational numbers
  return i;
}

template <>
z_interval_t lower_half_line(const z_interval_t &i, bool /*is_signed*/) {
  return i.lower_half_line();
}

template <>
q_interval_t lower_half_line(const q_interval_t &i, bool /*is_signed*/) {
  return i.lower_half_line();
}

template <>
z_interval_t upper_half_line(const z_interval_t &i, bool /*is_signed*/) {
  return i.upper_half_line();
}

template <>
q_interval_t upper_half_line(const q_interval_t &i, bool /*is_signed*/) {
  return i.upper_half_line();
}
} // namespace linear_interval_solver_impl

// Default instantiations
template class bound<z_number>;
template class bound<q_number>;
template class interval<z_number>;
template class interval<q_number>;

} // end namespace ikos
