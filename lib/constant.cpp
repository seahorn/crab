#include <crab/domains/constant_impl.hpp>
#include <crab/numbers/bignums.hpp>

namespace crab {
namespace domains {

template <typename Number>
static constant<Number> sRem(const constant<Number> &x,
                             const constant<Number> &y) {
  if (y.is_constant() && y.get_constant() == 0) {
    return constant<Number>::bottom();
  }
  if (x.is_constant() && y.is_constant()) {
    return constant<Number>(x.get_constant() % y.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
static constant<Number> bitwiseAnd(const constant<Number> &x,
                                   const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    return constant<Number>(x.get_constant() & y.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
static constant<Number> bitwiseOr(const constant<Number> &x,
                                  const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    return constant<Number>(x.get_constant() | y.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
static constant<Number> bitwiseXor(const constant<Number> &x,
                                   const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    return constant<Number>(x.get_constant() ^ y.get_constant());
  } else {
    return constant<Number>::top();
  }
}

template <typename Number>
static constant<Number> bitwiseShl(const constant<Number> &x,
                                   const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    if (y.get_constant() >= 0) {
      return constant<Number>(x.get_constant() << y.get_constant());
    }
  }
  return constant<Number>::top();
}

template <typename Number>
static constant<Number> bitwiseLShr(const constant<Number> &x,
                                    const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    // if get_contant() is non-negative then LShr = AShr.
    if (x.get_constant() >= 0) {
      if (y.get_constant() >= 0) {
        return constant<Number>(x.get_constant() >> y.get_constant());
      }
    }
  }
  return constant<Number>::top();
}

template <typename Number>
static constant<Number> bitwiseAShr(const constant<Number> &x,
                                    const constant<Number> &y) {
  if (x.is_constant() && y.is_constant()) {
    if (y.get_constant() >= 0) {
      return constant<Number>(x.get_constant() >> y.get_constant());
    }
  }
  return constant<Number>::top();
}

/// These operations only instantiated for integer types.

template <>
constant<ikos::z_number>
constant<ikos::z_number>::SRem(const constant<ikos::z_number> &o) const {
  return sRem(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseAnd(const constant<ikos::z_number> &o) const {
  return bitwiseAnd(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseOr(const constant<ikos::z_number> &o) const {
  return bitwiseOr(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseXor(const constant<ikos::z_number> &o) const {
  return bitwiseXor(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseShl(const constant<ikos::z_number> &o) const {
  return bitwiseShl(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseLShr(const constant<ikos::z_number> &o) const {
  return bitwiseLShr(*this, o);
}
template <>
constant<ikos::z_number>
constant<ikos::z_number>::BitwiseAShr(const constant<ikos::z_number> &o) const {
  return bitwiseAShr(*this, o);
}

template <>
constant<uint64_t> constant<uint64_t>::SRem(const constant<uint64_t> &o) const {
  return sRem(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseAnd(const constant<uint64_t> &o) const {
  return bitwiseAnd(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseOr(const constant<uint64_t> &o) const {
  return bitwiseOr(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseXor(const constant<uint64_t> &o) const {
  return bitwiseXor(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseShl(const constant<uint64_t> &o) const {
  return bitwiseShl(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseLShr(const constant<uint64_t> &o) const {
  return bitwiseLShr(*this, o);
}
template <>
constant<uint64_t>
constant<uint64_t>::BitwiseAShr(const constant<uint64_t> &o) const {
  return bitwiseAShr(*this, o);
}

template <>
constant<int64_t> constant<int64_t>::SRem(const constant<int64_t> &o) const {
  return sRem(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseAnd(const constant<int64_t> &o) const {
  return bitwiseAnd(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseOr(const constant<int64_t> &o) const {
  return bitwiseOr(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseXor(const constant<int64_t> &o) const {
  return bitwiseXor(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseShl(const constant<int64_t> &o) const {
  return bitwiseShl(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseLShr(const constant<int64_t> &o) const {
  return bitwiseLShr(*this, o);
}
template <>
constant<int64_t>
constant<int64_t>::BitwiseAShr(const constant<int64_t> &o) const {
  return bitwiseAShr(*this, o);
}

// Default instantiations
template class constant<ikos::z_number>;
template class constant<uint64_t>;
template class constant<int64_t>;
template class constant<ikos::q_number>;

} // namespace domains
} // namespace crab
