#pragma once

#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp>
#include <boost/functional/hash_fwd.hpp> // for hash_combine

namespace crab {
  
enum variable_type_kind {
  // scalar types
  BOOL_TYPE,
  INT_TYPE,
  REAL_TYPE,
  REF_TYPE,
  // array types
  ARR_BOOL_TYPE,
  ARR_INT_TYPE,
  ARR_REAL_TYPE,
  // region types: a region can contain either a scalar or an array of
  // a non-reference type
  REG_UNKNOWN_TYPE, // any REG_X unifies with this
  REG_BOOL_TYPE,
  REG_INT_TYPE,
  REG_REAL_TYPE,
  REG_REF_TYPE,
  //
  REG_ARR_BOOL_TYPE,
  REG_ARR_INT_TYPE,
  REG_ARR_REAL_TYPE,

  // unknown type
  UNK_TYPE
};


/*
 * Variable type used by the abstract domains and linear expressions.
 *
 * We try to avoid having a type as a generic class that has a
 * subclass for each subtype.
 */
class variable_type {
  variable_type_kind m_kind;
  /* 
   * Important: the bitwidth (m_bitwidth) is only used if the type is
   * an integer or Boolean. Booleans always have a bitwidth of
   * 1. Reals do not have a bitwidth associated with. For references,
   * we might want to use a bitwidth in the future but for now, we
   * don't. Arrays of integers do not have a bitwidth associated with
   * because arrays can be indexed using different bitwidths so each
   * array store and load must say how many bytes is being
   * accessed. In the future, we might want to use bitwidth for
   * regions of integers but for now, we don't.
  */
  unsigned m_bitwidth; /* in bits */
public:
  variable_type(variable_type_kind kind, unsigned width = 0):
    m_kind(kind), m_bitwidth(width) {
    if (m_kind == INT_TYPE && m_bitwidth <= 1) {
      CRAB_ERROR("Cannot create integer variable without bitwidth");
    }
    
    if (m_kind == BOOL_TYPE) {
      m_bitwidth = 1;
    }
  }
  
  variable_type(const variable_type &o) = default;
  variable_type(variable_type &&o) = default;
  variable_type &operator=(const variable_type &o) = default;
  variable_type &operator=(variable_type &&o) = default;  
  bool operator==(const variable_type &o) const {
    if (m_kind == INT_TYPE) {
      return m_kind == o.m_kind && m_bitwidth == o.m_bitwidth;
    } else {
      return m_kind == o.m_kind;
    }
  }
  bool operator!=(const variable_type &o) const {
    return (!(operator==(o)));
  }

  size_t hash() const {
    size_t res = std::hash<size_t>{}(static_cast<size_t>(m_kind));
    boost::hash_combine(res, std::hash<size_t>{}(static_cast<size_t>(m_bitwidth)));
    return res;
  }

  bool is_typed() const { return m_kind != UNK_TYPE;}
  
  //// scalars
  bool is_bool() const { return m_kind == BOOL_TYPE; }
  bool is_integer() const { return m_kind == INT_TYPE; }
  bool is_integer(unsigned bitwidth) const { return m_kind == INT_TYPE && m_bitwidth == bitwidth; }
  unsigned get_integer_bitwidth() const { assert(m_kind == INT_TYPE); return m_bitwidth;}
  bool is_real() const { return m_kind == REAL_TYPE; }
  bool is_reference() const { return m_kind == REF_TYPE; }
  //// arrays
  bool is_bool_array() const { return m_kind == ARR_BOOL_TYPE; }  
  bool is_integer_array() const { return m_kind == ARR_INT_TYPE; }
  bool is_real_array() const { return m_kind == ARR_REAL_TYPE; }
  bool is_array() const { return is_integer_array() || is_bool_array() || is_real_array();}
  //// regions
  bool is_unknown_region() const { return m_kind == REG_UNKNOWN_TYPE; }
  bool is_bool_region() const { return m_kind == REG_BOOL_TYPE; }
  bool is_integer_region() const { return m_kind == REG_INT_TYPE; }
  bool is_real_region() const { return m_kind == REG_REAL_TYPE; }
  bool is_reference_region() const { return m_kind == REG_REF_TYPE; }
  bool is_bool_array_region() const { return m_kind == REG_ARR_BOOL_TYPE; }
  bool is_int_array_region() const { return m_kind == REG_ARR_INT_TYPE; }
  bool is_real_array_region() const { return m_kind == REG_ARR_REAL_TYPE; }
  bool is_region() const { return is_unknown_region() || is_scalar_region() || is_array_region();}
  bool is_scalar_region() const { return is_bool_region() || is_integer_region() || is_real_region() ||
                                         is_reference_region();}
  bool is_array_region() const { return is_bool_array_region() || is_int_array_region() || is_real_array_region(); }
  
  void write(crab_os &o) const {
    switch (m_kind) {
    case BOOL_TYPE:
      o << "bool";
      break;
    case INT_TYPE:
      o << "int" << m_bitwidth;
      break;
    case REAL_TYPE:
      o << "real";
      break;
    case REF_TYPE:
      o << "ref";
      break;
    case ARR_BOOL_TYPE:
      o << "arr(bool)";
      break;
    case ARR_INT_TYPE:
      o << "arr(int)";
      break;
    case ARR_REAL_TYPE:
      o << "arr(real)";
      break;
    case REG_UNKNOWN_TYPE:
      o << "region(unknown)";
      break;
    case REG_BOOL_TYPE:
      o << "region(bool)";
      break;
    case REG_INT_TYPE:
      o << "region(int)";
      break;
    case REG_REAL_TYPE:
      o << "region(real)";
      break;
    case REG_REF_TYPE:
      o << "region(ref)";
      break;
    case REG_ARR_BOOL_TYPE:
      o << "region(arr(bool))";
      break;
    case REG_ARR_INT_TYPE:
      o << "region(arr(int))";
      break;
    case REG_ARR_REAL_TYPE:
      o << "region(arr(real))";
      break;
    default:
      o << "unknown";
      break;
    }
  }
  
  friend inline crab_os &operator<<(crab_os &o, const variable_type &vt) {
    vt.write(o);
    return o;
  }
};
  
template <typename Number, typename VariableName>
class variable : public indexable {
  // XXX: template parameter Number is required even if the class
  // does not use it.  This allows, e.g., linear_constraint to
  // deduce the kind of Number from constraints like x < y.

public:

  using variable_t = variable<Number, VariableName>;
  using bitwidth_t = unsigned;
  using type_t = variable_type;
  using number_t = Number;
  using varname_t = VariableName;

private:
  VariableName _n;
  variable_type m_ty;

public:
  /* ========== Begin internal API  ============= */
  /* Call this constructor only from abstract domains */  
  explicit variable(const VariableName &n)
    : _n(n), m_ty(crab::UNK_TYPE, 0) {}

  /* Call this constructor only from abstract domains */
  variable(const VariableName &n, variable_type_kind ty_kind)
    : _n(n), m_ty(ty_kind, 0) {}

  /* Call this constructor only from abstract domains */  
  variable(const VariableName &n, variable_type_kind ty_kind, bitwidth_t width)
    : _n(n), m_ty(ty_kind, width) {}
  /* ========== End internal API  =============== */
public:

  variable(const VariableName &n, variable_type ty)
    : _n(n), m_ty(ty) {}

  variable(const variable_t &o) : _n(o._n), m_ty(o.m_ty) {}

  variable(variable_t &&o)
    : _n(std::move(o._n)), m_ty(std::move(o.m_ty)) {}

  variable_t &operator=(const variable_t &o) {
    if (this != &o) {
      _n = o._n;
      m_ty = o.m_ty;
    }
    return *this;
  }

  type_t get_type() const { return m_ty; }

  const VariableName &name() const { return _n; }

  virtual ikos::index_t index() const override { return _n.index(); }

  std::size_t hash() const {
    // casting to size_t may overflow but it shouldn't affect
    // correctness.
    return std::hash<size_t>{}(static_cast<size_t>(_n.index()));
  }

  bool operator==(const variable_t &o) const { return index() == o.index(); }

  bool operator!=(const variable_t &o) const { return (!(operator==(o))); }

  bool operator<(const variable_t &o) const { return index() < o.index(); }

  virtual void write(crab::crab_os &o) const override {
    o << _n;
    CRAB_LOG("crab-print-types", o << ":" << get_type(););
  }

  void dump(crab::crab_os &o) const {
    o << _n << ":" << get_type();
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const variable<Number, VariableName> &v) {
    v.write(o);
    return o;
  }

}; // class variable

/** specialization for boost::hash_combine **/
template <typename Number, typename VariableName>
inline size_t hash_value(const variable<Number, VariableName> &v) {
  return v.hash();
}

} // end namespace crab

/** specialization for std::hash for variables **/
namespace std {
template <typename Number, typename VariableName>
struct hash<crab::variable<Number, VariableName>> {
  using variable_t = crab::variable<Number, VariableName>;
  size_t operator()(const variable_t &v) const { return v.hash(); }
};
} // end namespace std
