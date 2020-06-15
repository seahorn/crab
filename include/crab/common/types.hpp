#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <functional>
#include <iosfwd>
#include <memory>

/* Basic type definitions */

namespace ikos {
// Numerical type for indexed objects
typedef uint64_t index_t;
} // end namespace ikos

namespace crab {
// XXX: we try to avoid having a type as a generic class that has a
// subclass for each subtype.
enum variable_type {
  BOOL_TYPE,
  INT_TYPE,
  REAL_TYPE,
  REF_TYPE,
  
  ARR_BOOL_TYPE,
  ARR_INT_TYPE,
  ARR_REAL_TYPE,
  
  UNK_TYPE  
};

inline crab_os &operator<<(crab_os &o, variable_type t) {
  switch (t) {
  case BOOL_TYPE:
    o << "bool";
    break;
  case INT_TYPE:
    o << "int";
    break;
  case REAL_TYPE:
    o << "real";
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
  case REF_TYPE:
    o << "ref";
    break;    
  default:
    o << "unknown";
    break;
  }
  return o;
}

typedef enum {
  BINOP_ADD,
  BINOP_SUB,
  BINOP_MUL,
  BINOP_SDIV,
  BINOP_UDIV,
  BINOP_SREM,
  BINOP_UREM,
  BINOP_AND,
  BINOP_OR,
  BINOP_XOR,
  BINOP_SHL,
  BINOP_LSHR,
  BINOP_ASHR,
  BINOP_FUNCTION
} binary_operation_t;

typedef enum { BINOP_BAND, BINOP_BOR, BINOP_BXOR } bool_binary_operation_t;

typedef enum { CAST_TRUNC, CAST_SEXT, CAST_ZEXT } cast_operation_t;

inline crab::crab_os &operator<<(crab::crab_os &o, binary_operation_t op) {
  switch (op) {
  case BINOP_ADD:
    o << "+";
    break;
  case BINOP_SUB:
    o << "-";
    break;
  case BINOP_MUL:
    o << "*";
    break;
  case BINOP_SDIV:
    o << "/";
    break;
  case BINOP_UDIV:
    o << "/_u";
    break;
  case BINOP_SREM:
    o << "%";
    break;
  case BINOP_UREM:
    o << "%_u";
    break;
  case BINOP_AND:
    o << "&";
    break;
  case BINOP_OR:
    o << "|";
    break;
  case BINOP_XOR:
    o << "^";
    break;
  case BINOP_SHL:
    o << "<<";
    break;
  case BINOP_LSHR:
    o << ">>_l";
    break;
  case BINOP_ASHR:
    o << ">>_r";
    break;
  case BINOP_FUNCTION:
    o << "uf";
    break;
  default:
    CRAB_ERROR("unexpected binary operation ", op);
  }
  return o;
}

inline crab::crab_os &operator<<(crab::crab_os &o, bool_binary_operation_t op) {
  switch (op) {
  case BINOP_BAND:
    o << "&";
    break;
  case BINOP_BOR:
    o << "|";
    break;
  case BINOP_BXOR:
    o << "^";
    break;
  default:
    CRAB_ERROR("unexpected boolean binary operation ", op);
  }
  return o;
}

inline crab::crab_os &operator<<(crab::crab_os &o, cast_operation_t op) {
  switch (op) {
  case CAST_TRUNC:
    o << "trunc";
    break;
  case CAST_SEXT:
    o << "sext";
    break;
  case CAST_ZEXT:
    o << "zext";
    break;
  default:
    CRAB_ERROR("unexpected cast operation", op);
  }
  return o;
}

template <typename T> inline boost::optional<T> conv_op(binary_operation_t op);

template <typename T>
inline boost::optional<T> conv_op(bool_binary_operation_t op);

template <typename T> inline boost::optional<T> conv_op(cast_operation_t op);

} // end namespace crab

namespace ikos {
// Interface for writeable objects
class writeable {
public:
  virtual void write(crab::crab_os &o) = 0;
  virtual ~writeable() {}
}; // class writeable

inline crab::crab_os &operator<<(crab::crab_os &o, writeable &x) {
  x.write(o);
  return o;
}

// Container for typed variables used by the crab abstract domains
// and linear_constraints.
template <typename Number, typename VariableName> class variable {
  // XXX: template parameter Number is required even if the class
  // does not use it.  This allows, e.g., linear_constraint to
  // deduce the kind of Number from constraints like x < y.

public:
  typedef variable<Number, VariableName> variable_t;
  typedef unsigned bitwidth_t;
  typedef crab::variable_type type_t;
  typedef Number number_t;
  typedef VariableName varname_t;

private:
  VariableName _n;
  type_t _type;
  bitwidth_t _width;

public:
  /**
   * DO NOT USE this constructor to create a CFG since all CFG
   * statements must be strongly typed.  This constructor is
   * intended to be used only abstract domains to generate temporary
   * variables.
   **/
  explicit variable(const VariableName &n)
      : _n(n), _type(crab::UNK_TYPE), _width(0) {}

public:
  variable(const VariableName &n, type_t type)
      : _n(n), _type(type), _width(0) {}

  variable(const VariableName &n, type_t type, bitwidth_t width)
      : _n(n), _type(type), _width(width) {}

  variable(const variable_t &o) : _n(o._n), _type(o._type), _width(o._width) {}

  variable(variable_t &&o)
      : _n(std::move(o._n)), _type(std::move(o._type)),
        _width(std::move(o._width)) {}

  variable_t &operator=(const variable_t &o) {
    if (this != &o) {
      _n = o._n;
      _type = o._type;
      _width = o._width;
    }
    return *this;
  }

  bool is_typed() const { return _type != crab::UNK_TYPE; }

  bool is_int_type() const { return _type == crab::INT_TYPE; }

  bool is_bool_type() const { return _type == crab::BOOL_TYPE; }

  bool is_real_type() const { return _type == crab::REAL_TYPE; }  

  bool is_ref_type() const { return _type == crab::REF_TYPE; }
  
  bool is_array_type() const {
    return _type >= crab::ARR_BOOL_TYPE && _type <= crab::ARR_REAL_TYPE;
  }

  type_t get_type() const { return _type; }

  bool has_bitwidth() const { return _width > 0; }

  bitwidth_t get_bitwidth() const { return _width; }

  const VariableName &name() const { return _n; }

  // Cannot be const because from VariableName we might want to
  // access to its variable factory to create e.g., new
  // VariableName's.
  VariableName &name() { return _n; }

  ikos::index_t index() const { return _n.index(); }

  std::size_t hash() const {
    // casting to size_t may overflow but it shouldn't affect
    // correctness.
    return std::hash<size_t>{}(static_cast<size_t>(_n.index()));
  }

  bool operator==(const variable_t &o) const {
    return index() == o.index();
  }

  bool operator!=(const variable_t &o) const { return (!(operator==(o))); }

  bool operator<(const variable_t &o) const {
    return index() < o.index();
  }

  void write(crab::crab_os &o) const {
    o << _n;
    CRAB_LOG("crab-print-types",
	     o << ":" << get_type();
	     switch (get_type()) {
	     case crab::INT_TYPE:
	     case crab::ARR_INT_TYPE:
	       o << ":" << get_bitwidth();
	       break;
	     default:;
	     });
  }

  void dump(crab::crab_os &o) const {
    o << _n << ":" << get_type() << ":" << get_bitwidth();
  }

}; // class variable

/* used by boost::hash_combine */
template <typename Number, typename VariableName>
inline size_t hash_value(const variable<Number, VariableName> &v) {
  return v.hash();
}

template <typename Number, typename VariableName>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const variable<Number, VariableName> &v) {
  v.write(o);
  return o;
}
} // end namespace ikos

namespace crab {
namespace domains {
using namespace ikos;
}
} // namespace crab

#ifdef BOOST_NO_EXCEPTIONS
namespace boost {
template <class E> void throw_exception(E const &e) { std::exit(1); }
} // namespace boost
#endif

/** specialization for std::hash for variables **/
namespace std {
template <typename Number, typename VariableName>
struct hash<ikos::variable<Number, VariableName>> {
  using variable_t = ikos::variable<Number, VariableName>;
  size_t operator()(const variable_t &v) const { return v.hash(); }
};
} // namespace std
