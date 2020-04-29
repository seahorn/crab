#pragma once

/*
   A simple flat 3-valued boolean domain and a reduced product of this
   flat bool domain with an arbitrary numerical domain
*/

#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/separate_domains.hpp>

namespace crab {

namespace domains {

class boolean_value {
  /*
            Top
            / \
           /   \
        True  False
           \   /
            \ /
           Bottom
  */
  typedef enum { False = 0x0, True = 0x1, Bottom = 0x2, Top = 0x3 } kind_t;

  kind_t _value;

  boolean_value(kind_t v) : _value(v){};

public:
  boolean_value() : _value(Top) {}

  static boolean_value bottom() { return boolean_value(Bottom); }

  static boolean_value top() { return boolean_value(Top); }

  static boolean_value get_true() { return boolean_value(True); }

  static boolean_value get_false() { return boolean_value(False); }

  boolean_value(const boolean_value &other) : _value(other._value) {}

  boolean_value &operator=(const boolean_value &other) {
    if (this != &other)
      _value = other._value;
    return *this;
  }

  bool is_bottom() { return (_value == Bottom); }

  bool is_top() { return (_value == Top); }

  bool is_true() const { return (_value == True); }

  bool is_false() const { return (_value == False); }

  bool operator<=(boolean_value other) {

    if (_value == Bottom || other._value == Top)
      return true;
    else if (_value == Top)
      return (other._value == Top);
    else if (_value == True)
      return ((other._value == True) || (other._value == Top));
    else if (_value == False)
      return ((other._value == False) || (other._value == Top));
    // this should be unreachable
    return false;
  }

  bool operator==(boolean_value other) { return (_value == other._value); }

  boolean_value operator|(boolean_value other) {
    if (is_bottom())
      return other;
    if (other.is_bottom())
      return *this;
    if (is_top() || other.is_top())
      return top();
    if (_value == other._value)
      return *this;

    // othewise true | false or false | true ==> top
    return top();
  }

  boolean_value operator&(boolean_value other) {
    if (is_bottom())
      return *this;
    if (other.is_bottom())
      return other;
    if (is_top())
      return other;
    if (other.is_top())
      return *this;
    if (_value == other._value)
      return *this;

    // othewise true & false or false & true ==> bottom
    return bottom();
  }

  // the lattice satisfy ACC so join is the widening
  boolean_value operator||(boolean_value other) {
    return this->operator|(other);
  }

  // the lattice satisfy DCC so meet is the narrowing
  boolean_value operator&&(boolean_value other) {
    return this->operator&(other);
  }

  // Boolean operations
  /*
            And  Or  X0r
       0 0   0   0    0
       0 1   0   1    1
       1 0   0   1    1
       1 1   1   1    0
       0 *   0   *    *
       * 0   0   *    *
       1 *   *   1    *
       * 1   *   1    *
       * *   *   *    *
  */

  boolean_value And(boolean_value other) {
    if (is_bottom() || other.is_bottom())
      return bottom();

    if (!is_top() && !other.is_top())
      return boolean_value(static_cast<kind_t>(static_cast<int>(_value) &
                                               static_cast<int>(other._value)));

    int x = static_cast<int>(_value);
    int y = static_cast<int>(other._value);
    if (x == 0 || y == 0)
      return get_false();
    else
      return top();
  }

  boolean_value Or(boolean_value other) {
    if (is_bottom() || other.is_bottom())
      return bottom();

    if (!is_top() && !other.is_top())
      return boolean_value(static_cast<kind_t>(static_cast<int>(_value) |
                                               static_cast<int>(other._value)));

    int x = static_cast<int>(_value);
    int y = static_cast<int>(other._value);
    if (x == 1 || y == 1)
      return get_true();
    else
      return top();
  }

  boolean_value Xor(boolean_value other) {
    if (is_bottom() || other.is_bottom())
      return bottom();

    if (!is_top() && !other.is_top())
      return boolean_value(static_cast<kind_t>(static_cast<int>(_value) ^
                                               static_cast<int>(other._value)));
    else
      return top();
  }

  boolean_value Negate() {
    if (is_bottom())
      return bottom();
    if (_value == True)
      return get_false();
    if (_value == False)
      return get_true();
    return top();
  }

  void write(crab_os &o) const {
    switch (_value) {
    case Bottom:
      o << "_|_";
      break;
    case Top:
      o << "*";
      break;
    case True:
      o << "true";
      break;
    default: /*False*/
      o << "false";
    }
  }

}; // end class boolean_value

inline crab_os &operator<<(crab_os &o, const boolean_value &v) {
  v.write(o);
  return o;
}

// A simple abstract domain for booleans
template <typename Number, typename VariableName>
class flat_boolean_domain final
    : public abstract_domain<flat_boolean_domain<Number, VariableName>> {
  typedef flat_boolean_domain<Number, VariableName> flat_boolean_domain_t;

public:
  typedef abstract_domain<flat_boolean_domain_t> abstract_domain_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;

  typedef interval<number_t> interval_t;
  typedef boolean_value bool_t;
  typedef separate_domain<variable_t, boolean_value> separate_domain_t;
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  flat_boolean_domain(separate_domain_t env) : _env(env) {}

public:
  void set_to_top() {
    flat_boolean_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    flat_boolean_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_domain() : _env(separate_domain_t::top()) {}

  flat_boolean_domain(const flat_boolean_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  flat_boolean_domain_t &operator=(const flat_boolean_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o)
      _env = o._env;
    return *this;
  }

  iterator begin() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return _env.begin();
  }

  iterator end() {
    if (is_bottom())
      CRAB_ERROR("Cannot return iterator from bottom");
    return _env.end();
  }

  bool is_bottom() { return _env.is_bottom(); }

  bool is_top() { return _env.is_top(); }

  bool operator<=(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    return (_env <= o._env);
  }

  flat_boolean_domain_t operator|(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    flat_boolean_domain_t res(_env | o._env);
    CRAB_LOG("flat-boolean", crab::outs() << "After join " << *this << " and "
                                          << o << "=" << res << "\n";);
    return res;
  }

  void operator|=(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    CRAB_LOG("flat-boolean",
             crab::outs() << "After join " << *this << " and " << o << "=");
    _env = _env | o._env;
    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  flat_boolean_domain_t operator&(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");

    flat_boolean_domain_t res(_env & o._env);
    CRAB_LOG("flat-boolean", crab::outs() << "After meet " << *this << " and "
                                          << o << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t operator||(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    flat_boolean_domain_t res(_env || o._env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t
  widening_thresholds(flat_boolean_domain_t o,
                      const iterators::thresholds<number_t> &) {

    flat_boolean_domain_t res(_env || o._env);
    CRAB_LOG("flat-boolean", crab::outs()
                                 << "After widening " << *this << " and " << o
                                 << "=" << res << "\n");
    return res;
  }

  flat_boolean_domain_t operator&&(flat_boolean_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    return (_env && o._env);
  }

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
    if (!is_bottom())
      _env -= v;
  }

  // flat_boolean_domains_api

  // XXX: the flat boolean domain cannot reason about linear
  // constraints so we assign top to x.
  void assign_bool_cst(variable_t x, linear_constraint_t cst) {
    _env -= x;
    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << x << ":=" << bx << "\n");
  }

  void assign_bool_var(variable_t x, variable_t y, bool is_not_y) {
    crab::CrabStats::count(getDomainName() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_var");
    _env.set(x, (is_not_y ? _env[y].Negate() : _env[y]));
    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << "After " << x << ":=";
             if (is_not_y) crab::outs() << "not(" << y << ")";
             else crab::outs() << y;
             crab::outs() << " --->" << x << "=" << bx << "\n");
  }

  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply_binary_bool");

    switch (op) {
    case OP_BAND:
      _env.set(x, _env[y].And(_env[z]));
      break;
    case OP_BOR:
      _env.set(x, _env[y].Or(_env[z]));
      break;
    case OP_BXOR:
      _env.set(x, _env[y].Xor(_env[z]));
      break;
    default:
      CRAB_ERROR("Unknown boolean operator");
    }

    CRAB_LOG("flat-boolean", auto bx = _env[x];
             crab::outs() << "After " << x << ":=" << y << " " << op << " " << z
                          << " --->" << x << "=" << bx << "\n");
  }

  void assume_bool(variable_t x, bool is_negated) {
    crab::CrabStats::count(getDomainName() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".assume_bool");

    if (!is_negated)
      _env.set(x, _env[x] & boolean_value::get_true());
    else
      _env.set(x, _env[x] & boolean_value::get_false());

    CRAB_LOG("flat-boolean", auto bx = _env[x];
             if (!is_negated) crab::outs()
             << "After assume(" << x << ") --> " << x << "=" << bx << "\n";
             else crab::outs() << "After assume(not(" << x << ")) --> " << x
                               << "=" << bx << "\n";);
  }

  // XXX: these methods are not actually part of boolean_operators
  // api but they are used by flat_boolean_numerical_domain and
  // domain_traits.

  void set_bool(variable_t x, boolean_value v) { _env.set(x, v); }

  boolean_value get_bool(variable_t x) { return _env[x]; }

  // backward boolean operators
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                flat_boolean_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign_bool_cst");
    if (is_bottom())
      return;

    /* nothing to do: flat_boolean_domain ignores this */
    _env -= lhs;
  }

  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                flat_boolean_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign_bool_var");

    if (is_bottom())
      return;
    /** TODO  **/
    /*
       assume(lhs == rhs);
       assume(lhs == not(rhs))
    */
    _env -= lhs;
    *this = *this & inv;
  }

  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  flat_boolean_domain_t inv) {
    crab::CrabStats::count(getDomainName() +
                           ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".backward_apply_binary_bool");

    if (is_bottom())
      return;

    /** TODO **/
    /*
       if x is true and op=AND then y=true and z=true
       if x is false and op=OR then y=false and z=false
    */
    _env -= x;
    *this = *this & inv;
  }

  /*
     Begin unimplemented operations

     flat_boolean_domain implements only boolean operations.  The
     implementation of the rest of operations is empty because they
     should never be called.
  */

  // numerical_domains_api
  // XXX: needed for making a reduced product with a numerical domain
  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {}
  void apply(operation_t op, variable_t x, variable_t y, number_t k) {}
  void assign(variable_t x, linear_expression_t e) {}
  void backward_assign(variable_t x, linear_expression_t e,
                       flat_boolean_domain_t invariant) {}
  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      flat_boolean_domain_t invariant) {}
  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      flat_boolean_domain_t invariant) {}
  void operator+=(linear_constraint_system_t csts) {}
  void operator+=(linear_constraint_t cst) {}
  // not part of the numerical_domains api but it should be
  void set(variable_t x, interval_t intv) {}
  interval_t operator[](variable_t x) { return interval_t::top(); }

  // int_cast_operators_api and bitwise_operators_api
  // XXX: needed for making a reduced product with a numerical domain
  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {}
  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
  }
  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t z) {}

  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           flat_boolean_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           flat_boolean_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            flat_boolean_domain_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            flat_boolean_domain_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  flat_boolean_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  flat_boolean_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             flat_boolean_domain_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  static std::string getDomainName() { return "Boolean"; }

  void write(crab_os &o) {
    crab::CrabStats::count(getDomainName() + ".count.write");
    crab::ScopedCrabStats __st__(getDomainName() + ".write");

    _env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    if (is_bottom())
      return linear_constraint_t::get_false();

    if (is_top())
      return linear_constraint_t::get_true();

    linear_constraint_system_t res;
    for (auto kv : _env) {
      variable_t v(kv.first);
      if (kv.second.is_true()) {
        res += linear_constraint_t(v == number_t(1));
      } else if (kv.second.is_false()) {
        res += linear_constraint_t(v == number_t(0));
      } else {
        res += linear_constraint_t(v >= number_t(0));
        res += linear_constraint_t(v >= number_t(1));
      }
    }
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    _env.rename(from, to);    
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  flat_boolean_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }

    for (variable_t v : variables) {
      this->operator-=(v);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    flat_boolean_domain_t res;
    for (variable_t v : variables) {
      res.set_bool(v, get_bool(v));
    }
    std::swap(*this, res);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set_bool(new_x, get_bool(x));
  }

  void normalize() {}

  void minimize() {}

}; // class flat_boolean_domain

template <typename Number, typename VariableName>
struct abstract_domain_traits<flat_boolean_domain<Number, VariableName>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

// Simple reduced product of the flat boolean domain with an
// arbitrary numerical domain.
// The reduction happens in two situations:
//    (1) when bvar := linear_constraint from numerical to boolean
//    (2) when assume_bool(bvar) from boolean to numerical
// The step (2) is quite weak.
//
// For instance, code like this won't trigger any propagation so
// we won't know the value of x:
//
//   havoc(x);
//   b1 = (x >= 0);
//   y  = zext b1 to i32
//   assume (y >= 1);
//
// Instead, the code should be rewritten like this:
//   havoc(x);
//   b1 = (x >= 0);
//   assume_bool(b1);
//
template <typename NumDom>
class flat_boolean_numerical_domain final
    : public abstract_domain<flat_boolean_numerical_domain<NumDom>> {
  typedef flat_boolean_numerical_domain<NumDom> bool_num_domain_t;
  typedef abstract_domain<bool_num_domain_t> abstract_domain_t;

public:
  typedef typename NumDom::number_t number_t;
  typedef typename NumDom::varname_t varname_t;
  typedef flat_boolean_domain<number_t, varname_t> bool_domain_t;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;

  typedef interval<number_t> interval_t;
  typedef bound<number_t> bound_t;

private:
  // This lattice is the dual of a discrete lattice where
  // elements are linear constraints.
  class lin_cst_set_domain : public writeable {

    typedef discrete_domain<linear_constraint_t> set_t;
    set_t m_set;

  public:
    typedef typename set_t::iterator iterator;

    lin_cst_set_domain(set_t s) : m_set(s) {}
    lin_cst_set_domain() : m_set(set_t::bottom()) /*top by default*/ {}
    lin_cst_set_domain(const lin_cst_set_domain &other) : m_set(other.m_set) {}

    static lin_cst_set_domain bottom() { return set_t::top(); }
    static lin_cst_set_domain top() { return set_t::bottom(); }

    bool is_top() { return m_set.is_bottom(); }
    bool is_bottom() { return m_set.is_top(); }

    bool operator<=(lin_cst_set_domain other) {
      if (other.is_top() || is_bottom())
        return true;
      else
        return other.m_set <= m_set;
    }

    bool operator==(lin_cst_set_domain other) {
      return (*this <= other && other <= *this);
    }

    void operator|=(lin_cst_set_domain other) { m_set = m_set & other.m_set; }

    lin_cst_set_domain operator|(lin_cst_set_domain other) {
      return lin_cst_set_domain(m_set & other.m_set);
    }

    lin_cst_set_domain operator&(lin_cst_set_domain other) {
      return lin_cst_set_domain(m_set | other.m_set);
    }

    lin_cst_set_domain operator||(lin_cst_set_domain other) {
      return this->operator|(other);
    }

    lin_cst_set_domain operator&&(lin_cst_set_domain other) {
      return this->operator&(other);
    }

    lin_cst_set_domain &operator+=(linear_constraint_t c) {
      m_set += c;
      return *this;
    }
    lin_cst_set_domain &operator-=(linear_constraint_t c) {
      m_set -= c;
      return *this;
    }

    std::size_t size() { return m_set.size(); }
    iterator begin() { return m_set.begin(); }
    iterator end() { return m_set.end(); }
    void write(crab::crab_os &o) {
      if (is_bottom())
        o << "_|_";
      else if (is_top())
        o << "top";
      else
        m_set.write(o);
    }
  };

  typedef domain_product2<number_t, varname_t, bool_domain_t, NumDom>
      domain_product2_t;

  // For performing reduction from the boolean domain to the
  // numerical one.
  // Map bool variables to sets of constraints such that if the
  // bool variable is true then the conjunction of the constraints
  // must be satisfiable.
  typedef separate_domain<variable_t, lin_cst_set_domain> var_lincons_map_t;
  /** we need to keep track which constraints still hold at the
      time the reduction from boolean variables to numerical ones
      is done. For instance,
      a := x > y;
      // unchanged = {x,y}
      if (*)
        x := 0;
        // unchanged = {y}
      else
        // unchanged = {x,y}

      // unchanged = {y}
      assume(a);
      // we cannot say here that x>y holds since x might have been modified
  **/

  // Simple wrapper for performing a must-forward dataflow
  // analysis that keeps track whether a variable has been
  // modified from any path since the variable was defined up to
  // the current location.
  class invariance_domain : public writeable {

  public:
    typedef typename linear_constraint_t::variable_t variable_t;
    typedef typename linear_constraint_t::variable_set_t variable_set_t;

  private:
    typedef discrete_domain<variable_t> set_t;
    set_t m_set;

  public:
    typedef typename set_t::iterator iterator;

    invariance_domain(set_t s) : m_set(s) {}
    invariance_domain(variable_set_t s) : m_set(set_t::bottom()) {
      for (auto v : s) {
        m_set += v;
      }
    }
    invariance_domain(variable_t v) : m_set(set_t::bottom()) { m_set += v; }

    invariance_domain() : m_set(set_t::bottom()) /*top by default*/ {}
    invariance_domain(const invariance_domain &other) : m_set(other.m_set) {}

    static invariance_domain bottom() { return set_t::top(); }
    static invariance_domain top() { return set_t::bottom(); }

    bool is_top() { return m_set.is_bottom(); }
    bool is_bottom() { return m_set.is_top(); }

    /*
             top  = everything might change
          {x,y,z} = nothing changes
              {x} = everything might change except x

                 top
                / | \
             {x} {y} {z}
             / \ / \ / \
           {x,y} {x,z} {y,z}
              \   |   /
               {x,y,z}
                  |
                 bot
     */
    bool operator<=(invariance_domain other) {
      if (other.is_top() || is_bottom())
        return true;
      else
        return other.m_set <= m_set;
    }

    void operator|=(invariance_domain other) { m_set = m_set & other.m_set; }

    invariance_domain operator|(invariance_domain other) {
      return invariance_domain(m_set & other.m_set);
    }

    invariance_domain operator&(invariance_domain other) {
      return invariance_domain(m_set | other.m_set);
    }

    invariance_domain operator||(invariance_domain other) {
      return this->operator|(other);
    }

    invariance_domain operator&&(invariance_domain other) {
      return this->operator&(other);
    }

    invariance_domain &operator+=(variable_t v) {
      m_set += v;
      return *this;
    }
    invariance_domain &operator-=(variable_t v) {
      m_set -= v;
      return *this;
    }

    bool operator[](const variable_t &v) {
      invariance_domain d(v);
      return (d <= *this);
    }
    std::size_t size() { return m_set.size(); }
    iterator begin() { return m_set.begin(); }
    iterator end() { return m_set.end(); }
    void write(crab::crab_os &o) {
      if (is_bottom())
        o << "_|_";
      else if (is_top())
        o << "top";
      else
        m_set.write(o);
    }
  };

  domain_product2_t _product;
  var_lincons_map_t _var_to_csts;
  invariance_domain _unchanged_vars;

  flat_boolean_numerical_domain(domain_product2_t &&product,
                                var_lincons_map_t &&var_to_csts,
                                invariance_domain &&unchanged_vars)
      : _product(std::move(product)), _var_to_csts(std::move(var_to_csts)),
        _unchanged_vars(std::move(unchanged_vars)) {}

public:
  void set_to_top() {
    bool_num_domain_t abs(domain_product2_t::top(), var_lincons_map_t::top(),
                          invariance_domain::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    bool_num_domain_t abs(domain_product2_t::bottom(),
                          var_lincons_map_t::bottom(),
                          invariance_domain::bottom());
    std::swap(*this, abs);
  }

  flat_boolean_numerical_domain()
      : _product(), _var_to_csts(), _unchanged_vars() {}

  flat_boolean_numerical_domain(const bool_num_domain_t &other)
      : _product(other._product), _var_to_csts(other._var_to_csts),
        _unchanged_vars(other._unchanged_vars) {}

  flat_boolean_numerical_domain(const bool_num_domain_t &&other)
    : _product(std::move(other._product)),
      _var_to_csts(std::move(other._var_to_csts)),
      _unchanged_vars(std::move(other._unchanged_vars)) {}
  
  bool_num_domain_t &operator=(const bool_num_domain_t &other) {
    if (this != &other) {
      _product = other._product;
      _var_to_csts = other._var_to_csts;
      _unchanged_vars = other._unchanged_vars;
    }
    return *this;
  }

  bool_num_domain_t &operator=(const bool_num_domain_t &&other) {
    if (this != &other) {
      _product = std::move(other._product);
      _var_to_csts = std::move(other._var_to_csts);
      _unchanged_vars = std::move(other._unchanged_vars);
    }
    return *this;
  }
  
  bool is_bottom() { return _product.is_bottom(); }

  bool is_top() { return _product.is_top(); }

  bool_domain_t &first() { return _product.first(); }

  NumDom &second() { return _product.second(); }

  bool operator<=(bool_num_domain_t other) {
    return _product <= other._product;
  }

  bool operator==(bool_num_domain_t other) {
    return _product == other._product;
  }

  void operator|=(bool_num_domain_t other) {
    _product |= other._product;
    _var_to_csts = _var_to_csts | other._var_to_csts;
    _unchanged_vars = _unchanged_vars | other._unchanged_vars;
  }

  bool_num_domain_t operator|(bool_num_domain_t other) {
    return bool_num_domain_t(_product | other._product,
                             _var_to_csts | other._var_to_csts,
                             _unchanged_vars | other._unchanged_vars);
  }

  bool_num_domain_t operator&(bool_num_domain_t other) {
    return bool_num_domain_t(_product & other._product,
                             _var_to_csts & other._var_to_csts,
                             _unchanged_vars & other._unchanged_vars);
  }

  bool_num_domain_t operator||(bool_num_domain_t other) {
    return bool_num_domain_t(_product || other._product,
                             _var_to_csts || other._var_to_csts,
                             _unchanged_vars || other._unchanged_vars);
  }

  bool_num_domain_t
  widening_thresholds(bool_num_domain_t other,
                      const iterators::thresholds<number_t> &ts) {
    return bool_num_domain_t(_product.widening_thresholds(other._product, ts),
                             _var_to_csts || other._var_to_csts,
                             _unchanged_vars || other._unchanged_vars);
  }

  bool_num_domain_t operator&&(bool_num_domain_t other) {
    return bool_num_domain_t(_product && other._product,
                             _var_to_csts && other._var_to_csts,
                             _unchanged_vars && other._unchanged_vars);
  }

  // numerical_domains_api

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    _product.apply(op, x, y, z);
    _unchanged_vars -= variable_t(x);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    _product.apply(op, x, y, k);
    _unchanged_vars -= variable_t(x);
  }

  void assign(variable_t x, linear_expression_t e) {
    _product.assign(x, e);
    _unchanged_vars -= variable_t(x);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       bool_num_domain_t invariant) override {
    _product.backward_assign(x, e, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      bool_num_domain_t invariant) override {
    _product.backward_apply(op, x, y, z, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      bool_num_domain_t invariant) override {
    _product.backward_apply(op, x, y, z, invariant._product);
    _unchanged_vars -= variable_t(x);
  }

  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraint");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraint");

    if (csts.is_true() || csts.is_false()) {
      return;
    }

    bool all_non_boolean = true;
    for (const linear_constraint_t &cst : csts) {
      auto const &variables = cst.variables();
      if (std::any_of(variables.begin(), variables.end(),
                      [](const variable_t &v) { return v.is_bool_type(); })) {
        all_non_boolean = false;
        break;
      }
    }

    if (all_non_boolean) {
      /// -- Common case: all constraints are non-boolean.
      _product.second() += csts;
    } else {
      /// -- We have both boolean and non-boolean constraints.

      // Normalization ensures that inequality pairs x <= k and x >= k
      // are replaced with a single equality x = k, where x always
      // appears as positive.
      linear_constraint_system_t norm_csts = csts.normalize();
      linear_constraint_system_t non_boolean_csts;
      for (const linear_constraint_t &cst : norm_csts) {
        auto const &variables = cst.variables();
        if (cst.is_equality() &&
            std::all_of(variables.begin(), variables.end(),
                        [](const variable_t &v) { return v.is_bool_type(); })) {
          // boolean component
          const linear_expression_t exp = cst.expression();
          if (exp.size() == 1) {
            number_t coef = (*(exp.begin())).first;
            variable_t var = (*(exp.begin())).second;
            number_t k = exp.constant();
            if (coef == number_t(1)) {
              if (k == number_t(0)) {
                _product.first().assume_bool(var, true /*is_negated*/);
              } else if (k == number_t(1)) {
                _product.first().assume_bool(var, false /*is_negated*/);
              }
            }
          }
        } else {
          // numerical component
          non_boolean_csts += cst;
        }
      }
      _product.second() += non_boolean_csts;
    }

    #if 0
    // update unchanged_vars
    for (auto v : csts.variables()) {
      _unchanged_vars -= v;
    }
    #endif 
  }

  void set(variable_t x, interval_t intv) {
    // domain_product2 does not define set method
    _product.second().set(x, intv); // only on the numerical domain
    _unchanged_vars -= variable_t(x);
  }

  interval_t operator[](variable_t v) {
    // domain_product2 does not define [] method
    boolean_value bv = _product.first().get_bool(v);
    interval_t isecond = _product.second()[v];

    if (bv.is_bottom() || isecond.is_bottom())
      return interval_t::bottom();

    if (bv.is_true())
      return interval_t(number_t(1)) & isecond;
    else if (bv.is_false())
      return interval_t(number_t(0)) & isecond;
    else
      return isecond;
  }

  void operator-=(variable_t v) {
    _product -= v;
    _var_to_csts -= v;
    _unchanged_vars -= variable_t(v);
  }

  // boolean_operators

  void assign_bool_cst(variable_t x, linear_constraint_t cst) {
    crab::CrabStats::count(getDomainName() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_cst");

    /// Reduction from the numerical domain to the flat boolean
    /// domain

    if (_product.is_bottom())
      return;

    _product.assign_bool_cst(x, cst);

    NumDom inv1(_product.second());
    inv1 += cst;
    if (inv1.is_bottom()) {
      // -- definitely false
      _product.first().set_bool(x, boolean_value::get_false());
    } else {
      NumDom inv2(_product.second());
      inv2 += cst.negate();
      if (inv2.is_bottom()) {
        // -- definitely true
        _product.first().set_bool(x, boolean_value::get_true());
      } else {
#if 0
	    // XXX: before we give up we convert into intervals and
	    // check again if the negated constraint is bottom.  This
	    // is useful because e.g., apron domains completely ignore
	    // disequations.
	    // 
	    // JN: I disable this code because AFIK all domains reason
	    // now about disequalities, included apron/elina domains.
	    interval_domain<number_t,varname_t> inv3;
	    inv3 += cst.negate();	    
	    for (auto c: _product.second().to_linear_constraint_system())
	    { inv3 += c;}
	    if (inv3.is_bottom()) {
	      // -- definitely true	  
	      _product.first().set_bool(x, boolean_value::get_true());
	    } else {
	      // -- inconclusive
	      _product.first().set_bool(x, boolean_value::top());
	    }
#else
        // -- inconclusive
        _product.first().set_bool(x, boolean_value::top());
#endif
      }
    }
    _var_to_csts.set(x, lin_cst_set_domain(cst));
    // We assume all variables in cst are unchanged unless the
    // opposite is proven
    for (auto v : cst.variables()) {
      _unchanged_vars += v;
    }

    CRAB_LOG("flat-boolean", auto bx = _product.first().get_bool(x);
             crab::outs() << "*** Reduction numerical --> boolean\n "
                          << "\t" << x << " := (" << cst << ")\n"
                          << "\t" << x << " := " << bx << "\n"
                          << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);
  }

  void assign_bool_var(variable_t x, variable_t y, bool is_not_y) {
    crab::CrabStats::count(getDomainName() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_var");

    if (is_bottom())
      return;

    _product.assign_bool_var(x, y, is_not_y);

    if (!is_not_y)
      _var_to_csts.set(x, _var_to_csts[y]);
    else {
      auto csts = _var_to_csts[y];
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        _var_to_csts.set(x, lin_cst_set_domain(cst.negate()));
        return;
      }
      // we do not negate multiple conjunctions because it would
      // become a disjunction so we give up
      if (csts.size() > 1)
        _var_to_csts -= x;
    }

    CRAB_LOG("flat-boolean",
             crab::outs() << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);
  }

  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply_binary_bool");

    if (is_bottom())
      return;

    _product.apply_binary_bool(op, x, y, z);

    // // --- for reduction from boolean to the numerical domain
    // if (op == OP_BAND) {
    //   _var_to_csts.set
    //     (x, _var_to_csts [y] & _var_to_csts [z]);
    //   return;
    // }

    // // we almost lose precision with or and xor except if one of
    // // the operands is false
    // if (op == OP_BOR || op == OP_BXOR) {
    //   if (_product.first().get_bool(y).is_false()) {
    //     _var_to_csts.set(x, _var_to_csts [z]);
    //     return;
    //   }
    //   if (_product.first().get_bool(z).is_false()) {
    //     _var_to_csts.set(x, _var_to_csts [y]);
    //     return;
    //   }
    // }

    /// otherwise we give up
    _var_to_csts -= x;
  }

  void assume_bool(variable_t x, bool is_negated) {
    crab::CrabStats::count(getDomainName() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".assume_bool");

    if (is_bottom())
      return;

    _product.assume_bool(x, is_negated);

    CRAB_LOG("flat-boolean",
             crab::outs() << "*** Reduction boolean --> numerical\n"
                          << "\tassume" << (is_negated ? "(not " : "(") << x
                          << ")\n"
                          << "\tINV=" << _product.second() << "\n"
                          << "\tunchanged vars=" << _unchanged_vars << "\n"
                          << "\tconstraints for reduction=" << _var_to_csts
                          << "\n";);

    if (_var_to_csts[x].is_top() || _var_to_csts[x].is_bottom())
      return;

    // Perform reduction from the flat boolean domain to the
    // numerical domain.
    if (!is_negated) {
      for (auto cst : _var_to_csts[x]) {
        // -- we only apply reduction if we know that all the
        // constraint variables have not been modified since they
        // were added into _var_to_csts.
        if (_unchanged_vars <= cst.variables()) {
          _product.second() += cst;
          CRAB_LOG("flat-boolean", crab::outs() << "\t"
                                                << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean", crab::outs()
                                       << "\t"
                                       << "Cannot applied " << cst << "\n";);
        }
      }
    } else {
      // we only perform reduction if there is only one constraint
      auto csts = _var_to_csts[x];
      if (csts.size() == 1) {
        auto cst = *(csts.begin());
        if (_unchanged_vars <= cst.variables()) {
          _product.second() += cst.negate();
          CRAB_LOG("flat-boolean", crab::outs() << "\t"
                                                << "Applied " << cst << "\n";);
        } else {
          CRAB_LOG("flat-boolean", crab::outs()
                                       << "\t"
                                       << "Cannot apply " << cst << "\n";);
        }
      }
    }

    CRAB_LOG("flat-boolean",
             crab::outs() << "After reduction=" << _product.second() << "\n";);
  }

  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                bool_num_domain_t inv) {
    /** TODO **/
    /*
       if lhs is true than assume(rhs)
       if lhs is false then assume(not rhs)
    */
    /** TODO: this can be done better **/
    _var_to_csts -= lhs;
  }

  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                bool_num_domain_t inv) {
    _product.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv._product);
    /** TODO: this can be done better **/
    _var_to_csts -= lhs;
  }

  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  bool_num_domain_t inv) {
    _product.backward_apply_binary_bool(op, x, y, z, inv._product);
    /** TODO: this can be done better **/
    _var_to_csts -= x;
  }

  // cast_operators_api

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    CRAB_LOG("flat-boolean", crab::outs() << src << ":" << src.get_bitwidth()
                                          << " " << op << " " << dst << ":"
                                          << dst.get_bitwidth() << " with "
                                          << *this << "\n");

    if (op == OP_TRUNC && (src.get_bitwidth() > 1 && dst.get_bitwidth() == 1)) {
      // -- int to bool:
      // assume that zero is false and non-zero is true
      interval_t i_src = _product.second()[src];
      interval_t zero = interval_t(number_t(0));
      if (i_src == zero) {
        _product.first().set_bool(dst, boolean_value::get_false());
      } else if (!(zero <= i_src)) {
        _product.first().set_bool(dst, boolean_value::get_true());
      } else {
        _product.first().set_bool(dst, boolean_value::top());
      }
    } else if ((op == OP_ZEXT || op == OP_SEXT) &&
               (src.get_bitwidth() == 1 && dst.get_bitwidth() > 1)) {
      // -- bool to int:
      // if OP_SEXT then true is -1 and false is zero
      // if OP_ZEXT then true is 1 and false is zero
      boolean_value b_src = _product.first().get_bool(src);
      if (b_src.is_true()) {
        _product.second().assign(dst,
                                 linear_expression_t(op == OP_SEXT ? -1 : 1));
      } else if (b_src.is_false()) {
        _product.second().assign(dst, linear_expression_t(0));
      } else {
        _product.second() -= dst;
        if (op == OP_SEXT) {
          _product.second() += linear_constraint_t(variable_t(dst) >= -1);
          _product.second() += linear_constraint_t(variable_t(dst) <= 0);
        } else {
          _product.second() += linear_constraint_t(variable_t(dst) >= 0);
          _product.second() += linear_constraint_t(variable_t(dst) <= 1);
        }
      }
      _unchanged_vars -= variable_t(dst);
    } else {
      _product.apply(op, dst, src);
      _unchanged_vars -= variable_t(dst);
    }

    CRAB_LOG("flat-boolean", crab::outs() << *this << "\n");
  }

  // bitwise_operators_api

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    _product.apply(op, x, y, z);
    _unchanged_vars -= variable_t(x);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    _product.apply(op, x, y, k);
    _unchanged_vars -= variable_t(x);
  }

  // array_operators_api

  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {
    _product.array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {
    _product.array_load(lhs, a, elem_size, i);
    if (a.get_type() == ARR_INT_TYPE || a.get_type() == ARR_REAL_TYPE)
      _unchanged_vars -= variable_t(lhs);
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    _product.array_store(a, elem_size, i, val, is_strong_update);
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t elem_size, linear_expression_t i,
                           linear_expression_t val,
                           bool is_strong_update) override {
    _product.array_store(a_new, a_old, elem_size, i, val, is_strong_update);
  }

  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t v) override {
    _product.array_store_range(a, elem_size, i, j, v);
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t v) override {
    _product.array_store_range(a_new, a_old, elem_size, i, j, v);
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    _product.array_assign(lhs, rhs);
  }

  // backward array operations

  virtual void backward_array_init(variable_t a, linear_expression_t elem_size,
                                   linear_expression_t lb_idx,
                                   linear_expression_t ub_idx,
                                   linear_expression_t val,
                                   bool_num_domain_t invariant) override {
    _product.backward_array_init(a, elem_size, lb_idx, ub_idx, val,
                                 invariant._product);
  }

  virtual void backward_array_load(variable_t lhs, variable_t a,
                                   linear_expression_t elem_size,
                                   linear_expression_t i,
                                   bool_num_domain_t invariant) override {
    _product.backward_array_load(lhs, a, elem_size, i, invariant._product);
    if (a.get_type() == ARR_INT_TYPE || a.get_type() == ARR_REAL_TYPE) {
      _unchanged_vars -= variable_t(lhs);
    }
  }

  virtual void backward_array_store(variable_t a, linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    bool_num_domain_t invariant) override {
    _product.backward_array_store(a, elem_size, i, val, is_strong_update,
                                  invariant._product);
  }

  virtual void backward_array_store(variable_t a_new, variable_t a_old,
                                    linear_expression_t elem_size,
                                    linear_expression_t i,
                                    linear_expression_t val,
                                    bool is_strong_update,
                                    bool_num_domain_t invariant) override {
    _product.backward_array_store(a_new, a_old, elem_size, i, val,
                                  is_strong_update, invariant._product);
  }

  virtual void
  backward_array_store_range(variable_t a, linear_expression_t elem_size,
                             linear_expression_t i, linear_expression_t j,
                             linear_expression_t v,
                             bool_num_domain_t invariant) override {
    _product.backward_array_store_range(a, elem_size, i, j, v,
                                        invariant._product);
  }

  virtual void backward_array_store_range(
      variable_t a_new, variable_t a_old, linear_expression_t elem_size,
      linear_expression_t i, linear_expression_t j, linear_expression_t v,
      bool_num_domain_t invariant) override {
    _product.backward_array_store_range(a_new, a_old, elem_size, i, j, v,
                                        invariant._product);
  }

  virtual void backward_array_assign(variable_t lhs, variable_t rhs,
                                     bool_num_domain_t invariant) override {
    _product.backward_array_assign(lhs, rhs, invariant._product);
  }

  // pointer_operators_api
  virtual void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) override {
    _product.pointer_load(lhs, rhs, elem_size);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) override {
    _product.pointer_store(lhs, rhs, elem_size);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) override {
    _product.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    _product.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    _product.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    _product.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(pointer_constraint_t cst) override {
    _product.pointer_assume(cst);
  }

  virtual void pointer_assert(pointer_constraint_t cst) override {
    _product.pointer_assert(cst);
  }

  void write(crab_os &o) { _product.write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    linear_constraint_system_t res;
    res += _product.first().to_linear_constraint_system();
    res += _product.second().to_linear_constraint_system();
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    disjunctive_linear_constraint_system_t res;
    res += _product.first().to_disjunctive_linear_constraint_system();
    res += _product.second().to_disjunctive_linear_constraint_system();
    return res;
  }

  static std::string getDomainName() {
    return domain_product2_t::getDomainName();
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    _product.rename(from, to);
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    _product.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  bool_num_domain_t invariant) override {
    _product.backward_intrinsic(name, inputs, outputs, invariant._product);    
  }
  /* end intrinsics operations */
  
  void normalize() { _product.normalize(); }

  void minimize() { _product.minimize(); }

  void forget(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    _product.forget(variables);
    for (variable_t v : variables) {
      _var_to_csts -= v;
      _unchanged_vars -= v;
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    if (variables.empty()) {
      set_to_top();
      return;
    }

    _product.project(variables);

    var_lincons_map_t new_var_to_csts;
    invariance_domain new_unchanged_vars;
    for (variable_t v : variables) {
      new_var_to_csts.set(v, _var_to_csts[v]);
      if (_unchanged_vars[v]) {
        new_unchanged_vars += v;
      }
    }
    std::swap(_var_to_csts, new_var_to_csts);
    std::swap(_unchanged_vars, new_unchanged_vars);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    _product.expand(x, new_x);
    _var_to_csts.set(new_x, _var_to_csts[x]);
    if (_unchanged_vars[variable_t(x)]) {
      _unchanged_vars += variable_t(new_x);
    }
  }
}; // class flat_boolean_numerical_domain

template <typename Num>
struct abstract_domain_traits<flat_boolean_numerical_domain<Num>> {
  typedef typename Num::number_t number_t;
  typedef typename Num::varname_t varname_t;
};

template <typename NumDom>
class checker_domain_traits<flat_boolean_numerical_domain<NumDom>> {
public:
  typedef flat_boolean_numerical_domain<NumDom> this_type;
  typedef typename this_type::linear_constraint_t linear_constraint_t;
  typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    NumDom &lhs_dom = lhs.second();
    return checker_domain_traits<NumDom>::entail(lhs_dom, rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    NumDom &rhs_dom = rhs.second();
    return checker_domain_traits<NumDom>::entail(lhs, rhs_dom);
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    NumDom &lhs_dom = lhs.second();
    return checker_domain_traits<NumDom>::entail(lhs_dom, rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    NumDom &dom = inv.second();
    return checker_domain_traits<NumDom>::intersect(dom, cst);
  }
};

} // end namespace domains
} // end namespace crab
