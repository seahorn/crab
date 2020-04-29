#pragma once

/* A flat lattice for nullity */

#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>

namespace crab {

namespace domains {

class nullity_value {
  /*
            Top
            / \
           /   \
        Null  NonNull
           \   /
            \ /
           Bottom
  */
  typedef enum {
    Bottom = 0x0, /*unused*/
    Null = 0x1,
    NonNull = 0x2,
    Top = 0x3
  } kind_t;

  kind_t _value;

  nullity_value(kind_t v) : _value(v){};

public:
  nullity_value() : _value(Top) {}

  static nullity_value bottom() { return nullity_value(Bottom); }

  static nullity_value top() { return nullity_value(Top); }

  static nullity_value non_null() { return nullity_value(NonNull); }

  static nullity_value null() { return nullity_value(Null); }

  nullity_value(const nullity_value &other) : _value(other._value) {}

  nullity_value &operator=(const nullity_value &other) {
    if (this != &other) {
      _value = other._value;
    }
    return *this;
  }

  bool is_bottom() { return (_value == Bottom); }

  bool is_top() { return (_value == Top); }

  bool is_non_null() const { return (_value == NonNull); }

  bool is_null() const { return (_value == Null); }

  bool operator<=(nullity_value other) {

    if (_value == Bottom || other._value == Top)
      return true;
    else if (_value == Top)
      return (other._value == Top);
    else if (_value == NonNull)
      return (other._value >= NonNull);
    else if (_value == Null)
      return ((other._value == Null) || (other._value == Top));

    // this should be unreachable
    return false;
  }

  bool operator==(nullity_value other) { return (_value == other._value); }

  nullity_value operator|(nullity_value other) {
    return nullity_value(static_cast<kind_t>(static_cast<int>(this->_value) |
                                             static_cast<int>(other._value)));
  }

  // the lattice satisfy ACC so join is the widening
  nullity_value operator||(nullity_value other) {
    return this->operator|(other);
  }

  nullity_value operator&(nullity_value other) {
    return nullity_value(static_cast<kind_t>(static_cast<int>(this->_value) &
                                             static_cast<int>(other._value)));
  }

  // the lattice satisfy DCC so meet is the narrowing
  nullity_value operator&&(nullity_value other) {
    return this->operator&(other);
  }

  void write(crab_os &o) const {
    switch (_value) {
    case Bottom:
      o << "_|_";
      break;
    case Top:
      o << "T";
      break;
    case NonNull:
      o << "NN";
      break;
    default: /*Null*/
      o << "N";
    }
  }

}; // end class nullity_value

inline crab_os &operator<<(crab_os &o, const nullity_value &v) {
  v.write(o);
  return o;
}

// Abstract domain for nullity
template <typename Number, typename VariableName>
class nullity_domain final
    : public abstract_domain<nullity_domain<Number, VariableName>> {

  typedef nullity_domain<Number, VariableName> nullity_domain_t;
  typedef abstract_domain<nullity_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef interval<number_t> interval_t;

private:
  typedef separate_domain<variable_t, nullity_value> separate_domain_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  nullity_domain(separate_domain_t env) : _env(env) {}

public:
  void set_to_top() {
    nullity_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    nullity_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  nullity_domain() : _env(separate_domain_t::top()) {}

  nullity_domain(const nullity_domain_t &e) : _env(e._env) {}

  nullity_domain_t &operator=(const nullity_domain_t &o) {
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

  bool operator<=(nullity_domain_t o) { return (_env <= o._env); }

  nullity_domain_t operator|(nullity_domain_t o) {
    nullity_domain_t res(_env | o._env);
    CRAB_LOG("nullity", crab::outs() << "After join " << *this << " and " << o
                                     << "=" << res << "\n";);
    return res;
  }

  void operator|=(nullity_domain_t o) {
    CRAB_LOG("nullity", crab::outs()
                            << "After join " << *this << " and " << o << "=");
    _env = _env | o._env;
    CRAB_LOG("nullity", crab::outs() << *this << "\n");
  }

  nullity_domain_t operator&(nullity_domain_t o) {
    nullity_domain_t res(_env & o._env);
    CRAB_LOG("nullity", crab::outs() << "After meet " << *this << " and " << o
                                     << "=" << res << "\n");
    return res;
  }

  nullity_domain_t operator||(nullity_domain_t o) {
    nullity_domain_t res(_env || o._env);

    CRAB_LOG("nullity", crab::outs() << "After widening " << *this << " and "
                                     << o << "=" << res << "\n");
    return res;
  }

  nullity_domain_t
  widening_thresholds(nullity_domain_t o,
                      const iterators::thresholds<number_t> &) {
    nullity_domain_t res(_env || o._env);

    CRAB_LOG("nullity", crab::outs() << "After widening " << *this << " and "
                                     << o << "=" << res << "\n");
    return res;
  }

  nullity_domain_t operator&&(nullity_domain_t o) { return (_env && o._env); }

  void set_nullity(variable_t v, nullity_value n) {
    if (!is_bottom())
      _env.set(v, n);
  }

  void set_nullity(variable_t x, variable_t y) {
    if (!is_bottom())
      _env.set(x, _env[y]);

    CRAB_LOG("nullity", crab::outs() << "After " << x << ":=" << y << "="
                                     << *this << "\n");
  }

  nullity_value get_nullity(variable_t v) { return _env[v]; }

  void operator-=(variable_t v) {
    if (!is_bottom())
      _env -= v;
  }

  void equality(variable_t p, variable_t q) {
    if (is_bottom())
      return;

    // if (p == q) ...
    nullity_value p_meet_q = _env[p] & _env[q];
    _env.set(p, p_meet_q);
    _env.set(q, p_meet_q);

    CRAB_LOG("nullity", crab::outs() << "After " << p << "==" << q << "="
                                     << *this << "\n");
  }

  void equality(variable_t p, nullity_value v) {
    if (!is_bottom())
      _env.set(p, _env[p] & v);

    CRAB_LOG("nullity", crab::outs() << "After " << p << "==" << v << "="
                                     << *this << "\n");
  }

  void disequality(variable_t p, variable_t q) {
    if (is_bottom())
      return;

    // if (p != q) ...
    if (_env[p].is_null() && _env[q].is_null()) {
      set_to_bottom();
    } else if (_env[p].is_top() && _env[q].is_null()) {
      _env.set(p, nullity_value::non_null()); // refine p
    } else if (_env[q].is_top() && _env[p].is_null()) {
      _env.set(q, nullity_value::non_null()); // refine q
    }

    CRAB_LOG("nullity", crab::outs() << "After " << p << "!=" << q << "="
                                     << *this << "\n");
  }

  void disequality(variable_t p, nullity_value v) {
    if (is_bottom())
      return;

    if (_env[p].is_null() && v.is_null()) {
      set_to_bottom();
    } else if (_env[p].is_top() && v.is_null()) { // refine p
      _env.set(p, nullity_value::non_null());
    }

    CRAB_LOG("nullity", crab::outs() << "After " << p << "!=" << v << "="
                                     << *this << "\n");
  }

  /*
     Begin unimplemented operations

     The nullity domain only implements pointer operations. The
     implementation of the rest of operations (e.g., numerical,
     boolean, array, etc) is empty because they should not be
     called.
  */
  // arithmetic operations
  // XXX: needed for making a reduced product with a numerical domain
  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {}
  void apply(operation_t op, variable_t x, variable_t y, number_t k) {}
  void assign(variable_t x, linear_expression_t e) {}
  void backward_assign(variable_t x, linear_expression_t e,
                       nullity_domain_t invariant) {}
  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      nullity_domain_t invariant) {}
  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      nullity_domain_t invariant) {}
  void operator+=(linear_constraint_system_t csts) {}
  void operator+=(linear_constraint_t cst) {}
  // not part of the numerical_domains api but it should be
  void set(variable_t x, interval_t intv) {}
  interval_t operator[](variable_t x) { return interval_t::top(); }

  // int cast operations and bitwise operations
  // XXX: needed for making a reduced product with a numerical domain
  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {}
  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
  }
  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t z) {}

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                nullity_domain_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                nullity_domain_t invariant) {}
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  nullity_domain_t invariant) {}
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
                           nullity_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           nullity_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, nullity_domain_t invariant) {
  }
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, nullity_domain_t invariant) {
  }
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  nullity_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  nullity_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             nullity_domain_t invariant) {}
  /* End unimplemented operations */

  // pointer operations
  void pointer_load(variable_t /*lhs*/, variable_t rhs, linear_expression_t elem_size) {
    // XXX: assume after the load the rhs must be non-null otherwise
    // the program failed.
    equality(rhs, nullity_value::non_null());
  }

  void pointer_store(variable_t lhs, variable_t /*rhs*/, linear_expression_t elem_size) {
    // XXX: assume after the store the lhs must be non-null
    // otherwise the program failed.
    equality(lhs, nullity_value::non_null());
  }

  void pointer_assign(variable_t lhs, variable_t rhs,
                      linear_expression_t /*offset*/) {
    set_nullity(lhs, rhs);
    CRAB_LOG("nullity", crab::outs() << "After " << lhs << ":=" << rhs << "="
                                     << *this << "\n");
  }

  void pointer_mk_obj(variable_t lhs, ikos::index_t /*address*/) {
    set_nullity(lhs, nullity_value::non_null());
    CRAB_LOG("nullity", crab::outs() << "After " << lhs << ":= mk_object()"
                                     << "=" << *this << "\n");
  }

  void pointer_function(variable_t lhs, varname_t /*func*/) {
    set_nullity(lhs, nullity_value::non_null());
  }

  void pointer_mk_null(variable_t lhs) {
    set_nullity(lhs, nullity_value::null());
    CRAB_LOG("nullity", crab::outs() << "After " << lhs << ":= NULL"
                                     << "=" << *this << "\n");
  }

  void pointer_assume(pointer_constraint_t cst) {
    if (cst.is_tautology())
      return;

    if (cst.is_contradiction()) {
      set_to_bottom();
      return;
    }

    if (cst.is_unary()) {
      if (cst.is_equality())
        equality(cst.lhs(), nullity_value::null());
      else // cst.is_disequality();
        disequality(cst.lhs(), nullity_value::null());
    } else {
      assert(cst.is_binary());
      if (cst.is_equality())
        equality(cst.lhs(), cst.rhs());
      else // cst.is_disequality();
        disequality(cst.lhs(), cst.rhs());
    }
  }

  void pointer_assert(pointer_constraint_t cst) {
    CRAB_WARN("nullity pointer_assert not implemented");
  }

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    nullity_domain_t res;
    for (variable_t v : variables) {
      res.set_nullity(v, get_nullity(v));
    }
    std::swap(*this, res);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set_nullity(new_x, get_nullity(x));
  }

  void normalize() {}

  void minimize() {}

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
			  nullity_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  linear_constraint_system_t to_linear_constraint_system() {
    if (is_bottom())
      return linear_constraint_t::get_false();
    else
      return linear_constraint_t::get_true();
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

  static std::string getDomainName() { return "Nullity"; }

  void write(crab_os &o) { _env.write(o); }

}; // class nullity_domain

template <typename Number, typename VariableName>
struct abstract_domain_traits<nullity_domain<Number, VariableName>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

} // end namespace domains
} // end namespace crab
