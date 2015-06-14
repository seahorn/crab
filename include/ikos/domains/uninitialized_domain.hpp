/*******************************************************************************
 * Simple lattice for uninitialized variables.
 ******************************************************************************/

#ifndef IKOS_UNINITIALIZED_VARIABLES_HPP
#define IKOS_UNINITIALIZED_VARIABLES_HPP

#include <ikos/common/types.hpp>
#include <ikos/domains/separate_domains.hpp>
#include <ikos/domains/domain_products.hpp>

namespace ikos {

// A simple lattice for uninitialization
class uninitializedValue : public ikos::writeable {
  typedef enum {
    Bottom = 0x0, /*unused*/
    Initialized = 0x1,
    Uninitialized = 0x2,
    MayUninitialized = 0x3
  } kind_t;

  kind_t _value;

  uninitializedValue(kind_t v) : _value(v){};

public:
  uninitializedValue() : _value(MayUninitialized) {}

public:
  static uninitializedValue bottom() { return uninitializedValue(Bottom); }

  static uninitializedValue top() {
    return uninitializedValue(MayUninitialized);
  }

  static uninitializedValue initialized() {
    return uninitializedValue(Initialized);
  }

  static uninitializedValue uninitialized() {
    return uninitializedValue(Uninitialized);
  }

public:
  uninitializedValue(const uninitializedValue& other)
      : writeable(), _value(other._value) {}

  uninitializedValue& operator=(uninitializedValue other) {
    _value = other._value;
    return *this;
  }

  bool is_bottom() { return (_value == Bottom); }

  bool is_top() { return (_value == MayUninitialized); }

  bool is_initialized() const { return (_value == Initialized); }

  bool is_uninitialized() const { return (_value == Uninitialized); }

  bool operator<=(uninitializedValue other) {
    if (_value == Bottom)
      return true;
    else if (_value == MayUninitialized)
      return (other._value == MayUninitialized);
    else if (_value == Initialized)
      return ((other._value == _value) ||
              (other._value == MayUninitialized));
    else if (_value == Uninitialized)
      return (_value <= other._value);
    else {
      assert(false && "unreachable");
      return false;
    }
  }

  bool operator==(uninitializedValue other) {
    return (_value == other._value);
  }

  uninitializedValue operator|(uninitializedValue other) {
    return uninitializedValue(static_cast< kind_t >(
        static_cast< int >(_value) | static_cast< int >(other._value)));
  }

  uninitializedValue operator||(uninitializedValue other) {
    return this->operator|(other);
  }

  uninitializedValue operator&(uninitializedValue other) {
    return uninitializedValue(static_cast< kind_t >(
        static_cast< int >(_value) & static_cast< int >(other._value)));
  }

  uninitializedValue operator&&(uninitializedValue other) {
    return this->operator&(other);
  }

  ostream& write(ostream& o) {
    switch (_value) {
      case Bottom: {
        o << "_|_";
        break;
      }
      case MayUninitialized: {
        o << "T";
        break;
      }
      case Initialized: {
        o << "I";
        break;
      }
      case Uninitialized: {
        o << "U";
        break;
      }
    }
    return o;
  }
}; // end class uninitializedValue

// An abstract domain for reasoning about uninitialized variables.
template < typename VariableName >
class uninitialized_domain : public ikos::writeable {
  template < typename Any1, typename Any2, typename Any3 >
  friend class uninitialized_array_domain;

public:
  typedef uninitialized_domain< VariableName > uninitialized_domain_t;

private:
  typedef separate_domain< VariableName, uninitializedValue > separate_domain_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

  uninitialized_domain(separate_domain_t env) : _env(env) {}

  struct mkVal
      : public std::unary_function< VariableName, uninitializedValue > {
    separate_domain_t _env;
    mkVal(separate_domain_t env) : _env(env) {}
    uninitializedValue operator()(VariableName v) { return _env[v]; }
  };

public:
  static uninitialized_domain_t top() {
    return uninitialized_domain(separate_domain_t::top());
  }

  static uninitialized_domain_t bottom() {
    return uninitialized_domain(separate_domain_t::bottom());
  }

public:
  uninitialized_domain() : _env(separate_domain_t::top()) {}

  uninitialized_domain(const uninitialized_domain_t& e)
      : writeable(), _env(e._env) {}

  uninitialized_domain_t& operator=(uninitialized_domain_t e) {
    _env = e._env;
    return *this;
  }

  iterator begin() { return _env.begin(); }

  iterator end() { return _env.end(); }

  bool is_bottom() { return _env.is_bottom(); }

  bool is_top() { return _env.is_top(); }

  bool operator<=(uninitialized_domain_t e) { return (_env <= e._env); }

  bool operator==(uninitialized_domain_t e) {
    return ((_env <= e._env) && (e._env <= _env));
  }

  uninitialized_domain_t operator|(uninitialized_domain_t e) {
    return (_env | e._env);
  }

  uninitialized_domain_t operator&(uninitialized_domain_t e) {
    return (_env & e._env);
  }

  uninitialized_domain_t operator||(uninitialized_domain_t e) {
    return (_env || e._env);
  }

  uninitialized_domain_t operator&&(uninitialized_domain_t e) {
    return (_env && e._env);
  }

  void set(VariableName v, uninitializedValue e) {
    if (is_bottom()) return;

    _env.set(v, e);
  }

  void assign(VariableName x, VariableName y) {
    if (is_bottom()) return;

    _env.set(x, _env[y]);
  }

  template<typename Iterator>
  void assign(VariableName x, Iterator It, Iterator End) {
    if (is_bottom()) return;
    
    std::vector< uninitializedValue > out;
    mkVal f(_env);
    std::transform(It, End, back_inserter(out), f);
    assign(x, out);
  }

  // if all elements of ys are initialized so does x. If some ys is
  // uninitialized so does x. Otherwise, x is top.
  void assign(VariableName x, std::vector< uninitializedValue > ys) {
    if (is_bottom()) return;
      
    bool all_init = true;
    bool some_uninit = false;
    for (auto abs_val: ys) {
      all_init &= abs_val.is_initialized();
      some_uninit |= abs_val.is_uninitialized();
    }
    if (all_init)
      _env.set(x, uninitializedValue::initialized());
    else if (some_uninit)
      _env.set(x, uninitializedValue::uninitialized());
    else
      _env.set(x, uninitializedValue::top());
  }

  void apply(operation_t /*op*/, VariableName x, VariableName y, 
             uninitializedValue z) {
    if (is_bottom()) return;

    apply(x, _env[y], z);
  }

  void apply(operation_t /*op*/, VariableName x, VariableName y, VariableName z) {
    if (is_bottom()) return;

    apply(x, _env[y], _env[z]);
  }

  void apply(VariableName x, VariableName y, VariableName z) {
    if (is_bottom()) return;

    apply(x, _env[y], _env[z]);
  }

  void apply(VariableName x, uninitializedValue y, uninitializedValue z) {
    if (is_bottom()) return;
      
    if (y.is_uninitialized() || z.is_uninitialized())
      _env.set(x, uninitializedValue::uninitialized());
    else
      _env.set(x, y | z);
  }

  void operator-=(VariableName v) { _env -= v; }

  uninitializedValue operator[](VariableName v) {
    return _env[v]; 
  }

  ostream& write(ostream& o) { 
    _env.write(o); 
    return o;
  }

}; // class uninitialized_domain

} // end ikos namespace

#endif 
