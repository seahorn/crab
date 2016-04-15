#ifndef CRAB_COMMON_HPP
#define CRAB_COMMON_HPP

#include <stdint.h>
#include <string>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <errno.h>
#include <boost/optional.hpp>

/* Basic type definitions */


template<typename... ArgTypes>
inline void ___print___(ArgTypes... args)
{
  // trick to expand variadic argument pack without recursion
  using expand_variadic_pack = int[];
  // first zero is to prevent empty braced-init-list
  // void() is to prevent overloaded operator, messing things up
  // trick is to use the side effect of list-initializer to call a function
  // on every argument.
  // (void) is to suppress "statement has no effect" warnings
    (void)expand_variadic_pack{0, ((std::cerr << args), void(), 0)... };
}

#define CRAB_ERROR(...)              \
  do {                               \
    ___print___(__VA_ARGS__);        \
    std::cerr << "\n";               \
    std::exit (EXIT_FAILURE);        \
  } while (0)

#define CRAB_WARN(...)               \
  do {                               \
    std::cerr << "WARNING:";         \
    ___print___(__VA_ARGS__);        \
    std::cerr << "\n";               \
  } while (0)


namespace crab {

   typedef enum { 
     BINOP_ADD, BINOP_SUB, BINOP_MUL, 
     BINOP_SDIV, BINOP_UDIV, BINOP_SREM, BINOP_UREM,
     BINOP_AND, BINOP_OR, BINOP_XOR, BINOP_SHL, BINOP_LSHR, BINOP_ASHR
   } binary_operation_t;

   inline std::ostream& operator<<(std::ostream&o, binary_operation_t op) {
     switch (op) {
       case BINOP_ADD: o << "+"; break;
       case BINOP_SUB: o << "-"; break;
       case BINOP_MUL: o << "*"; break;
       case BINOP_SDIV: o << "/"; break;
       case BINOP_UDIV: o << "/_u"; break;
       case BINOP_SREM: o << "%"; break;
       case BINOP_UREM: o << "%_u"; break;
       case BINOP_AND: o << "&"; break;
       case BINOP_OR: o << "|"; break;
       case BINOP_XOR: o << "^"; break;
       case BINOP_SHL: o << "<<"; break;
       case BINOP_LSHR: o << ">>_l"; break;
       case BINOP_ASHR: o << ">>_r"; break;
       default: CRAB_ERROR("unreachable");
      }
     return o;
   }

  template<typename T>
  inline boost::optional<T> convOp (binary_operation_t op); 


  // toy language for pointer constraints
  typedef enum  { PTR_EQUALITY, PTR_DISEQUALITY } ptr_cst_kind_t;  

  template<typename VariableName>
  class pointer_constraint {
   public:

    typedef pointer_constraint <VariableName> ptr_cst_t;
    typedef boost::optional<VariableName> opt_var_t;

   private:

    opt_var_t _lhs; 
    opt_var_t _rhs;
    ptr_cst_kind_t _kind;

    /* 
       We can only express constraints of the form p {==,!=} q where
       p,q can be either null or a variable.

       !lhs, !rhs, kind == PTR_EQUALITY    -> null == null -> true
       !lhs, !rhs, kind == PTR_DISEQUALITY -> null != null -> false
       !lhs, rhs, kind == PTR_EQUALITY     -> *rhs == null
       lhs, !rhs, kind == PTR_EQUALITY     -> *lhs == null
       !lhs, rhs, kind == PTR_DISEQUALITY  -> *rhs != null
       lhs, !rhs, kind == PTR_DISEQUALITY  -> *lhs != nul
       lhs, rhs, kind == PTR_EQUALITY      -> *lhs == *rhs
       lhs, rhs, kind == PTR_DISEQUALITY   -> *lhs != *rhs
    */

    pointer_constraint (opt_var_t lhs, opt_var_t rhs, ptr_cst_kind_t kind):
        _lhs (lhs), _rhs (rhs), _kind (kind) {
      if (!_lhs && _rhs){ // normalize
        std::swap (_lhs,_rhs);
      }
    }

   public:

    pointer_constraint () : _kind (PTR_EQUALITY) { }
    
    // return true iff null != null
    bool is_contradiction () const {
      return (!_lhs && !_rhs && _kind == PTR_DISEQUALITY);
    }

    // return true iff null == null
    bool is_tautology () const {
      return (!_lhs && !_rhs && _kind == PTR_EQUALITY);
    }

    bool is_equality () const {
      return (_kind == PTR_EQUALITY);
    }

    bool is_disequality () const {
      return (_kind == PTR_DISEQUALITY);
    }

    // return true iff  p == null or p != null
    bool is_unary () const {
      return (_lhs && !_rhs);
    }

    // return true iff p == q or p != q
    bool is_binary () const {
      return (_lhs && _rhs);
    }

    VariableName lhs () const {
      if (!_lhs) CRAB_ERROR ("pointer constraint lhs is null");
      return *_lhs;
    }

    VariableName rhs () const {
      if (!_rhs) CRAB_ERROR ("pointer constraint rhs is null");
      return *_rhs;
    }

    static ptr_cst_t mk_true () {
      return ptr_cst_t ();
    }

    static ptr_cst_t mk_false () {
      return ptr_cst_t (opt_var_t (), opt_var_t (), PTR_DISEQUALITY);
    }

    static ptr_cst_t mk_eq_null (VariableName v)  {
      return ptr_cst_t (opt_var_t (v), opt_var_t (), PTR_EQUALITY);
    }

    static ptr_cst_t mk_diseq_null (VariableName v) {
      return ptr_cst_t (opt_var_t (v), opt_var_t (), PTR_DISEQUALITY);
    }

    static ptr_cst_t mk_eq (VariableName v1, VariableName v2)  {
      return ptr_cst_t (opt_var_t (v1), opt_var_t (v2), PTR_EQUALITY);
    }

    static ptr_cst_t mk_diseq (VariableName v1, VariableName v2)  {
      return ptr_cst_t (opt_var_t (v1), opt_var_t (v2), PTR_DISEQUALITY);
    }

    void write (std::ostream& o) const {
      if (is_contradiction () ) {
        o << "false";
      } else if (is_tautology ()) {
        o << "true";
      } else {
        assert (_lhs);

        o << lhs ();
        
        if (_kind == PTR_EQUALITY) 
          o << " == ";
        else 
          o << " != ";
        
        if (!_rhs)
          o << "NULL";
        else 
          o <<  rhs ();
      }
    }
    
    friend std::ostream& operator<<(std::ostream& o, const ptr_cst_kind_t& k) {
      if (k == PTR_EQUALITY ) 
        o << " == ";
      else 
        o << " != ";
      return o;
    }
    
    friend std::ostream& operator<<(std::ostream& o, const ptr_cst_t& cst) {
      cst.write (o);
      return o;
    }

  };
} // end namespace
 
namespace ikos {
  // Numerical type for indexed objects
  typedef uint64_t index_t;

  // Interface for writeable objects
  class writeable {
   public:
    virtual void write(std::ostream& o) = 0;
    virtual ~writeable() {}
  }; // class writeable

  inline std::ostream& operator<<(std::ostream& o, writeable& x) {
    x.write(o);
    return o;
  }

  // Container for typed variables
  template< typename Type, typename VariableName >
  class variable: public writeable {
  
   public:
    typedef variable< Type, VariableName > variable_t;
    typedef typename VariableName::index_t index_t;
    
   private:
    VariableName _n;
    boost::optional<Type> _ty;
    
   public:
    
    variable(const VariableName& n): writeable (), _n(n) { }

    variable(VariableName&& o): _n(std::move(o)) { }
    
    variable(const VariableName& n, const Type& ty): 
        writeable (), _n(n), _ty (ty) { }
    
    variable(const variable_t& v): 
        writeable(), _n(v._n), _ty (v._ty) { }
    
    variable_t& operator=(const variable_t &o) {
      if (this != &o) {
        this->_n = o._n;
        this->_ty = o._ty;
      }
      return *this;
    }
    
    const VariableName& name() const { return _n; }
    
    boost::optional<Type> type() const { return _ty; }
    
    index_t index() const { return _n.index(); }
    
    // bool operator==(const variable_t& o) const {
    //   return _n.index () == o._n.index ();
    // }
    
    bool operator<(const variable_t& o) const {
      // ignore type for comparing variables
      return _n.index () < o._n.index ();
    }
    
    void write(std::ostream& o) { 
      if (_ty)
        o << _n << ":" << *_ty; 
      else
        o << _n;
    }
    
    friend index_t hash_value (const variable_t& v) {
      // ignore type for computing hash value
      return v.index ();
    }
  
    friend std::ostream& operator<<(std::ostream& o, variable_t& v) {
      v.write(o);
      return o;
    }
    
  }; // class variable

} // end namespace 

#endif // CRAB_COMMON_HPP
