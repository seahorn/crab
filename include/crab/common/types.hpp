/*******************************************************************************
 *
 * Basic type definitions.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011-2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#ifndef IKOS_COMMON_HPP
#define IKOS_COMMON_HPP

#include <stdint.h>
#include <string>
#include <iostream>
#include <stdarg.h>
#include <errno.h>
#include <boost/noncopyable.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/slist.hpp>
#include <boost/iterator/iterator_facade.hpp>

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

// Default definition but it can be redefined by dbg.hpp
#define CRAB_DEBUG(...)

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

} // end namespace crab

namespace ikos 
{
  
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

template < typename Element >
class collection : public writeable {
public:
  typedef collection< Element > collection_t;

private:
  typedef boost::container::slist< Element > slist_t;
  typedef boost::shared_ptr< slist_t > slist_ptr;

private:
  slist_ptr _slist;

public:
  class iterator
      : public boost::iterator_facade< iterator,
                                       const Element,
                                       boost::forward_traversal_tag > {
    friend class boost::iterator_core_access;

  private:
    typename slist_t::const_iterator _it;
    slist_ptr _l;

  public:
    iterator(slist_ptr l, bool b) : _it(b ? l->begin() : l->end()), _l(l) {}

  private:
    void increment() { ++(this->_it); }

    bool equal(const iterator& other) const {
      return this->_l == other._l && this->_it == other._it;
    }

    const Element& dereference() const {
      if (this->_it != this->_l->end()) {
        return *(this->_it);
      } else {
        CRAB_ERROR("Collection: trying to dereference an empty iterator");
      }
    }

  }; // class iterator

public:
  collection() : _slist(slist_ptr(new slist_t)) {}

  collection(const collection_t& c) : writeable(), _slist(c._slist) {}

  collection_t& operator=(const collection_t& c) {
    this->_slist = c._slist;
    return *this;
  }

  collection_t& operator+=(const Element& e) {
    this->_slist->push_front(e);
    return *this;
  }

  collection_t& operator+=(const collection_t& c) {
    for (iterator it = c.begin(); it != c.end(); ++it) {
      this->_slist->push_front(*it);
    }
    return *this;
  }

  collection_t operator+(const collection_t& c) const {
    collection_t r;
    r.operator+=(c);
    r.operator+=(*this);
    return r;
  }

  iterator begin() const { return iterator(this->_slist, true); }

  iterator end() const { return iterator(this->_slist, false); }

  std::size_t size() const { return this->_slist->size(); }

  void write(std::ostream& o) {
    o << "{";
    for (iterator it = this->begin(); it != this->end();) {
      Element e = *it;
      o << e;
      ++it;
      if (it != end()) {
        o << "; ";
      }
    }
    o << "}";
  }

}; // class collection

} // namespace ikos

#endif // IKOS_COMMON_HPP
