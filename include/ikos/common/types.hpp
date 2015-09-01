/*******************************************************************************
 *
 * Basic type definitions.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 * Contributors: Jorge A. Navas (jorge.a.navaslaserna@nasa.gov)
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
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/slist.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace ikos {

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

#define IKOS_ERROR(...)              \
    do {                             \
      ___print___(__VA_ARGS__);      \
      std::cerr << "\n";             \
      std::exit (EXIT_FAILURE);      \
    } while (0)

// Default definition but it can be redefined by dbg.hpp
#define IKOS_DEBUG(...)

} // end namespace

namespace ikos 
{
  
  // Numerical type for indexed objects
  typedef uint64_t index_t;
  
  // Interface for writeable objects
  class writeable {

  public:
    virtual std::ostream& write(std::ostream& o) = 0;
    
    virtual ~writeable() { }

    friend std::ostream& operator<<(std::ostream& o, writeable& x) {
      x.write (o);
      return o;
    }
    
  }; // class writeable



  // Exception for IKOS internal errors
  class error: public writeable {
    
  protected:
    std::string msg;
    
  private:
    error();
    
  public:
    error(std::string msg_): msg(msg_) { }
    
    std::string message() {
      return msg;
    }

    std::ostream& write(std::ostream& o) {
      o << message();
      return o;
    }
    
  }; // class error


  // Container for typed variables
  template< typename Type, typename VariableName >
  class variable: public writeable {
    
   public:
    typedef variable< Type, VariableName > variable_t;
    typedef typename VariableName::index_t index_t;

   private:
    VariableName _n;
    
   public:
    variable(VariableName n): writeable (), _n(n) { }

    variable(const variable_t& v): writeable(), _n(v._n) { }
    
    variable_t& operator=(const variable_t &o) {
      if (this != &o) {
        this->_n = o._n;
      }
      return *this;
    }

    VariableName name() const {
      return _n;
    }
    
    index_t index() const {
      return _n.index();
    }

    bool operator==(const variable_t& o) const {
      return _n.index () == o._n.index ();
    }

    bool operator<(const variable_t& o) const {
      return _n.index () < o._n.index ();
    }

    std::ostream& write(std::ostream& o) {
      o << _n;
      return o;
    }

    friend index_t hash_value (variable_t  v) {
      return v.index ();
    }

    friend std::ostream& operator<<(std::ostream& o, variable_t  v) {
      v.write(o);
      return o;
    }

  }; // class variable

  template< typename Element >
  class collection: public writeable {

  public:
    typedef collection< Element > collection_t;

  private:
    typedef boost::container::slist< Element > slist_t;
    typedef boost::shared_ptr< slist_t > slist_ptr;
    
  private:
    slist_ptr _slist;

  public:
    class iterator: public boost::iterator_facade< iterator
						   , Element
						   , boost::forward_traversal_tag
						   , Element
						   > {
      
      friend class boost::iterator_core_access;
      
    private:
      typename slist_t::iterator _it;
      slist_ptr _l;
      
    public:
      iterator(slist_ptr l, bool b): _it(b ? l->begin() : l->end()), _l(l) { }
      
    private:
      void increment() { 
	++(this->_it);
      }
      
      bool equal(const iterator& other) const {
	return (this->_l == other._l && this->_it == other._it);
      }
      
      Element dereference() const {
	if (this->_it != this->_l->end()) {
	  return *(this->_it);
	} else {
	  throw error("Collection: trying to dereference an empty iterator");
	}
      }
      
    }; // class iterator

  public:
    collection(): _slist(slist_ptr(new slist_t)) { }

    collection(const collection_t& c): writeable(), _slist(c._slist) { }
    
    collection_t& operator=(collection_t c) {
      this->_slist = c._slist;
      return *this;
    }

    collection_t& operator+=(Element e) {
      this->_slist->push_front(e);
      return *this;
    }

    collection_t& operator+=(collection_t c) {
      for (iterator it = c.begin(); it != c.end(); ++it) {
	this->_slist->push_front(*it);
      }
      return *this;
    }

    collection_t operator+(collection_t c) {
      collection_t r;
      r.operator+=(c);
      r.operator+=(*this);
      return r;
    }

    iterator begin() {
      return iterator(this->_slist, true);
    }

    iterator end() {
      return iterator(this->_slist, false);
    }

    std::size_t size() {
      return this->_slist->size();
    }
    
    std::ostream& write(std::ostream& o) {
      o << "{";
      for (iterator it = this->begin(); it != this->end(); ) {
	Element e = *it;
	o << e;
	++it;
	if (it != end()) {
	  o << "; ";
	}
      }
      o << "}";
      return o;
    }

  }; // class collection

  // Enumeration type for basic arithmetic operations
  typedef enum {
    OP_ADDITION,
    OP_SUBTRACTION,
    OP_MULTIPLICATION,
    OP_DIVISION
  } operation_t;

  inline std::ostream& operator<<(std::ostream&o, operation_t op) {
    switch (op) {
      case OP_ADDITION: o << "+"; break;
      case OP_SUBTRACTION: o << "-"; break;
      case OP_MULTIPLICATION: o << "*"; break;
      default: o << "/"; break;
    }
    return o;
  }

} // namespace ikos

#endif // IKOS_COMMON_HPP
