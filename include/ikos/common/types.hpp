/*******************************************************************************
 * Common types and classes used in IKOS.
 ******************************************************************************/


#ifndef IKOS_COMMON_HPP
#define IKOS_COMMON_HPP

#include <stdint.h>
#include <string>
#include <iostream>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/slist.hpp>
#include <boost/iterator/iterator_facade.hpp>

#define IKOS_ERROR(msg)                             \
    do {                                            \
      std::cout << msg << "\n";                     \
      exit (EXIT_FAILURE);                          \
    } while (0)

namespace ikos 
{
  

  // Numerical type for indexed objects
  typedef uint64_t index_t;
  
  // Interface for writeable objects
  class writeable {

  public:
    virtual std::ostream& write(std::ostream& o) = 0;

    virtual ~writeable() { }

  }; // class writeable

  inline std::ostream& operator<<(std::ostream& o, writeable& x) {
    return (x.write(o));
  }


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


  // Container data structure for typed variables
  template< typename Type, typename VariableName >
  class variable: public writeable {
    
  public:
    typedef variable< Type, VariableName > variable_t;

  public:
    VariableName _n;
    
  public:
    variable(VariableName n): _n(n) { }

    variable(const variable_t& v): writeable(), _n(v._n) { }
    
    variable_t& operator=(variable_t v) {
      this->_n = v._n;
      return *this;
    }

    VariableName name() const {
      return _n;
    }
    
    index_t index() {
      return _n.index();
    }

    bool operator<(const variable_t& v) const {
      variable_t v1 = const_cast< variable_t& >(*this);
      variable_t v2 = const_cast< variable_t& >(v);
      return v1._n.index() < v2._n.index();
    }

    std::ostream& write(std::ostream& o) {
      o << _n;
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
  
} // namespace ikos

#endif // IKOS_COMMON_HPP
