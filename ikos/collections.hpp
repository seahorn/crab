/*******************************************************************************
 * Basic bag-like collection based on a singly linked list.
 ******************************************************************************/


#ifndef IKOS_COLLECTIONS_HPP
#define IKOS_COLLECTIONS_HPP

#include <iostream>
#include <boost/container/slist.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <ikos/common.hpp>

namespace ikos {

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

} // namespace ikos

#endif // IKOS_COLLECTIONS_HPP
