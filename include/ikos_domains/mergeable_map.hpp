/*******************************************************************************
 * Implementation of a mergeable map based on patricia trees
 ******************************************************************************/

#ifndef MERGEABLE_MAP_HPP
#define MERGEABLE_MAP_HPP

#include <iostream>
#include <boost/optional.hpp>
#include <ikos_domains/patricia_trees.hpp>
#include <ikos_domains/common.hpp>

namespace ikos {

template < typename Key, typename Value >
class mergeable_map: public writeable {

 private:
  typedef patricia_tree< Key, Value > patricia_tree_t;
  typedef typename patricia_tree_t::binary_op_t binary_op_t;
  
 public:
  typedef mergeable_map< Key, Value > mergeable_map_t;
  typedef typename patricia_tree_t::iterator iterator;
  
 private:
  patricia_tree_t _tree;
  
 private:

  class union_op: public binary_op_t {
    
    boost::optional< Value > apply(Value x, Value y) 
    {
      if (x == y)
        return boost::optional< Value >(x);
      else
        throw error("mergeable_map: merging a key with two different values");
    };

    bool default_is_absorbing() 
    {
      return false;
    }

  }; // class union_op
  
 private:

  static patricia_tree_t do_union(patricia_tree_t t1, patricia_tree_t t2) {
    union_op o;
    t1.merge_with(t2, o);
    return t1;
  }
  
 private:
  mergeable_map(patricia_tree_t t): _tree(t) { }
  
 public:
  mergeable_map(): _tree(patricia_tree_t()) { }
  
  mergeable_map(const mergeable_map_t& e): writeable(), _tree(e._tree) { }
  
  mergeable_map_t& operator=(mergeable_map_t e) {
    this->_tree = e._tree;
    return *this;
  }
  
  iterator begin() {
    return this->_tree.begin();
  }
  
  iterator end() {
    return this->_tree.end();
  }
  
  std::size_t size(){
    return this->_tree.size();
  }
  
  mergeable_map_t operator|(mergeable_map_t e) {
    mergeable_map_t u(do_union(this->_tree, e._tree));
    return u;
  }
  
  void set(Key k, Value v) {
    this->_tree.insert(k, v);
  }
  
  mergeable_map_t& operator-=(Key k) {
    this->_tree.remove(k);
    return *this;
  }
  
  boost::optional<Value> operator[](Key k) {
    return this->_tree.lookup(k);
  }

  void clear() {
    this->_tree = patricia_tree_t();
  }
  
  std::ostream& write(std::ostream& o) {
    o << "{";
    for (typename patricia_tree_t::iterator it = this->_tree.begin(); it != this->_tree.end(); ) {
      it->first.write(o);
      o << " -> ";
      o << it->second;
      ++it;
      if (it != this->_tree.end()) {
        o << "; ";
      }
    }
    o << "}";
    return o;
  }    
}; // class mergeable_map

} // end namespace 

#endif /*MERGEABLE_MAP_HPP*/
