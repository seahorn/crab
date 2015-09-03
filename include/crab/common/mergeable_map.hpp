/*******************************************************************************
 * Implementation of a mergeable map based on patricia trees
 ******************************************************************************/

#ifndef MERGEABLE_MAP_HPP
#define MERGEABLE_MAP_HPP

#include <boost/optional.hpp>

#include <crab/common/types.hpp>
#include <crab/domains/patricia_trees.hpp>

using namespace ikos;

namespace crab {

   namespace domains {

      template<typename Key, typename Value>
      class merge_op_check_equal: public patricia_tree< Key, Value >::binary_op_t {
        boost::optional< Value > apply(Value x, Value y) {
          if (x == y)
            return x;
          else
            CRAB_ERROR("mergeable_map: merging a key with two different values");
        };
        bool default_is_absorbing() { return false; }
      }; 
   
      template<typename Key, typename Value>
      class merge_op_first: public patricia_tree< Key, Value >::binary_op_t {
        boost::optional< Value > apply(Value x, Value y)  {
          return x;
        };
        bool default_is_absorbing() { return false; } 
      }; 
   
      template<typename Key, typename Value>
      class merge_op_second: public patricia_tree< Key, Value >::binary_op_t {  
        boost::optional< Value > apply(Value x, Value y)  {
          return y;
        };
        bool default_is_absorbing() { return false; }
      }; 
   
     template < typename Key, typename Value, 
                typename MergeOp = merge_op_check_equal <Key, Value> >
     class mergeable_map: public writeable {
       
      private:
       typedef patricia_tree< Key, Value > patricia_tree_t;
       typedef typename patricia_tree_t::binary_op_t binary_op_t;
  
      public:
       typedef mergeable_map< Key, Value > mergeable_map_t;
       typedef typename patricia_tree_t::iterator iterator;
       
      private:
       patricia_tree_t _tree;
       
       static patricia_tree_t do_union(patricia_tree_t t1, patricia_tree_t t2) {
         MergeOp o;
         t1.merge_with(t2, o);
         return t1;
       }
       
       mergeable_map(patricia_tree_t t): _tree(t) { }
       
      public:
       
       mergeable_map(): _tree(patricia_tree_t()) { }
  
       mergeable_map(const mergeable_map_t& e): writeable(), _tree(e._tree) { }
       
       mergeable_map_t& operator=(mergeable_map_t e) {
         _tree = e._tree;
         return *this;
       }
       
       iterator begin() {
         return _tree.begin();
       }
       
       iterator end() {
         return _tree.end();
       }
       
       std::size_t size(){
         return _tree.size();
       }
       
       mergeable_map_t operator|(mergeable_map_t e) {
         mergeable_map_t u(do_union(_tree, e._tree));
         return u;
       }
       
       void set(Key k, Value v) {
         _tree.insert(k, v);
       }
       
       mergeable_map_t& operator-=(Key k) {
         _tree.remove(k);
         return *this;
       }
       
       boost::optional<Value> operator[](Key k) {
         return _tree.lookup(k);
       }
       
       void clear() {
         _tree = patricia_tree_t();
       }
       
       void write(std::ostream& o) {
         o << "{";
         for (auto it = _tree.begin(); it != _tree.end(); ) {
           Key k = it->first;
           k.write(o);
           o << " -> ";
           Value v = it->second;
           o << v;
           ++it;
           if (it != _tree.end()) {
             o << "; ";
           }
         }
         o << "}";
       }    
     }; // class mergeable_map

   } // end namespace domains
} // end namespace crab

#endif /*MERGEABLE_MAP_HPP*/
