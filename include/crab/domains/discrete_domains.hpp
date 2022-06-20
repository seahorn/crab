#pragma once

/** Several implementations of set-based domains:
 * 
 * - discrete_domain uses patricia trees to implement sets. As a
 *   constraint, it requires the set element to be indexable. Given a
 *   set of elements E, the discrete domain is a finite lattice that
 *   consists of all subsets of E. For instance if E={x,y,z}:
 * 
 *              {x,y,z}    least precise element
 *              /  |   \
 *          {x,y} {x,z} {y,z} 
 *           \   /    /  /
 *            \ / \  / \/
 *            {x} {y} {z}
 *              \  | / 
 *                { }      most precise element
 *
 * - set_domain has the same functionality than discrete_domain but it
 *   uses std::set. It should be used only when the set element cannot
 *   be indexable because set_domain is slower than discrete_domain.
 * 
 * - dual_set_domain: it takes a set-based domain (discrete_domain or
 *   set_domain) and invert it. That is, 
 *
 *                { }      least precise element
 *               / |  \
 *            {x} {y} {z}
 *            / \ / \ /	\
 *          {x,y} {x,z} {y,z}
 *             \   |   /
 *              {x,y,z}    most precise element
 *
 * - discrete_pair_domain is a set of pairs (key,value) where key must
 *   be indexable.
 **/
#include <crab/domains/patricia_trees.hpp>
#include <crab/support/debug.hpp>

#include "boost/range/iterator_range_core.hpp"
#include "boost/optional.hpp"

namespace ikos {

/*******************************************************************************
 *
 * An implementation of discrete domains based on Patricia trees.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
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
  
template <typename Element> class discrete_domain {

private:
  using ptset_t = patricia_tree_set<Element>;
  using discrete_domain_t = discrete_domain<Element>;

public:
  using element_t = Element;
  using iterator = typename ptset_t::iterator;
  
private:
  bool m_is_top;
  ptset_t m_set;

private:
  discrete_domain(bool is_top) : m_is_top(is_top) {}
  discrete_domain(ptset_t set) : m_is_top(false), m_set(set) {}

public:
  // Return an empty set
  static discrete_domain_t bottom() { return discrete_domain_t(false); }
  // Return a set with all variables
  static discrete_domain_t top() { return discrete_domain_t(true); }

public:
  // Default constructor creates an empty set rather than top
  discrete_domain() : m_is_top(false) {}

  discrete_domain(const discrete_domain_t &other) = default;
  discrete_domain(discrete_domain_t &&other) = default;
  discrete_domain_t &operator=(const discrete_domain_t &other) = default;
  discrete_domain_t &operator=(discrete_domain_t &&other) = default;  

  discrete_domain(Element s) : m_is_top(false), m_set(s) {}

  template <typename Iterator>
  discrete_domain(Iterator eIt, Iterator eEt) : m_is_top(false) {
    for (auto e : boost::make_iterator_range(eIt, eEt)) {
      m_set += e;
    }
  }

  bool is_top() const { return m_is_top; }

  bool is_bottom() const { return (!m_is_top && m_set.empty()); }

  bool operator<=(const discrete_domain_t &other) const {
    return other.m_is_top || (!m_is_top && m_set <= other.m_set);
  }

  bool operator==(const discrete_domain_t &other) const {
    return (m_is_top && other.m_is_top) || (m_set == other.m_set);
  }

  void operator|=(const discrete_domain_t &other) { *this = *this | other; }

  discrete_domain_t operator|(const discrete_domain_t &other) const {
    if (m_is_top || other.m_is_top) {
      return discrete_domain_t(true);
    } else {
      return discrete_domain_t(m_set | other.m_set);
    }
  }

  discrete_domain_t operator&(const discrete_domain_t &other) const {
    if (is_bottom() || other.is_bottom()) {
      return discrete_domain_t::bottom();
    } else if (is_top()) {
      return other;
    } else if (other.is_top()) {
      return *this;
    } else {
      return discrete_domain_t(m_set & other.m_set);
    }
  }

  discrete_domain_t operator||(const discrete_domain_t &other) const {
    return operator|(other);
  }

  discrete_domain_t operator&&(const discrete_domain_t &other) const {
    return operator&(other);
  }

  discrete_domain_t &operator+=(Element s) {
    if (!m_is_top) {
      m_set += s;
    }
    return *this;
  }

  template <typename Range> discrete_domain_t &operator+=(Range es) {
    if (!m_is_top) {
      for (auto e : es) {
        m_set += e;
      }
    }
    return *this;
  }

  discrete_domain_t operator+(Element s) {
    discrete_domain_t r(*this);
    r.operator+=(s);
    return r;
  }

  template <typename Range> discrete_domain_t operator+(Range es) {
    discrete_domain_t r(*this);
    r.operator+=(es);
    return r;
  }

  discrete_domain_t &operator-=(Element s) {
    if (!m_is_top) {
      m_set -= s;
    }
    return *this;
  }

  template <typename Range> discrete_domain_t &operator-=(Range es) {
    if (!m_is_top) {
      for (auto e : es) {
        m_set -= e;
      }
    }
    return *this;
  }

  discrete_domain_t operator-(Element s) {
    discrete_domain_t r(*this);
    r.operator-=(s);
    return r;
  }

  template <typename Range> discrete_domain_t operator-(Range es) {
    discrete_domain_t r(*this);
    r.operator-=(es);
    return r;
  }

  bool contain(Element e) {
    if (is_bottom()) {
      return false;
    } else if (is_top()) {
      return true;
    } else {
      return m_set[e];
    }
  }
  
  void rename(const std::vector<Element> &from, const std::vector<Element> &to) {
    if (is_top() || is_bottom()) {
      return;
    }
    if (from.size() != to.size()) {
      CRAB_ERROR("discrete_domain::rename with input vectors of different sizes");
    }

    for(unsigned i=0, sz=from.size(); i<sz; ++i) {
      if (from[i] == to[i]) {
	continue;
      }
      if (contain(from[i])) {
	this->operator-=(from[i]);
	this->operator+=(to[i]);
      }
    }
  }
  
  std::size_t size() const {
    if (m_is_top) {
      assert(false);
      CRAB_ERROR("Size for discrete domain TOP is undefined");
    } else {
      return m_set.size();
    }
  }

  iterator begin() const {
    if (m_is_top) {
      assert(false);      
      CRAB_ERROR("Iterator for discrete domain TOP is undefined");
    } else {
      return m_set.begin();
    }
  }

  iterator end() const {
    if (m_is_top) {
      assert(false);      
      CRAB_ERROR("Iterator for discrete domain TOP is undefined");
    } else {
      return m_set.end();
    }
  }

  void write(crab::crab_os &o) const {
    if (m_is_top) {
      o << "{...}";
    } else if (m_set.empty()) {
      o << "{}";
    } else {
      o << m_set;
    }
  }

}; // class discrete_domain

template <typename Elem>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const discrete_domain<Elem> &d) {
  d.write(o);
  return o;
}

} // namespace ikos


namespace crab {
namespace domains {
  
//===================================================================//  
//        The NOSA license does not apply to this code  
//===================================================================//  


/**
 * This domain is semantically equivalent to discrete_domain but it
 * uses std::set instead of patricia tries as underlying
 * datastructure. Because of this, this domain is slower than
 * discrete_domain so the only reason to use it is when Element cannot
 * be a subclass of indexable.
**/

template <class Element, class Compare>
class set_domain {
private:
  using set_t = std::set<Element, Compare>;
  using set_domain_t = set_domain<Element, Compare>;

public:
  using element_t = Element;  
  using iterator = typename set_t::const_iterator;

private:
  bool m_is_top;
  set_t m_set;

  set_domain(bool is_top) : m_is_top(is_top) {}
  set_domain(set_t set) : m_is_top(false), m_set(set) {}

public:
  // Default constructor creates an empty set rather than top
  set_domain() : m_is_top(false) {}
  set_domain(Element s) : m_is_top(false) {
    m_set.insert(s);
  }  
  set_domain(const set_domain_t &other) = default;
  set_domain(set_domain_t &&other) = default;
  set_domain_t &operator=(const set_domain_t &other) = default;
  set_domain_t &operator=(set_domain_t &&other) = default;  

  // Return an empty set
  static set_domain_t bottom() { return set_domain_t(false); }
  // Return a set with all variables
  static set_domain_t top() { return set_domain_t(true); }
  
  bool is_top() const { return m_is_top; }
  bool is_bottom() const { return (!m_is_top && m_set.empty()); }

  bool operator<=(const set_domain_t &other) const {
    return other.m_is_top ||
      (!m_is_top && std::includes(other.m_set.begin(), other.m_set.end(),
				  m_set.begin(), m_set.end(),
				  [](const Element &e1, const Element &e2) {
				    Compare cmp;
				    return cmp(e1,e2);
				  }));
  }

  bool operator==(const set_domain_t &other) const {
    return (m_is_top && other.m_is_top) || (m_set == other.m_set);
  }

  void operator|=(const set_domain_t &other) {
    if (is_top()) {
      return;
    } else if (other.is_top()) {
      *this = other;
    } else {
      m_set.insert(other.begin(), other.end());
    }
  }

  set_domain_t operator|(const set_domain_t &other) const {
    if (is_top() || other.is_top()) {
      return set_domain_t::top();
    } else {
      set_domain_t res(*this);
      res.m_set.insert(other.begin(), other.end());
      return res;
    }
  }

  set_domain_t operator&(const set_domain_t &other) const {
    if (is_bottom() || other.is_bottom()) {
      return set_domain_t::bottom();
    } else if (is_top()) {
      return other;
    } else if (other.is_top()) {
      return *this;
    } else {
      set_t s;
      std::set_intersection(m_set.begin(), m_set.end(),
			    other.m_set.begin(), other.m_set.end(),
			    std::inserter(s, s.end()),
			    [](const Element &e1, const Element &e2) {
			      Compare cmp;
			      return cmp(e1,e2);
			    });
      return set_domain_t(s);
    }
  }

  set_domain_t operator||(const set_domain_t &other) const {
    return operator|(other);
  }

  set_domain_t operator&&(const set_domain_t &other) const {
    return operator&(other);
  }

  set_domain_t &operator+=(Element s) {
    if (!is_top()) {
      m_set.insert(s);
    }
    return *this;
  }

  set_domain_t operator+(Element s) {
    set_domain_t r(*this);
    r.operator+=(s);
    return r;
  }

  set_domain_t &operator-=(Element s) {
    if (!is_top()) {
      m_set.erase(s);
    }
    return *this;
  }

  set_domain_t operator-(Element s) {
    set_domain_t r(*this);
    r.operator-=(s);
    return r;
  }

  bool contain(Element e) {
    if (is_bottom()) {
      return false;
    } else if (is_top()) {
      return true;
    } else {
      return m_set.count(e) > 0;
    }
  }
  
  void rename(const std::vector<Element> &from, const std::vector<Element> &to) {
    if (is_top() || is_bottom()) {
      return;
    }
    if (from.size() != to.size()) {
      CRAB_ERROR("set domain::rename with input vectors of different sizes");
    }

    for(unsigned i=0, sz=from.size(); i<sz; ++i) {
      if (from[i] == to[i]) {
	continue;
      }
      if (contain(from[i])) {
	this->operator-=(from[i]);
	this->operator+=(to[i]);
      }
    }
  }
  
  std::size_t size() const {
    if (is_top()) {
      assert(false);
      CRAB_ERROR("Size for set domain TOP is undefined");
    } else {
      return m_set.size();
    }
  }

  iterator begin() const {
    if (is_top()) {
      assert(false);      
      CRAB_ERROR("Iterator for set domain TOP is undefined");
    } else {
      return m_set.begin();
    }
  }

  iterator end() const {
    if (is_top()) {
      assert(false);      
      CRAB_ERROR("Iterator for set domain TOP is undefined");
    } else {
      return m_set.end();
    }
  }

  void write(crab::crab_os &o) const {
    if (is_top()) {
      o << "{...}";
    } else if (is_bottom()) {
      o << "{}";
    } else {
      o << "{";
      for (auto it = m_set.begin(), et = m_set.end(); it!=et; ) {
	o << *it;
	++it;
	if (it != et) {
	  o << ",";
	}
      }
      o << "}";
    }
  }

}; // class set_domain

template<class Element, class Compare>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const set_domain<Element,Compare> &d) {
  d.write(o);
  return o;
}

// Dual of set_domain/discrete_domain: the larger is the set, the more
// precise is the domain.
template<class Set>
class dual_set_domain {
  using dual_set_domain_t = dual_set_domain<Set>;
       
public:
  using set_domain_t = Set;   
  using element_t = typename set_domain_t::element_t;
  using iterator = typename set_domain_t::iterator;

private:
    set_domain_t m_set;
public:
  
  dual_set_domain()
    : m_set(set_domain_t::bottom()) /*top by default*/ {}
  dual_set_domain(set_domain_t s)
    : m_set(s) {}

  dual_set_domain(const dual_set_domain_t &other) = default;
  dual_set_domain(dual_set_domain_t &&other) = default;  
  dual_set_domain_t&operator=(const dual_set_domain_t &other) = default;
  dual_set_domain_t&operator=(dual_set_domain_t &&other) = default;  

  static dual_set_domain_t bottom() { return set_domain_t::top(); }
  static dual_set_domain_t top() { return set_domain_t::bottom(); }
  bool is_top() const { return m_set.is_bottom(); }
  bool is_bottom() const { return m_set.is_top(); }

  bool operator<=(const dual_set_domain_t &other) const {
    if (other.is_top() || is_bottom()) {
      return true;
    } else {
      return other.m_set <= m_set;
    }
  }

  bool operator==(const dual_set_domain_t &other) const {
    return (*this <= other && other <= *this);
  }

  void operator|=(const dual_set_domain_t &other) {
    m_set = m_set & other.m_set;
  }

  dual_set_domain_t operator|(const dual_set_domain_t &other) const {
    return dual_set_domain_t(m_set & other.m_set);
  }

  dual_set_domain_t operator&(const dual_set_domain_t &other) const {
    return dual_set_domain_t(m_set | other.m_set);
  }

  dual_set_domain_t operator||(const dual_set_domain_t &other) const {
    return this->operator|(other);
  }

  dual_set_domain_t operator&&(const dual_set_domain_t &other) const {
    return this->operator&(other);
  }

  dual_set_domain_t &operator+=(const element_t &c) {
    m_set += c;
    return *this;
  }
  dual_set_domain_t &operator-=(const element_t &c) {
    m_set -= c;
    return *this;
  }

  bool at(const element_t &e) const{
    dual_set_domain_t s(e);
    return (s <= *this);
  }
  
  std::size_t size() { return m_set.size(); }

  iterator begin() const  { return m_set.begin(); }
  iterator end() const { return m_set.end(); }
  
  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      m_set.write(o);
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
				   const dual_set_domain_t &dom) {
    dom.write(o);
    return o;
  }
}; // end dual_set_domain


/*
 * Represent sets of pairs (Key,Value).
 * 
 * The pair must consist of a Key (must inherit from indexable class)
 * and a Value that can be a generic abstract domain. 
 * 
 * Bottom means empty set rather than failure.
 */ 
template <typename Key, typename Value> class discrete_pair_domain {

private:
  using patricia_tree_t = ikos::patricia_tree<Key, Value>;
  using unary_op_t = typename patricia_tree_t::unary_op_t;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

public:
  using discrete_pair_domain_t = discrete_pair_domain<Key, Value>;
  using iterator = typename patricia_tree_t::iterator;
  using key_type = Key;
  using value_type = Value;

private:
  bool m_is_top;
  patricia_tree_t m_tree;

  static patricia_tree_t apply_operation(binary_op_t &o,
                                         const patricia_tree_t &t1,
                                         const patricia_tree_t &t2,
                                         bool &is_bottom) {
    patricia_tree_t res(t1);
    is_bottom = res.merge_with(t2, o);
    return res;
  }

  discrete_pair_domain(patricia_tree_t t) : m_is_top(false), m_tree(t) {}

  discrete_pair_domain(bool b) : m_is_top(b) {}

  class join_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator|(y);
      if (z.is_top()) {
        return {false, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    }
    virtual bool default_is_absorbing() override { return false; }
  }; // class join_op

  class meet_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator&(y);
      if (z.is_bottom()) {
        return {true, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };
    virtual bool default_is_absorbing() override { return true; }
  }; // class meet_op

  class domain_po : public partial_order_t {
    virtual bool leq(const Value &x, const Value &y) override { return x.operator<=(y); }
    virtual bool default_is_top() override { return false; }
  }; // class domain_po

public:
  static discrete_pair_domain_t top() {
    return discrete_pair_domain_t(true);
  }

  static discrete_pair_domain_t bottom() {
    return discrete_pair_domain_t(false);
  }
  
  discrete_pair_domain() : m_is_top(false), m_tree(patricia_tree_t()) {}


  discrete_pair_domain(const discrete_pair_domain_t &o) = default;
  discrete_pair_domain(discrete_pair_domain_t &&o) = default;
  discrete_pair_domain_t &operator=(const discrete_pair_domain_t &o)  = default;
  discrete_pair_domain_t &operator=(discrete_pair_domain_t &&o)  = default;
  
  iterator begin() const {
    if (is_top()) {
      CRAB_ERROR("discrete_pair_domain: trying to invoke iterator on top");
    } else {
      return m_tree.begin();
    }
  }

  iterator end() const {
    if (is_top()) {
      CRAB_ERROR("discrete_pair_domain: trying to invoke iterator on top");
    } else {
      return m_tree.end();
    }
  }

  bool is_top() const { return m_is_top; }

  bool is_bottom() const { return (!is_top() && m_tree.empty()); }

  bool operator<=(const discrete_pair_domain_t &o) const {
    domain_po po;
    return (o.is_top() || (!is_top() && (m_tree.leq(o.m_tree, po))));
  }

  discrete_pair_domain_t
  operator|(const discrete_pair_domain_t &o) const {
    if (is_top() || o.is_top()) {
      return discrete_pair_domain_t::top();
    } else {
      join_op op;
      bool is_bottom;
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree, is_bottom);
      return discrete_pair_domain_t(std::move(res));
    }
  }

  discrete_pair_domain_t
  operator&(const discrete_pair_domain_t &o) const {
    if (is_top()) {
      return o;
    } else if (o.is_top()) {
      return *this;
    } else {
      meet_op op;
      bool is_bottom;
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree, is_bottom);
      if (is_bottom) {
        return discrete_pair_domain_t::bottom();
      } else {
        return discrete_pair_domain_t(std::move(res));
      }
    }
  }

  void set(Key k, Value v) {
    if (!is_top()) {
      m_tree.insert(k, v);
    }
  }

  discrete_pair_domain_t &operator-=(Key k) {
    if (!is_top()) {
      m_tree.remove(k);
    }
    return *this;
  }

  Value operator[](Key k) {
    if (is_top())
      return Value::top();
    else {
      boost::optional<Value> v = m_tree.lookup(k);
      if (v) {
        return *v;
      } else {
        return Value::bottom();
      }
    }
  }
  
  void write(crab::crab_os &o) const {
    if (is_top()) {
      o << "{...}";
    } else if (is_bottom()) {
      o << "{}";
    } else {
      o << "{";
      for (typename patricia_tree_t::iterator it = m_tree.begin();
           it != m_tree.end();) {
        Key k = it->first;
        k.write(o);
        o << " -> ";
        Value v = it->second;
        v.write(o);
        ++it;
        if (it != m_tree.end()) {
          o << "; ";
        }
      }
      o << "}";
    }
  }
}; // class discrete_pair_domain

template <typename Key, typename Value>
inline crab_os &operator<<(crab_os &o, const discrete_pair_domain<Key, Value> &d) {
  d.write(o);
  return o;
}

} //end namespace domains
} //end namespace crab
