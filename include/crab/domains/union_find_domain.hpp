#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <set>

/**
 * A standard union-find equipped with lattice operations (inclusion,
 * lub, meet, etc).
 **/

namespace crab {
namespace domains {

template <class Domain> class equivalence_class {
private:
  std::size_t m_rank;
  Domain m_val;

public:
  explicit equivalence_class(Domain val) : m_rank(0), m_val(val) {}

  std::size_t &get_rank() { return m_rank; }

  std::size_t get_rank() const { return m_rank; }

  Domain &get_domain() { return m_val; }

  const Domain &get_domain() const { return m_val; }

}; // end class equivalence_class

/*
  Partition a set of elements into equivalence classes. Apart from
  lattice operations, the domain provides standard union-find
  operations:

  - make(v)
  - union(v1, v2)
  - find(v)

  In addition, each equivalence class has associated an abstract state
  of type Domain.

  The method "top()" returns an union-find in which all elements are
  in the same equivalence class. Since the domain does not know the
  universe of elements, the method "is_top()" only succeeds if the
  union_find_domain state was created by calling directly "top()".
  The method "bottom()" represents an unreachable/failure abstract
  state. Do not use "top()" if you intend to generate an empty
  union-find. Instead, use the constructor without parameters to
  produce an empty union-find.
*/

// Important: when two equivalence classes are merged, there are two
// options: either apply union or intersection semantics.
template<class Domain>  
struct uf_merge_semantics_t {
  virtual bool is_bottom_absorbing() const = 0;
  virtual Domain apply(const Domain &d1, const Domain &d2) const = 0;
};

template<class Domain>    
struct uf_union_semantics_t: public uf_merge_semantics_t<Domain> {
  virtual bool is_bottom_absorbing() const override {
    return false;
  }
  virtual Domain apply(const Domain &d1, const Domain &d2) const override {
    return d1 | d2;
  }
};
  
template<class Domain>    
struct uf_intersection_semantics_t: public uf_merge_semantics_t<Domain> {
  virtual bool is_bottom_absorbing() const override {
    return true;
  }  
  virtual Domain apply(const Domain &d1, const Domain &d2) const override {
    return d1 & d2;
  }
};
  
template<class Element, class Domain,
	 // By default, we apply union semantics.
	 class MergeSemantics = uf_union_semantics_t<Domain>>
class union_find_domain {
public:
  using element_t = Element;
  using domain_t = Domain;
  enum class lattice_val { bottom, top, neither_top_nor_bot };    
private:
  using union_find_domain_t = union_find_domain<Element, Domain, MergeSemantics>;
  using parents_map_t = std::unordered_map<element_t, element_t>;
public:
  using element_set_t = std::set<element_t>;  
  using equivalence_class_elems_t = std::unordered_map<element_t, element_set_t>;
  using equivalence_class_t = equivalence_class<Domain>;  
private:  
  using classes_map_t = std::unordered_map<element_t, equivalence_class_t>;

  // immediate link to the parent
  parents_map_t m_parents;
  // each root has its equivalence class
  classes_map_t m_classes;
  lattice_val m_val;

  void clear() {
    m_parents.clear();
    m_classes.clear();
  }
 
  // Return the representative after merging all elements in elems.
  // Pre-condition: if dom != none then forall v \in elems:: contains(v)
  boost::optional<element_t>
  merge_elems(const element_set_t &elems,
             boost::optional<domain_t> dom = boost::none) {
    auto it = elems.begin();
    auto et = elems.end();
    if (it == et) {
      return boost::none;
    }
    element_t v = *it;
    if (dom && !contains(v)) {
      make(v, *dom);
      if (is_bottom()) {
	return boost::none;
      }
    }
    ++it;

    if (it == et) {
      return (contains(v) ? boost::optional<element_t>(find(v)) : boost::none);
    }

    boost::optional<element_t> repr;
    for (; it != et; ++it) {
      if (dom && !contains(*it)) {
        make(*it, *dom);
	if (is_bottom()) {
	  return boost::none;
	}	
      }
      repr = join(v, *it);      
      if (!repr /*is_bottom()*/) {
	return boost::none;
      }      
    }
    return repr;
  }

  void print(crab_os &o) const {
    o << "({";
    for (auto it = m_parents.begin(), et = m_parents.end(); it != et;) {
      element_t key = it->first, value = it->second;
      o << key << " -> " << value;
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "}, {";
    for (auto it = m_classes.begin(), et = m_classes.end(); it != et;) {
      element_t rep = it->first;
      equivalence_class_t ec = it->second;
      o << rep << " -> " << ec.get_domain();
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "})";
  }

  void get_all_members(element_set_t &out) const {
    for (auto &kv: m_parents) {
      out.insert(kv.first);
    }
  }

  union_find_domain_t join_or_widening(const union_find_domain_t &o, bool is_join) const {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() | o.is_top()) {
      union_find_domain_t res;
      return res;
    } else {
      CRAB_LOG("union-find",
               crab::outs() << (is_join ? "Join" : "Widen")
	                    << "\n\t" << *this << "\n\t" << o << "\n";);
      union_find_domain_t left(*this);
      union_find_domain_t right(o);

      // Keep only common elements      
      element_set_t left_elems;
      left.get_all_members(left_elems);
      for (auto v: left_elems) {
	if (!right.contains(v)) {
	  left.forget(v);
	}
      }
      element_set_t right_elems;
      right.get_all_members(right_elems);      
      for (auto v: right_elems) {
	if (!left.contains(v)) {
	  right.forget(v);
	}
      }
      
      equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
      equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();

      // Merge equivalence classes from the left to the right while
      // joining/widening the domains associated to the equivalence
      // classes.
      //
      // The merging on the right is needed so that right_dom is updated.
      for (auto &kv : left_equiv_classes) {
        boost::optional<element_t> right_repr =
	  right.merge_elems(kv.second);
        if (!right_repr) {
	  // this shouldn't happen
	  CRAB_ERROR("unexpected situation in join_or_widening 1");
	}
        domain_t &left_dom = left.get_equiv_class(kv.first).get_domain();
        const domain_t &right_dom =
            right.get_equiv_class(*right_repr).get_domain();
        left_dom = (is_join ? (left_dom | right_dom): (left_dom || right_dom));
      }

      // Merge equivalence classes from the right to the left while
      // joining/widening the domains associated to the equivalence
      // classes.
      for (auto &kv : right_equiv_classes) {
        boost::optional<element_t> left_repr =
	  left.merge_elems(kv.second);
        if (!left_repr) {
	   // this shouldn't happen
	  CRAB_ERROR("unexpected situation in join_or_widening 2");
	}
        domain_t &left_dom = left.get_equiv_class(*left_repr).get_domain();
        const domain_t &right_dom =
            right.get_equiv_class(kv.first).get_domain();
	left_dom = (is_join ? (left_dom | right_dom): (left_dom || right_dom));	
      }
      CRAB_LOG("union-find", crab::outs() << "Res=" << left << "\n";);
      return left;
    }
  }

  union_find_domain_t meet_or_narrowing(const union_find_domain_t &o, bool is_meet) const {
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    } else {
      /* TOIMPROVE: this operation is imprecise because we still merge
       * the equivalence classes.
       */
      union_find_domain_t left(*this);
      union_find_domain_t right(o);
      equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
      equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();

      // Merge equivalence classes from the right to the left while
      // applying meet/narrowing
      //
      // The merging on the right is needed so that right_dom is updated.
      for (auto &kv : left_equiv_classes) {
        domain_t &left_dom = left.get_equiv_class(kv.first).get_domain();
        // add elements on the right if they do not exist
        boost::optional<element_t> right_repr =
	  right.merge_elems(kv.second, left_dom);
        if (!right_repr) {
	  // this shouldn't happen
	  CRAB_ERROR("unexpected situation in meet_or_narrowing 1");
	}
        const domain_t &right_dom =
            right.get_equiv_class(*right_repr).get_domain();
        left_dom = (is_meet ? left_dom & right_dom: left_dom && right_dom);
        if (left_dom.is_bottom()) {
          left.set_to_bottom();
          return left;
        }
      }

      // Merge equivalence classes from the right to the left while
      // applying meet/narrowing
      for (auto &kv : right_equiv_classes) {
        const domain_t &right_dom =
            right.get_equiv_class(kv.first).get_domain();
        boost::optional<element_t> left_repr =
	  left.merge_elems(kv.second, right_dom);
        if (!left_repr) {
	   // this shouldn't happen
	  CRAB_ERROR("unexpected situation in meet_or_narrowing 2");
	}
        domain_t &left_dom = left.get_equiv_class(*left_repr).get_domain();
	left_dom = (is_meet ? left_dom & right_dom: left_dom && right_dom);
        if (left_dom.is_bottom()) {
          left.set_to_bottom();
          return left;
        }
      }
      return left;
    }
  }
 
  struct get_domain_from_ec:
    public std::unary_function<typename classes_map_t::value_type, Domain> {    
    const Domain& operator()(const typename classes_map_t::value_type &kv) const {
      return kv.second.get_domain();
    }
    Domain& operator()(typename classes_map_t::value_type &kv) const {
      return kv.second.get_domain();
    }     
  };  
public:
 
  using const_domain_iterator =
    boost::transform_iterator<get_domain_from_ec, typename classes_map_t::const_iterator>;
  using domain_iterator =
    boost::transform_iterator<get_domain_from_ec, typename classes_map_t::iterator>;
  using const_domain_range = boost::iterator_range<const_domain_iterator>;
  using domain_range = boost::iterator_range<domain_iterator>;

  // empty union-find 
  union_find_domain(lattice_val val = lattice_val::neither_top_nor_bot)
      : m_val(val) {}
  union_find_domain(const union_find_domain_t &o) = default;
  union_find_domain(union_find_domain_t &&o) = default;
  union_find_domain_t &operator=(const union_find_domain_t &o) = default;
  union_find_domain_t &operator=(union_find_domain_t &&o) = default;

  /* =============================================================
   * Standard union-find operations: make, find, and union plus some
   * extra helpers.
   * =============================================================
   */

  bool is_empty() const {
    return m_parents.empty();
  }
    
  // Pre-condition: !contains(v)
  void make(const element_t &v, Domain val) {
    if (is_bottom()) {
      CRAB_ERROR("calling union_find_domain::make on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("calling union_find_domain::make on top");
    }
    if (contains(v)) {
      CRAB_ERROR("element already exists when called union_find_domain::make");
    }

    MergeSemantics op;
    if (op.is_bottom_absorbing() && val.is_bottom()) {
      set_to_bottom();
    } else {
      m_parents.insert({v, v});
      m_classes.insert({v, equivalence_class_t(val)});
    }
  }

  // Pre-condition: contains(v)
  // NOTE: it is not "const" because it does path-compression
  element_t &find(const element_t &v) {
    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      assert(false);
      CRAB_ERROR("called union_find_domain::find on a non-existing element ",
                 v, " in ", *this);
    }
    element_t &parent = it->second;
    if (parent == v) {
      return parent;
    } else {
      return parent = find(parent);
    }
  }

  // Pre-condition: contains(v)
  // find operation without path-compression
  const element_t &find(const element_t &v) const {
    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      assert(false);
      CRAB_ERROR("called union_find_domain::find on a non-existing elem ",
                 v, " in ", *this);
    }
    const element_t &parent = it->second;
    if (parent == v) {
      return parent;
    } else {
      return find(parent);
    }
  }  

  // Pre-condition: contains(x) && contains(y)
  // Post-condition: merge x and y. If the result is not bottom then
  // it returns the representative of the new equivalence class,
  // otherwise, none.
  boost::optional<element_t> join(const element_t &x, const element_t &y) {
    element_t rep_x = find(x);
    element_t rep_y = find(y);
    MergeSemantics merge_op;
    if (rep_x != rep_y) {
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      equivalence_class_t &ec_y = m_classes.at(rep_y);
      domain_t &dom_x = ec_x.get_domain();
      domain_t &dom_y = ec_y.get_domain();
      if (ec_x.get_rank() > ec_y.get_rank()) {
        dom_x = merge_op.apply(dom_x, dom_y);
	if (merge_op.is_bottom_absorbing() && dom_x.is_bottom()) {
	  set_to_bottom();
	  return boost::none;
	}								   
        m_parents.at(rep_y) = rep_x;
        m_classes.erase(rep_y);
        return rep_x;
      } else {
        dom_y = merge_op.apply(dom_y, dom_x);
	if (merge_op.is_bottom_absorbing() && dom_y.is_bottom()) {
	  set_to_bottom();
	  return boost::none;
	} 
        m_parents.at(rep_x) = rep_y;
        if (ec_x.get_rank() == ec_y.get_rank()) {
          ec_y.get_rank()++;
        }
        m_classes.erase(rep_x);
        return rep_y;
      }
    } else {
      return rep_x;
    }
  }

  bool contains(const element_t &v) const {
    return m_parents.find(v) != m_parents.end();
  }

  // Pre-condition: contains(v)
  equivalence_class_t &get_equiv_class(const element_t &v) {
    return m_classes.at(find(v));
  }

  // Pre-condition: contains(v). No path-compression.
  const equivalence_class_t &get_equiv_class(const element_t &v) const {
    return m_classes.at(find(v));
  }

  void remove_equiv_class(const element_t &v) {
    if (!contains(v)) {
      return;
    }
    element_t rep_v = find(v);
    m_classes.erase(rep_v);

    equivalence_class_elems_t map = equiv_classes_elems();
    auto it = map.find(rep_v);
    if (it != map.end()) {
      for (auto v : it->second) {
        m_parents.erase(v);
      }
    }
  }

  /* Iterate over the domains associated to each equivalence class */
  const_domain_iterator begin_domains() const {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: begin_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: begin_domains() on top");
    }    
    return boost::make_transform_iterator(m_classes.begin(),
					  get_domain_from_ec());
  }
  
  const_domain_iterator end_domains() const {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: end_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: end_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.end(),
					  get_domain_from_ec());
  }

  domain_iterator begin_domains() {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: begin_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: begin_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.begin(),
					  get_domain_from_ec());
  }

  domain_iterator end_domains() {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: end_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: end_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.end(),
					  get_domain_from_ec());
  }    

  domain_range domains() {
    domain_iterator it = begin_domains();
    domain_iterator et = end_domains();
    return boost::make_iterator_range(it, et);
  }

  const_domain_range domains() const {
    const_domain_iterator it = begin_domains();
    const_domain_iterator et = end_domains();    
    return boost::make_iterator_range(it, et);
  }

  // Build a map from representative to a set with all the elements
  // in the equivalence class.
  equivalence_class_elems_t equiv_classes_elems() {
    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t rep = find(kv.second);
      auto it = res.find(rep);
      if (it == res.end()) {
        element_set_t s;
        s.insert(kv.first);
        res.insert({rep, s});
      } else {
        element_set_t &s = it->second;
        s.insert(kv.first);
      }
    }
    return res;
  }

  // Without path-compression
  equivalence_class_elems_t equiv_classes_elems() const {
    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t rep = find(kv.second);
      auto it = res.find(rep);
      if (it == res.end()) {
        element_set_t s;
        s.insert(kv.first);
        res.insert({rep, s});
      } else {
        element_set_t &s = it->second;
        s.insert(kv.first);
      }
    }
    return res;
  }  
  
  /* =============================================================
   *                     Abstract domain API
   * =============================================================
   */

  static union_find_domain_t top() {
    return union_find_domain_t(lattice_val::top);
  }

  static union_find_domain_t bottom() {
    return union_find_domain_t(lattice_val::bottom);
  }

  bool is_bottom() const { return (m_val == lattice_val::bottom); }

  bool is_top() const { return (m_val == lattice_val::top); }

  void set_to_top() {
    clear();
    m_val = lattice_val::top;
  }

  void set_to_bottom() {
    clear();
    m_val = lattice_val::bottom;
  }

  void set(const element_t &x, domain_t dom) {
    if (is_bottom()) { 
      return;
    }

    MergeSemantics op;
    if (op.is_bottom_absorbing() && dom.is_bottom()) {
      set_to_bottom();
      return;
    }

    if (!contains(x)) {
      // Create a singleton equivalence class
      make(x, dom);
    } else {
      // Modify the abstract state of the whole equivalence class
      element_t rep_x = find(x);
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      ec_x.get_domain() = dom;
    }
  }

  // Return null if !contains(x)
  domain_t *get(const element_t &x) {
    if (is_bottom()) {
      CRAB_ERROR("called union_find_domain::operator[] on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("called union_find_domain::operator[] on top");
    }
    if (!contains(x)) {
      return nullptr;
    }

    element_t rep_x = find(x);
    equivalence_class_t &ec_x = m_classes.at(rep_x);
    return &(ec_x.get_domain());
  }

  // Return null if !contains(x)
  const domain_t *get(const element_t &x) const {
    if (is_bottom()) {
      CRAB_ERROR("called union_find_domain::operator[] on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("called union_find_domain::operator[] on top");
    }
    if (!contains(x)) {
      return nullptr;
    }

    element_t rep_x = find(x);
    const equivalence_class_t &ec_x = m_classes.at(rep_x);
    return &(ec_x.get_domain());
  }  

  // Return true if *this is a refined partitioning of o
  //    forall x,y \in elems(o) :: o.find(x) != o.find(y) =>
  //                              this.find(x) != this.find(y)
  bool operator<=(const union_find_domain_t &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    } else {
      union_find_domain_t left(*this);
      union_find_domain_t right(o);
      equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
      equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();
      for (auto &kv : right_equiv_classes) {
        element_t right_repr = kv.first;
        element_set_t &right_elems = kv.second;
        for (const element_t &right_v : right_elems) {
          if (!left.contains(right_v)) {
            return false;
          }
          element_t left_repr = left.find(right_v);
          element_set_t &left_elems = left_equiv_classes[left_repr];

	  // !(left_elems <= right_elems)
	  if (!std::includes(right_elems.begin(), right_elems.end(),
			     left_elems.begin(), left_elems.end())) {	    
            return false;
          }
          const domain_t &left_dom =
              left.get_equiv_class(left_repr).get_domain();
          const domain_t &right_dom =
              right.get_equiv_class(right_repr).get_domain();
          if (!(left_dom <= right_dom)) {
            return false;
          }
        }
      }
      return true;
    }
  }

  union_find_domain_t operator|(const union_find_domain_t &o) const {
    return join_or_widening(o, true /*is_join*/);
  }

  union_find_domain_t operator&(const union_find_domain_t &o) const {
    return meet_or_narrowing(o, true /*is_meet*/);
  }

  union_find_domain_t operator||(const union_find_domain_t &o) const {
    return join_or_widening(o, false /*is_join*/);    
  }
  
  union_find_domain_t operator&&(const union_find_domain_t &o) const {
    return meet_or_narrowing(o, false /*is_meet*/);
  }

  // Add y into the equivalence class of x
  void add(const element_t &x, const element_t &y) {
    if (is_bottom() || is_top()) {
      return;
    }
    if (!contains(x)) {
      return;
    }
    if (contains(y)) {
      forget(y);
    }
    element_t rep_x = find(x);
    m_parents.insert({y, rep_x});
  }

  // The domain attached to the affected equivalence class is not
  // modified.  
  void forget(const element_t &v) {
    if (is_bottom() || is_top()) {
      return;
    }

    if (contains(v)) {
      if (m_parents.at(v) != v) {
        // v is not a representative
        element_t rep_v = find(v);
        for (auto &kv : m_parents) {
          if (kv.second == v) {
            m_parents.at(kv.first) = rep_v;
          }
        }
      } else {
        // v is the representative of the equivalence class
        boost::optional<element_t> new_rep;
        for (auto &kv : m_parents) {
          if (kv.first != v && kv.second == v) {
            if (!new_rep) {
              // choose the first one as representative
              new_rep = kv.first;
              m_classes.insert({*new_rep, m_classes.at(v)});
            }
            m_parents.at(kv.first) = *new_rep;
          }
        }
        m_classes.erase(v);
      }
      m_parents.erase(v);
    }
  }

  // The domains attached to the affected equivalence classes are not
  // modified
  void project(const std::vector<element_t> &elements) {
    if (is_bottom() || is_top()) {
      return;
    }

    // First, collect all elements to be forgotten
    std::vector<element_t> sorted_elements(elements);
    std::sort(sorted_elements.begin(), sorted_elements.end());
    std::vector<element_t> elements_to_forget;
    elements_to_forget.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      auto lower = std::lower_bound(sorted_elements.begin(),
                                    sorted_elements.end(), kv.first);
      if (lower == sorted_elements.end() || kv.first < *lower) { // not found
        elements_to_forget.push_back(kv.first);
      }
    }
    // Forget elements
    for (auto &v : elements_to_forget) {
      forget(v);
    }
  }

  // The domain attached to the affected equivalence class is not
  // modified
  void rename(const std::vector<element_t> &old_elements,
              const std::vector<element_t> &new_elements) {
    if (is_top() || is_bottom()) {
      return;
    }
    if (old_elements.size() != new_elements.size()) {
      CRAB_ERROR(
          "union_find_domain::rename with input vectors of different sizes");
    }

    auto rename_elem = [this](const element_t &x, const element_t &y) {
      if (contains(y)) {
        CRAB_ERROR("union_find_domain::rename assumes that ", y,
                   " does not exist");
      }
      if (contains(x)) {
        // tuples where x is the left operand
        std::vector<typename parents_map_t::value_type> left_x;
        for (auto &kv : m_parents) {
          if (kv.first == x) {
            left_x.push_back(kv);
          } else if (kv.second == x) {
            kv.second = y;
          }
        }
        for (auto &kv : left_x) {
          m_parents.erase(kv.first);
          m_parents.insert({y, (kv.second == x ? y : kv.second)});
        }
        auto it = m_classes.find(x);
        if (it != m_classes.end()) {
          equivalence_class_t equiv_class = it->second;
          m_classes.erase(it);
          m_classes.insert({y, std::move(equiv_class)});
        }
      }
    };

    for (unsigned i = 0, sz = old_elements.size(); i < sz; ++i) {
      rename_elem(old_elements[i], new_elements[i]);
    }
  }

  void write(crab_os &o) const {
    if (is_top()) {
      o << "{}";
    } else if (is_bottom()) {
      o << "_|_";
    } else {

      // Sort everything to print always the same thing
      auto sorted_equiv_classes = [](const equivalence_class_elems_t &unsorted_map) {
	     std::vector<std::pair<element_t, std::vector<element_t>>> sorted_eq_classes;
	     for (auto &kv: unsorted_map) {
	       std::vector<element_t> sorted_eq_class(kv.second.begin(), kv.second.end());
	       std::sort(sorted_eq_class.begin(), sorted_eq_class.end());
	       sorted_eq_classes.emplace_back(std::make_pair(kv.first, sorted_eq_class));
	     }
	     std::sort(sorted_eq_classes.begin(), sorted_eq_classes.end(),
		       [](const std::pair<element_t, std::vector<element_t>> &p1,
			  const std::pair<element_t, std::vector<element_t>> &p2) {
			 return p1.first < p2.first;
		       });
	     return sorted_eq_classes;
      };

      auto print_elems_vector = [&o](const std::vector<element_t>& elems) {
				 o << "{";
				 for (auto it = elems.begin(), et = elems.end(); it!=et;) {
				   o << *it;
				   ++it;
				   if (it != et) {
				     o << ",";
				   }
				 }
				 o << "}";
			       };
      
      union_find_domain_t tmp(*this);
      CRAB_LOG("union-find-print", tmp.print(o); o << "\n";);
      equivalence_class_elems_t ec_elems = tmp.equiv_classes_elems();
      auto sorted_ec_elems = sorted_equiv_classes(ec_elems);
      o << "{";
      for (auto it = sorted_ec_elems.begin(), et = sorted_ec_elems.end(); it != et;) {
	auto &p = *it;
        element_t &rep = tmp.find(p.first);
	print_elems_vector(p.second);
	o << "=>" << tmp.m_classes.at(rep).get_domain();
        ++it;
        if (it != et) {
          o << ",";
        }
      }
      o << "}";
    }
  }

  friend class crab::crab_os &operator<<(crab::crab_os &o,
                                         const union_find_domain_t &dom) {
    dom.write(o);
    return o;
  }
}; // end class union_find_domain

} // end namespace domains
} // end namespace crab
