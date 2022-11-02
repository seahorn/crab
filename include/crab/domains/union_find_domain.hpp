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
  std::shared_ptr<Domain> m_val;

  // Copy-on-write: call always this function before get_absval() if
  // m_val might be modified.
  void detach_absval() { m_val.reset(new Domain(*m_val));}
  
public:

  explicit equivalence_class(std::shared_ptr<Domain> val)
    : m_rank(0), m_val(val) {}

  std::size_t &get_rank() { return m_rank; }

  std::size_t get_rank() const { return m_rank; }

  std::shared_ptr<Domain> detach_and_get_absval() {
    if (m_val.use_count() > 1) {
      detach_absval();
    }
    return m_val;
  }
  
  std::shared_ptr<const Domain> get_absval() const {
    return m_val;
  }

  void set_absval(std::shared_ptr<Domain> val) { m_val = val;}
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
struct uf_merge_semantics {
  virtual ~uf_merge_semantics() {} 
  virtual bool is_bottom_absorbing() const = 0;
  virtual void apply(Domain &d1, const Domain &d2) const = 0;
};

template<class Domain>    
struct uf_union_semantics: public uf_merge_semantics<Domain> {
  virtual bool is_bottom_absorbing() const override {
    return false;
  }
  virtual void apply(Domain &d1, const Domain &d2) const override {
    d1 |= d2;
  }
};
  
template<class Domain>    
struct uf_intersection_semantics: public uf_merge_semantics<Domain> {
  virtual bool is_bottom_absorbing() const override {
    return true;
  }  
  virtual void apply(Domain &d1, const Domain &d2) const override {
    d1 &= d2;
  }
};

  
// If Element and Domain are completely disjoint then this operation
// is non-op. However, Element can be a variable and Domain can be,
// for instance, an abstract domain that maps variables to abstract
// states. In that case, if we forget a variable from an equivalence
// class then we might want also to forget the variable from the
// attached abstract domain. So this is something that should be
// decided by the client of the union-find (see template parameters
// ForgetElementInDomain and RenameElementInDomain).
template<class Element, class Domain>  
struct uf_forget_element_in_domain {
  void operator()(Domain &dom, const Element&e) const {}
};
template<class Element, class Domain>  
struct uf_rename_element_in_domain {
  void operator()(Domain &dom, const Element&e, const Element &new_e) const {}
};
  
template<class Element, class Domain,
	 class MergeSemantics = uf_union_semantics<Domain>,
	 class ForgetElementInDomain = uf_forget_element_in_domain<Element, Domain>,
	 class RenameElementInDomain = uf_rename_element_in_domain<Element, Domain>>
class union_find_domain {
public:
  using element_t = Element;
  using domain_t = Domain;
  enum class lattice_val { bottom, top, neither_top_nor_bot };    
private:
  using union_find_domain_t = union_find_domain<Element,
						Domain,
						MergeSemantics,
						ForgetElementInDomain,
						RenameElementInDomain>;
  using parents_map_t = std::unordered_map<element_t, element_t>;
public:
  using element_set_t = std::vector<element_t>; 
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
	      std::shared_ptr<domain_t> absval = nullptr) {
    crab::CrabStats::count("union_find.count.merge_elements");
    crab::ScopedCrabStats __st__("union_find.merge_elements");
    
    auto it = elems.begin();
    auto et = elems.end();
    if (it == et) {
      return boost::none;
    }
    element_t v = *it;
    if (absval && !contains(v)) {
      make(v, absval);
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
      if (absval && !contains(*it)) {
        make(*it, absval);
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
      const equivalence_class_t &ec = it->second;
      o << rep << " -> " << *(ec.get_absval());
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "})";
  }

  void get_all_members(std::vector<element_t> &out) const {
    out.reserve(m_parents.size());
    for (auto &kv: m_parents) {
      out.push_back(kv.first);
    }
  }

  // 1. Ensure that left and right have the same equivalence classes.
  // 2. Apply component-wise the join or the widening.
  union_find_domain_t join_or_widening(const union_find_domain_t &o, bool is_join) const {
    crab::CrabStats::count("union_find.count.join_or_widening");
    crab::ScopedCrabStats __st__("union_find.join_or_widening");
    
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() || o.is_top()) {
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

      // Merge equivalence classes from the left to the right
      for (auto &kv : left_equiv_classes) {
        boost::optional<element_t> right_repr =
	  right.merge_elems(kv.second);
        if (!right_repr) {
	  // this shouldn't happen
	  CRAB_ERROR("unexpected situation in join_or_widening 1");
	}
      }
      left_equiv_classes.clear();
      
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
	auto &left_ec = left.get_equiv_class(*left_repr);
	std::shared_ptr<domain_t> left_absval = left_ec.detach_and_get_absval();
	std::shared_ptr<const domain_t> right_absval =
	  right.get_equiv_class(kv.first).get_absval();
	if (is_join) {
	  *left_absval |= *right_absval;
	} else {
	  *left_absval = *left_absval || *right_absval;
	}
      }
      CRAB_LOG("union-find", crab::outs() << "Res=" << left << "\n";);
      return left;
    }
  }

  // 1. Ensure that left and right have the same equivalence classes.
  // 2. Apply component-wise the meet or the narrowing
  union_find_domain_t meet_or_narrowing(const union_find_domain_t &o, bool is_meet) const {
    crab::CrabStats::count("union_find.count.meet_or_narrowing");
    crab::ScopedCrabStats __st__("union_find.meet_or_narrowing");
    
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
      // The merging on the right is needed so that right_absval is updated.
      for (auto &kv : left_equiv_classes) {
	std::shared_ptr<domain_t> left_absval =
	  left.get_equiv_class(kv.first).detach_and_get_absval();
        // add elements on the right if they do not exist
        boost::optional<element_t> right_repr =
	  right.merge_elems(kv.second, left_absval);
        if (!right_repr) {
	  // this shouldn't happen
	  CRAB_ERROR("unexpected situation in meet_or_narrowing 1");
	}
      }
      left_equiv_classes.clear();
      
      // Merge equivalence classes from the right to the left while
      // applying meet/narrowing
      for (auto &kv : right_equiv_classes) {
	std::shared_ptr<domain_t> right_absval =
            right.get_equiv_class(kv.first).detach_and_get_absval();
        boost::optional<element_t> left_repr =
	  left.merge_elems(kv.second, right_absval);
        if (!left_repr) {
	   // this shouldn't happen
	  CRAB_ERROR("unexpected situation in meet_or_narrowing 2");
	}
	auto &left_ec = left.get_equiv_class(*left_repr);
	std::shared_ptr<domain_t> left_absval = left_ec.detach_and_get_absval();
	if (is_meet) {
	  // TODO: operator &= might not be defined if domain_t is not
	  // derived from abstract_domain
	  (*left_absval) = (*left_absval) & (*right_absval);	  	  
	} else {
	  (*left_absval) = (*left_absval) && (*right_absval);	  
	}
        if (left_absval->is_bottom()) {
          left.set_to_bottom();
          return left;
        }
      }
      return left;
    }
  }

  struct get_absval_from_ec:
    public std::unary_function<typename classes_map_t::value_type, Domain> {

    std::shared_ptr<const Domain> operator()(const typename classes_map_t::value_type &kv) const {
      return kv.second.get_absval();
    }
    
    std::shared_ptr<Domain> operator()(typename classes_map_t::value_type &kv) const {
      return kv.second.detach_and_get_absval();
    }     
  };  
public:
 
  using const_domain_iterator =
    boost::transform_iterator<get_absval_from_ec, typename classes_map_t::const_iterator>;
  using domain_iterator =
    boost::transform_iterator<get_absval_from_ec, typename classes_map_t::iterator>;
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
  void make(const element_t &v, std::shared_ptr<Domain> val) {
    crab::CrabStats::count("union_find.count.make");
    crab::ScopedCrabStats __st__("union_find.make");

    if (crab::CrabSanityCheckFlag) {
      if (is_bottom()) {
	CRAB_ERROR("calling union_find_domain::make on bottom");
      }
      if (is_top()) {
	CRAB_ERROR("calling union_find_domain::make on top");
      }
      if (contains(v)) {
	CRAB_ERROR("element already exists when called union_find_domain::make");
      }
    }

    MergeSemantics op;
    if (op.is_bottom_absorbing() && val->is_bottom()) {
      set_to_bottom();
    } else {
      m_parents.insert({v, v});
      m_classes.insert({v, equivalence_class_t(val)});
    }
  }
  
  // Pre-condition: contains(v)
  // NOTE: it is not "const" because it does path-compression
  element_t &find(const element_t &v) {
    crab::CrabStats::count("union_find.count.find");
    crab::ScopedCrabStats __st__("union_find.find");
    
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
    crab::CrabStats::count("union_find.count.find");
    crab::ScopedCrabStats __st__("union_find.find");
    
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
    crab::CrabStats::count("union_find.count.union");
    crab::ScopedCrabStats __st__("union_find.union");
    
    element_t rep_x = find(x);
    element_t rep_y = find(y);
    MergeSemantics merge_op;
    if (rep_x != rep_y) {
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      equivalence_class_t &ec_y = m_classes.at(rep_y);
      if (ec_x.get_rank() > ec_y.get_rank()) {
	std::shared_ptr<domain_t> dom_x = ec_x.detach_and_get_absval();
	std::shared_ptr<const domain_t> dom_y = ec_y.get_absval();
        merge_op.apply(*dom_x, *dom_y);
	if (merge_op.is_bottom_absorbing() && dom_x->is_bottom()) {
	  set_to_bottom();
	  return boost::none;
	}								   
        m_parents.at(rep_y) = rep_x;
        m_classes.erase(rep_y);
        return rep_x;
      } else {
	std::shared_ptr<const domain_t> dom_x = ec_x.get_absval();
	std::shared_ptr<domain_t> dom_y = ec_y.detach_and_get_absval();
        merge_op.apply(*dom_y, *dom_x);
	if (merge_op.is_bottom_absorbing() && dom_y->is_bottom()) {
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
					  get_absval_from_ec());
  }
  
  const_domain_iterator end_domains() const {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: end_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: end_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.end(),
					  get_absval_from_ec());
  }

  domain_iterator begin_domains() {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: begin_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: begin_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.begin(),
					  get_absval_from_ec());
  }

  domain_iterator end_domains() {
    if (is_bottom()) {
      CRAB_ERROR("union-find domain: end_domains() on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("union-find domain: end_domains() on top");
    }     
    return boost::make_transform_iterator(m_classes.end(),
					  get_absval_from_ec());
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
    crab::CrabStats::count("union_find.count.equiv_classes_elems");
    crab::ScopedCrabStats __st__("union_find.equiv_classes_elems");
    
    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t rep = find(kv.second);
      element_set_t &s = res[rep];
      auto it = std::upper_bound(s.begin(), s.end(), kv.first);
      s.insert(it, kv.first);
    }
    return res;    
  }

  // Without path-compression
  equivalence_class_elems_t equiv_classes_elems() const {
    crab::CrabStats::count("union_find.count.equiv_classes_elems");
    crab::ScopedCrabStats __st__("union_find.equiv_classes_elems");
    
    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t rep = find(kv.second);
      element_set_t &s = res[rep];
      auto it = std::upper_bound(s.begin(), s.end(), kv.first);
      s.insert(it, kv.first);
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

  void set(const element_t &x, domain_t absval) {
    if (is_bottom()) { 
      return;
    }

    MergeSemantics op;
    if (op.is_bottom_absorbing() && absval.is_bottom()) {
      set_to_bottom();
      return;
    }

    std::shared_ptr<domain_t> absval_ptr = std::make_shared<domain_t>(absval);
    
    if (!contains(x)) {
      // Create a singleton equivalence class
      make(x, absval_ptr);
    } else {
      // Modify the abstract state of the whole equivalence class
      element_t rep_x = find(x);
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      ec_x.set_absval(absval_ptr);
    }
  }

  // Return null if !contains(x)
  std::shared_ptr<domain_t> get(const element_t &x) {
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
    return ec_x.detach_and_get_absval();
  }

  // Return null if !contains(x)
  std::shared_ptr<const domain_t> get(const element_t &x) const {
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
    return ec_x.get_absval();
  }  

  // Return true if *this is a refined partitioning of o
  //    forall x,y \in elems(o) :: o.find(x) != o.find(y) =>
  //                              this.find(x) != this.find(y)
  bool operator<=(const union_find_domain_t &o) const {
    crab::CrabStats::count("union_find.count.less_or_equal");
    crab::ScopedCrabStats __st__("union_find.less_or_equal");
    
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    } else {
      const union_find_domain_t &left = *this;
      const union_find_domain_t &right = o;
      equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
      equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();
      for (auto &kv : right_equiv_classes) {
        element_t right_repr = kv.first;
        const element_set_t &right_elems = kv.second;
        for (const element_t &right_v : right_elems) {
          if (!left.contains(right_v)) {
            return false;
          }
          element_t left_repr = left.find(right_v);
          const element_set_t &left_elems = left_equiv_classes[left_repr];

	  if (!std::includes(right_elems.begin(), right_elems.end(),
			     left_elems.begin(), left_elems.end())) {	    
            return false;
          }
	  std::shared_ptr<const domain_t> left_absval =
              left.get_equiv_class(left_repr).get_absval();
	  std::shared_ptr<const domain_t> right_absval =
              right.get_equiv_class(right_repr).get_absval();
          if (!(*left_absval <= *right_absval)) {
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

  // Remove v from its equivalence class. Moreover, it applies
  // ForgetElementInDomain on the domain attached to the equivalence
  // class.
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
	equivalence_class_t &ec = m_classes.at(rep_v);
	ForgetElementInDomain{}(*(ec.detach_and_get_absval()), v);
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
	if (new_rep) {
	  equivalence_class_t &ec = m_classes.at(*new_rep);
	  ForgetElementInDomain{}(*(ec.detach_and_get_absval()), v);
	} else {
	  // do nothing: v is the only element in its equivalence
	  // class
	}
      }
      m_parents.erase(v);
    }
  }

  
  // Remove from all equivalence classes any variable that is not in
  // elements. 
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

  // Rename old_elements in the affected equivalence classes with
  // new_elements. Moreover, it applies RenameElementInDomain on the
  // domain attached to the equivalence class.
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
	// remove the key-value entry where key==x
	// rename key-value entries   where value==x
	boost::optional<element_t> parent_x;
	for (auto it=m_parents.begin(), et=m_parents.end();it!=et;) {
	  if ((*it).first == x) {
	    parent_x = (*it).second;
	    it = m_parents.erase(it);
	  } else if ((*it).second == x) {
	    (*it).second = y;
	    ++it;
	  } else {
	    ++it;
	  }
	}

	if (parent_x) {
	  m_parents.insert({y, (*parent_x == x ? y: *parent_x)});
	  auto it = m_classes.find(x);
	  if (it != m_classes.end()) {
	    equivalence_class_t ec = it->second;
	    RenameElementInDomain{}(*(ec.detach_and_get_absval()), x, y); 
	    m_classes.erase(it);
	    m_classes.insert({y, std::move(ec)});
	  }
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
	o << "=>" << *(tmp.m_classes.at(rep).get_absval());
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
