#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/constant.hpp>
#include <crab/domains/symbolic_variable_eq_domain.hpp>
#include <crab/domains/term/term_operators.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

namespace symbolic_variable_equiality_domain_impl {
class symbolic_var {
  using this_domain_t = symbolic_var;

public:
  using var_id = uint32_t;

  var_id var;

  symbolic_var(var_id _var) : var(_var) {}
  symbolic_var(const this_domain_t &o) = default;
  symbolic_var(this_domain_t &&o) = default;
  this_domain_t &operator=(const this_domain_t &o) = default;
  this_domain_t &operator=(this_domain_t &&o) = default;

  void write(crab_os &o) const { o << "#var" << var; }

  bool operator==(this_domain_t o) const { return var == o.var; }

  operator var_id() const { return var; }

  friend class crab::crab_os &operator<<(crab::crab_os &o,
                                         const this_domain_t &dom) {
    dom.write(o);
    return o;
  }
};

struct symb_var_hash {
  std::size_t operator()(const symbolic_var &k) const {
    return std::hash<uint32_t>()(k.var);
  }
};

template <class InputIt1, class InputIt2>
bool set_is_absolute_complement(InputIt1 first1, InputIt1 last1,
                                InputIt2 first2, InputIt2 last2) {
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2)
      ++first1;
    else {
      if (!(*first2 < *first1))
        return false; // *first1 and *first2 are equivalent.
      ++first2;
    }
  }
  return true;
}

template <class SYMBVAR> SYMBVAR make_fresh_var_symbol() {
  return SYMBVAR(term::term_op_val_generator_t::get_next_val());
}

template <typename Key, typename Value>
void print_unordered_map(crab::crab_os &o,
                         std::unordered_map<Key, Value> const &m) {
  o << "{";
  for (auto it = m.begin(), et = m.end(); it != et;) {
    o << it->first << " => " << it->second;
    ++it;
    if (it != et) {
      o << ", ";
    }
  }
  o << "}";
}
} // namespace symbolic_variable_equiality_domain_impl

/// @brief An abstract domain represents equalities used in analyses such as
/// allocation-site abstraction or domain reduction. In short, giving two
/// elements are known to hold equal values in concrete semantics if
/// the domain captures that the elements hold equal symbolic variable.
/// @tparam BaseDomain
template <class BaseDomain> class symbolic_variable_equiality_domain {
public:
  // typedefs for BaseDomain
  using linear_constraint_system_t =
      typename BaseDomain::linear_constraint_system_t;
  using linear_constraint_t = typename BaseDomain::linear_constraint_t;
  using linear_expression_t = typename BaseDomain::linear_expression_t;
  using reference_constraint_t = typename BaseDomain::reference_constraint_t;
  using variable_or_constant_t = typename BaseDomain::variable_or_constant_t;
  using variable_t = typename BaseDomain::variable_t;
  using variable_vector_t = typename BaseDomain::variable_vector_t;
  using varname_t = typename BaseDomain::varname_t;

  // typedefs for equality domain
  using element_t = variable_t;
  using element_set_t = std::vector<element_t>;
  using domain_t = class symbolic_variable_equiality_domain_impl::symbolic_var;
  using var_id = typename domain_t::var_id;
  using parents_map_t = std::unordered_map<element_t, element_t>;
  using equivalence_class_t = equivalence_class<domain_t>;
  using equivalence_class_elems_t =
      std::unordered_map<element_t, element_set_t>;

private:
  using this_domain_t = symbolic_variable_equiality_domain<BaseDomain>;
  using classes_map_t = std::unordered_map<element_t, equivalence_class_t>;
  enum class lattice_val { bottom, top, neither_top_nor_bot };

  // a map that stores a reference to its immediate representative
  // The map is a many to one hash map
  parents_map_t m_parents;
  // a map that keeps a representative to its corresponding symbolic variable
  // each representative has its equivalence class
  classes_map_t m_classes;

  lattice_val m_val;
  // For example,
  // The disjoin set for a state of symbolic_variable_equiality_domain<int>:
  // { 1, 2, 5 } |-> #var1
  // { 4, 3 } |-> #var2
  // { 6 } |-> #var3
  // can be represented as:
  //                #var1        #var2          #var3
  //  m_classes       ▲            ▲              ▲
  // ─ ─ ─ ─ ─ ─ ─ ─ ─│─ ─ ─ ─ ─ ─ ┼ ─ ─ ─ ─ ─ ─ ─│─ ─ ─ ─ ─ ─ ─ ─
  //  m_parents      ┌─┐          ┌┴┐            ┌┴┐
  //             ┌──▶│1│◀─┐       │4│◀─┐         │6│
  //             │   └─┘  │       └─┘  │         └─┘
  //            ┌─┐      ┌─┐          ┌─┐
  //            │2│      │5│          │3│
  //            └─┘      └─┘          └─┘
  // Note that, { 6 } |-> #var3 will be removed after normalization since
  // it does not caputre any relations between elements.

  /// @brief a helper method to empty the disjoin set
  void clear() {
    m_parents.clear();
    m_classes.clear();
  }

  /// @brief Build a map from representative to a ordered set with all the
  /// elements in the equivalence class.
  /// @return a map computed as brief described.
  equivalence_class_elems_t equiv_classes_elems() {
    crab::CrabStats::count("union_find.count.equiv_classes_elems");
    crab::ScopedCrabStats __st__("union_find.equiv_classes_elems");

    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t &rep = find(kv.second);
      element_set_t &s = res[rep];
      auto it = std::upper_bound(s.begin(), s.end(), kv.first);
      s.insert(it, kv.first);
    }
    return res;
  }

  /// @brief Build a map from representative to a ordered set with all the
  /// elements in the equivalence class **without path-compression**
  /// @return a map computed as brief described.
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

  element_set_t get_all_members() const {
    element_set_t out;
    out.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      out.push_back(kv.first);
    }
    return out;
  }

  element_set_t get_all_members_from_an_equiv_class(const element_t &e) const {
    element_t e_rep = find(e);
    element_set_t out;
    out.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      element_t rep = find(kv.second);
      if (rep == e_rep) {
        auto it = std::upper_bound(out.begin(), out.end(), kv.first);
        out.insert(it, kv.first);
      }
    }
    return out;
  }

  /// @brief a helper method for domain operation join
  /// @param left one abstract state
  /// @param right one abstract state
  /// @return a new abstract state saved the joined result
  /// TODO: This operation is not an optimal solution. This computation does
  /// not follow the feature of union-find data structure.
  /// The solution is computing the intersection of two equivalence classes
  /// to find least common elements. It only computes a set of
  /// equivalence classes where each class is subset of two classes from the
  /// given two abstract states. That is, the result of the join only remains
  /// equalities appeared both in left and in right states.
  /// Formally, forall cls_x \in left. forall cls_y \in right ::
  ///           cls_x ^ cls_y
  this_domain_t join(const this_domain_t &left,
                     const this_domain_t &right) const {
    equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();
    equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
    this_domain_t res;
    var_id current_id = 0;
    for (auto &left_kv : left_equiv_classes) {
      const element_set_t &left_set = left_kv.second;
      for (auto &right_kv : right_equiv_classes) {
        const element_set_t &right_set = right_kv.second;
        element_set_t out_set;
        out_set.reserve(std::min(left_set.size(), right_set.size()));
        CRAB_LOG("symb-var-join", crab::outs() << "intersect left: ";
                 print_elems_vector(crab::outs(), left_set);
                 crab::outs() << ", right: ";
                 print_elems_vector(crab::outs(), right_set);
                 crab::outs() << "\n";);
        // set_intersection is a linear time operation
        std::set_intersection(left_set.begin(), left_set.end(),
                              right_set.begin(), right_set.end(),
                              std::back_inserter(out_set));
        CRAB_LOG("symb-var-join", crab::outs() << "output: ";
                 print_elems_vector(crab::outs(), out_set);
                 crab::outs() << "\n";);
        if (out_set.size() > 1) {
          // select one element as representative
          const element_t &rep = out_set[0];
          std::shared_ptr<domain_t> absval_ptr =
              std::make_shared<domain_t>(domain_t(current_id));
          equivalence_class_t new_val = equivalence_class_t(absval_ptr);
          res.m_classes.insert({rep, new_val});
          for (auto &v : out_set) {
            res.m_parents.insert({v, rep});
          }
          current_id++;
        }
      }
    }
    res.normalize();
    return res;
  }

  /// @brief a helper method for domain operation meet
  /// @param left one abstract state
  /// @param right one abstract state
  /// @return a new abstract state saved the meet result
  /// TODO: This operation is not an optimal solution. This computation does
  /// not follow the feature of union-find data structure.
  /// The solution is merging two equivalence classes if two equivalence classes
  /// shared common elements. As the result, it computed a set of a set of
  /// equivalence classes where each class is a super set of two or more classes
  /// from the given two abstract states. That is, the result of the meet will
  /// keep equalities appeared either in left or in right states.
  /// Formally, forall cls_x \in left. forall cls_y \in right ::
  ///           cls_x ^ cls_y != \empty => cls_x union cls_y
  this_domain_t meet(const this_domain_t &left,
                     const this_domain_t &right) const {
    equivalence_class_elems_t left_equiv_classes = left.equiv_classes_elems();
    equivalence_class_elems_t right_equiv_classes = right.equiv_classes_elems();
    this_domain_t res;
    var_id current_id = 0;
    auto meet_helper = [&current_id,
                        &res](const equivalence_class_elems_t &first_classes,
                              const equivalence_class_elems_t &second_classes) {
      for (auto &left_kv : first_classes) {
        const element_set_t &left_set = left_kv.second;
        element_set_t out_set = left_set;
        out_set.reserve(left_set.size());
        for (auto &right_kv : second_classes) {
          const element_set_t &right_set = right_kv.second;
          // we merge two equivalence classes if they shared common elements
          bool is_absolute_complement =
              symbolic_variable_equiality_domain_impl::
                  set_is_absolute_complement(left_set.begin(), left_set.end(),
                                             right_set.begin(),
                                             right_set.end());
          if (!is_absolute_complement) {
            std::set_union(left_set.begin(), left_set.end(), right_set.begin(),
                           right_set.end(), std::back_inserter(out_set));
          }
        }
        if (out_set.size() > 1) {
          // select one element as representative
          const element_t &rep = out_set[0];
          std::shared_ptr<domain_t> absval_ptr =
              std::make_shared<domain_t>(domain_t(current_id));
          equivalence_class_t new_val = equivalence_class_t(absval_ptr);
          res.m_classes.insert({rep, new_val});
          for (auto &v : out_set) {
            res.m_parents.insert({v, rep});
          }
          current_id++;
        }
      }
    };
    meet_helper(left_equiv_classes, right_equiv_classes);
    meet_helper(right_equiv_classes, left_equiv_classes);
    res.normalize();
    return res;
  }

  /// @brief create an equivalent class
  /// @param v the representative element for the new class
  /// @param val the domain value
  void make_set(const element_t &v, std::shared_ptr<domain_t> val) {
    if (is_bottom()) {
      CRAB_ERROR(domain_name(), "::", __func__, " make on bottom");
    }
    if (is_top()) {
      set_neither_top_or_bottom();
    }
    if (contains(v)) {
      CRAB_ERROR(domain_name(), "::", __func__, " the new element ", v,
                 " is already consisted in ", *this);
    }

    m_parents.insert({v, v});
    m_classes.insert({v, equivalence_class_t(val)});
  }

  void print_elems_vector(crab::crab_os &o,
                          const std::vector<element_t> &elems) const {
    o << "[";
    for (auto it = elems.begin(), et = elems.end(); it != et;) {
      o << *it;
      ++it;
      if (it != et) {
        o << ",";
      }
    }
    o << "]";
  }

  void print_equiv_classes(crab_os &o,
                           const equivalence_class_elems_t &equiv_classes,
                           bool verbose = false) const {
    o << "{";
    for (auto it = equiv_classes.begin(), et = equiv_classes.end(); it != et;) {
      print_elems_vector(o, it->second);
      if (!verbose) {
        o << "=>" << *(m_classes.at(it->first).get_absval());
      }
      ++it;
      if (it != et) {
        o << ",";
      }
    }
    o << "}";
  }

  void print_classes_vals(crab_os &o) const {
    o << "{";
    for (auto it = m_classes.begin(), et = m_classes.end(); it != et;) {
      o << it->first << "=>" << *(it->second.get_absval());
      ++it;
      if (it != et) {
        o << ",";
      }
    }
    o << "}";
  }

public:
  /**------------------ Begin union find APIs ------------------**/
  /// @brief Check whether domain constains element v in some class
  /// @param v an element
  /// @return true if the elment exists; otherwise, false.
  bool contains(const element_t &v) const {
    return m_parents.find(v) != m_parents.end();
  }

  /// @brief Check whether two elements in the same class even if they may not
  /// exist
  /// @param x an element
  /// @param y an element differs from x
  /// @return true if they are in the same class; otherwise, return false
  bool equals(const element_t &x, const element_t &y) const {
    if (!contains(x) || !contains(y)) {
      return false;
    }
    const element_t &rep_x = find(x);
    const element_t &rep_y = find(y);
    return rep_x == rep_y;
  }

  /// @brief find the representative with path-compression
  /// @param v an element in some set
  /// @attention v must constain in current domain; otherwise, an error returns
  /// @return returns the representative of the set that contains the element v
  /// NOTE: it is not "const" because it does path-compression
  element_t &find(const element_t &v) {
    crab::CrabStats::count("union_find.count.find");
    crab::ScopedCrabStats __st__("union_find.find");

    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      CRAB_ERROR(domain_name(), "::", __func__, " on a non-existing elem ", v,
                 " in ", *this);
    }
    element_t &parent = it->second;
    if (v == parent) {
      return parent;
    } else {
      return parent = find(parent);
    }
  }

  /// @brief find the representative without path-compression
  /// @param v an element in some set
  /// @attention v must constain in current domain; otherwise, an error returns
  /// @return returns the representative of the set that contains the element v
  element_t find(const element_t &v) const {
    crab::CrabStats::count("union_find.count.find");
    crab::ScopedCrabStats __st__("union_find.find");

    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      CRAB_ERROR(domain_name(), "::", __func__, " on a non-existing elem ", v,
                 " in ", *this);
    }
    const element_t &parent = it->second;
    if (v == parent) {
      return parent;
    } else {
      return find(parent);
    }
  }

  /// @brief union two sets
  /// @param x one set in domain
  /// @param y one set in domain
  void union_sets(const element_t &x, const element_t &y) {
    const element_t &rep_x = find(x);
    const element_t &rep_y = find(y);
    if (rep_x == rep_y) {
      return; // x and y are already in the same set
    }
    equivalence_class_t &ec_x = m_classes.at(rep_x);
    equivalence_class_t &ec_y = m_classes.at(rep_y);
    // Check the rank for two sets.
    // The rank is the depth of the tree
    if (ec_x.get_rank() > ec_y.get_rank()) {
      // merge elements in set rep_y to set rep_x
      std::shared_ptr<domain_t> dom_x = ec_x.detach_and_get_absval();
      std::shared_ptr<const domain_t> dom_y = ec_y.get_absval();
      m_parents.at(rep_y) = rep_x;
      m_classes.erase(rep_y);
    } else {
      // merge elements in set rep_x to set rep_y
      std::shared_ptr<const domain_t> dom_x = ec_x.get_absval();
      std::shared_ptr<domain_t> dom_y = ec_y.detach_and_get_absval();
      m_parents.at(rep_x) = rep_y;
      if (ec_x.get_rank() == ec_y.get_rank()) {
        ec_y.get_rank()++;
      }
      m_classes.erase(rep_x);
    }
  }

  /// @brief set a domain value to x's class
  /// @param x an element
  /// @param absval a domain value
  /// @note  if x does not exist, create an class with element x and absval
  ///        if x \in some cls, update current domain value by absval
  void set(const element_t &x, domain_t absval) {
    if (is_bottom()) {
      return;
    }

    std::shared_ptr<domain_t> absval_ptr = std::make_shared<domain_t>(absval);

    if (!contains(x)) {
      make_set(x, absval_ptr);
    } else {
      // Modify the abstract state of the whole equivalence class
      element_t rep_x = find(x);
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      ec_x.set_absval(absval_ptr);
    }
  }

  /// @brief get the domain value stored in current equivalent class
  /// @param v an element in some set
  /// @return null if !contains(x), otherwise, returns a shared pointer to value
  std::shared_ptr<domain_t> get(const element_t &x) {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top() || !contains(x)) {
      return nullptr;
    }

    element_t rep_x = find(x);
    equivalence_class_t &ec_x = m_classes.at(rep_x);
    return ec_x.detach_and_get_absval();
  }

  /// @brief get the equivalent class by giving a element
  /// @param x an element
  /// @return if the element exists, return corresponding element class
  boost::optional<element_set_t> get_variables(const element_t &x) const {
    if (contains(x)) {
      return get_all_members_from_an_equiv_class(x);
    } else {
      return boost::none;
    }
  }

  /// @brief get the equivalent class by giving a domain value
  /// @param dom a domain value such as #var1
  /// @return if the value exists, return corresponding element class
  boost::optional<element_set_t> get_variables(const domain_t &dom) const {
    for (auto &kv : m_classes) {
      const element_t &key = kv.first;
      const equivalence_class_t &ec = kv.second;
      auto dom_ptr = ec.get_absval();
      if (dom_ptr && dom == *dom_ptr) {
        return get_all_members_from_an_equiv_class(key);
      }
    }
    return boost::none;
  }

  /// @brief get the domain value stored in current equivalent class
  /// @param v an element in some set
  /// @return null if !contains(x), otherwise, returns a shared pointer to value
  std::shared_ptr<const domain_t> get(const element_t &x) const {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top() || !contains(x)) {
      return nullptr;
    }

    element_t rep_x = find(x);
    const equivalence_class_t &ec_x = m_classes.at(rep_x);
    return ec_x.get_absval();
  }

  /// @brief Add y into the equivalence class of x
  /// @param x an element may exists
  /// @param y an element may exists
  /// @note  if x does not exist, nothing changed.
  ///        if y \in some cls, forget it from cls, add y into x's class
  void add(const element_t &x, const element_t &y) {
    if (is_bottom() || is_top() || x == y) {
      return;
    }
    if (!contains(x)) {
      set(x, symbolic_variable_equiality_domain_impl::make_fresh_var_symbol<
                 symbolic_variable_equiality_domain_impl::symbolic_var>());
    }
    if (contains(y)) {
      *this -= y;
    }
    element_t rep_x = find(x);
    m_parents.insert({y, rep_x});
  }
  /**------------------ End union find APIs ------------------**/

  /**------------------ Begin domain APIs ------------------**/
  // empty union-find
  symbolic_variable_equiality_domain(
      lattice_val val = lattice_val::neither_top_nor_bot)
      : m_val(val) {}
  symbolic_variable_equiality_domain(const this_domain_t &o)
      : m_parents(o.m_parents), m_classes(o.m_classes), m_val(o.m_val) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  symbolic_variable_equiality_domain(this_domain_t &&o)
      : m_parents(std::move(o.m_parents)), m_classes(std::move(o.m_classes)),
        m_val(o.m_val) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  this_domain_t &operator=(const this_domain_t &o) {
    if (this != &o) {
      m_parents = o.m_parents;
      m_classes = o.m_classes;
      m_val = o.m_val;
    }
    return *this;
  }
  this_domain_t &operator=(this_domain_t &&o) {
    if (this != &o) {
      m_parents = std::move(o.m_parents);
      m_classes = std::move(o.m_classes);
      m_val = o.m_val;
    }
    return *this;
  }

  bool is_bottom() const { return m_val == lattice_val::bottom; }

  bool is_top() const { return m_val == lattice_val::top; }

  void set_to_top() {
    clear();
    m_val = lattice_val::top;
  }

  void set_to_bottom() {
    clear();
    m_val = lattice_val::bottom;
  }

  void set_neither_top_or_bottom() { m_val = lattice_val::neither_top_nor_bot; }

  bool operator<=(const this_domain_t &e) const {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    if (is_bottom() || e.is_top()) {
      // _|_ <= e || this <= top
      return true;
    } else if (is_top() || e.is_bottom()) {
      // top >= e || this >= _|_
      return false;
    }
    // Return true if *this is a refined partitioning of o
    //    forall cls_x \in this. exists cls_y \in e ::
    //      cls_y <= cls_x
    equivalence_class_elems_t left_equiv_classes =
        (*this).equiv_classes_elems();
    equivalence_class_elems_t right_equiv_classes = e.equiv_classes_elems();
    bool res = true;
    for (auto &right_kv : right_equiv_classes) {
      const element_set_t &right_set = right_kv.second;
      CRAB_LOG("symb-var-equiv-classes",
               print_elems_vector(crab::outs(), right_set););
      assert(right_set.size() > 0);
      bool tmp = false;
      for (auto &left_kv : left_equiv_classes) {
        const element_set_t &left_set = left_kv.second;
        CRAB_LOG("symb-var-equiv-classes",
                 print_elems_vector(crab::outs(), left_set););
        // subset relations implies the dual of domain operation \sqsubseteq
        tmp |= std::includes(left_set.begin(), left_set.end(),
                             right_set.begin(), right_set.end());
      }
      res &= tmp;
    }
    return res;
  }

  this_domain_t operator|(const this_domain_t &e) const {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    if (is_bottom()) {
      return e;
    } else if (e.is_bottom()) {
      return *this;
    } else if (is_top() || e.is_top()) {
      this_domain_t res;
      return res;
    } else {
      return join(*this, e);
    }
  }

  void operator|=(const this_domain_t &e) {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    if (is_bottom()) {
      *this = e;
    } else if (e.is_bottom()) {
      return;
    } else if (is_top() || e.is_top()) {
      set_to_top();
    } else {
      // TODO: improve performance
      *this = join(*this, e);
    }
  }

  this_domain_t operator&(const this_domain_t &e) const {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    if (is_bottom() || e.is_top()) {
      return *this;
    } else if (e.is_bottom() || is_top()) {
      return e;
    } else {
      return meet(*this, e);
    }
  }

  this_domain_t operator||(const this_domain_t &e) const {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return *this;
  }

  this_domain_t operator&&(const this_domain_t &e) const {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return *this;
  }

  /// @brief expand what x equals to to y. This is equivalent to add
  /// @param x the orginal variable in some class
  /// @param y
  void expand(const element_t &x, const element_t &y) {
    if (is_bottom() || is_top() || !contains(x)) {
      return;
    }
    add(x, y);
  }

  /// @brief assign equality, this is equivalent to add y to x's equiv cls
  /// @param x an element
  /// @param y an element
  void assign(const element_t &x, const element_t &y) {
    if (!is_bottom()) {
      add(x, y);
    }
  }

  /// @brief assign bool lhs := rhs, not used
  /// @param lhs an element
  /// @param rhs an element
  /// @param is_not_rhs if true the assignment is lhs := not(rhs)
  void assign_bool_var(const element_t &lhs, const element_t &rhs,
                       bool is_not_rhs) {
    // WARN: this should not happen
    return;
  }

  /// @brief remove an element from a class
  /// @param v an element
  /// @note  if v \in cls \land v is representative, select a new representative
  ///        before removing v.
  void operator-=(const element_t &v) {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    if (is_bottom() || is_top()) {
      return;
    }

    if (contains(v)) {
      if (m_parents.at(v) != v) {
        // v is not a representative
        element_t rep_v = find(v);
        for (auto &kv : m_parents) {
          // set all elements referring v to v's parent
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
    normalize();
  }

  /// @brief alternative operation to forget a set of elements
  /// @param elements a set of elements
  void forget(const std::vector<element_t> &elements) {
    if (is_bottom() || is_top()) {
      return;
    }

    for (auto v : elements) {
      *this -= v;
    }
  }

  /// @brief keep equivalence classes only for element in vector elements
  /// @param elements a vector of elements required to keep
  /// @note  After project operation, the state only keeps equivalence classes
  ///        for element in vector elements
  void project(const std::vector<element_t> &elements) {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_LOG("symb-var-eq", crab::outs()
                                << "Before projection " << *this << "\n";);
    // group elements based on current equivalence classes
    equivalence_class_elems_t elems;

    this_domain_t res;

    for (auto &v : elements) {
      if (contains(v)) {
        element_t rep = find(v);
        element_set_t &s = elems[rep];
        auto it = std::upper_bound(s.begin(), s.end(), v);
        s.insert(it, v);
      }
    }

    for (auto &kv : elems) {
      const element_t &rep = kv.first;
      const element_set_t &s = kv.second;
      bool keep_rep = false;
      for (auto &e : s) {
        if (e == rep) {
          keep_rep = true;
          break;
        }
      }
      element_t new_rep = keep_rep ? rep : s[0];
      for (auto &e : s) {
        res.m_parents.insert({e, new_rep});
      }
      equivalence_class_t ec = m_classes.at(rep);
      res.m_classes.insert({new_rep, std::move(ec)});
    }
    std::swap(res, *this);
    CRAB_LOG("symb-var-eq", crab::outs()
                                << "After projection " << *this << "\n";);
    normalize();
  }

  void rename(const std::vector<element_t> &old_elements,
              const std::vector<element_t> &new_elements) {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom()) {
      return;
    }
    if (old_elements.size() != new_elements.size()) {
      CRAB_ERROR(domain_name(),
                 "::rename with input vectors of different sizes");
    }
    for (unsigned i = 0, size = old_elements.size(); i < size; ++i) {
      const element_t &old_v = old_elements[i];
      const element_t &new_v = new_elements[i];
      if (!contains(old_v)) {
        CRAB_ERROR(domain_name(), "::rename assumes that ", old_v, " exists");
      }
      if (contains(new_v)) {
        CRAB_ERROR(domain_name(), "::rename assumes that ", new_v,
                   " does not exist");
      }
      boost::optional<element_t> parent_old_v;
      for (auto it = m_parents.begin(), et = m_parents.end(); it != et;) {
        if ((*it).first == old_v) {
          // find it as key, remove it
          parent_old_v = (*it).second;
          it = m_parents.erase(it);
        } else if ((*it).second == old_v) {
          // find it as value
          (*it).second = new_v;
          ++it;
        } else {
          ++it;
        }
        if (parent_old_v) {
          // insert a key-value pair for new_v
          m_parents.insert(
              {new_v, (*parent_old_v == old_v ? new_v : *parent_old_v)});
          auto it = m_classes.find(old_v);
          if (it != m_classes.end()) {
            // update m_classes
            equivalence_class_t ec = it->second;
            m_classes.erase(it);
            m_classes.insert({new_v, std::move(ec)});
          }
        }
      }
    }
  }

  /// @brief A normalization function used for other domain operations
  /// @note  perform normalization to eliminate any class such as:
  ///          {v3}=>#var5
  void normalize() {
    // FIXME: change this
    return;

    CRAB_LOG("symb-var-eq", crab::outs()
                                << "Before normalization " << *this << "\n";);
    equivalence_class_elems_t equiv_classes = equiv_classes_elems();
    for (auto &kv : equiv_classes) {
      const element_t &rep = kv.first;
      const element_set_t &s = kv.second;
      if (s.size() == 1) {
        // the equivalence class only remain one element
        m_classes.erase(rep);
        m_parents.erase(rep);
      }
    }
    if (m_classes.size() == 0) {
      set_to_top();
    }
    CRAB_LOG("symb-var-eq", crab::outs()
                                << "After normalization " << *this << "\n";);
  }

  void write(crab_os &o) const {
    if (is_top()) {
      o << "{}";
    } else if (is_bottom()) {
      o << "_|_";
    } else {
      equivalence_class_elems_t equiv_classes = equiv_classes_elems();
      CRAB_LOG("symb-var-eq-print", o << "("
                                      << "EquivClass=";
               print_equiv_classes(o, equiv_classes, true); o << ","
                                                              << "DomainVal=";
               print_classes_vals(o); o << ")"; return;);
      print_equiv_classes(o, equiv_classes);
    }
  }

  friend class crab::crab_os &operator<<(crab::crab_os &o,
                                         const this_domain_t &dom) {
    dom.write(o);
    return o;
  }

  friend class crab::crab_os &operator<<(crab::crab_os &o,
                                         const parents_map_t &map) {
    symbolic_variable_equiality_domain_impl::print_unordered_map(o, map);
    return o;
  }

  std::string domain_name() const { return "SymbolicVarEqDomain"; }
  /**------------------ End domain APIs ------------------**/
  // WARN: a special function to check equivalence classes over two domain
  // values
  bool is_included(const element_t &v, const this_domain_t &e) const {
    if (!contains(v) || !e.contains(v)) {
      return false;
    }
    const element_set_t left_s = get_all_members_from_an_equiv_class(v);
    const element_set_t right_s = e.get_all_members_from_an_equiv_class(v);
    return std::includes(left_s.begin(), left_s.end(), right_s.begin(),
                         right_s.end()) ||
           std::includes(right_s.begin(), right_s.end(), left_s.begin(),
                         left_s.end());
  }
};

} // end namespace domains
} // end namespace crab