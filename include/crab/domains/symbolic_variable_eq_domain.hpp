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
  using var_id_t = uint32_t;

  var_id_t var;

  symbolic_var(var_id_t _var) : var(_var) {}
  symbolic_var(const this_domain_t &o) = default;
  symbolic_var(this_domain_t &&o) = default;
  this_domain_t &operator=(const this_domain_t &o) = default;
  this_domain_t &operator=(this_domain_t &&o) = default;

  void write(crab_os &o) const { o << "#var" << var; }

  bool operator==(this_domain_t o) const { return var == o.var; }

  operator var_id_t() const { return var; }

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

/// @brief a copy-on-write wrapper for the domain value
/// @tparam Domain the type of equality domain
template <class Domain> class equivalence_class {
private:
  std::size_t m_size;
  std::shared_ptr<Domain> m_val;

  // Copy-on-write: call always this function before get_absval() if
  // m_val might be modified.
  void detach_absval() { m_val.reset(new Domain(*m_val)); }

public:
  explicit equivalence_class(std::shared_ptr<Domain> val)
      : m_size(0), m_val(val) {} // avoid copy initilization

  std::size_t &get_size() { return m_size; }

  std::size_t get_size() const { return m_size; }

  std::shared_ptr<Domain> detach_and_get_absval() {
    if (m_val.use_count() > 1) {
      detach_absval();
    }
    return m_val;
  }

  std::shared_ptr<const Domain> get_absval() const { return m_val; }

  void set_absval(std::shared_ptr<Domain> val) { m_val = val; }
}; // end class equivalence_class

/// @brief helper function for logging and debugging map
/// @tparam Key the type of value for keys
/// @tparam Value the type of value for values
/// @param o crab os stream
/// @param m unorder map object
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

class SVEQDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};

#define SVEQ_DOMAIN_SCOPED_STATS(NAME) CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define SVEQ_COUNT_STATS(NAME) CRAB_DOMAIN_COUNT_STATS(NAME, 0)

/// @brief An abstract domain represents equalities used in analyses such as
/// allocation-site abstraction or domain reduction. In short, giving two
/// elements are known to hold equal values in concrete semantics if
/// the domain captures that the elements hold equal symbolic variable.
/// @tparam BaseDomain the type of the domain value represents the base domain in the object domain
/// @tparam DomainParams domain parameter for inter-procedural analysis
template <class BaseDomain, class DomainParams = SVEQDefaultParams> 
class symbolic_variable_equiality_domain final
  : public abstract_domain_api<symbolic_variable_equiality_domain<BaseDomain, DomainParams>> {
public:
  using symb_eq_domain_t = symbolic_variable_equiality_domain<BaseDomain, DomainParams>;
  using abstract_domain_t = abstract_domain_api<symb_eq_domain_t>;

  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using number_t = typename BaseDomain::number_t;

  // typedefs for equality domain
  using element_t = variable_t;
  using element_set_t = std::vector<element_t>;
  using domain_t = class symbolic_variable_equiality_domain_impl::symbolic_var;
  using var_id_t = typename domain_t::var_id_t;
  using parents_map_t = std::unordered_map<element_t, element_t>;
  using equivalence_class_t =
      symbolic_variable_equiality_domain_impl::equivalence_class<domain_t>;
  using equivalence_class_elems_t =
      std::unordered_map<element_t, element_set_t>;

private:
  using this_domain_t = symb_eq_domain_t;
  using classes_map_t = std::unordered_map<element_t, equivalence_class_t>;
  enum class lattice_val { bottom, top, neither_top_nor_bot };

  // a map that stores a variable to its immediate representative
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
  // The implementation is path compression based. So there is no case that
  // a variable maps to another variable whic is not representative in m_parents
  // (i.e. no tree like representation for each equivalent class)

  /// @brief a helper method to empty the disjoin set
  void clear() {
    m_parents.clear();
    m_classes.clear();
  }

  /// @brief Build a map from representative to a ordered set with all the
  /// elements in the equivalence class.
  /// @return a map computed as brief described.
  equivalence_class_elems_t equiv_classes_elems() {
    SVEQ_DOMAIN_SCOPED_STATS(".classes");

    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t &rep = kv.second; // already path compressed
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
    SVEQ_DOMAIN_SCOPED_STATS(".classes");

    equivalence_class_elems_t res;
    for (auto &kv : m_parents) {
      element_t rep = kv.second; // already path compressed
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
      element_t rep = kv.second;
      if (rep == e_rep) {
        auto it = std::upper_bound(out.begin(), out.end(), kv.first);
        out.insert(it, kv.first);
      }
    }
    return out; // no need std::move since Return Value Optimization (RVO) is enabled
  }

  /// @brief a helper method for domain operation join
  /// @param left one abstract state
  /// @param right one abstract state
  /// @return a new abstract state saved the joined result
  /// @details Computing the intersection of two equivalence classes
  /// to find least common elements. It only computes a set of
  /// equivalence classes where each class is subset of two classes from the
  /// given two abstract states. That is, the result of the join only remains
  /// equalities appeared both in left and in right states.
  /// Formally, forall cls_x \in left. forall cls_y \in right ::
  ///           cls_x ^ cls_y
  /// 
  /// This operation is running on Quadratic Time -- O(n^2).
  /// The computation follows the structure of unordered_map.
  /// The operation requires to construct a new unordered_map.
  this_domain_t join(const this_domain_t &left,
                         const this_domain_t &right) const {
    this_domain_t res;
    var_id_t current_id = 0;
    // res.dump(); 
    for (auto it = left.m_parents.begin(); it != left.m_parents.end(); it++) {
      // for each equality k == v
      const element_t &k = it->first;
      const element_t &v = it->second;
      if (k == v) {
        continue;
      }
      auto it_k = right.m_parents.find(k);
      auto it_v = right.m_parents.find(v);
      // check if k == v exists in another map
      if (it_k != right.m_parents.end() && it_v != right.m_parents.end() &&
          it_k->second == it_v->second) {
        res.m_parents.insert({v, k}); // insert new pair <v, k>
        if (!res.contains(k)) {
          std::shared_ptr<domain_t> absval_ptr =
              std::make_shared<domain_t>(current_id);
          res.make_set(k, absval_ptr);
          current_id++;
        }
      }
      for (auto it2 = std::next(it); it2 != left.m_parents.end(); it2++) {
        // for each k2 == v2
        const element_t &k2 = it2->first;
        const element_t &v2 = it2->second;
        if (v != v2 || (k == k2 && v == v2)) {
          // k2 is not in k's class or k2, v2 is the same as k, v
          // skip
          continue;
        }
        // check if k2 == k exists in another map
        auto it_k2 = right.m_parents.find(k2);
        if (it_k != right.m_parents.end() && it_k2 != right.m_parents.end() &&
            it_k2->second == it_k->second) {
          res.m_parents.insert({k2, k}); // insert new pair <k2, k>
          if (!res.contains(k)) {
            std::shared_ptr<domain_t> absval_ptr =
                std::make_shared<domain_t>(current_id);
            res.make_set(k, absval_ptr);
            current_id++;
          }
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
  /// The solution is merging two equivalence classes if two equivalence classes
  /// shared common elements. As the result, it computed a set of a set of
  /// equivalence classes where each class is a super set of two or more classes
  /// from the given two abstract states. That is, the result of the meet will
  /// keep equalities appeared either in left or in right states.
  /// Bt definition, our meet is an over-approximation meet of concrete.
  /// E.g. there is no bottom state after meet.
  /// Formally, forall cls_x \in left. forall cls_y \in right ::
  ///           cls_x ^ cls_y != \empty => cls_x union cls_y
  this_domain_t meet(const this_domain_t &left,
                     const this_domain_t &right) const {
    this_domain_t res;
    var_id_t current_id = 0;
    for (auto it = left.m_parents.begin(); it != left.m_parents.end(); it++) {
      // for each equality k == v
      const element_t &k = it->first;
      const element_t &v = it->second;
      // k.dump(crab::outs());
      // res.dump();
      // skip if k already in some class on the res state
      // why exists? k is already merged from the right state, so skip it.
      if (res.m_parents.find(k) != res.m_parents.end()) {
        continue;
      }
      bool is_k_added = false;
      // if some of k's equivlence member k_p from the left state exists on the
      // res state, this mean we also need to add equality k == k_p on the res
      element_set_t k_cls = left.get_all_members_from_an_equiv_class(k);
      for (auto &k_p : k_cls) {
        if (res.m_parents.find(k_p) != res.m_parents.end()) {
          res.add(k_p, k);
          is_k_added = true;
          break;
        }
      }
      if (!is_k_added) {
        // insert k and v into res map since k == v will keep
        res.m_parents.insert({k, v});
        if (!res.contains(v)) {
          std::shared_ptr<domain_t> absval_ptr =
              std::make_shared<domain_t>(current_id);
          res.make_set(v, absval_ptr);
          current_id++;
        }
      }

      // check if k does not in some class on the right state
      // , we do not need to merge any classes on the right state
      auto it2 = right.m_parents.find(k);
      if (it2 == right.m_parents.end()) {
        continue;
      }
      const element_t &v2 = it2->second;
      for (auto it3 = right.m_parents.begin(); it3 != right.m_parents.end(); it3++) {
        // for each k3 == v3
        const element_t &k3 = it3->first;
        const element_t &v3 = it3->second;
        if (v2 == v3 && k3 != k) {
          // inside right's equivalence class including k, any variables equal to k will keep
          res.m_parents.insert({k3, res.m_parents.find(k)->second});
        }
      }
    }

    // Insert equalities remained on the right only
    for (auto it = right.m_parents.begin(); it != right.m_parents.end(); it++) {
      // for each equality k == v
      const element_t &k = it->first;
      const element_t &v = it->second;
      if (res.m_parents.find(k) == res.m_parents.end()) {
        res.m_parents.insert({k, v});
        if (!res.contains(v)) {
          std::shared_ptr<domain_t> absval_ptr =
              std::make_shared<domain_t>(current_id);
          res.make_set(v, absval_ptr);
          current_id++;
        }
      }
    }
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

  /// @brief find the representative without path-compression
  /// @param v an element in some set
  /// @attention v must constain in current domain; otherwise, an error returns
  /// @return returns the representative of the set that contains the element v
  element_t find(const element_t &v) const {
    SVEQ_DOMAIN_SCOPED_STATS(".find");

    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      CRAB_ERROR(domain_name(), "::", __func__, " on a non-existing elem ", v,
                 " in ", *this);
    }
    return it->second;
  }

  boost::optional<element_t> find_opt(const element_t &v) const {
    SVEQ_DOMAIN_SCOPED_STATS(".find");

    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      return boost::none;
    }
    return it->second;
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

    std::shared_ptr<domain_t> absval_ptr =
        std::make_shared<domain_t>(std::move(absval));

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
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top() || !contains(x)) {
      return boost::none;
    }
    return get_all_members_from_an_equiv_class(x);
  }

  /// @brief get the equivalent class by giving a domain value
  /// @param dom a domain value such as #var1
  /// @return if the value exists, return corresponding element class
  boost::optional<element_set_t> get_variables(const domain_t &dom) const {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top()) {
      return boost::none;
    }
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

  boost::optional<element_set_t> get_all_variables() const {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top()) {
      return boost::none;
    }
    return get_all_members();
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
  symbolic_variable_equiality_domain(const this_domain_t &o) = default;

  symbolic_variable_equiality_domain(this_domain_t &&o) = default;
  this_domain_t &operator=(const this_domain_t &o) = default;
  this_domain_t &operator=(this_domain_t &&o) = default;

  this_domain_t make_bottom() const override {
    this_domain_t res(lattice_val::bottom);
    return res;
  }

  this_domain_t make_top() const override {
    this_domain_t res(lattice_val::top);
    return res;
  }

  bool is_bottom() const override { return m_val == lattice_val::bottom; }

  bool is_top() const override { return m_val == lattice_val::top; }

  void set_to_top() override {
    clear();
    m_val = lattice_val::top;
  }

  void set_to_bottom() override {
    clear();
    m_val = lattice_val::bottom;
  }

  void set_neither_top_or_bottom() { m_val = lattice_val::neither_top_nor_bot; }

  bool operator<=(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".leq");
    if (is_bottom() || e.is_top()) {
      // _|_ <= e || this <= top
      return true;
    } else if (is_top() || e.is_bottom()) {
      // top >= e || this >= _|_
      return false;
    }
    // Return true if *this is a refined partitioning of o
    //    \forall cls_x \in *this. \exists cls_y \in e ::
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

  this_domain_t operator|(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".join");
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

  void operator|=(const this_domain_t &e) override {
    SVEQ_DOMAIN_SCOPED_STATS(".join");
    if (is_bottom()) {
      if (!e.is_bottom()) {
        *this = e;
      }
      return;
    } else if (e.is_bottom()) {
      return;
    } else if (is_top() || e.is_top()) {
      set_to_top();
      return;
    } else {
      *this = std::move(join(*this, e));
      return;
    }
  }

  this_domain_t operator&(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".meet");
    if (is_bottom() || e.is_top()) {
      return *this;
    } else if (e.is_bottom() || is_top()) {
      return e;
    } else {
      return meet(*this, e);
    }
  }

  void operator&=(const this_domain_t &e) override {
    SVEQ_DOMAIN_SCOPED_STATS(".meet");
    if (is_bottom() || e.is_top()) {
      return;
    } else if (e.is_bottom() || is_top()) {
      *this = e;
      return;
    } else {
      *this = std::move(meet(*this, e));
      return;
    }
  }

  this_domain_t operator||(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".widening");
    return *this;
  }

  this_domain_t operator&&(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".narrowing");
    return *this;
  }

  this_domain_t widening_thresholds(const this_domain_t &abs, const thresholds<number_t> &ts) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".widening");
    return *this;
  }

  /// @brief expand what x equals to to y. This is equivalent to \c add(x,y)
  /// @param x the orginal variable in some class
  /// @param y a new variable that has same equalities as \p x
  void expand(const element_t &x, const element_t &y) override {
    if (is_bottom() || is_top() || !contains(x)) {
      return;
    }
    add(x, y);
  }

  NUMERICAL_OPERATIONS_NOT_IMPLEMENTED(this_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(this_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(this_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(this_domain_t)

  /// @brief remove an element from a class
  /// @param v an element
  /// @note  if v \in cls and v is the representative, select a new representative
  ///        before removing v.
  void operator-=(const element_t &v) override {
    SVEQ_DOMAIN_SCOPED_STATS(".forget");
    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_LOG("symb-var-eq", crab::outs() << "Forgetting " << v << ": ";
             dump(););

    if (contains(v)) {
      // v is not a representative, jusrt remove v
      // if not, pick a new one and update all members
      if (m_parents.at(v) == v) {
        // v is the representative of the equivalence class
        boost::optional<element_t> new_rep;
        for (auto &kv : m_parents) {
          // search element which is not v but in v's class
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
    CRAB_LOG("symb-var-eq", crab::outs() << "After forget "; dump(););
  }

  /// @brief alternative operation to forget a set of elements
  /// @param elements a set of elements
  void forget(const std::vector<element_t> &elements) override {
    if (is_bottom() || is_top()) {
      return;
    }

    // TODO: if necessary, provide a dirtect implementation without too much
    // replacing representative
    for (auto v : elements) {
      *this -= v;
    }
  }

  /// @brief keep equivalence classes only for element in vector elements
  /// @param elements a vector of elements required to keep
  /// @note  After project operation, the state only keeps equivalence classes
  ///        for element in vector elements
  void project(const std::vector<element_t> &elements) override {
    SVEQ_DOMAIN_SCOPED_STATS(".project");

    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_LOG("symb-var-eq", crab::outs() << "Projecting ";
             print_elems_vector(crab::outs(), elements); crab::outs() << "\n";
             dump(););
    // group elements based on current equivalence classes
    equivalence_class_elems_t elems;

    this_domain_t res;

    for (auto &v : elements) {
      if (contains(v)) {
        element_t rep = find(v);
        element_set_t &s = elems[rep];
        s.push_back(v);
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
      res.m_classes.insert({new_rep, ec});
    }
    std::swap(res, *this);
    normalize();
    CRAB_LOG("symb-var-eq", crab::outs() << "After projection "; dump(););
  }

  void rename(const std::vector<element_t> &old_elements,
              const std::vector<element_t> &new_elements) override {
    SVEQ_DOMAIN_SCOPED_STATS(".rename");

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

  // Reduce the size of the abstract domain representation.
  void minimize() override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  /// @brief A normalization function used for other domain operations
  /// @note  perform normalization to eliminate any class such as:
  ///          {v3}=>#var5
  void normalize() override {
    SVEQ_DOMAIN_SCOPED_STATS(".normalize");

    CRAB_LOG("symb-var-eq", crab::outs()
                                << "Before normalization " << *this << "\n";);
    // TODO: could implement to check each size of class instead of creating
    // equivalence_class_elems_t map
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

  void write(crab_os &o) const override {
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

  void dump() const {
    if (is_top()) {
      crab::outs() << "{}";
    } else if (is_bottom()) {
      crab::outs() << "_|_";
    } else {
      equivalence_class_elems_t equiv_classes = equiv_classes_elems();
      crab::outs() << "(EquivClass=";
      print_equiv_classes(crab::outs(), equiv_classes, true);
      crab::outs() << ", DomainVal=";
      print_classes_vals(crab::outs());
      crab::outs() << ")";
      crab::outs() << "HashTable=";
      symbolic_variable_equiality_domain_impl::print_unordered_map<element_t,
                                                                   element_t>(
          crab::outs(), m_parents);
    }
    crab::outs() << "\n";
  }

  std::string domain_name() const override { return "EqDomain"; }

  linear_constraint_system_t to_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  disjunctive_linear_constraint_system_t to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs, const variable_vector_t &outputs) override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_domain_t &invariant) override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void callee_entry(const callsite_info<variable_t> &callsite, const this_domain_t &caller) override {
    SVEQ_DOMAIN_SCOPED_STATS(".callee_entry");
    inter_abstract_operations<this_domain_t, DomainParams::implement_inter_transformers>::
      callee_entry(callsite, caller, *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite, const this_domain_t &callee) override {
    SVEQ_DOMAIN_SCOPED_STATS(".caller_cont");
    inter_abstract_operations<this_domain_t, DomainParams::implement_inter_transformers>::
      callee_entry(callsite, callee, *this);
  }
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

template <typename BaseDomain, typename DomainParams>
struct abstract_domain_traits<symbolic_variable_equiality_domain<BaseDomain, DomainParams>> {
  using number_t = typename BaseDomain::number_t;
  using varname_t = typename BaseDomain::varname_t;
};
} // end namespace domains
} // end namespace crab