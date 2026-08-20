#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/constant.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

namespace crab {
namespace domains {

namespace symbolic_variable_equality_domain_impl {
// Identifier of an equivalence class. Two variables are known equal when their
// classes carry the same id. Tagging classes with ids (rather than a
// representative variable) is a design choice: it allows an equality relation
// to be split across several domain values, where setting the same id in two
// values records an equality between their elements. Because ids may travel
// between values, they must be globally unique -- an id may appear in two
// places only because it was deliberately copied there. A per-value allocator
// cannot provide this (two values can allocate the same id independently),
// hence the single shared counter below.
using class_id_t = uint32_t;

// The single, monotonically-increasing source of class ids, shared by all
// domain instances. Assumes single-threaded use.
inline class_id_t &class_id_counter() {
  static class_id_t counter = 0;
  return counter;
}

// Produce a fresh class id, strictly greater than every id produced before
// and therefore different from every id in use anywhere.
inline class_id_t fresh_class_id() {
  class_id_t &counter = class_id_counter();
  if (counter == std::numeric_limits<class_id_t>::max()) {
    CRAB_ERROR("symbolic_variable_equality_domain: class id counter overflow");
  }
  return counter++;
}

// True iff `id` was produced by fresh_class_id() at some point.
inline bool is_valid_class_id(class_id_t id) { return id < class_id_counter(); }
} // namespace symbolic_variable_equality_domain_impl

class SVEQDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
  // when non-zero, domain operations validate their internal representation
  // (see check_lattice_val). Off by default so normal analyses pay nothing.
  enum { check_lattice_val = 0 };
  // when non-zero, operations drop singleton classes (which carry no equality)
  // so a singleton-only state collapses to top. On by default.
  enum { normalize = 1 };
};

// Like SVEQDefaultParams but with internal consistency checks enabled. Use for
// unit testing, or pass it when debugging a specific analysis.
class SVEQCheckedParams {
public:
  enum { implement_inter_transformers = 0 };
  enum { check_lattice_val = 1 };
  enum { normalize = 1 };
};

// Like SVEQDefaultParams but keeps singleton classes (disables normalize). Use
// where a lone element's class must be retained, e.g. object_domain's field /
// register equalities.
class SVEQNoNormalizeParams {
public:
  enum { implement_inter_transformers = 0 };
  enum { check_lattice_val = 0 };
  enum { normalize = 0 };
};

#define SVEQ_DOMAIN_SCOPED_STATS(NAME) CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define SVEQ_COUNT_STATS(NAME) CRAB_DOMAIN_COUNT_STATS(NAME, 0)

/// @brief An abstract domain that represents equalities, used in analyses
/// such as allocation-site abstraction or domain reduction. In short, two
/// elements are known to hold equal values in the concrete semantics if
/// the domain captures that the elements are assigned the same symbolic
/// variable.
/// @tparam Number numeric type of the program variables
/// @tparam VariableName name type of the program variables
/// @tparam DomainParams domain parameters (normalize, checks, inter-procedural)
template <typename Number, typename VariableName,
          typename DomainParams = SVEQDefaultParams>
class symbolic_variable_equality_domain final
    : public abstract_domain_api<symbolic_variable_equality_domain<
          Number, VariableName, DomainParams>> {
public:
  using symb_eq_domain_t =
      symbolic_variable_equality_domain<Number, VariableName, DomainParams>;
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
  using number_t = Number;
  using varname_t = VariableName;

  // typedefs for equality domain
  using element_t = variable_t;
  using element_set_t = std::vector<element_t>;
  using class_id_t = symbolic_variable_equality_domain_impl::class_id_t;
  using parents_map_t = std::unordered_map<element_t, element_t>;
  using equivalence_class_elems_t =
      std::unordered_map<element_t, element_set_t>;

private:
  using this_domain_t = symb_eq_domain_t;
  using classes_map_t = std::unordered_map<element_t, class_id_t>;
  enum class lattice_val { bottom, top, neither_top_nor_bot };

  // a map that stores a variable to its immediate representative
  // The map is a many to one hash map
  parents_map_t m_parents;
  // a map that keeps a representative to its corresponding symbolic variable
  // each representative has its equivalence class
  classes_map_t m_classes;

  lattice_val m_val;
  // For example,
  // The disjoint set for a state of symbolic_variable_equality_domain<int>:
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
  // it does not capture any relations between elements.
  // The implementation is path compression based. So there is no case that
  // a variable maps to another variable which is not a representative in
  // m_parents (i.e. no tree like representation for each equivalent class)

  /// @brief a helper method to empty the disjoint set
  void clear() {
    m_parents.clear();
    m_classes.clear();
  }

  void try_make_set_by_raw_val(const element_t &v, class_id_t id) {
    if (!contains(v)) {
      make_set(v, id);
    }
  }

  // Return the representative of the class currently tagged with `id`, if any.
  // Backs the "same id => same class" rule: set()/add() merge into the existing
  // class instead of duplicating the id.
  boost::optional<element_t> rep_with_class_id(class_id_t id) const {
    for (const auto &kv : m_classes) {
      if (kv.second == id) {
        return kv.first;
      }
    }
    return boost::none;
  }

  /// @brief Build a map from representative to an ordered set with all the
  /// elements in the equivalence class.
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
    return out;
  }

  /// @brief join helper: keep only the equalities present in BOTH operands.
  /// @param left,right the operands.
  /// @pre left and right are both neither top nor bottom; operator| handles
  ///      those cases before reaching this helper.
  /// @return a new state in which x == y holds iff it holds in left AND right.
  /// @details For every pair (x,y) that is an equality in `left`
  ///          (left.find(x) == left.find(y)), keep it iff it is also an
  ///          equality in `right`. A pair equal in both is necessarily equal in
  ///          left, so iterating left's equalities and filtering by right is
  ///          sufficient. Quadratic in the number of variables.
  this_domain_t join(const this_domain_t &left,
                     const this_domain_t &right) const {
    this_domain_t res;
    CRAB_LOG("symb-var-eq-join", crab::outs() << "Join "; left.dump();
             crab::outs() << " and "; right.dump(););
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
        res.try_make_set_by_raw_val(k, fresh_class_id());
        res.m_parents.insert({v, k}); // insert new pair <v, k>
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
          res.try_make_set_by_raw_val(k, fresh_class_id());
          res.m_parents.insert({k2, k}); // insert new pair <k2, k>
        }
      }
    }
    res.normalize();
    res.check_lattice_val();
    CRAB_LOG("symb-var-eq-join", res.dump(); crab::outs() << "\n");
    return res;
  }

  /// @brief meet helper: keep every equality from EITHER operand (their union).
  /// @param left,right the operands.
  /// @pre left and right are both neither top nor bottom; operator& handles
  ///      those cases before reaching this helper.
  /// @return a new state in which x == y holds iff it holds in left OR right
  ///         (the union of both equivalence relations; never bottom).
  /// @details Start from a copy of left, then for each equality k == v in right
  ///          (i.e. k != v in right.m_parents) union the classes of k and v.
  ///          Symbols of the result are not significant, so merged classes keep
  ///          left's representative/symbol.
  this_domain_t meet(const this_domain_t &left,
                     const this_domain_t &right) const {
    this_domain_t res = left;
    CRAB_LOG("symb-var-eq-meet", crab::outs() << "Meet "; left.dump();
             crab::outs() << " and "; right.dump(););
    for (auto &kv : right.m_parents) {
      const element_t &k = kv.first;
      const element_t &v = kv.second;
      if (k != v) { // k == v is an equality in right (v is k's representative)
        res.add(k, v);
      }
    }
    res.normalize();
    res.check_lattice_val();
    CRAB_LOG("symb-var-eq-meet", res.dump(); crab::outs() << "\n");
    return res;
  }

  /// @brief create a new equivalence class {v} tagged with class id `id`
  /// @param v the representative element for the new class
  /// @param id the class id
  void make_set(const element_t &v, class_id_t id) {
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
    if (DomainParams::check_lattice_val && rep_with_class_id(id)) {
      // each equivalence class must carry a distinct id
      CRAB_ERROR(domain_name(), "::", __func__, " class id #id", id,
                 " is already used by another class in ", *this);
    }

    m_parents.insert({v, v});
    m_classes.insert({v, id});
  }

  // Consistency check of the internal representation. Enabled via
  // DomainParams::check_lattice_val (e.g. SVEQCheckedParams); compiles to a
  // no-op otherwise. Validates the invariants linking m_val, m_parents and
  // m_classes.
  void check_lattice_val() const {
    if (!DomainParams::check_lattice_val) {
      return;
    }
    switch (m_val) {
    case lattice_val::bottom:
      if (!m_parents.empty() || !m_classes.empty()) {
        CRAB_ERROR(domain_name(),
                   "::check_lattice_val: bottom must have empty maps");
      }
      break;
    case lattice_val::top:
      if (!m_parents.empty() || !m_classes.empty()) {
        CRAB_ERROR(domain_name(),
                   "::check_lattice_val: top must have empty maps");
      }
      break;
    case lattice_val::neither_top_nor_bot:
      if (m_parents.empty() && m_classes.empty()) {
        CRAB_ERROR(domain_name(), "::check_lattice_val: neither-top-nor-bottom "
                                  "must not be empty (it should be top)");
      }
      if (!m_parents.empty() && m_classes.empty()) {
        CRAB_ERROR(domain_name(), "::check_lattice_val: m_parents is non-empty "
                                  "but there are no equivalence classes");
      }
      // every element points to a representative that owns an equivalence class
      // and every representative is self-parented
      for (auto &kv : m_parents) {
        const element_t &rep = kv.second;
        auto pit = m_parents.find(rep);
        if (pit == m_parents.end() || !(pit->second == rep)) {
          CRAB_ERROR(domain_name(), "::check_lattice_val: representative ", rep,
                     " of ", kv.first, " is not self-parented in m_parents");
        }
        if (m_classes.find(rep) == m_classes.end()) {
          CRAB_ERROR(domain_name(), "::check_lattice_val: representative ", rep,
                     " has no equivalence class in m_classes");
        }
      }
      // every equivalence class is keyed by a self-parented representative and
      // carries an id that the factory actually produced
      for (auto &kv : m_classes) {
        const element_t &rep = kv.first;
        auto pit = m_parents.find(rep);
        if (pit == m_parents.end() || !(pit->second == rep)) {
          CRAB_ERROR(domain_name(), "::check_lattice_val: class key ", rep,
                     " is not a representative in m_parents");
        }
        if (!symbolic_variable_equality_domain_impl::is_valid_class_id(
                kv.second)) {
          CRAB_ERROR(domain_name(), "::check_lattice_val: class id #id",
                     kv.second, " was never produced by fresh_class_id()");
        }
      }
      break;
    }
  }

  // True iff this state encodes no equalities, i.e. it is semantically top:
  // either flagged top, or every class is a singleton (each element is its own
  // representative). Lets leq treat a not-yet-normalized singleton-only state
  // the same as top.
  bool has_no_equalities() const {
    if (is_top()) {
      return true;
    }
    if (is_bottom()) {
      return false;
    }
    for (auto &kv : m_parents) {
      if (!(kv.first == kv.second)) {
        return false;
      }
    }
    return true;
  }

  // Print each element of `r` via `fn`, wrapped in open/close and joined by
  // sep.
  template <typename Range, typename PrintFn>
  void print_separated(crab_os &o, const Range &r, PrintFn fn, char open,
                       char close, const char *sep = ",") const {
    o << open;
    bool first = true;
    for (const auto &elem : r) {
      if (!first) {
        o << sep;
      }
      first = false;
      fn(elem);
    }
    o << close;
  }

  void print_elems_vector(crab::crab_os &o,
                          const std::vector<element_t> &elems) const {
    print_separated(
        o, elems, [&](const element_t &e) { o << e; }, '[', ']');
  }

  void print_equiv_classes(crab_os &o,
                           const equivalence_class_elems_t &equiv_classes,
                           bool verbose = false) const {
    print_separated(
        o, equiv_classes,
        [&](const std::pair<const element_t, element_set_t> &kv) {
          print_elems_vector(o, kv.second);
          if (!verbose) {
            o << "=>#id" << m_classes.at(kv.first);
          }
        },
        '{', '}');
  }

  void print_classes_vals(crab_os &o) const {
    print_separated(
        o, m_classes,
        [&](const std::pair<const element_t, class_id_t> &kv) {
          o << kv.first << "=>#id" << kv.second;
        },
        '{', '}');
  }

public:
  /**------------------ Begin union find APIs ------------------**/
  /// @brief Check whether domain contains element v in some class
  /// @param v an element
  /// @return true if the element exists; otherwise, false.
  bool contains(const element_t &v) const {
    return m_parents.find(v) != m_parents.end();
  }

  /// @brief Check whether two elements are in the same class, even if they
  /// may not exist
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
  /// @attention v must be contained in the current domain; otherwise, an
  /// error is reported
  /// @return returns the representative of the set that contains the element v
  element_t find(const element_t &v) const {
    if (auto rep = find_opt(v)) {
      return *rep;
    }
    CRAB_ERROR(domain_name(), "::", __func__, " on a non-existing elem ", v,
               " in ", *this);
  }

  boost::optional<element_t> find_opt(const element_t &v) const {
    SVEQ_DOMAIN_SCOPED_STATS(".find");

    auto it = m_parents.find(v);
    if (it == m_parents.end()) {
      return boost::none;
    }
    return it->second;
  }

  /// @brief Produce a fresh class id, different from every id ever produced.
  /// @note  Ids passed to set() must originate from here (either directly or
  ///        read back from a class via get_class_id()).
  static class_id_t fresh_class_id() {
    return symbolic_variable_equality_domain_impl::fresh_class_id();
  }

  /// @brief tag x's class with class id \p id
  /// @param x an element
  /// @param id a class id previously produced by fresh_class_id()
  /// @note  if x is new and no class already holds \p id, create a class {x}
  ///        with \p id;
  ///        if some class already holds \p id, x's class is merged into it --
  ///        i.e. sharing an id asserts an equality;
  ///        otherwise relabel x's class to \p id.
  void set(const element_t &x, class_id_t id) {
    if (is_bottom()) {
      return;
    }
    if (!symbolic_variable_equality_domain_impl::is_valid_class_id(id)) {
      // an id the factory never produced could later be handed out by
      // fresh_class_id(), silently merging unrelated classes
      CRAB_ERROR(domain_name(), "::set: class id #id", id,
                 " was never produced by fresh_class_id()");
    }

    if (!contains(x)) {
      if (auto holder = rep_with_class_id(id)) {
        add(*holder, x); // id already in use: x joins that class
      } else {
        make_set(x, id);
      }
    } else {
      element_t rep_x = find(x);
      if (auto holder = rep_with_class_id(id)) {
        if (!(*holder == rep_x)) {
          // another class already owns this id: relabeling x's class to it
          // asserts the two are equal, so union them instead of duplicating
          add(*holder, x);
        }
        // else: x's class already holds this id -> nothing to do
      } else {
        // id unused: relabel x's whole class
        m_classes.at(rep_x) = id;
      }
    }
    check_lattice_val();
  }

  /// @brief get the equivalent class by giving an element
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

  /// @brief get the members of the class tagged with class id \p id
  /// @param id a class id such as #id1
  /// @return if a class has that id, return its members
  boost::optional<element_set_t> get_variables(class_id_t id) const {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top()) {
      return boost::none;
    }
    if (auto rep = rep_with_class_id(id)) {
      return get_all_members_from_an_equiv_class(*rep);
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
    element_set_t out;
    out.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      out.push_back(kv.first);
    }
    return out;
  }

  /// @brief get the class id of x's equivalence class
  /// @param x an element in some set
  /// @return none if x is not in a class, otherwise its class id
  boost::optional<class_id_t> get_class_id(const element_t &x) const {
    if (is_bottom()) {
      CRAB_ERROR("called ", domain_name(), "::", __func__, " on bottom");
    }
    if (is_top() || !contains(x)) {
      return boost::none;
    }
    return m_classes.at(find(x));
  }

  /// @brief Assert x == y by unioning their equivalence classes.
  /// @param x,y elements; either may be new (a fresh class is created first).
  /// @note  If y already belongs to a class, that whole class is merged in, so
  ///        y's existing equalities are preserved (this is a union, not a
  ///        move). The merged class keeps x's representative and symbol. Adding
  ///        an equality to top creates a class (top is not absorbing).
  void add(const element_t &x, const element_t &y) {
    if (is_bottom() || x == y) {
      return;
    }
    if (!contains(x)) {
      // Create an isolated fresh class for x. The factory is monotonic and
      // set() only accepts ids the factory produced, so a fresh id is greater
      // than every id in use and cannot land on an existing class.
      try_make_set_by_raw_val(x, fresh_class_id());
    }
    element_t rep_x = find(x);
    if (!contains(y)) {
      m_parents.insert({y, rep_x}); // y is new: place it in x's class
    } else {
      element_t rep_y = find(y);
      if (rep_x != rep_y) {
        // union: repoint every member of y's class to x's representative,
        // then drop y's now-empty class entry
        for (auto &kv : m_parents) {
          if (kv.second == rep_y) {
            kv.second = rep_x;
          }
        }
        m_classes.erase(rep_y);
      }
    }
    check_lattice_val();
  }
  /**------------------ End union find APIs ------------------**/

  /**------------------ Begin domain APIs ------------------**/
  // A default-constructed domain is top (the empty union-find: no equalities
  // known). It transitions to neither_top_nor_bot once an equality is recorded
  // (via make_set, reached from set/add).
  symbolic_variable_equality_domain(lattice_val val = lattice_val::top)
      : m_val(val) {}
  symbolic_variable_equality_domain(const this_domain_t &o) = default;

  symbolic_variable_equality_domain(this_domain_t &&o) = default;
  this_domain_t &operator=(const this_domain_t &o) = default;
  this_domain_t &operator=(this_domain_t &&o) = default;

  this_domain_t make_bottom() const override {
    this_domain_t res(lattice_val::bottom);
    return res;
  }

  DEFAULT_MAKE_PROJECTION(this_domain_t)
  DEFAULT_MAKE_FORGET(this_domain_t)

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
    check_lattice_val();
    e.check_lattice_val();
    if (is_bottom() || e.is_top()) {
      // _|_ <= e || this <= top
      return true;
    } else if (e.is_bottom()) {
      // this (non-_|_) is not <= _|_
      return false;
    } else if (is_top()) {
      // top has no equalities, so top <= e iff e has none either (e may be
      // top-like without being flagged top, e.g. a not-yet-normalized
      // singleton)
      return e.has_no_equalities();
    }
    CRAB_LOG("symb-var-eq-leq", crab::outs() << "Inclusion test: "; dump();
             crab::outs() << " and "; e.dump(););
    // equals
    // Return true if *this is a refined partitioning of e
    //    \forall cls_y \in e. \exists cls_x \in *this ::
    //      the set cls_y is a subset of cls_x
    // e.g. {[v3,v4]=>#var0,[v1,v2]=>#var1} \leq {[v1,v2]=>#var0}
    //      {[v3,v4]=>#var0,[v1,v2]=>#var1} \leq {[v1,v2]=>#var0, [v5]=>#var2}
    //      {[v3,v4]=>#var0,[v1,v2]=>#var1} \leq {[v1,v2]=>#var0, [v3]=>#var2}
    // or alternative, \for each equality \in e, find equality \in *this
    bool res = true;
    for (auto it = e.m_parents.begin(); it != e.m_parents.end(); it++) {
      const element_t &k1 = it->first;
      const element_t &v1 = it->second;
      for (auto it2 = std::next(it); it2 != e.m_parents.end(); it2++) {
        const element_t &k2 = it2->first;
        const element_t &v2 = it2->second;
        if (v1 == v2) { // k1 == k2
          auto itk1 = m_parents.find(k1);
          auto itk2 = m_parents.find(k2);
          res &= (itk1 != m_parents.end() && itk2 != m_parents.end() &&
                  itk1->second == itk2->second);
          if (!res) {
            break;
          }
        }
      }
    }
    CRAB_LOG("symb-var-eq-leq", crab::outs()
                                    << (res ? "true" : "false") << "\n");
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
    return this->operator|(e);
  }

  this_domain_t operator&&(const this_domain_t &e) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".narrowing");
    return this->operator&(e);
  }

  this_domain_t
  widening_thresholds(const this_domain_t &abs,
                      const thresholds<number_t> &ts) const override {
    SVEQ_DOMAIN_SCOPED_STATS(".widening");
    return this->operator|(abs);
  }

  /// @brief expand what x is equal to into y. This is equivalent to \c add(x,y)
  /// @param x the original variable in some class
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
  /// @note  if v \in cls and v is the representative, select a new
  /// representative
  ///        before removing v.
  void operator-=(const element_t &v) override {
    SVEQ_DOMAIN_SCOPED_STATS(".forget");
    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_LOG("symb-var-eq", crab::outs() << "Forgetting " << v << ": ";
             dump(););

    if (contains(v)) {
      // v is not a representative, just remove v
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
    if (m_parents.empty()) {
      // no equalities left: collapse back to top
      set_to_top();
    }
    normalize();
    check_lattice_val();
    CRAB_LOG("symb-var-eq", crab::outs() << "After forget "; dump(););
  }

  /// @brief alternative operation to forget a set of elements
  /// @param elements a set of elements
  void forget(const std::vector<element_t> &elements) override {
    if (is_bottom() || is_top()) {
      return;
    }

    // TODO: if necessary, provide a direct implementation without too much
    // replacing representative
    for (auto v : elements) {
      *this -= v;
    }
  }

  /// @brief keep equivalence classes only for elements in vector elements
  /// @param elements a vector of elements required to keep
  /// @note  After project operation, the state only keeps equivalence classes
  ///        for elements in vector elements
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
      res.m_classes.insert({new_rep, m_classes.at(rep)});
    }
    if (!res.m_parents.empty()) {
      // res was default-constructed as top; mark it as a real state
      res.set_neither_top_or_bottom();
    }
    std::swap(res, *this);
    normalize();
    check_lattice_val();
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
    CRAB_LOG("symb-var-eq", crab::outs() << "Renaming ";
             print_elems_vector(crab::outs(), old_elements);
             crab::outs() << " to ";
             print_elems_vector(crab::outs(), new_elements);
             crab::outs() << "\n"; dump(););

    // Build the substitution old_v |-> new_v for the elements that actually
    // exist. rename() assumes each new element is fresh (does not exist yet).
    parents_map_t subst;
    for (unsigned i = 0, size = old_elements.size(); i < size; ++i) {
      const element_t &old_v = old_elements[i];
      const element_t &new_v = new_elements[i];
      if (!contains(old_v)) {
        continue;
      }
      if (contains(new_v)) {
        CRAB_ERROR(domain_name(), "::rename assumes that ", new_v,
                   " does not exist");
      }
      subst.insert({old_v, new_v});
    }
    if (subst.empty()) {
      return;
    }

    // Apply the substitution to a single element.
    auto rename_elem = [&subst](const element_t &e) -> element_t {
      auto it = subst.find(e);
      return it == subst.end() ? e : it->second;
    };

    // An element appears in m_parents both as a key (every element) and, when
    // it is a representative, as a value (its members point to it). Both sides
    // may be renamed, so rewrite the whole map in one pass.
    parents_map_t new_parents;
    new_parents.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      new_parents.insert({rename_elem(kv.first), rename_elem(kv.second)});
    }
    std::swap(m_parents, new_parents);

    // m_classes is keyed by representative, so only the key may be renamed.
    classes_map_t new_classes;
    new_classes.reserve(m_classes.size());
    for (auto &kv : m_classes) {
      new_classes.insert({rename_elem(kv.first), std::move(kv.second)});
    }
    std::swap(m_classes, new_classes);

    check_lattice_val();
    CRAB_LOG("symb-var-eq", crab::outs() << "After renaming "; dump(););
  }

  // Reduce the size of the abstract domain representation.
  void minimize() override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  /// @brief A normalization function to produce a standard form
  /// @note  Drops singleton classes such as {v3}=>#var5: a class with a single
  /// element conveys no equality, so removing it brings the state closer to top
  /// (and a singleton-only state becomes top). Enabled by
  /// DomainParams::normalize, which is on by default. It can be turned off (see
  /// SVEQNoNormalizeParams) when the domain keeps equalities between fields and
  /// registers and singleton classes must be retained.
  void normalize() override {
    SVEQ_DOMAIN_SCOPED_STATS(".normalize");
    if (!DomainParams::normalize) {
      return;
    }
    if (is_top() || is_bottom()) {
      return;
    }
    // Count members per representative.
    std::unordered_map<element_t, unsigned> class_size;
    for (auto &kv : m_parents) {
      class_size[kv.second]++;
    }
    // Drop singletons (a representative whose only member is itself).
    for (auto it = m_parents.begin(); it != m_parents.end();) {
      if (class_size[it->second] == 1) {
        m_classes.erase(it->second);
        it = m_parents.erase(it);
      } else {
        ++it;
      }
    }
    if (m_parents.empty()) {
      set_to_top();
    }
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
      print_separated(
          crab::outs(), m_parents,
          [&](const std::pair<const element_t, element_t> &kv) {
            crab::outs() << kv.first << " => " << kv.second;
          },
          '{', '}', ", ");
    }
    crab::outs() << "\n";
  }

  std::string domain_name() const override { return "EqDomain"; }

  // Concretize the state as a conjunction of variable equalities: for each
  // equivalence class, every non-representative member m yields m == rep.
  linear_constraint_system_t to_linear_constraint_system() const override {
    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }
    if (is_top()) {
      return csts; // empty conjunction == true
    }
    for (auto &kv : m_parents) {
      const element_t &member = kv.first;
      const element_t &rep = kv.second;
      if (member == rep) {
        continue; // representative / singleton carries no equality
      }
      csts += linear_constraint_t(linear_expression_t(member) -
                                      linear_expression_t(rep),
                                  linear_constraint_t::EQUALITY);
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_domain_t &invariant) override {
    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const this_domain_t &caller) override {
    SVEQ_DOMAIN_SCOPED_STATS(".callee_entry");
    inter_abstract_operations<
        this_domain_t,
        DomainParams::implement_inter_transformers>::callee_entry(callsite,
                                                                  caller,
                                                                  *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const this_domain_t &callee) override {
    SVEQ_DOMAIN_SCOPED_STATS(".caller_cont");
    inter_abstract_operations<this_domain_t,
                              DomainParams::implement_inter_transformers>::
        caller_continuation(callsite, callee, *this);
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

template <typename Number, typename VariableName, typename DomainParams>
struct abstract_domain_traits<
    symbolic_variable_equality_domain<Number, VariableName, DomainParams>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // end namespace domains
} // end namespace crab