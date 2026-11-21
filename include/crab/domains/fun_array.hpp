#pragma once

/**
 * FunArray: array content analysis by segmentation (P. Cousot, R. Cousot,
 * F. Logozzo. A Parametric Segmentation Functor for Fully Automatic and
 * Scalable Array Content Analysis. POPL 2011).
 *
 * Arrays are partitioned into groups, and each group is described by one
 * segmentation
 *
 *   {B0} P0 {B1}[?] P1 ... {Bn}[?] Pn {+oo}
 *
 * whose bounds are sets of equal symbolic expressions (x + k). The value Pk
 * of a segment [Bk, Bk+1) of a group with arrays A, B, ... is a relational
 * state (zones by default) that stands for
 *
 *   forall j in [Bk, Bk+1).  phi(A[j], B[j], ..., j, scalars)
 *
 * with one quantified index j shared by all the arrays of the group. This
 * is the shape of every invariant the domain can express; for instance:
 *
 *   A[j] == 0          every element of the segment is 0
 *   A[j] == j + 2      every element is its index plus 2
 *   A[j] <= m          every element is at most the scalar m
 *   A[j] == B[j]       A and B are equal on the segment
 *
 * Relations between elements at different indices, such as sortedness
 * (A[j1] <= A[j2] for j1 < j2) or A[j1] < B[j2] for all j1, j2, are not
 * expressible; the latter can still be derived when a scalar separates both
 * arrays (A[j] < m and B[j] > m).
 *
 * Every array starts in its own group; two groups merge
 * when a store relates them (A[i] := t where a value of B's group mentions
 * t, e.g. after t := B[i]) or on array_assign. Scalars are kept in
 * ScalarDom.
 *
 * Assumptions:
 *  - Indices count elements, not bytes: A[i+1] is the element after A[i].
 *    An array accessed with an element size different from its first one is
 *    forgotten.
 *  - Indices are non-negative: every load and store at index i adds i >= 0
 *    to the scalar state.
 *  - Only arrays of integers are tracked; other arrays are ignored.
 *  - array_init and array_store_range forget the array's contents.
 **/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_mixins.hpp>
#include <crab/domains/fun_array/expression_set_bound.hpp>
#include <crab/domains/fun_array/offset_expression.hpp>
#include <crab/domains/fun_array/segmentation.hpp>
#include <crab/domains/fun_array/zone_segment_values.hpp>
#include <crab/domains/inter_abstract_operations.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>

#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace crab {
namespace domains {

class FunArrayDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};

namespace fun_array_impl {
// union_find_domain hooks: an array leaving or renamed in a group is
// forgotten or renamed in the group's segmentation
template <class Seg> struct forget_array_in_group {
  void operator()(Seg &seg, const typename Seg::variable_t &a) const {
    seg.forget_array(a);
  }
};
template <class Seg> struct rename_array_in_group {
  void operator()(Seg &seg, const typename Seg::variable_t &a,
                  const typename Seg::variable_t &b) const {
    seg.rename({a}, {b});
  }
};
} // end namespace fun_array_impl

template <class ScalarDom,
          class BoundDom = fun_array_impl::expression_set_bound<
              fun_array_impl::offset_expression<ScalarDom>, ScalarDom>,
          class ValueDom = fun_array_impl::zone_segment_values<ScalarDom>,
          class Params = FunArrayDefaultParams>
class fun_array final
    : public abstract_domain_base<
          fun_array<ScalarDom, BoundDom, ValueDom, Params>,
          region_and_reference_operations_not_implemented,
          default_inter_operations> {
public:
  using number_t = typename ScalarDom::number_t;
  using varname_t = typename ScalarDom::varname_t;

private:
  using fun_array_t = fun_array<ScalarDom, BoundDom, ValueDom, Params>;
  using abstract_domain_t = abstract_domain_api<fun_array_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using interval_t = ikos::interval<number_t>;
  using segmentation_t = fun_array_impl::segmentation<BoundDom, ValueDom>;
  using expression_t = typename BoundDom::expression_t;

  static_assert(std::is_same<decltype(ValueDom::make_top_for(
                                 std::declval<const variable_t &>())),
                             ValueDom>::value,
                "ValueDom must have static ValueDom "
                "make_top_for(const variable_t &)");

  // read by the default_inter_operations mixin
  using params_t = Params;

private:
  // The groups: arrays are the elements, a group's segmentation is its
  // class state. Merging two groups meets their segmentations: they
  // describe different arrays, so the group's state is their conjunction.
  // Arrays that are not in m_groups are {0} T {+oo}.
  using groups_t = union_find_domain<
      variable_t, segmentation_t, uf_intersection_semantics<segmentation_t>,
      fun_array_impl::forget_array_in_group<segmentation_t>,
      fun_array_impl::rename_array_in_group<segmentation_t>>;
  using elem_size_map_t = std::map<variable_t, uint64_t>;

  ScalarDom m_scalar;
  groups_t m_groups;
  // bytes per element of each array in m_groups, from its first access
  elem_size_map_t m_elem_size;

  fun_array(ScalarDom scalar, groups_t groups, elem_size_map_t elem_size)
      : m_scalar(std::move(scalar)), m_groups(std::move(groups)),
        m_elem_size(std::move(elem_size)) {
    normalize_arrays();
  }

  static segmentation_t top_segmentation(const variable_t &a) {
    return segmentation_t(ValueDom::make_top_for(a));
  }

  static bool is_tracked_array(const variable_t &a) {
    auto ty = a.get_type();
    return ty.is_integer_array() || ty.is_math_integer_array();
  }

  // The arrays in m_groups
  std::vector<variable_t> arrays() const {
    std::vector<variable_t> res;
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      res.insert(res.end(), kv.second.begin(), kv.second.end());
    }
    return res;
  }

  // A bottom segmentation makes the whole state bottom. A group of one
  // array whose segmentation is top is removed: it is the same as the array
  // not being there. Groups of two or more arrays are kept even when top,
  // so that groups only grow (termination).
  void normalize_arrays() {
    if (m_scalar.is_bottom() || m_groups.is_bottom()) {
      set_to_bottom();
      return;
    }
    std::vector<variable_t> drop;
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      const segmentation_t &seg = shared_group_of(kv.first);
      if (seg.is_bottom()) {
        set_to_bottom();
        return;
      }
      if (kv.second.size() == 1 && seg.is_top()) {
        drop.push_back(kv.second[0]);
      }
    }
    for (auto const &a : drop) {
      m_groups.forget(a);
      m_elem_size.erase(a);
    }
  }

  // Convert elem_size to an uint64_t (as array_smashing does)
  uint64_t check_and_get_elem_size(const linear_expression_t &elem_size) const {
    interval_t r(elem_size.constant());
    for (auto const &kv : elem_size) {
      r += interval_t(kv.first) * m_scalar.at(kv.second);
    }
    if (boost::optional<number_t> n_bytes = r.singleton()) {
      if (static_cast<int64_t>(*n_bytes) > 0) {
        return (uint64_t) static_cast<int64_t>(*n_bytes);
      }
    }
    CRAB_ERROR("array domains expect constant array element sizes ",
               "between 1 and ", std::numeric_limits<uint64_t>::max(),
               ". Found ", elem_size);
  }

  // Make sure a is in m_groups for an access of size bytes: an array not
  // accessed before gets its own group {0} T {+oo}; an array accessed with
  // another size is forgotten and gets a new group. An array in a group but
  // with no size yet (from array_assign) takes this size.
  void prepare_access(const variable_t &a, uint64_t size) {
    auto known = elem_size_of(m_elem_size, a);
    if (m_groups.contains(a) && (!known || *known == size)) {
      m_elem_size[a] = size;
      return;
    }
    if (m_groups.contains(a)) {
      m_groups.forget(a);
    }
    m_groups.set(a, top_segmentation(a));
    m_elem_size[a] = size;
  }

  // The segmentation of a's group, to modify (copied first if shared).
  // Pre-condition: m_groups.contains(a)
  segmentation_t &group_of(const variable_t &a) { return *m_groups.get(a); }

  // The segmentation of a's group, to read (never copied).
  // Pre-condition: m_groups.contains(a)
  const segmentation_t &shared_group_of(const variable_t &a) const {
    return *m_groups.get(a);
  }

  // The element size of array a, if known
  static boost::optional<uint64_t> elem_size_of(const elem_size_map_t &sizes,
                                                const variable_t &a) {
    auto it = sizes.find(a);
    if (it == sizes.end()) {
      return boost::none;
    }
    return it->second;
  }

  // Apply f to every group's segmentation that mentions x, leaving the
  // others shared
  template <class F> void for_each_mentioning(const variable_t &x, F f) {
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      if (shared_group_of(kv.first).mentions(x)) {
        f(group_of(kv.first));
      }
    }
  }

  // Apply f to every group's segmentation
  template <class F> void for_each_group(F f) {
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      f(group_of(kv.first));
    }
  }

  // x is assigned: drop it from every segmentation
  void forget_in_arrays(const variable_t &x) {
    for_each_mentioning(x, [&x](segmentation_t &seg) { seg.forget(x); });
  }

  // Indices are non-negative
  void assume_non_negative(const linear_expression_t &i) {
    m_scalar += linear_constraint_t(-i, linear_constraint_t::INEQUALITY);
  }

  // x := e on the segmentations, before m_scalar sees it: bounds are
  // rewritten and the assignment is replayed in the values that mention x
  // (segmentation::assign)
  void assign_in_arrays(const variable_t &x, const linear_expression_t &e) {
    if (auto ee = expression_t::make(e)) {
      // x may be added to bounds containing e: every group is concerned
      for_each_group([&](segmentation_t &seg) { seg.assign(x, e); });
    } else {
      for_each_mentioning(x,
                          [&](segmentation_t &seg) { seg.assign(x, e); });
    }
  }

  // The bound equalities and orderings of every segmentation, as constraints
  // over scalars (to_linear_constraint_system)
  linear_constraint_system_t bound_constraints() const {
    linear_constraint_system_t res;
    if (is_bottom()) {
      return res;
    }
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      const segmentation_t &seg = shared_group_of(kv.first);
      res += seg.template bound_equalities<linear_constraint_system_t>();
      res += seg.template bound_orderings<linear_constraint_system_t>();
    }
    return res;
  }

  // Reduce the scalar state from the bounds of every segmentation: their
  // equalities go to the scalar state (e.g., i == n for a bound {i n}). The
  // converse of segmentation::reduce_bounds_from_scalar.
  //
  // Not called on the result of widening. In X(k+1) = X(k) || Y(k),
  // widening guarantees that the iterates stabilize for any right operands
  // Y(k), which may be reduced (they come from join), but only if each left
  // operand is the previous result unchanged. Widening drops information to
  // stabilize; reducing its result could give some of it back, and the
  // iterates might never stabilize.
  void reduce_scalar_from_bounds() {
    if (m_scalar.is_bottom() || m_groups.is_bottom()) {
      return;
    }
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      m_scalar += shared_group_of(kv.first)
                      .template bound_equalities<linear_constraint_system_t>();
    }
  }

  // Merge the groups of a and b (both in m_groups)
  void merge_groups(const variable_t &a, const variable_t &b) {
    m_groups.join(a, b);
    if (m_groups.is_bottom()) {
      set_to_bottom();
    }
  }

  // Prepare two states for a binary lattice operation:
  //  - arrays with different element sizes on each side are forgotten on
  //    both;
  //  - when saturate is true (join, <=), every array's group segmentation on
  //    each side gains the other side's bound expressions that its own
  //    scalar state proves equal to its bounds (segmentation::saturate); an
  //    array on one side only is given {0} T {+oo} on the other first, so
  //    that it can gain them too.
  // The union-find never sees the scalar states, hence this pre-pass.
  static void prepare(fun_array_t &l, fun_array_t &r, bool saturate) {
    for (auto const &a : l.arrays()) {
      auto ls = elem_size_of(l.m_elem_size, a);
      auto rs = elem_size_of(r.m_elem_size, a);
      if (ls && rs && *ls != *rs) {
        l.m_groups.forget(a);
        l.m_elem_size.erase(a);
        r.m_groups.forget(a);
        r.m_elem_size.erase(a);
      }
    }
    if (!saturate) {
      return;
    }
    auto materialize = [](fun_array_t &to, const fun_array_t &from) {
      for (auto const &a : from.arrays()) {
        if (!to.m_groups.contains(a)) {
          to.m_groups.set(a, top_segmentation(a));
          if (auto s = elem_size_of(from.m_elem_size, a)) {
            to.m_elem_size[a] = *s;
          }
        }
      }
    };
    materialize(r, l);
    materialize(l, r);
    for (auto const &a : l.arrays()) {
      segmentation_t ls = l.group_of(a);
      l.group_of(a).saturate(r.group_of(a), l.m_scalar);
      r.group_of(a).saturate(ls, r.m_scalar);
    }
  }

  // The element sizes of the arrays in groups, taken from l or r
  static elem_size_map_t elem_sizes(const groups_t &groups,
                                    const fun_array_t &l,
                                    const fun_array_t &r) {
    elem_size_map_t res;
    if (groups.is_bottom()) {
      return res;
    }
    for (auto const &kv : groups.equiv_classes_elems()) {
      for (auto const &a : kv.second) {
        if (auto s = elem_size_of(l.m_elem_size, a)) {
          res[a] = *s;
        } else if (auto s = elem_size_of(r.m_elem_size, a)) {
          res[a] = *s;
        }
      }
    }
    return res;
  }

public:
  fun_array() {}

  fun_array(const fun_array_t &o) = default;
  fun_array(fun_array_t &&o) = default;
  fun_array_t &operator=(const fun_array_t &o) = default;
  fun_array_t &operator=(fun_array_t &&o) = default;

  bool is_asc_phase() const override { return m_scalar.is_asc_phase(); }
  void set_phase(bool is_ascending) override {
    m_scalar.set_phase(is_ascending);
  }

  fun_array_t make_top() const override { return fun_array_t(); }

  fun_array_t make_bottom() const override {
    fun_array_t res;
    res.set_to_bottom();
    return res;
  }

  void set_to_top() override {
    m_scalar.set_to_top();
    m_groups = groups_t();
    m_elem_size.clear();
  }

  void set_to_bottom() override {
    m_scalar.set_to_bottom();
    m_groups = groups_t();
    m_elem_size.clear();
  }

  bool is_bottom() const override { return m_scalar.is_bottom(); }

  bool is_top() const override {
    return m_scalar.is_top() && m_groups.is_empty();
  }

  // As for join, the segmentations are first saturated (prepare): e.g., with
  // i == 0, an untouched {0} T {+oo} is {0 i} T {+oo}, which is included in
  // {0} 0 {i}? T {+oo}. Then union_find_domain::operator<=, which is false
  // unless this side's groups refine o's.
  bool operator<=(const fun_array_t &o) const override {
    if (is_bottom()) {
      return true;
    }
    if (o.is_bottom()) {
      return false;
    }
    if (!(m_scalar <= o.m_scalar)) {
      return false;
    }
    for (auto const &a : o.arrays()) {
      auto ls = elem_size_of(m_elem_size, a);
      auto os = elem_size_of(o.m_elem_size, a);
      if (ls && os && *ls != *os) {
        return false;
      }
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, true);
    return l.m_groups <= o.m_groups;
  }

  void operator|=(const fun_array_t &o) override { *this = *this | o; }

  // Join: the segmentations are saturated (prepare), then the union-find
  // aligns the groups of both sides (merging them by meet) and joins the
  // segmentations group by group. Arrays on one side only are forgotten.
  fun_array_t operator|(const fun_array_t &o) const override {
    if (is_bottom()) {
      return o;
    }
    if (o.is_bottom()) {
      return *this;
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, true);
    groups_t groups = l.m_groups | r.m_groups;
    elem_size_map_t sizes = elem_sizes(groups, l, r);
    fun_array_t res(m_scalar | o.m_scalar, std::move(groups),
                    std::move(sizes));
    res.reduce_scalar_from_bounds();
    return res;
  }

  void operator&=(const fun_array_t &o) override { *this = *this & o; }

  fun_array_t operator&(const fun_array_t &o) const override {
    if (is_bottom() || o.is_bottom()) {
      return make_bottom();
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, false);
    groups_t groups = l.m_groups & r.m_groups;
    elem_size_map_t sizes = elem_sizes(groups, l, r);
    fun_array_t res(m_scalar & o.m_scalar, std::move(groups),
                    std::move(sizes));
    res.reduce_scalar_from_bounds();
    return res;
  }

  // Widening; no saturation, so that bounds can only lose expressions
  fun_array_t operator||(const fun_array_t &o) const override {
    if (is_bottom()) {
      return o;
    }
    if (o.is_bottom()) {
      return *this;
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, false);
    groups_t groups = l.m_groups || r.m_groups;
    elem_size_map_t sizes = elem_sizes(groups, l, r);
    return fun_array_t(m_scalar || o.m_scalar, std::move(groups),
                       std::move(sizes));
  }

  fun_array_t
  widening_thresholds(const fun_array_t &o,
                      const thresholds<number_t> &ts) const override {
    if (is_bottom()) {
      return o;
    }
    if (o.is_bottom()) {
      return *this;
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, false);
    groups_t groups = l.m_groups || r.m_groups;
    elem_size_map_t sizes = elem_sizes(groups, l, r);
    return fun_array_t(m_scalar.widening_thresholds(o.m_scalar, ts),
                       std::move(groups), std::move(sizes));
  }

  // Narrowing: the segmentations are met; the fixpoint runs a bounded
  // number of descending iterations
  fun_array_t operator&&(const fun_array_t &o) const override {
    if (is_bottom() || o.is_bottom()) {
      return make_bottom();
    }
    fun_array_t l(*this), r(o);
    prepare(l, r, false);
    groups_t groups = l.m_groups & r.m_groups;
    elem_size_map_t sizes = elem_sizes(groups, l, r);
    return fun_array_t(m_scalar && o.m_scalar, std::move(groups),
                       std::move(sizes));
  }

  interval_t operator[](const variable_t &v) override { return m_scalar[v]; }

  interval_t at(const variable_t &v) const override { return m_scalar.at(v); }

  // Forgetting an array removes it from its group: it is forgotten in the
  // group's values (segmentation::forget_array) and the other arrays keep
  // their facts.
  void operator-=(const variable_t &v) override {
    if (v.get_type().is_array()) {
      m_groups.forget(v);
      m_elem_size.erase(v);
    } else {
      forget_in_arrays(v);
      m_scalar -= v;
    }
  }

  void forget(const variable_vector_t &vars) override {
    for (auto const &v : vars) {
      *this -= v;
    }
  }

  void project(const variable_vector_t &vars) override {
    variable_vector_t scalars, arrays;
    for (auto const &v : vars) {
      if (v.get_type().is_array()) {
        arrays.push_back(v);
      } else {
        scalars.push_back(v);
      }
    }
    m_groups.project(arrays);
    for (auto it = m_elem_size.begin(); it != m_elem_size.end();) {
      if (!m_groups.contains(it->first)) {
        it = m_elem_size.erase(it);
      } else {
        ++it;
      }
    }
    for_each_group([&](segmentation_t &seg) { seg.project(scalars); });
    m_scalar.project(scalars);
    normalize_arrays();
  }

  // An expanded array gets its own group: a copy of x's group without the
  // other arrays
  void expand(const variable_t &x, const variable_t &new_x) override {
    if (x.get_type() != new_x.get_type()) {
      CRAB_ERROR(domain_name(), "::expand must preserve the same type");
    }
    if (x.get_type().is_array()) {
      if (m_groups.contains(x)) {
        segmentation_t seg = group_of(x);
        for (auto const &kv : m_groups.equiv_classes_elems()) {
          for (auto const &b : kv.second) {
            if (!(b == x) && m_groups.find(b) == m_groups.find(x)) {
              seg.forget_array(b);
            }
          }
        }
        seg.rename({x}, {new_x});
        *this -= new_x;
        m_groups.set(new_x, seg);
        if (auto s = elem_size_of(m_elem_size, x)) {
          m_elem_size[new_x] = *s;
        }
      }
    } else {
      m_scalar.expand(x, new_x);
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename expects vectors same sizes");
    }
    variable_vector_t scalars_from, scalars_to, arrays_from, arrays_to;
    for (unsigned k = 0, n = from.size(); k < n; ++k) {
      if (from[k].get_type() != to[k].get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve the same type");
      }
      if (from[k].get_type().is_array()) {
        arrays_from.push_back(from[k]);
        arrays_to.push_back(to[k]);
      } else {
        scalars_from.push_back(from[k]);
        scalars_to.push_back(to[k]);
      }
    }
    // arrays: membership, element sizes and their placeholders in values
    m_groups.rename(arrays_from, arrays_to);
    elem_size_map_t sizes;
    for (auto const &kv : m_elem_size) {
      variable_t name = kv.first;
      for (unsigned k = 0, n = arrays_from.size(); k < n; ++k) {
        if (arrays_from[k] == kv.first) {
          name = arrays_to[k];
        }
      }
      sizes[name] = kv.second;
    }
    m_elem_size = std::move(sizes);
    // scalars: bounds and values
    if (!scalars_from.empty()) {
      for_each_group([&](segmentation_t &seg) {
        seg.rename(scalars_from, scalars_to);
      });
    }
    m_scalar.rename(scalars_from, scalars_to);
  }

  void normalize() override { m_scalar.normalize(); }
  void minimize() override { m_scalar.minimize(); }

  // A guard refines the scalar state, is replayed in the segment values that
  // mention its variables, and reduces the bounds (segmentation::assume,
  // segmentation::reduce_bounds_from_scalar). E.g., in the else branch of
  // `t := A[i]; if (t > m) ..`, t <= m turns A==t on [i, i+1) into A<=m.
  void operator+=(const linear_constraint_system_t &csts) override {
    m_scalar += csts;
    if (m_scalar.is_bottom()) {
      set_to_bottom();
      return;
    }
    for_each_group([&](segmentation_t &seg) {
      seg.assume(csts);
      seg.reduce_bounds_from_scalar(m_scalar);
    });
    normalize_arrays();
    reduce_scalar_from_bounds();
  }

  bool entails(const linear_constraint_t &cst) const override {
    return m_scalar.entails(cst);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    assign_in_arrays(x, e);
    m_scalar.assign(x, e);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    forget_in_arrays(x);
    m_scalar.weak_assign(x, e);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    // operations with a linear result are replayed as assignments
    if (op == OP_ADDITION) {
      assign_in_arrays(x, linear_expression_t(y) + z);
    } else if (op == OP_SUBTRACTION) {
      assign_in_arrays(x, linear_expression_t(y) - z);
    } else if (op == OP_MULTIPLICATION) {
      assign_in_arrays(x, linear_expression_t(z, y));
    } else {
      forget_in_arrays(x);
    }
    m_scalar.apply(op, x, y, z);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (op == OP_ADDITION) {
      assign_in_arrays(x, linear_expression_t(y) + z);
    } else if (op == OP_SUBTRACTION) {
      assign_in_arrays(x, linear_expression_t(y) - z);
    } else {
      forget_in_arrays(x);
    }
    m_scalar.apply(op, x, y, z);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    forget_in_arrays(lhs);
    m_scalar.select(lhs, cond, e1, e2);
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    forget_in_arrays(dst);
    m_scalar.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    forget_in_arrays(x);
    m_scalar.apply(op, x, y, z);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    forget_in_arrays(x);
    m_scalar.apply(op, x, y, k);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const fun_array_t &inv) override {
    CRAB_WARN("backward_assign in fun_array domain not implemented");
    *this -= x;
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const fun_array_t &inv) override {
    CRAB_WARN("backward_apply in fun_array domain not implemented");
    *this -= x;
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const fun_array_t &inv) override {
    CRAB_WARN("backward_apply in fun_array domain not implemented");
    *this -= x;
  }

  // boolean operations: booleans never appear in bounds, but segment
  // values may mention them
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    forget_in_arrays(lhs);
    m_scalar.assign_bool_cst(lhs, rhs);
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    forget_in_arrays(lhs);
    m_scalar.assign_bool_ref_cst(lhs, rhs);
  }

  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    forget_in_arrays(lhs);
    m_scalar.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {
    forget_in_arrays(lhs);
    m_scalar.weak_assign_bool_cst(lhs, rhs);
  }

  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {
    forget_in_arrays(lhs);
    m_scalar.weak_assign_bool_var(lhs, rhs, is_not_rhs);
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    forget_in_arrays(x);
    m_scalar.apply_binary_bool(op, x, y, z);
  }

  void assume_bool(const variable_t &v, bool is_negated) override {
    m_scalar.assume_bool(v, is_negated);
    normalize_arrays();
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    forget_in_arrays(lhs);
    m_scalar.select_bool(lhs, cond, b1, b2);
  }

  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const fun_array_t &inv) override {
    CRAB_WARN("backward_assign_bool_cst in fun_array domain not implemented");
    *this -= lhs;
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const fun_array_t &inv) override {
    CRAB_WARN(
        "backward_assign_bool_ref_cst in fun_array domain not implemented");
    *this -= lhs;
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const fun_array_t &inv) override {
    CRAB_WARN("backward_assign_bool_var in fun_array domain not implemented");
    *this -= lhs;
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const fun_array_t &inv) override {
    CRAB_WARN(
        "backward_apply_binary_bool in fun_array domain not implemented");
    *this -= x;
  }

  // array operations

  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    *this -= a;
  }

  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    if (is_bottom()) {
      return;
    }
    uint64_t size = check_and_get_elem_size(elem_size);
    if (!is_tracked_array(a)) {
      *this -= lhs;
      return;
    }
    assume_non_negative(i);
    if (m_scalar.is_bottom()) {
      set_to_bottom();
      return;
    }
    prepare_access(a, size);
    if (i[lhs] != number_t(0)) {
      // lhs := A[lhs + ..]: the index needs the old lhs, so the element is
      // loaded into a fresh temporary first, then assigned to lhs
      auto &vfac = const_cast<varname_t *>(&(lhs.name()))->get_var_factory();
      variable_t tmp(vfac.make_temporary_varname(), lhs.get_type());
      group_of(a).load(a, tmp, i, m_scalar);
      normalize_arrays();
      if (is_bottom()) {
        return;
      }
      assign(lhs, linear_expression_t(tmp));
      *this -= tmp;
    } else {
      *this -= lhs;
      group_of(a).load(a, lhs, i, m_scalar);
    }
    normalize_arrays();
    reduce_scalar_from_bounds();
  }

  // A[i] := val. If a value of another group mentions a variable of val
  // (e.g. after t := B[i], B's group has B == t), the two groups merge first,
  // so that the store relates A and B (A == B on [i, i+1)).
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i,
                   const linear_expression_t &val,
                   bool is_strong_update) override {
    if (is_bottom()) {
      return;
    }
    uint64_t size = check_and_get_elem_size(elem_size);
    if (!is_tracked_array(a)) {
      return;
    }
    assume_non_negative(i);
    if (m_scalar.is_bottom()) {
      set_to_bottom();
      return;
    }
    prepare_access(a, size);
    for (auto const &v : val.variables()) {
      for (auto const &kv : m_groups.equiv_classes_elems()) {
        if (m_groups.find(kv.first) != m_groups.find(a) &&
            shared_group_of(kv.first).mentions(v)) {
          merge_groups(a, kv.first);
          if (is_bottom()) {
            return;
          }
        }
      }
    }
    // is_strong_update is the client's guarantee that the store may be
    // modelled as strong even when the domain cannot tell (array_smashing
    // needs it). It is safe to ignore here: whether the update is strong
    // depends on the index (see segmentation::store).
    group_of(a).store(a, i, val, m_scalar);
    normalize_arrays();
    reduce_scalar_from_bounds();
  }

  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {
    *this -= a;
  }

  // The whole array lhs becomes rhs: lhs joins rhs's group and every value
  // of the group gets lhs == rhs
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    if (is_bottom() || lhs == rhs) {
      return;
    }
    *this -= lhs;
    if (!m_groups.contains(rhs)) {
      m_groups.set(rhs, top_segmentation(rhs));
    }
    m_groups.add(rhs, lhs);
    group_of(rhs).assign_array(lhs, rhs);
    auto it = m_elem_size.find(rhs);
    if (it != m_elem_size.end()) {
      m_elem_size[lhs] = it->second;
    }
    normalize_arrays();
  }

  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const fun_array_t &invariant) override {
    CRAB_WARN("backward_array_init in fun_array domain not implemented");
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const fun_array_t &invariant) override {
    CRAB_WARN("backward_array_load in fun_array domain not implemented");
    *this -= lhs;
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const fun_array_t &invariant) override {
    CRAB_WARN("backward_array_store in fun_array domain not implemented");
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const fun_array_t &invariant) override {
    CRAB_WARN(
        "backward_array_store_range in fun_array domain not implemented");
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const fun_array_t &invariant) override {
    CRAB_WARN("backward_array_assign in fun_array domain not implemented");
  }

  // The scalar state, plus what the bounds of every segmentation say about
  // scalars: equalities within a bound ({i n} gives i == n) and orderings
  // between consecutive bounds ({i} .. {n} gives i < n, or i <= n with a ?).
  // Segment values are not exported: they are universally quantified over
  // array elements, which linear constraints over scalars cannot express.
  linear_constraint_system_t to_linear_constraint_system() const override {
    linear_constraint_system_t res = m_scalar.to_linear_constraint_system();
    res += bound_constraints();
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto disj_csts = m_scalar.to_disjunctive_linear_constraint_system();
    if (disj_csts.is_false()) {
      return disj_csts;
    }
    linear_constraint_system_t bounds = bound_constraints();
    disjunctive_linear_constraint_system_t res;
    for (auto csts : disj_csts) {
      csts += bounds;
      res += csts;
    }
    return res;
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    for (auto const &v : outputs) {
      forget_in_arrays(v);
    }
    m_scalar.intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const fun_array_t &invariant) override {
    CRAB_WARN("backward_intrinsic in fun_array domain not implemented");
    forget(outputs);
  }

  // The segmentation of a's group, for tests and clients: {0} T {+oo} if a
  // was not accessed
  segmentation_t segmentation_of(const variable_t &a) const {
    if (is_bottom() || !m_groups.contains(a)) {
      return top_segmentation(a);
    }
    return *m_groups.get(a);
  }

  // Are a and b in the same group?
  bool same_group(const variable_t &a, const variable_t &b) const {
    return !is_bottom() && m_groups.contains(a) && m_groups.contains(b) &&
           m_groups.find(a) == m_groups.find(b);
  }

  // scalar state, then one line per group: {A, B}: {0} A==B {i}? T {+oo}
  void write(crab_os &o) const override {
    o << m_scalar;
    if (is_bottom()) {
      return;
    }
    for (auto const &kv : m_groups.equiv_classes_elems()) {
      o << "\n{";
      for (unsigned k = 0; k < kv.second.size(); ++k) {
        o << (k > 0 ? ", " : "") << kv.second[k];
      }
      o << "}: " << *m_groups.get(kv.first);
    }
  }

  std::string domain_name() const override { return "FunArray"; }
};

template <class ScalarDom, class BoundDom, class ValueDom, class Params>
struct abstract_domain_traits<fun_array<ScalarDom, BoundDom, ValueDom, Params>> {
  using number_t = typename ScalarDom::number_t;
  using varname_t = typename ScalarDom::varname_t;
};

} // end namespace domains
} // end namespace crab
