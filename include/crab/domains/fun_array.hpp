#pragma once

/**
 * FunArray: array content analysis by segmentation (P. Cousot, R. Cousot,
 * F. Logozzo. A Parametric Segmentation Functor for Fully Automatic and
 * Scalable Array Content Analysis. POPL 2011).
 *
 * Each array is described by a segmentation
 *
 *   {B0} P0 {B1}[?] P1 ... {Bn}[?] Pn {+oo}
 *
 * whose bounds are sets of equal symbolic expressions (x + k) and whose
 * values describe every element of their segment. Scalars are kept in
 * ScalarDom. Every array currently has its own segmentation.
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
#include <crab/support/debug.hpp>

#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace crab {
namespace domains {

class FunArrayDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};

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
  struct array_state {
    segmentation_t seg;
    // bytes per element, from the first access
    uint64_t elem_size;
  };
  // Arrays that are not in m_arrays are {0} T {+oo}
  using array_map_t = std::map<variable_t, array_state>;

  ScalarDom m_scalar;
  array_map_t m_arrays;

  fun_array(ScalarDom scalar, array_map_t arrays)
      : m_scalar(std::move(scalar)), m_arrays(std::move(arrays)) {
    normalize_arrays();
  }

  static segmentation_t top_segmentation(const variable_t &a) {
    return segmentation_t(ValueDom::make_top_for(a));
  }

  static bool is_tracked_array(const variable_t &a) {
    auto ty = a.get_type();
    return ty.is_integer_array() || ty.is_math_integer_array();
  }

  // A bottom segmentation makes the whole state bottom; a top one is
  // removed from the map
  void normalize_arrays() {
    if (m_scalar.is_bottom()) {
      m_arrays.clear();
      return;
    }
    for (auto it = m_arrays.begin(); it != m_arrays.end();) {
      if (it->second.seg.is_bottom()) {
        set_to_bottom();
        return;
      }
      if (it->second.seg.is_top()) {
        it = m_arrays.erase(it);
      } else {
        ++it;
      }
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

  // The state of array a for an access of size bytes: a's segmentation, or
  // {0} T {+oo} if a was not accessed before or with another size
  array_state &array_for_access(const variable_t &a, uint64_t size) {
    auto it = m_arrays.find(a);
    if (it != m_arrays.end() && it->second.elem_size == size) {
      return it->second;
    }
    if (it != m_arrays.end()) {
      m_arrays.erase(it);
    }
    return m_arrays.emplace(a, array_state{top_segmentation(a), size})
        .first->second;
  }

  // x is assigned: drop it from every segmentation
  void forget_in_arrays(const variable_t &x) {
    for (auto &kv : m_arrays) {
      kv.second.seg.forget(x);
    }
  }

  // Indices are non-negative
  void assume_non_negative(const linear_expression_t &i) {
    m_scalar += linear_constraint_t(-i, linear_constraint_t::INEQUALITY);
  }

  // x := e on the segmentations, before m_scalar sees it: bounds are
  // rewritten and the assignment is replayed in the values that mention x
  // (segmentation::assign)
  void assign_in_arrays(const variable_t &x, const linear_expression_t &e) {
    for (auto &kv : m_arrays) {
      kv.second.seg.assign(x, e);
    }
  }

  // Reduce the scalar state from the bounds of every segmentation: their
  // equalities go to the scalar state (e.g., i == n for a bound {i n}). The
  // converse of segmentation::reduce_bounds_from_scalar. Not after widening,
  // so that a widened scalar state is not refined back.
  void reduce_scalar_from_bounds() {
    if (m_scalar.is_bottom()) {
      return;
    }
    for (auto const &kv : m_arrays) {
      m_scalar += kv.second.seg.template bound_equalities<
          linear_constraint_system_t>();
    }
  }

  // Combine the arrays of two states. Arrays missing on one side are
  // {0} T {+oo} there; op gets both segmentations, each saturated with the
  // other when saturate is true.
  template <class Op>
  static array_map_t combine_arrays(const fun_array_t &l, const fun_array_t &r,
                                    bool saturate, Op op) {
    std::set<variable_t> names;
    for (auto const &kv : l.m_arrays) {
      names.insert(kv.first);
    }
    for (auto const &kv : r.m_arrays) {
      names.insert(kv.first);
    }
    array_map_t res;
    for (auto const &a : names) {
      auto li = l.m_arrays.find(a);
      auto ri = r.m_arrays.find(a);
      if (li != l.m_arrays.end() && ri != r.m_arrays.end() &&
          li->second.elem_size != ri->second.elem_size) {
        continue;
      }
      uint64_t size =
          li != l.m_arrays.end() ? li->second.elem_size : ri->second.elem_size;
      segmentation_t ls =
          li != l.m_arrays.end() ? li->second.seg : top_segmentation(a);
      segmentation_t rs =
          ri != r.m_arrays.end() ? ri->second.seg : top_segmentation(a);
      if (saturate) {
        segmentation_t ls0 = ls;
        ls.saturate(rs, l.m_scalar);
        rs.saturate(ls0, r.m_scalar);
      }
      res.emplace(a, array_state{op(ls, rs), size});
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
    m_arrays.clear();
  }

  void set_to_bottom() override {
    m_scalar.set_to_bottom();
    m_arrays.clear();
  }

  bool is_bottom() const override { return m_scalar.is_bottom(); }

  bool is_top() const override { return m_scalar.is_top() && m_arrays.empty(); }

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
    // As for join, this side's segmentation first gains o's bound
    // expressions that this side's scalar state proves equal to its bounds:
    // e.g., with i == 0, an untouched {0} T {+oo} is {0 i} T {+oo}, which is
    // included in {0} 0 {i}? T {+oo}.
    for (auto const &kv : o.m_arrays) {
      auto it = m_arrays.find(kv.first);
      if (it != m_arrays.end() && it->second.elem_size != kv.second.elem_size) {
        return false;
      }
      segmentation_t seg =
          it != m_arrays.end() ? it->second.seg : top_segmentation(kv.first);
      seg.saturate(kv.second.seg, m_scalar);
      if (!(seg <= kv.second.seg)) {
        return false;
      }
    }
    return true;
  }

  void operator|=(const fun_array_t &o) override { *this = *this | o; }

  // Join; before joining two segmentations, each gets the other's bound
  // expressions that its own scalar state proves equal to its bounds
  // (segmentation::saturate)
  fun_array_t operator|(const fun_array_t &o) const override {
    if (is_bottom()) {
      return o;
    }
    if (o.is_bottom()) {
      return *this;
    }
    fun_array_t res(
        m_scalar | o.m_scalar,
        combine_arrays(*this, o, true,
                       [](const segmentation_t &x, const segmentation_t &y) {
                         return x | y;
                       }));
    res.reduce_scalar_from_bounds();
    return res;
  }

  void operator&=(const fun_array_t &o) override { *this = *this & o; }

  fun_array_t operator&(const fun_array_t &o) const override {
    if (is_bottom() || o.is_bottom()) {
      return make_bottom();
    }
    fun_array_t res(
        m_scalar & o.m_scalar,
        combine_arrays(*this, o, false,
                       [](const segmentation_t &x, const segmentation_t &y) {
                         return x & y;
                       }));
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
    return fun_array_t(
        m_scalar || o.m_scalar,
        combine_arrays(*this, o, false,
                       [](const segmentation_t &x, const segmentation_t &y) {
                         return x || y;
                       }));
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
    return fun_array_t(
        m_scalar.widening_thresholds(o.m_scalar, ts),
        combine_arrays(*this, o, false,
                       [](const segmentation_t &x, const segmentation_t &y) {
                         return x || y;
                       }));
  }

  // Narrowing: the segmentations are met; the fixpoint runs a bounded
  // number of descending iterations
  fun_array_t operator&&(const fun_array_t &o) const override {
    if (is_bottom() || o.is_bottom()) {
      return make_bottom();
    }
    return fun_array_t(
        m_scalar && o.m_scalar,
        combine_arrays(*this, o, false,
                       [](const segmentation_t &x, const segmentation_t &y) {
                         return x & y;
                       }));
  }

  interval_t operator[](const variable_t &v) override { return m_scalar[v]; }

  interval_t at(const variable_t &v) const override { return m_scalar.at(v); }

  void operator-=(const variable_t &v) override {
    if (v.get_type().is_array()) {
      // Correct only while every array has its own segmentation: an array
      // variable then appears only in the values of its own segmentation,
      // so erasing it removes every mention. Once arrays share a
      // segmentation, this must instead forget the array in its group's
      // values (segment_value_api::forget_array) and keep the group for
      // the other arrays.
      m_arrays.erase(v);
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

  // As operator-=, dropping the arrays not in vars by erasing their
  // segmentations is correct only while every array has its own
  // segmentation.
  void project(const variable_vector_t &vars) override {
    variable_vector_t scalars;
    std::set<variable_t> arrays;
    for (auto const &v : vars) {
      if (v.get_type().is_array()) {
        arrays.insert(v);
      } else {
        scalars.push_back(v);
      }
    }
    for (auto it = m_arrays.begin(); it != m_arrays.end();) {
      if (arrays.count(it->first) == 0) {
        it = m_arrays.erase(it);
      } else {
        it->second.seg.project(scalars);
        ++it;
      }
    }
    m_scalar.project(scalars);
    normalize_arrays();
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    if (x.get_type() != new_x.get_type()) {
      CRAB_ERROR(domain_name(), "::expand must preserve the same type");
    }
    if (x.get_type().is_array()) {
      auto it = m_arrays.find(x);
      if (it != m_arrays.end()) {
        array_state st = it->second;
        st.seg.rename({x}, {new_x});
        m_arrays.erase(new_x);
        m_arrays.emplace(new_x, std::move(st));
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
    variable_vector_t scalars_from, scalars_to;
    for (unsigned k = 0, n = from.size(); k < n; ++k) {
      if (from[k].get_type() != to[k].get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve the same type");
      }
      if (!from[k].get_type().is_array()) {
        scalars_from.push_back(from[k]);
        scalars_to.push_back(to[k]);
      }
    }
    array_map_t arrays;
    for (auto &kv : m_arrays) {
      kv.second.seg.rename(from, to);
      variable_t name = kv.first;
      for (unsigned k = 0, n = from.size(); k < n; ++k) {
        if (from[k] == kv.first) {
          name = to[k];
        }
      }
      arrays.emplace(name, std::move(kv.second));
    }
    m_arrays = std::move(arrays);
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
    for (auto &kv : m_arrays) {
      kv.second.seg.assume(csts);
      kv.second.seg.reduce_bounds_from_scalar(m_scalar);
    }
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
    m_arrays.erase(a);
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
    array_state &st = array_for_access(a, size);
    if (i[lhs] != number_t(0)) {
      // lhs := A[lhs + ..]: the index needs the old lhs, the segment would
      // need the new one, so the load is weak
      *this -= lhs;
      st.seg.weak_load(a, lhs, m_scalar);
    } else {
      *this -= lhs;
      st.seg.load(a, lhs, i, m_scalar);
    }
    normalize_arrays();
    reduce_scalar_from_bounds();
  }

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
    // is_strong_update is ignored, as array_adaptive does for arrays it does
    // not smash: whether the update is strong depends on the index (see
    // segmentation::store)
    array_state &st = array_for_access(a, size);
    st.seg.store(a, i, val, m_scalar);
    normalize_arrays();
    reduce_scalar_from_bounds();
  }

  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {
    m_arrays.erase(a);
  }

  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    auto it = m_arrays.find(rhs);
    if (it == m_arrays.end()) {
      m_arrays.erase(lhs);
      return;
    }
    array_state st = it->second;
    st.seg.rename({rhs}, {lhs});
    m_arrays.erase(lhs);
    m_arrays.emplace(lhs, std::move(st));
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

  // Only the scalar state is exported for now
  linear_constraint_system_t to_linear_constraint_system() const override {
    return m_scalar.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return m_scalar.to_disjunctive_linear_constraint_system();
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

  // The segmentation of array a, for tests and clients: {0} T {+oo} if a
  // was not accessed
  segmentation_t segmentation_of(const variable_t &a) const {
    auto it = m_arrays.find(a);
    return it != m_arrays.end() ? it->second.seg : top_segmentation(a);
  }

  // scalar state, then one line per array: A: {0} 0 {i}? T {+oo}
  void write(crab_os &o) const override {
    o << m_scalar;
    for (auto const &kv : m_arrays) {
      o << "\n" << kv.first << ": " << kv.second.seg;
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
