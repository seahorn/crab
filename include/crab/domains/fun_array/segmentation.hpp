#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <cassert>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * An array segmentation (POPL'11, §11):
 *
 *   {B0} P0 {B1}[?] P1 ... {Bn}[?] Pn {+oo}
 *
 * The bounds B0 < B1 < ... < Bn are sets of expressions that are equal in
 * the concrete (BoundDom). The value Pk (ValueDom) describes every element
 * with index in [Bk, Bk+1), where Bn+1 = +oo. A ? after Bk means that the
 * segment ending at Bk may be empty; without it the segment is not empty.
 *
 * Indices are assumed non-negative: B0 always contains 0, and elements at
 * negative indices are not described. There is no array length, so the last
 * segment runs to +oo, which is implicit. It is never empty.
 *
 * The lattice operations first unify the two segmentations so that both have
 * the same bounds, and then combine their values segment by segment.
 **/
template <class BoundDom, class ValueDom> class segmentation {
public:
  using segmentation_t = segmentation<BoundDom, ValueDom>;
  using bound_t = BoundDom;
  using value_t = ValueDom;
  using expression_t = typename BoundDom::expression_t;
  using variable_t = typename BoundDom::variable_t;
  using linear_expression_t = typename BoundDom::linear_expression_t;
  using number_t = typename BoundDom::number_t;

private:
  bool m_is_bottom;
  // B0 ... Bn
  std::vector<BoundDom> m_bounds;
  // m_maybe_empty[k]: the segment ending at m_bounds[k] may be empty.
  // m_maybe_empty[0] is false: no segment ends at B0.
  std::vector<bool> m_maybe_empty;
  // P0 ... Pn: m_values[k] starts at m_bounds[k]
  std::vector<ValueDom> m_values;

  // Two segmentations with the same bounds
  struct unified {
    std::vector<BoundDom> bounds;
    std::vector<bool> left_maybe_empty;
    std::vector<bool> right_maybe_empty;
    std::vector<ValueDom> left_values;
    std::vector<ValueDom> right_values;
    // Unification split or dropped a bound of the right segmentation, or
    // merged two of its segments: the right result may be weaker than the
    // right input.
    bool right_weakened = false;
  };

  enum class unify_op { join, widening, meet, leq };

  static BoundDom zero_bound() {
    return BoundDom(*expression_t::make(linear_expression_t(number_t(0))));
  }

  // Some index lies in the segment m_values[k]
  bool is_non_empty(std::size_t k) const {
    return k + 1 == m_bounds.size() || !m_maybe_empty[k + 1];
  }

  // Remove m_bounds[k], k > 0: its two segments are merged, joining their
  // values. The merged segment is not empty if the second one was not, so
  // the flag of the second one is kept.
  //
  // E.g., dropping {i} (k = 1) from
  //   {0} A==0 {i}? A==1 {n} A==5 {+oo}
  // gives
  //   {0} 0<=A<=1 {n} A==5 {+oo}
  // where [0, n) is not empty because [i, n) was not.
  void drop_bound(std::size_t k) {
    assert(k > 0);
    m_values[k - 1] = m_values[k - 1] | m_values[k];
    m_values.erase(m_values.begin() + k);
    m_maybe_empty.erase(m_maybe_empty.begin() + k);
    m_bounds.erase(m_bounds.begin() + k);
  }

  // Bottom if some segment that is not empty has a bottom value.
  //
  // E.g., {0} _|_ {i} T {+oo} is bottom: [0, i) has at least one element
  // and no value satisfies _|_. {0} _|_ {i}? T {+oo} is not: [0, i) may be
  // empty (i == 0). The last segment [Bn, +oo) is never empty.
  void normalize() {
    if (m_is_bottom) {
      return;
    }
    for (std::size_t k = 0, n = m_values.size(); k < n; ++k) {
      if (is_non_empty(k) && m_values[k].is_bottom()) {
        set_to_bottom();
        return;
      }
    }
  }

  // Purification (Clousot's ArraySegmentation.Purify): drop the expressions
  // whose variable is not in vars; constants are kept, so B0 keeps 0. A
  // bound left empty is removed and its two segments are merged.
  //
  // E.g., with vars = {b} (the variables of the other segmentation),
  //   {0 t} P0 {a} P1 {b, a+1} P2 {+oo}
  // becomes
  //   {0} P0|P1 {b} P2 {+oo}
  // t, a and a+1 are dropped; {a} is left empty, so it is removed and its
  // segments are joined. This is POPL'11's Example 10: without it, the
  // one-bound lookahead of unify could lose {b}.
  segmentation_t purify(const std::set<variable_t> &vars) const {
    segmentation_t res(*this);
    for (std::size_t k = 0; k < res.m_bounds.size();) {
      BoundDom b;
      for (auto const &e : res.m_bounds[k].expressions()) {
        auto x = e.variable();
        if (!x || vars.count(*x) > 0) {
          b.add(e);
        }
      }
      if (b.is_empty()) {
        res.drop_bound(k);
      } else {
        res.m_bounds[k] = b;
        ++k;
      }
    }
    return res;
  }

  // Both bounds have a constant and a's is smaller (SegmentLimit.IsAtLeftOf).
  //
  // E.g., {1, i} is at the left of {2}; {i} is not at the left of {2}, nor
  // is {2} at the left of {i}: without a constant on both sides there is no
  // answer.
  static bool is_at_left_of(const BoundDom &a, const BoundDom &b) {
    auto ka = a.find_constant();
    auto kb = b.find_constant();
    return ka && kb && *ka < *kb;
  }

  /**
   * Segmentation unification, following Clousot's
   * ArraySegmentation.TryUnifySegments. The inputs are scanned left to
   * right, a being the current bound of l and b the current bound of r; the
   * cases are described in the loop. N is the neutral element of the side
   * being split: the new segments it fills are empty in the concrete, so any
   * value is sound for them. Dropping a bound merges its two segments.
   *
   * Both B0 contain 0, so the first bounds always share 0 and are handled by
   * cases 3 to 6; the result's B0 contains 0 as well.
   *
   * Unlike Clousot, when one input runs out of bounds before the other, the
   * remaining bounds of the other are dropped, except for meet, which adds
   * them to the exhausted input as ? bounds (cases 1 and 2).
   *
   * E.g. (POPL'11, Example 8, the first iteration of an initialization
   * loop), for join, with N = _|_:
   *
   *   l: {0 i} T {+oo}                r: {0 i-1} 0 {1 i} T {+oo}
   *
   *   case 6, {0 i} and {0 i-1} share {0}:
   *   l: {0} N {i}? T                 r: {0} N {i-1}? 0 {1 i} T
   *   case 7b, {i} is in r's next bound {1 i}: drop {i-1} from r
   *   l: {0} N {i}? T                 r: {0} N|0 {1 i} T
   *   case 4, {i} is a strict subset of {1 i}: split r
   *   l: {0} N {i}? T                 r: {0} 0 {i} N {1}? T
   *   case 1, l has no bounds left: drop {1} from r
   *   l: {0} N {i}? T {+oo}           r: {0} 0 {i} N|T {+oo}
   *
   * The join then gives {0} 0 {i}? T {+oo}.
   **/
  static unified unify(unify_op op, segmentation_t l, segmentation_t r,
                       const ValueDom &neutral_l, const ValueDom &neutral_r) {
    unified u;

    std::vector<BoundDom> &lb = l.m_bounds;
    std::vector<BoundDom> &rb = r.m_bounds;
    std::vector<bool> &lq = l.m_maybe_empty;
    std::vector<bool> &rq = r.m_maybe_empty;
    const std::vector<ValueDom> &lv = l.m_values;
    const std::vector<ValueDom> &rv = r.m_values;
    const std::size_t nl = lb.size();
    const std::size_t nr = rb.size();
    std::size_t i = 0, j = 0;
    assert(nl > 0 && nr > 0 && !(lb[0] & rb[0]).is_empty());

    // Add a bound to both results. ql and qr are the flags of the segments
    // ending at it; vl and vr the values of the segments starting at it.
    auto keep = [&u](const BoundDom &bound, bool ql, bool qr,
                     const ValueDom &vl, const ValueDom &vr) {
      u.bounds.push_back(bound);
      u.left_maybe_empty.push_back(ql);
      u.right_maybe_empty.push_back(qr);
      u.left_values.push_back(vl);
      u.right_values.push_back(vr);
    };
    // Drop the current bound of one side: the segment starting at it is
    // merged into the last segment of that side's result.
    auto drop_left = [&]() {
      u.left_values.back() = u.left_values.back() | lv[i];
      ++i;
    };
    auto drop_right = [&]() {
      u.right_values.back() = u.right_values.back() | rv[j];
      ++j;
      u.right_weakened = true;
    };

    while (i < nl || j < nr) {
      // Case 1: l has no bounds left, r has.
      //   meet: add b to l as a ? bound; l's last segment is split, both
      //         parts keep its value.
      //   join, widening, leq: drop b, merging r's segments around it.
      if (i == nl) {
        if (op == unify_op::meet) {
          const ValueDom last = u.left_values.back();
          keep(rb[j], true, rq[j], last, rv[j]);
          ++j;
        } else {
          drop_right();
        }
        continue;
      }
      // Case 2: r has no bounds left, l has. Symmetric to case 1.
      if (j == nr) {
        if (op == unify_op::meet) {
          const ValueDom last = u.right_values.back();
          keep(lb[i], lq[i], true, lv[i], last);
          ++i;
        } else {
          drop_left();
        }
        continue;
      }

      const BoundDom a = lb[i];
      const BoundDom b = rb[j];
      const bool a_in_b = a <= b;
      const bool b_in_a = b <= a;

      // Case 3: same bound.
      //   l: {a} P     r: {a} Q     =>  keep {a}
      if (a_in_b && b_in_a) {
        keep(a, lq[i], rq[j], lv[i], rv[j]);
        ++i;
        ++j;
        continue;
      }
      // Case 4: a is a strict subset of b. Split r.
      //   l: {a} P     r: {a b'} Q  =>  r: {a} N {b'}? Q
      // The segment [a, b') is empty in r, since a == b' there. r's bound
      // b' is compared next against l's next bound.
      if (a_in_b) {
        keep(a, lq[i], rq[j], lv[i], neutral_r);
        rb[j] = b - a;
        rq[j] = true;
        ++i;
        u.right_weakened = true;
        continue;
      }
      // Case 5: b is a strict subset of a. Split l, symmetric to case 4.
      //   l: {b a'} P  r: {b} Q     =>  l: {b} N {a'}? P
      if (b_in_a) {
        keep(b, lq[i], rq[j], neutral_l, rv[j]);
        lb[i] = a - b;
        lq[i] = true;
        ++j;
        continue;
      }

      // Case 6: a and b intersect in c, and neither contains the other.
      // Keep c and split both.
      //   l: {c a'} P  r: {c b'} Q  =>  l: {c} N {a'}? P   r: {c} N {b'}? Q
      // a' and b' are compared next.
      const BoundDom c = a & b;
      if (!c.is_empty()) {
        keep(c, lq[i], rq[j], neutral_l, neutral_r);
        lb[i] = a - c;
        lq[i] = true;
        rb[j] = b - c;
        rq[j] = true;
        u.right_weakened = true;
        continue;
      }

      // Case 7: a and b are disjoint. Drop one or both, with a lookahead of
      // one bound.
      if (j + 1 < nr && !(a & rb[j + 1]).is_empty()) {
        // 7a: a meets r's next bound, so b is extra in r.
        //   l: {a} P     r: {b} Q {a ..} R  =>  drop b
        drop_right();
      } else if (i + 1 < nl && !(lb[i + 1] & b).is_empty()) {
        // 7b: b meets l's next bound, so a is extra in l.
        //   l: {a} P {b ..} R  r: {b} Q     =>  drop a
        drop_left();
      } else if (is_at_left_of(a, b)) {
        // 7c: by their constants, a < b: drop a, the one r went past.
        drop_left();
      } else if (is_at_left_of(b, a)) {
        // 7d: by their constants, b < a: drop b.
        drop_right();
      } else if (a.find_constant()) {
        // 7e: only a has a constant: keep it, drop b.
        drop_right();
      } else if (b.find_constant()) {
        // 7f: only b has a constant: keep it, drop a.
        drop_left();
      } else {
        // 7g: nothing to go by: drop both.
        drop_left();
        drop_right();
      }
    }
    return u;
  }

  // The variables mentioned by the bounds
  std::set<variable_t> variables() const {
    std::set<variable_t> res;
    for (auto const &b : m_bounds) {
      for (auto const &x : b.variables()) {
        res.insert(x);
      }
    }
    return res;
  }

  template <class Combine>
  static segmentation_t combine(const unified &u, bool maybe_empty_or,
                                Combine combine_values) {
    std::vector<bool> maybe_empty;
    for (std::size_t k = 0, n = u.bounds.size(); k < n; ++k) {
      maybe_empty.push_back(maybe_empty_or ? (u.left_maybe_empty[k] ||
                                              u.right_maybe_empty[k])
                                           : (u.left_maybe_empty[k] &&
                                              u.right_maybe_empty[k]));
    }
    std::vector<ValueDom> values;
    for (std::size_t k = 0, n = u.left_values.size(); k < n; ++k) {
      values.push_back(combine_values(u.left_values[k], u.right_values[k]));
    }
    segmentation_t res(u.bounds, std::move(maybe_empty), std::move(values));
    return res;
  }

public:
  // {0} v {+oo}
  explicit segmentation(const ValueDom &v) : m_is_bottom(false) {
    m_bounds.push_back(zero_bound());
    m_maybe_empty.push_back(false);
    m_values.push_back(v);
    normalize();
  }

  // {bounds[0]} values[0] {bounds[1]}[?] values[1] ... {bounds[n]}[?]
  // values[n] {+oo}, where maybe_empty[k] is the ? after bounds[k]. The
  // bounds must be non-empty and in increasing order, bounds[0] must contain
  // 0 and maybe_empty[0] must be false.
  segmentation(std::vector<BoundDom> bounds, std::vector<bool> maybe_empty,
               std::vector<ValueDom> values)
      : m_is_bottom(false), m_bounds(std::move(bounds)),
        m_maybe_empty(std::move(maybe_empty)), m_values(std::move(values)) {
    if (m_bounds.empty() || m_maybe_empty.size() != m_bounds.size() ||
        m_values.size() != m_bounds.size()) {
      CRAB_ERROR("segmentation: expected n > 0 bounds, n flags and n values");
    }
    if (!(zero_bound() <= m_bounds[0])) {
      CRAB_ERROR("segmentation: the first bound must contain 0");
    }
    if (m_maybe_empty[0]) {
      CRAB_ERROR("segmentation: no segment ends at the first bound");
    }
    for (auto const &b : m_bounds) {
      if (b.is_empty()) {
        CRAB_ERROR("segmentation: empty bound");
      }
    }
    normalize();
  }

  segmentation_t make_bottom() const {
    segmentation_t res(*this);
    res.set_to_bottom();
    return res;
  }

  segmentation_t make_top() const {
    return segmentation_t(m_values[0].make_top());
  }

  // m_values keeps one value, to build new values from
  void set_to_bottom() {
    m_is_bottom = true;
    m_bounds.assign(1, zero_bound());
    m_maybe_empty.assign(1, false);
    m_values.resize(1, m_values[0]);
    m_values[0].set_to_bottom();
  }

  bool is_bottom() const { return m_is_bottom; }

  // {0} T {+oo}
  bool is_top() const {
    return !m_is_bottom && m_bounds.size() == 1 &&
           m_bounds[0].expressions().size() == 1 && m_values[0].is_top();
  }

  // Number of bounds, B0 included
  std::size_t size() const { return m_bounds.size(); }
  const BoundDom &bound(std::size_t k) const { return m_bounds[k]; }
  // The segment ending at bound(k) may be empty
  bool maybe_empty(std::size_t k) const { return m_maybe_empty[k]; }
  // value(k) starts at bound(k)
  const ValueDom &value(std::size_t k) const { return m_values[k]; }

  // Clousot's ArraySegmentation.LessEqual. False if unification weakens o,
  // since then the comparison would be against something weaker than o.
  //
  // E.g.,
  //   {0} 0 {i n} T {+oo}  <=  {0} 0 {i}? T {+oo}
  // holds: the left also knows i == n and that [0, i) is not empty. The
  // converse does not. And
  //   {0} 0 {a} 1 {b} 2 {+oo}  <=  {0} 0 {a b} 2 {+oo}
  // is false: the right requires a == b. Unification splits the right into
  // {0} 0 {a} N {b}? 2, which would compare as included, so a weakened
  // right answers false.
  bool operator<=(const segmentation_t &o) const {
    if (m_is_bottom) {
      return true;
    }
    if (o.m_is_bottom) {
      return false;
    }
    unified u = unify(unify_op::leq, *this, o, m_values[0].make_bottom(),
                      m_values[0].make_top());
    if (u.right_weakened) {
      return false;
    }
    for (std::size_t k = 0, n = u.bounds.size(); k < n; ++k) {
      if (u.left_maybe_empty[k] && !u.right_maybe_empty[k]) {
        return false;
      }
    }
    for (std::size_t k = 0, n = u.left_values.size(); k < n; ++k) {
      if (!(u.left_values[k] <= u.right_values[k])) {
        return false;
      }
    }
    return true;
  }

  // Clousot's ArraySegmentation.Join: purify, unify with bottom as neutral
  // element, then join values; a segment may be empty if it may be empty
  // on either side.
  //
  // E.g.,
  //   {0 i} T {+oo}  |  {0 i-1} 0 {1 i} T {+oo}  =  {0} 0 {i}? T {+oo}
  // (see unify for the steps), and
  //   {0} A==1 {i n} T {+oo}  |  {0} A==2 {i m} T {+oo}
  //     =  {0} 1<=A<=2 {i} T {+oo}
  // where n and m, each on one side only, are dropped by purification.
  segmentation_t operator|(const segmentation_t &o) const {
    if (m_is_bottom) {
      return o;
    }
    if (o.m_is_bottom) {
      return *this;
    }
    ValueDom bot = m_values[0].make_bottom();
    unified u = unify(unify_op::join, purify(o.variables()),
                      o.purify(variables()), bot, bot);
    return combine(u, true, [](const ValueDom &x, const ValueDom &y) {
      return x | y;
    });
  }

  // Clousot's ArraySegmentation.Widening: as join, then widen values
  // segment by segment.
  //
  // E.g.,
  //   {0} A==0 {i}? T {+oo}  ||  {0} 0<=A<=1 {i}? T {+oo}
  //     =  {0} A>=0 {i}? T {+oo}
  // The bounds of the result are subsets of the left's bounds, so a
  // sequence of widenings cannot keep adding bounds.
  segmentation_t operator||(const segmentation_t &o) const {
    if (m_is_bottom) {
      return o;
    }
    if (o.m_is_bottom) {
      return *this;
    }
    ValueDom bot = m_values[0].make_bottom();
    unified u = unify(unify_op::widening, purify(o.variables()),
                      o.purify(variables()), bot, bot);
    return combine(u, true, [](const ValueDom &x, const ValueDom &y) {
      return x || y;
    });
  }

  // Clousot's ArraySegmentation.Meet: unify with bottom as neutral element,
  // then meet values; a segment may be empty only if it may be empty on
  // both sides.
  //
  // E.g.,
  //   {0} T {+oo}  &  {0} 0 {i}? T {+oo}  =  {0} 0 {i}? T {+oo}
  // where the left, having no bounds left, gets {i} (case 1 of unify), and
  //   {0} A==1 {1} T {+oo}  &  {0} A==1 {2} T {+oo}  =  {0} A==1 {2} T {+oo}
  // where {1} is dropped from the left (by their constants, 1 < 2) and {2}
  // is added to it.
  segmentation_t operator&(const segmentation_t &o) const {
    if (m_is_bottom || o.m_is_bottom) {
      return make_bottom();
    }
    ValueDom bot = m_values[0].make_bottom();
    unified u = unify(unify_op::meet, *this, o, bot, bot);
    return combine(u, false, [](const ValueDom &x, const ValueDom &y) {
      return x & y;
    });
  }

  // {B0} P0 {B1} P1 {B2}? P2 ... {+oo}
  void write(crab_os &o) const {
    if (m_is_bottom) {
      o << "_|_";
      return;
    }
    for (std::size_t k = 0, n = m_bounds.size(); k < n; ++k) {
      if (k > 0) {
        o << " ";
      }
      o << m_bounds[k];
      if (m_maybe_empty[k]) {
        o << "?";
      }
      o << " " << m_values[k];
    }
    o << " {+oo}";
  }
};

template <class BoundDom, class ValueDom>
inline crab_os &operator<<(crab_os &o,
                           const segmentation<BoundDom, ValueDom> &s) {
  s.write(o);
  return o;
}

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
