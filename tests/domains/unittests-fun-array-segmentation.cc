// Tests for fun_array's segmentation: unification and the lattice
// operations (join, widening, meet, <=), with offset_expression bounds and
// zone_segment_values.
//
// The paper examples are from P. Cousot, R. Cousot, F. Logozzo, "A
// Parametric Segmentation Functor for Fully Automatic and Scalable Array
// Content Analysis", POPL 2011, Sect. 11.4 and 12.2. Segmentations here have
// no array length: the paper's last bound ({n}, {A.Length}) becomes the
// implicit +oo end, and the first bound always contains 0.
//
// Nothing is printed to standard output unless --print-invariants is given:
// the Boost.Test log goes to stderr (tests/run_tests.sh diffs stdout against
// a golden file) and ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_fun_array_segmentation
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../crab_lang.hpp"
#include <crab/domains/fun_array/expression_set_bound.hpp>
#include <crab/domains/fun_array/offset_expression.hpp>
#include <crab/domains/fun_array/segmentation.hpp>
#include <crab/domains/fun_array/zone_segment_values.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/support/os.hpp>

#include <initializer_list>
#include <string>
#include <vector>

using namespace crab::domains;
using namespace crab::domains::fun_array_impl;
using namespace ikos;

namespace {
using varname_t = crab::cfg_impl::varname_t;
using scalar_t = split_dbm_domain<z_number, varname_t>;
using expr_t = offset_expression<scalar_t>;
using bound_t = expression_set_bound<expr_t, scalar_t>;
using value_t = zone_segment_values<scalar_t>;
using seg_t = segmentation<bound_t, value_t>;
using var_t = scalar_t::variable_t;
using lin_t = scalar_t::linear_expression_t;

// Set by --print-invariants
bool print_invariants = false;

void print(const std::string &name, const seg_t &s) {
  if (print_invariants) {
    crab::outs() << name << " = " << s << "\n";
  }
}

bool same_bound(const bound_t &x, const bound_t &y) {
  return x <= y && y <= x;
}

bool same_value(const value_t &x, const value_t &y) {
  return x <= y && y <= x;
}

// s is {bounds[0]} values[0] {bounds[1]}[?] values[1] ... {+oo}, where
// maybe_empty[k] is the ? after bounds[k]
void check(const seg_t &s, const std::vector<bound_t> &bounds,
           const std::vector<bool> &maybe_empty,
           const std::vector<value_t> &values) {
  BOOST_TEST_REQUIRE(!s.is_bottom());
  BOOST_TEST_REQUIRE(s.size() == bounds.size());
  for (std::size_t k = 0; k < bounds.size(); ++k) {
    BOOST_TEST(same_bound(s.bound(k), bounds[k]), "bound " << k);
    BOOST_TEST(s.maybe_empty(k) == maybe_empty[k], "? after bound " << k);
    BOOST_TEST(same_value(s.value(k), values[k]), "value " << k);
  }
}

struct fixture {
  crab::var_factory_impl::str_variable_factory vfac;
  var_t A, i, n, a, b, c;
  value_t top;

  fixture()
      : A(vfac["A"], crab::ARR_INT_TYPE, 32), i(vfac["i"], crab::INT_TYPE, 32),
        n(vfac["n"], crab::INT_TYPE, 32), a(vfac["a"], crab::INT_TYPE, 32),
        b(vfac["b"], crab::INT_TYPE, 32), c(vfac["c"], crab::INT_TYPE, 32),
        top(value_t::make_index_placeholder(A)) {}

  static lin_t k(int x) { return lin_t(z_number(x)); }

  static bound_t B(std::initializer_list<lin_t> es) {
    bound_t res;
    for (auto const &e : es) {
      res.add(*expr_t::make(e));
    }
    return res;
  }

  // Every element of the segment is x
  value_t eq(int x) const {
    value_t v(top);
    v.store(A, k(x), scalar_t());
    return v;
  }
  value_t T() const { return top; }
  value_t bot() const { return top.make_bottom(); }
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(paper_examples, fixture)

// Example 8 (Fig. 2, first iteration of the initialization loop): the loop
// entry {0 i} T {n} and the back edge {0 i-1} 0 {1 i} T {n}? join into
// {0} 0 {i}? T {n}?.
BOOST_AUTO_TEST_CASE(example_8_initialization_loop) {
  seg_t entry({B({k(0), i})}, {false}, {T()});
  seg_t back({B({k(0), lin_t(i) - 1}), B({k(1), i})}, {false, false},
             {eq(0), T()});
  seg_t j1 = entry | back;
  seg_t j2 = back | entry;
  print("example 8: entry", entry);
  print("example 8: back edge", back);
  print("example 8: entry | back edge", j1);

  check(j1, {B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  check(j2, {B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  BOOST_TEST((entry <= j1));
  BOOST_TEST((back <= j1));
  BOOST_TEST(!(j1 <= entry));
  BOOST_TEST(!(j1 <= back));

  // The invariant is stable under widening with the next iterate
  seg_t w = j1 || (entry | back);
  check(w, {B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
}

// Example 5: when unifying {0}...{a b} with {0}...{a b c}, the paper notes
// that {0}...{a} and {0}...{b} are both minimal solutions. The common
// expressions a and b are kept together; c, which only the right side has,
// is dropped.
BOOST_AUTO_TEST_CASE(example_5_common_expressions_kept) {
  seg_t l({B({k(0)}), B({a, b})}, {false, false}, {eq(1), eq(2)});
  seg_t r({B({k(0)}), B({a, b, c})}, {false, false}, {eq(1), eq(3)});
  seg_t j = l | r;
  print("example 5: join", j);

  check(j, {B({k(0)}), B({a, b})}, {false, false},
        {eq(1), eq(2) | eq(3)});
}

// Example 6: when unifying {0}...{a}...{b}...{c} with {0}...{b}...{a}...{c},
// {0}...{a}...{c} and {0}...{b}...{c} are both maximal solutions. The paper's
// algorithm returns {0}...{c}; the one-bound lookahead of Clousot's
// TryUnifySegments finds a in the right side's next bound, drops b there,
// and returns the maximal solution {0}...{a}...{c}.
BOOST_AUTO_TEST_CASE(example_6_lookahead) {
  value_t p0 = eq(0), p1 = eq(1), p2 = eq(2), p3 = eq(3);
  value_t q0 = eq(10), q1 = eq(11), q2 = eq(12), q3 = eq(13);
  seg_t l({B({k(0)}), B({a}), B({b}), B({c})}, {false, false, false, false},
          {p0, p1, p2, p3});
  seg_t r({B({k(0)}), B({b}), B({a}), B({c})}, {false, false, false, false},
          {q0, q1, q2, q3});
  seg_t j = l | r;
  print("example 6: join", j);

  check(j, {B({k(0)}), B({a}), B({c})}, {false, false, false},
        {p0 | q0 | q1, p1 | p2 | q2, p3 | q3});
}

// Example 7: the coarsest common refinement of {0 i-1} P1 {i} and
// {0 i-2} P2 {i} needs the bound {i-2} below {0}, which is undefined when
// i < 2. Without scalar information, unification keeps {0} and {i} and
// drops i-1 and i-2, merging the segments between them.
BOOST_AUTO_TEST_CASE(example_7_no_bounds_below_zero) {
  seg_t l({B({k(0), lin_t(i) - 1}), B({i})}, {false, false}, {eq(1), T()});
  seg_t r({B({k(0), lin_t(i) - 2}), B({i})}, {false, false}, {eq(2), T()});
  seg_t j = l | r;
  print("example 7: join", j);

  check(j, {B({k(0)}), B({i})}, {false, false}, {eq(1) | eq(2), T()});
}

// Example 10 (purification): unifying {0}...{a}...{b}...{A.Length} with
// {0}...{b}...{A.Length} keeps {b}: a appears only on the left, so the left
// bound {a} is removed first.
BOOST_AUTO_TEST_CASE(example_10_purification) {
  value_t p0 = eq(0), p1 = eq(1), p2 = eq(2);
  value_t q0 = eq(10), q1 = eq(11);
  seg_t l({B({k(0)}), B({a}), B({b})}, {false, false, false}, {p0, p1, p2});
  seg_t r({B({k(0)}), B({b})}, {false, false}, {q0, q1});
  seg_t j = l | r;
  print("example 10: join", j);

  check(j, {B({k(0)}), B({b})}, {false, false}, {p0 | p1 | q0, p2 | q1});
}

// Example 11: Clousot materializes constant bounds, unifying
// {0}...{25}...{26} with {0}...{20}...{30} into
// {0}...{20}...{25}...{26}...{30}. That refinement is not implemented: the
// disjoint constant bounds are dropped, the smaller one first.
BOOST_AUTO_TEST_CASE(example_11_constant_refinement_not_implemented) {
  value_t p0 = eq(0), p1 = eq(1), p2 = eq(2);
  value_t q0 = eq(10), q1 = eq(11), q2 = eq(12);
  seg_t l({B({k(0)}), B({k(25)}), B({k(26)})}, {false, false, false},
          {p0, p1, p2});
  seg_t r({B({k(0)}), B({k(20)}), B({k(30)})}, {false, false, false},
          {q0, q1, q2});
  seg_t j = l | r;
  print("example 11: join", j);

  check(j, {B({k(0)})}, {false}, {p0 | p1 | p2 | q0 | q1 | q2});
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(lattice, fixture)

// {0} T {+oo} is top, and the identity of meet
BOOST_AUTO_TEST_CASE(untouched_array) {
  seg_t untouched(T());
  seg_t s({B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  print("untouched", untouched);
  print("untouched & s", untouched & s);

  BOOST_TEST(untouched.is_top());
  BOOST_TEST(!s.is_top());
  BOOST_TEST((s <= untouched));
  BOOST_TEST(!(untouched <= s));
  BOOST_TEST((untouched | s).is_top());
  check(untouched & s, {B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
}

// A bound in the first bound is information: {0 i} says i == 0
BOOST_AUTO_TEST_CASE(first_bound_with_variable_is_not_top) {
  seg_t s({B({k(0), i})}, {false}, {T()});
  BOOST_TEST(!s.is_top());
  BOOST_TEST(!(seg_t(T()) <= s));
}

BOOST_AUTO_TEST_CASE(leq_reflexive_and_join_upper_bound) {
  seg_t s1({B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  seg_t s2({B({k(0)}), B({i, n})}, {false, false}, {eq(0), T()});
  seg_t j = s1 | s2;
  print("s1 | s2", j);

  BOOST_TEST((s1 <= s1));
  BOOST_TEST((s2 <= s2));
  BOOST_TEST((s1 <= j));
  BOOST_TEST((s2 <= j));
  // s2 is more precise: i == n and the segment [0, i) is not empty
  BOOST_TEST((s2 <= s1));
  BOOST_TEST(!(s1 <= s2));
}

// The right side requires a == b; the left side does not, so it is not
// included, although unifying them splits the right side into
// {a} N {b}? Q, which would compare as included.
BOOST_AUTO_TEST_CASE(leq_false_when_right_is_split) {
  seg_t l({B({k(0)}), B({a}), B({b})}, {false, false, false},
          {eq(0), eq(1), eq(2)});
  seg_t r({B({k(0)}), B({a, b})}, {false, false}, {eq(0), eq(2)});
  BOOST_TEST(!(l <= r));
  // r has the segment [a, b) empty while l requires it non-empty
  BOOST_TEST(!(r <= l));
  // with [a, b) possibly empty in l, r is included in l
  seg_t lq({B({k(0)}), B({a}), B({b})}, {false, false, true},
           {eq(0), eq(1), eq(2)});
  BOOST_TEST((r <= lq));
}

// Meet keeps the bounds of the side that has more of them
BOOST_AUTO_TEST_CASE(meet_adds_bounds) {
  seg_t s5({B({k(0)}), B({k(1)})}, {false, false}, {eq(1), T()});
  seg_t s6({B({k(0)}), B({k(2)})}, {false, false}, {eq(1), T()});
  seg_t m = s5 & s6;
  print("s5 & s6", m);

  check(m, {B({k(0)}), B({k(2)})}, {false, false}, {eq(1), T()});
}

BOOST_AUTO_TEST_CASE(widening_widens_values) {
  seg_t s1({B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  seg_t s2({B({k(0)}), B({i})}, {false, true}, {eq(0) | eq(1), T()});
  seg_t w = s1 || s2;
  print("s1 || s2", w);

  BOOST_TEST_REQUIRE(w.size() == 2u);
  BOOST_TEST((s1 <= w));
  BOOST_TEST((s2 <= w));
  BOOST_TEST(!(w <= s2));
}

// A segment that is not empty and has a bottom value makes the whole
// segmentation bottom; a segment that may be empty does not.
BOOST_AUTO_TEST_CASE(bottom_segments) {
  seg_t non_empty({B({k(0)}), B({i})}, {false, false}, {bot(), T()});
  seg_t maybe_empty({B({k(0)}), B({i})}, {false, true}, {bot(), T()});
  seg_t last({B({k(0)}), B({i})}, {false, true}, {T(), bot()});
  BOOST_TEST(non_empty.is_bottom());
  BOOST_TEST(!maybe_empty.is_bottom());
  BOOST_TEST(last.is_bottom());

  seg_t s({B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  check(non_empty | s, {B({k(0)}), B({i})}, {false, true}, {eq(0), T()});
  BOOST_TEST((non_empty <= s));
  BOOST_TEST((non_empty & s).is_bottom());
}

BOOST_AUTO_TEST_SUITE_END()

// Keep stdout empty unless --print-invariants is given: the golden-output
// harness (tests/run_tests.sh) diffs the standard output of every test
// binary, so the whole Boost.Test log is routed to stderr. ctest passes
// --disable-warnings, which Boost.Test would reject, so it is translated
// into the crab flag it stands for.
int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int idx = 0; idx < argc; ++idx) {
    std::string arg(argv[idx]);
    if (arg == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    if (arg == "--print-invariants") {
      print_invariants = true;
      continue;
    }
    args.push_back(argv[idx]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
