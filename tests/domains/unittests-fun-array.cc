// Tests for single operations of the fun_array domain with one array per
// group: loads, stores, scalar assignments and joins. Whole loops run
// through the fixpoint are in fun-array.cc.
//
// Nothing is printed to standard output unless --print-invariants is given:
// the Boost.Test log goes to stderr (tests/run_tests.sh diffs stdout against
// a golden file) and ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_fun_array
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../crab_lang.hpp"
#include <crab/domains/fun_array.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/support/os.hpp>

#include <initializer_list>
#include <string>
#include <vector>

using namespace crab::cfg_impl;
using namespace crab::domains;
using namespace ikos;

namespace {
using scalar_t = split_dbm_domain<z_number, varname_t>;
using dom_t = fun_array<scalar_t>;
using seg_t = dom_t::segmentation_t;
using bound_t = seg_t::bound_t;
using expr_t = dom_t::expression_t;
using var_t = dom_t::variable_t;
using lin_t = dom_t::linear_expression_t;
using cst_t = dom_t::linear_constraint_t;

// Set by --print-invariants
bool print_invariants = false;

void print(const std::string &name, const dom_t &d) {
  if (print_invariants) {
    crab::outs() << name << " = " << d << "\n";
  }
}

lin_t k(int x) { return lin_t(z_number(x)); }

bound_t B(std::initializer_list<lin_t> es) {
  bound_t res;
  for (auto const &e : es) {
    res.add(*expr_t::make(e));
  }
  return res;
}

bool same_bound(const bound_t &x, const bound_t &y) {
  return x <= y && y <= x;
}

// The bounds and ? marks of s are bounds and maybe_empty
void check_bounds(const seg_t &s, const std::vector<bound_t> &bounds,
                  const std::vector<bool> &maybe_empty) {
  BOOST_TEST_REQUIRE(!s.is_bottom());
  BOOST_TEST_REQUIRE(s.size() == bounds.size());
  for (std::size_t n = 0; n < bounds.size(); ++n) {
    BOOST_TEST(same_bound(s.bound(n), bounds[n]), "bound " << n);
    BOOST_TEST(s.maybe_empty(n) == maybe_empty[n], "? after bound " << n);
  }
}

struct fixture {
  variable_factory_t vfac;
  var_t A, i, j, n, x;
  fixture()
      : A(vfac["A"], crab::ARR_INT_TYPE, 32), i(vfac["i"], crab::INT_TYPE, 32),
        j(vfac["j"], crab::INT_TYPE, 32), n(vfac["n"], crab::INT_TYPE, 32),
        x(vfac["x"], crab::INT_TYPE, 32) {}
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(transfer_functions, fixture)

BOOST_AUTO_TEST_CASE(store_then_load) {
  dom_t d;
  d.array_store(A, k(4), i, k(5), true);
  d.array_load(x, A, k(4), i);
  print("A[i] := 5; x := A[i]", d);
  BOOST_TEST(d.entails(cst_t(lin_t(x) - 5, cst_t::EQUALITY)));
  // the access assumed i >= 0
  BOOST_TEST(d.entails(cst_t(-lin_t(i), cst_t::INEQUALITY)));
  check_bounds(d.segmentation_of(A), {B({k(0)}), B({i}), B({lin_t(i) + 1})},
               {false, true, false});
}

// A store at an index equal to an existing bound reuses it
BOOST_AUTO_TEST_CASE(store_at_existing_bound) {
  dom_t d;
  d.assign(i, k(0));
  d.array_store(A, k(4), i, k(0), true);
  print("i := 0; A[i] := 0", d);
  check_bounds(d.segmentation_of(A), {B({k(0), i}), B({k(1), lin_t(i) + 1})},
               {false, false});
}

// i := i + 1 rewrites the bounds; i := n replaces them
BOOST_AUTO_TEST_CASE(scalar_assignments) {
  dom_t d;
  d.array_store(A, k(4), i, k(0), true);
  d.apply(OP_ADDITION, i, i, z_number(1));
  print("A[i] := 0; i := i + 1", d);
  check_bounds(d.segmentation_of(A), {B({k(0)}), B({lin_t(i) - 1}), B({i})},
               {false, true, false});

  d.assign(i, n);
  print("i := n", d);
  check_bounds(d.segmentation_of(A), {B({k(0)})}, {false});

  dom_t e;
  e.array_store(A, k(4), n, k(0), true);
  e.assign(j, n);
  print("A[n] := 0; j := n", e);
  // j joins the bound holding n, not the one holding n+1
  check_bounds(e.segmentation_of(A),
               {B({k(0)}), B({n, j}), B({lin_t(n) + 1})},
               {false, true, false});
}

BOOST_AUTO_TEST_CASE(negative_index_is_unreachable) {
  dom_t d;
  d += cst_t(lin_t(i) + 1, cst_t::INEQUALITY); // i <= -1
  d.array_store(A, k(4), i, k(0), true);
  BOOST_TEST(d.is_bottom());
}

// An access with another element size forgets the array
BOOST_AUTO_TEST_CASE(element_size_change_forgets) {
  dom_t d;
  d.array_store(A, k(4), i, k(5), true);
  d.array_load(x, A, k(8), i);
  BOOST_TEST(!d.entails(cst_t(lin_t(x) - 5, cst_t::EQUALITY)));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(reduction, fixture)

// A guard reaches the segment values that mention its variables
BOOST_AUTO_TEST_CASE(guard_reaches_values) {
  var_t m(vfac["m"], crab::INT_TYPE, 32);
  dom_t d;
  d.array_load(x, A, k(4), i); // [i, i+1) gets A == x
  d += cst_t(lin_t(x) - m, cst_t::INEQUALITY); // x <= m
  d.assign(x, k(7)); // x no longer the element
  print("x := A[i]; x <= m; x := 7", d);
  // A[i] <= m survives x := 7
  dom_t t(d);
  t.array_load(n, A, k(4), i);
  BOOST_TEST(t.entails(cst_t(lin_t(n) - m, cst_t::INEQUALITY)));
}

// An assignment is replayed in the values that mention its variable,
// instead of forgetting it
BOOST_AUTO_TEST_CASE(assignment_reaches_values) {
  var_t m(vfac["m"], crab::INT_TYPE, 32);
  dom_t d;
  d.array_load(m, A, k(4), i);              // [i, i+1): A == m
  d += cst_t(lin_t(m) - x + 1, cst_t::INEQUALITY); // m <= x - 1
  d.assign(m, x);                           // [i, i+1): A < m
  print("m := A[i]; m < x; m := x", d);
  dom_t t(d);
  t.array_load(n, A, k(4), i);
  BOOST_TEST(t.entails(cst_t(lin_t(n) - m + 1, cst_t::INEQUALITY)));
}

// POPL'11, Fig. 5, point 6: with i == 0, {0} T {i}? 0 {n}? becomes
// {0 i} 0 {n}?: the segment [0, i) is empty.
BOOST_AUTO_TEST_CASE(guard_merges_bounds_proven_equal) {
  dom_t d;
  d += cst_t(-lin_t(i), cst_t::INEQUALITY);      // 0 <= i
  d += cst_t(lin_t(i) - n + 1, cst_t::INEQUALITY); // i < n
  d.array_store(A, k(4), i, k(0), true);          // {0} T {i}? 0 {i+1} T
  d += cst_t(lin_t(i), cst_t::EQUALITY);          // i == 0
  print("{0} T {i}? 0 {i+1} T with i == 0", d);
  // only the bounds of the empty segment merge; {i+1} does not gain 1
  check_bounds(d.segmentation_of(A), {B({k(0), i}), B({lin_t(i) + 1})},
               {false, false});
}

BOOST_AUTO_TEST_CASE(guard_removes_maybe_empty) {
  dom_t d;
  d.array_store(A, k(4), i, k(0), true); // {0} T {i}? 0 {i+1} T
  d += cst_t(-lin_t(i) + 1, cst_t::INEQUALITY); // i >= 1
  check_bounds(d.segmentation_of(A), {B({k(0)}), B({i}), B({lin_t(i) + 1})},
               {false, false, false});
}

// A segment that cannot be empty between bounds proven equal is
// unreachable
BOOST_AUTO_TEST_CASE(guard_detects_unreachable) {
  dom_t d;
  d.array_store(A, k(4), i, k(0), true);
  d += cst_t(-lin_t(i) + 1, cst_t::INEQUALITY); // i >= 1: [0, i) not empty
  BOOST_TEST(!d.is_bottom());
  d.assign(j, i);                                 // {0} T {i j} 0 ..
  d += cst_t(lin_t(j), cst_t::INEQUALITY);        // j <= 0, so i == 0
  BOOST_TEST(d.is_bottom());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(array_groups, fixture)

// x := A[idx] and y := Bs[idx] in a copy of d entail x == y
bool equal_elements(const dom_t &d, const var_t &x, const var_t &y,
                    const var_t &A, const var_t &Bs, const var_t &idx) {
  dom_t t(d);
  t.array_load(x, A, k(4), idx);
  t.array_load(y, Bs, k(4), idx);
  return t.entails(cst_t(lin_t(x) - y, cst_t::EQUALITY));
}

// t := B[i]; A[i] := t puts A and B in one group, with A == B on [i, i+1)
BOOST_AUTO_TEST_CASE(store_merges_groups) {
  var_t Bs(vfac["B"], crab::ARR_INT_TYPE, 32);
  var_t t(vfac["t"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32);
  dom_t d;
  d.array_load(t, Bs, k(4), i);
  BOOST_TEST(!d.same_group(A, Bs));
  d.array_store(A, k(4), i, t, true);
  print("t := B[i]; A[i] := t", d);
  BOOST_TEST(d.same_group(A, Bs));
  d.assign(t, k(0)); // t is no longer the element
  BOOST_TEST(equal_elements(d, x, y, A, Bs, i));
}

// A store of a value no other group mentions does not merge groups
BOOST_AUTO_TEST_CASE(unrelated_store_does_not_merge) {
  var_t Bs(vfac["B"], crab::ARR_INT_TYPE, 32);
  dom_t d;
  d.array_store(Bs, k(4), i, k(1), true);
  d.array_store(A, k(4), i, k(1), true);
  BOOST_TEST(!d.same_group(A, Bs));
}

// Forgetting A keeps B's facts in their shared group
BOOST_AUTO_TEST_CASE(forget_keeps_other_arrays) {
  var_t Bs(vfac["B"], crab::ARR_INT_TYPE, 32);
  var_t t(vfac["t"], crab::INT_TYPE, 32);
  dom_t d;
  d.array_store(Bs, k(4), i, k(5), true);
  d.array_load(t, Bs, k(4), i);
  d.array_store(A, k(4), i, t, true);
  BOOST_TEST(d.same_group(A, Bs));
  d -= A;
  print("B[i] := 5; t := B[i]; A[i] := t; forget A", d);
  dom_t e(d);
  e.array_load(x, Bs, k(4), i);
  BOOST_TEST(e.entails(cst_t(lin_t(x) - 5, cst_t::EQUALITY)));
}

// A := B relates A and B everywhere
BOOST_AUTO_TEST_CASE(array_assign_relates) {
  var_t Bs(vfac["B"], crab::ARR_INT_TYPE, 32);
  var_t y(vfac["y"], crab::INT_TYPE, 32);
  dom_t d;
  d.array_store(Bs, k(4), i, k(3), true);
  d.array_assign(A, Bs);
  print("B[i] := 3; A := B", d);
  BOOST_TEST(d.same_group(A, Bs));
  dom_t e(d);
  e.array_load(x, A, k(4), i);
  BOOST_TEST(e.entails(cst_t(lin_t(x) - 3, cst_t::EQUALITY)));
  // and at any index, not only where B is known
  BOOST_TEST(equal_elements(d, x, y, A, Bs, n));
}

// A join first gives both sides the same groups: on the right A and B are
// in separate groups, but both elements at i are 1, so the merged group has
// A == B there, as on the left
BOOST_AUTO_TEST_CASE(join_aligns_groups) {
  var_t Bs(vfac["B"], crab::ARR_INT_TYPE, 32);
  var_t t(vfac["t"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32);
  dom_t l;
  l.array_load(t, Bs, k(4), i);
  l.array_store(A, k(4), i, t, true);
  l -= t;
  dom_t r;
  r.array_store(Bs, k(4), i, k(1), true);
  r.array_store(A, k(4), i, k(1), true);
  dom_t j = l | r;
  print("join: left", l);
  print("join: right", r);
  print("join", j);
  BOOST_TEST(j.same_group(A, Bs));
  BOOST_TEST(equal_elements(j, x, y, A, Bs, i));
  BOOST_TEST((l <= j));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(lattice, fixture)

// The first iteration of `i := 0; while (..) { A[i] := 0; i := i + 1; }`:
// the loop entry has not accessed A yet, so its segmentation gains {i}
// from the back edge (its scalar state proves i == 0) before the join.
BOOST_AUTO_TEST_CASE(join_saturates_untouched_array) {
  dom_t entry;
  entry.assign(i, k(0));
  dom_t back(entry);
  back.array_store(A, k(4), i, k(0), true);
  back.apply(OP_ADDITION, i, i, z_number(1));
  dom_t joined = entry | back;
  print("entry", entry);
  print("back edge", back);
  print("entry | back edge", joined);

  check_bounds(joined.segmentation_of(A), {B({k(0)}), B({i})}, {false, true});
  BOOST_TEST((entry <= joined));
  BOOST_TEST((back <= joined));
}

BOOST_AUTO_TEST_SUITE_END()

// Keep stdout empty unless --print-invariants is given: the golden-output
// harness (tests/run_tests.sh) diffs the standard output of every test
// binary, so the whole Boost.Test log is routed to stderr. ctest passes
// --disable-warnings, which Boost.Test would reject, so it is translated
// into the crab flag it stands for.
int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int a = 0; a < argc; ++a) {
    std::string arg(argv[a]);
    if (arg == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    if (arg == "--print-invariants") {
      print_invariants = true;
      continue;
    }
    args.push_back(argv[a]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
