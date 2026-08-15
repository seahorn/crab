// Unit tests for the tDBM domain, organized as Boost.Test suites:
//   transfer_functions      assign/assume/join/meet/widening/forget/rename/
//                           project and entailment over the {2,3,4} template
//   incremental_saturation  facts the lazy variant derives while constraints
//                           are added (no explicit normalize())
//   wide_template           the TVPI C-string example over {10,255}
//   normalization           lazy vs eager Params, idempotence, bottom
//   reduction               the saturation rules one derivation at a time
//   regressions             one case per soundness bug found in review, plus
//                           coverage for every public API (template {1,2,3,4};
//                           coefficient 1 is deliberate: ghost(v,1) == v)
//   division_template_gaps  division fallbacks the {1,2,3,4} template cannot
//                           reach (template {1,6})
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_tvpi
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../../common.hpp"

#include <crab/support/debug.hpp>

#include <string>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using test_domain_t = z_tvpi_dbm_domain_t;
using eager_dom_t =
    crab::domains::tvpi_dbm_domain<z_sdbm_domain_t,
                                   crab::domains::TVPIDBMNormalizeParams>;

namespace {

// The global coefficient template, scoped to a fixture (values must be given
// sorted; restored on scope exit).
struct scoped_template {
  std::vector<unsigned> m_saved;
  explicit scoped_template(std::initializer_list<unsigned> t)
      : m_saved(crab_domain_params_man::get().coefficients()) {
    crab_domain_params_man::get().coefficients().assign(t);
  }
  ~scoped_template() { crab_domain_params_man::get().coefficients() = m_saved; }
};

// Per-case fixture: a fresh variable factory, the usual variables, and the
// coefficient template of the suite.
struct tvpi_fixture {
  scoped_template tpl;
  variable_factory_t vfac;
  z_var x, y, z, w, i, j, k, n, o, u, v, c, s64;
  explicit tvpi_fixture(std::initializer_list<unsigned> t)
      : tpl(t), x(vfac["x"], crab::INT_TYPE, 32),
        y(vfac["y"], crab::INT_TYPE, 32), z(vfac["z"], crab::INT_TYPE, 32),
        w(vfac["w"], crab::INT_TYPE, 32), i(vfac["i"], crab::INT_TYPE, 32),
        j(vfac["j"], crab::INT_TYPE, 32), k(vfac["k"], crab::INT_TYPE, 32),
        n(vfac["n"], crab::INT_TYPE, 32), o(vfac["o"], crab::INT_TYPE, 32),
        u(vfac["u"], crab::INT_TYPE, 32), v(vfac["v"], crab::INT_TYPE, 32),
        c(vfac["c"], crab::INT_TYPE, 32), s64(vfac["s64"], crab::INT_TYPE, 64) {
  }
};

struct template_234 : tvpi_fixture {
  template_234() : tvpi_fixture({2, 3, 4}) {}
};
struct template_1234 : tvpi_fixture {
  template_1234() : tvpi_fixture({1, 2, 3, 4}) {}
};
struct template_10_255 : tvpi_fixture {
  template_10_255() : tvpi_fixture({10, 255}) {}
};
struct template_16 : tvpi_fixture {
  template_16() : tvpi_fixture({1, 6}) {}
};

// Assume lo <= var <= hi.
template <typename Dom>
void set_range(Dom &d, const z_var &var, int lo, int hi) {
  d += (var >= z_number(lo));
  d += (var <= z_number(hi));
}

// Fundamental lattice laws for one join and one meet.
void check_lattice_laws(const test_domain_t &dom1, const test_domain_t &dom2) {
  test_domain_t dom3 = dom1 | dom2;
  BOOST_TEST((dom1 <= dom3), "join soundness: dom1 <= dom1|dom2");
  BOOST_TEST((dom2 <= dom3), "join soundness: dom2 <= dom1|dom2");
  test_domain_t dom4 = dom1 & dom2;
  BOOST_TEST((dom4 <= dom1), "meet soundness: dom1&dom2 <= dom1");
  BOOST_TEST((dom4 <= dom2), "meet soundness: dom1&dom2 <= dom2");
}

// Is `target` (an inequality e <= 0) present, verbatim after ghost
// translation, in dom's constraint system?
bool materialized(const test_domain_t &dom, const z_lin_cst_t &target) {
  for (auto const &cst : dom.to_linear_constraint_system()) {
    if (!cst.is_inequality()) {
      continue;
    }
    auto diff = cst.expression() - target.expression();
    if (diff.is_constant() && diff.constant() == 0) {
      return true;
    }
  }
  return false;
}

bool materialized_after_normalize(const test_domain_t &dom,
                                  const z_lin_cst_t &target) {
  test_domain_t tmp(dom);
  tmp.normalize();
  return materialized(tmp, target);
}

} // namespace

//===----------------------------------------------------------------------===//
// Transfer functions over the {2,3,4} template.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(transfer_functions, template_234)

BOOST_AUTO_TEST_CASE(assign_exact_values) {
  // x=1, y=2x, z=3x+7, n=2x+2y+5, k=5, o=-2k-5
  test_domain_t dom1;
  dom1.assign(x, z_number(1));
  dom1.assign(y, x * z_number(2));
  dom1.assign(z, x * z_number(3) + z_number(7));
  dom1.assign(n, x * z_number(2) + y * z_number(2) + z_number(5));
  dom1.assign(o, z_number(-2) * k - z_number(5));
  dom1.assign(k, z_number(0) * o + z_number(5));
  BOOST_TEST(!dom1.is_bottom());
  BOOST_TEST(dom1.entails(x == z_number(1)));
  BOOST_TEST(dom1.entails(y == z_number(2)));
  BOOST_TEST(dom1.entails(z == z_number(10)));
  BOOST_TEST(dom1.entails(k == z_number(5)));
}

BOOST_AUTO_TEST_CASE(assume_tvpi_equalities) {
  test_domain_t dom1;
  dom1 += (x == z_number(1));
  dom1 += (y == x * z_number(2));
  dom1 += (z == x * z_number(3) + z_number(7));
  dom1 += (x * z_number(5) + z_number(6) * n == z_number(4));
  dom1 += (z * z_number(3) + y * z_number(6) + k * z_number(9) == z_number(15));
  BOOST_TEST(!dom1.is_bottom());
  // Coefficient 3 is in the template so 3z+6y+9k=15 scales to z+2y+3k=5;
  // with z=10 and y=2, k must equal -3.
  BOOST_TEST(dom1.entails(k == z_number(-3)));
}

BOOST_AUTO_TEST_CASE(exact_values_meet_is_bottom) {
  test_domain_t dom1, dom2;
  // dom1: x=1, y=2x, z=3x     dom2: x=2, y=3x, z=4x
  dom1.assign(x, z_number(1));
  dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
  dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));
  dom2 += (x == z_number(2));
  dom2 += (y == x * z_number(3));
  dom2 += (z == x + x + x + x);
  check_lattice_laws(dom1, dom2);
  BOOST_TEST((dom1 & dom2).is_bottom(), "meet must be bottom (x=1 vs x=2)");

  test_domain_t dom5(dom1);
  dom5 -= x;
  BOOST_TEST(!dom5.is_bottom());
  BOOST_TEST(dom5.entails(y == z_number(2)), "y survives forgetting x");
  BOOST_TEST(dom5.entails(z == z_number(3)), "z survives forgetting x");

  test_domain_t dom6(dom1);
  dom6.rename({x}, {i});
  BOOST_TEST(!dom6.is_bottom());
  BOOST_TEST(dom6.entails(i == z_number(1)), "i == 1 after renaming x to i");

  test_domain_t dom7(dom1);
  dom7.project({y, z});
  BOOST_TEST(!dom7.is_bottom());
  BOOST_TEST(dom7.entails(y == z_number(2)));
  BOOST_TEST(dom7.entails(z == z_number(3)));

  test_domain_t dom8(dom1);
  dom8.project({x});
  BOOST_TEST(!dom8.is_bottom());
  BOOST_TEST(dom8.entails(x == z_number(1)));
}

BOOST_AUTO_TEST_CASE(ranges_meet_not_bottom) {
  test_domain_t dom1, dom2;
  // dom1: x>=1, y=2x, z=3x     dom2: x<=20, i=3x, j=4x
  dom1 += (x >= z_number(1));
  dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
  dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));
  dom2 += (x <= z_number(20));
  dom2.assign(i, z_number(3) * x);
  dom2.assign(j, z_number(4) * x);
  check_lattice_laws(dom1, dom2);
  BOOST_TEST(!(dom1 & dom2).is_bottom(), "meet keeps x in [1,20]");

  test_domain_t dom5(dom1);
  dom5 -= x;
  BOOST_TEST(!dom5.is_bottom());
  BOOST_TEST(dom5.entails(y >= z_number(2)));
  BOOST_TEST(dom5.entails(z >= z_number(3)));

  test_domain_t dom6(dom1);
  dom6.rename({x}, {i});
  BOOST_TEST(!dom6.is_bottom());
  BOOST_TEST(dom6.entails(i >= z_number(1)));

  test_domain_t dom7(dom1);
  dom7.project({y, z});
  BOOST_TEST(!dom7.is_bottom());
  BOOST_TEST(dom7.entails(y >= z_number(2)));
  BOOST_TEST(dom7.entails(z >= z_number(3)));

  test_domain_t dom8(dom1);
  dom8.project({x});
  BOOST_TEST(!dom8.is_bottom());
  BOOST_TEST(dom8.entails(x >= z_number(1)));
}

BOOST_AUTO_TEST_CASE(widening_loop_counter) {
  test_domain_t dom1, dom2;
  dom1 += (x == z_number(0));
  dom1 += (i == z_number(0));
  BOOST_TEST(!dom1.is_bottom());
  dom2 = dom1;
  dom2.intrinsic("loop_counter", {i}, {});
  dom2.apply(OP_ADDITION, i, i, z_number(1));
  dom2.apply(OP_ADDITION, x, x, z_number(3));
  BOOST_TEST(!dom2.is_bottom());
  check_lattice_laws(dom1, dom2);
  test_domain_t dom5 = dom1 || (dom1 | dom2);
  BOOST_TEST(!dom5.is_bottom());
  BOOST_TEST((dom1 <= dom5), "widening result covers its left argument");
}

BOOST_AUTO_TEST_CASE(entailment_via_resultant) {
  // dom1: x>=5, y=4x+2, z=2x+3y, so y = 4x+2 <= 14x+6 = z for x >= 5.
  test_domain_t dom1;
  dom1 += (x >= z_number(5));
  dom1.assign(y, z_number(4) * x + z_number(2));
  dom1.assign(z, z_number(2) * x + z_number(3) * y);
  BOOST_TEST(dom1.entails(y <= z));

  // dom2: x in [1,10], y=4x+17, i in [0,x), z=4i+3:
  // i <= x-1 gives 4i+3 <= 4x-1 <= 4x+17 = y.
  test_domain_t dom2;
  dom2 += (x >= z_number(1));
  dom2 += (x <= z_number(10));
  dom2.assign(y, z_number(4) * x + z_number(17));
  dom2 += (i >= z_number(0));
  dom2 += (i <= x - z_number(1));
  dom2.assign(z, z_number(4) * i + z_number(3));
  BOOST_TEST(dom2.entails(z <= y));
}

BOOST_AUTO_TEST_CASE(lazy_leq_incompleteness) {
  // With lazy normalization the domain cannot derive x<=1 from
  // {x-y<=4, 2y-x<=-3, 3x-y<=5} without an explicit normalize(), so the
  // semantically valid inclusions below answer false.  This pins the
  // documented behavior of operator<= on unclosed lazy operands.
  test_domain_t dom1;
  dom1 += (x - y <= z_number(4));
  dom1 += (z_number(2) * y - x <= z_number(-3));
  dom1 += (z_number(3) * x - y <= z_number(5));
  BOOST_TEST(!dom1.is_bottom());
  test_domain_t dom2;
  dom2 += (y <= z_number(0));
  dom2 += (x <= z_number(2));
  dom2 += (y - x <= z_number(2));
  BOOST_TEST(!dom2.is_bottom());
  check_lattice_laws(dom1, dom2);
  test_domain_t dom4;
  dom4 += (y <= z_number(0));
  dom4 += (x <= z_number(2));
  dom4 += (z_number(2) * y - x <= z_number(2));
  BOOST_TEST(!(dom1 <= dom4), "lazy leq misses the saturated inclusion");
  BOOST_TEST(!(dom2 <= dom4), "lazy leq misses the saturated inclusion");
}

BOOST_AUTO_TEST_CASE(leq_refuted) {
  // Neither operand is semantically included in dom4 (e.g. (x=-1, y=1)
  // satisfies dom1 but violates y <= 0).
  test_domain_t dom1;
  dom1 += (z_number(2) * x - y <= z_number(5));
  dom1 += (z_number(3) * y - x <= z_number(7));
  dom1 += (-z_number(3) * x + y <= z_number(4));
  BOOST_TEST(!dom1.is_bottom());
  test_domain_t dom2;
  dom2 += (z_number(3) * x - z_number(2) * y <= z_number(8));
  dom2 += (z_number(4) * y - x <= z_number(10));
  dom2 += (z_number(2) * x - z_number(3) * y <= z_number(-6));
  BOOST_TEST(!dom2.is_bottom());
  check_lattice_laws(dom1, dom2);
  test_domain_t dom4;
  dom4 += (y <= z_number(0));
  dom4 += (x <= z_number(2));
  dom4 += (z_number(2) * y - x <= z_number(2));
  BOOST_TEST(!(dom1 <= dom4));
  BOOST_TEST(!(dom2 <= dom4));
}

BOOST_AUTO_TEST_CASE(array_access_bound) {
  z_var isz(vfac["t"], crab::INT_TYPE, 32);
  z_var len(vfac["l"], crab::INT_TYPE, 32);
  z_var tsz(vfac["s"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  z_var offset(vfac["off"], crab::INT_TYPE, 32);
  test_domain_t dom1;
  dom1 += (isz == z_number(4));
  set_range(dom1, len, 1, 10);
  dom1.apply(OP_MULTIPLICATION, tmp, isz, len);
  dom1 += (tmp <= tsz); // 4*len <= tsz
  dom1 += (i >= z_number(0));
  dom1 += (i <= len - 1);
  dom1.apply(OP_MULTIPLICATION, x, isz, i); // x = 4*i
  dom1.apply(OP_ADDITION, offset, x, isz);  // offset = 4*i + 4
  // offset = 4i+4 <= 4(len-1)+4 = 4len <= tsz
  BOOST_TEST(dom1.entails(offset <= tsz));
}

BOOST_AUTO_TEST_CASE(derived_bound_through_scaled_terms) {
  test_domain_t dom1;
  set_range(dom1, x, 1, 10);
  dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
  dom1 += (y <= z);
  dom1 += (i >= z_number(0));
  dom1 += (i <= x - z_number(1));
  dom1.apply(OP_MULTIPLICATION, n, i, z_number(2));
  // n = 2i <= 2(x-1) < 2x = y <= z
  BOOST_TEST(dom1.entails(n <= z));
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// Facts the lazy variant derives incrementally, while constraints are added.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(incremental_saturation, template_234)

BOOST_AUTO_TEST_CASE(bounds_prove_leq) {
  test_domain_t dom;
  set_range(dom, i, 0, 6);
  set_range(dom, j, 4, 10);
  dom += (x == z_number(3) * j); // x in [12, 30]
  dom += (y == z_number(2) * i); // y in [0, 12]
  BOOST_TEST(!dom.is_bottom());
  // max(y) = 12 = min(x): provable through bounds alone.
  BOOST_TEST(dom.entails(y <= x));
}

BOOST_AUTO_TEST_CASE(chain_needs_incremental_reduce) {
  test_domain_t dom;
  set_range(dom, i, 0, 6);
  set_range(dom, j, 0, 10);
  dom += (j > i);
  dom += (x == z_number(3) * j);
  dom += (y == z_number(2) * i);
  BOOST_TEST(!dom.is_bottom());
  // j >= i+1 gives 3j >= 3i+3 > 2i.  The lazy domain proves it because
  // incremental saturation derives ghost(i,2) - ghost(j,2) <= -2 when
  // j > i is added; historically this needed a full TvpiReduce.
  BOOST_TEST(dom.entails(y <= x));
}

BOOST_AUTO_TEST_CASE(partial_proofs) {
  test_domain_t dom;
  set_range(dom, i, 1, 4);
  set_range(dom, j, 1, 4);
  dom += (k >= z_number(0));
  dom += (k <= i - 1);
  dom += (x == z_number(3) * i);
  dom += (o == z_number(2) + k);
  dom += (y == z_number(2) * j);
  BOOST_TEST(!dom.is_bottom());
  // k <= i-1 gives k+2 <= i+1 <= 3i for i >= 1.
  BOOST_TEST(dom.entails(o <= x));
  dom += (z == o + y); // z = k + 2 + 2j
  // z can exceed x (i=1, k=0, j=4 gives z=10, x=3): must not be entailed.
  BOOST_TEST(!dom.entails(z <= x));
}

BOOST_AUTO_TEST_CASE(join_of_exact_values) {
  test_domain_t dom1, dom2;
  dom1 += (i == z_number(4));
  dom2 += (i == z_number(4));
  dom1 += (x == 2 * i); // x = 8
  dom2 += (x == 3 * i); // x = 12
  BOOST_TEST(!dom1.is_bottom());
  BOOST_TEST(!dom2.is_bottom());
  auto dom3 = dom1 | dom2;
  BOOST_TEST(!dom3.is_bottom());
  BOOST_TEST(dom3.entails(x >= 8));
  BOOST_TEST(dom3.entails(x <= 12));
}

BOOST_AUTO_TEST_CASE(widening_law) {
  test_domain_t dom1, dom2, dom3;
  dom1 += (x == z_number(0));
  dom1 += (i == z_number(0));
  BOOST_TEST(!dom1.is_bottom());
  dom2 = dom1;
  dom2.intrinsic("loop_counter", {i}, {});
  dom2.apply(OP_ADDITION, i, i, z_number(1));
  BOOST_TEST(!dom2.is_bottom());
  dom3 = dom2;
  dom2.apply(OP_ADDITION, x, x, z_number(3)); // x=3, i=1
  dom3.apply(OP_ADDITION, x, x, z_number(2)); // x=2, i=1
  BOOST_TEST(!dom3.is_bottom());
  test_domain_t dom4 = dom2 | dom3; // x in [2,3], i=1
  BOOST_TEST(!dom4.is_bottom());
  test_domain_t dom5 = dom1 || (dom1 | dom4); // x >= 0, i >= 0
  BOOST_TEST(!dom5.is_bottom());
  BOOST_TEST((dom1 <= dom5), "widening result covers its left argument");
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// The TVPI paper's C-string example over the {10,255} template.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(wide_template, template_10_255)

BOOST_AUTO_TEST_CASE(cstring_join_keeps_tvpi_differences) {
  // char s[32] = "the string"; i iterates until s[i] == 0:
  //   dom1: i in [0,9],  c in [1,255]   (inside the string)
  //   dom2: i == 10,     c == 0         (the terminator)
  //   dom3: i > 10                      (unreachable given i <= 10)
  test_domain_t dom1;
  set_range(dom1, i, 0, 9);
  set_range(dom1, c, 1, 255);
  BOOST_TEST(!dom1.is_bottom());
  test_domain_t dom2;
  dom2 += (i == z_number(10));
  dom2 += (c == z_number(0));
  BOOST_TEST(!dom2.is_bottom());
  test_domain_t dom3;
  dom3 += (i > z_number(10));
  set_range(dom3, c, 0, 255);
  BOOST_TEST(!dom3.is_bottom());
  test_domain_t input;
  set_range(input, i, 0, 10);
  BOOST_TEST((input & dom3).is_bottom(), "i <= 10 and i > 10 is unsat");

  test_domain_t output = input & dom1;
  output |= input & dom2;
  output |= input & dom3; // joins bottom: no effect

  BOOST_TEST(!output.is_bottom());
  BOOST_TEST(output.entails(i >= z_number(0)));
  BOOST_TEST(output.entails(i <= z_number(10)));
  BOOST_TEST(output.entails(c >= z_number(0)));
  BOOST_TEST(output.entails(c <= z_number(255)));
  // The TVPI difference constraints survive the join through the ghost
  // dimensions (the paper's sum form 255i + c <= 2550 is not a difference
  // and is out of reach by design):
  BOOST_TEST(output.entails(z_number(255) * i - c <= z_number(2550)),
             "255i - c <= 2550 via ghost(i,255) - ghost(c,1)");
  BOOST_TEST(output.entails(i - z_number(10) * c <= z_number(10)),
             "i - 10c <= 10 via ghost(i,1) - ghost(c,10)");
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// Lazy vs eager normalization.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(normalization, template_234)

BOOST_AUTO_TEST_CASE(both_params_prove_basic_fact) {
  // x = 4*i, n = i: both variants must prove x == 4*n (DbmClosure suffices).
  test_domain_t dom_def;
  set_range(dom_def, i, 0, 10);
  dom_def += (x == z_number(4) * i);
  dom_def.assign(n, i);
  BOOST_TEST(dom_def.entails(x == z_number(4) * n));

  eager_dom_t dom_eager;
  set_range(dom_eager, i, 0, 10);
  dom_eager += (x == z_number(4) * i);
  dom_eager.assign(n, i);
  BOOST_TEST(dom_eager.entails(x == z_number(4) * n));
}

BOOST_AUTO_TEST_CASE(manual_normalize_matches_eager) {
  // x = 4*i, y = 2*i: y <= x needs the cross-coefficient constraint.
  test_domain_t d_raw;
  set_range(d_raw, i, 0, 10);
  d_raw += (x == z_number(4) * i);
  d_raw += (y == z_number(2) * i);
  test_domain_t d_manual = d_raw;
  d_manual.normalize();
  bool r_manual = d_manual.entails(y <= x);
  BOOST_TEST(r_manual, "explicit normalize() must prove y <= x");

  eager_dom_t d_auto;
  set_range(d_auto, i, 0, 10);
  d_auto += (x == z_number(4) * i);
  d_auto += (y == z_number(2) * i);
  BOOST_TEST(d_auto.entails(y <= x) == r_manual,
             "the eager variant matches manual normalization");
}

BOOST_AUTO_TEST_CASE(array_bound_eager) {
  z_var isz(vfac["isz"], crab::INT_TYPE, 32);
  z_var len(vfac["len"], crab::INT_TYPE, 32);
  z_var tsz(vfac["tsz"], crab::INT_TYPE, 32);
  z_var offset(vfac["offset"], crab::INT_TYPE, 32);
  eager_dom_t dom;
  dom += (isz == z_number(4));
  set_range(dom, len, 1, 10);
  dom += (tsz >= z_number(4) * len);
  dom += (i >= z_number(0));
  dom += (i <= len - 1);
  dom += (x == z_number(4) * i);
  dom += (offset == x + z_number(4));
  BOOST_TEST(dom.entails(offset <= tsz));
}

BOOST_AUTO_TEST_CASE(normalize_idempotent) {
  eager_dom_t dom1;
  set_range(dom1, i, 0, 6);
  set_range(dom1, j, 4, 10);
  dom1 += (x == z_number(3) * j);
  dom1 += (y == z_number(2) * i);
  eager_dom_t dom2 = dom1;
  dom2.normalize();
  BOOST_TEST((dom1 <= dom2));
  BOOST_TEST((dom2 <= dom1));
}

BOOST_AUTO_TEST_CASE(bottom_detection) {
  eager_dom_t dom;
  dom += (i >= z_number(5));
  dom += (i <= z_number(3));
  BOOST_TEST(dom.is_bottom());
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// The saturation rules, one derivation at a time.  "incremental" is the fact
// being materialized right after the additions, with no normalize().
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(reduction, template_234)

BOOST_AUTO_TEST_CASE(purple_new_edge_first_argument) {
  //   3y - z <= 0  and  [new] 2x - 2y <= 0   gives  3x - z <= 0
  test_domain_t d;
  d += (3 * y - z <= 0);
  d += (2 * x - 2 * y <= 0);
  z_lin_cst_t t(3 * x - z <= 0);
  BOOST_TEST(materialized(d, t), "incremental must derive 3x - z <= 0");
  BOOST_TEST(materialized_after_normalize(d, t));
}

BOOST_AUTO_TEST_CASE(orange_new_edge_second_argument) {
  //   w - 3x <= 0  and  [new] 2x - 2y <= 0   gives  w - 3y <= 0
  test_domain_t d;
  d += (w - 3 * x <= 0);
  d += (2 * x - 2 * y <= 0);
  z_lin_cst_t t(w - 3 * y <= 0);
  BOOST_TEST(materialized(d, t), "incremental must derive w - 3y <= 0");
  BOOST_TEST(materialized_after_normalize(d, t));
}

BOOST_AUTO_TEST_CASE(green_with_aligned_seed) {
  //   w - 2x <= 0, 3y - z <= 0  and  [new] 2x - 2y <= 0  gives 3w - 2z <= 0
  test_domain_t d;
  d += (w - 2 * x <= 0);
  d += (3 * y - z <= 0);
  d += (2 * x - 2 * y <= 0); // the trigger
  z_lin_cst_t t(3 * w - 2 * z <= 0);
  BOOST_TEST(materialized(d, t), "seeded Green must derive 3w - 2z <= 0");
  BOOST_TEST(materialized_after_normalize(d, t));
}

BOOST_AUTO_TEST_CASE(green_interior) {
  //   w - 3x <= 0, 2y - z <= 0  and  [new] 2x - 2y <= 0  gives 2w - 3z <= 0
  test_domain_t d;
  d += (w - 3 * x <= 0);
  d += (2 * y - z <= 0);
  d += (2 * x - 2 * y <= 0); // the trigger
  z_lin_cst_t t(2 * w - 3 * z <= 0);
  BOOST_TEST(materialized(d, t), "Green must derive 2w - 3z <= 0");
  BOOST_TEST(materialized_after_normalize(d, t));
}

BOOST_AUTO_TEST_CASE(v0_bound_with_floor) {
  //   y <= 5  and  [new] 2x - 3y <= 0   gives  2x <= 15, so x <= 7
  test_domain_t d;
  d += (y <= z_number(5));
  d += (2 * x - 3 * y <= 0);
  auto ub_incr = d[x].ub().number();
  BOOST_TEST((ub_incr && *ub_incr == 7), "incremental v0 bound floors 15/2");
  test_domain_t f(d);
  f.normalize();
  auto ub_full = f[x].ub().number();
  BOOST_TEST((ub_full && *ub_full == 7));
}

BOOST_AUTO_TEST_CASE(bound_insertion_gap) {
  // DOCUMENTED GAP: bound insertions do not trigger the incremental pass
  // (a per-bound alignment sweep would cost O(|T|^2 N) on the most common
  // constraint kind).  normalize() derives the fact.  If the incremental
  // check below starts passing, the gap was closed: revisit this test.
  test_domain_t d;
  d += (2 * y - x <= 0);
  d += (x <= z_number(7));
  auto ub_incr = d[y].ub().number();
  BOOST_TEST(!(ub_incr && *ub_incr <= 3), "expected incremental gap");
  test_domain_t f(d);
  f.normalize();
  auto ub_full = f[y].ub().number();
  BOOST_TEST((ub_full && *ub_full <= 3), "full reduce must derive y <= 3");
}

BOOST_AUTO_TEST_CASE(insertion_order_robustness) {
  test_domain_t d;
  d += (2 * x - 2 * y <= 0);
  d += (3 * y - z <= 0);
  d += (w - 2 * x <= 0); // trigger arrives last, from the other side
  z_lin_cst_t t(3 * w - 2 * z <= 0);
  // Order-dependence of the incremental pass is a documented property, so
  // only the full reduce is pinned here.
  BOOST_TEST(materialized_after_normalize(d, t),
             "full reduce must be insertion-order insensitive");
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// Regressions (one case per soundness bug found in review) and coverage for
// the APIs no other suite exercises.  The template deliberately contains 1:
// ghost(v,1) is v itself, so every ghost loop that failed to skip
// coefficient 1 operated on the program variable (bug A4).
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(regressions, template_1234)

BOOST_AUTO_TEST_CASE(a1_self_referential_assign) {
  // The ghost refinement used to rewrite x := x + 2*y into the meet
  // x == x + ghost(y,2); the two x terms cancel and the state collapsed
  // via ghost(y,2) == 0.
  test_domain_t dom;
  dom.assign(y, z_number(1));
  dom.assign(x, z_number(5));
  dom.assign(x, x + z_number(2) * y);
  BOOST_TEST(!dom.is_bottom(), "x := x + 2*y must not collapse to bottom");
  BOOST_TEST(dom.entails(x == z_number(7)));
  BOOST_TEST(dom.entails(y == z_number(1)));
  BOOST_TEST(!dom.entails(z_number(2) * y <= z_number(0)));

  test_domain_t dom2; // negative variant
  dom2.assign(y, z_number(1));
  dom2.assign(x, z_number(5));
  dom2.assign(x, x - z_number(2) * y);
  BOOST_TEST(!dom2.is_bottom());
  BOOST_TEST(dom2.entails(x == z_number(3)));

  test_domain_t dom3; // control: non-unit self-reference keeps working
  dom3.assign(x, z_number(3));
  dom3.assign(x, z_number(2) * x);
  BOOST_TEST(!dom3.is_bottom());
  BOOST_TEST(dom3.entails(x == z_number(6)));

  test_domain_t dom4; // control: all-unit rhs (identity rewrite path)
  dom4.assign(y, z_number(1));
  dom4.assign(x, z_number(5));
  dom4.assign(x, x + y);
  BOOST_TEST(!dom4.is_bottom());
  BOOST_TEST(dom4.entails(x == z_number(6)));
}

BOOST_AUTO_TEST_CASE(a3_aliased_var_apply) {
  // The singleton reduction used to run AFTER the base operation and then
  // re-enter the scalar apply: the op was applied twice (x aliases y), or
  // the singleton was read from the post-state (x aliases z).
  test_domain_t dom;
  dom.assign(x, z_number(5));
  dom.assign(z, z_number(2));
  dom.apply(OP_ADDITION, x, x, z);
  BOOST_TEST(dom.entails(x == z_number(7)), "x == 7, not 9: op applied once");

  test_domain_t dom2;
  dom2.assign(y, z_number(3));
  dom2.assign(x, z_number(2));
  dom2.apply(OP_SUBTRACTION, x, y, x);
  BOOST_TEST(dom2.entails(x == z_number(1)), "singleton of x read pre-state");

  test_domain_t dom3;
  dom3.assign(y, z_number(2));
  dom3.assign(x, z_number(5));
  dom3.apply(OP_MULTIPLICATION, x, y, x); // y singleton, commutative
  BOOST_TEST(dom3.entails(x == z_number(10)), "x == 10, not 20");

  test_domain_t dom4; // precision lock-in: non-aliased singleton reduction
  set_range(dom4, w, 0, 10);
  dom4.assign(z, z_number(2));
  dom4.apply(OP_MULTIPLICATION, x, w, z);
  BOOST_TEST(dom4.entails(x == z_number(2) * w), "x == 2*w survives");
  BOOST_TEST(dom4.entails(x >= z_number(0)));
  BOOST_TEST(dom4.entails(x <= z_number(20)));
}

BOOST_AUTO_TEST_CASE(a4_coefficient_one_in_template) {
  // Ghost loops that did not skip c == 1 operated on the variable itself:
  // rename and expand passed duplicate pairs into the base domain (fatal),
  // and the apply loops rewrote the variable in place.
  test_domain_t dom;
  dom.assign(x, z_number(2));
  dom.apply(OP_MULTIPLICATION, y, x, z_number(3));
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(dom.entails(y == z_number(6)));

  test_domain_t dp(dom);
  dp.project({x, y});
  BOOST_TEST(dp.entails(x == z_number(2)));
  BOOST_TEST(dp.entails(y == z_number(6)));

  test_domain_t dr(dom);
  dr.rename({x}, {i});
  BOOST_TEST(dr.entails(i == z_number(2)));

  test_domain_t de(dom);
  de.expand(x, j);
  BOOST_TEST(de.entails(j == z_number(2)));
  BOOST_TEST(de.entails(x == z_number(2)));

  test_domain_t da;
  da.assign(u, z_number(4));
  da.assign(w, z_number(5));
  da.apply(OP_ADDITION, k, u, w); // the var-RHS loop must skip c == 1 too
  BOOST_TEST(da.entails(k == z_number(9)));
}

BOOST_AUTO_TEST_CASE(reuse_after_set_to_top) {
  // The incremental saturation hook must still fire on a reset value
  // (m_reduce_active must not be stuck: RAII guard regression).
  test_domain_t dom;
  dom.assign(x, z_number(1));
  dom.assign(y, z_number(2));
  dom.normalize();
  dom.set_to_top();
  BOOST_TEST(dom.is_top());
  dom += (z_number(3) * y - z <= z_number(0));
  dom += (z_number(2) * x - z_number(2) * y <= z_number(0));
  BOOST_TEST(materialized(dom, z_number(3) * x - z <= z_number(0)),
             "incremental saturation derives 3x - z <= 0 without normalize");
  BOOST_TEST(!dom.entails(x == z_number(1)), "old facts do not resurrect");
}

BOOST_AUTO_TEST_CASE(lattice_constants_and_var_lifecycle) {
  test_domain_t dom;
  set_range(dom, x, 1, 3);
  BOOST_TEST(!dom.is_top());
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(dom.make_top().is_top());
  BOOST_TEST(dom.make_bottom().is_bottom());
  BOOST_TEST(dom.make_bottom().at(x).is_bottom());
  const test_domain_t &cdom = dom;
  auto itv = cdom.at(x);
  BOOST_TEST((itv.lb().number() && *itv.lb().number() == z_number(1)));
  BOOST_TEST((itv.ub().number() && *itv.ub().number() == z_number(3)));
  dom.set_to_bottom();
  BOOST_TEST(dom.is_bottom());
  dom.set_to_top();
  BOOST_TEST(dom.is_top());
  dom -= x; // forget on top: no crash, stays top
  BOOST_TEST(dom.is_top());
  dom += (x == z_number(4)); // reusable after the resets
  BOOST_TEST(dom.entails(x == z_number(4)));
}

BOOST_AUTO_TEST_CASE(weak_assign_hull) {
  test_domain_t dom;
  dom.assign(x, z_number(1));
  dom.assign(y, z_number(5));
  dom.weak_assign(y, z_number(2) * x); // y may stay 5 or become 2
  BOOST_TEST(!dom.entails(y == z_number(5)));
  BOOST_TEST(!dom.entails(y == z_number(2)));
  BOOST_TEST(dom.entails(y >= z_number(2)));
  BOOST_TEST(dom.entails(y <= z_number(5)));
  BOOST_TEST(dom.entails(x == z_number(1)));
}

BOOST_AUTO_TEST_CASE(select_default) {
  test_domain_t dom;
  set_range(dom, x, 0, 10);
  dom.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
  BOOST_TEST(dom.entails(y >= z_number(0)));
  BOOST_TEST(dom.entails(y <= z_number(30)));

  test_domain_t dom2;
  dom2.assign(x, z_number(2)); // condition entailed
  dom2.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
  BOOST_TEST(dom2.entails(y == z_number(4)));

  test_domain_t dom3;
  dom3.assign(x, z_number(7)); // condition refuted
  dom3.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
  BOOST_TEST(dom3.entails(y == z_number(21)));
}

BOOST_AUTO_TEST_CASE(int_conversions) {
  test_domain_t dom;
  dom.assign(y, z_number(7));
  dom.apply(crab::domains::OP_SEXT, s64, y);
  BOOST_TEST(dom.entails(s64 == z_number(7)));
  dom.apply(OP_MULTIPLICATION, w, y, z_number(2)); // w == 2*y
  BOOST_TEST(dom.entails(w == z_number(2) * y));
  dom.apply(crab::domains::OP_TRUNC, w, s64); // redefines w: ghosts must go
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(!dom.entails(w == z_number(2) * y), "stale w == 2*y dropped");
}

BOOST_AUTO_TEST_CASE(bitwise_conservative_no_stale) {
  test_domain_t dom;
  dom.assign(y, z_number(12));
  dom.assign(z, z_number(10));
  dom.apply(OP_MULTIPLICATION, v, y, z_number(2)); // v == 2*y
  BOOST_TEST(dom.entails(v == z_number(2) * y));
  dom.apply(crab::domains::OP_AND, x, y, z); // variable overload
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(dom.entails(x >= z_number(0)));
  BOOST_TEST(dom.entails(x <= z_number(12)));
  dom.apply(crab::domains::OP_AND, w, y, z_number(4)); // scalar overload
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(dom.entails(w >= z_number(0)));
  BOOST_TEST(dom.entails(w <= z_number(12)));
  dom.apply(crab::domains::OP_AND, v, y, z_number(3)); // redefines v
  BOOST_TEST(!dom.entails(v == z_number(2) * y), "stale v == 2*y dropped");
}

BOOST_AUTO_TEST_CASE(forget_vector_and_meet_assign) {
  test_domain_t dom;
  dom.assign(x, z_number(1));
  dom.assign(y, z_number(2));
  dom.assign(z, z_number(3));
  dom.forget({x, y});
  BOOST_TEST(dom.entails(z == z_number(3)));
  BOOST_TEST(dom.at(x).is_top());
  BOOST_TEST(dom.at(y).is_top());
  dom.forget({}); // no-op
  BOOST_TEST(dom.entails(z == z_number(3)));

  test_domain_t d1, d2;
  set_range(d1, x, 0, 10);
  set_range(d2, x, 5, 20);
  d1 &= d2;
  BOOST_TEST(d1.entails(x >= z_number(5)));
  BOOST_TEST(d1.entails(x <= z_number(10)));
  test_domain_t d3, d4;
  d3 += (x == z_number(1));
  d4 += (x == z_number(2));
  d3 &= d4;
  BOOST_TEST(d3.is_bottom(), "meet of disjoint values is bottom");
  test_domain_t d5;
  d5 += (x == z_number(1));
  d5 &= d5.make_top();
  BOOST_TEST(d5.entails(x == z_number(1)), "meet with top is identity");
  d5 &= d5.make_bottom();
  BOOST_TEST(d5.is_bottom());
}

BOOST_AUTO_TEST_CASE(to_disjunctive_constraints) {
  test_domain_t dom;
  set_range(dom, x, 1, 3);
  dom.apply(OP_MULTIPLICATION, y, x, z_number(2));
  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_TEST(!dcsts.is_false());
  BOOST_TEST(!dcsts.is_true());
  BOOST_TEST(
      dom.make_bottom().to_disjunctive_linear_constraint_system().is_false());
  BOOST_TEST(
      dom.make_top().to_disjunctive_linear_constraint_system().is_true());
}

BOOST_AUTO_TEST_CASE(nonsingleton_aliasing) {
  // The fixed paths under general intervals: try_scalar_apply does not fire
  // and the plain var-RHS path runs, with aliasing handled by the base
  // domain's pre-state reads.
  test_domain_t dom;
  set_range(dom, x, 1, 5);
  set_range(dom, z, 1, 2);
  dom.apply(OP_ADDITION, x, x, z); // aliased, both operands intervals
  BOOST_TEST(dom.entails(x >= z_number(2)));
  BOOST_TEST(dom.entails(x <= z_number(7)));
  BOOST_TEST(!dom.entails(x <= z_number(6)), "upper bound is tight");

  test_domain_t dom2; // relational input, relational output
  set_range(dom2, w, 0, 10);
  dom2.assign(y, w + z_number(3)); // y - w == 3
  set_range(dom2, z, 1, 2);
  dom2.apply(OP_ADDITION, x, y, z); // all non-singleton
  BOOST_TEST(dom2.entails(x - w >= z_number(4)));
  BOOST_TEST(dom2.entails(x - w <= z_number(5)));
  BOOST_TEST(dom2.entails(x >= z_number(4)));

  test_domain_t dom3; // the A1 shape with an interval operand
  set_range(dom3, y, 1, 3);
  dom3.assign(x, z_number(5));
  dom3.assign(x, x + z_number(2) * y);
  BOOST_TEST(!dom3.is_bottom());
  BOOST_TEST(dom3.entails(x >= z_number(7)));
  BOOST_TEST(dom3.entails(x <= z_number(11)));
}

BOOST_AUTO_TEST_CASE(weak_self_reference) {
  // The weak path skips the ghost refinement, so the A1 bug cannot fire
  // here; this pins that weak semantics stay sound under aliasing.
  test_domain_t dom;
  set_range(dom, y, 1, 3);
  dom.assign(x, z_number(5));
  dom.weak_assign(x, x + z_number(2) * y); // x may stay 5 or become [7,11]
  BOOST_TEST(!dom.is_bottom());
  BOOST_TEST(dom.entails(x >= z_number(5)));
  BOOST_TEST(dom.entails(x <= z_number(11)));
  BOOST_TEST(!dom.entails(x == z_number(5)));
}

BOOST_AUTO_TEST_CASE(a6_stale_ghost_hook_window) {
  // The incremental saturation hook in the scalar MUL special used to run
  // BEFORE x's ghost family was re-established: eliminate_x combined the
  // still-stale edge w - ghost(x,2) <= 5 (about x_old) with the new
  // x == 2*y (about x_new), deriving w <= 4*y + 5 = 9 against w >= 20.
  // x and w must be genuine intervals: with singletons the base domain
  // subsumes the relational edge into bounds and the hazard is masked.
  test_domain_t dom;
  set_range(dom, x, 8, 10);
  set_range(dom, w, 20, 25);
  dom += (w - z_number(2) * x <= z_number(5)); // satisfiable: x=10, w=25
  dom.assign(y, z_number(1));
  dom.apply(OP_MULTIPLICATION, x, y, z_number(2)); // x := 2*y == 2
  BOOST_TEST(!dom.is_bottom(), "reachable state must not collapse");
  BOOST_TEST(dom.entails(x == z_number(2)));
  BOOST_TEST(dom.entails(w >= z_number(20)), "old bound not misused");
  BOOST_TEST(!dom.entails(w <= z_number(9)),
             "no derivation may conflate x_old with x_new");
}

BOOST_AUTO_TEST_CASE(weak_assign_relational) {
  // weak y := 2*x from y == 10 joins {y == 10} with {y == 2*x}; both
  // branches satisfy 2*x <= y, so the relational bound must survive.
  test_domain_t dom;
  set_range(dom, x, 1, 3);
  dom.assign(y, z_number(10));
  dom.weak_assign(y, z_number(2) * x);
  BOOST_TEST(!dom.entails(y == z_number(10)));
  BOOST_TEST(!dom.entails(y == z_number(2) * x));
  BOOST_TEST(dom.entails(y >= z_number(2)));
  BOOST_TEST(dom.entails(y <= z_number(10)));
  BOOST_TEST(dom.entails(z_number(2) * x <= y),
             "2*x <= y holds in both branches of the weak join");
}

BOOST_AUTO_TEST_CASE(ghost_family_exactness) {
  { // assign, non-unit self-reference plus a second term:
    // x==5, y==3, x := 2*x + y: x_new = 13 and ghost(x,2) = 26.
    test_domain_t d;
    d.assign(x, z_number(5));
    d.assign(y, z_number(3));
    d.assign(x, z_number(2) * x + y);
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x == z_number(13)));
    BOOST_TEST(d.entails(z_number(2) * x - y == z_number(23)),
               "2*x - y == 23 through ghost(x,2) = 26");
  }
  { // scalar MUL self-reference over an interval, with a bystander:
    // x in [2,5], w := x, x := x*3: x in [6,15], ghost(x,2) = 2*x.
    test_domain_t d;
    set_range(d, x, 2, 5);
    d.assign(w, x);
    d.apply(OP_MULTIPLICATION, x, x, z_number(3));
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x >= z_number(6)));
    BOOST_TEST(d.entails(x <= z_number(15)));
    BOOST_TEST(!d.entails(x <= z_number(14)), "upper bound is tight");
    BOOST_TEST(d.entails(w >= z_number(2)));
    BOOST_TEST(d.entails(w <= z_number(5)));
  }
  { // var-RHS ADD, aliased, both operands intervals.
    test_domain_t d;
    set_range(d, x, 1, 3);
    set_range(d, z, 1, 2);
    d.apply(OP_ADDITION, x, x, z);
    BOOST_TEST(d.entails(x >= z_number(2)));
    BOOST_TEST(d.entails(x <= z_number(5)));
  }
  { // chained composition: identity rewrite, then non-unit self-reference,
    // then mixed unit and non-unit: x = 1 -> 2 -> 4 -> 8, ghost(x,2) = 16.
    test_domain_t d;
    d.assign(x, z_number(1));
    d.assign(y, z_number(2));
    d.assign(x, x + z_number(1));
    d.assign(x, z_number(2) * x);
    d.assign(x, x + z_number(2) * y);
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x == z_number(8)));
    BOOST_TEST(d.entails(z_number(2) * x - y == z_number(14)),
               "2*x - y == 14 through ghost(x,2) = 16");
  }
}

BOOST_AUTO_TEST_CASE(eager_oracle_replay) {
  // The eager variant normalizes after every operation, meeting each ghost
  // family member against the others: a stale or wrong ghost value tends to
  // collapse the state.  The exactness sequences must replay bottom-free.
  {
    eager_dom_t d;
    d.assign(x, z_number(5));
    d.assign(y, z_number(3));
    d.assign(x, z_number(2) * x + y);
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x == z_number(13)));
  }
  {
    eager_dom_t d;
    set_range(d, x, 2, 5);
    d.assign(w, x);
    d.apply(OP_MULTIPLICATION, x, x, z_number(3));
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x >= z_number(6)));
    BOOST_TEST(d.entails(x <= z_number(15)));
    BOOST_TEST(d.entails(w <= z_number(5)));
  }
  {
    eager_dom_t d;
    d.assign(x, z_number(1));
    d.assign(y, z_number(2));
    d.assign(x, x + z_number(1));
    d.assign(x, z_number(2) * x);
    d.assign(x, x + z_number(2) * y);
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x == z_number(8)));
  }
  {
    eager_dom_t d; // A1 and A3 shapes under the eager variant
    d.assign(y, z_number(1));
    d.assign(x, z_number(5));
    d.assign(x, x + z_number(2) * y);
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(x == z_number(7)));
    eager_dom_t d2;
    d2.assign(x, z_number(5));
    d2.assign(z, z_number(2));
    d2.apply(OP_ADDITION, x, x, z);
    BOOST_TEST(d2.entails(x == z_number(7)));
  }
}

BOOST_AUTO_TEST_CASE(a7_aliased_division) {
  // The DIV special used to write ghost(x,|z|) := x AFTER the base apply
  // (post-state), and the rewrite loop then divided that again: y == 7,
  // y := y/2 left every ghost of y wrong (ghost(y,2) == 1 while 2*y == 6).
  // Ghosts are now refreshed from their definition instead.
  test_domain_t d;
  d.assign(y, z_number(7));
  d.apply(OP_SDIV, y, y, z_number(2));
  BOOST_TEST(d.entails(y == z_number(3)));
  BOOST_TEST(!d.entails(z_number(2) * y <= z_number(5)),
             "2*y <= 5 must not be entailed (2*y == 6)");
  d.normalize();
  BOOST_TEST(!d.is_bottom());

  eager_dom_t d2;
  d2.assign(y, z_number(7));
  d2.apply(OP_SDIV, y, y, z_number(2));
  BOOST_TEST(!d2.is_bottom());
  BOOST_TEST(d2.entails(y == z_number(3)));
}

BOOST_AUTO_TEST_CASE(division_bands) {
  // x := y / z with a = |z| satisfies y = sign(z)*a*x + r, |r| <= a-1,
  // sign(r) = sign(y); the band replaces the old, false equality
  // a*x == sign(z)*y.
  { // one-sided band and both ghost pairs, y in [10,21]
    test_domain_t d;
    set_range(d, y, 10, 21);
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(d.entails(x >= z_number(5)));
    BOOST_TEST(d.entails(x <= z_number(10)));
    BOOST_TEST(!d.entails(z_number(2) * x >= y), "y may be odd");
    BOOST_TEST(d.entails(z_number(2) * x <= y));
    BOOST_TEST(d.entails(y - z_number(2) * x <= z_number(1)));
    BOOST_TEST(d.entails(z_number(2) * y - z_number(4) * x <= z_number(3)),
               "band at the ghost pair (2,4)");
  }
  { // odd singleton: the old equality claimed 2*x == 21
    test_domain_t d;
    d.assign(y, z_number(21));
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(d.entails(x == z_number(10)));
    BOOST_TEST(!d.entails(z_number(2) * x >= z_number(21)), "2*x == 20");
  }
  { // a does not divide c: the interval fallback must be exact
    test_domain_t d;
    d.assign(y, z_number(7));
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(d.entails(x == z_number(3)));
    BOOST_TEST(!d.entails(z_number(3) * x >= z_number(10)),
               "3*x == 9; the old ghost said 10");
    BOOST_TEST(d.entails(z_number(3) * x >= z_number(8)));
  }
  { // negative divisor, odd singleton: the old ghost said 2*x == -21
    test_domain_t d;
    d.assign(y, z_number(21));
    d.apply(OP_SDIV, x, y, z_number(-2));
    BOOST_TEST(d.entails(x == z_number(-10)));
    BOOST_TEST(!d.entails(z_number(2) * x <= z_number(-21)), "2*x == -20");
    BOOST_TEST(d.entails(z_number(2) * x <= z_number(-19)));
  }
  { // unknown sign of y: two-sided band, and no false 2*x <= y
    test_domain_t d;
    set_range(d, y, -5, 9);
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(!d.entails(z_number(2) * x <= y),
               "y = -5 gives x = -2 and 2x = -4 > -5");
    BOOST_TEST(d.entails(y - z_number(2) * x <= z_number(1)));
    BOOST_TEST(d.entails(z_number(2) * x - y <= z_number(1)));
  }
  { // UDIV over entailed-nonnegative y behaves like SDIV
    test_domain_t d;
    set_range(d, y, 10, 21);
    d.apply(OP_UDIV, x, y, z_number(2));
    BOOST_TEST(d.entails(z_number(2) * x <= y));
    BOOST_TEST(!d.entails(z_number(2) * x >= y));
  }
  { // eager replay: the band must survive (and feed) full normalization
    eager_dom_t d;
    set_range(d, y, 10, 21);
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(z_number(2) * x <= y));
  }
}

BOOST_AUTO_TEST_CASE(division_branch_coverage) {
  { // SDIV, y <= 0, z > 0: the band flips one-sided the other way
    test_domain_t d;
    set_range(d, y, -21, -10);
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(d.entails(x >= z_number(-10)));
    BOOST_TEST(d.entails(x <= z_number(-5)));
    BOOST_TEST(d.entails(y - z_number(2) * x <= z_number(0)), "y <= 2*x");
    BOOST_TEST(d.entails(z_number(2) * x - y <= z_number(1)));
    BOOST_TEST(!d.entails(z_number(2) * x <= y),
               "y = -21 gives x = -10 and 2x = -20 > -21");

    test_domain_t ds; // odd singleton
    ds.assign(y, z_number(-21));
    ds.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(ds.entails(x == z_number(-10)));
    BOOST_TEST(!ds.entails(z_number(2) * x <= z_number(-21)), "2*x == -20");
  }
  { // SDIV, y <= 0, z < 0: sum band, interval content only
    test_domain_t d;
    d.assign(y, z_number(-21));
    d.apply(OP_SDIV, x, y, z_number(-2));
    BOOST_TEST(d.entails(x == z_number(10)));
    BOOST_TEST(!d.entails(z_number(2) * x >= z_number(21)), "2*x == 20");
  }
  { // UDIV with a possibly-negative dividend: no band may be added
    test_domain_t d;
    set_range(d, y, -5, 9);
    d.apply(OP_UDIV, x, y, z_number(2));
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(!d.entails(y - z_number(2) * x <= z_number(1)),
               "no band for sign-unknown udiv");
  }
  { // UDIV with a negative divisor: no band either
    test_domain_t d;
    set_range(d, y, 10, 21);
    d.apply(OP_UDIV, x, y, z_number(-2));
    BOOST_TEST(!d.is_bottom());
  }
  { // var-RHS DIV with a non-singleton divisor: stale facts must drop
    test_domain_t d;
    d.assign(w, z_number(10));
    d.apply(OP_MULTIPLICATION, v, w, z_number(2)); // v == 2*w
    BOOST_TEST(d.entails(v == z_number(2) * w));
    set_range(d, z, 2, 3);
    d.apply(OP_SDIV, v, w, z); // redefines v; divisor not singleton
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(!d.entails(v == z_number(2) * w), "stale v == 2*w dropped");
    BOOST_TEST(d.entails(v >= z_number(3)));
    BOOST_TEST(d.entails(v <= z_number(5)));
  }
  { // eager replay of the mirrored band
    eager_dom_t d;
    set_range(d, y, -21, -10);
    d.apply(OP_SDIV, x, y, z_number(2));
    BOOST_TEST(!d.is_bottom());
    BOOST_TEST(d.entails(y - z_number(2) * x <= z_number(0)));
  }
}

BOOST_AUTO_TEST_SUITE_END()

//===----------------------------------------------------------------------===//
// Division fallbacks the {1,2,3,4} template cannot reach.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(division_template_gaps, template_16)

BOOST_AUTO_TEST_CASE(quotient_outside_template) {
  // a | c but q = c/a = 3 is not in {1,6}: interval fallback only.
  test_domain_t d;
  set_range(d, y, 10, 21);
  d.apply(OP_SDIV, x, y, z_number(2));
  BOOST_TEST(d.entails(x >= z_number(5)));
  BOOST_TEST(d.entails(x <= z_number(10)));
  BOOST_TEST(d.entails(z_number(6) * x >= z_number(29)),
             "fallback ghost(x,6) = [30,60]");
  BOOST_TEST(!d.entails(y - z_number(2) * x <= z_number(1)),
             "no band: neither coefficient 2 nor q = 3 is tracked");
}

BOOST_AUTO_TEST_CASE(divisor_outside_template) {
  // a = 5 is not in {1,6}: every coefficient hits the a-does-not-divide-c
  // fallback.
  test_domain_t d;
  set_range(d, y, 10, 21);
  d.apply(OP_SDIV, x, y, z_number(5));
  BOOST_TEST(d.entails(x >= z_number(2)));
  BOOST_TEST(d.entails(x <= z_number(4)));
  BOOST_TEST(d.entails(z_number(6) * x <= z_number(25)),
             "fallback ghost(x,6) = [12,24]");
}

BOOST_AUTO_TEST_SUITE_END()

// Keep stdout empty: the golden-output harness (tests/run_tests.sh) diffs
// the standard output of every test binary, so the whole Boost.Test log is
// routed to stderr.  ctest passes --disable-warnings, which Boost.Test would
// reject, so it is translated into the crab flag it stands for.
int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int idx = 0; idx < argc; ++idx) {
    if (std::string(argv[idx]) == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
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
