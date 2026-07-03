// Functional audit of the two saturation procedures:
//   - incremental_tvpi_reduce (runs inside operator+= per TVPI constraint)
//   - tvpi_reduce             (runs inside normalize())
// Each case is a concrete Resultant derivation.  "incr" checks that the
// constraint is MATERIALIZED in the lazy state right after the additions
// (no normalize); "full" checks it after an explicit normalize().
#include "../../common.hpp"
#include "../../program_options.hpp"

#include <cassert>

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using test_domain_t = z_tvpi_dbm_domain_t;

// Is `target` (an inequality e <= 0) present, verbatim after ghost
// translation, in dom's constraint system?
static bool materialized(test_domain_t &dom, const z_lin_cst_t &target) {
  for (auto const &c : dom.to_linear_constraint_system()) {
    if (!c.is_inequality()) {
      continue;
    }
    auto diff = c.expression() - target.expression();
    if (diff.is_constant() && diff.constant() == 0) {
      return true;
    }
  }
  return false;
}

static bool materialized_after_normalize(const test_domain_t &dom,
                                         const z_lin_cst_t &target) {
  test_domain_t tmp(dom);
  tmp.normalize();
  return materialized(tmp, target);
}

static void report(const char *name, bool incr, bool full) {
  crab::outs() << name << ": incremental=" << (incr ? "yes" : "NO")
               << " full-reduce=" << (full ? "yes" : "NO") << "\n";
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  auto &coeffs = crab_domain_params_man::get().coefficients();
  coeffs.insert(coeffs.end(), {2, 3, 4});

  variable_factory_t vfac;
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  { // Case A — Purple (Case 1): new edge as FIRST Resultant argument.
    //   3y - z <= 0  ∧  [new] 2x - 2y <= 0
    //   λ1=3, λ2=2:  6x - 2z <= 0  --scale-->  3x - z <= 0
    test_domain_t d;
    d += (3 * y - z <= 0);
    d += (2 * x - 2 * y <= 0);
    z_lin_cst_t t(3 * x - z <= 0);
    bool incr = materialized(d, t), full = materialized_after_normalize(d, t);
    report("A purple  3x-z<=0   ", incr, full);
    assert(incr && "A: incremental must derive 3x-z<=0");
    assert(full && "A: full reduce must derive 3x-z<=0");
  }

  { // Case B — Orange (Case 2): new edge as SECOND Resultant argument.
    //   w - 3x <= 0  ∧  [new] 2x - 2y <= 0
    //   λ1=2, λ2=3:  2w - 6y <= 0  --scale-->  w - 3y <= 0
    test_domain_t d;
    d += (w - 3 * x <= 0);
    d += (2 * x - 2 * y <= 0);
    z_lin_cst_t t(w - 3 * y <= 0);
    bool incr = materialized(d, t), full = materialized_after_normalize(d, t);
    report("B orange  w-3y<=0   ", incr, full);
    assert(incr && "B: incremental must derive w-3y<=0");
    assert(full && "B: full reduce must derive w-3y<=0");
  }

  { // Case C — Green with the q==a SEED: 3-hop through both x and y where
    // the first hop is already aligned (base closure derives w-2y<=0, then
    // the seeded W_y entry aligns 2 vs 3):
    //   w - 2x <= 0  ∧  [new] 2x - 2y <= 0  ∧  3y - z <= 0
    //   ==>  3w - 2z <= 0
    test_domain_t d;
    d += (w - 2 * x <= 0);
    d += (3 * y - z <= 0);
    d += (2 * x - 2 * y <= 0); // the trigger
    z_lin_cst_t t(3 * w - 2 * z <= 0);
    bool incr = materialized(d, t), full = materialized_after_normalize(d, t);
    report("C green/s 3w-2z<=0  ", incr, full);
    assert(incr && "C: incremental (seeded Green) must derive 3w-2z<=0");
    assert(full && "C: full reduce must derive 3w-2z<=0");
  }

  { // Case D — Green, interior sub-case (q!=a and d!=b): Orange derives
    // w-3y<=0 into W_y, the second Resultant aligns 3 vs 2:
    //   w - 3x <= 0  ∧  [new] 2x - 2y <= 0  ∧  2y - z <= 0
    //   ==>  2w - 3z <= 0
    test_domain_t d;
    d += (w - 3 * x <= 0);
    d += (2 * y - z <= 0);
    d += (2 * x - 2 * y <= 0); // the trigger
    z_lin_cst_t t(2 * w - 3 * z <= 0);
    bool incr = materialized(d, t), full = materialized_after_normalize(d, t);
    report("D green   2w-3z<=0  ", incr, full);
    assert(incr && "D: incremental Green must derive 2w-3z<=0");
    assert(full && "D: full reduce must derive 2w-3z<=0");
  }

  { // Case E — v0 pattern on insertion, with floor rounding:
    //   y <= 5  ∧  [new] 2x - 3y <= 0   ==>  2x <= 15  ==>  x <= 7
    test_domain_t d;
    d += (y <= z_number(5));
    d += (2 * x - 3 * y <= 0);
    auto ub_incr = d[x].ub().number();
    test_domain_t f(d);
    f.normalize();
    auto ub_full = f[x].ub().number();
    bool incr = ub_incr && *ub_incr == 7, full = ub_full && *ub_full == 7;
    report("E v0-ub   x<=7      ", incr, full);
    assert(incr && "E: incremental v0 must derive x<=7 (floor 15/2)");
    assert(full && "E: full reduce must keep x<=7");
  }

  { // Case F — DOCUMENTED GAP: bound insertions do not trigger the
    // incremental pass (deliberate: a per-bound alignment sweep would cost
    // O(|T|²·N) on the most common constraint kind).  normalize() derives it.
    //   2y - x <= 0  ∧  [new] x <= 7   ==>  2y <= 7  ==>  y <= 3
    test_domain_t d;
    d += (2 * y - x <= 0);
    d += (x <= z_number(7));
    auto ub_incr = d[y].ub().number();
    test_domain_t f(d);
    f.normalize();
    auto ub_full = f[y].ub().number();
    bool incr = ub_incr && *ub_incr <= 3, full = ub_full && *ub_full <= 3;
    report("F bound   y<=3      ", incr, full);
    assert(!incr && "F: expected gap — revisit if this starts passing");
    assert(full && "F: full reduce must derive y<=3");
  }

  { // Case G — insertion-order robustness for Case C's fact.
    test_domain_t d;
    d += (2 * x - 2 * y <= 0);
    d += (3 * y - z <= 0);
    d += (w - 2 * x <= 0); // trigger arrives last, from the other side
    z_lin_cst_t t(3 * w - 2 * z <= 0);
    bool incr = materialized(d, t), full = materialized_after_normalize(d, t);
    report("G order   3w-2z<=0  ", incr, full);
    // No assert on incr: order-dependence of the incremental pass is a
    // documented property; full reduce must be order-insensitive.
    assert(full && "G: full reduce must derive 3w-2z<=0 in any order");
  }

  crab::outs() << "unittests-tvpi-reduce: done\n";
  return 0;
}
