// Regression tests for the aliasing-audit bugs (lhs occurring in the rhs)
// plus coverage for the tvpi_dbm APIs no other test exercises:
// weak_assign, select, expand, int_conv/bitwise apply, forget(vector),
// operator&=, at(), set_to_top/set_to_bottom, make_top/make_bottom, is_top,
// to_disjunctive_linear_constraint_system, and the eager Params variant.
//
// The coefficient template deliberately contains 1: get_ghost_var(v,1) == v,
// so every ghost loop that fails to skip coefficient 1 operates on the
// program variable itself (bug A4).
//
// NOTE: the test binaries build with -DNDEBUG, so plain assert() is inert;
// TCHECK prints ok/FAIL (byte-diffable) and drives the exit code.
//
// Extending this file:
//  - add a new numbered case block at the END (idx auto-increments — the
//    banner TEXT, not the number, is the stable identifier);
//  - assert through TCHECK, and print any state whose exact facts matter:
//    the byte-diffed stdout is the real oracle, ctest only sees exit codes;
//  - a case needing its own coefficient template uses scoped_template
//    (global and sorted; restored on scope exit);
//  - after an intentional behavior change, re-baseline only by reviewing
//    the output diff hunk-by-hunk (equal-or-tighter), never blindly.
#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using test_domain_t = z_tvpi_dbm_domain_t;
using z_sdbm_t = z_sdbm_domain_t;
using eager_dom_t =
    crab::domains::tvpi_dbm_domain<z_sdbm_t,
                                   crab::domains::TVPIDBMNormalizeParams>;

static unsigned idx = 1;
static unsigned failures = 0;

#define TCHECK(COND, MSG)                                                      \
  do {                                                                         \
    bool _c = (COND);                                                          \
    crab::outs() << (_c ? "  ok:   " : "  FAIL: ") << MSG << "\n";             \
    if (!_c)                                                                   \
      failures++;                                                              \
  } while (0)

// Assume lo <= v <= hi.
template <typename Dom>
static void set_range(Dom &d, const z_var &v, int lo, int hi) {
  d += (v >= z_number(lo));
  d += (v <= z_number(hi));
}

// Temporarily replace the global coefficient template (values must be given
// sorted); the previous template is restored on scope exit.  Per-case idiom
// from tests/domains/fixedtvpi.cc, made scope-safe.
struct scoped_template {
  std::vector<unsigned> m_saved;
  explicit scoped_template(std::initializer_list<unsigned> t)
      : m_saved(crab_domain_params_man::get().coefficients()) {
    crab_domain_params_man::get().coefficients().assign(t);
  }
  ~scoped_template() { crab_domain_params_man::get().coefficients() = m_saved; }
};

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

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  auto &coeffs = crab_domain_params_man::get().coefficients();
  coeffs.insert(coeffs.end(), {1, 2, 3, 4});

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var u(vfac["u"], crab::INT_TYPE, 32);
  z_var v(vfac["v"], crab::INT_TYPE, 32);
  z_var s64(vfac["s64"], crab::INT_TYPE, 64);

  { // case 1 [bug A1]: assignment whose rhs mentions the lhs with unit
    // coefficient next to a non-unit term.  The ghost refinement used to
    // rewrite x := x + 2*y into the meet x == x + ghost(y,2), where the two
    // x's cancel and the state collapses via ghost(y,2) == 0.
    crab::outs() << "---- case " << idx++ << " (A1: x := x + 2*y) ----\n";
    test_domain_t dom;
    dom.assign(y, z_number(1));
    dom.assign(x, z_number(5));
    dom.assign(x, x + z_number(2) * y);
    crab::outs() << "after x := x + 2*y: " << dom << "\n";
    TCHECK(!dom.is_bottom(), "x := x + 2*y must not collapse to bottom");
    TCHECK(dom.entails(x == z_number(7)), "x == 7");
    TCHECK(dom.entails(y == z_number(1)), "y == 1");
    TCHECK(!dom.entails(z_number(2) * y <= z_number(0)),
           "2*y <= 0 not implied");

    test_domain_t dom2;
    dom2.assign(y, z_number(1));
    dom2.assign(x, z_number(5));
    dom2.assign(x, x - z_number(2) * y); // negative variant
    TCHECK(!dom2.is_bottom(), "x := x - 2*y must not collapse to bottom");
    TCHECK(dom2.entails(x == z_number(3)), "x == 3");

    test_domain_t dom3; // control: non-unit self-reference keeps working
    dom3.assign(x, z_number(3));
    dom3.assign(x, z_number(2) * x);
    TCHECK(!dom3.is_bottom(), "x := 2*x must not collapse to bottom");
    TCHECK(dom3.entails(x == z_number(6)), "x == 6");

    test_domain_t dom4; // control: all-unit rhs (identity rewrite path)
    dom4.assign(y, z_number(1));
    dom4.assign(x, z_number(5));
    dom4.assign(x, x + y);
    TCHECK(!dom4.is_bottom(), "x := x + y must not collapse to bottom");
    TCHECK(dom4.entails(x == z_number(6)), "x == 6");
  }

  { // case 2 [bug A3]: var-RHS apply where the lhs aliases an operand and
    // the other operand is a singleton.  eval_apply used to re-enter the
    // scalar apply AFTER the base op had run: the op was applied twice
    // (x==y), or the singleton was read from the post-state (x==z).
    crab::outs() << "---- case " << idx++ << " (A3: aliased var apply) ----\n";
    test_domain_t dom;
    dom.assign(x, z_number(5));
    dom.assign(z, z_number(2));
    dom.apply(OP_ADDITION, x, x, z); // x := x + z
    crab::outs() << "after x := x + z: " << dom << "\n";
    TCHECK(dom.entails(x == z_number(7)), "x == 7 (not 9: op applied once)");

    test_domain_t dom2;
    dom2.assign(y, z_number(3));
    dom2.assign(x, z_number(2));
    dom2.apply(OP_SUBTRACTION, x, y, x); // x := y - x
    TCHECK(dom2.entails(x == z_number(1)),
           "x == 1 (singleton of x read pre-state)");

    test_domain_t dom3;
    dom3.assign(y, z_number(2));
    dom3.assign(x, z_number(5));
    dom3.apply(OP_MULTIPLICATION, x, y, x); // x := y * x, y singleton
    TCHECK(dom3.entails(x == z_number(10)), "x == 10 (not 20)");

    test_domain_t dom4; // precision lock-in: non-aliased singleton reduction
    set_range(dom4, w, 0, 10);
    dom4.assign(z, z_number(2));
    dom4.apply(OP_MULTIPLICATION, x, w, z); // x := w * z == w * 2
    crab::outs() << "after x := w * 2: " << dom4 << "\n";
    TCHECK(dom4.entails(x == z_number(2) * w), "x == 2*w (TVPI relation)");
    TCHECK(dom4.entails(x >= z_number(0)), "x >= 0");
    TCHECK(dom4.entails(x <= z_number(20)), "x <= 20");
  }

  { // case 3 [bug A4]: coefficient 1 in the template set.  Ghost loops that
    // do not skip c == 1 operate on the program variable itself: apply's
    // rewrite loop rewrites x in place, and project/rename/expand pass
    // duplicate (v,v) pairs into the base domain (CRAB_ERROR pre-fix).
    crab::outs() << "---- case " << idx++ << " (A4: coefficient 1) ----\n";
    test_domain_t dom;
    dom.assign(x, z_number(2));
    dom.apply(OP_MULTIPLICATION, y, x, z_number(3)); // y := 3*x
    crab::outs() << "after y := 3*x: " << dom << "\n";
    TCHECK(!dom.is_bottom(), "not bottom");
    TCHECK(dom.entails(y == z_number(6)), "y == 6");

    test_domain_t dp(dom);
    dp.project({x, y});
    TCHECK(dp.entails(x == z_number(2)) && dp.entails(y == z_number(6)),
           "project keeps x == 2 and y == 6");

    test_domain_t dr(dom);
    dr.rename({x}, {i}); // pre-fix: CRAB_ERROR (duplicate pair, exits)
    TCHECK(dr.entails(i == z_number(2)), "rename x -> i keeps i == 2");

    test_domain_t de(dom);
    de.expand(x, j); // pre-fix: CRAB_ERROR (j "already exists", exits)
    TCHECK(de.entails(j == z_number(2)) && de.entails(x == z_number(2)),
           "expand x -> j: both equal 2");

    test_domain_t da;
    da.assign(u, z_number(4));
    da.assign(w, z_number(5));
    da.apply(OP_ADDITION, k, u, w); // var-RHS loop must skip c == 1 too
    TCHECK(da.entails(k == z_number(9)), "k == 9");
  }

  { // case 4 [A2-class]: reuse after set_to_top; the incremental saturation
    // hook must still fire on the reset value (m_reduce_active not stuck).
    crab::outs() << "---- case " << idx++ << " (reuse after set_to_top) ----\n";
    test_domain_t dom;
    dom.assign(x, z_number(1));
    dom.assign(y, z_number(2));
    dom.normalize();
    dom.set_to_top();
    TCHECK(dom.is_top(), "is_top after set_to_top");
    dom += (z_number(3) * y - z <= z_number(0));
    dom += (z_number(2) * x - z_number(2) * y <= z_number(0));
    TCHECK(materialized(dom, z_number(3) * x - z <= z_number(0)),
           "incremental saturation derives 3x - z <= 0");
    TCHECK(!dom.entails(x == z_number(1)), "old facts do not resurrect");
  }

  { // case 5 [A5]: lattice constants and m_vars lifecycle.
    crab::outs() << "---- case " << idx++ << " (lattice constants) ----\n";
    test_domain_t dom;
    set_range(dom, x, 1, 3);
    TCHECK(!dom.is_top() && !dom.is_bottom(), "x in [1,3]: neither extreme");
    TCHECK(dom.make_top().is_top(), "make_top().is_top()");
    TCHECK(dom.make_bottom().is_bottom(), "make_bottom().is_bottom()");
    TCHECK(dom.make_bottom().at(x).is_bottom(), "at(x) on bottom is bottom");
    const test_domain_t &cdom = dom;
    auto itv = cdom.at(x);
    TCHECK(itv.lb().number() && *itv.lb().number() == z_number(1), "lb == 1");
    TCHECK(itv.ub().number() && *itv.ub().number() == z_number(3), "ub == 3");
    dom.set_to_bottom();
    TCHECK(dom.is_bottom(), "set_to_bottom");
    dom.set_to_top();
    TCHECK(dom.is_top(), "set_to_top");
    dom -= x; // forget on top: no crash, stays top
    TCHECK(dom.is_top(), "forget on top stays top");
    dom += (x == z_number(4)); // reuse after reset
    TCHECK(dom.entails(x == z_number(4)), "reusable after reset: x == 4");
  }

  { // case 6: weak_assign — may keep the old value, so only the hull.
    crab::outs() << "---- case " << idx++ << " (weak_assign) ----\n";
    test_domain_t dom;
    dom.assign(x, z_number(1));
    dom.assign(y, z_number(5));
    dom.weak_assign(y, z_number(2) * x);
    crab::outs() << "after weak y := 2*x: " << dom << "\n";
    TCHECK(!dom.entails(y == z_number(5)), "y == 5 no longer certain");
    TCHECK(!dom.entails(y == z_number(2)), "y == 2 not certain either");
    TCHECK(dom.entails(y >= z_number(2)), "y >= 2");
    TCHECK(dom.entails(y <= z_number(5)), "y <= 5");
    TCHECK(dom.entails(x == z_number(1)), "x untouched");
  }

  { // case 7: select (DEFAULT_SELECT).
    crab::outs() << "---- case " << idx++ << " (select) ----\n";
    test_domain_t dom;
    set_range(dom, x, 0, 10);
    dom.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
    crab::outs() << "select both-feasible: " << dom << "\n";
    TCHECK(dom.entails(y >= z_number(0)), "y >= 0");
    TCHECK(dom.entails(y <= z_number(30)), "y <= 30");

    test_domain_t dom2;
    dom2.assign(x, z_number(2)); // condition entailed
    dom2.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
    TCHECK(dom2.entails(y == z_number(4)), "cond true: y == 4");

    test_domain_t dom3;
    dom3.assign(x, z_number(7)); // condition refuted
    dom3.select(y, x <= z_number(5), z_number(2) * x, z_number(3) * x);
    TCHECK(dom3.entails(y == z_number(21)), "cond false: y == 21");
  }

  { // case 8: int_conv — SEXT is exact; TRUNC drops the dst's ghosts.
    crab::outs() << "---- case " << idx++ << " (int_conv) ----\n";
    test_domain_t dom;
    dom.assign(y, z_number(7));
    dom.apply(OP_SEXT, s64, y);
    TCHECK(dom.entails(s64 == z_number(7)), "sext: s64 == 7");
    dom.apply(OP_MULTIPLICATION, w, y, z_number(2)); // w == 2*y
    TCHECK(dom.entails(w == z_number(2) * y), "w == 2*y before trunc");
    dom.apply(OP_TRUNC, w, s64); // redefines w; its ghosts must go
    crab::outs() << "after trunc: " << dom << "\n";
    TCHECK(!dom.is_bottom(), "not bottom after trunc");
    TCHECK(!dom.entails(w == z_number(2) * y), "stale w == 2*y dropped");
  }

  { // case 9: bitwise, both overloads — conservative, but never stale.
    crab::outs() << "---- case " << idx++ << " (bitwise) ----\n";
    test_domain_t dom;
    dom.assign(y, z_number(12));
    dom.assign(z, z_number(10));
    dom.apply(OP_MULTIPLICATION, v, y, z_number(2)); // v == 2*y
    TCHECK(dom.entails(v == z_number(2) * y), "v == 2*y established");
    dom.apply(OP_AND, x, y, z); // var overload
    TCHECK(!dom.is_bottom(), "AND var: not bottom");
    TCHECK(dom.entails(x >= z_number(0)), "AND var: x >= 0");
    TCHECK(dom.entails(x <= z_number(12)), "AND var: x <= 12");
    dom.apply(OP_AND, w, y, z_number(4)); // scalar overload
    TCHECK(!dom.is_bottom(), "AND scalar: not bottom");
    TCHECK(dom.entails(w >= z_number(0)), "AND scalar: w >= 0");
    TCHECK(dom.entails(w <= z_number(12)), "AND scalar: w <= 12");
    dom.apply(OP_AND, v, y, z_number(3)); // redefine v
    TCHECK(!dom.entails(v == z_number(2) * y), "stale v == 2*y dropped");
  }

  { // case 10: forget(vector) and operator&=.
    crab::outs() << "---- case " << idx++ << " (forget/meet-assign) ----\n";
    test_domain_t dom;
    dom.assign(x, z_number(1));
    dom.assign(y, z_number(2));
    dom.assign(z, z_number(3));
    dom.forget({x, y});
    TCHECK(dom.entails(z == z_number(3)), "z survives forget({x,y})");
    TCHECK(dom.at(x).is_top(), "x forgotten");
    TCHECK(dom.at(y).is_top(), "y forgotten");
    dom.forget({}); // no-op
    TCHECK(dom.entails(z == z_number(3)), "empty forget is a no-op");

    test_domain_t d1, d2;
    set_range(d1, x, 0, 10);
    set_range(d2, x, 5, 20);
    d1 &= d2;
    TCHECK(d1.entails(x >= z_number(5)) && d1.entails(x <= z_number(10)),
           "meet-assign: x in [5,10]");
    test_domain_t d3, d4;
    d3 += (x == z_number(1));
    d4 += (x == z_number(2));
    d3 &= d4;
    TCHECK(d3.is_bottom(), "meet-assign of disjoint values is bottom");
    test_domain_t d5;
    d5 += (x == z_number(1));
    d5 &= d5.make_top();
    TCHECK(d5.entails(x == z_number(1)), "meet with top is identity");
    d5 &= d5.make_bottom();
    TCHECK(d5.is_bottom(), "meet with bottom is bottom");
  }

  { // case 11: to_disjunctive_linear_constraint_system.
    crab::outs() << "---- case " << idx++ << " (to_disjunctive) ----\n";
    test_domain_t dom;
    set_range(dom, x, 1, 3);
    dom.apply(OP_MULTIPLICATION, y, x, z_number(2));
    auto dcsts = dom.to_disjunctive_linear_constraint_system();
    crab::outs() << "disjunctive form: " << dcsts << "\n";
    TCHECK(!dcsts.is_false() && !dcsts.is_true(),
           "non-trivial state yields non-trivial system");
    TCHECK(
        dom.make_bottom().to_disjunctive_linear_constraint_system().is_false(),
        "bottom yields false");
    TCHECK(dom.make_top().to_disjunctive_linear_constraint_system().is_true(),
           "top yields true");
  }

  { // case 12: the eager variant over the same aliased paths.
    crab::outs() << "---- case " << idx++ << " (eager variant) ----\n";
    eager_dom_t dom;
    dom.assign(y, z_number(1));
    dom.assign(x, z_number(5));
    dom.assign(x, x + z_number(2) * y);
    TCHECK(!dom.is_bottom(), "eager: x := x + 2*y not bottom");
    TCHECK(dom.entails(x == z_number(7)), "eager: x == 7");

    eager_dom_t dom2;
    dom2.assign(x, z_number(5));
    dom2.assign(z, z_number(2));
    dom2.apply(OP_ADDITION, x, x, z);
    TCHECK(dom2.entails(x == z_number(7)), "eager: x := x + z gives 7");
  }

  { // case 13: the fixed paths under GENERAL intervals (no singleton, so
    // try_scalar_apply does not fire and the plain var-RHS path runs) and
    // relational inputs.  Aliasing must still be handled — here by the base
    // domain's own pre-state reads.
    crab::outs() << "---- case " << idx++ << " (non-singleton aliasing) ----\n";
    test_domain_t dom;
    set_range(dom, x, 1, 5);
    set_range(dom, z, 1, 2);
    dom.apply(OP_ADDITION, x, x, z); // aliased, both operands intervals
    crab::outs() << "after x := x + z: " << dom << "\n";
    TCHECK(dom.entails(x >= z_number(2)), "x >= 2");
    TCHECK(dom.entails(x <= z_number(7)), "x <= 7");
    TCHECK(!dom.entails(x <= z_number(6)), "x <= 6 not entailed (tight)");

    test_domain_t dom2; // relational input, relational output
    set_range(dom2, w, 0, 10);
    dom2.assign(y, w + z_number(3)); // y - w == 3
    set_range(dom2, z, 1, 2);
    dom2.apply(OP_ADDITION, x, y, z); // x := y + z, all non-singleton
    crab::outs() << "after x := y + z: " << dom2 << "\n";
    TCHECK(dom2.entails(x - w >= z_number(4)), "x - w >= 4 survives");
    TCHECK(dom2.entails(x - w <= z_number(5)), "x - w <= 5 survives");
    TCHECK(dom2.entails(x >= z_number(4)), "x >= 4");

    test_domain_t dom3; // A1 shape with an interval operand
    set_range(dom3, y, 1, 3);
    dom3.assign(x, z_number(5));
    dom3.assign(x, x + z_number(2) * y); // pre-fix: bottom, exactly as with
    TCHECK(!dom3.is_bottom(), "x := x + 2*y, y in [1,3]: not bottom");
    TCHECK(dom3.entails(x >= z_number(7)), "x >= 7");
    TCHECK(dom3.entails(x <= z_number(11)), "x <= 11");
  }

  { // case 14: weak_assign with a self-referential rhs (the weak path skips
    // the ghost refinement, so bug A1 cannot fire here; this pins that the
    // weak semantics stay sound under aliasing).
    crab::outs() << "---- case " << idx++ << " (weak self-reference) ----\n";
    test_domain_t dom;
    set_range(dom, y, 1, 3);
    dom.assign(x, z_number(5));
    dom.weak_assign(x, x + z_number(2) * y); // x may stay 5 or become [7,11]
    crab::outs() << "after weak x := x + 2*y: " << dom << "\n";
    TCHECK(!dom.is_bottom(), "not bottom");
    TCHECK(dom.entails(x >= z_number(5)), "x >= 5");
    TCHECK(dom.entails(x <= z_number(11)), "x <= 11");
    TCHECK(!dom.entails(x == z_number(5)), "x == 5 no longer certain");
  }

  { // case 15 [bug A6]: mid-transfer stale-ghost window.  The incremental
    // saturation hook in the scalar MUL special used to run BEFORE x's ghost
    // family was re-established: eliminate_x combined the still-stale edge
    // w - ghost(x,2) <= 5 (about x_old = 10) with the new x == 2*y (about
    // x_new = 2), deriving w <= 4*y + 5 = 9 against w == 25 -> bottom.
    // NOTE: x and w must be genuine intervals — with singletons split_dbm
    // subsumes the relational edge into the bounds and the hazard is masked.
    crab::outs() << "---- case " << idx++
                 << " (stale-ghost hook window) ----\n";
    test_domain_t dom;
    set_range(dom, x, 8, 10);
    set_range(dom, w, 20, 25);
    dom += (w - z_number(2) * x <= z_number(5)); // satisfiable: x=10, w=25
    dom.assign(y, z_number(1));
    dom.apply(OP_MULTIPLICATION, x, y, z_number(2)); // x := 2*y == 2
    crab::outs() << "after x := 2*y: " << dom << "\n";
    TCHECK(!dom.is_bottom(), "reachable state must not collapse to bottom");
    TCHECK(dom.entails(x == z_number(2)), "x == 2");
    TCHECK(dom.entails(w >= z_number(20)), "w >= 20 (old bound not misused)");
    TCHECK(!dom.entails(w <= z_number(9)),
           "no derivation may conflate x_old with x_new");
  }

  { // case 16: weak_assign, relational content.  weak y := 2*x from y == 10
    // joins {y == 10} with {y == 2*x}; both branches satisfy 2*x <= y, so
    // the relational lower bound must survive the weak join.
    crab::outs() << "---- case " << idx++ << " (weak_assign relational) ----\n";
    test_domain_t dom;
    set_range(dom, x, 1, 3);
    dom.assign(y, z_number(10));
    dom.weak_assign(y, z_number(2) * x);
    crab::outs() << "after weak y := 2*x: " << dom << "\n";
    TCHECK(!dom.entails(y == z_number(10)), "y == 10 no longer certain");
    TCHECK(!dom.entails(y == z_number(2) * x), "y == 2*x not certain either");
    TCHECK(dom.entails(y >= z_number(2)), "y >= 2");
    TCHECK(dom.entails(y <= z_number(10)), "y <= 10");
    TCHECK(dom.entails(z_number(2) * x <= y),
           "2*x <= y holds in both branches of the weak join");
  }

  { // case 17: ghost-family exactness after each rewritten transfer
    // function.  Expected ghost values are computed by hand in the comments;
    // the printed states lock them into the baseline.
    crab::outs() << "---- case " << idx++ << " (family exactness) ----\n";

    { // (a) assign, non-unit self-reference plus a second term.
      // x==5, y==3, x := 2*x + y: x_new = 13.
      //   c=2: ghost(x,2) := ghost(x,4) + ghost(y,2) = 20 + 6 = 26 = 2*13 ok
      //   c=3: needs ghost(x,6) — not in template -> dropped
      //   c=4: needs ghost(x,8) — not in template -> dropped
      test_domain_t d;
      d.assign(x, z_number(5));
      d.assign(y, z_number(3));
      d.assign(x, z_number(2) * x + y);
      crab::outs() << "(a) " << d << "\n";
      TCHECK(!d.is_bottom(), "a: not bottom");
      TCHECK(d.entails(x == z_number(13)), "a: x == 13");
      TCHECK(d.entails(z_number(2) * x - y == z_number(23)),
             "a: 2*x - y == 23 (through ghost(x,2) = 26)");
    }
    { // (b) scalar MUL self-reference over an interval, with a bystander
      // relation.  x in [2,5], w := x, x := x*3: x_new in [6,15];
      //   c=2: ghost(x,2) := ghost(x,2)_old * 3 = [12,30] = 2*x_new ok
      //   (same shape for c=3: [18,45], c=4: [24,60])
      test_domain_t d;
      set_range(d, x, 2, 5);
      d.assign(w, x);
      d.apply(OP_MULTIPLICATION, x, x, z_number(3));
      crab::outs() << "(b) " << d << "\n";
      TCHECK(!d.is_bottom(), "b: not bottom");
      TCHECK(d.entails(x >= z_number(6)), "b: x >= 6");
      TCHECK(d.entails(x <= z_number(15)), "b: x <= 15");
      TCHECK(!d.entails(x <= z_number(14)), "b: x <= 14 not entailed (tight)");
      TCHECK(d.entails(w >= z_number(2)) && d.entails(w <= z_number(5)),
             "b: bystander w keeps [2,5]");
    }
    { // (c) var-RHS ADD, aliased, both operands intervals.
      // x in [1,3], z in [1,2], x := x + z: x_new in [2,5];
      //   c: ghost(x,c) := ghost(x,c)_old + ghost(z,c) = c*[2,5] ok
      test_domain_t d;
      set_range(d, x, 1, 3);
      set_range(d, z, 1, 2);
      d.apply(OP_ADDITION, x, x, z);
      crab::outs() << "(c) " << d << "\n";
      TCHECK(d.entails(x >= z_number(2)), "c: x >= 2");
      TCHECK(d.entails(x <= z_number(5)), "c: x <= 5");
    }
    { // (d) chained composition: identity-rewrite path, then non-unit
      // self-reference, then mixed unit/non-unit self-reference.
      // x=1, y=2: x := x+1 -> 2; x := 2*x -> 4; x := x + 2*y -> 8;
      //   final ghost(x,2) = ghost(x,2)_pre + ghost(y,4) = 8 + 8 = 16 ok
      //   ghost(x,3), ghost(x,4): dropped along the way (out-of-template
      //   reads), never stale.
      test_domain_t d;
      d.assign(x, z_number(1));
      d.assign(y, z_number(2));
      d.assign(x, x + z_number(1));
      d.assign(x, z_number(2) * x);
      d.assign(x, x + z_number(2) * y);
      crab::outs() << "(d) " << d << "\n";
      TCHECK(!d.is_bottom(), "d: not bottom");
      TCHECK(d.entails(x == z_number(8)), "d: x == 8");
      TCHECK(d.entails(z_number(2) * x - y == z_number(14)),
             "d: 2*x - y == 14 (through ghost(x,2) = 16)");
    }
  }

  { // case 18: the eager variant as an independent oracle.  Its normalize()
    // runs full saturation after every operation, meeting each family
    // member against x and the others — a ghost holding a stale or wrong
    // value tends to collapse the state to bottom (Scaling rule turns a
    // wrong ghost bound into a contradicting bound on x).  Replaying the
    // case-17 sequences eagerly must reproduce the same facts, bottom-free.
    crab::outs() << "---- case " << idx++ << " (eager oracle replay) ----\n";
    {
      eager_dom_t d;
      d.assign(x, z_number(5));
      d.assign(y, z_number(3));
      d.assign(x, z_number(2) * x + y);
      TCHECK(!d.is_bottom(), "a: not bottom");
      TCHECK(d.entails(x == z_number(13)), "a: x == 13");
    }
    {
      eager_dom_t d;
      set_range(d, x, 2, 5);
      d.assign(w, x);
      d.apply(OP_MULTIPLICATION, x, x, z_number(3));
      TCHECK(!d.is_bottom(), "b: not bottom");
      TCHECK(d.entails(x >= z_number(6)) && d.entails(x <= z_number(15)),
             "b: x in [6,15]");
      TCHECK(d.entails(w <= z_number(5)), "b: w <= 5");
    }
    {
      eager_dom_t d;
      d.assign(x, z_number(1));
      d.assign(y, z_number(2));
      d.assign(x, x + z_number(1));
      d.assign(x, z_number(2) * x);
      d.assign(x, x + z_number(2) * y);
      TCHECK(!d.is_bottom(), "d: not bottom");
      TCHECK(d.entails(x == z_number(8)), "d: x == 8");
    }
  }

  { // case 19 [bug A7]: aliased division x := x / z.  The DIV special used
    // to write ghost(x,|z|) := x AFTER the base apply (post-state), and the
    // rewrite loop then divided that again: y == 7, y := y/2 left every
    // ghost of y wrong (ghost(y,2) == 1 while 2*y == 6), and entailment
    // queries through the ghosts answered unsoundly.
    crab::outs() << "---- case " << idx++ << " (aliased division) ----\n";
    test_domain_t d;
    d.assign(y, z_number(7));
    d.apply(OP_SDIV, y, y, z_number(2));
    crab::outs() << "after y := y / 2: " << d << "\n";
    TCHECK(d.entails(y == z_number(3)), "y == 3");
    TCHECK(!d.entails(z_number(2) * y <= z_number(5)),
           "2*y <= 5 must not be entailed (2*y == 6)");
    d.normalize();
    TCHECK(!d.is_bottom(), "normalize must not collapse a reachable state");

    eager_dom_t d2; // eager oracle: normalize runs right after the apply
    d2.assign(y, z_number(7));
    d2.apply(OP_SDIV, y, y, z_number(2));
    TCHECK(!d2.is_bottom(), "eager: not bottom");
    TCHECK(d2.entails(y == z_number(3)), "eager: y == 3");
  }

  { // case 20 [DIV fix]: sound division bands.  x := y / z with a = |z|
    // satisfies y = sign(z)*a*x + r, |r| <= a-1, sign(r) = sign(y); the old
    // code asserted the false equality a*x == sign(z)*y instead.
    crab::outs() << "---- case " << idx++ << " (sound division bands) ----\n";
    { // P1/P2: y in [10,21], x := y/2 — one-sided band, both ghost pairs.
      test_domain_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_SDIV, x, y, z_number(2));
      crab::outs() << "P1 " << d << "\n";
      TCHECK(d.entails(x >= z_number(5)) && d.entails(x <= z_number(10)),
             "P1: x in [5,10]");
      TCHECK(!d.entails(z_number(2) * x >= y),
             "P1: 2*x >= y NOT entailed (y may be odd)");
      TCHECK(d.entails(z_number(2) * x <= y), "P1: 2*x <= y (band edge)");
      TCHECK(d.entails(y - z_number(2) * x <= z_number(1)),
             "P1: y - 2*x <= 1 (band edge)");
      TCHECK(d.entails(z_number(2) * y - z_number(4) * x <= z_number(3)),
             "P2: 2*y - 4*x <= 3 (band at ghost pair (2,4))");
    }
    { // P1s: odd singleton — the old equality claimed 2*x == 21.
      test_domain_t d;
      d.assign(y, z_number(21));
      d.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(d.entails(x == z_number(10)), "P1s: x == 10");
      TCHECK(!d.entails(z_number(2) * x >= z_number(21)),
             "P1s: 2*x >= 21 NOT entailed (2*x == 20)");
    }
    { // P3: a does not divide c — interval fallback must be exact.
      test_domain_t d;
      d.assign(y, z_number(7));
      d.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(d.entails(x == z_number(3)), "P3: x == 3");
      TCHECK(!d.entails(z_number(3) * x >= z_number(10)),
             "P3: 3*x >= 10 NOT entailed (3*x == 9, old ghost said 10)");
      TCHECK(d.entails(z_number(3) * x >= z_number(8)), "P3: 3*x >= 8");
    }
    { // P4: negative divisor, odd singleton — old ghost said 2*x == -21.
      test_domain_t d;
      d.assign(y, z_number(21));
      d.apply(OP_SDIV, x, y, z_number(-2));
      TCHECK(d.entails(x == z_number(-10)), "P4: x == -10");
      TCHECK(!d.entails(z_number(2) * x <= z_number(-21)),
             "P4: 2*x <= -21 NOT entailed (2*x == -20)");
      TCHECK(d.entails(z_number(2) * x <= z_number(-19)), "P4: 2*x <= -19");
    }
    { // P5: unknown sign of y — two-sided band, and no false 2*x <= y.
      test_domain_t d;
      set_range(d, y, -5, 9);
      d.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(!d.entails(z_number(2) * x <= y),
             "P5: 2*x <= y NOT entailed (y=-5 gives x=-2, 2x=-4 > -5)");
      TCHECK(d.entails(y - z_number(2) * x <= z_number(1)) &&
                 d.entails(z_number(2) * x - y <= z_number(1)),
             "P5: two-sided band |y - 2*x| <= 1");
    }
    { // P6: UDIV over entailed-nonnegative y behaves like SDIV.
      test_domain_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_UDIV, x, y, z_number(2));
      TCHECK(d.entails(z_number(2) * x <= y), "P6: udiv 2*x <= y");
      TCHECK(!d.entails(z_number(2) * x >= y),
             "P6: udiv 2*x >= y NOT entailed");
    }
    { // eager replay: the band must survive (and feed) full normalization.
      eager_dom_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(!d.is_bottom(), "eager: not bottom");
      TCHECK(d.entails(z_number(2) * x <= y), "eager: band survives");
    }
  }

  { // case 21: remaining DIV branches — negative dividends, UDIV gating,
    // and the var-RHS drop.
    crab::outs() << "---- case " << idx++ << " (DIV branch coverage) ----\n";
    { // (a) SDIV, y <= 0, z > 0: band flips one-sided the other way,
      // y - 2*x in [-1, 0].
      test_domain_t d;
      set_range(d, y, -21, -10);
      d.apply(OP_SDIV, x, y, z_number(2));
      crab::outs() << "(a) " << d << "\n";
      TCHECK(d.entails(x >= z_number(-10)) && d.entails(x <= z_number(-5)),
             "a: x in [-10,-5]");
      TCHECK(d.entails(y - z_number(2) * x <= z_number(0)), "a: y <= 2*x");
      TCHECK(d.entails(z_number(2) * x - y <= z_number(1)), "a: 2*x - y <= 1");
      TCHECK(!d.entails(z_number(2) * x <= y),
             "a: 2*x <= y NOT entailed (y=-21 gives x=-10, 2x=-20 > -21)");

      test_domain_t ds; // odd singleton
      ds.assign(y, z_number(-21));
      ds.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(ds.entails(x == z_number(-10)), "a: singleton x == -10");
      TCHECK(!ds.entails(z_number(2) * x <= z_number(-21)),
             "a: 2*x <= -21 NOT entailed (2*x == -20)");
    }
    { // (b) SDIV, y <= 0, z < 0: sum band, interval content only.
      test_domain_t d;
      d.assign(y, z_number(-21));
      d.apply(OP_SDIV, x, y, z_number(-2));
      TCHECK(d.entails(x == z_number(10)), "b: x == 10");
      TCHECK(!d.entails(z_number(2) * x >= z_number(21)),
             "b: 2*x >= 21 NOT entailed (2*x == 20)");
    }
    { // (c) UDIV with a possibly-negative dividend: NO band may be added
      // (unsigned reinterpretation), only the interval fallback.
      test_domain_t d;
      set_range(d, y, -5, 9);
      d.apply(OP_UDIV, x, y, z_number(2));
      TCHECK(!d.is_bottom(), "c: not bottom");
      TCHECK(!d.entails(y - z_number(2) * x <= z_number(1)),
             "c: no band for sign-unknown udiv");
    }
    { // (d) UDIV with a negative divisor: no band either.
      test_domain_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_UDIV, x, y, z_number(-2));
      TCHECK(!d.is_bottom(), "d: not bottom");
    }
    { // (e) var-RHS DIV, non-singleton divisor: x's stale facts must drop.
      test_domain_t d;
      d.assign(w, z_number(10));
      d.apply(OP_MULTIPLICATION, v, w, z_number(2)); // v == 2*w
      TCHECK(d.entails(v == z_number(2) * w), "e: v == 2*w established");
      set_range(d, z, 2, 3);
      d.apply(OP_SDIV, v, w, z); // redefines v; divisor not singleton
      TCHECK(!d.is_bottom(), "e: not bottom");
      TCHECK(!d.entails(v == z_number(2) * w), "e: stale v == 2*w dropped");
      TCHECK(d.entails(v >= z_number(3)) && d.entails(v <= z_number(5)),
             "e: v in [3,5] (10/3 .. 10/2)");
    }
    { // (f) eager replay of (a).
      eager_dom_t d;
      set_range(d, y, -21, -10);
      d.apply(OP_SDIV, x, y, z_number(2));
      TCHECK(!d.is_bottom(), "f: eager not bottom");
      TCHECK(d.entails(y - z_number(2) * x <= z_number(0)),
             "f: eager band survives");
    }
  }

  { // case 22: template-gap branches, which {1,2,3,4} cannot reach —
    // switch to T = {1,6} (fixedtvpi per-case idiom, restored below).
    crab::outs() << "---- case " << idx++ << " (DIV template gaps) ----\n";
    scoped_template tmpl({1, 6});
    { // (a) a | c but q = c/a outside the template: c=6, a=2 -> q=3 not in
      // {1,6}: interval fallback only, no band at any pair.
      test_domain_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_SDIV, x, y, z_number(2));
      crab::outs() << "(a) " << d << "\n";
      TCHECK(d.entails(x >= z_number(5)) && d.entails(x <= z_number(10)),
             "a: x in [5,10]");
      TCHECK(d.entails(z_number(6) * x >= z_number(29)),
             "a: 6*x >= 29 (fallback ghost(x,6) = [30,60])");
      TCHECK(!d.entails(y - z_number(2) * x <= z_number(1)),
             "a: no band (neither coefficient 2 nor q=3 tracked)");
    }
    { // (b) a itself outside the template (z=5): every c hits a∤c.
      test_domain_t d;
      set_range(d, y, 10, 21);
      d.apply(OP_SDIV, x, y, z_number(5));
      TCHECK(d.entails(x >= z_number(2)) && d.entails(x <= z_number(4)),
             "b: x in [2,4]");
      TCHECK(d.entails(z_number(6) * x <= z_number(25)),
             "b: 6*x <= 25 (fallback ghost(x,6) = [12,24])");
    }
  } // scoped_template restores {1,2,3,4}

  crab::outs() << "unittests-tvpi4: " << (failures ? "FAILED" : "passed")
               << " (" << failures << " failing checks)\n";
  return failures ? 1 : 0;
}
