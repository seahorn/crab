#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/domains/symbolic_variable_eq_domain.hpp>

#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

// Instantiate with checks enabled so the unit test always validates the
// domain's internal representation in addition to its observable results.
using test_domain_t = crab::domains::symbolic_variable_equality_domain<
    z_number, varname_t, crab::domains::SVEQCheckedParams>;
using value_domain_t = symbolic_variable_equality_domain_impl::class_id_t;

// The *meaning* of an EqDom state is a partition of variables into equivalence
// classes (variables sharing a class are known equal). We describe expected
// results as partitions; variables not listed in any class are singletons
// (equal to nothing).
//
// This test is self-contained: it never asks the domain "are x and y equal?"
// (that query, equals(), is itself part of the code under test). Instead it
// builds ground-truth reference domains from partitions and checks the computed
// results using only the public lattice interface (<=, |, &). The whole
// partition check reduces to entailment of single equalities:
//   "dom entails a == b"  <=>  dom <= make_dom({{a, b}})
// and <= is in turn pinned down by the asymmetric leq expectations in each case,
// so a degenerate <= (always-true / always-false) cannot pass silently.
using partition_t = std::vector<std::vector<z_var>>;

namespace {
unsigned g_checks = 0;
unsigned g_failures = 0;

// Trusted (pure test-side) oracle: are a and b in the same class of `parts`?
bool expected_equal(const partition_t &parts, const z_var &a, const z_var &b) {
  if (a == b) {
    return true;
  }
  for (const auto &cls : parts) {
    bool has_a = false, has_b = false;
    for (const auto &v : cls) {
      has_a = has_a || (v == a);
      has_b = has_b || (v == b);
    }
    if (has_a && has_b) {
      return true;
    }
  }
  return false;
}

// Build a reference domain denoting exactly `parts`. A fresh class id per
// class avoids set()'s merge-on-equal-id behaviour.
test_domain_t make_dom(const partition_t &parts) {
  test_domain_t d;
  for (const auto &cls : parts) {
    if (cls.size() < 2) {
      continue; // singletons carry no equality
    }
    d.set(cls[0], test_domain_t::fresh_class_id());
    for (unsigned i = 1; i < cls.size(); ++i) {
      d.add(cls[0], cls[i]);
    }
  }
  return d;
}

void print_partition(crab::crab_os &o, const partition_t &parts) {
  o << "{";
  bool first = true;
  for (const auto &cls : parts) {
    if (cls.size() < 2) {
      continue;
    }
    if (!first) {
      o << ", ";
    }
    first = false;
    o << "[";
    for (unsigned i = 0; i < cls.size(); ++i) {
      if (i) {
        o << ",";
      }
      o << cls[i];
    }
    o << "]";
  }
  o << "}";
}

// Verify that `dom` entails exactly the equalities in `expected` over
// `universe`, using only <= against minimal reference domains. On any
// disagreement, report the computed state, the expected partition, and each
// offending pair with the reason.
void check_partition(const std::string &what, const test_domain_t &dom,
                     const std::vector<z_var> &universe,
                     const partition_t &expected) {
  ++g_checks;
  struct mismatch_t {
    z_var a, b;
    bool expected_eq, actual_eq;
  };
  std::vector<mismatch_t> mismatches;
  for (unsigned i = 0; i < universe.size(); ++i) {
    for (unsigned j = i + 1; j < universe.size(); ++j) {
      const z_var &a = universe[i];
      const z_var &b = universe[j];
      bool exp = expected_equal(expected, a, b);
      bool act = dom <= make_dom({{a, b}}); // does dom entail a == b ?
      if (exp != act) {
        mismatches.push_back({a, b, exp, act});
      }
    }
  }
  if (mismatches.empty()) {
    crab::outs() << "[PASS] " << what << "\n";
    return;
  }
  ++g_failures;
  crab::outs() << "[FAIL] " << what << "\n";
  crab::outs() << "       computed : " << dom << "\n";
  crab::outs() << "       expected : ";
  print_partition(crab::outs(), expected);
  crab::outs() << "\n       reason   :\n";
  for (const auto &m : mismatches) {
    crab::outs() << "         - expected " << m.a
                 << (m.expected_eq ? " == " : " != ") << m.b
                 << " but domain " << (m.actual_eq ? "entails " : "does not entail ")
                 << m.a << " == " << m.b << "\n";
  }
}

void check_bool(const std::string &what, bool actual, bool expected) {
  ++g_checks;
  if (actual == expected) {
    crab::outs() << "[PASS] " << what << "\n";
    return;
  }
  ++g_failures;
  crab::outs() << "[FAIL] " << what << "\n";
  crab::outs() << "       computed : " << (actual ? "true" : "false") << "\n";
  crab::outs() << "       expected : " << (expected ? "true" : "false") << "\n";
}

// Exercise and verify <=, join and meet for a pair of domains:
//  - the inputs are built as intended (catches accidental set()-merges),
//  - leq holds in the expected directions,
//  - join/meet produce the expected partitions,
//  - the universal lattice laws hold (join is an upper bound, meet a lower one).
void check_latticeops(const std::string &name, const test_domain_t &dom1,
                      const test_domain_t &dom2,
                      const std::vector<z_var> &universe, const partition_t &p1,
                      const partition_t &p2, bool exp_leq_12, bool exp_leq_21,
                      const partition_t &exp_join, const partition_t &exp_meet) {
  crab::outs() << "==== " << name << " ====\n";
  check_partition(name + ": dom1 built as intended", dom1, universe, p1);
  check_partition(name + ": dom2 built as intended", dom2, universe, p2);

  check_bool(name + ": dom1 <= dom2", dom1 <= dom2, exp_leq_12);
  check_bool(name + ": dom2 <= dom1", dom2 <= dom1, exp_leq_21);

  test_domain_t join = dom1 | dom2;
  check_partition(name + ": join (dom1 | dom2)", join, universe, exp_join);
  check_bool(name + ": dom1 <= join (upper bound)", dom1 <= join, true);
  check_bool(name + ": dom2 <= join (upper bound)", dom2 <= join, true);

  test_domain_t meet = dom1 & dom2;
  check_partition(name + ": meet (dom1 & dom2)", meet, universe, exp_meet);
  check_bool(name + ": meet <= dom1 (lower bound)", meet <= dom1, true);
  check_bool(name + ": meet <= dom2 (lower bound)", meet <= dom2, true);
}
} // namespace

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_var v1(vfac["v1"], crab::INT_TYPE, 32);
  z_var v2(vfac["v2"], crab::INT_TYPE, 32);
  z_var v3(vfac["v3"], crab::INT_TYPE, 32);
  z_var v4(vfac["v4"], crab::INT_TYPE, 32);
  z_var v5(vfac["v5"], crab::INT_TYPE, 32);
  z_var v6(vfac["v6"], crab::INT_TYPE, 32);
  z_var v7(vfac["v7"], crab::INT_TYPE, 32);
  z_var v8(vfac["v8"], crab::INT_TYPE, 32);
  z_var v9(vfac["v9"], crab::INT_TYPE, 32);
  z_var v10(vfac["v10"], crab::INT_TYPE, 32);
  z_var v11(vfac["v11"], crab::INT_TYPE, 32);
  z_var v12(vfac["v12"], crab::INT_TYPE, 32);

  // NOTE on set(): set(x, id) merges x into any existing class that already
  // holds the same id. So reusing an idom value across two set() calls in the
  // same domain merges those classes -- this is intended (it is how
  // object_domain establishes equalities), and a few cases below rely on it.
  // Ids must come from fresh_class_id(); invented ids are rejected by set().
  value_domain_t idom1 = test_domain_t::fresh_class_id();
  value_domain_t idom2 = test_domain_t::fresh_class_id();
  value_domain_t idom3 = test_domain_t::fresh_class_id();

  { // disjoint inputs: no shared equalities
    test_domain_t dom1, dom2;
    dom1.set(v1, idom1);
    dom1.add(v1, v2); // dom1 : {v1,v2}
    dom2.set(v4, idom2);
    dom2.add(v4, v5); // dom2 : {v4,v5}
    check_latticeops("case0 (disjoint)", dom1, dom2, {v1, v2, v4, v5},
                     /*p1*/ {{v1, v2}}, /*p2*/ {{v4, v5}},
                     /*leq12*/ false, /*leq21*/ false,
                     /*join*/ {}, /*meet*/ {{v1, v2}, {v4, v5}});
  }

  { // inputs share variable v2 but in different classes
    test_domain_t dom1, dom2;
    dom1.set(v1, idom1);
    dom1.add(v1, v2); // dom1 : {v1,v2}
    dom2.set(v4, idom2);
    dom2.add(v4, v2); // dom2 : {v2,v4}
    check_latticeops("case1 (share v2)", dom1, dom2, {v1, v2, v4},
                     {{v1, v2}}, {{v2, v4}}, false, false,
                     /*join*/ {}, /*meet*/ {{v1, v2, v4}});
  }

  { // set()-merge: reusing idom1 for v2 and v5 merges them into one class
    test_domain_t dom1, dom2;
    dom1.set(v2, idom1);
    dom1.add(v2, v3);
    dom1.set(v5, idom1); // merges v5 into {v2,v3} (same symbolic value)
    dom1.add(v5, v1);    // dom1 : {v1,v2,v3,v5}
    dom2.set(v1, idom2);
    dom2.add(v1, v2);    // dom2 : {v1,v2}
    check_latticeops("case2 (set-merge in dom1)", dom1, dom2, {v1, v2, v3, v5},
                     {{v1, v2, v3, v5}}, {{v1, v2}},
                     /*leq12*/ true, /*leq21*/ false,
                     /*join*/ {{v1, v2}}, /*meet*/ {{v1, v2, v3, v5}});
  }

  { // refinement: dom1 has all of dom2's equalities and more
    test_domain_t dom1, dom2;
    dom1.set(v1, idom1);
    dom1.add(v1, v2);
    dom1.set(v3, idom3);
    dom1.add(v3, v4); // dom1 : {v1,v2},{v3,v4}
    dom2.set(v1, idom3);
    dom2.add(v1, v2); // dom2 : {v1,v2}
    check_latticeops("case3 (refinement)", dom1, dom2, {v1, v2, v3, v4},
                     {{v1, v2}, {v3, v4}}, {{v1, v2}},
                     /*leq12*/ true, /*leq21*/ false,
                     /*join*/ {{v1, v2}}, /*meet*/ {{v1, v2}, {v3, v4}});
  }

  { // single class on each side + forget/rename/project
    test_domain_t dom1, dom2;
    dom1.set(v1, idom1);
    dom1.add(v1, v2);
    dom1.add(v2, v3);
    dom1.add(v3, v4); // dom1 : {v1,v2,v3,v4}
    dom2.set(v2, idom2);
    dom2.add(v2, v1);
    dom2.add(v1, v3); // dom2 : {v1,v2,v3}
    check_latticeops("case4 (chain)", dom1, dom2, {v1, v2, v3, v4},
                     {{v1, v2, v3, v4}}, {{v1, v2, v3}},
                     /*leq12*/ true, /*leq21*/ false,
                     /*join*/ {{v1, v2, v3}}, /*meet*/ {{v1, v2, v3, v4}});

    test_domain_t forgotten(dom1);
    forgotten -= v2; // drop v2 from {v1,v2,v3,v4}
    check_partition("case4: forget v2", forgotten, {v1, v2, v3, v4},
                    {{v1, v3, v4}});

    test_domain_t renamed(dom1);
    renamed.rename({v1, v2, v3, v4}, {v5, v6, v7, v8});
    check_partition("case4: rename {v1..v4} -> {v5..v8}", renamed,
                    {v1, v2, v3, v4, v5, v6, v7, v8}, {{v5, v6, v7, v8}});

    test_domain_t projected(dom1);
    projected.project({v1, v3});
    check_partition("case4: project on {v1,v3}", projected, {v1, v2, v3, v4},
                    {{v1, v3}});
  }

  { // two classes on each side, fully crossing + forget/rename/project
    test_domain_t dom1, dom2;
    dom1.set(v1, idom2);
    dom1.add(v1, v2);
    dom1.set(v3, idom3);
    dom1.add(v3, v4); // dom1 : {v1,v2},{v3,v4}
    dom2.set(v2, idom3);
    dom2.add(v2, v3);
    dom2.set(v1, idom1);
    dom2.add(v1, v4); // dom2 : {v2,v3},{v1,v4}
    check_latticeops("case5 (crossing)", dom1, dom2, {v1, v2, v3, v4},
                     {{v1, v2}, {v3, v4}}, {{v2, v3}, {v1, v4}}, false, false,
                     /*join*/ {}, /*meet*/ {{v1, v2, v3, v4}});

    test_domain_t forgotten(dom1);
    forgotten -= v2; // {v1,v2} loses v2 -> v1 becomes a singleton
    check_partition("case5: forget v2", forgotten, {v1, v2, v3, v4},
                    {{v3, v4}});

    test_domain_t renamed(dom1);
    renamed.rename({v1, v2, v3, v4}, {v5, v6, v7, v8});
    check_partition("case5: rename {v1..v4} -> {v5..v8}", renamed,
                    {v1, v2, v3, v4, v5, v6, v7, v8}, {{v5, v6}, {v7, v8}});

    test_domain_t projected(dom1);
    projected.project({v1, v3}); // v1,v3 are in different classes -> no equality
    check_partition("case5: project on {v1,v3}", projected, {v1, v2, v3, v4},
                    {});
  }

  { // three classes on each side, heavily crossing
    test_domain_t dom1, dom2;
    dom1.set(v6, idom2);
    dom1.add(v6, v1);
    dom1.add(v1, v2);
    dom1.add(v6, v8);
    dom1.add(v8, v11);
    dom1.set(v3, idom3);
    dom1.add(v3, v7);
    dom1.add(v7, v12);
    dom1.set(v4, idom1);
    dom1.add(v4, v9);
    dom1.add(v4, v10);
    dom1.add(v4, v5);
    // dom1 : {v1,v2,v6,v8,v11},{v3,v7,v12},{v4,v5,v9,v10}
    dom2.set(v5, idom3);
    dom2.add(v5, v12);
    dom2.add(v5, v4);
    dom2.add(v4, v10);
    dom2.add(v12, v2);
    dom2.set(v6, idom1);
    dom2.add(v6, v8);
    dom2.add(v6, v7);
    dom2.add(v7, v9);
    dom2.set(v3, idom2);
    dom2.add(v3, v11);
    dom2.add(v11, v1);
    // dom2 : {v2,v4,v5,v10,v12},{v6,v7,v8,v9},{v1,v3,v11}
    std::vector<z_var> all = {v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12};
    check_latticeops(
        "case6 (heavy crossing)", dom1, dom2, all,
        {{v1, v2, v6, v8, v11}, {v3, v7, v12}, {v4, v5, v9, v10}},
        {{v2, v4, v5, v10, v12}, {v6, v7, v8, v9}, {v1, v3, v11}}, false, false,
        /*join*/ {{v6, v8}, {v1, v11}, {v4, v5, v10}},
        /*meet*/
        {{v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12}});

    test_domain_t forgotten(dom1);
    forgotten -= v2;
    check_partition("case6: forget v2", forgotten, all,
                    {{v1, v6, v8, v11}, {v3, v7, v12}, {v4, v5, v9, v10}});

    test_domain_t proj1(dom1);
    proj1.project({v1, v3}); // different classes in dom1
    check_partition("case6: project dom1 on {v1,v3}", proj1, all, {});

    test_domain_t proj2(dom2);
    proj2.project({v1, v3}); // same class {v1,v3,v11} in dom2
    check_partition("case6: project dom2 on {v1,v3}", proj2, all, {{v1, v3}});
  }

  { // set()-merge: reusing idom3 for v10 and v11 collapses dom2 into one class
    test_domain_t dom1, dom2;
    dom1.set(v12, idom2);
    dom1.add(v12, v1);
    dom1.add(v12, v2);
    dom1.add(v12, v3);
    dom1.add(v12, v4);
    dom1.set(v11, idom3);
    dom1.add(v11, v10); // dom1 : {v1,v2,v3,v4,v12},{v10,v11}
    dom2.set(v10, idom3);
    dom2.add(v10, v2);
    dom2.add(v10, v3);
    dom2.set(v11, idom3); // merges v11 into {v2,v3,v10} (same symbolic value)
    dom2.add(v11, v4);
    dom2.add(v4, v12); // dom2 : {v2,v3,v4,v10,v11,v12}
    std::vector<z_var> univ = {v1, v2, v3, v4, v10, v11, v12};
    check_latticeops("case7 (set-merge in dom2)", dom1, dom2, univ,
                     {{v1, v2, v3, v4, v12}, {v10, v11}},
                     {{v2, v3, v4, v10, v11, v12}}, false, false,
                     /*join*/ {{v2, v3, v4, v12}, {v10, v11}},
                     /*meet*/ {{v1, v2, v3, v4, v10, v11, v12}});

    test_domain_t forgotten(dom1);
    forgotten -= v10; // {v10,v11} loses v10 -> v11 becomes a singleton
    check_partition("case7: forget v10", forgotten, univ,
                    {{v1, v2, v3, v4, v12}});

    test_domain_t proj1(dom1);
    proj1.project({v1, v3}); // same class in dom1
    check_partition("case7: project dom1 on {v1,v3}", proj1, univ, {{v1, v3}});

    test_domain_t proj2(dom2);
    proj2.project({v1, v3}); // v1 not in dom2 -> no equality
    check_partition("case7: project dom2 on {v1,v3}", proj2, univ, {});
  }

  { // object-domain-like: dom2 refines into dom1's first class only
    test_domain_t dom1, dom2;
    dom1.set(v1, idom1);
    dom1.add(v1, v2);
    dom1.set(v3, idom3);
    dom1.add(v3, v4); // dom1 : {v1,v2},{v3,v4}
    dom2.set(v1, idom3);
    dom2.add(v1, v2); // dom2 : {v1,v2}
    check_latticeops("case8 (object-like)", dom1, dom2, {v1, v2, v3, v4},
                     {{v1, v2}, {v3, v4}}, {{v1, v2}},
                     /*leq12*/ true, /*leq21*/ false,
                     /*join*/ {{v1, v2}}, /*meet*/ {{v1, v2}, {v3, v4}});
  }

  { // singletons only: both states carry no equalities (equiv. to top)
    test_domain_t dom1, dom2;
    dom1.set(v4, idom3); // dom1 : {v4} (singleton)
    dom2.set(v3, idom2); // dom2 : {v3} (singleton)
    check_latticeops("case9 (singletons)", dom1, dom2, {v3, v4},
                     /*p1*/ {}, /*p2*/ {},
                     /*leq12*/ true, /*leq21*/ true,
                     /*join*/ {}, /*meet*/ {});
  }

  { // adding an equality to top must create a class (top is not absorbing)
    crab::outs() << "==== case10 (add on top) ====\n";
    test_domain_t dom; // default-constructed == top
    check_bool("case10: fresh domain is top", dom.is_top(), true);
    check_partition("case10: fresh domain has no equalities", dom, {v1, v2},
                    {});
    dom.add(v1, v2);
    check_bool("case10: no longer top after add", dom.is_top(), false);
    check_partition("case10: top.add(v1,v2)", dom, {v1, v2}, {{v1, v2}});
  }

  { // normalize drops singleton classes (no equality), collapsing to top
    crab::outs() << "==== case11 (normalize) ====\n";
    test_domain_t dom;
    dom.set(v1, idom1); // a lone class {v1} carries no equality
    check_bool("case11: singleton present before normalize", dom.is_top(),
               false);
    dom.normalize();
    check_bool("case11: singleton dropped by normalize (now top)", dom.is_top(),
               true);
    check_partition("case11: normalized state has no equalities", dom, {v1, v2},
                    {});
  }

  { // to_linear_constraint_system concretizes class equalities
    crab::outs() << "==== case12 (to_linear_constraint_system) ====\n";
    test_domain_t top;
    check_bool("case12: top yields a true constraint system",
               top.to_linear_constraint_system().is_true(), true);

    test_domain_t dom;
    dom.add(v1, v2);
    dom.add(v2, v3); // one class {v1,v2,v3}: 3 members, 1 rep -> 2 equalities
    auto csts = dom.to_linear_constraint_system();
    unsigned n_eq = 0;
    for (auto &c : csts) {
      if (c.is_equality()) {
        n_eq++;
      }
    }
    check_bool("case12: {v1,v2,v3} yields 2 constraints", csts.size() == 2,
               true);
    check_bool("case12: all emitted constraints are equalities",
               n_eq == csts.size(), true);
  }

  { // an id read from one domain value can be shared into another value to
    // record a cross-value equality; fresh ids never land on it by accident
    crab::outs() << "==== case13 (share an id across domain values) ====\n";
    test_domain_t regs, flds;
    regs.add(v1, v2); // {v1,v2} with a factory-produced id
    value_domain_t id = *regs.get_class_id(v1);
    flds.set(v9, id);   // deliberately share the id into another value
    flds.add(v10, v11); // fresh ids in flds cannot collide with `id`
    check_partition("case13: sharing links only the intended class", flds,
                    {v9, v10, v11}, {{v10, v11}});
    check_bool("case13: shared id matches across the two values",
               *flds.get_class_id(v9) == id, true);
  }

  { // set() on an existing var to an id already owned by another class must
    // MERGE the two classes, not leave two classes sharing one id
    crab::outs() << "==== case14 (set merges on an in-use id) ====\n";
    value_domain_t l1 = test_domain_t::fresh_class_id();
    value_domain_t l2 = test_domain_t::fresh_class_id();
    test_domain_t dom;
    dom.set(v1, l1);
    dom.add(v1, v2); // {v1,v2} tagged l1
    dom.set(v3, l2);
    dom.add(v3, v4); // {v3,v4} tagged l2
    dom.set(v3, l1); // relabel v3's class to an in-use id
    check_partition("case14: set to an in-use id merges the classes", dom,
                    {v1, v2, v3, v4}, {{v1, v2, v3, v4}});
  }

  crab::outs() << "\n[SUMMARY] " << (g_checks - g_failures) << "/" << g_checks
               << " checks passed\n";
  if (g_failures > 0) {
    crab::outs() << "[SUMMARY] " << g_failures << " check(s) FAILED\n";
    return 1;
  }
  crab::outs() << "[SUMMARY] all checks passed\n";
  return 0;
}
