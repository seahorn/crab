// Unit tests for exporting abstract states as linear constraint systems,
// organized as Boost.Test suites:
//
//   flat_bool_export     to_linear_constraint_system on the 3-valued boolean
//                        domain: top is the empty system, bottom carries a
//                        contradiction, known values become 0/1 equalities
//   product_export       to_disjunctive_linear_constraint_system on
//                        flat_boolean_numerical_domain: a product denotes the
//                        *conjunction* of its components
//   array_export         to_disjunctive_linear_constraint_system on the array
//                        wrappers, which must propagate a bottom base domain
//                        rather than iterating it
//
// The product and array suites are regression tests: exporting a state whose
// boolean component was top used to terminate the process ("cannot add true"),
// and a state with both boolean and numerical information used to be exported
// as their disjunction instead of their conjunction.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_lincst_export
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../common.hpp"

#include <crab/support/debug.hpp>

#include <set>
#include <string>
#include <vector>

using namespace crab::cfg_impl;
using namespace crab::domain_impl;

using flat_bool_domain_t =
    crab::domains::flat_boolean_domain<ikos::z_number, varname_t>;
using product_domain_t = z_bool_interval_domain_t;
using z_lin_cst_sys_t = ikos::linear_constraint_system<ikos::z_number, varname_t>;
using z_disj_lin_cst_sys_t =
    ikos::disjunctive_linear_constraint_system<ikos::z_number, varname_t>;

namespace {

/** The constraints of a system, as strings, so they can be compared by set. */
std::set<std::string> constraints_of(const z_lin_cst_sys_t &csts) {
  std::set<std::string> res;
  for (auto const &c : csts) {
    crab::crab_string_os os;
    os << c;
    res.insert(os.str());
  }
  return res;
}

/** The disjuncts of a disjunctive system. Must not be called on bottom. */
std::vector<std::set<std::string>>
disjuncts_of(const z_disj_lin_cst_sys_t &dcsts) {
  std::vector<std::set<std::string>> res;
  for (auto const &csts : dcsts) {
    res.push_back(constraints_of(csts));
  }
  return res;
}

} // namespace

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(flat_bool_export)

// Top is the empty system, matching the convention of the numerical domains
// (e.g. split_dbm): is_true() is "no constraints". Returning a system that
// *contains* a tautology instead would not be recognized as top, and would add
// a redundant `true` to any conjunction it took part in.
BOOST_AUTO_TEST_CASE(top_exports_the_empty_system) {
  flat_bool_domain_t dom;
  auto csts = dom.to_linear_constraint_system();
  BOOST_CHECK_MESSAGE(csts.is_true(), "top must export as the empty system");
  BOOST_CHECK_EQUAL(csts.size(), 0u);
}

// Bottom, by contrast, *is* represented by a constraint, so that is_false()
// can detect it.
BOOST_AUTO_TEST_CASE(bottom_exports_a_contradiction) {
  flat_bool_domain_t dom;
  dom.set_to_bottom();
  auto csts = dom.to_linear_constraint_system();
  BOOST_CHECK(csts.is_false());
}

BOOST_AUTO_TEST_CASE(known_booleans_export_as_zero_one_equalities) {
  variable_factory_t vfac;
  z_var b(vfac["b"], crab::BOOL_TYPE);
  z_var c(vfac["c"], crab::BOOL_TYPE);

  flat_bool_domain_t dom;
  dom.assume_bool(b, false /*not negated*/); // b is true
  dom.assume_bool(c, true /*negated*/);      // c is false

  auto csts = constraints_of(dom.to_linear_constraint_system());
  BOOST_CHECK_EQUAL(csts.count("b = 1"), 1u);
  BOOST_CHECK_EQUAL(csts.count("c = 0"), 1u);
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(product_export)

BOOST_AUTO_TEST_CASE(top_exports_as_top) {
  product_domain_t dom;
  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_CHECK(!dcsts.is_false());
  BOOST_CHECK_MESSAGE(dcsts.is_true(), "a top product must export as top");
}

BOOST_AUTO_TEST_CASE(bottom_exports_as_bottom) {
  product_domain_t dom;
  dom.set_to_bottom();
  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_CHECK_MESSAGE(dcsts.is_false(),
                      "a bottom product must export as bottom");
}

// Regression: with a top boolean component, exporting used to reach
// `operator+=(true)` on a disjunctive system, which terminates the process.
// Most program points have no boolean information, so this was the common case.
BOOST_AUTO_TEST_CASE(numerical_only_survives_a_top_boolean_component) {
  variable_factory_t vfac;
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  product_domain_t dom;
  dom += z_lin_cst_t(z_lin_exp_t(y) >= ikos::z_number(1));
  dom += z_lin_cst_t(z_lin_exp_t(y) <= ikos::z_number(10));

  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_REQUIRE(!dcsts.is_false());
  BOOST_CHECK_MESSAGE(!dcsts.is_true(),
                      "numerical information must not be exported as top");

  auto disjuncts = disjuncts_of(dcsts);
  BOOST_REQUIRE_EQUAL(disjuncts.size(), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("-y <= -1"), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("y <= 10"), 1u);
  // No stray tautology from the top boolean component.
  BOOST_CHECK_EQUAL(disjuncts[0].count("true"), 0u);
}

BOOST_AUTO_TEST_CASE(boolean_only_is_exported) {
  variable_factory_t vfac;
  z_var b(vfac["b"], crab::BOOL_TYPE);

  product_domain_t dom;
  dom.assume_bool(b, false /*not negated*/);

  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_REQUIRE(!dcsts.is_false());
  BOOST_REQUIRE(!dcsts.is_true());
  auto disjuncts = disjuncts_of(dcsts);
  BOOST_REQUIRE_EQUAL(disjuncts.size(), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("b = 1"), 1u);
}

// Regression, and the reason this file exists: a product denotes the
// conjunction of its components. Adding each component to a disjunctive system
// in turn computes their *disjunction* instead, which is far weaker -- here it
// would yield `(b = 1) OR (1 <= y <= 10)` as two disjuncts rather than one
// disjunct holding both facts.
BOOST_AUTO_TEST_CASE(components_are_conjoined_not_disjoined) {
  variable_factory_t vfac;
  z_var b(vfac["b"], crab::BOOL_TYPE);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  product_domain_t dom;
  dom.assume_bool(b, false /*not negated*/);
  dom += z_lin_cst_t(z_lin_exp_t(y) >= ikos::z_number(1));
  dom += z_lin_cst_t(z_lin_exp_t(y) <= ikos::z_number(10));

  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_REQUIRE(!dcsts.is_false());
  BOOST_REQUIRE(!dcsts.is_true());

  auto disjuncts = disjuncts_of(dcsts);
  BOOST_REQUIRE_MESSAGE(disjuncts.size() == 1u,
                        "the two components must be conjoined into a single "
                        "disjunct, not unioned into two");
  BOOST_CHECK_EQUAL(disjuncts[0].count("b = 1"), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("-y <= -1"), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("y <= 10"), 1u);
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(array_export)

// Regression: the array wrappers range-iterate the base domain's disjunctive
// system, and begin() raises an error on a bottom system. This was unreachable
// only while the product never reported bottom.
BOOST_AUTO_TEST_CASE(array_adaptive_propagates_bottom) {
  z_aa_bool_int_t dom;
  dom.set_to_bottom();
  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_CHECK(dcsts.is_false());
}

BOOST_AUTO_TEST_CASE(array_smashing_propagates_bottom) {
  z_as_bool_num_t dom;
  dom.set_to_bottom();
  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_CHECK(dcsts.is_false());
}

BOOST_AUTO_TEST_CASE(array_adaptive_exports_scalar_facts) {
  variable_factory_t vfac;
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  z_aa_bool_int_t dom;
  dom += z_lin_cst_t(z_lin_exp_t(y) >= ikos::z_number(1));
  dom += z_lin_cst_t(z_lin_exp_t(y) <= ikos::z_number(10));

  auto dcsts = dom.to_disjunctive_linear_constraint_system();
  BOOST_REQUIRE(!dcsts.is_false());
  BOOST_REQUIRE(!dcsts.is_true());
  auto disjuncts = disjuncts_of(dcsts);
  BOOST_REQUIRE_EQUAL(disjuncts.size(), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("-y <= -1"), 1u);
  BOOST_CHECK_EQUAL(disjuncts[0].count("y <= 10"), 1u);
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
