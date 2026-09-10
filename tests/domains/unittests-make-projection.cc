// Unit tests for the functional project/forget API
// (make_projection/make_forget), organized as Boost.Test suites:
//   functional_law   for every domain: make_projection(V) equals
//                    {copy; project(V)} by mutual inclusion (same for
//                    make_forget), over a relational state and the edge
//                    cases the implementations special-case: bottom, top,
//                    empty V, untracked variables, V == all tracked.
//   semantics        direct expectations on the specialized domains: a
//                    kept relation survives the projection, a dropped
//                    variable becomes unconstrained.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_make_projection
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>

#include "../crab_dom.hpp"
#include <crab/domains/sign_domain.hpp>

using namespace crab::domain_impl;
using namespace ikos;

using z_sign_domain_t = crab::domains::sign_domain<z_number, varname_t>;

namespace {
// The domains under test: every specialized implementation that compiles
// in this configuration plus representatives of the default macro
// (intervals & friends) and of the forwarding wrappers (bool+num, ric).
using test_domains = boost::mpl::list<z_interval_domain_t,
                                      z_constant_domain_t, z_sign_domain_t,
                                      z_ric_domain_t, z_dis_interval_domain_t,
                                      z_sdbm_domain_t, z_soct_domain_t,
                                      z_term_domain_t, z_bool_num_domain_t
#ifdef HAVE_ELINA
                                      ,
                                      z_oct_elina_domain_t,
                                      z_pk_elina_domain_t
#endif
                                      >;

struct fixture {
  variable_factory_t vfac;
  z_var x, y, z, w;
  fixture()
      : x(vfac["x"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32),
        z(vfac["z"], crab::INT_TYPE, 32), w(vfac["w"], crab::INT_TYPE, 32) {}

  // relational state over {x, y, z}; w stays untracked
  template <typename Dom> Dom mk_state() const {
    Dom d;
    d += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
    d += z_lin_cst_t(z_lin_exp_t(x) <= z_number(9));
    d += z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2));
    d += z_lin_cst_t(z_lin_exp_t(z) <= z_lin_exp_t(y));
    return d;
  }
};

template <typename Dom>
bool same_meaning(const Dom &a, const Dom &b) {
  return a <= b && b <= a;
}
} // namespace

BOOST_FIXTURE_TEST_SUITE(functional_law, fixture)

BOOST_AUTO_TEST_CASE_TEMPLATE(equals_copy_then_inplace, Dom, test_domains) {
  std::vector<std::vector<z_var>> vsets = {
      {x},       // strict subset
      {x, y},    // subset keeping a relation
      {x, y, z}, // all tracked
      {},        // empty
      {w},       // untracked only
      {x, w},    // mixed tracked/untracked
  };
  for (auto const &vs : vsets) {
    { // projection law on a relational state
      Dom d = mk_state<Dom>();
      Dom ref(d);
      ref.project(vs);
      BOOST_TEST((same_meaning(d.make_projection(vs), ref)));
    }
    { // forget law on a relational state
      Dom d = mk_state<Dom>();
      Dom ref(d);
      ref.forget(vs);
      BOOST_TEST((same_meaning(d.make_forget(vs), ref)));
    }
    { // bottom
      Dom d = mk_state<Dom>();
      d += z_lin_cst_t(z_lin_exp_t(x) >= z_number(100)); // contradiction
      Dom rp(d);
      rp.project(vs);
      Dom rf(d);
      rf.forget(vs);
      BOOST_TEST((same_meaning(d.make_projection(vs), rp)));
      BOOST_TEST((same_meaning(d.make_forget(vs), rf)));
    }
    { // top
      Dom d;
      Dom rp(d);
      rp.project(vs);
      Dom rf(d);
      rf.forget(vs);
      BOOST_TEST((same_meaning(d.make_projection(vs), rp)));
      BOOST_TEST((same_meaning(d.make_forget(vs), rf)));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(semantics, fixture)

// split_dbm: the projection keeps the relation among kept variables --
// including one that only holds THROUGH a dropped variable -- and the
// forgotten variable becomes unconstrained.
BOOST_AUTO_TEST_CASE(sdbm_projection_keeps_relations) {
  z_sdbm_domain_t d;
  d += z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2));
  d += z_lin_cst_t(z_lin_exp_t(z) == z_lin_exp_t(y) + z_number(3));
  // x--z relation flows through the dropped y
  z_sdbm_domain_t p = d.make_projection({x, z});
  BOOST_TEST(
      (p.entails(z_lin_cst_t(z_lin_exp_t(z) == z_lin_exp_t(x) + z_number(5)))));
  BOOST_TEST((p.at(y).is_top()));
  z_sdbm_domain_t f = d.make_forget({y});
  BOOST_TEST(
      (f.entails(z_lin_cst_t(z_lin_exp_t(z) == z_lin_exp_t(x) + z_number(5)))));
  BOOST_TEST((f.at(y).is_top()));
  // the source state is untouched (functional, not in-place)
  BOOST_TEST(
      (d.entails(z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2)))));
}

// split_oct: variable bounds live in the v+/v- pair encoding; they must
// survive the projection.
BOOST_AUTO_TEST_CASE(soct_projection_keeps_bounds_and_relations) {
  z_soct_domain_t d;
  d += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
  d += z_lin_cst_t(z_lin_exp_t(x) <= z_number(9));
  d += z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2));
  z_soct_domain_t p = d.make_projection({x, y});
  BOOST_TEST((p.entails(z_lin_cst_t(z_lin_exp_t(x) >= z_number(1)))));
  BOOST_TEST((p.entails(z_lin_cst_t(z_lin_exp_t(x) <= z_number(9)))));
  BOOST_TEST(
      (p.entails(z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2)))));
  z_soct_domain_t q = d.make_projection({y});
  BOOST_TEST((q.entails(z_lin_cst_t(z_lin_exp_t(y) >= z_number(3)))));
  BOOST_TEST((q.at(x).is_top()));
}

#ifdef HAVE_ELINA
// elina: the functional forms must leave the source untouched and agree
// with the in-place results (the wrap-the-fresh-state path).
BOOST_AUTO_TEST_CASE(elina_functional_leaves_source_untouched) {
  z_oct_elina_domain_t d;
  d += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
  d += z_lin_cst_t(z_lin_exp_t(y) == z_lin_exp_t(x) + z_number(2));
  z_oct_elina_domain_t p = d.make_projection({y});
  BOOST_TEST((p.entails(z_lin_cst_t(z_lin_exp_t(y) >= z_number(3)))));
  BOOST_TEST((p.at(x).is_top()));
  z_oct_elina_domain_t f = d.make_forget({x});
  BOOST_TEST((f.at(x).is_top()));
  BOOST_TEST((d.entails(z_lin_cst_t(z_lin_exp_t(x) >= z_number(1)))));
}
#endif

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
