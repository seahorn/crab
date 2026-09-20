// Tests for the cast dispatch in flat_boolean_numerical_domain.
//
// The dispatch used to compare bitwidths. A BOOL_TYPE variable reports width
// 1, and so did an INT_TYPE variable of bitwidth 1, so `trunc i32 -> i1` took
// the int-to-bool branch. That branch records a boolean value for the
// destination and never writes the numerical factor, so an i1 *integer*
// destination kept its previous value instead of being assigned -- an unsound
// stale read. The guards now test types, which tells the two apart.
//
// The regression test for that case lived here and constructed an
// INT_TYPE of bitwidth 1 directly. variable_type now rejects that width, so
// the case is prevented at construction and can no longer be exercised
// in-process: CRAB_ERROR calls std::exit rather than throwing. The executable
// evidence for the fix is the parent commit, where the test failed before the
// dispatch change and passed after it.
//
// What remains pins the behaviour the type-based guards must not disturb, and
// keeps covering them for the destination types that are still constructible.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_flat_bool_cast
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../crab_dom.hpp"

using namespace crab::domain_impl;
using namespace ikos;

namespace {
struct fixture {
  variable_factory_t vfac;
  z_var x; // i32, the wide side of a cast
  z_var b; // a boolean
  fixture()
      : x(vfac["x"], crab::INT_TYPE, 32), b(vfac["b"], crab::BOOL_TYPE) {}
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(flat_bool_cast, fixture)

// A boolean destination keeps its precise handling: a source known to be zero
// makes the boolean false.
BOOST_AUTO_TEST_CASE(trunc_to_bool_zero_source_is_false) {
  z_bool_num_domain_t d;
  d.assign(x, z_number(0));
  d.apply(crab::domains::OP_TRUNC, b, x);

  z_bool_num_domain_t assume_b(d);
  assume_b.assume_bool(b, false /*is_negated*/);
  BOOST_TEST((assume_b.is_bottom()));

  z_bool_num_domain_t assume_not_b(d);
  assume_not_b.assume_bool(b, true /*is_negated*/);
  BOOST_TEST((!assume_not_b.is_bottom()));
}

// ... and a source known to be non-zero makes it true.
BOOST_AUTO_TEST_CASE(trunc_to_bool_nonzero_source_is_true) {
  z_bool_num_domain_t d;
  d.assign(x, z_number(5));
  d.apply(crab::domains::OP_TRUNC, b, x);

  z_bool_num_domain_t assume_not_b(d);
  assume_not_b.assume_bool(b, true /*is_negated*/);
  BOOST_TEST((assume_not_b.is_bottom()));

  z_bool_num_domain_t assume_b(d);
  assume_b.assume_bool(b, false /*is_negated*/);
  BOOST_TEST((!assume_b.is_bottom()));
}

// bool -> int with the boolean's value known: zext gives 1, sext gives -1.
BOOST_AUTO_TEST_CASE(zext_and_sext_from_a_known_bool) {
  {
    z_bool_num_domain_t d;
    d.assign_bool_cst(b, z_lin_cst_t::get_true());
    d.apply(crab::domains::OP_ZEXT, x, b);
    BOOST_TEST((d.at(x) == z_interval_t(z_number(1))));
  }
  {
    z_bool_num_domain_t d;
    d.assign_bool_cst(b, z_lin_cst_t::get_true());
    d.apply(crab::domains::OP_SEXT, x, b);
    BOOST_TEST((d.at(x) == z_interval_t(z_number(-1))));
  }
  {
    z_bool_num_domain_t d;
    d.assign_bool_cst(b, z_lin_cst_t::get_false());
    d.apply(crab::domains::OP_ZEXT, x, b);
    BOOST_TEST((d.at(x) == z_interval_t(z_number(0))));
  }
}

// bool -> int with the boolean unknown delegates to the numerical domain,
// where int_cast_domain_traits refines a zext of a boolean to [0,1]. This is
// the path a loop-carried boolean takes, and it is the one the type-based
// guard has to keep reaching.
BOOST_AUTO_TEST_CASE(zext_from_an_unknown_bool_is_refined_to_zero_one) {
  z_bool_num_domain_t d;
  d -= b; // havoc
  d.apply(crab::domains::OP_ZEXT, x, b);
  BOOST_TEST((d.at(x) == z_interval_t(z_number(0), z_number(1))));
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
