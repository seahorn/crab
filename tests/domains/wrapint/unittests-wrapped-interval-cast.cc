// Regression tests for truncation in wrapped_interval_domain.
//
// The OP_TRUNC case of apply(int_conv_operation_t, ...) declared a second
// `dst_i` that shadowed the one the function stores, so the truncated value
// was computed into a local and discarded. `set(dst, dst_i)` then wrote the
// outer, default-constructed wrapped_interval -- which is [0,7] at 3 bits and
// satisfies is_top(), so the binding was erased instead of assigned. Every
// truncation of a known range produced top for the destination: sound, but a
// total loss of precision in the one domain whose whole purpose is modelling
// machine arithmetic. The zext/sext case was never affected.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_wrapped_interval_cast
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../../crab_dom.hpp"

using namespace crab::domain_impl;
using namespace ikos;

namespace {
struct fixture {
  variable_factory_t vfac;
  z_var x;  // i32 source
  z_var y;  // i8 destination
  fixture()
      : x(vfac["x"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 8) {}

  // x constrained to [lo, hi], then y := trunc(x)
  z_wrapped_interval_domain_t trunc_from(int lo, int hi) const {
    z_wrapped_interval_domain_t d;
    d += z_lin_cst_t(z_lin_exp_t(x) >= z_number(lo));
    d += z_lin_cst_t(z_lin_exp_t(x) <= z_number(hi));
    d.apply(crab::domains::OP_TRUNC, y, x);
    return d;
  }
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(wrapped_interval_cast, fixture)

// A range that fits in the destination and does not cross the signed limit
// transfers exactly. Before the fix the destination was top.
BOOST_AUTO_TEST_CASE(trunc_keeps_a_range_that_fits) {
  z_wrapped_interval_domain_t d = trunc_from(10, 20);
  BOOST_TEST((d.at(y) == z_interval_t(z_number(10), z_number(20))));
}

// A value that does not fit wraps, rather than being dropped: 300 mod 2^8.
BOOST_AUTO_TEST_CASE(trunc_wraps_a_value_that_does_not_fit) {
  z_wrapped_interval_domain_t d = trunc_from(300, 300);
  BOOST_TEST((d.at(y) == z_interval_t(z_number(44))));
}

// Truncation must still lose nothing that matters when the result spans the
// signed boundary: [100,200] as an i8 is the wrapped interval 100..200, whose
// signed reading crosses 127, so at() is top -- but the domain still holds the
// information, and a later meet sees it.
BOOST_AUTO_TEST_CASE(trunc_across_the_signed_limit_keeps_information) {
  z_wrapped_interval_domain_t d = trunc_from(100, 200);
  BOOST_TEST((d.at(y).is_top()));
  // signed values of {100..200} are {100..127} U {-128..-56}; meeting with
  // y <= 3 leaves only the negative part. Before the fix y was unconstrained
  // and the meet gave [-128, 3].
  d += z_lin_cst_t(z_lin_exp_t(y) <= z_number(3));
  BOOST_TEST((d.at(y) == z_interval_t(z_number(-128), z_number(-56))));
}

// The sext/zext path never had the shadowing bug; pin it so a future edit to
// this function cannot silently break it the same way.
BOOST_AUTO_TEST_CASE(zext_still_widens) {
  z_wrapped_interval_domain_t d;
  d += z_lin_cst_t(z_lin_exp_t(y) >= z_number(10));
  d += z_lin_cst_t(z_lin_exp_t(y) <= z_number(20));
  d.apply(crab::domains::OP_ZEXT, x, y);
  BOOST_TEST((d.at(x) == z_interval_t(z_number(10), z_number(20))));
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
