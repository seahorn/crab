// Tests for MATH_INT_TYPE and ARR_MATH_INT_TYPE: integers, and arrays of
// integers, that carry no representation width.
//
// The type-level properties (predicates, equality, hashing, printing) are
// checked directly. The domain-level behaviour covers the two casts that are
// legal for a mathematical integer -- trunc to bool and zext from bool -- and
// an array round trip, which is the case that motivated a separate array kind:
// the array domains build each cell's ghost variable from the array's type, so
// an array of mathematical integers must be distinguishable from an array of
// iN.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr (tests/run_tests.sh diffs stdout against a golden file) and
// ctest checks the exit code.
#define BOOST_TEST_MODULE unittests_math_int_type
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../crab_dom.hpp"
#include <crab/support/os.hpp>

using namespace crab::domain_impl;
using namespace ikos;

namespace {
crab::variable_type math_int() { return crab::variable_type(crab::MATH_INT_TYPE); }
crab::variable_type int_ty(unsigned w) {
  return crab::variable_type(crab::INT_TYPE, w);
}

std::string to_str(const crab::variable_type &ty) {
  crab::crab_string_os os;
  os << ty;
  return os.str();
}

struct fixture {
  variable_factory_t vfac;
  z_var m, n;  // mathematical integers
  z_var b;     // boolean
  z_var x;     // i32
  fixture()
      : m(vfac["m"], math_int()), n(vfac["n"], math_int()),
        b(vfac["b"], crab::BOOL_TYPE), x(vfac["x"], crab::INT_TYPE, 32) {}
};
} // namespace

BOOST_AUTO_TEST_SUITE(math_int_type)

// is_integer() is inclusive so that the numerical domains treat a
// mathematical integer as the numeric scalar it is; is_fixed_width_integer()
// is the exclusive test for code that needs a width.
BOOST_AUTO_TEST_CASE(predicates) {
  BOOST_TEST((math_int().is_integer()));
  BOOST_TEST((math_int().is_math_integer()));
  BOOST_TEST((!math_int().is_fixed_width_integer()));
  BOOST_TEST((math_int().is_scalar()));
  BOOST_TEST((!math_int().is_real()));
  BOOST_TEST((!math_int().is_bool()));
  // is_integer(w) asks for a specific width, so it is false for every w
  BOOST_TEST((!math_int().is_integer(32)));
  BOOST_TEST((!math_int().is_integer(0)));

  BOOST_TEST((int_ty(32).is_integer()));
  BOOST_TEST((int_ty(32).is_fixed_width_integer()));
  BOOST_TEST((!int_ty(32).is_math_integer()));
}

// Arrays go the other way: is_integer_array() stays exclusive, because its
// callers mean "elements have a width".
BOOST_AUTO_TEST_CASE(array_predicates) {
  crab::variable_type ma(crab::ARR_MATH_INT_TYPE);
  crab::variable_type ia(crab::ARR_INT_TYPE);
  BOOST_TEST((ma.is_array()));
  BOOST_TEST((ma.is_math_integer_array()));
  BOOST_TEST((!ma.is_integer_array()));
  BOOST_TEST((ia.is_integer_array()));
  BOOST_TEST((!ia.is_math_integer_array()));
  BOOST_TEST((ia.is_array()));
}

// Distinct from every iN and from real; equal only to itself.
BOOST_AUTO_TEST_CASE(equality_and_hashing) {
  BOOST_TEST((math_int() == math_int()));
  BOOST_TEST((math_int().hash() == math_int().hash()));
  BOOST_TEST((math_int() != int_ty(8)));
  BOOST_TEST((math_int() != int_ty(32)));
  BOOST_TEST((math_int() != int_ty(64)));
  BOOST_TEST((math_int() != crab::variable_type(crab::REAL_TYPE)));
  BOOST_TEST((math_int() != crab::variable_type(crab::BOOL_TYPE)));
  BOOST_TEST((crab::variable_type(crab::ARR_MATH_INT_TYPE) !=
              crab::variable_type(crab::ARR_INT_TYPE)));
}

// A fixed-width integer prints its width, so a bare "int" is unambiguous.
BOOST_AUTO_TEST_CASE(printing) {
  BOOST_TEST((to_str(math_int()) == "int"));
  BOOST_TEST((to_str(int_ty(32)) == "int32"));
  BOOST_TEST((to_str(crab::variable_type(crab::ARR_MATH_INT_TYPE)) ==
              "arr(mathint)"));
  BOOST_TEST((to_str(crab::variable_type(crab::ARR_INT_TYPE)) == "arr(int)"));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(math_int_domains, fixture)

// Ordinary arithmetic is unaffected: every domain already models iN as an
// unbounded integer, so a mathematical integer behaves the same way.
BOOST_AUTO_TEST_CASE(arithmetic_works) {
  z_bool_num_domain_t d;
  d.assign(m, z_number(10));
  d.apply(crab::domains::OP_ADDITION, n, m, z_number(5));
  BOOST_TEST((d.at(n) == z_interval_t(z_number(15))));
}

// Constraints over a mathematical integer. This goes through
// linear_interval_solver, whose get_bitwidth() reports 0 for a variable with
// no width; the default mk_interval ignores that, so the constraint applies
// normally. Only wrapped intervals specialize mk_interval, and they reject
// mathematical integers outright.
BOOST_AUTO_TEST_CASE(constraints_work) {
  z_bool_num_domain_t d;
  d += z_lin_cst_t(z_lin_exp_t(m) >= z_number(3));
  d += z_lin_cst_t(z_lin_exp_t(m) <= z_number(9));
  BOOST_TEST((d.at(m) == z_interval_t(z_number(3), z_number(9))));
  // and a relational constraint between two of them
  d += z_lin_cst_t(z_lin_exp_t(n) == z_lin_exp_t(m) + z_number(1));
  BOOST_TEST((d.at(n) == z_interval_t(z_number(4), z_number(10))));
}

// zext bool -> math int: the 0/1 embedding, with the boolean known and
// unknown. The unknown case goes through int_cast_domain_traits, whose zext
// refinement must still fire for a Boolean source.
BOOST_AUTO_TEST_CASE(bool_to_math_int) {
  {
    z_bool_num_domain_t d;
    d.assign_bool_cst(b, z_lin_cst_t::get_true());
    d.apply(crab::domains::OP_ZEXT, m, b);
    BOOST_TEST((d.at(m) == z_interval_t(z_number(1))));
  }
  {
    z_bool_num_domain_t d;
    d.assign_bool_cst(b, z_lin_cst_t::get_false());
    d.apply(crab::domains::OP_ZEXT, m, b);
    BOOST_TEST((d.at(m) == z_interval_t(z_number(0))));
  }
  {
    z_bool_num_domain_t d;
    d -= b; // unknown boolean
    d.apply(crab::domains::OP_ZEXT, m, b);
    BOOST_TEST((d.at(m) == z_interval_t(z_number(0), z_number(1))));
  }
}

// trunc math int -> bool: zero is false, non-zero is true.
BOOST_AUTO_TEST_CASE(math_int_to_bool) {
  {
    z_bool_num_domain_t d;
    d.assign(m, z_number(0));
    d.apply(crab::domains::OP_TRUNC, b, m);
    z_bool_num_domain_t assume_b(d);
    assume_b.assume_bool(b, false /*is_negated*/);
    BOOST_TEST((assume_b.is_bottom()));
  }
  {
    z_bool_num_domain_t d;
    d.assign(m, z_number(42));
    d.apply(crab::domains::OP_TRUNC, b, m);
    z_bool_num_domain_t assume_not_b(d);
    assume_not_b.assume_bool(b, true /*is_negated*/);
    BOOST_TEST((assume_not_b.is_bottom()));
  }
}

// The zext refinement must not invent a bound for a source that has no width.
// For an i32 source it still adds dst <= 2^32-1.
BOOST_AUTO_TEST_CASE(zext_refinement_only_for_fixed_width_sources) {
  {
    z_bool_num_domain_t d;
    d -= x; // unknown i32
    d.apply(crab::domains::OP_ZEXT, m, x);
    BOOST_TEST((d.at(m).ub().is_finite()));
  }
  {
    z_bool_num_domain_t d;
    d -= n; // unknown math int
    d.apply(crab::domains::OP_ZEXT, m, n);
    BOOST_TEST((d.at(m).is_top()));
  }
}

BOOST_AUTO_TEST_SUITE_END()

// Arrays of mathematical integers: the element size still says how many bytes
// an access covers, but the value stored has no width, so the cell ghost must
// not claim one. Before ARR_MATH_INT_TYPE existed, a load built an iN ghost
// from the element size and array_adaptive rejected the assignment as a type
// mismatch.
BOOST_AUTO_TEST_CASE(math_int_array_round_trip) {
  variable_factory_t vfac;
  z_var a(vfac["a"], crab::variable_type(crab::ARR_MATH_INT_TYPE));
  z_var v(vfac["v"], crab::variable_type(crab::MATH_INT_TYPE));
  z_var w(vfac["w"], crab::variable_type(crab::MATH_INT_TYPE));

  z_aa_bool_int_t d;
  d.assign(v, z_number(7));
  d.array_store(a, z_number(4) /*elem_size*/, z_number(0) /*index*/, v,
                true /*strong update*/);
  d.array_load(w, a, z_number(4), z_number(0));
  BOOST_TEST((d.at(w) == z_interval_t(z_number(7))));
}

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
