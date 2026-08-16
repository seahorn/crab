// End-to-end tests for the tDBM domain: each program is analyzed by the
// intra-procedural fixpoint and its __CRAB_assert statements are then
// discharged by the assertion checker.  Every assertion carries a debug_info
// id, so the tests pin the verdict of each individual assertion as well as
// the per-program totals.
//
// Nothing is printed to standard output: the whole Boost.Test log is routed
// to stderr and ctest checks the exit code.
#define BOOST_TEST_MODULE tvpi_cfg
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../../common.hpp"

#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>
#include <crab/support/debug.hpp>

#include <memory>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::checker;
using namespace ikos;

namespace {

// Run the forward analyzer (same fixpoint knobs the old run_and_check calls
// used: widening delay 2, 1 descending iteration, jump set size 20) and
// discharge the program's assertions.  Returns the checker's database.
checks_db analyze_and_check(z_cfg_t &cfg) {
  using analyzer_t =
      crab::analyzer::intra_fwd_analyzer<z_cfg_ref_t, z_tvpi_dbm_domain_t>;
  z_cfg_ref_t cfg_ref(cfg);
  z_tvpi_dbm_domain_t init;
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = 2;
  fixpo_params.get_descending_iterations() = 1;
  fixpo_params.get_max_thresholds() = 20;
  analyzer_t analyzer(cfg_ref, init.make_top(), nullptr /*no liveness*/,
                      fixpo_params);
  typename analyzer_t::assumption_map_t assumptions;
  analyzer.run(cfg.entry(), init, assumptions);

  using checker_t = intra_checker<analyzer_t>;
  using prop_t = assert_property_checker<analyzer_t>;
  typename checker_t::prop_checker_ptr prop(new prop_t(0 /*verbose*/));
  checker_t checker(analyzer, {prop});
  checker.run();
  return checker.get_all_checks();
}

// Every check recorded for assertion @p id is CRAB_SAFE (and at least one
// check was recorded: an unreachable-and-unchecked assertion must not count
// as proven).
bool assertion_safe(const checks_db &db, int64_t id) {
  crab::cfg::debug_info di(id);
  if (!db.has_checks(di)) {
    return false;
  }
  for (check_kind k : db.get_checks(di)) {
    if (k != check_kind::CRAB_SAFE && k != check_kind::CRAB_UNREACH) {
      return false;
    }
  }
  return true;
}

// The global coefficient template, scoped to a fixture.
struct scoped_template {
  std::vector<unsigned> m_saved;
  explicit scoped_template(std::initializer_list<unsigned> t)
      : m_saved(crab_domain_params_man::get().coefficients()) {
    crab_domain_params_man::get().coefficients().assign(t);
  }
  ~scoped_template() { crab_domain_params_man::get().coefficients() = m_saved; }
};

struct template_234 {
  scoped_template tpl;
  variable_factory_t vfac;
  template_234() : tpl({2, 3, 4}) {}
};

//===----------------------------------------------------------------------===//
// The analyzed programs.  Assertion ids: <prog><index>.
//===----------------------------------------------------------------------===//

// int N = nd_int(); assume(N >= 1);
// i = x = y = 0;
// while (i < N) { i++; x += 4; y += 8; }
// assert(x == 4*N);  // 11: needs cross-coefficient saturation
// assert(y == 8*N);  // 12: needs cross-coefficient saturation
std::unique_ptr<z_cfg_t> prog1(variable_factory_t &vfac) {
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var n(vfac["N"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_header;
  loop_exit >> exit;
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  entry.assign(y, 0);
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x == 4 * n, debug_info(11));
  loop_exit.assertion(y == 8 * n, debug_info(12));
  loop_body.intrinsic("loop_counter", {}, {i});
  loop_body.add(i, i, 1);
  loop_body.add(x, x, 4);
  loop_body.add(y, y, 8);
  return cfg;
}

// int N = nd_int(); assume(N >= 1);
// i = x = 0;
// while (i < N) { i++; if (*) x += 2; else x += 3; }
// assert(x >= 2*N);  // 21: needs cross-coefficient saturation
// assert(x <= 3*N);  // 22: needs cross-coefficient saturation
std::unique_ptr<z_cfg_t> prog2(variable_factory_t &vfac) {
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_body_then = cfg->insert("loop_body_then");
  z_basic_block_t &loop_body_else = cfg->insert("loop_body_else");
  z_basic_block_t &loop_body_tail = cfg->insert("loop_body_tail");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_body_then;
  loop_body >> loop_body_else;
  loop_body_then >> loop_body_tail;
  loop_body_else >> loop_body_tail;
  loop_body_tail >> loop_header;
  loop_exit >> exit;
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x >= 2 * n, debug_info(21));
  loop_exit.assertion(x <= 3 * n, debug_info(22));
  loop_body.intrinsic("loop_counter", {}, {i});
  loop_body.add(i, i, 1);
  loop_body_then.add(x, x, 2);
  loop_body_else.add(x, x, 3);
  return cfg;
}

// int isz = 4;
// int len = nd_int(); assume(1 <= len <= 10);
// int tsz = nd_int(); assume(tsz >= len * isz);
// for (i = 0; i < len; i++) {
//   idx = i * isz; offset = idx + isz;
//   assert(offset <= tsz);  // 31
// }
// assert(tsz >= len * isz);  // 32
std::unique_ptr<z_cfg_t> prog3(variable_factory_t &vfac) {
  z_var isz(vfac["isz"], crab::INT_TYPE, 32);
  z_var len(vfac["len"], crab::INT_TYPE, 32);
  z_var tsz(vfac["tsz"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var idx(vfac["idx"], crab::INT_TYPE, 32);
  z_var offset(vfac["offset"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> loop_header;
  loop_header >> loop_body;
  loop_body >> loop_header;
  loop_header >> loop_exit;
  loop_exit >> exit;
  entry.assign(isz, 4);
  entry.havoc(len);
  entry.assume(len >= 1);
  entry.assume(len <= 10);
  entry.havoc(tsz);
  entry.mul(tmp, len, isz);
  entry.assume(z_lin_exp_t(tsz) >= tmp);
  entry.assign(i, 0);
  loop_body.assume(z_lin_exp_t(i) < len);
  loop_body.mul(idx, i, isz);
  loop_body.add(offset, idx, isz);
  loop_body.assertion(z_lin_exp_t(offset) <= tsz, debug_info(31));
  loop_body.add(i, i, 1);
  loop_exit.assume(z_lin_exp_t(i) >= len);
  loop_exit.mul(tmp2, len, isz);
  loop_exit.assertion(z_lin_exp_t(tsz) >= tmp2, debug_info(32));
  return cfg;
}

// int x = nd_int(); y = 2*x;
// int z = nd_int(); k = 2*z;
// assume(z - x >= 3);
// assert(y == 2*x);     // 41
// assert(k - y >= 6);   // 42: 2x - k <= -6, so y - k <= -6
// assert(2*z - y >= 6); // 43: y - 2z <= -6
std::unique_ptr<z_cfg_t> prog4(variable_factory_t &vfac) {
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header_1 = cfg->insert("header_1");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> header_1;
  header_1 >> exit;
  entry.havoc(x);
  entry.mul(y, x, 2);
  entry.havoc(z);
  entry.mul(k, z, 2);
  header_1.assume(z - x >= 3);
  exit.assertion(y == 2 * x, debug_info(41));
  exit.assertion(k - y >= 6, debug_info(42));
  exit.assertion(2 * z - y >= 6, debug_info(43));
  return cfg;
}

// assume(0 <= x <= 4); assume(y == 2);
// z = x * y;            // variable * variable, y a singleton
// int m = nd_int(); assume(m >= z);
// j = x - 1; k = j * y;
// assert(k + y <= m);   // 51
std::unique_ptr<z_cfg_t> prog5(variable_factory_t &vfac) {
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var m(vfac["m"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb_1 = cfg->insert("bb_1");
  z_basic_block_t &bb_2 = cfg->insert("bb_2");
  z_basic_block_t &bb_3 = cfg->insert("bb_3");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> bb_1;
  bb_1 >> bb_2;
  bb_2 >> bb_3;
  bb_3 >> exit;
  entry.havoc(x);
  entry.assume(x >= 0);
  entry.assume(x <= 4);
  entry.havoc(y);
  entry.assume(y == 2);
  bb_1.havoc(z);
  bb_1.mul(z, x, y);
  bb_2.havoc(m);
  bb_2.assume(m >= z);
  bb_3.sub(j, x, 1);
  bb_3.mul(k, j, y);
  exit.assertion(k + y <= m, debug_info(51));
  return cfg;
}

// y = 2*n; x = y*1; z = x; m = y/1; z' = m;
// assert(z == 2*n);   // 61: equality preserved through *1 and copies
// assert(z' == 2*n);  // 62: equality preserved through /1 and copies
std::unique_ptr<z_cfg_t> prog6(variable_factory_t &vfac) {
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var z_p(vfac["z'"], crab::INT_TYPE, 32);
  z_var m(vfac["m"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("entry", "exit"));
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header_1 = cfg->insert("header_1");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> header_1;
  header_1 >> exit;
  entry.havoc(n);
  entry.havoc(x);
  entry.havoc(z);
  entry.havoc(m);
  entry.mul(y, n, 2);
  header_1.mul(x, y, 1);
  header_1.assign(z, x);
  header_1.div(m, y, 1);
  header_1.assign(z_p, m);
  exit.assertion(z == 2 * n, debug_info(61));
  exit.assertion(z_p == 2 * n, debug_info(62));
  return cfg;
}

} // namespace

//===----------------------------------------------------------------------===//
// Verdicts.  All assertions in these programs are expected to be proven.
//===----------------------------------------------------------------------===//

BOOST_FIXTURE_TEST_SUITE(tvpi_cfg_checker, template_234)

BOOST_AUTO_TEST_CASE(prog1_loop_scaling) {
  auto cfg = prog1(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 2u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 11), "x == 4*N after the loop");
  BOOST_TEST(assertion_safe(db, 12), "y == 8*N after the loop");
}

BOOST_AUTO_TEST_CASE(prog2_loop_two_strides) {
  auto cfg = prog2(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 2u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 21), "x >= 2*N after the loop");
  BOOST_TEST(assertion_safe(db, 22), "x <= 3*N after the loop");
}

BOOST_AUTO_TEST_CASE(prog3_array_bounds) {
  auto cfg = prog3(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 2u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 31), "offset <= tsz inside the loop");
  BOOST_TEST(assertion_safe(db, 32), "tsz >= len*isz after the loop");
}

BOOST_AUTO_TEST_CASE(prog4_scaled_differences) {
  auto cfg = prog4(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 3u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 41), "y == 2*x");
  BOOST_TEST(assertion_safe(db, 42), "k - y >= 6");
  BOOST_TEST(assertion_safe(db, 43), "2*z - y >= 6");
}

BOOST_AUTO_TEST_CASE(prog5_var_var_multiplication) {
  auto cfg = prog5(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 1u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 51), "k + y <= m");
}

BOOST_AUTO_TEST_CASE(prog6_unit_mul_div) {
  auto cfg = prog6(vfac);
  checks_db db = analyze_and_check(*cfg);
  BOOST_TEST(db.get_total_safe() == 2u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_safe(db, 61), "z == 2*n through x := y*1");
  BOOST_TEST(assertion_safe(db, 62), "z' == 2*n through m := y/1");
}

BOOST_AUTO_TEST_SUITE_END()

// Keep stdout empty: the golden-output harness (tests/run_tests.sh) diffs
// the standard output of every test binary, so the whole Boost.Test log is
// routed to stderr.  ctest passes --disable-warnings, which Boost.Test would
// reject, so it is translated into the crab flag it stands for.
int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int i = 0; i < argc; ++i) {
    if (std::string(argv[i]) == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    args.push_back(argv[i]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
