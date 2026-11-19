// Tests for the fun_array domain on whole programs: each test builds a CFG
// and runs the Crab fixpoint on it. Single operations are tested in
// unittests-fun-array.cc.
//
// The loops are Fig. 2 (initialization) and Fig. 3 (partial
// initialization) of P. Cousot, R. Cousot, F. Logozzo, "A Parametric
// Segmentation Functor for Fully Automatic and Scalable Array Content
// Analysis", POPL 2011. Segmentations here have no array length: the
// paper's last bound ({A.Length}) becomes the implicit +oo end.
//
// Nothing is printed to standard output unless --print-invariants is given:
// the Boost.Test log goes to stderr (tests/run_tests.sh diffs stdout against
// a golden file) and ctest checks the exit code.
#define BOOST_TEST_MODULE fun_array_fixpoint
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "../crab_lang.hpp"
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/domains/fun_array.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/support/os.hpp>

#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domains;
using namespace ikos;

namespace {
using scalar_t = split_dbm_domain<z_number, varname_t>;
using dom_t    = fun_array<scalar_t>;
using seg_t    = dom_t::segmentation_t;
using bound_t  = seg_t::bound_t;
using expr_t   = dom_t::expression_t;
using var_t    = dom_t::variable_t;
using lin_t    = dom_t::linear_expression_t;
using cst_t    = dom_t::linear_constraint_t;

// Set by --print-invariants
bool print_invariants = false;

void print(const std::string &name, const dom_t &d) {
  if (print_invariants) {
    crab::outs() << name << " = " << d << "\n";
  }
}

lin_t k(int x) { return lin_t(z_number(x)); }

bound_t B(std::initializer_list<lin_t> es) {
  bound_t res;
  for (auto const &e : es) {
    res.add(*expr_t::make(e));
  }
  return res;
}

bool same_bound(const bound_t &x, const bound_t &y) {
  return x <= y && y <= x;
}

// The bounds and ? marks of s are bounds and maybe_empty
void check_bounds(const seg_t &s, const std::vector<bound_t> &bounds,
                  const std::vector<bool> &maybe_empty) {
  BOOST_TEST_REQUIRE(!s.is_bottom());
  BOOST_TEST_REQUIRE(s.size() == bounds.size());
  for (std::size_t n = 0; n < bounds.size(); ++n) {
    BOOST_TEST(same_bound(s.bound(n), bounds[n]), "bound " << n);
    BOOST_TEST(s.maybe_empty(n) == maybe_empty[n], "? after bound " << n);
  }
}

// x := A[idx] in a copy of d, with the constraints csts on idx, entails
// x == val
bool loads(const dom_t &d, const var_t &x, const var_t &A, const var_t &idx,
           std::vector<cst_t> csts, int val) {
  dom_t t(d);
  for (auto const &c : csts) {
    t += c;
  }
  t.array_load(x, A, k(4), idx);
  return t.entails(cst_t(lin_t(x) - val, cst_t::EQUALITY));
}

struct fixture {
  variable_factory_t vfac;
  var_t A, C, i, j, n, t, x, idx;
  fixture()
      : A(vfac["A"], crab::ARR_INT_TYPE, 32), C(vfac["C"], crab::ARR_INT_TYPE, 32),
        i(vfac["i"], crab::INT_TYPE, 32), j(vfac["j"], crab::INT_TYPE, 32),
        n(vfac["n"], crab::INT_TYPE, 32), t(vfac["t"], crab::INT_TYPE, 32),
        x(vfac["x"], crab::INT_TYPE, 32), idx(vfac["idx"], crab::INT_TYPE, 32) {}

  // The invariants of each block: "b" at its entry, "b.post" at its exit
  std::map<std::string, dom_t> analyze(z_cfg_t &cfg) {
    using analyzer_t = crab::analyzer::intra_fwd_analyzer<z_cfg_ref_t, dom_t>;
    z_cfg_ref_t cfg_ref(cfg);
    crab::fixpoint_parameters params;
    analyzer_t a(cfg_ref, dom_t(), nullptr, params);
    typename analyzer_t::assumption_map_t assumptions;
    a.run(cfg.entry(), dom_t(), assumptions);
    std::map<std::string, dom_t> res;
    for (auto const &bb : cfg) {
      res.emplace(bb.label(), a.get_pre(bb.label()));
      res.emplace(bb.label() + ".post", a.get_post(bb.label()));
    }
    return res;
  }
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(paper_figures, fixture)

// Fig. 2:  i := 0; while (i < n) { A[i] := 0; i := i + 1; }
// At the exit A is {0} 0 {i}? T, so A[k] == 0 for every 0 <= k < n.
BOOST_AUTO_TEST_CASE(fig2_initialization) {
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit");
  BB(cfg, entry);
  BB(cfg, head);
  BB(cfg, body);
  BB(cfg, exit);
  entry >> head;
  head >> body;
  body >> head;
  head >> exit;
  entry.assign(i, 0);
  body.assume(lin_t(i) <= lin_t(n) - 1);
  body.array_store(A, i, 0, 4);
  body.add(i, i, 1);
  exit.assume(lin_t(i) >= lin_t(n));

  auto inv = analyze(*cfg);
  print("fig 2 head", inv.at("head"));
  print("fig 2 exit", inv.at("exit.post"));

  check_bounds(inv.at("head").segmentation_of(A), {B({k(0)}), B({i})},
               {false, true});
  const dom_t &out = inv.at("exit.post");
  check_bounds(out.segmentation_of(A), {B({k(0)}), B({i})}, {false, true});
  BOOST_TEST(loads(out, x, A, idx,
                   {cst_t(-lin_t(idx), cst_t::INEQUALITY),
                    cst_t(lin_t(idx) - n + 1, cst_t::INEQUALITY)},
                   0));
  // nothing is known at or after i
  BOOST_TEST(!loads(out, x, A, idx,
                    {cst_t(lin_t(i) - idx, cst_t::INEQUALITY)}, 0));
}

// Fig. 3:  i := 0; j := 0;
//          while (i < n) { t := A[i]; if (*) { C[j] := 1; j := j + 1; }
//                          i := i + 1; }
// At the exit C is {0} 1 {j}? T, so C[k] == 1 for every 0 <= k < j.
BOOST_AUTO_TEST_CASE(fig3_partial_initialization) {
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit");
  BB(cfg, entry);
  BB(cfg, head);
  BB(cfg, body);
  BB(cfg, then_);
  BB(cfg, else_);
  BB(cfg, inc);
  BB(cfg, exit);
  entry >> head;
  head >> body;
  body >> then_;
  body >> else_;
  then_ >> inc;
  else_ >> inc;
  inc >> head;
  head >> exit;
  entry.assign(i, 0);
  entry.assign(j, 0);
  body.assume(lin_t(i) <= lin_t(n) - 1);
  body.array_load(t, A, i, 4);
  then_.array_store(C, j, 1, 4);
  then_.add(j, j, 1);
  inc.add(i, i, 1);
  exit.assume(lin_t(i) >= lin_t(n));

  auto inv = analyze(*cfg);
  print("fig 3 head", inv.at("head"));
  print("fig 3 exit", inv.at("exit.post"));

  const dom_t &out = inv.at("exit.post");
  check_bounds(out.segmentation_of(C), {B({k(0)}), B({j})}, {false, true});
  BOOST_TEST(loads(out, x, C, idx,
                   {cst_t(-lin_t(idx), cst_t::INEQUALITY),
                    cst_t(lin_t(idx) - j + 1, cst_t::INEQUALITY)},
                   1));
}

BOOST_AUTO_TEST_SUITE_END()

// Keep stdout empty unless --print-invariants is given: the golden-output
// harness (tests/run_tests.sh) diffs the standard output of every test
// binary, so the whole Boost.Test log is routed to stderr. ctest passes
// --disable-warnings, which Boost.Test would reject, so it is translated
// into the crab flag it stands for.
int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int a = 0; a < argc; ++a) {
    std::string arg(argv[a]);
    if (arg == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    if (arg == "--print-invariants") {
      print_invariants = true;
      continue;
    }
    args.push_back(argv[a]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
