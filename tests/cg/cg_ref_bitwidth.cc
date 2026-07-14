/** Regression test for a variable_type hash/equality inconsistency.
 *
 * variable_type::operator== ignores m_bitwidth for every kind except
 * INT_TYPE and REG_INT_TYPE (see crab/types/variable.hpp). Its hash(),
 * however, used to fold m_bitwidth in unconditionally, breaking the
 * invariant a == b => hash(a) == hash(b) for references.
 *
 * This is invisible until callsite_or_fdecl hashing depends on the
 * argument/return variable_type hashes (as it does since the fix to
 * callsite_or_fdecl::compute_hash). A callsite whose reference argument
 * was built with one bitwidth (e.g. clam's mkRefVar uses 32) and a callee
 * declaration whose reference formal was built with another (e.g. the
 * default 0) are equal under operator== but hash into different buckets,
 * so the call graph fails to resolve the callee and aborts the analysis
 * with "Function not found for callsite".
 **/

#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/cg/cg.hpp>

#include <memory>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::cg;
using namespace crab::domain_impl;

using call_graph_t = call_graph<z_cfg_ref_t>;

// callee(ref p) -> (ref r): the reference formals are built with the
// default bitwidth (0), matching crab's own type helpers.
static std::unique_ptr<z_cfg_t> callee(variable_factory_t &vfac) {
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var r(vfac["r"], crab::REF_TYPE);
  function_decl<z_number, varname_t> decl("callee", {p}, {r});
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);
  BB(cfg, entry);
  BB(cfg, exit);
  entry >> exit;
  return cfg;
}

// main() calls callee passing reference actuals built with bitwidth 32,
// as clam's crabLitFactory::mkRefVar() does. REF_TYPE(0) and REF_TYPE(32)
// are equal under variable_type::operator== but hashed differently before
// the fix.
static std::unique_ptr<z_cfg_t> mk_main(variable_factory_t &vfac) {
  z_var a(vfac["a"], crab::REF_TYPE, 32);
  z_var b(vfac["b"], crab::REF_TYPE, 32);
  std::vector<z_var> inputs, outputs;
  function_decl<z_number, varname_t> decl("main", inputs, outputs);
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit", decl);
  BB(cfg, entry);
  BB(cfg, exit);
  entry >> exit;
  // callsite(func, lhs/outputs, args/inputs)
  entry.callsite("callee", {b}, {a});
  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool /*stats_enabled*/) -> int {
    variable_factory_t vfac;

    // (1) Directly exercise the invariant that regressed: two reference
    //     types that compare equal must hash equally.
    z_var ref0(vfac["ref0"], crab::REF_TYPE);       // bitwidth 0
    z_var ref32(vfac["ref32"], crab::REF_TYPE, 32); // bitwidth 32
    const auto t0 = ref0.get_type();
    const auto t32 = ref32.get_type();
    if (!(t0 == t32)) {
      crab::outs() << "FAIL: reference types must be equal regardless of "
                      "bitwidth\n";
      return 1;
    }
    if (t0.hash() != t32.hash()) {
      crab::outs() << "FAIL: variable_type::hash inconsistent with operator== "
                      "for REF_TYPE\n";
      return 1;
    }

    // (2) End-to-end: the call graph must resolve a callsite whose reference
    //     argument bitwidth differs from the callee declaration's. Before the
    //     fix, type_check() aborts with "Function not found for callsite".
    auto f_callee = callee(vfac);
    auto f_main = mk_main(vfac);
    std::vector<z_cfg_ref_t> cfgs({*f_callee, *f_main});
    call_graph_t cg(cfgs);
    cg.type_check();

    crab::outs() << "OK: callee resolved despite reference bitwidth mismatch\n";
    return 0;
  });
}
