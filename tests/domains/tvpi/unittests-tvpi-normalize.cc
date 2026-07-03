/**
 * Tests for TVPIDBMNormalizeParams.
 *
 * Each case exercises TVPIDBMNormalizeParams (normalize=1) which auto-calls
 * normalize() after every transfer function, keeping the domain in closed form.
 */
#include "../../common.hpp"
#include "../../program_options.hpp"

#include <cassert>

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using z_sdbm_t = z_sdbm_domain_t;
using default_dom_t =
    crab::domains::tvpi_dbm_domain<z_sdbm_t,
                                   crab::domains::TVPIDBMDefaultParams>;
using normalize_dom_t =
    crab::domains::tvpi_dbm_domain<z_sdbm_t,
                                   crab::domains::TVPIDBMNormalizeParams>;

static unsigned idx = 1;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  auto &coeffs = crab_domain_params_man::get().coefficients();
  coeffs.insert(coeffs.end(), {2, 3, 4});

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);

  // ------------------------------------------------------------------
  // Case 1: basic TVPI fact — x = 4*i, i = n  =>  x = 4*n
  //   DbmClosure suffices; both params should prove this.
  // ------------------------------------------------------------------
  {
    crab::outs() << "\n---- case " << idx << " (basic TVPI, both prove) ----\n";

    default_dom_t dom_def;
    dom_def += (i >= z_number(0));
    dom_def += (i <= z_number(10));
    dom_def += (x == z_number(4) * i);
    dom_def.assign(n, i);
    bool r_def = dom_def.entails(x == z_number(4) * n);

    normalize_dom_t dom_norm;
    dom_norm += (i >= z_number(0));
    dom_norm += (i <= z_number(10));
    dom_norm += (x == z_number(4) * i);
    dom_norm.assign(n, i);
    bool r_norm = dom_norm.entails(x == z_number(4) * n);

    crab::outs() << "  default   entails(x == 4*n): "
                 << (r_def ? "true" : "false") << "\n";
    crab::outs() << "  normalize entails(x == 4*n): "
                 << (r_norm ? "true" : "false") << "\n";
    assert(r_def && "case1: default should prove x==4*n");
    assert(r_norm && "case1: normalize should prove x==4*n");
    idx++;
  }

  // ------------------------------------------------------------------
  // Case 2: auto-normalize vs manual.
  //   x = 4*i, y = 2*i, i >= 0  =>  y <= x  (2i <= 4i, trivially).
  //   - default without normalize():  may miss cross-coeff constraint
  //   - default + manual normalize(): should prove it
  //   - TVPIDBMNormalizeParams (auto): same result as manual normalize
  // ------------------------------------------------------------------
  {
    crab::outs() << "\n---- case " << idx
                 << " (auto-normalize vs manual) ----\n";

    // Default: no auto-normalize
    default_dom_t d_raw;
    d_raw += (i >= z_number(0));
    d_raw += (i <= z_number(10));
    d_raw += (x == z_number(4) * i);
    d_raw += (y == z_number(2) * i);
    bool r_raw = d_raw.entails(y <= x);

    // Default + explicit normalize
    default_dom_t d_manual = d_raw;
    d_manual.normalize();
    bool r_manual = d_manual.entails(y <= x);

    // Normalize params: auto-normalize after each +=
    normalize_dom_t n_auto;
    n_auto += (i >= z_number(0));
    n_auto += (i <= z_number(10));
    n_auto += (x == z_number(4) * i);
    n_auto += (y == z_number(2) * i);
    bool r_auto = n_auto.entails(y <= x);

    crab::outs() << "  default  (no normalize)     entails(y<=x): "
                 << (r_raw ? "true" : "false") << "\n";
    crab::outs() << "  default  (manual normalize) entails(y<=x): "
                 << (r_manual ? "true" : "false") << "\n";
    crab::outs() << "  normalize (auto)             entails(y<=x): "
                 << (r_auto ? "true" : "false") << "\n";

    assert(r_manual && "case2: default+manual_normalize must prove y<=x");
    assert(r_auto == r_manual &&
           "case2: auto-normalize must match manual normalize");
    idx++;
  }

  // ------------------------------------------------------------------
  // Case 3: prog3-style array bound.
  //   isz=4, x = 4*i, offset = x+4, tsz >= 4*len, i < len  =>  offset <= tsz
  // ------------------------------------------------------------------
  {
    crab::outs() << "\n---- case " << idx << " (array bound) ----\n";

    z_var isz(vfac["isz"], crab::INT_TYPE, 32);
    z_var len(vfac["len"], crab::INT_TYPE, 32);
    z_var tsz(vfac["tsz"], crab::INT_TYPE, 32);
    z_var offset(vfac["offset"], crab::INT_TYPE, 32);

    normalize_dom_t dom;
    dom += (isz == z_number(4));
    dom += (len >= z_number(1));
    dom += (len <= z_number(10));
    dom += (tsz >= z_number(4) * len);
    dom += (i >= z_number(0));
    dom += (i <= len - 1);
    dom += (x == z_number(4) * i);
    dom += (offset == x + z_number(4));

    bool r = dom.entails(offset <= tsz);
    crab::outs() << "  normalize entails(offset <= tsz): "
                 << (r ? "true" : "false") << "\n";
    assert(r && "case3: normalize should prove offset <= tsz");
    idx++;
  }

  // ------------------------------------------------------------------
  // Case 4: normalize() is idempotent.
  //   After auto-normalize, calling normalize() again gives the same domain.
  // ------------------------------------------------------------------
  {
    crab::outs() << "\n---- case " << idx << " (normalize idempotent) ----\n";

    normalize_dom_t dom1;
    dom1 += (i >= z_number(0));
    dom1 += (i <= z_number(6));
    dom1 += (j >= z_number(4));
    dom1 += (j <= z_number(10));
    dom1 += (x == z_number(3) * j);
    dom1 += (y == z_number(2) * i);

    normalize_dom_t dom2 = dom1;
    dom2.normalize();

    bool leq1 = (dom1 <= dom2);
    bool leq2 = (dom2 <= dom1);
    crab::outs() << "  dom1 <= dom1.normalize(): " << (leq1 ? "true" : "false")
                 << "\n";
    crab::outs() << "  dom1.normalize() <= dom1: " << (leq2 ? "true" : "false")
                 << "\n";
    assert(leq1 && leq2 && "case4: normalize must be idempotent");
    idx++;
  }

  // ------------------------------------------------------------------
  // Case 5: bottom detection.
  // ------------------------------------------------------------------
  {
    crab::outs() << "\n---- case " << idx << " (bottom detection) ----\n";

    normalize_dom_t dom;
    dom += (i >= z_number(5));
    dom += (i <= z_number(3));

    bool is_bot = dom.is_bottom();
    crab::outs() << "  after 5<=i<=3, is_bottom: "
                 << (is_bot ? "true" : "false") << "\n";
    assert(is_bot && "case5: domain must detect bottom");
    idx++;
  }

  return 0;
}
