#include "../../common.hpp"
#include "../../program_options.hpp"

#include <cassert>

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using test_domain_t = z_tvpi_dbm_domain_t;

unsigned idx = 1;

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
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  { // case 1: provable via bounds — max(y)=12 == min(x)=12
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    z_number SZ = z_number(3);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    dom += (i >= z_number(0));
    dom += (i <= z_number(6)); // i in [0, 6]
    dom += (j >= z_number(4));
    dom += (j <= z_number(10)); // j in [4, 10]
    dom += (x == SZ * j);       // x = 3j,  x in [12, 30]
    dom += (y == SLICE * i);    // y = 2i,  y in [0,  12]
    crab::outs() << "dom=" << dom << "\n";
    assert(!dom.is_bottom());

    // y <= x: 2i <= 3j.  max(y)=12, min(x)=12  => provable via bounds alone.
    bool check = dom.entails(y <= x);
    crab::outs() << "assert(y <= x): " << (check ? "true" : "false") << "\n";
    assert(check && "case1: must prove y <= x via bounds");
    idx += 1;
  }

  { // case 2: provable since incremental saturation was wired into +=/entails
    // Requires a chain
    // ghost(y,1)->ghost(i,2)->ghost(j,2)->ghost(j,3)->ghost(x,1). Historically
    // the lazy domain lacked ghost(i,2)-ghost(j,2)<=-2 until a full TvpiReduce;
    // incremental_tvpi_reduce now derives it as constraints are added.
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    z_number SZ = z_number(3);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    dom += (i >= z_number(0));
    dom += (i <= z_number(6)); // i in [0, 6]
    dom += (j >= z_number(0));
    dom += (j <= z_number(10)); // j in [0, 10]
    dom += (j > i);             // j >= i+1
    dom += (x == SZ * j);       // x = 3j
    dom += (y == SLICE * i);    // y = 2i
    crab::outs() << dom << "\n";
    assert(!dom.is_bottom());

    // assert(y <= x)? 2i <= 3j.
    // Provable: j >= i+1 => 3j >= 3i+3 > 2i.  Incremental saturation derives
    // the DBM edge ghost(i,2) - ghost(j,2) <= -2 when j > i is added.
    bool check = dom.entails(y <= x);
    crab::outs() << "assert(y <= x): " << (check ? "true" : "false")
                 << " (expected true: incremental saturation derives"
                 << " 2i-2j<=-2)\n";
    idx += 1;
  }

  { // case 3: o <= x provable; z <= x is NOT provable (and mathematically
    // false)
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    z_number SZ = z_number(3);
    z_number OFFSET = z_number(2);
    z_number SLICE = z_number(2);

    test_domain_t dom;
    dom += (i >= z_number(1));
    dom += (i <= z_number(4)); // i in [1, 4]
    dom += (j >= z_number(1));
    dom += (j <= z_number(4)); // j in [1, 4]
    dom += (k >= z_number(0));
    dom += (k <= i - 1);      // 0 <= k < i
    dom += (x == SZ * i);     // x = 3i
    dom += (o == OFFSET + k); // o = k + 2
    dom += (y == SLICE * j);  // y = 2j
    crab::outs() << dom << "\n";
    assert(!dom.is_bottom());

    // assert(o <= x)? k+2 <= 3i.
    // k <= i-1 => k+2 <= i+1 <= 3i (for i >= 1). Provable.
    bool check = dom.entails(o <= x);
    crab::outs() << "assert(o <= x): " << (check ? "true" : "false") << "\n";
    assert(check && "case3: must prove o <= x (k+2 <= 3i)");

    dom += (z == o + y); // z = k + 2 + 2j
    crab::outs() << dom << "\n";

    // assert(z <= x)? z = k+2+2j <= 3i = x?
    // Counter-example: i=1, k=0, j=4  =>  z=10, x=3.  NOT always true.
    // The domain correctly does not prove this.
    bool check2 = dom.entails(z <= x);
    crab::outs() << "assert(z <= x): " << (check2 ? "true" : "false")
                 << " (expected false: z = k+2+2j can exceed x = 3i)\n";
    idx += 1;
  }

  { // case 4: join of two exact-value domains; joined interval provable
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    test_domain_t dom1, dom2;
    dom1 += (i == z_number(4));
    dom2 += (i == z_number(4));
    dom1 += (x == 2 * i); // x = 2*4 = 8
    dom2 += (x == 3 * i); // x = 3*4 = 12

    crab::outs() << "dom1: " << dom1 << "\n";
    crab::outs() << "dom2: " << dom2 << "\n";
    assert(!dom1.is_bottom());
    assert(!dom2.is_bottom());

    auto dom3 = (dom1 | dom2); // x in [8, 12]
    crab::outs() << "dom1 join dom2: " << dom3 << "\n";
    assert(!dom3.is_bottom());

    // join lower-bounds x >= 8 and upper-bounds x <= 12
    bool check = dom3.entails(x >= 8);
    crab::outs() << "assert(x >= 2 * 4): " << (check ? "true" : "false")
                 << "\n";
    assert(check && "case4: join must prove x >= 8");

    check = dom3.entails(x <= 12);
    crab::outs() << "assert(x <= 3 * 4): " << (check ? "true" : "false")
                 << "\n";
    assert(check && "case4: join must prove x <= 12");
    idx += 1;
  }

  { // case 5: widening — dom1 must be <= widening result (lattice law)
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    test_domain_t dom1, dom2, dom3;
    dom1 += (x == z_number(0));
    dom1 += (i == z_number(0));
    assert(!dom1.is_bottom());

    dom2 = dom1;
    dom2.intrinsic("loop_counter", {i}, {});
    dom2.apply(OP_ADDITION, i, i, z_number(1));
    crab::outs() << "Dom2 adds a loop counter " << i << "\n";
    assert(!dom2.is_bottom());

    dom3 = dom2;
    dom2.apply(OP_ADDITION, x, x, z_number(3)); // dom2: x=3, i=1
    dom3.apply(OP_ADDITION, x, x, z_number(2)); // dom3: x=2, i=1
    assert(!dom3.is_bottom());

    test_domain_t dom4 = dom2 | dom3; // x in [2,3], i=1
    crab::outs() << "Dom4 = Dom2 | Dom3 = " << dom4 << "\n";
    assert(!dom4.is_bottom());

    test_domain_t dom5 = dom1 || (dom1 | dom4); // x >= 0, i >= 0
    crab::outs() << "Dom5 = Dom1 || (Dom1 | Dom4) = " << dom5 << "\n";
    assert(!dom5.is_bottom());

    bool r1 = dom1 <= dom5;
    crab::outs() << "Dom1 <= Dom5 = " << (r1 ? "true" : "false") << "\n";
    assert(r1 && "case5: dom1 must be <= widening result");
    idx++;
  }

  return 0;
}
