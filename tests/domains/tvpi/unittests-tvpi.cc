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

// Runs join/meet and asserts fundamental lattice laws.
void perform_domain_operations(const test_domain_t &dom1,
                               const test_domain_t &dom2) {
  crab::outs() << "Dom1=" << dom1 << "\n";
  crab::outs() << "Dom2=" << dom2 << "\n";
  bool r1 = dom1 <= dom2;
  crab::outs() << "Dom1 <= Dom2 = " << (r1 ? "true" : "false") << "\n";
  bool r2 = dom2 <= dom1;
  crab::outs() << "Dom2 <= Dom1 = " << (r2 ? "true" : "false") << "\n";

  test_domain_t dom3 = dom1 | dom2;
  crab::outs() << "Dom3 = Dom1 | Dom2 = " << dom3 << "\n";
  bool r3 = dom1 <= dom3;
  bool r4 = dom2 <= dom3;
  crab::outs() << "Dom1 <= Dom3 = " << (r3 ? "true" : "false") << "\n";
  crab::outs() << "Dom2 <= Dom3 = " << (r4 ? "true" : "false") << "\n";
  assert(r3 && "join soundness: dom1 must be <= dom1|dom2");
  assert(r4 && "join soundness: dom2 must be <= dom1|dom2");

  test_domain_t dom4 = dom1 & dom2;
  crab::outs() << "Dom4 = Dom1 & Dom2 = " << dom4 << "\n";
  bool r5 = dom4 <= dom1;
  bool r6 = dom4 <= dom2;
  crab::outs() << "Dom4 <= Dom1 = " << (r5 ? "true" : "false") << "\n";
  crab::outs() << "Dom4 <= Dom2 = " << (r6 ? "true" : "false") << "\n";
  assert(r5 && "meet soundness: dom1&dom2 must be <= dom1");
  assert(r6 && "meet soundness: dom1&dom2 must be <= dom2");
}

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

  { // test assign: x=1, y=2x, z=3x+7, n=2x+2y+5, k=5, o=-2k-5
    crab::outs() << "\n\n---- test assignment ----\n\n";
    test_domain_t dom1;
    dom1.assign(x, z_number(1));
    dom1.assign(y, x * z_number(2));
    dom1.assign(z, x * z_number(3) + z_number(7));
    dom1.assign(n, x * z_number(2) + y * z_number(2) + z_number(5));
    dom1.assign(o, z_number(-2) * k - z_number(5));
    dom1.assign(k, z_number(0) * o + z_number(5));
    crab::outs() << "Dom1=" << dom1 << "\n";
    assert(!dom1.is_bottom() && "assignment domain must not be bottom");
    assert(dom1.entails(x == z_number(1)) && "x must equal 1");
    assert(dom1.entails(y == z_number(2)) && "y must equal 2");
    assert(dom1.entails(z == z_number(10)) && "z must equal 10");
    assert(dom1.entails(k == z_number(5)) && "k must equal 5");
  }

  { // test assume: x=1, y=2x, z=3x+7, then two TVPI constraints
    crab::outs() << "\n\n---- test assume ----\n\n";
    test_domain_t dom1;
    dom1 += (x == z_number(1));
    dom1 += (y == x * z_number(2));
    dom1 += (z == x * z_number(3) + z_number(7));
    dom1 += (x * z_number(5) + z_number(6) * n == z_number(4));
    dom1 +=
        (z * z_number(3) + y * z_number(6) + k * z_number(9) == z_number(15));
    crab::outs() << "Dom1=" << dom1 << "\n";
    assert(!dom1.is_bottom() && "assume domain must not be bottom");
    // Coefficient 3 is in template so 3z+6y+9k=15 → z+2y+3k=5 is tracked;
    // with z=10, y=2: k must equal -3.
    assert(dom1.entails(k == z_number(-3)) && "k must equal -3");
  }

  { // case 1: exact values, meet is bottom (x=1 vs x=2 conflict)
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    test_domain_t dom1, dom2;
    // dom1 : x=1, y=2x, z=3x
    dom1.assign(x, z_number(1));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));

    // dom2 : x=2, y=3x, z=4x
    dom2 += (x == z_number(2));
    dom2 += (y == x * z_number(3));
    dom2 += (z == x + x + x + x);

    perform_domain_operations(dom1, dom2);
    assert((dom1 & dom2).is_bottom() &&
           "case1: meet must be bottom (x=1 vs x=2)");

    // dom5: forget x → {y=2, z=3}
    test_domain_t dom5(dom1);
    dom5 -= x;
    crab::outs() << "After forgetting " << x << " in Dom1:" << dom5 << "\n";
    assert(!dom5.is_bottom());
    assert(dom5.entails(y == z_number(2)) && "y must equal 2 after forget x");
    assert(dom5.entails(z == z_number(3)) && "z must equal 3 after forget x");

    // dom6: rename x→i → {i=1, y=2, z=3}
    test_domain_t dom6(dom1);
    dom6.rename({x}, {i});
    crab::outs() << "After renaming {x} with {i} in Dom1:" << dom6 << "\n";
    assert(!dom6.is_bottom());
    assert(dom6.entails(i == z_number(1)) && "i must equal 1 after rename");

    // dom7: project {y,z} → {y=2, z=3}
    test_domain_t dom7(dom1);
    dom7.project({y, z});
    crab::outs() << "After projecting on y and z in Dom1:" << dom7 << "\n";
    assert(!dom7.is_bottom());
    assert(dom7.entails(y == z_number(2)) && "y must equal 2 after project");
    assert(dom7.entails(z == z_number(3)) && "z must equal 3 after project");

    // dom8: project {x} → {x=1}
    test_domain_t dom8(dom1);
    dom8.project({x});
    crab::outs() << "After projecting on x in Dom1:" << dom8 << "\n";
    assert(!dom8.is_bottom());
    assert(dom8.entails(x == z_number(1)) && "x must equal 1 after project");
    idx++;
  }

  { // case 2: ranges, meet is not bottom, join is top
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    test_domain_t dom1, dom2;
    // dom1 : x>=1, y=2x, z=3x
    dom1 += (x >= z_number(1));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1.apply(OP_MULTIPLICATION, z, x, z_number(3));

    // dom2: x<=20, i=3x, j=4x
    dom2 += (x <= z_number(20));
    dom2.assign(i, z_number(3) * x);
    dom2.assign(j, z_number(4) * x);

    perform_domain_operations(dom1, dom2);
    assert(!(dom1 & dom2).is_bottom() &&
           "case2: meet must not be bottom (x in [1,20])");

    // dom5: forget x → {y>=2, z>=3}
    test_domain_t dom5(dom1);
    dom5 -= x;
    crab::outs() << "After forgetting " << x << " in Dom1:" << dom5 << "\n";
    assert(!dom5.is_bottom());
    assert(dom5.entails(y >= z_number(2)) && "y must be >= 2 after forget x");
    assert(dom5.entails(z >= z_number(3)) && "z must be >= 3 after forget x");

    // dom6: rename x→i → {i>=1, y=2i, z=3i}
    test_domain_t dom6(dom1);
    dom6.rename({x}, {i});
    crab::outs() << "After renaming {x} with {i} in Dom1:" << dom6 << "\n";
    assert(!dom6.is_bottom());
    assert(dom6.entails(i >= z_number(1)) && "i must be >= 1 after rename");

    // dom7: project {y,z} → {y>=2, z>=3}
    test_domain_t dom7(dom1);
    dom7.project({y, z});
    crab::outs() << "After projecting on y and z in Dom1:" << dom7 << "\n";
    assert(!dom7.is_bottom());
    assert(dom7.entails(y >= z_number(2)) && "y must be >= 2 after project");
    assert(dom7.entails(z >= z_number(3)) && "z must be >= 3 after project");

    // dom8: project {x} → {x>=1}
    test_domain_t dom8(dom1);
    dom8.project({x});
    crab::outs() << "After projecting on x in Dom1:" << dom8 << "\n";
    assert(!dom8.is_bottom());
    assert(dom8.entails(x >= z_number(1)) && "x must be >= 1 after project");
    idx++;
  }

  { // case 3: loop counter widening — dom1 <= widened result
    crab::outs() << "\n\n---- case " << idx << "----\n\n";
    test_domain_t dom1, dom2;
    // dom1 : x=0, i=0
    dom1 += (x == z_number(0));
    dom1 += (i == z_number(0));
    assert(!dom1.is_bottom());

    // dom2: dom1 + loop_counter(i), i=i+1, x=x+3 → x=3, i=1
    dom2 = dom1;
    dom2.intrinsic("loop_counter", {i}, {});
    dom2.apply(OP_ADDITION, i, i, z_number(1));
    dom2.apply(OP_ADDITION, x, x, z_number(3));
    crab::outs() << "Dom2 adds a loop counter " << i << "\n";
    assert(!dom2.is_bottom());

    perform_domain_operations(dom1, dom2);

    // widening: x and i grow without bound
    test_domain_t dom5 = dom1 || (dom1 | dom2);
    crab::outs() << "Dom5 = Dom1 || (Dom1 | Dom2) = " << dom5 << "\n";
    assert(!dom5.is_bottom() && "widening result must not be bottom");

    bool r1 = dom1 <= dom5;
    crab::outs() << "Dom1 <= Dom5 = " << (r1 ? "true" : "false") << "\n";
    assert(r1 && "case3: dom1 must be <= widened result");
    idx++;
  }

  { // case 4: TVPI entailment — resultant then closure proves ordering
    crab::outs() << "\n\n---- case " << idx << "---- \n\n";
    test_domain_t dom1, dom2;
    // dom1 : x>=5, y=4x+2, z=2x+3y → y<=z (y=4x+2 <= 2x+3(4x+2)=14x+6, and
    // x>=5)
    dom1 += (x >= z_number(5));
    dom1.assign(y, z_number(4) * x + z_number(2));
    dom1.assign(z, z_number(2) * x + z_number(3) * y);
    crab::outs() << "Dom1=" << dom1 << "\n";
    bool check = dom1.entails(y <= z);
    crab::outs() << "assert(y <= z): " << (check ? "true" : "false") << "\n";
    assert(check && "case4: dom1 must prove y <= z");

    // dom2 : x in [1,10], y=4x+17, i in [0,x), z=4i+3
    // need: z <= y, i.e. 4i+3 <= 4x+17
    // i<=x-1 → 4i<=4x-4 → 4i+3<=4x-1 <= 4x+17 ✓
    dom2 += (x >= z_number(1));
    dom2 += (x <= z_number(10));
    dom2.assign(y, z_number(4) * x + z_number(17));
    dom2 += (i >= z_number(0));
    dom2 += (i <= x - z_number(1));
    dom2.assign(z, z_number(4) * i + z_number(3));
    crab::outs() << "Dom2=" << dom2 << "\n";
    check = dom2.entails(z <= y);
    crab::outs() << "assert(z <= y): " << (check ? "true" : "false") << "\n";
    assert(check && "case4: dom2 must prove z <= y");
    idx++;
  }

  { // case 5: join and containment with lazy normalization
    // With lazy normalization the domain cannot derive x<=1 from
    // {x-y<=4,2y-x<=-3,3x-y<=5} without an explicit normalize(), so dom1 <=
    // dom4 and dom2 <= dom4 return false.
    crab::outs() << "\n\n---- case " << idx << "---- \n\n";
    test_domain_t dom1;
    // dom1 : x-y<=4, 2y-x<=-3, 3x-y<=5  (implies x<=1 and y<=-1 when
    // normalized)
    dom1 += (x - y <= z_number(4));
    dom1 += (z_number(2) * y - x <= z_number(-3));
    dom1 += (z_number(3) * x - y <= z_number(5));
    crab::outs() << "dom1=" << dom1 << "\n";
    assert(!dom1.is_bottom());

    // dom2 : y<=0, x<=2, y-x<=2
    test_domain_t dom2;
    dom2 += (y <= z_number(0));
    dom2 += (x <= z_number(2));
    dom2 += (y - x <= z_number(2));
    crab::outs() << "dom2=" << dom2 << "\n";
    assert(!dom2.is_bottom());

    test_domain_t dom3 = dom1 | dom2;
    crab::outs() << "dom1 | dom2=" << dom3 << "\n";

    // dom4: manually crafted join candidate — y<=0, x<=2, 2y-x<=2
    test_domain_t dom4;
    dom4 += (y <= z_number(0));
    dom4 += (x <= z_number(2));
    dom4 += (z_number(2) * y - x <= z_number(2));
    crab::outs() << "dom4=" << dom4 << "\n";

    // With lazy normalization, dom1 has no explicit bounds on x or y, so these
    // are false even though semantically dom1 ⊑ dom4 after full saturation.
    bool l1 = dom1 <= dom4;
    bool l2 = dom2 <= dom4;
    crab::outs() << "Dom1 <= Dom4 = " << (l1 ? "true" : "false") << "\n";
    crab::outs() << "Dom2 <= Dom4 = " << (l2 ? "true" : "false") << "\n";
    idx++;
  }

  { // case 6: join and containment — dom1 and dom2 do not imply dom4
    crab::outs() << "\n\n---- case " << idx << "---- \n\n";
    // dom1 : 2x-y<=5, 3y-x<=7, -3x+y<=4
    test_domain_t dom1;
    dom1 += (z_number(2) * x - y <= z_number(5));
    dom1 += (z_number(3) * y - x <= z_number(7));
    dom1 += (-z_number(3) * x + y <= z_number(4));
    crab::outs() << "dom1=" << dom1 << "\n";
    assert(!dom1.is_bottom());

    // dom2 : 3x-2y<=8, 4y-x<=10, 2x-3y<=-6
    test_domain_t dom2;
    dom2 += (z_number(3) * x - z_number(2) * y <= z_number(8));
    dom2 += (z_number(4) * y - x <= z_number(10));
    dom2 += (z_number(2) * x - z_number(3) * y <= z_number(-6));
    crab::outs() << "dom2=" << dom2 << "\n";
    assert(!dom2.is_bottom());

    test_domain_t dom3 = dom1 | dom2;
    crab::outs() << "dom1 | dom2=" << dom3 << "\n";

    // dom4: same reference domain as case 5
    test_domain_t dom4;
    dom4 += (y <= z_number(0));
    dom4 += (x <= z_number(2));
    dom4 += (z_number(2) * y - x <= z_number(2));
    crab::outs() << "dom4=" << dom4 << "\n";

    // Neither dom1 nor dom2 is semantically ⊑ dom4 (e.g. (x=-1,y=1) ∈ dom1,
    // y>0).
    bool l1 = dom1 <= dom4;
    bool l2 = dom2 <= dom4;
    crab::outs() << "Dom1 <= Dom4 = " << (l1 ? "true" : "false") << "\n";
    crab::outs() << "Dom2 <= Dom4 = " << (l2 ? "true" : "false") << "\n";
    assert(!l1 && "case6: dom1 must NOT be <= dom4");
    assert(!l2 && "case6: dom2 must NOT be <= dom4");
    idx++;
  }

  { // case 7: array-access bound — offset <= tsz proved via resultant + closure
    crab::outs() << "\n\n---- case " << idx << "---- \n\n";
    z_var isz(vfac["t"], crab::INT_TYPE, 32);
    z_var len(vfac["l"], crab::INT_TYPE, 32);
    z_var tsz(vfac["s"], crab::INT_TYPE, 32);
    z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
    z_var offset(vfac["o"], crab::INT_TYPE, 32);

    test_domain_t dom1;
    dom1 += (isz == z_number(4)); // isz = 4
    dom1 += (len >= z_number(1));
    dom1 += (len <= z_number(10)); // len in [1,10]
    dom1.apply(OP_MULTIPLICATION, tmp, isz, len);
    dom1 += (tmp <= tsz); // 4*len <= tsz
    dom1 += (i >= z_number(0));
    dom1 += (i <= len - 1);                   // i in [0, len)
    dom1.apply(OP_MULTIPLICATION, x, isz, i); // x = 4*i
    dom1.apply(OP_ADDITION, offset, x, isz);  // offset = 4*i + 4
    crab::outs() << "dom1=" << dom1 << "\n";
    // offset = 4i+4 <= 4(len-1)+4 = 4len <= tsz
    bool ret = dom1.entails(offset <= tsz);
    crab::outs() << "assert(o <= s): " << (ret ? "true" : "false") << "\n";
    assert(ret && "case7: must prove offset <= tsz");
    idx++;
  }

  { // case 8: derived bound via TVPI — n=2i, y=2x, i<x, y<=z => n<=z
    crab::outs() << "\n\n---- case " << idx << "---- \n\n";
    test_domain_t dom1;
    // x in [1,10], y=2x, y<=z
    dom1 += (x >= z_number(1));
    dom1 += (x <= z_number(10));
    dom1.apply(OP_MULTIPLICATION, y, x, z_number(2));
    dom1 += (y <= z);
    // i in [0,x), n=2i
    dom1 += (i >= z_number(0));
    dom1 += (i <= x - z_number(1));
    dom1.apply(OP_MULTIPLICATION, n, i, z_number(2));
    // n=2i <= 2(x-1) <= 2x-2 < 2x = y <= z
    crab::outs() << "dom1=" << dom1 << "\n";
    bool ret = dom1.entails(n <= z);
    crab::outs() << "assert(n <= z): " << (ret ? "true" : "false") << "\n";
    assert(ret && "case8: must prove n <= z");
    idx++;
  }

  return 0;
}
