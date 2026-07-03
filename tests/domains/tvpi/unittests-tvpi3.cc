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
  coeffs.insert(coeffs.end(), {10, 255});

  variable_factory_t vfac;
  z_var c(vfac["c"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var o(vfac["o"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  { // tvpi paper example for C string
    // We convert into a simple program with no string but need string property
    /*
        char s[32] = "the string";
        int i = 0;
        while (true) {
            c = s[i]; <------ output
            if (c==0) break;
            i = i + 1;
        };
     */
    // dom1: i \in [0, 9], c \in [1, 255]
    test_domain_t dom1;
    dom1 += (i >= z_number(0));
    dom1 += (i <= z_number(9)); // i \in [0, 10)
    dom1 += (c <= z_number(255));
    dom1 += (c >= z_number(1)); // c \in [1, 255]
    assert(!dom1.is_bottom());

    // dom2: i = 10, c = 0
    test_domain_t dom2;
    dom2 += (i == z_number(10));
    dom2 += (c == z_number(0));
    assert(!dom2.is_bottom());

    // dom3: i > 10, c \in [0, 255]
    test_domain_t dom3;
    dom3 += (i > z_number(10));
    dom3 += (c <= z_number(255));
    dom3 += (c >= z_number(0));
    assert(!dom3.is_bottom());

    // input: i \in [0, 10]
    test_domain_t input;
    input += (i <= z_number(10));
    input += (i >= z_number(0));

    // input & dom3 = bottom: i<=10 AND i>10 is unsatisfiable (integers)
    assert((input & dom3).is_bottom() &&
           "input & dom3 must be bottom: i<=10 AND i>10");

    test_domain_t output = input & dom1;
    output |= input & dom2;
    output |= input & dom3; // adds bottom, no effect on join

    crab::outs() << "output: " << output << "\n";

    // TVPI domain captures (from the paper):
    // i \in [0, 10], c \in [0, 255]
    // 255i + c <= 2550
    // -i - 10c <= -10

    // for us, we can only compute the range if using zones
    // besides, coefficients cannot be extended unless we know extreme points
    // for new convex hull.
    // dom1: i \in [0, 9],       c \in [1, 255]
    // dom2: i = 10, c = 0
    // dom1 join dom2:
    //  i \in [0, 10], c \in [0, 255]
    //
    // The domain tracks ghost(v,k) for every coefficient k of the fixed
    // template set, for every original variable v it sees.
    // So i<=10 constrains ghost(i,255)<=2550 and c>=0 constrains ghost(c,1)>=0.
    // DbmClosure through v0 then derives: ghost(i,255) - ghost(c,1) <= 2550,
    // i.e. 255i - c <= 2550.  Similarly i - 10c <= 10.
    //
    // 255i + c <= 2550 (TVPI paper form): NOT provable — two positive terms,
    // not a difference constraint.  The difference form 255i - c <= 2550 IS
    // provable (see above).

    assert(!output.is_bottom());
    assert(output.entails(i >= z_number(0)));
    assert(output.entails(i <= z_number(10)));
    assert(output.entails(c >= z_number(0)));
    assert(output.entails(c <= z_number(255)));
    // TVPI difference constraints provable after join via DbmClosure on ghosts:
    //   255*10 - 0 = 2550  (dom2 extreme point drives the upper bound)
    assert(output.entails(z_number(255) * i - c <= z_number(2550)) &&
           "255i - c <= 2550 via ghost(i,255) - ghost(c,1)");
    //   10 - 10*0 = 10  (dom2: i=10, c=0)
    assert(output.entails(i - z_number(10) * c <= z_number(10)) &&
           "i - 10c <= 10 via ghost(i,1) - ghost(c,10)");
  }
  return 0;
}
