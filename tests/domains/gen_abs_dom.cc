#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);

  z_sdbm_domain_t sdbm_top;
  z_abs_domain_t dom1(sdbm_top);
  dom1 += (x >= y);
  dom1 += (z >= 10);
  // This also allowed because the constructor is not explicit
  z_abs_domain_t dom2 = sdbm_top;
  dom2 += (x >= y);
  dom2 += (z >= 5);

  auto dom3 = dom1 | dom2;
  crab::outs() << "Join(" << dom1 << "," << dom2 << ")=" << dom3 << "\n";

  auto dom4 = dom1 & dom2;
  crab::outs() << "Meet(" << dom1 << "," << dom2 << ")=" << dom4 << "\n";

  auto dom5 = dom1 || dom2;
  crab::outs() << "Widen(" << dom1 << "," << dom2 << ")=" << dom5 << "\n";

  z_abs_domain_t dom6(dom1);
  dom6.apply(OP_ADDITION, w, x, 5);
  crab::outs() << w << ":=" << x << "+ 5 in " << dom1 << "=" << dom6 << "\n";

  //// This should not compile:
  // z_abs_domain_t dom4;
  // crab::outs() << dom4;

  return 0;
}
