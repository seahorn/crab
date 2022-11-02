#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main(int argc, char **argv) {
#ifdef HAVE_ELINA
  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var v(vfac["v"], crab::INT_TYPE, 32);

  z_pk_elina_domain_t inv;
  crab::outs() << " ===== " << inv.domain_name() << " ====\n";
  inv += z_lin_cst_t(x <= 0);
  inv += z_lin_cst_t(y == 5);
  crab::outs() << "Before renaming " << inv << "\n";
  
  {
    z_pk_elina_domain_t tmp(inv);
    tmp.rename({x, y}, {w, z});
    crab::outs() << "After rename x with w and y with z=" << tmp << "\n";
  }
  {
    z_pk_elina_domain_t tmp(inv);
    tmp.rename({x}, {w});
    crab::outs() << "After rename x with w=" << tmp << "\n";
  }
  // {
  //   z_pk_elina_domain_t tmp(inv);
  //   tmp.rename({x}, {y});
  //   crab::outs() << "After rename x with y=" << tmp << "\n";
  //   // This should raise an error:
  //   //  CRAB ERROR: rename expects y to be unconstrained
  // }
#endif

  return 0;
}
