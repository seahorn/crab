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

  { // linear constraints normalization
    using linear_constraint_system_t =
        ikos::linear_constraint_system<ikos::z_number, varname_t>;

    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var w(vfac["w"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);

    linear_constraint_system_t csts;
    csts += z_lin_cst_t(x <= 8);
    csts += z_lin_cst_t(y >= 5);
    csts += z_lin_cst_t(y <= 6);
    csts += z_lin_cst_t(x >= 8);
    csts += z_lin_cst_t(v >= 9888);
    csts += z_lin_cst_t(z == 10);
    csts += z_lin_cst_t(w <= 0);
    csts += z_lin_cst_t(v <= 9888);
    crab::outs() << "Before normalize constraints " << csts << "\n";
    crab::outs() << "After normalize constraints " << csts.normalize() << "\n";
  }

  return 0;
}
