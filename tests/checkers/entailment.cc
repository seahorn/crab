#include "../common.hpp"
#include "../program_options.hpp"
#include "crab/domains/abstract_domain_specialized_traits.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

template <typename Dom> void check_entailment(Dom inv, const z_lin_cst_t &cst) {
  bool r = inv.entails(cst);
  if (r) {
    crab::outs() << inv << " entails " << cst << "\n";
  } else {
    crab::outs() << inv << " does not entail " << cst << "\n";
  }
}

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  z_sdbm_domain_t inv1;
  inv1 += crab::cfg_impl::z_lin_cst_t(x >= 5);
  inv1 += crab::cfg_impl::z_lin_cst_t(x <= 10);
  crab::cfg_impl::z_lin_cst_t c1(x == 7);
  crab::cfg_impl::z_lin_cst_t c2(x == 11);
  check_entailment(inv1, c1);
  check_entailment(inv1, c2);

  z_sdbm_domain_t inv2;
  inv2 += crab::cfg_impl::z_lin_cst_t(x == 5);
  inv2 += crab::cfg_impl::z_lin_cst_t(y >= 0);
  inv2 += crab::cfg_impl::z_lin_cst_t(y <= 10);
  crab::cfg_impl::z_lin_cst_t c3(x == 5);
  crab::cfg_impl::z_lin_cst_t c4(x == 7);
  crab::cfg_impl::z_lin_cst_t c5(y != 42);
  check_entailment(inv2, c3);
  check_entailment(inv2, c4);
  check_entailment(inv2, c5);

  return 0;
}
