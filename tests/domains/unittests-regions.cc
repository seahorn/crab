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
  
  { // join
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);    
    z_var rgn1(vfac["region_0"], crab::REG_INT_TYPE, 32);
    z_var_or_num_t n34_32(z_number(34), crab::variable_type(crab::INT_TYPE, 32));
    
    z_rgn_int_t inv1, inv2, inv3;
    inv1 += (x >= z_number(5));
    inv1 += (y >= z_number(10));
    inv1.region_init(rgn1);
    
    inv2 += (x >= z_number(5));
    inv2.region_init(rgn1);
    inv2.ref_make(ref, rgn1);
    inv2.ref_store(ref, rgn1, n34_32);


    crab::outs() << "Join of\n\t" << inv1  << "\nand\n\t" << inv2 << " =\n\t";
    inv3 = inv1 | inv2;
    crab::outs() << inv3 << "\n";
  }
  
  { // join
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);    
    z_var rgn1(vfac["region_0"], crab::REG_INT_TYPE, 32);

    z_var_or_num_t n34_32(z_number(34), crab::variable_type(crab::INT_TYPE, 32));
    z_rgn_int_t inv1, inv2;
    inv1 += (x >= z_number(5));
    inv1 += (y >= z_number(10));
    inv1.region_init(rgn1);
    
    inv2 += (x >= z_number(5));
    inv2.region_init(rgn1);
    inv2.ref_make(ref, rgn1);
    inv2.ref_store(ref, rgn1, n34_32);


    crab::outs() << "Join of\n\t" << inv1  << "\nand\n\t" << inv2 << " =\n\t";
    inv1 |= inv2;
    crab::outs() << inv1 << "\n";
  }

  return 0;
}
