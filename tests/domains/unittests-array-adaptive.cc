#include "../program_options.hpp"
#include "../common.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main (int argc, char** argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
    return 0;
  }

  crab::outs() << "Unit tests for array adaptive domain\n\n";
  
  { // join
    variable_factory_t vfac;
    unsigned elem_sz = 4;
    z_var a(vfac["A"], crab::ARR_INT_TYPE);
    z_var b(vfac["B"], crab::ARR_INT_TYPE);
    z_var c(vfac["C"], crab::ARR_INT_TYPE);            
    z_aa_bool_int_t inv1,inv2,inv3;
    inv1.array_store(a, elem_sz, 4, 5, true);
    inv1.array_store(b, elem_sz, 0, 66, true);    
    inv2.array_store(a, elem_sz, 4, 10, true);
    inv2.array_store(c, elem_sz, 8, 32, true);                
    inv3 = inv1 | inv2;
    crab::outs() << inv1 << " | " << inv2 << "=" << inv3 << "\n";
    crab::outs() << inv1 << " |= " << inv2 << "=";
    inv1 |= inv2;
    crab::outs() << inv1 << "\n";
  }

  { // meet
    variable_factory_t vfac;
    unsigned elem_sz = 4;    
    z_var a(vfac["A"], crab::ARR_INT_TYPE);
    z_var b(vfac["B"], crab::ARR_INT_TYPE);
    z_var c(vfac["C"], crab::ARR_INT_TYPE);        
    z_aa_bool_int_t inv1,inv2,inv3;
    inv1.array_store(a, elem_sz, 4, 5, true);
    inv1.array_store(b, elem_sz, 0, 66, true);        
    inv2.array_store(a, elem_sz, 4, 10, true);
    inv2.array_store(c, elem_sz, 8, 32, true);            
    inv3 = inv1 & inv2;
    crab::outs() << inv1 << " & " << inv2 << "=" << inv3 << "\n";
  }

  { // meet
    variable_factory_t vfac;
    unsigned elem_sz = 4;    
    z_var a(vfac["A"], crab::ARR_INT_TYPE);
    z_var b(vfac["B"], crab::ARR_INT_TYPE);
    z_var c(vfac["C"], crab::ARR_INT_TYPE);        
    z_aa_bool_int_t inv1,inv2,inv3;
    inv1.array_store(a, elem_sz, 4, 5, true);
    inv1.array_store(b, elem_sz, 0, 66, true);        
    inv2.array_store(a, elem_sz, 4, 5, true);
    inv2.array_store(c, elem_sz, 8, 32, true);            
    inv3 = inv1 & inv2;
    crab::outs() << inv1 << " & " << inv2 << "=" << inv3 << "\n";
  }
  
  return 0;
}
