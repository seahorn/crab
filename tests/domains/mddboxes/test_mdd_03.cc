#include "crab/config.h"

#ifndef HAVE_MDD
int main(int argc, char**argv) {
  return 0;
}
#else
#include "../../common.hpp"
using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

void test_abstract_operations() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);  
  
  {
    z_mdd_boxes_domain_t inv1;
    inv1 += (z != 0);
    inv1.assign(x, 0); // <------ this is the problem
    z_mdd_boxes_domain_t inv2(inv1);
    z_mdd_boxes_domain_t inv3(inv1);
    inv2.assign(y, 5);
    inv3.assign(y, 0);    
    z_mdd_boxes_domain_t inv4 = inv2 | inv3;
    inv4.assign(x, y);
    crab::outs() << inv4 << "\n";
  }

}

int main (int argc, char** argv ) {
  test_abstract_operations();
  return 0;
}
#endif
  
