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
    z_mdd_boxes_domain_t inv = z_mdd_boxes_domain_t::top();    
    inv.set(x, z_interval_t(2, 5));
    inv.set(y, z_interval_t(5, 10));    
    z_lin_t e(x + y + 3);
    inv.assign(z, e);
    crab::outs() << "z := x + y + 3 ---> " << inv << "\n";
  }

}

int main (int argc, char** argv ) {
  //variable_factory_t vfac;

  test_abstract_operations();
  return 0;
}
#endif
  
