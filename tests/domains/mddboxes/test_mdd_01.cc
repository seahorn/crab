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
  z_var w(vfac["w"], crab::INT_TYPE, 32);        
  z_var z(vfac["z"], crab::INT_TYPE, 32);      
  
  {
    z_mdd_boxes_domain_t inv = z_mdd_boxes_domain_t::top();    
    inv.set(x, z_interval_t(2, 5));
    crab::outs() << "x := [2,5] ---> " << inv << "\n";        
    inv.set(y, z_interval_t(5, 10));
    crab::outs() << "y := [5,10] ---> " << inv << "\n";    
    z_lin_t e(x + y + 3);
    inv.assign(z, e);
    crab::outs() << "z := x+y+3 ---> " << inv << "\n";
    inv.apply(OP_ADDITION, w, x, 3);
    crab::outs() << "w := x + 3 ---> " << inv << "\n";
    inv.apply(OP_SUBTRACTION, w, y, 3);
    crab::outs() << "w := y - 3 ---> " << inv << "\n";    
    inv.apply(OP_MULTIPLICATION, w, y, 3);
    crab::outs() << "w := y * 3 ---> " << inv << "\n";        
  }

  {
    z_mdd_boxes_domain_t inv = z_mdd_boxes_domain_t::top();    
    inv.set(x, z_interval_t(2, 5));
    crab::outs() << "x := [2,5] ---> " << inv << "\n";        
    inv.set(y, z_interval_t(5, 10));
    crab::outs() << "y := [5,10] ---> " << inv << "\n";

    inv.set(z, z_interval_t(8, 10));
    crab::outs() << "z := [8,10] ---> " << inv << "\n";
    
    z_lin_cst_t c1(y <= x);
    inv += c1;
    crab::outs() << "add constraint " << c1 << " ---> " << inv << "\n";

    z_lin_cst_t c2(z <= x);
    inv += c2;
    crab::outs() << "add constraint " << c2 << " ---> " << inv << "\n";
    
  }
  
  {
    z_mdd_boxes_domain_t inv1 = z_mdd_boxes_domain_t::top();
    z_mdd_boxes_domain_t inv2 = z_mdd_boxes_domain_t::top();        
    inv1.set(x, z_interval_t(2, 5));
    inv1.set(y, z_interval_t(10,20));
    crab::outs() << "INV X=" << inv1 << "\n";
    
    inv2.set(x, z_interval_t(7, 10));
    inv2.set(y, z_interval_t(10,20));
    crab::outs() << "INV 2=" << inv2 << "\n";

    z_mdd_boxes_domain_t inv3 = inv1 | inv2;
    crab::outs() << "INV 1 | INV 2 =" << inv3 << "\n";
  }
  
}

int main (int argc, char** argv ) {
  test_abstract_operations();
  return 0;
}
#endif
  
