/* This test case does not show a bug but an imprecision of the
   widening operator with rationals */

#include "../common.hpp"
#include "../program_options.hpp"
using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

#ifdef HAVE_LDD
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);      
  z_var z(vfac["z"], crab::INT_TYPE, 32);  

   {
     crab::outs() << "-------------------\n";     
    z_boxes_domain_t inv;
    inv += z_lin_cst_t(y >= 0);
    inv += z_lin_cst_t(y <= 9);
    crab::outs() << inv << "\n";        
    inv.apply(OP_ADDITION, z, y, 24);
    crab::outs() << z <<  ":=" << y << "+ 24:\n";
    crab::outs() << "expected= z = [24, 33]\n" ;
    crab::outs() << "result=" << inv << "\n";
   }  
   {
     crab::outs() << "-------------------\n";     
     z_boxes_domain_t inv;
     inv += z_lin_cst_t(y >= 0);
     inv += z_lin_cst_t(y <= 9);
     crab::outs() << inv << "\n";          
     inv.apply(OP_MULTIPLICATION, z, y, 24);
     crab::outs() << z <<  ":=" << y << "* 24:\n";     
     crab::outs() << "expected= z = [0, 216]\n" ;
     crab::outs() << "result=" << inv << "\n";
   }

   {
     crab::outs() << "-------------------\n";     
     z_boxes_domain_t inv;
     inv += z_lin_cst_t(y >= 0);
     inv += z_lin_cst_t(y <= 9);
     inv += z_lin_cst_t(x == 24);
     crab::outs() << inv << "\n";          
     inv.apply(OP_MULTIPLICATION, z, y, x);
     crab::outs() << z <<  ":=" << y << "*" << x << ":\n";          
     crab::outs() << "expected= z = [0, 216]\n" ;
     crab::outs() << "result=" << inv << "\n";
   }

   {
     crab::outs() << "-------------------\n";     
     z_boxes_domain_t inv;
     inv += z_lin_cst_t(y >= -1);
     inv += z_lin_cst_t(y <= 9);
     crab::outs() << inv << "\n";          
     inv.apply(OP_MULTIPLICATION, z, y, 24);
     crab::outs() << z <<  ":=" << y << "* 24:\n";               
     crab::outs() << "expected= z = [-24, 216]\n" ;
     crab::outs() << "result=" << inv << "\n";
   }
   
   {
     crab::outs() << "-------------------\n";
     z_boxes_domain_t inv;
     inv += z_lin_cst_t(y >= 0);
     inv += z_lin_cst_t(y <= 80);
     crab::outs() << inv << "\n";          
     inv.apply(OP_SDIV, z, y, 4);
     crab::outs() << z <<  ":=" << y << "/ 4:\n";          
     crab::outs() << "expected= z = [0, 20]\n" ;
     crab::outs() << "result=" << inv << "\n";
   }

   {
     crab::outs() << "-------------------\n";     
     z_boxes_domain_t inv;
     inv += z_lin_cst_t(y >= 0);
     inv += z_lin_cst_t(y <= 80);
     inv += z_lin_cst_t(x == 4);
     crab::outs() << inv << "\n";     
     inv.apply(OP_SDIV, z, y, x);
     crab::outs() << z <<  ":=" << y << "/" << x << ":\n";          
     crab::outs() << "expected= z = [0, 20]\n" ;
     crab::outs() << "result=" << inv << "\n";
   }
#endif
  return 0;
}
