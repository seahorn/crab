#include "../program_options.hpp"
#include "../common.hpp"
#include <crab/domains/numerical_packing.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

namespace {
  
  using z_packing_domain_t = numerical_packing_domain<z_sdbm_domain_t>;
}

int main(int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  {
    z_packing_domain_t dom;
    dom += (z_lin_exp_t(x) == z_lin_exp_t(y));
    crab::outs() << "EXPECTED={({x,y},{x=y})} ACTUAL=" << dom << "\n";
  }
  
  {
    z_packing_domain_t dom;
    dom.assign(x, z_number(5));
    dom.assign(y, z_number(10));
    //crab::outs() << dom << "\n";
    dom += (z_lin_exp_t(x) == z_lin_exp_t(y));
    crab::outs() << "EXPECTED=_|_ ACTUAL=" << dom << "\n";
  }

  {
    z_packing_domain_t dom1,dom2;
    crab::outs() << "EXPECTED=top ACTUAL=" << dom1 << "\n";
    dom1.assign(x, z_number(5));
    dom2.assign(x, z_number(10));
    dom1 |= dom2;  
    dom1.assign(y, z_number(10));
    dom1.assign(z, z_number(0));  
    crab::outs() << "EXPECTED={({x},{x->[5,10]}),({y},{y->10}),({z},{z->0})} ACTUAL=" <<  dom1 << "\n";
    dom1 += (z_lin_exp_t(x) == z_lin_exp_t(y));
    crab::outs() << "EXPECTED={({x,y},{x->10, y->10}), ({z},{z->0})} ACTUAL=" <<  dom1 << "\n";
    dom1 -= y;
    crab::outs() << "EXPECTED={({x},{x->10}), ({z},{z->0})} ACTUAL=" <<  dom1 << "\n";
  }  

 {
    z_packing_domain_t dom;
    dom.assign(w, z_number(5));
    dom.assign(z, z_number(10));
    crab::outs() << "EXPECTED={({w},{w->5}),({z},{z->10})} ACTUAL=" << dom << "\n";
    dom += (z_lin_exp_t(x) == z_lin_exp_t(y));
    crab::outs() << "EXPECTED={({w},{w->5}),({z},{z->10}), ({x,y},{x=y})} ACTUAL=" << dom << "\n";
  }
 
  return 0;
}
