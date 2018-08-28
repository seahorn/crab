#include "../program_options.hpp"
#include "../common.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  { // expand 
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var A(vfac["A"], crab::INT_TYPE, 32);
    z_var x_prime(vfac["x\'"], crab::INT_TYPE, 32);

    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    // for all i. A[i] >= 0
    dbm += (A >= z_number (0));
    // x = A[..];
    dbm.expand (A, x_prime);
    dbm.assign (x, x_prime);
    dbm -= x_prime;
    // if (x <= 0)
    dbm += (x <= z_number (0));
    crab::outs() << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));    
    dbm += (z_lin_t(x) == z_lin_t(y));
    dbm += (z_lin_t(y) == z_lin_t(z));
    crab::outs() << "Before x != 0: " << dbm << "\n";    
    dbm += (z_lin_t(x) != z_number(0));
    crab::outs() << "After x != 0: " << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var u(vfac["u"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);        
    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));    
    dbm += (z_lin_t(x) == z_lin_t(y));
    dbm += (z_lin_t(u) == z_lin_t(v));    
    dbm += (z_lin_t(y) == z_lin_t(z));
    crab::outs() << "Before x != 10: " << dbm << "\n";    
    dbm += (z_lin_t(x) != z_number(10));
    crab::outs() << "After x != 10: " << dbm << "\n";
  }
  
  return 0;
}
