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
  {  
    z_var x(vfac["x"], crab::INT_TYPE);
    z_var A(vfac["A"], crab::INT_TYPE);
    z_var x_prime(vfac["x\'"], crab::INT_TYPE);

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
  return 0;
}
