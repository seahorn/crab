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
    varname_t x = vfac["x"];
    varname_t A = vfac["A"];
    varname_t x_prime = vfac["x\'"];

    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    // for all i. A[i] >= 0
    dbm += (z_var(A) >= z_number (0));
    // x = A[..];
    dbm.expand (A, x_prime);
    dbm.assign (x, z_var (x_prime));
    dbm -= x_prime;
    // if (x <= 0)
    dbm += (z_var(x) <= z_number (0));
    crab::outs() << dbm << "\n";
  }
  return 0;
}
