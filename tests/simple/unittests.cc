#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

typedef linear_constraint<z_number, varname_t> linear_constraint_t;
typedef linear_expression<z_number, varname_t> linear_expression_t;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  {  
    varname_t x = vfac["x"];
    varname_t A = vfac["A"];
    varname_t x_prime = vfac["x\'"];

    dbm_domain_t dbm = dbm_domain_t::top ();
    // for all i. A[i] >= 0
    dbm += linear_constraint_t ( linear_expression_t (A) >= z_number (0));
    // x = A[..];
    dbm.expand (A, x_prime);
    dbm.assign (x, z_var (x_prime));
    dbm -= x_prime;
    // if (x <= 0)
    dbm += linear_constraint_t ( linear_expression_t (x) <= z_number (0));
    crab::outs() << dbm << "\n";
  }
  return 0;
}
