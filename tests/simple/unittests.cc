#include <ikos/tests/Cfg_impl.hpp>
#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/intervals_congruences.hpp>                      
#include <ikos/domains/octagons.hpp>                      
#include <ikos/domains/dbm.hpp>                      

using namespace std;

namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains
  typedef interval_domain< z_number, varname_t > interval_domain_t;
  typedef interval_congruence_domain< z_number, varname_t > ric_t;
  typedef DBM<z_number, varname_t> dbm_domain_t;
  typedef octagon< z_number, varname_t > octagon_domain_t;


} // end namespace

using namespace cfg_impl;
using namespace domain_impl;

typedef linear_constraint<z_number, varname_t> linear_constraint_t;
typedef linear_expression<z_number, varname_t> linear_expression_t;

int main (int argc, char** argv )
{
  VariableFactory vfac;



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
    cout << dbm << endl;
  }
       
  return 0;
}
