#include <ikos/tests/Cfg_impl.hpp>
#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/intervals_congruences.hpp>                      
#include <ikos/domains/octagons.hpp>                      
#include <ikos/domains/dbm.hpp>                      
#include <ikos/domains/term_equiv.hpp>

using namespace std;

namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains
  typedef interval_domain< z_number, varname_t > interval_domain_t;
  typedef ikos::term::TDomInfo<z_number, varname_t, interval_domain_t> dom_info_t;
  typedef anti_unif<dom_info_t>::anti_unif_t term_domain_t;
} // end namespace

using namespace cfg_impl;
using namespace domain_impl;

int main (int argc, char** argv )
{
  VariableFactory vfac;

  term_domain_t dom_left = term_domain_t::top ();
  term_domain_t dom_right = term_domain_t::top ();


  varname_t x = vfac["x"];
  varname_t y = vfac["y"];
  
  dom_left.assign(y, z_number(8));
  dom_left.apply(OP_MULTIPLICATION, x, y, z_number(5));

  dom_right.apply(OP_MULTIPLICATION, x, y, z_number(5));

  term_domain_t l_join_r = dom_left | dom_right;

  std::cout << dom_left << " | " << dom_right << " = " << l_join_r << std::endl;
//  std::cout << dom_left;
  return 0;
}
