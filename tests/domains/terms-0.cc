#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  z_term_domain_t dom_left = z_term_domain_t::top ();
  z_term_domain_t dom_right = z_term_domain_t::top ();


  z_var x(vfac["x"]);
  z_var y(vfac["y"]);
  
  dom_left.assign(y, z_number(8));
  dom_left.apply(OP_MULTIPLICATION, x, y, z_number(5));

  dom_right.apply(OP_MULTIPLICATION, x, y, z_number(5));

  z_term_domain_t l_join_r = dom_left | dom_right;

  crab::outs() << dom_left << " | " << dom_right << " = " << l_join_r << "\n";
//  crab::outs() << dom_left;
  return 0;
}
