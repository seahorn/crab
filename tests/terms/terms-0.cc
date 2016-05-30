#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  term_domain_t dom_left = term_domain_t::top ();
  term_domain_t dom_right = term_domain_t::top ();


  varname_t x = vfac["x"];
  varname_t y = vfac["y"];
  
  dom_left.assign(y, z_number(8));
  dom_left.apply(OP_MULTIPLICATION, x, y, z_number(5));

  dom_right.apply(OP_MULTIPLICATION, x, y, z_number(5));

  term_domain_t l_join_r = dom_left | dom_right;

  crab::outs() << dom_left << " | " << dom_right << " = " << l_join_r << "\n";
//  crab::outs() << dom_left;
  return 0;
}
