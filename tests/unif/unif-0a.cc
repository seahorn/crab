#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;


int main (int argc, char** argv )
{
  VariableFactory vfac;

  varname_t x = vfac["x"];
  varname_t y = vfac["y"];

  term_domain_t dom_left = term_domain_t::top ();
  dom_left.assign(x, z_number(1));

  term_domain_t dom_right = dom_left;
  dom_right.apply(OP_ADDITION, x, x, z_number(1));

  term_domain_t l_join_r = dom_left | dom_right;

  std::cout << dom_left << " | " << dom_right << " = " << l_join_r << std::endl;
//  std::cout << dom_left;
  return 0;
}
