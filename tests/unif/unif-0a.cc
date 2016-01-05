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

  {
    term_domain_t dom_left = term_domain_t::top ();
    dom_left.assign(x, z_number(1));

    term_domain_t dom_right = dom_left;
    dom_right.apply(OP_ADDITION, x, x, z_number(1));
    
    term_domain_t l_join_r = dom_left | dom_right;
    
    std::cout << dom_left << " | " << dom_right << " = " << l_join_r << std::endl;
    //  std::cout << dom_left;
  }

  {
    varname_t z = vfac["z"];

    term_domain_t inv = term_domain_t::top ();
    inv.assign (x, 5);
    inv.assign (y, 5);
    inv.assign (z, 9);
    std::vector <varname_t> vs;
    vs.push_back (x);
    vs.push_back (y);
    cout << "Before project " << inv << "\n";
    crab::domain_traits::project (inv, vs.begin (), vs.end ());
    cout << "After project " << inv << "\n";


  }
  return 0;
}
