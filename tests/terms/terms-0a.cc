#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  varname_t x = vfac["x"];
  varname_t y = vfac["y"];

  {
    term_domain_t dom_left = term_domain_t::top ();
    dom_left.assign(x, z_number(1));

    term_domain_t dom_right = dom_left;
    dom_right.apply(OP_ADDITION, x, x, z_number(1));
    
    term_domain_t l_join_r = dom_left | dom_right;
    
    crab::outs() << dom_left << " | " << dom_right << " = " << l_join_r << "\n";
    //  crab::outs() << dom_left;
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
    crab::outs() << "Before project " << inv << "\n";
    crab::domains::domain_traits<term_domain_t>::project (inv, vs.begin (), vs.end ());
    crab::outs() << "After project " << inv << "\n";


  }
  return 0;
}
