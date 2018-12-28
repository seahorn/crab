#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::cfg;
using namespace crab::domain_impl;


int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;

  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  {
    z_term_domain_t dom_left = z_term_domain_t::top ();
    dom_left.assign(x, z_number(1));

    z_term_domain_t dom_right = dom_left;
    dom_right.apply(OP_ADDITION, x, x, z_number(1));
    
    z_term_domain_t l_join_r = dom_left | dom_right;
    
    crab::outs() << dom_left << " | " << dom_right << " = " << l_join_r << "\n";
    //  crab::outs() << dom_left;
  }

  {
    z_var z(vfac["z"], crab::INT_TYPE, 32);

    z_term_domain_t inv = z_term_domain_t::top ();
    inv.assign (x, 5);
    inv.assign (y, 5);
    inv.assign (z, 9);
    std::vector<z_var> vs;
    vs.push_back (x);
    vs.push_back (y);
    crab::outs() << "Before project " << inv << "\n";
    inv.project(vs);
    crab::outs() << "After project " << inv << "\n";


  }
  return 0;
}
