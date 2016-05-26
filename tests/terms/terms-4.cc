#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  typedef interval< z_number> interval_t;

  varname_t x = vfac["x"];
  varname_t y = vfac["y"];
  varname_t w = vfac["w"];
  varname_t z = vfac["z"];

  {

    term_domain_t dom_left = term_domain_t::top ();
    term_domain_t dom_right = term_domain_t::top ();

    // ({w = a0, x = a0, y = '+' (a0,a1), z = a1}, { x=5, w=5, z=3, y=8 })
    dom_left.assign (x, 5);
    dom_left.assign (w, z_lin_t (x));
    dom_left.assign (z, 3);
    dom_left.apply(OP_ADDITION, y, x, z);
    
    // ({w = b0,  x = '+' (b0,b1), y = b0, z = b1}, {y=8, w=8,z=2,x=10 })
    dom_right.assign (y, 8); 
    dom_right.assign (w, z_lin_t (y));
    dom_right.assign (z, 2); 
    dom_right.apply(OP_ADDITION, x, w, z);

    // meet = ({x=y=w= '+' (c0,c1), z=c2},{_|_}) = _|_
    crab::outs() << "Meet" << "\n" << dom_left << " \n " << dom_right << "\n"; 
    term_domain_t l_meet_r = dom_left & dom_right;
    crab::outs() << "Result=" << l_meet_r << "\n";
  }

  {
    term_domain_t dom_left = term_domain_t::top ();
    term_domain_t dom_right = term_domain_t::top ();

    // ({w = a0, x = a0, y = '+' (a0,a1), z = a1}, {x=[5,8],w=[5,8],z=[1,10],y=[6,18]})
    dom_left.set (x, interval_t (5,8));
    dom_left.assign (w, z_lin_t (x));
    dom_left.set (z, interval_t (1,10));
    dom_left.apply(OP_ADDITION, y, x, z);
    
    // ({w = b0,  x = '+' (b0,b1), y = b0, z = b1}, {y=[2,7],w=[2,7],z=[3,5],x=[5,12]})
    dom_right.set (y, interval_t (2,7)); 
    dom_right.assign (w, z_lin_t (y));
    dom_right.set (z, interval_t (3,5)); 
    dom_right.apply(OP_ADDITION, x, w, z);

    // meet = ({x=y=w= '+' (c0,c1), z=c2},{x=[5,8], y=[6,7], z=[3,5], w=[5,7]})
    crab::outs() << "Meet" << "\n" << dom_left << " \n " << dom_right << "\n"; 
    term_domain_t l_meet_r = dom_left & dom_right;
    crab::outs() << "Result=" << l_meet_r << "\n";
  }

  {
    term_domain_t dom = term_domain_t::top ();
    varname_t zero = vfac["v0"];
    varname_t one = vfac["v1"];

    dom.set (zero, interval_t (0,0));
    dom.set (one, interval_t (1,1));

    dom.apply(OP_ADDITION, x, one, zero);
    dom.apply(OP_ADDITION, y, zero, one);
    z_lin_cst_t c1 (z_lin_t(x) == z_lin_t(y));
    crab::outs() << "Added " << c1 << "\n";
    crab::outs() << dom << "\n";
    dom += c1;
    crab::outs() << "Result=" << dom << "\n";
    z_lin_cst_t c2 (z_lin_t(x) != z_lin_t(y));
    crab::outs() << "Added " << c2 << "\n";
    dom += c2;
    crab::outs() << "Result=" << dom << "\n";
  }

  return 0;
}
