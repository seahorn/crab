#include "../program_options.hpp"
#include "../common.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  { // expand 
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var A(vfac["A"], crab::INT_TYPE, 32);
    z_var x_prime(vfac["x\'"], crab::INT_TYPE, 32);

    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    // for all i. A[i] >= 0
    dbm += (A >= z_number (0));
    // x = A[..];
    dbm.expand (A, x_prime);
    dbm.assign (x, x_prime);
    dbm -= x_prime;
    // if (x <= 0)
    dbm += (x <= z_number (0));
    crab::outs() << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));    
    dbm += (z_lin_t(x) == z_lin_t(y));
    dbm += (z_lin_t(y) == z_lin_t(z));
    crab::outs() << "Before x != 0: " << dbm << "\n";    
    dbm += (z_lin_t(x) != z_number(0));
    crab::outs() << "After x != 0: " << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var u(vfac["u"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);        
    z_dbm_domain_t dbm = z_dbm_domain_t::top ();
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));    
    dbm += (z_lin_t(x) == z_lin_t(y));
    dbm += (z_lin_t(u) == z_lin_t(v));    
    dbm += (z_lin_t(y) == z_lin_t(z));
    crab::outs() << "Before x != 10: " << dbm << "\n";    
    dbm += (z_lin_t(x) != z_number(10));
    crab::outs() << "After x != 10: " << dbm << "\n";
  }

  { // meet
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm1 = z_dbm_domain_t::top ();
    dbm1 += (x == z_lin_t(y));
    dbm1 += (x == z_number(1));
    crab::outs() << "DBM1=" << dbm1 << "\n";
    z_dbm_domain_t dbm2 = z_dbm_domain_t::top ();
    dbm2 += (y >= z_number(1));
    crab::outs() << "DBM2=" << dbm2 << "\n";    
    z_dbm_domain_t dbm3 = dbm2 & dbm1;
    crab::outs() << "DBM1 & DBM2=" << dbm3 << "\n";
  }

  { // linear constraints normalization (should not be here)
    typedef ikos::linear_constraint_system<ikos::z_number, varname_t>
      linear_constraint_system_t;

    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var w(vfac["w"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);    
    
    linear_constraint_system_t csts;
    csts += z_lin_cst_t(x <= 8);
    csts += z_lin_cst_t(y >= 5);
    csts += z_lin_cst_t(y <= 6);
    csts += z_lin_cst_t(x >= 8);
    csts += z_lin_cst_t(v >= 9888);    
    csts += z_lin_cst_t(z == 10);
    csts += z_lin_cst_t(w <= 0);
    csts += z_lin_cst_t(v <= 9888);        
    crab::outs() << "Before normalize constraints " << csts << "\n";
    crab::outs() << "After normalize constraints " << csts.normalize() << "\n";    
  }
  
  { // overflow cases with zones domain
    typedef SplitDBM<ikos::z_number,
		     varname_t,
		     DBM_impl::DefaultParams<ikos::z_number,
					     DBM_impl::GraphRep::ss>>
      SplitDBM_t;
    
    z_var x(vfac["x"], crab::INT_TYPE, 32);

    SplitDBM_t d1, d2;    
    z_lin_cst_t c1(x == z_number("-9223372036854775808"));
    crab::outs () << "Adding constraint 1 " << c1 << "\n";
    d1 += c1;
    crab::outs() << "Dom1=" << d1 << "\n";
    auto csts1 = d1.to_linear_constraint_system();
    crab::outs() << "Csts1=" << csts1 << "\n";
    /////////////////////////
    z_lin_cst_t c2(x == z_number("-9223372036854775807"));
    crab::outs () << "Adding constraint 2 " << c2 << "\n";
    d2 += c2;
    crab::outs() << "Dom2=" << d2 << "\n";
    auto csts2 = d2.to_linear_constraint_system();
    crab::outs() << "Csts2=" << csts2 << "\n";
    
  }
  
  return 0;
}
