#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  { // expand
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var A(vfac["A"], crab::INT_TYPE, 32);
    z_var x_prime(vfac["x\'"], crab::INT_TYPE, 32);

    z_dbm_domain_t dbm;
    // for all i. A[i] >= 0
    dbm += (A >= z_number(0));
    // x = A[..];
    dbm.expand(A, x_prime);
    dbm.assign(x, x_prime);
    dbm -= x_prime;
    // if (x <= 0)
    dbm += (x <= z_number(0));
    crab::outs() << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm;
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));
    dbm += (z_lin_exp_t(x) == z_lin_exp_t(y));
    dbm += (z_lin_exp_t(y) == z_lin_exp_t(z));
    crab::outs() << "Before x != 0: " << dbm << "\n";
    dbm += (z_lin_exp_t(x) != z_number(0));
    crab::outs() << "After x != 0: " << dbm << "\n";
  }

  { // disequalities
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var u(vfac["u"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm;
    dbm += (x >= z_number(0));
    dbm += (x <= z_number(10));
    dbm += (z_lin_exp_t(x) == z_lin_exp_t(y));
    dbm += (z_lin_exp_t(u) == z_lin_exp_t(v));
    dbm += (z_lin_exp_t(y) == z_lin_exp_t(z));
    crab::outs() << "Before x != 10: " << dbm << "\n";
    dbm += (z_lin_exp_t(x) != z_number(10));
    crab::outs() << "After x != 10: " << dbm << "\n";
  }

  { // meet
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_dbm_domain_t dbm1;
    dbm1 += (x == z_lin_exp_t(y));
    dbm1 += (x == z_number(1));
    crab::outs() << "DBM1=" << dbm1 << "\n";
    z_dbm_domain_t dbm2;
    dbm2 += (y >= z_number(1));
    crab::outs() << "DBM2=" << dbm2 << "\n";
    z_dbm_domain_t dbm3 = dbm2 & dbm1;
    crab::outs() << "DBM1 & DBM2=" << dbm3 << "\n";
  }

  { // linear constraints normalization (should not be here)
    using linear_constraint_system_t =
        ikos::linear_constraint_system<ikos::z_number, varname_t>;

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
    using split_dbm_domain_t = split_dbm_domain<
        ikos::z_number, varname_t,
        DBM_impl::DefaultParams<ikos::z_number, DBM_impl::GraphRep::ss>>;

    z_var x(vfac["x"], crab::INT_TYPE, 32);

    split_dbm_domain_t d1, d2;
    z_lin_cst_t c1(x == z_number("-9223372036854775808"));
    crab::outs() << "Adding constraint 1 " << c1 << "\n";
    d1 += c1;
    crab::outs() << "Dom1=" << d1 << "\n";
    auto csts1 = d1.to_linear_constraint_system();
    crab::outs() << "Csts1=" << csts1 << "\n";
    /////////////////////////
    z_lin_cst_t c2(x == z_number("-9223372036854775807"));
    crab::outs() << "Adding constraint 2 " << c2 << "\n";
    d2 += c2;
    crab::outs() << "Dom2=" << d2 << "\n";
    auto csts2 = d2.to_linear_constraint_system();
    crab::outs() << "Csts2=" << csts2 << "\n";
  }
  {
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);

    z_sdbm_domain_t inv;
    inv += (y >= 0);
    inv.apply(OP_AND, x, y, 9223372036854775296);
    crab::outs() << inv << "\n";
  }

  {  // Tests the other graph representations: at least making sure
     // the code compiles.
    using ss_graph_t = DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::ss>;
    using pt_graph_t = DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::pt>;
    using ht_graph_t = DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::ht>;        
      
    using ss_dbm = split_dbm_domain<z_number, varname_t, ss_graph_t>;
    using pt_dbm = split_dbm_domain<z_number, varname_t, pt_graph_t>;
    using ht_dbm = split_dbm_domain<z_number, varname_t, ht_graph_t>;    

    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);

    z_lin_cst_t cst(z_lin_exp_t(x) == z_lin_exp_t(y));
    
    ss_dbm ss_inv_1, ss_inv_2, ss_inv_3, ss_inv_4;
    ss_inv_1.assign(x, 5);
    ss_inv_1 += cst;
    ss_inv_2.assign(x, 10);
    ss_inv_2 += cst;
    ss_inv_3 = ss_inv_1 | ss_inv_2;
    ss_inv_4 = ss_inv_1 & ss_inv_2;    
    crab::outs() << "Join of " << ss_inv_1 << " and " << ss_inv_2 << "=" << ss_inv_3 << "\n";
    crab::outs() << "Meet of " << ss_inv_1 << " and " << ss_inv_2 << "=" << ss_inv_4 << "\n";    

    pt_dbm pt_inv_1, pt_inv_2, pt_inv_3, pt_inv_4;
    pt_inv_1.assign(x, 5);
    pt_inv_1 += cst;
    pt_inv_2.assign(x, 10);
    pt_inv_2 += cst;
    pt_inv_3 = pt_inv_1 | pt_inv_2;
    pt_inv_4 = pt_inv_1 & pt_inv_2;    
    crab::outs() << "Join of " << pt_inv_1 << " and " << pt_inv_2 << "=" << pt_inv_3 << "\n";
    crab::outs() << "Meet of " << pt_inv_1 << " and " << pt_inv_2 << "=" << pt_inv_4 << "\n";    

    ht_dbm ht_inv_1, ht_inv_2, ht_inv_3, ht_inv_4;
    ht_inv_1.assign(x, 5);
    ht_inv_1 += cst;
    ht_inv_2.assign(x, 10);
    ht_inv_2 += cst;
    ht_inv_3 = ht_inv_1 | ht_inv_2;
    ht_inv_4 = ht_inv_1 & ht_inv_2;    
    crab::outs() << "Join of " << ht_inv_1 << " and " << ht_inv_2 << "=" << ht_inv_3 << "\n";
    crab::outs() << "Meet of " << ht_inv_1 << " and " << ht_inv_2 << "=" << ht_inv_4 << "\n";    
    
  }

  
  return 0;
}
