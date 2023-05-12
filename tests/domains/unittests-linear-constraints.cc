#include "../common.hpp"
#include "../program_options.hpp"

#include "crab/domains/discrete_domains.hpp"

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
  using linear_constraint_system_t =
    ikos::linear_constraint_system<ikos::z_number, varname_t>;

  { // linear constraints normalization
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

  { 
    struct linear_constraint_compare {
      bool operator()(const z_lin_cst_t &c1, const z_lin_cst_t &c2) const {
	return c1.lexicographical_compare(c2);
      }
    };
    using set_domain_t = set_domain<z_lin_cst_t, linear_constraint_compare>;
    
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);
    z_var w(vfac["w"], crab::INT_TYPE, 32);
    z_var v(vfac["v"], crab::INT_TYPE, 32);

    set_domain_t dom1 = set_domain_t::bottom();
    dom1 += z_lin_cst_t(x <= 8);
    crab::outs() << "After adding x <= 8: " << dom1 << "\n";
    dom1 += z_lin_cst_t(y >= 5);
    crab::outs() << "After adding y >= 5: " << dom1 << "\n";    
    dom1 += z_lin_cst_t(y <= 6);
    crab::outs() << "After adding y <= 6: " << dom1 << "\n";        
    dom1 += z_lin_cst_t(x >= 8);
    crab::outs() << "After adding x >= 8: " << dom1 << "\n";
    dom1 += z_lin_cst_t(x <= 8);
    crab::outs() << "After adding x <= 8: " << dom1 << "\n";
    crab::outs() << "--------------\n";

    set_domain_t dom2 = set_domain_t::bottom();
    dom2 += z_lin_cst_t(x <= 8);
    dom2 += z_lin_cst_t(y >= 5);
    
    crab::outs() << "Checking if " << dom1 << " <= " << dom2 << "=" << (dom1 <= dom2) << "\n";
    crab::outs() << "Checking if " << dom2 << " <= " << dom1 << "=" << (dom2 <= dom1) << "\n";    
    crab::outs() << "--------------\n";

    set_domain_t dom3 = dom1 & dom2;
    crab::outs() << dom1  << " & " << dom2 << "=" << dom3 << "\n";
    crab::outs() << "--------------\n";


    dom2 += z_lin_cst_t(z >= 4);
    set_domain_t dom4 = dom1 | dom2;
    crab::outs() << dom1  << " | " << dom2 << "=" << dom4 << "\n";
    crab::outs() << "--------------\n";

    set_domain_t dom5 = set_domain_t::top();
    set_domain_t dom6 = dom1 | dom5;
    crab::outs() << dom1  << " | " << dom5 << "=" << dom6 << "\n";
    crab::outs() << "--------------\n";
    
  }

  {
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    linear_constraint_system_t csts;
    csts += z_lin_cst_t(x >= 0);
    csts += z_lin_cst_t::get_false();
    crab::outs() << csts << "\n";
    crab::outs() << "is false=" << csts.is_false() << "\n";
    crab::outs() << "--------------\n";    
  }
  return 0;
}
