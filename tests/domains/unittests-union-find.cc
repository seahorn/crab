#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/domains/union_find_domain.hpp>

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
  
  z_var v1(vfac["v1"], crab::INT_TYPE, 32);
  z_var v2(vfac["v2"], crab::INT_TYPE, 32);
  z_var v3(vfac["v3"], crab::INT_TYPE, 32);
  z_var v4(vfac["v4"], crab::INT_TYPE, 32);
  z_var v5(vfac["v5"], crab::INT_TYPE, 32);
  z_var v6(vfac["v6"], crab::INT_TYPE, 32);
  z_var v7(vfac["v7"], crab::INT_TYPE, 32);
  z_var v8(vfac["v8"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  
  { // test all operations

    using union_find_domain_t =
      crab::domains::union_find_domain<z_var, z_interval_domain_t>;    
    z_interval_domain_t idom1,idom2; // some abstract domain
    union_find_domain_t dom1, dom2, dom3;

    idom1.assign(x, 5);
    idom1.assign(y, 10);
    idom2.assign(x,10);
    idom2.assign(y, 20);
    
    dom1.set(v1, idom1);
    dom1.set(v2, idom1);
    dom1.set(v3, idom1);
    dom1.set(v4, idom1);
    crab::outs() << "Dom1=" << dom1 << "\n";    
    dom1.join(v1,v2);
    crab::outs() << "After join " << v1 << " and " << v2 << "=" << dom1 << "\n";
    dom2.set(v1, idom2);
    dom2.set(v2, idom2);
    dom2.set(v3, idom2);
    dom2.set(v4, idom2);
    crab::outs() << "Dom2=" << dom2 << "\n";    
    dom2.join(v3,v4);
    crab::outs() << "After join " << v3 << " and " << v4 << "=" << dom2 << "\n";    
    bool r1 = dom1 <= dom2;
    bool r2 = dom2 <= dom1;
    crab::outs() << "Dom1 <= Dom2 = " << r1 << "\n";
    crab::outs() << "Dom2 <= Dom1 = " << r2 << "\n";
    dom3 = dom1 | dom2;
    crab::outs() << "Dom3 = Dom1 | Dom2 = " << dom3 << "\n";
    bool r3 = dom1 <= dom3;
    bool r4 = dom2 <= dom3;
    crab::outs() << "Dom1 <= Dom3 = " << r3 << "\n";
    crab::outs() << "Dom2 <= Dom3 = " << r4 << "\n";
    union_find_domain_t dom4 = dom1 & dom2;
    crab::outs() << "Dom4 = Dom1 & Dom2 = " << dom4 << "\n";
    
    union_find_domain_t dom5(dom3);
    dom5.forget(v2);
    crab::outs() << "After forgetting " << v2 << ":" << dom5 << "\n";

    union_find_domain_t dom6(dom3);
    dom6.rename({v1,v2,v3,v4}, {v5,v6,v7,v8});
    crab::outs() << "After renaming {v1,v2,v3,v4} with {v5,v6,v7,v8} in Dom3:" << dom6 << "\n";

    union_find_domain_t dom7(dom3);
    dom7.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom3:" << dom7 << "\n";
  }

  // Test meet
  {
    using union_find_domain_t =
      crab::domains::union_find_domain<z_var, boolean_value>;    
    union_find_domain_t dom1, dom2, dom3;
    dom1.set(v1, boolean_value::get_false());
    dom1.set(v2, boolean_value::get_false());
    dom1.set(v3, boolean_value::get_false());
    dom1.set(v4, boolean_value::get_false());
    dom1.join(v1,v2);
    dom1.join(v3,v4);
    crab::outs() << "Dom1=" << dom1 << "\n";
    /////
    dom2.set(v1, boolean_value::get_false());
    dom2.set(v2, boolean_value::get_false());
    dom2.set(v3, boolean_value::get_false());
    dom2.set(v4, boolean_value::get_false());
    dom2.join(v2,v3);
    crab::outs() << "Dom2=" << dom2 << "\n";        
    dom3 = dom1 & dom2;
    crab::outs() << "Dom1 & Dom2 = " << dom3 << "\n";
  }

  {
    using union_find_domain_t =
      crab::domains::union_find_domain<z_var, boolean_value>;        
    union_find_domain_t dom1, dom2, dom3;
    dom1.set(v1, boolean_value::get_false());
    dom1.set(v2, boolean_value::get_false());
    dom1.set(v3, boolean_value::get_false());
    dom1.set(v4, boolean_value::get_false());
    dom1.join(v1,v2);
    dom1.join(v3,v4);
    crab::outs() << "Dom1=" << dom1 << "\n";
    /////
    dom2.set(v5, boolean_value::get_false());
    dom2.set(v6, boolean_value::get_false());
    dom2.set(v7, boolean_value::get_false());
    dom2.set(v8, boolean_value::get_false());
    dom2.join(v6,v7);
    crab::outs() << "Dom2=" << dom2 << "\n";        
    dom3 = dom1 & dom2;
    crab::outs() << "Dom1 & Dom2 = " << dom3 << "\n";
  }


  {
    using union_find_domain_t =
      crab::domains::union_find_domain<z_var, boolean_value>;        
    union_find_domain_t dom1, dom2, dom3;
    dom1.set(v1, boolean_value::get_false());
    dom1.set(v2, boolean_value::get_false());
    dom1.set(v3, boolean_value::get_false());
    dom1.set(v4, boolean_value::get_false());
    dom1.join(v1,v2);
    dom1.join(v3,v4);
    crab::outs() << "Dom1=" << dom1 << "\n";
    /////
    dom2.set(v3, boolean_value::get_false());
    dom2.set(v4, boolean_value::get_false());
    dom2.set(v7, boolean_value::get_false());
    dom2.set(v8, boolean_value::get_false());
    dom2.join(v7,v8);
    crab::outs() << "Dom2=" << dom2 << "\n";        
    dom3 = dom1 & dom2;
    crab::outs() << "Dom1 & Dom2 = " << dom3 << "\n";
  }
  
  

  
  return 0;
}
