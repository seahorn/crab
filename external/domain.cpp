#include <crab/config.h>
#include <crab/types/linear_constraints.hpp>
#include <crab/types/varname_factory.hpp>
#include <crab/types/variable.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

// Abstract domains
#include <crab/domains/intervals.hpp>
#include <crab/domains/split_dbm.hpp>

using namespace crab;
using namespace ikos;

// A variable factory based on strings
using variable_factory_t = var_factory_impl::str_variable_factory;
namespace crab {
template<>
class variable_name_traits<std::string> {
public:
  static std::string to_string(std::string varname) {
    return varname;
  }
};
} //end namespace crab
// Expressions
using varname_t = typename variable_factory_t::varname_t;
using var_t = variable<z_number, varname_t>;
using lin_exp_t = linear_expression<z_number, varname_t>;
using lin_cst_t = linear_constraint<z_number, varname_t> ;
using lin_cst_sys_t = linear_constraint_system<z_number, varname_t> ;

///////// Begin Crab Abstract Domains /////////////
using interval_domain_t = interval_domain<z_number,varname_t>;
using zones_domain_t = domains::split_dbm_domain<z_number, varname_t>;
///////// End Crab Abstract Domains /////////////

int main(int argc, char**argv) {

  variable_factory_t vfac;
  var_t x(vfac["x"], INT_TYPE, 32);
  var_t y(vfac["y"], INT_TYPE, 32);
  var_t z(vfac["z"], INT_TYPE, 32);

  {
    outs() << "Example using intervals\n";
    interval_domain_t inv1, inv2;
    inv1.assign(x, 5);
    inv1.assign(y, 10);
    outs() << "inv1=" << inv1 << "\n";  
    inv2.assign(x, 10);
    inv2.assign(y, 20);
    outs() << "inv2=" << inv1 << "\n";    
    interval_domain_t inv3 = inv1 | inv2;
    inv3.apply(domains::OP_ADDITION, z, x, y);
    outs() << "inv1 | inv2 = " << inv3 << "\n";
  }

  {
    outs() << "Example using zones\n";
    zones_domain_t inv1, inv2;
    inv1.assign(x, 5);
    inv1.assign(y, 10);
    outs() << "inv1=" << inv1 << "\n";  
    inv2.assign(x, 10);
    inv2.assign(y, 20);
    outs() << "inv2=" << inv1 << "\n";    
    zones_domain_t inv3 = inv1 | inv2;
    inv3.apply(domains::OP_ADDITION, z, x, y);
    outs() << "inv1 | inv2 = " << inv3 << "\n";
  }
  
  return 0;
}
