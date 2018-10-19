// Crab CFG stuff
#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/domains/linear_constraints.hpp>

// Abstract domains
#include <crab/domains/intervals.hpp>
#include <crab/domains/apron_domains.hpp>

using namespace crab;
using namespace crab::domains;
using namespace crab::cfg;
using namespace ikos;

///////// Begin Crab CFG /////////////
// A variable factory based on strings
typedef var_factory_impl::str_variable_factory variable_factory_t;
typedef typename variable_factory_t::varname_t varname_t;
// CFG basic block labels
typedef std::string basic_block_label_t;
/// To define CFG over integers
typedef cfg<basic_block_label_t, varname_t, z_number> cfg_t;
typedef cfg_ref<cfg_t> cfg_ref_t;
typedef cfg_t::basic_block_t basic_block_t;
typedef variable<z_number, varname_t> var_t;
typedef linear_expression<z_number, varname_t> lin_exp_t;
typedef linear_constraint<z_number, varname_t> lin_cst_t;
typedef linear_constraint_system<z_number, varname_t> lin_cst_sys_t;
///////// End Crab CFG /////////////

///////// Begin Crab Abstract Domains /////////////
typedef interval_domain<z_number,varname_t> interval_domain_t;
typedef apron_domain<z_number,varname_t, apron_domain_id_t::APRON_PK> pk_domain_t;
///////// End Crab Abstract Domains /////////////

int main(int argc, char**argv) {

  variable_factory_t vfac;
  var_t x(vfac["x"], INT_TYPE);
  var_t y(vfac["y"], INT_TYPE);
  var_t z(vfac["z"], INT_TYPE);

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
    inv3.apply(OP_ADDITION, z, x, y);
    outs() << "inv1 | inv2 = " << inv3 << "\n";
  }

  {
    outs() << "Example using polyhedra\n";
    pk_domain_t inv1, inv2;
    inv1.assign(x, 5);
    inv1.assign(y, 10);
    outs() << "inv1=" << inv1 << "\n";  
    inv2.assign(x, 10);
    inv2.assign(y, 20);
    outs() << "inv2=" << inv1 << "\n";    
    pk_domain_t inv3 = inv1 | inv2;
    inv3.apply(OP_ADDITION, z, x, y);
    outs() << "inv1 | inv2 = " << inv3 << "\n";
  }
  
  return 0;
}
