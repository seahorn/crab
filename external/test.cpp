// To define a Crab CFG
#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/domains/linear_constraints.hpp>

// To define abstract domains
#include <crab/domains/intervals.hpp>

///////// Begin Crab CFG /////////////
// A variable factory based on strings
typedef crab::cfg::var_factory_impl::str_variable_factory variable_factory_t;
typedef typename variable_factory_t::varname_t varname_t;
// CFG basic block labels
typedef std::string basic_block_label_t;

namespace crab{
namespace cfg_impl{  
template<> inline std::string get_label_str(std::string e) 
{ return e; }
}
}

/// To define CFG over integers
typedef crab::cfg::cfg<basic_block_label_t, varname_t, ikos::z_number> cfg_t;
typedef crab::cfg::cfg_ref<cfg_t> cfg_ref_t;
typedef cfg_t::basic_block_t basic_block_t;
typedef ikos::variable<ikos::z_number, varname_t> var_t;
typedef ikos::linear_expression<ikos::z_number, varname_t> lin_exp_t;
typedef ikos::linear_constraint<ikos::z_number, varname_t> lin_cst_t;
typedef ikos::linear_constraint_system<ikos::z_number, varname_t> lin_cst_sys_t;
///////// End Crab CFG /////////////

///////// Begin Crab Abstract Domains /////////////
typedef crab::domains::interval_domain<ikos::z_number,varname_t> interval_domain_t;
///////// End Crab Abstract Domains /////////////

int main(int argc, char**argv) {

  variable_factory_t vfac;
  var_t x(vfac["x"], crab::INT_TYPE, 32);
  var_t y(vfac["y"], crab::INT_TYPE, 32);

  interval_domain_t inv;
  inv.set(x, 5);
  inv.set(y, 10);
  crab::outs () << inv << "\n";
  
  return 0;
}
