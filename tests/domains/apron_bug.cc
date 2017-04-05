#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

namespace crab
{
   namespace cfg_impl
   {
     typedef Cfg< basic_block_label_t, varname_t, q_number> q_cfg_t;
     typedef q_cfg_t::basic_block_t                         q_basic_block_t;
     
     typedef cfg_ref<q_cfg_t>                               q_cfg_ref_t;
     typedef cfg_rev<q_cfg_ref_t>                           q_cfg_rev_t;

   }
   namespace domain_impl
   {
    
     using namespace crab::cfg_impl;
     using namespace crab::domains; 
     using namespace ikos;
     
     typedef variable< q_number, varname_t> q_var;
     typedef linear_expression<q_number, varname_t> q_lin_t;
     typedef linear_constraint<q_number, varname_t> q_lin_cst_t;
     typedef linear_constraint_system<q_number, varname_t> q_lin_cst_sys_t;
     typedef interval<q_number> q_interval_t;
     typedef bound<q_number> q_bound_t;

     // Numerical domains
     typedef interval_domain<q_number, varname_t> q_interval_domain_t;     
     typedef apron_domain<q_number, varname_t, apron_domain_id_t::APRON_PK> q_pk_apron_domain_t;
     typedef boxes_domain< q_number, varname_t > q_boxes_domain_t;     
   }
}

/* Example of how to infer invariants from the above CFG */
int main (int argc, char** argv ) {

#ifdef HAVE_APRON
  variable_factory_t vfac;
  {
    q_pk_apron_domain_t inv = q_pk_apron_domain_t::top ();
    q_lin_cst_sys_t csts;
    csts += (q_lin_t (vfac ["x"]) >= q_number(1));
    csts += (q_lin_t (vfac ["x"]) <= q_number(1));
    crab::outs () << "INITIALLY=" << inv << "\n";
    crab::outs () << "ADDING CONSTRAINTS=" << csts << "\n";
    inv += csts;
    crab::outs () << "EXPECTED={x = 1}\n";
    crab::outs () << "RESULT=" << inv << "\n";
  }
  crab::outs () << "------------------------------------\n";
  {
    q_pk_apron_domain_t inv = q_pk_apron_domain_t::top ();
    q_lin_cst_sys_t csts;
    csts += (q_lin_t (vfac ["x"]) >= q_number(0.5));
    csts += (q_lin_t (vfac ["x"]) <= q_number(0.5));
    crab::outs () << "INITIALLY=" << inv << "\n";
    crab::outs () << "ADDING CONSTRAINTS=" << csts << "\n";
    inv += csts;
    crab::outs () << "EXPECTED={x = 1/2}\n";
    crab::outs () << "RESULT=" << inv << "\n";
  }
#endif  
  return 0;
}

