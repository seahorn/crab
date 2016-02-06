#ifndef __TESTS_COMMON__
#define __TESTS_COMMON__

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/cfg/Cfg.hpp>
#include <crab/cfg/CfgBgl.hpp> 
#include <crab/cfg/VarFactory.hpp>
#include <crab/analysis/FwdAnalyzer.hpp>
#include <crab/analysis/Pointer.hpp>
#include <crab/analysis/Liveness.hpp>

#include <crab/domains/linear_constraints.hpp> 
#include <crab/domains/intervals.hpp>                      
#include <crab/domains/dbm.hpp>                      
#include <crab/domains/split_dbm.hpp>
#include <crab/domains/var_packing_naive_dbm.hpp>
#include <crab/domains/boxes.hpp>                      
#include <crab/domains/apron_domains.hpp>                      
#include <crab/domains/dis_intervals.hpp>
#include <crab/domains/term_equiv.hpp>
#include <crab/domains/array_graph.hpp>                      
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/combined_domains.hpp>                      

namespace crab {

  namespace cfg_impl {

    using namespace cfg;
    using namespace std;

    template<> inline std::string get_label_str(std::string e) 
    { return e; }

    // A variable factory based on strings
    typedef cfg::var_factory_impl::StrVariableFactory VariableFactory;
    typedef typename VariableFactory::varname_t varname_t;

    // CFG
    typedef variable< z_number, varname_t >      z_var;
    typedef std::string                          basic_block_label_t;
    typedef Cfg< basic_block_label_t, varname_t> cfg_t;
    typedef cfg_t::basic_block_t                 basic_block_t;
  }

  namespace domain_impl {
    using namespace crab::cfg_impl;
    using namespace crab::domains; 
    using namespace ikos;
    typedef linear_expression<z_number, varname_t> z_lin_t;
    typedef linear_constraint<z_number, varname_t> z_lin_cst_t;
    typedef linear_constraint_system<z_number, varname_t> z_lin_cst_sys_t;
    typedef pointer_constraint <varname_t> ptr_cst_t;
    typedef interval<z_number> z_interval_t;
    typedef bound<z_number> z_bound_t;
    // Numerical domains
    typedef interval_domain< z_number, varname_t > interval_domain_t;
    typedef numerical_congruence_domain< interval_domain_t> ric_domain_t;
    typedef DBM<z_number, varname_t> dbm_domain_t;
    typedef SplitDBM<z_number, varname_t> sdbm_domain_t;
    typedef var_packing_naive_dbm<z_number, varname_t> pdbm_domain_t;
    typedef boxes_domain< z_number, varname_t > boxes_domain_t;
    typedef dis_interval_domain<z_number, varname_t > dis_interval_domain_t;
    typedef apron_domain< z_number, varname_t, apron_domain_id_t::APRON_INT > box_apron_domain_t;
    typedef apron_domain< z_number, varname_t, apron_domain_id_t::APRON_OCT > oct_apron_domain_t;
    typedef apron_domain< z_number, varname_t, apron_domain_id_t::APRON_OPT_OCT > opt_oct_apron_domain_t;
    typedef apron_domain< z_number, varname_t, apron_domain_id_t::APRON_PK > pk_apron_domain_t;
    typedef term_domain<term::TDomInfo<z_number, varname_t, interval_domain_t> > term_domain_t;
    typedef term_domain<term::TDomInfo<z_number, varname_t, dbm_domain_t> > term_dbm_t;
    typedef term_domain<term::TDomInfo<z_number, varname_t, dis_interval_domain_t> > term_dis_int_t;
    typedef reduced_numerical_domain_product2<term_dis_int_t, sdbm_domain_t> num_domain_t; 
    // Array domains
    typedef array_graph_domain<dbm_domain_t, interval_domain_t> array_graph_domain_t;
    typedef array_smashing<dis_interval_domain_t> array_smashing_t;
  } 

}

#endif 
