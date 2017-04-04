#ifndef __TESTS_COMMON__
#define __TESTS_COMMON__

/* To be included by all the tests */

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cfg/cfg_bgl.hpp> 
#include <crab/cfg/var_factory.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/liveness.hpp>

#include <crab/domains/linear_constraints.hpp> 
#include <crab/domains/intervals.hpp>                      
#include <crab/domains/sparse_dbm.hpp>                      
#include <crab/domains/split_dbm.hpp>
#include <crab/domains/boxes.hpp>                      
#include <crab/domains/apron_domains.hpp>                      
#include <crab/domains/dis_intervals.hpp>
#include <crab/domains/term_equiv.hpp>
#include <crab/domains/array_sparse_graph.hpp>                      
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/nullity.hpp>
#include <crab/domains/flat_boolean_domain.hpp>                      
#include <crab/domains/combined_domains.hpp>                      

#include <boost/program_options.hpp>

namespace crab {

  namespace cfg_impl {

    using namespace cfg;
    using namespace std;

    template<> inline std::string get_label_str(std::string e) 
    { return e; }

    // A variable factory based on strings
    typedef cfg::var_factory_impl::str_variable_factory variable_factory_t;
    typedef typename variable_factory_t::varname_t varname_t;

    // CFG basic block labels
    typedef std::string basic_block_label_t;

    /// CFG over integers
    typedef Cfg< basic_block_label_t, varname_t, ikos::z_number> z_cfg_t;
    typedef cfg_ref<z_cfg_t> z_cfg_ref_t;
    typedef cfg_rev<z_cfg_ref_t> z_cfg_rev_t;
    typedef z_cfg_t::basic_block_t z_basic_block_t;
    /// CFG over rationals
    typedef Cfg< basic_block_label_t, varname_t, ikos::q_number> q_cfg_t;
    typedef q_cfg_t::basic_block_t q_basic_block_t;
    typedef cfg_ref<q_cfg_t> q_cfg_ref_t;
    typedef cfg_rev<q_cfg_ref_t> q_cfg_rev_t;
    
  }

  namespace domain_impl {
    using namespace crab::cfg_impl;
    using namespace crab::domains; 
    using namespace ikos;
    
    typedef variable<ikos::z_number, varname_t> z_var;
    typedef linear_expression<ikos::z_number, varname_t> z_lin_t;
    typedef linear_constraint<ikos::z_number, varname_t> z_lin_cst_t;
    typedef linear_constraint_system<ikos::z_number, varname_t> z_lin_cst_sys_t;

    typedef variable< q_number, varname_t> q_var;
    typedef linear_expression<q_number, varname_t> q_lin_t;
    typedef linear_constraint<q_number, varname_t> q_lin_cst_t;
    typedef linear_constraint_system<q_number, varname_t> q_lin_cst_sys_t;
    
    typedef pointer_constraint <varname_t> ptr_cst_t;
    typedef interval<ikos::z_number> z_interval_t;
    typedef interval<q_number> q_interval_t;
    
    // Numerical domains over integers
    typedef interval_domain< ikos::z_number, varname_t > z_interval_domain_t;
    typedef numerical_congruence_domain< z_interval_domain_t> z_ric_domain_t;
    typedef SpDBM_impl::DefaultParams<ikos::z_number, SpDBM_impl::GraphRep::adapt_ss> z_SparseDBMGraph;
    typedef SparseDBM<ikos::z_number, varname_t, z_SparseDBMGraph> z_dbm_domain_t;
    typedef SDBM_impl::DefaultParams<ikos::z_number, SDBM_impl::GraphRep::adapt_ss> z_SplitDBMGraph;
    typedef SplitDBM<ikos::z_number, varname_t, z_SplitDBMGraph> z_sdbm_domain_t;
    typedef boxes_domain<ikos::z_number, varname_t > z_boxes_domain_t;
    typedef dis_interval_domain<ikos::z_number, varname_t > z_dis_interval_domain_t;
    typedef apron_domain<ikos::z_number, varname_t, apron_domain_id_t::APRON_INT > z_box_apron_domain_t;
    typedef apron_domain<ikos::z_number, varname_t, apron_domain_id_t::APRON_OCT > z_oct_apron_domain_t;
    typedef apron_domain<ikos::z_number, varname_t, apron_domain_id_t::APRON_OPT_OCT > z_opt_oct_apron_domain_t;
    typedef apron_domain<ikos::z_number, varname_t, apron_domain_id_t::APRON_PK > z_pk_apron_domain_t;
    typedef term_domain<term::TDomInfo<ikos::z_number, varname_t, z_interval_domain_t> > z_term_domain_t;
    typedef term_domain<term::TDomInfo<ikos::z_number, varname_t, z_sdbm_domain_t> > z_term_dbm_t;
    typedef term_domain<term::TDomInfo<ikos::z_number, varname_t, z_dis_interval_domain_t> > z_term_dis_int_t;
    typedef reduced_numerical_domain_product2<z_term_dis_int_t, z_sdbm_domain_t> z_num_domain_t;
    // Pointer domains over integers
    typedef nullity_domain<ikos::z_number, varname_t > z_nullity_domain_t;
    // Numerical x pointer domains over integers
    typedef numerical_nullity_domain<z_sdbm_domain_t> z_num_null_domain_t;
    // Boolean-numerical domain over integers
    typedef flat_boolean_numerical_domain<z_dbm_domain_t> z_bool_num_domain_t;
    
    /// Numerical over rationals
    typedef interval_domain<ikos::q_number,varname_t > q_interval_domain_t;
    typedef apron_domain<ikos::q_number,varname_t,apron_domain_id_t::APRON_PK>
    q_pk_apron_domain_t;
    typedef boxes_domain< q_number, varname_t > q_boxes_domain_t;         
  } 
}

namespace {

  bool stats_enabled = false;                                                                         

  #define SET_TEST_OPTIONS(ARGC,ARGV)                                                                 \
  boost::program_options::options_description po("Test Options");                                     \
  po.add_options()                                                                                    \
      ("log",  boost::program_options::value<std::vector<string> >(), "Enable specified log level");  \
  po.add_options()                                                                                    \
      ("stats",boost::program_options::bool_switch(&stats_enabled), "Enable stats");                  \
  boost::program_options::options_description cmmdline_options;                                       \
  cmmdline_options.add(po);                                                                           \
  boost::program_options::variables_map vm;                                                           \
  boost::program_options::positional_options_description p;                                           \
  boost::program_options::store(boost::program_options::command_line_parser(ARGC, ARGV).              \
            options(cmmdline_options).                                                                \
            positional(p).                                                                            \
              run(), vm);                                                                             \
  boost::program_options::notify(vm);                                                                 \
  if (vm.count("log")) {                                                                              \
    vector<string> loggers = vm ["log"].as<vector<string> > ();                                       \
    for(unsigned int i=0; i<loggers.size (); i++)                                                     \
      crab::CrabEnableLog (loggers [i]);                                                              \
  }                                                                                                   
  #endif 

  template<typename CFG, typename Dom, typename IntraFwdAnalyzer>
  void intra_run_impl (CFG* cfg, 
		       crab::cfg_impl::variable_factory_t& vfac, 
		       bool run_liveness,
		       unsigned widening, 
		       unsigned narrowing, 
		       unsigned jump_set_size){
    typedef crab::cfg::cfg_ref<CFG> cfg_ref_t;
    
    crab::analyzer::liveness<cfg_ref_t> *live = nullptr;
    if (run_liveness) {
      crab::analyzer::liveness<cfg_ref_t> live_(*cfg);
      live_.exec ();
      live=&live_;
    }
    // Run fixpoint
    Dom inv = Dom::top ();        
    crab::outs() << "Invariants using " << inv.getDomainName () << "\n";
    IntraFwdAnalyzer a (*cfg, inv, live, widening, narrowing, jump_set_size);
			
    a.run ();
    
    // Print invariants
    for (auto &b : *cfg) {
      auto inv = a[b.label ()];
      crab::outs() << crab::cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
    }
    crab::outs() << "\n";
    if (stats_enabled) {
      crab::CrabStats::Print(crab::outs());
      crab::CrabStats::reset();
    }
  }

  // To run abstract domains defined over integers
  template<typename Dom>
  void run (crab::cfg_impl::z_cfg_t* cfg, 
            crab::cfg_impl::variable_factory_t& vfac, 
            bool run_liveness,
            unsigned widening, 
            unsigned narrowing, 
            unsigned jump_set_size)
  {
    using namespace crab::analyzer;
    typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom> intra_fwd_analyzer_t;
			       			       
    intra_run_impl<crab::cfg_impl::z_cfg_t, Dom, intra_fwd_analyzer_t>
      (cfg, vfac, run_liveness, widening, narrowing, jump_set_size);
  }

  // To run abstract domains defined over rationals
  template<typename Dom>
  void run (crab::cfg_impl::q_cfg_t* cfg, 
	    crab::cfg_impl::variable_factory_t& vfac, 
	    bool run_liveness,
	    unsigned widening, 
	    unsigned narrowing, 
	    unsigned jump_set_size)
  {
    using namespace crab::analyzer;
    typedef intra_fwd_analyzer<crab::cfg_impl::q_cfg_ref_t,Dom> intra_fwd_analyzer_t;
    intra_run_impl<crab::cfg_impl::q_cfg_t, Dom, intra_fwd_analyzer_t>
      (cfg, vfac, run_liveness, widening, narrowing, jump_set_size);
  }

   
} //end namespace

