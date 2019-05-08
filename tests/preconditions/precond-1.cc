#include "../program_options.hpp"
#include "../common.hpp"
#include <crab/analysis/bwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* CheckedC example: args.c */
z_cfg_t* prog (variable_factory_t &vfac)  {

  // Defining program variables
  z_var argc(vfac["argc"], crab::INT_TYPE, 32);
  z_var len(vfac["len"], crab::INT_TYPE, 32);
  
  // entry and exit block
  auto cfg = new z_cfg_t("bb0","bb6");
  /*
    if (argc > 2) 
       NumNodes = atoi(argv[2]);
    else 
       NumNodes = 4;
    
    if (argc > 1)
      nbody = atoi(argv[1]);
    else
      nbody = 32;
  */
  
  // adding blocks
  z_basic_block_t& bb0   = cfg->insert ("bb0");
  z_basic_block_t& bb1   = cfg->insert ("bb1");
  z_basic_block_t& bb2   = cfg->insert ("bb2");
  z_basic_block_t& bb3   = cfg->insert ("bb3");
  z_basic_block_t& bb4   = cfg->insert ("bb4");
  z_basic_block_t& bb5   = cfg->insert ("bb5");
  z_basic_block_t& bb6   = cfg->insert ("bb6");
  
  bb0 >> bb1;
  bb0 >> bb2;
  bb1 >> bb3;
  bb2 >> bb3;
  
  bb3 >> bb4;
  bb3 >> bb5;
  bb4 >> bb6;
  bb5 >> bb6;
  
  bb1.assume (argc >= 3);
  bb1.assertion(2 <= len - 1);
  
  bb4.assume (argc >= 2);
  bb4.assertion(1 <= len - 1);
  
  return cfg;
}


int main(int argc, char** argv) {
#if (defined(HAVE_APRON) ||defined(HAVE_ELINA))
  
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  z_cfg_t* cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  {
#ifdef HAVE_APRON
    typedef crab::analyzer::necessary_preconditions_fixpoint_iterator
    <z_cfg_ref_t,z_pk_apron_domain_t> analysis_t;  
    z_pk_apron_domain_t final_states = z_pk_apron_domain_t::bottom();
#endif
#ifdef HAVE_ELINA
    typedef crab::analyzer::necessary_preconditions_fixpoint_iterator
      <z_cfg_ref_t,z_pk_elina_domain_t> analysis_t;  
    z_pk_elina_domain_t final_states = z_pk_elina_domain_t::bottom();
#endif
    analysis_t analyzer(*cfg, nullptr, final_states, false /*error states*/);
    analyzer.run_backward();    
    crab::outs () << "Necessary preconditions from error states using Polyhedra:\n";
    // Print preconditions
    for (z_basic_block_t& bb : *cfg) {
      std::string bb_name = bb.label();
      auto inv = analyzer[bb_name];
      crab::outs () << bb_name << ":" << inv << "\n";
    }  
  }
  {
#ifdef HAVE_APRON
    typedef crab::analyzer::necessary_preconditions_fixpoint_iterator
    <z_cfg_ref_t,z_pk_apron_domain_t> analysis_t;  
    z_pk_apron_domain_t final_states = z_pk_apron_domain_t::top();
#endif
#ifdef HAVE_ELINA
    typedef crab::analyzer::necessary_preconditions_fixpoint_iterator
      <z_cfg_ref_t,z_pk_elina_domain_t> analysis_t;  
    z_pk_elina_domain_t final_states = z_pk_elina_domain_t::top();
#endif
    analysis_t analyzer(*cfg, nullptr, final_states, true /*good states*/);
    analyzer.run_backward();    
    crab::outs () << "Necessary preconditions from safe states using Polyhedra:\n";
    // Print preconditions
    for (z_basic_block_t& bb : *cfg) {
      std::string bb_name = bb.label();
      auto inv = analyzer[bb_name];
      crab::outs () << bb_name << ":" << inv << "\n";
    }  
  }    
  // free the CFG
  delete cfg;
#endif
  return 0;
}
