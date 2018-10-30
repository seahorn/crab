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

// Analysis
#include <crab/analysis/fwd_analyzer.hpp>

using namespace crab;
using namespace crab::domains;
using namespace crab::cfg;
using namespace crab::analyzer;
using namespace ikos;

///////// Begin Crab CFG /////////////
// A variable factory based on strings
typedef var_factory_impl::str_variable_factory variable_factory_t;
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

///////// Begin Analyses /////////////
typedef intra_fwd_analyzer<cfg_ref_t, interval_domain_t> interval_analysis_t;
typedef intra_fwd_analyzer<cfg_ref_t, pk_domain_t> pk_analysis_t;
///////// End Analyses /////////////
int main(int argc, char**argv) {

  variable_factory_t vfac;
  var_t x(vfac["x"], INT_TYPE, 32);
  var_t i(vfac["i"], INT_TYPE, 32);

  // Build a CFG
  cfg_t prog("entry", "exit");
  
  basic_block_t& entry = prog.insert("entry");
  basic_block_t& header = prog.insert("header");
  basic_block_t& exit = prog.insert("exit");
  basic_block_t& body_if = prog.insert("body_if");
  basic_block_t& body_if_tt = prog.insert("body_if_tt");
  basic_block_t& body_if_ff = prog.insert("body_if_ff");
  basic_block_t& body_tail = prog.insert("body_tail");

  entry >> header;
  entry >> exit;
  header >> body_if;
  body_if >> body_if_tt;
  body_if >> body_if_ff;
  body_if_tt >> body_tail;
  body_if_ff >> body_tail;
  body_tail >> header;

  entry.assign(x, 0);
  entry.assign(i, 0);
  body_if.assume(lin_cst_t(i <= 99));
  exit.assume(lin_cst_t(i >= 100));
  body_if_tt.add(x, x, 1);
  body_if_ff.add(x, x, 2);
  body_tail.add(i, i, 1);

  outs() << prog << "\n";

  {
    // Analysis of the CFG
    interval_analysis_t ianalyzer(prog);
    ianalyzer.run();
    
    // Print invariants
    outs () << "Invariants using intervals\n";
    for (auto &bb : prog) {
      std::string bb_name = bb.label();
      auto inv = ianalyzer.get_pre(bb_name);
      outs () << bb_name << ":" << inv << "\n";
    }
  }

  {
    // Analysis of the CFG
    pk_analysis_t ianalyzer(prog);
    ianalyzer.run();

    // Print invariants
    outs () << "Invariants using polyhedra\n";        
    for (auto &bb : prog) {
      std::string bb_name = bb.label();
      auto inv = ianalyzer.get_pre(bb_name);
      outs () << bb_name << ":" << inv << "\n";
    }
  }
  
  return 0;
}
