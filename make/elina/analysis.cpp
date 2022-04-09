#include <crab/config.h>
#include <crab/types/linear_constraints.hpp>
#include <crab/types/varname_factory.hpp>
#include <crab/types/variable.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

// Abstract domains
#include <crab/domains/intervals.hpp>
#include <crab/domains/elina_domains.hpp>
#include <crab/domains/generic_abstract_domain.hpp>

// Analysis
#include <crab/analysis/fwd_analyzer.hpp>

using namespace crab;
using namespace ikos;

// A variable factory based on strings
using variable_factory_t = crab::var_factory_impl::str_variable_factory;
using varname_t = typename variable_factory_t::varname_t;
namespace crab {
template<>
class variable_name_traits<std::string> {
public:
  static std::string to_string(std::string varname) {
    return varname;
  }
};
} // end namespace crab
// CFG basic block labels
using basic_block_label_t = std::string;
/// To define CFG over integers
using cfg_t = crab::cfg::cfg<basic_block_label_t, varname_t, z_number>;
using cfg_ref_t = crab::cfg::cfg_ref<cfg_t>;
using basic_block_t = cfg_t::basic_block_t;
using var_t = crab::variable<z_number, varname_t>;
using lin_exp_t = linear_expression<z_number, varname_t>;
using lin_cst_t = linear_constraint<z_number, varname_t> ;
using lin_cst_sys_t = linear_constraint_system<z_number, varname_t> ;
namespace crab {
template<>
class basic_block_traits<basic_block_t> {
public:
  using bb_label_t = typename basic_block_t::basic_block_label_t;  
  static std::string to_string(const bb_label_t &bbl) {
    return bbl;
  }
};
} // end namespace crab  
//// To define CFG over integers
using cfg_t = cfg::cfg<basic_block_label_t, varname_t, z_number>;
using cfg_ref_t = cfg::cfg_ref<cfg_t>;
using basic_block_t = cfg_t::basic_block_t ;

///////// Begin Crab Abstract Domains /////////////
using interval_domain_t = interval_domain<z_number,varname_t>;
using oct_domain_t = domains::elina_domain<z_number,varname_t,
					  domains::elina_domain_id_t::ELINA_OCT>;
using pk_domain_t = domains::elina_domain<z_number,varname_t,
					  domains::elina_domain_id_t::ELINA_PK>;
using abs_domain_t = domains::abstract_domain<var_t>;
///////// End Crab Abstract Domains /////////////

///////// Begin Analyses /////////////
using fwd_analyzer_t = analyzer::intra_fwd_analyzer<cfg_ref_t, abs_domain_t>;
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
  basic_block_t& body = prog.insert("body");
  basic_block_t& body_if_tt = prog.insert("body_if_tt");
  basic_block_t& body_if_ff = prog.insert("body_if_ff");
  basic_block_t& body_tail = prog.insert("body_tail");

  entry.add_succ(header);
  header.add_succ(body);
  header.add_succ(exit);  
  body.add_succ(body_if_tt);
  body.add_succ(body_if_ff);
  body_if_tt.add_succ(body_tail);
  body_if_ff.add_succ(body_tail);
  body_tail.add_succ(header);

  entry.assign(x, 0);
  entry.assign(i, 0);
  body.assume(lin_cst_t(i <= 9999));
  exit.assume(lin_cst_t(i >= 10000));
  body_if_tt.add(x, x, 1);
  body_if_ff.add(x, x, 3);
  body_tail.add(i, i, 1);

  outs() << prog << "\n";

  auto print_invariants = [&prog](const fwd_analyzer_t &analyzer) {
			   for (auto &bb : prog) {
			     std::string bb_name = bb.label();
			     auto inv = analyzer.get_pre(bb_name);
			     outs () << bb_name << ":" << inv << "\n";
			   }
			 };
  interval_domain_t top_intv;
  abs_domain_t init(top_intv);
  fwd_analyzer_t analyzer(prog, init);
  analyzer.run();
  outs () << "Invariants using intervals\n";
  print_invariants(analyzer);
  analyzer.clear();

  oct_domain_t top_oct;
  analyzer.get_abs_transformer().set_abs_value(top_oct);
  analyzer.run();
  outs () << "Invariants using Elina Octagons\n";
  print_invariants(analyzer);
  analyzer.clear();
  
  pk_domain_t top_pk;
  analyzer.get_abs_transformer().set_abs_value(top_pk);
  analyzer.run();
  outs () << "Invariants using Elina Polyhedra\n";
  print_invariants(analyzer);
  
  return 0;
}
