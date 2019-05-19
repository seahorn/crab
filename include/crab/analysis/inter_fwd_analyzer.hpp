#pragma once 

/* 
   A standard two-phase approach for summary-based,
   _context-insensitive_ inter-procedural analysis.
*/


#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/cfg/cfg.hpp>   // hasher of function declarations
#include <crab/cg/cg_bgl.hpp> // for sccg.hpp
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/inter_fwd_analyzer_ds.hpp>
#include <crab/analysis/dataflow/liveness.hpp>

#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>

#include <vector>
#include <set>
#include <string>

namespace crab {

namespace analyzer {

/**
 * Abstract transformer specialized for computing summaries.
 * class SumTable stores the summaries.
 **/
template<class SumTable> 
class bu_summ_abs_transformer final: 
    public intra_abs_transformer<typename SumTable::abs_domain_t> {
				      
public:

  typedef typename SumTable::abs_domain_t summ_abs_domain_t;
  typedef summ_abs_domain_t abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;    

private:

  typedef intra_abs_transformer<abs_dom_t> intra_abs_transform_t;
  typedef typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::variable_t variable_t;    
  typedef typename abs_dom_t::linear_expression_t linear_expression_t;    
    
  SumTable* m_sum_tbl;

public:

  using typename abs_transform_api_t::callsite_t;

  static void reuse_summary(abs_dom_t& caller, 
			    const callsite_t& cs,
			    const typename SumTable::Summary& summ) {
    // Before a summary is plug-in we rename it with unique variable
    // names so we avoid naming clashes in cases like for instance
    // summary variables have same names as lhs of callsites.
           
    CRAB_LOG("inter", 
	     crab::outs() << "    Reuse summary at " << cs << "\n";
	     crab::outs() << "    Summary:" << summ << "\n";); 

    std::set<variable_t> actuals, formals;
    // --- matching formal and actual parameters
    auto inputs = summ.get_renamed_inputs();
    unsigned i=0;
    // XXX: propagating down
    for (auto p : inputs) {
      auto a = cs.get_arg_name(i);
      if (!(a == p)) {
	CRAB_LOG("inter",
		 crab::outs() << "\t\tPropagate from caller to callee " 
		 << p << ":=" << a << "\n");
	inter_transformer_helpers<abs_dom_t>::unify(caller, p, a);
      }
      ++i;
      actuals.insert(a); formals.insert(p);
    }

    // --- meet caller's inv with summ
    auto sum_inv = summ.get_renamed_sum();
    caller = caller & sum_inv;
    CRAB_LOG("inter",
	     crab::outs() << "\t\tAfter meet caller and summary: "
	     <<  caller << "\n");

    // --- matching callsite's lhs and callee's return value 
    // XXX: propagate from the return values in the callee to the
    // lhs variables of the callsite in the caller.
    auto const &caller_vts = cs.get_lhs();
    auto const &callee_outs = summ.get_renamed_outputs();
    assert(caller_vts.size() == callee_outs.size());
	
    auto caller_it = caller_vts.begin();
    auto caller_et = caller_vts.end();
    auto callee_it = callee_outs.begin();
      
    // XXX: propagating up
    for (; caller_it != caller_et; ++caller_it, ++callee_it){
      auto vt = *caller_it;
      auto r = *callee_it;
      CRAB_LOG("inter",
	       crab::outs() << "\t\tPropagate from callee to caller " 
	       << vt << ":=" << r << "\n");
      inter_transformer_helpers<abs_dom_t>::unify(caller, vt, r);
      actuals.insert(vt); formals.insert(r);
    }
      
    // --- remove from caller only formal parameters so we can keep
    //     as much context from the caller as possible
    std::vector<variable_t> vs;
    boost::set_difference(formals, actuals, std::back_inserter(vs));
    caller.forget(vs);
      
    CRAB_LOG("inter", 
	     crab::outs() << "\t\tAfter forgetting formal parameters {";
	     for(auto v: vs) crab::outs() << v << ";";
	     crab::outs() << "}=" <<  caller << "\n");
  }

public:
    
  bu_summ_abs_transformer(abs_dom_t* inv, SumTable* sum_tbl)
    : intra_abs_transform_t(inv)
    , m_sum_tbl(sum_tbl) { }

  ~bu_summ_abs_transformer() = default;
    
  virtual void exec(callsite_t& cs) override {
    if (!m_sum_tbl) {
      CRAB_WARN("The summary table is empty: ignored analysis of callsite");
    } else {
      if (m_sum_tbl->hasSummary(cs)) {
	auto summ = m_sum_tbl->get(cs);
	reuse_summary(*(this->m_inv), cs, summ);
      }
      else {
	CRAB_LOG("inter",
		 crab::outs() << "\tSummary not found for " << cs << "\n");
	for (auto vt: cs.get_lhs()) {
	  *(this->m_inv) -= vt;  // havoc
	}
      }
    }
  }
}; 

/// Conversion between domains
template<typename Domain>
inline void convert_domains(Domain from, Domain& to) { to = from; }

template<typename Domain1, typename Domain2>
inline void convert_domains(Domain1 from, Domain2& to) {
  CRAB_LOG("inter",
	   crab::outs() << "Converting from "
	                << Domain1::getDomainName() << " to "
	                << Domain2::getDomainName() << " might lose precision if "
	                << Domain1::getDomainName() << " is more precise than "
	                << Domain2::getDomainName() << "\n");
	     
  for (auto cst : from.to_linear_constraint_system()) {
    to += cst;
  }
}

/**
 * Abstract transformer specialized for performing top-down forward
 * traversal while reusing summaries at the callsites.
 * class CallCtxTable stores the calling context.
 **/
template<class SumTable, class CallCtxTable> 
class td_summ_abs_transformer final: 
    public intra_abs_transformer<typename CallCtxTable::abs_domain_t> {

public:

  typedef typename SumTable::abs_domain_t     summ_abs_domain_t;
  typedef typename CallCtxTable::abs_domain_t call_abs_domain_t;
  typedef call_abs_domain_t abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;

private:
    
  typedef intra_abs_transformer<abs_dom_t> intra_abs_transform_t;
    
public:
    
  typedef typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
  using typename abs_transform_api_t::callsite_t;

private:
    
  SumTable* m_sum_tbl;
  CallCtxTable* m_call_tbl;
        
public:
    
  td_summ_abs_transformer(abs_dom_t* inv,
			  SumTable* sum_tbl, CallCtxTable* call_tbl)
    : intra_abs_transform_t(inv), 
      m_sum_tbl(sum_tbl),
      m_call_tbl(call_tbl) { }

  ~td_summ_abs_transformer() = default;
    
  virtual void exec(callsite_t& cs) override {
    if (!m_sum_tbl) {
      CRAB_WARN("The summary table is empty");
    } else if(m_sum_tbl->hasSummary(cs)) {
      auto summ = m_sum_tbl->get(cs);
        
      CRAB_LOG("inter", 
	       crab::outs() << "    Pluging caller context into callee\n"
	       << "    Summary: " << summ << "\n");
	
      ///////
      /// Generate the callee context and store it.
      ///////
      typedef typename abs_dom_t::variable_t variable_t;
	
      abs_dom_t callee_ctx_inv(*(this->m_inv));
      // --- matching formal and actual parameters
      // XXX: propagating down 	
      unsigned i=0;
      const std::vector<variable_t>& inputs = summ.get_inputs();
      for (variable_t p : inputs) {
	variable_t a = cs.get_arg_name(i);
	if (!(a == p)) {
	  inter_transformer_helpers<abs_dom_t>::unify(callee_ctx_inv, p, a);
	}
	++i;
      }

      // --- project only onto formal parameters
      callee_ctx_inv.project(inputs);
      // --- store the callee context
      CRAB_LOG("inter", 
	       crab::outs() << "\t\tCallee context stored: " 
	       << callee_ctx_inv << "\n");
	
      m_call_tbl->insert(cs, callee_ctx_inv);          
        
      /////
      // Generate the continuation at the caller
      /////
        
      // --- convert this->m_inv to the language of summ_abs_dom_t(summ)
      summ_abs_domain_t caller_ctx_inv;
      convert_domains(*(this->m_inv), caller_ctx_inv);
      // forget the variables of the lhs of the callsite, otherwise
      // caller_ctx_inv and m_inv may be inconsistent if the lhs
      // variables are constrained before the callsite (e.g., x=-5;
      // x := abs(x);)
      for (auto vt: cs.get_lhs()) {
	*(this->m_inv) -= vt;
      }
	
      CRAB_LOG("inter",
	       crab::outs() << "\t\tCaller context: " <<  caller_ctx_inv << "\n");
        
      // --- reuse summary to do the continuation
      bu_summ_abs_transformer<SumTable>::reuse_summary(caller_ctx_inv, cs, summ);
      CRAB_LOG("inter",
	       crab::outs() << "\t\tCaller context after plugin summary: " 
	       << caller_ctx_inv << "\n");
	
      // --- convert back inv to the language of abs_dom_t
      convert_domains(caller_ctx_inv, *(this->m_inv));
      CRAB_LOG("inter",
	       crab::outs() << "\t\tCaller continuation after " 
	       <<  cs << "=" <<  *(this->m_inv) << "\n");
      return;
    } else {
      // --- no summary found: do nothing
      if (CrabWarningFlag) {
	crab::errs() << "CRAB WARNING: no summary found for callsite " << cs << "\n";
	m_sum_tbl->write(crab::errs());
      }
    }
        
    // We could not reuse a summary so we just havoc lhs of the call
    // site.      
    for(auto vt: cs.get_lhs()) {
      *(this->m_inv) -= vt;
    }
  }
}; 

template<typename CallGraph, 
	 // abstract domain used for the bottom-up phase
	 typename BU_Dom, 
	 // abstract domain used for the top-down phase
	 typename TD_Dom>
class inter_fwd_analyzer: public boost::noncopyable {

  typedef typename CallGraph::node_t cg_node_t;
  typedef typename CallGraph::edge_t cg_edge_t;

public:

  typedef CallGraph cg_t;
  typedef typename cg_node_t::cfg_t cfg_t; 
  typedef typename cfg_t::varname_t varname_t;
  typedef typename cfg_t::number_t number_t;
  typedef typename cfg_t::variable_t variable_t;
  typedef liveness<cfg_t> liveness_t;     
  typedef boost::unordered_map <cfg_t, const liveness_t*> liveness_map_t;

private:

  typedef summary_table <cfg_t, BU_Dom> summ_tbl_t;
  typedef call_ctx_table <cfg_t, TD_Dom> call_tbl_t;

public:

  typedef bu_summ_abs_transformer<summ_tbl_t> bu_abs_tr;
  typedef td_summ_abs_transformer<summ_tbl_t, call_tbl_t> td_abs_tr;
  typedef fwd_analyzer<cfg_t, bu_abs_tr> bu_analyzer;
  typedef fwd_analyzer<cfg_t, td_abs_tr> td_analyzer;
  typedef typename summ_tbl_t::Summary summary_t;
  typedef boost::shared_ptr<summary_t> summary_ptr;
  // for checkers
  typedef TD_Dom abs_dom_t;
  typedef td_abs_tr abs_tr_t;
      
private:

  typedef boost::shared_ptr<td_analyzer> td_analyzer_ptr;
  typedef boost::unordered_map<std::size_t, td_analyzer_ptr> invariant_map_t;

  CallGraph m_cg;
  const liveness_map_t* m_live;
  invariant_map_t m_inv_map;
  summ_tbl_t m_summ_tbl;
  call_tbl_t m_call_tbl;
  unsigned int m_widening_delay;
  unsigned int m_descending_iters;
  size_t m_jump_set_size; // max size of the jump set (=0 if jump set disabled)
      
  const liveness_t* get_live(const cfg_t& c) {
    if (m_live) {
      auto it = m_live->find(c);
      if (it != m_live->end())
	return it->second;
    }
    return nullptr;
  }
      
public:
      
  inter_fwd_analyzer(CallGraph cg, const liveness_map_t* live,
		     // fixpoint parameters
		     unsigned int widening_delay=1,
		     unsigned int descending_iters=UINT_MAX,
		     size_t jump_set_size=0)
    : m_cg(cg), m_live(live),
      m_widening_delay(widening_delay), 
      m_descending_iters(descending_iters),
      m_jump_set_size(jump_set_size) {
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking call graph ... ";);
    crab::CrabStats::resume("CallGraph type checking");
    cg.type_check();
    crab::CrabStats::stop("CallGraph type checking");	
    CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);
  }
      
  void run(TD_Dom init = TD_Dom::top())  {

    CRAB_VERBOSE_IF(1, get_msg_stream() << "Started inter-procedural analysis\n";);
    CRAB_LOG("inter", 
	     m_cg.write(crab::outs()); crab::outs() << "\n");
                 
    crab::ScopedCrabStats __st__("Inter");

    bool has_noedges = true;
    for (auto const &n: boost::make_iterator_range(m_cg.nodes())) {
      if (m_cg.num_succs(n) > 0) {
	has_noedges = false;
	break;
      }
    } 
        
    if (has_noedges) {
      // -- Special case when the program has been inlined.
      CRAB_LOG("inter",
	       crab::outs() << "Callgraph has no edges so no summaries are computed.\n");

      CRAB_LOG("inter", 
	       m_cg.write(crab::outs());
	       crab::outs() << "\n";);
                   

      for (auto &v: boost::make_iterator_range(m_cg.nodes())) {
	crab::ScopedCrabStats __st__("Inter.TopDown");

	auto cfg = v.get_cfg();
	assert(cfg.has_func_decl());
	auto fdecl = cfg.get_func_decl();
	std::string fun_name = fdecl.get_func_name();
	if (fun_name != "main") continue;
            
	CRAB_LOG("inter",
		 crab::outs() << "++ Analyzing function "
		 << fdecl.get_func_name() << "\n");

	auto abs_tr = boost::make_shared<td_abs_tr>(&init, &m_summ_tbl, &m_call_tbl);
	auto a = boost::make_shared<td_analyzer>(cfg, nullptr, &*abs_tr, get_live(cfg),
						 m_widening_delay,
						 m_descending_iters,
						 m_jump_set_size);
						      
	a->run_forward();
	m_inv_map.insert({crab::cfg::cfg_hasher<cfg_t>::hash(fdecl), a});
      }
      return;
    }

    // -- General case 
    std::vector<cg_node_t> rev_order;
    graph_algo::scc_graph<CallGraph> Scc_g(m_cg);
    graph_algo::rev_topo_sort(Scc_g, rev_order);

    CRAB_VERBOSE_IF(1,get_msg_stream() << "== Bottom-up phase ...\n";);	
    for (auto n: rev_order) {
      crab::ScopedCrabStats __st__("Inter.BottomUp");
      std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members(n);
      for (auto m: scc_mems) {

	auto cfg = m.get_cfg();
	assert(cfg.has_func_decl());
	auto fdecl = cfg.get_func_decl();            
	std::string fun_name = fdecl.get_func_name();
	if (fun_name != "main" && cfg.has_exit()) {
	  CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Analyzing function "
			  << fdecl.get_func_name() << "\n";);
	  // --- run the analysis
	  auto init_inv = BU_Dom::top();
	  bu_abs_tr abs_tr(&init_inv, &m_summ_tbl);
	  bu_analyzer a(cfg, nullptr, &abs_tr, get_live(cfg),
			m_widening_delay, m_descending_iters, m_jump_set_size);
	  a.run_forward();
	      
	  // --- build the summary
	  std::vector<variable_t> formals, inputs, outputs;
	  formals.reserve(fdecl.get_num_inputs() + fdecl.get_num_outputs());
	  inputs.reserve(fdecl.get_num_inputs());
	  outputs.reserve(fdecl.get_num_outputs());

	  for (unsigned i=0; i < fdecl.get_num_inputs();i++) {
	    inputs.push_back(fdecl.get_input_name(i));
	    formals.push_back(fdecl.get_input_name(i));
	  }
	  for (unsigned i=0; i < fdecl.get_num_outputs();i++) {
	    outputs.push_back(fdecl.get_output_name(i));
	    formals.push_back(fdecl.get_output_name(i));
	  }
	      
	  // --- project onto formal parameters and return values
	  auto inv = a.get_post(cfg.exit());
	  //crab::CrabStats::count(BU_Dom::getDomainName() + ".count.project");
	  inv.project(formals);
	  m_summ_tbl.insert(fdecl, inv, inputs, outputs);
	}
      }
    } 

    CRAB_VERBOSE_IF(1, get_msg_stream() << "== Top-down phase ...\n";);
    bool is_root = true;
    for (auto n: boost::make_iterator_range(rev_order.rbegin(),
					    rev_order.rend())) {
      crab::ScopedCrabStats __st__("Inter.TopDown");
      std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members(n);
	  
      // The SCC is recursive if it has more than one element or
      // there is only one that calls directly to itself.
      bool is_recursive = (scc_mems.size() > 1) ||
	std::any_of(m_cg.succs(n).first, m_cg.succs(n).second,
		    [n](const cg_edge_t& e) {
		      return (n == e.dest());
		    });

      for (auto m: scc_mems) {
	auto cfg = m.get_cfg();
	assert(cfg.has_func_decl());
	auto fdecl = cfg.get_func_decl();
	CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Analyzing function " 
			<< fdecl.get_func_name() << "\n";);
	if (is_recursive) {
	  // If the SCC is recursive then what we have in
	  // m_call_tbl is incomplete and therefore it is unsound
	  // to use it. To remedy it, we insert another calling
	  // context with top value that approximates all the
	  // possible calling contexts during the recursive calls.
	  m_call_tbl.insert(fdecl, TD_Dom::top());
	}
	   
	auto init_inv = init;	    
	if (is_root) {
	  is_root = false;
	} else {
	  init_inv = m_call_tbl.get_call_ctx(fdecl);
	}
	    
	CRAB_LOG("inter",
		 crab::outs() << "    Starting analysis of "
		 << fdecl <<  " with " << init_inv << "\n");

	auto abs_tr = boost::make_shared<td_abs_tr>(&init_inv, &m_summ_tbl, &m_call_tbl);
	auto a = boost::make_shared<td_analyzer>(cfg, nullptr, &*abs_tr, get_live(cfg),
						 m_widening_delay,
						 m_descending_iters,
						 m_jump_set_size);
	a->run_forward();
	m_inv_map.insert(std::make_pair(crab::cfg::cfg_hasher<cfg_t>::hash(fdecl), a));
      }
    }
    CRAB_VERBOSE_IF(1,get_msg_stream() << "Finished inter-procedural analysis\n";);	
  }

  //! return the analyzed call graph
  CallGraph& get_call_graph() {
    return m_cg;
  }

  //! Return the invariants that hold at the entry of b in cfg
  TD_Dom get_pre(const cfg_t &cfg, typename cfg_t::basic_block_label_t b) const {
    assert(cfg.has_func_decl());
    auto fdecl = cfg.get_func_decl();
    auto const it = m_inv_map.find(crab::cfg::cfg_hasher<cfg_t>::hash(fdecl));
    if (it != m_inv_map.end()) {
      return it->second->get_pre(b);
    } else {
      return TD_Dom::top();
    }
  }
      
  //! Return the invariants that hold at the exit of b in cfg
  TD_Dom get_post(const cfg_t &cfg, typename cfg_t::basic_block_label_t b) const {
    assert(cfg.has_func_decl());
    auto fdecl = cfg.get_func_decl();
    auto const it = m_inv_map.find(crab::cfg::cfg_hasher<cfg_t>::hash(fdecl));
    if (it != m_inv_map.end()) {
      return it->second->get_post(b);
    } else {
      return TD_Dom::top();
    }
  }

  // clear all invariants
  void clear() {
    m_inv_map.clear();
  }
      
  //! Propagate inv through statements
  boost::shared_ptr<abs_tr_t> get_abs_transformer(TD_Dom* inv) {
    return boost::make_shared<abs_tr_t>(inv, &m_summ_tbl, &m_call_tbl);
  }

  //! Return true if there is a summary for a function
  bool has_summary(const cfg_t &cfg) const {
    assert(cfg.has_func_decl());
    auto fdecl = cfg.get_func_decl();
    return m_summ_tbl.hasSummary(fdecl);
  }

  //! Return the summary for a function
  summary_ptr get_summary(const cfg_t &cfg) const {
    assert(cfg.has_func_decl());
    auto fdecl = cfg.get_func_decl();
    if (m_summ_tbl.hasSummary(fdecl)) {
      summary_t summ = m_summ_tbl.get(fdecl);
      return boost::make_shared<summary_t>(summ);
    } else {
      CRAB_WARN("Summary not found");
      return nullptr;
    }
  }
        
}; 

} // end namespace
} // end namespace

