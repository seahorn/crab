#pragma once

/*
 * Standard top-down interprocedural analysis with memoization.
 *
 * The analysis is defined recursively. This is needed so that the
 * intra-procedural analysis can be called as a black-box without any
 * change. However, it might run out of stack on very deep call graphs
 * so special care must be taken.
 *
 * **Important**: the analysis assumes that the Crab CFG ensures that
 * an input parameter has only one use and that use is the right-side
 * operand of an assignment locate at the entry block of the
 * function. That is, that the following transformation took place
 * before the analysis starts:
 * 
 * foo(...,i1,...) {      foo(...,i1',...) {
 *                          // where := can be assign/expand/copy_region
 *                          i1 := i1';  
 *   BODY          ==>      BODY
 * }                      }
 */

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/inter/inter_analyzer_api.hpp>
#include <crab/analysis/inter/inter_params.hpp>
#include <crab/cg/cg_bgl.hpp> // for wto of callgraphs
#include <crab/fixpoint/wto.hpp>
#include <crab/domains/inter_abstract_operations_callsite_info.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>

#include <algorithm> // sort, set_difference
#include <climits>
#include <deque>
#include <memory> // unique_ptr
#include <set>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace crab {
namespace analyzer {

#define TimerInter "Inter"  
#define TimerInterCheckCache "Inter.lookup_summary"  
#define TimerInterStoreSum "Inter.store_summary"
#define TimerInterAnalyzeFunc "Inter.analyze_function"
const std::string CounterRecursiveCalls("Inter.recursive_callsites");
const std::string CounterReusedCalls("Inter.reused_callsites");
const std::string CounterAnalyzedCalls("Inter.analyzed_callsites");
const std::string CounterCalls("Inter.callsites");

#define TD_INTER_TIMER_STATS(NAME) CRAB_SCOPED_TIMER_STATS(NAME, 0)
#define TD_INTER_COUNT_STATS(NAME) CRAB_COUNT_STATS(NAME, 0)
  
namespace top_down_inter_impl {

/**
 *  This class represents the calling context of a function F. A
 *  calling context C consists of a summary (pair of pre/post
 *  conditions) that summarizes the analysis of F for a particular
 *  callsite. The calling context can also store optionally all the
 *  invariants inferred during the analysis of F.
 *
 *  Whether we keep the invariants for all F's block is an analysis
 *  parameter. For some clients, this information might be useful
 *  although it won't probably scale to keep them all in memory.
 **/
template <typename CFG, typename AbsDom, typename InvariantMap>
class calling_context {
public:
  using cfg_t = CFG;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;
  using fdecl_t = typename cfg_t::fdecl_t;
  using abs_dom_t = AbsDom;
  using invariant_map_t = InvariantMap;

private:
  static_assert(
      std::is_same<basic_block_label_t, typename InvariantMap::key_type>::value,
      "The key of an invariant map should be basic block label");
  static_assert(
      std::is_same<abs_dom_t, typename InvariantMap::mapped_type>::value,
      "The mapped value of an invariant map should be an abstract domain");

  using calling_context_t = calling_context<CFG, abs_dom_t, invariant_map_t>;

  // -- the summary: this is an abstraction of the call context

  // A summary is a pair (pre, post) of abstract states. pre's
  // variables only contains the input parameters of the function. pre
  // is used to check if a summary can be reused. post's variables
  // should contain both input and output parameters of a function. We
  // also include the inputs in post's variables so that the
  // underlying abstract domain can keep track of relationships
  // between input and outputs. Therefore, post is the actual
  // input-output relationship that summarizes the behavior of the
  // function.

  const fdecl_t &m_fdecl;
  abs_dom_t m_pre_summary;  // must be projected on m_fdecl inputs
  abs_dom_t m_post_summary; // must be projected on m_fdecl inputs and outputs

  // invariants that hold at the entry of each function's block
  invariant_map_t m_pre_invariants;
  // invariants that hold at the exit of each function's block
  invariant_map_t m_post_invariants;

  // if true then the calling context has not been joined yet with
  // other contexts.
  bool m_exact;
  // whether to keep all the invariants. Useful for printing
  // context-sensitive invariants but very expensive.
  bool m_keep_invariants;

  inline abs_dom_t make_bottom() const {
    return m_pre_summary.make_bottom();
  }

  static invariant_map_t join(invariant_map_t &m1, invariant_map_t &m2) {
    invariant_map_t out;
    for (auto &kv : m2) {
      auto it = m1.find(kv.first);
      if (it != m1.end()) {
        out.insert({kv.first, it->second | kv.second});
      } else {
        out.insert({kv.first, kv.second});
      }
    }
    return out;
  }

  // Private constructor used to join calling contexts while keeping
  // all invariants
  calling_context(const fdecl_t &fdecl,
		  abs_dom_t pre_summary, abs_dom_t post_summary,
                  invariant_map_t &&pre_invariants,
                  invariant_map_t &&post_invariants)
      : m_fdecl(fdecl), m_pre_summary(pre_summary), m_post_summary(post_summary),
        m_pre_invariants(std::move(pre_invariants)),
        m_post_invariants(std::move(post_invariants)), m_exact(false),
        m_keep_invariants(true) {}

  // Private constructor used to join calling contexts but without
  // keeping invariants
  calling_context(const fdecl_t &fdecl, abs_dom_t pre_summary,
		  abs_dom_t post_summary)
      : m_fdecl(fdecl), m_pre_summary(pre_summary),
	m_post_summary(post_summary), m_exact(false),
        m_keep_invariants(false) {}

public:
  calling_context(const fdecl_t &fdecl, abs_dom_t pre_summary, abs_dom_t post_summary,
                  bool keep_invariants, invariant_map_t &&pre_invariants,
                  invariant_map_t &&post_invariants)
      : m_fdecl(fdecl), m_pre_summary(pre_summary), m_post_summary(post_summary),
        m_pre_invariants(std::move(pre_invariants)),
        m_post_invariants(std::move(post_invariants)), m_exact(true),
        m_keep_invariants(keep_invariants) {

    if (!m_keep_invariants) {
      m_pre_invariants.clear();
      m_post_invariants.clear();
    }
  }

  calling_context(const calling_context_t &o) = delete;

  calling_context_t &operator=(const calling_context_t &o) = delete;

  const fdecl_t &get_fdecl() const { return m_fdecl; }

  const abs_dom_t &get_pre_summary() const { return m_pre_summary; }

  const abs_dom_t &get_post_summary() const { return m_post_summary; }

  // Check if d entails the summary precondition
  bool is_subsumed(const abs_dom_t &d, bool exact_check) const {
    const abs_dom_t &pre_summary = get_pre_summary();
    if (m_exact && exact_check) {
      return (d <= pre_summary && pre_summary <= d);
    } else {
      return (d <= pre_summary);
    }
  }

  // Join this and other
  std::unique_ptr<calling_context_t> join_with(calling_context_t &other) {
    if (!(get_fdecl() == other.get_fdecl())) {
      CRAB_ERROR("Cannot join calling contexts because they are from different "
                 "functions");
    }

    if (!m_keep_invariants) {
      return std::unique_ptr<calling_context_t>(new calling_context_t(
          m_fdecl, m_pre_summary | other.get_pre_summary(),
	  m_post_summary | other.get_post_summary()));

    } else {
      return std::unique_ptr<calling_context_t>(new calling_context_t(
          m_fdecl, m_pre_summary | other.get_pre_summary(),
	  m_post_summary | other.get_post_summary(),
          std::move(join(m_pre_invariants, other.m_pre_invariants)),
          std::move(join(m_post_invariants, other.m_post_invariants))));
    }
  }

  // invariants that hold at the entry of basic block bb
  abs_dom_t get_pre(const basic_block_label_t &bb) const {
    if (!m_keep_invariants) {
      abs_dom_t top;
      return top;
    } else {
      auto it = m_pre_invariants.find(bb);
      if (it == m_pre_invariants.end()) {
        return make_bottom(); // dead block under particular context
      } else {
        return it->second;
      }
    }
  }

  // invariants that hold at the exit of basic block bb
  abs_dom_t get_post(const basic_block_label_t &bb) const {
    if (!m_keep_invariants) {
      abs_dom_t top;
      return top;
    } else {
      auto it = m_post_invariants.find(bb);
      if (it == m_post_invariants.end()) {
        return make_bottom(); // dead block under particular context
      } else {
        return it->second;
      }
    }
  }
  
  void write(crab_os &o) const {
    const abs_dom_t &pre_summary = get_pre_summary();
    const abs_dom_t &post_summary = get_post_summary();
    o << "SUMMARY"
      << "\n\tPrecondition=" << pre_summary
      << "\n\tPostcondition=" << post_summary << "\n";
  }

  // for gdb debugging
  void dump() const { write(crab::outs()); }
};

/* Base class to tune the level of context-sensitivity */
template <typename CallingContext> class context_sensitivity_policy {
protected:
  unsigned int m_max_call_contexts;

public:
  using calling_context_ptr = std::unique_ptr<CallingContext>;
  using calling_context_ptr_deque = std::deque<calling_context_ptr>;

  context_sensitivity_policy(unsigned int max_call_contexts)
      : m_max_call_contexts(max_call_contexts) {}

  virtual ~context_sensitivity_policy() {}

  // Add one calling context cc in ccs
  virtual void add(calling_context_ptr_deque &ccs, calling_context_ptr cc) = 0;
};

/**
 * If limit of calling contexts is reached, then join the two oldest
 * calling contexts before adding the new calling context. Joining
 * calling contexts triggers the removal of subsumed calling contexts.
 *
 * XXX: this policy is chosen for now because it's pretty simple to
 * implement but we need to evaluate it more. We might need semantic
 * policies based on the abstract state. For instance, join two
 * calling contexts if the precision losses are small.
 **/
template <typename CallingContext>
class default_context_sensitivity_policy
    : public context_sensitivity_policy<CallingContext> {
  using context_sensitivity_policy_t =
      context_sensitivity_policy<CallingContext>;
  using abs_dom_t = typename CallingContext::abs_dom_t;

public:
  using calling_context_ptr =
      typename context_sensitivity_policy_t::calling_context_ptr;
  using calling_context_ptr_deque =
      typename context_sensitivity_policy_t::calling_context_ptr_deque;

  default_context_sensitivity_policy(unsigned int max_call_contexts)
      : context_sensitivity_policy_t(max_call_contexts) {}

  virtual void add(calling_context_ptr_deque &ccs,
                   calling_context_ptr cc) override {
    if (this->m_max_call_contexts == UINT_MAX) {
      ccs.push_back(std::move(cc));
      return;
    }

    CRAB_LOG(
        "inter3", crab::outs() << "All calling contexts before\n";
        for (unsigned i = 0, e = ccs.size(); i < e; ++i) { ccs[i]->dump(); });

    bool compress = false;
    unsigned num_of_ccs = ccs.size();
    if ((num_of_ccs >= 2) && (num_of_ccs > this->m_max_call_contexts)) {
      // -- join the two oldest contexts
      calling_context_ptr cc1 = std::move(ccs.front());
      ccs.pop_front();
      calling_context_ptr cc2 = std::move(ccs.front());
      ccs.pop_front();
      ccs.push_front(std::move(cc1->join_with(*cc2)));
      compress = true;
      CRAB_LOG("inter",
               crab::outs()
                   << "[INTER] joining two oldest calling contexts\n";);
    }
    ccs.push_back(std::move(cc));

    if (compress) {
      // -- remove redundant contexts
      assert(!ccs.empty());
      calling_context_ptr_deque new_ccs;
      const abs_dom_t &joined_pre_summary = ccs.front()->get_pre_summary();
      const abs_dom_t &joined_post_summary = ccs.front()->get_post_summary();
      auto it = ccs.begin();
      // the first one is the joined calling context so we keep it.
      new_ccs.push_back(std::move(*it));
      ++it;
      for (auto et = ccs.end(); it != et; ++it) {
        // discard any pre/post pair that is subsumed by the joined
        // calling context
        if (!((*it)->get_pre_summary() <= joined_pre_summary &&
              (*it)->get_post_summary() <= joined_post_summary)) {
          new_ccs.push_back(std::move(*it));
        }
      }
      std::swap(ccs, new_ccs);
    }

    CRAB_LOG(
        "inter3", crab::outs() << "All calling contexts after\n";
        for (unsigned i = 0, e = ccs.size(); i < e; ++i) { ccs[i]->dump(); });
  }
};

// Information needed during fixpoint of recursive calls
template <typename AbsDom> class func_fixpoint_map_entry {
public:
  using abs_dom_t = AbsDom;

private:
  abs_dom_t m_entry;
  abs_dom_t m_exit;

public:
  func_fixpoint_map_entry(abs_dom_t entry, abs_dom_t exit)
      : m_entry(entry), m_exit(exit) {}

  abs_dom_t get_entry() const { return m_entry; }

  void set_entry(abs_dom_t &&entry) { m_entry = std::move(entry); }

  abs_dom_t get_exit() const { return m_exit; }

  void set_exit(abs_dom_t &&exit) { m_exit = std::move(exit); }
};

template <typename CallGraph, typename AbsDom, typename InvariantMap>
class global_context {
  using this_type = global_context<CallGraph, AbsDom, InvariantMap>;

public:
  using callgraph_node_t = typename CallGraph::node_t;
  using cfg_t = typename callgraph_node_t::cfg_t;
  using variable_t = typename cfg_t::variable_t;
  using abs_dom_t = AbsDom;
  using invariant_map_t = InvariantMap;

  static_assert(std::is_same<typename cfg_t::basic_block_label_t,
                             typename InvariantMap::key_type>::value,
                "The key of an invariant map should be basic block label");
  static_assert(
      std::is_same<abs_dom_t, typename InvariantMap::mapped_type>::value,
      "The mapped value of an invariant map should be an abstract domain");

  using calling_context_t = calling_context<cfg_t, abs_dom_t, invariant_map_t>;
  using context_sensitivity_policy_t =
      context_sensitivity_policy<calling_context_t>;
  using default_context_sensitivity_policy_t =
      default_context_sensitivity_policy<calling_context_t>;
  using calling_context_ptr = std::unique_ptr<calling_context_t>;
  using calling_context_collection_t = std::deque<calling_context_ptr>;
  using calling_context_table_t =
      std::unordered_map<cfg_t, calling_context_collection_t>;
  using liveness_map_t =
      std::unordered_map<cfg_t, const live_and_dead_analysis<cfg_t> *>;
  using wto_cfg_t = ikos::wto<cfg_t>;
  // map cfgs to their WTO's
  using wto_cfg_map_t = std::unordered_map<cfg_t, const wto_cfg_t *>;
  using wto_cg_t = ikos::wto<crab::cg::call_graph_ref<CallGraph>>;
  using wto_cg_nesting_t = typename wto_cg_t::wto_nesting_t;
  // map callgraph entries to their WTO's
  using wto_cg_map_t =
      std::unordered_map<callgraph_node_t, std::unique_ptr<wto_cg_t>>;
  using checks_db_t = checker::checks_db;
  using global_invariant_map_t = std::unordered_map<cfg_t, invariant_map_t>;
  using func_fixpoint_map_entry_t = func_fixpoint_map_entry<abs_dom_t>;
  using func_fixpoint_map_t =
      std::unordered_map<callgraph_node_t, func_fixpoint_map_entry_t>;

  using variable_vector_t = std::vector<variable_t>;
  using region_equiv_class_t = std::vector<variable_vector_t>;
  using callsite_t = typename cfg_t::basic_block_t::callsite_t;
  using callsite_or_fdecl_t = crab::cfg::callsite_or_fdecl<cfg_t>;
  using callsite_info_t = domains::callsite_info<variable_t>;
  using callsite_info_map_t =
      std::unordered_map<const callsite_t *, callsite_info_t>;

private:
  /*
   * To deal with loops and recursive functions we compute WTOs (weak
   * topological orderings) of both the call graph and each
   * control-flow graph. In the case of the call graph, we compute a
   * separate WTO for each call graph entry.
   */

  // -- widening points of the callgraph (extracted from *all* entries' WTOs)
  std::set<callgraph_node_t> m_widening_set;
  // -- map each call site statement to static callsite info
  callsite_info_map_t m_callsite_info_map;
  // -- map each callgraph entry to its WTO
  wto_cg_map_t m_wto_cg_map;
  // -- to break cycles if no precise support for recursive functions
  // invariant: !m_call_stack.empty() &&  m_call_stack[0] is the current entry node.
  std::vector<callgraph_node_t> m_call_stack;
  // -- liveness symbols for each function
  const liveness_map_t *m_live_map;
  // -- wto for each function
  const wto_cfg_map_t *m_wto_cfg_map; // to avoid recomputing wto of cfgs
  // -- context-insensitive invariants (for external queries)
  //    populated only if keep_invariants is enabled.
  global_invariant_map_t m_pre_invariants;
  global_invariant_map_t m_post_invariants;
  // -- the policy for making tractable the number of calling contexts
  std::unique_ptr<context_sensitivity_policy_t> m_cs_policy;
  // -- for computing fixpoint of recursive functions
  func_fixpoint_map_t m_func_fixpoint_table;
  // -- enable checking interleaved with analysis
  bool m_enable_checker;
  checks_db_t m_checks_db;
  unsigned m_checker_verbosity;
  bool m_is_checking_phase;
  // -- all calling contexts
  calling_context_table_t m_cc_table;
  // -- keep context-sensitive invariants
  bool m_keep_cc_invariants;
  // -- keep context-insensitive invariants: used to populate
  //    m_pre_invariants and m_post_invariants.
  bool m_keep_invariants;
  // -- max number of calling contexts
  unsigned int m_max_call_contexts;
  // -- to analyze precisely recursive functions
  bool m_analyze_recursive_functions;
  // -- reuse summaries without losing precision
  bool m_exact_summary_reuse;
  // -- start the analysis only from main
  bool m_only_main_as_entry;
  // -- fixpoint parameters
  unsigned int m_widening_delay;
  unsigned int m_descending_iters;
  unsigned int m_thresholds_size;

  void join_with(global_invariant_map_t &global_table, cfg_t cfg,
                 invariant_map_t &other) {
    assert(m_keep_invariants);

    auto it = global_table.find(cfg);
    if (it != global_table.end()) {
      invariant_map_t &invariants = it->second;
      for (auto &kv : other) {
        auto it = invariants.find(kv.first);
        if (it != invariants.end()) {
          it->second |= kv.second;
        } else {
          invariants.insert({kv.first, kv.second});
        }
      }
    } else {
      invariant_map_t invariants;
      for (auto &kv : other) {
        invariants.insert({kv.first, kv.second});
      }
      global_table.insert(std::make_pair(cfg, std::move(invariants)));
    }
  }
  
  // Return the head of all WTO nested components of node. The result
  // is ordered. The last head belong to the innermost WTO nested
  // component and the first one is the outermost.
  boost::optional<wto_cg_nesting_t>
  get_all_wto_nested_heads(callgraph_node_t node) const {
    const callgraph_node_t &entry = get_current_entry();
    auto it = m_wto_cg_map.find(entry);
    if (it != m_wto_cg_map.end()) {
      boost::optional<wto_cg_nesting_t> res = it->second->nesting(node);
      CRAB_LOG("inter-callgraph-wto",
	       crab::outs() << "NESTING(" << node  << ")=";
	       if (res) {
		 crab::outs() << *res;
	       } else {
		 crab::outs() << "[]";
	       }
	       crab::outs() << "\n";);
      return res;
    } else {
      CRAB_ERROR("Not callgraph wto found for entry ",
		 entry.get_cfg().get_func_decl().get_func_name());
    }
  }

  // Apply fn for all nested components of root.
  class apply_fn_to_nested_wto_component_visitor
      : public wto_cg_t::wto_component_visitor_t {
  public:
    using wto_vertex_t = typename wto_cg_t::wto_vertex_t;
    using wto_cycle_t = typename wto_cg_t::wto_cycle_t;

  private:
    callgraph_node_t m_root;
    std::function<void(callgraph_node_t)> m_fn;
    bool m_active;

  public:
    apply_fn_to_nested_wto_component_visitor(
        callgraph_node_t root, std::function<void(callgraph_node_t)> fn)
        : m_root(root), m_fn(fn), m_active(false) {}

    virtual void visit(wto_cycle_t &cycle) override {
      if (cycle.head() == m_root) {
        m_active = true;
        for (auto &comp : cycle) {
          comp.accept(this);
        }
        m_active = false;
      } else {
        if (m_active) {
          m_fn(cycle.head());
        }
        for (auto &comp : cycle) {
          comp.accept(this);
        }
      }
    }
    virtual void visit(wto_vertex_t &vertex) override {
      if (m_active) {
        m_fn(vertex.node());
      }
    }
  };

public:
  global_context(const liveness_map_t *live_map,
                 const wto_cfg_map_t *wto_cfg_map, bool enable_checker,
                 unsigned checker_verbosity, bool keep_cc_invariants,
                 bool keep_invariants, unsigned int max_call_contexts,
                 bool analyze_recursive_functions, bool exact_summary_reuse,
                 bool only_main_as_entry, unsigned int widening_delay,
                 unsigned int descending_iters, unsigned int thresholds_size)
      : m_live_map(live_map), m_wto_cfg_map(wto_cfg_map),
        m_cs_policy(
            new default_context_sensitivity_policy_t(max_call_contexts)),
        m_enable_checker(enable_checker),
        m_checker_verbosity(checker_verbosity), m_is_checking_phase(false),
        m_keep_cc_invariants(keep_cc_invariants),
        m_keep_invariants(keep_invariants),
        m_max_call_contexts(max_call_contexts),
        m_analyze_recursive_functions(analyze_recursive_functions),
        m_exact_summary_reuse(exact_summary_reuse),
        m_only_main_as_entry(only_main_as_entry),
        m_widening_delay(widening_delay), m_descending_iters(descending_iters),
        m_thresholds_size(thresholds_size) {}

  global_context(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  std::set<callgraph_node_t> &get_widening_set() { return m_widening_set; }

  const callsite_info_map_t &get_callsite_info_map() const {
    return m_callsite_info_map;
  }

  callsite_info_map_t &get_callsite_info_map() { return m_callsite_info_map; }

  const std::set<callgraph_node_t> &get_widening_set() const {
    return m_widening_set;
  }

  callsite_info_t &get_call_site_info(const callsite_t &callsite) {
    auto it = m_callsite_info_map.find(&callsite);
    if (it != m_callsite_info_map.end()) {
      return it->second;
    }
    CRAB_ERROR("The callsite info for ", callsite,
               " should be computed before the anlaysis");
  }

  const callsite_info_t &get_call_site_info(const callsite_t &callsite) const {
    auto it = m_callsite_info_map.find(&callsite);
    if (it != m_callsite_info_map.end()) {
      return it->second;
    }
    CRAB_ERROR("The callsite info for ", callsite,
               " should be computed before the anlaysis");
  }

  // Return the current entry node of the analysis
  const callgraph_node_t& get_current_entry() const {
    assert(!m_call_stack.empty());
    return m_call_stack[0];
  }
  
  wto_cg_map_t &get_wto_cg_map() { return m_wto_cg_map; }

  const wto_cg_map_t &get_wto_cg_map() const { return m_wto_cg_map; }

  // Return true if the node is part of a non-empty wto nesting
  bool included_nested_wto_component(callgraph_node_t node) const {
    if (auto nesting_opt = get_all_wto_nested_heads(node)) {
      auto nesting = *nesting_opt;
      return nesting.begin() != nesting.end();
    }
    CRAB_ERROR(node.get_cfg().get_func_decl().get_func_name(),
               " should be part of some WTO component");
  }

  void
  apply_fn_to_nested_wto_component(callgraph_node_t root,
                                   std::function<void(callgraph_node_t)> fn) {
    apply_fn_to_nested_wto_component_visitor vis(root, fn);
    m_wto_cg_map[get_current_entry()]->accept(&vis);
  }

  bool find_call_stack(const callgraph_node_t &node) const {
    auto it = std::find(m_call_stack.begin(), m_call_stack.end(),
			node);
    return it!= m_call_stack.end();
  }

  void print_call_stack() const {
    crab::outs() << "callstack depth=" << m_call_stack.size()
		 << " content=[";
    for(unsigned i=0, sz=m_call_stack.size(); i<sz;) {
      crab::outs() << m_call_stack[i];
      ++i;
      if (i < sz) {
	crab::outs() << ",";
      }
    }
    crab::outs() << "]\n";
  }
  
  std::vector<callgraph_node_t> &get_call_stack() { return m_call_stack; }

  const liveness_map_t *get_live_map() const { return m_live_map; }

  const wto_cfg_map_t *get_wto_cfg_map() const { return m_wto_cfg_map; }

  context_sensitivity_policy_t &get_context_sensitivity_policy() {
    return *m_cs_policy;
  }

  bool run_checker() const { return m_enable_checker; }

  checks_db_t &get_checks_db() { return m_checks_db; }

  const checks_db_t &get_checks_db() const { return m_checks_db; }

  unsigned get_checker_verbosity() const { return m_checker_verbosity; }

  bool &get_is_checking_phase() { return m_is_checking_phase; }

  bool get_is_checking_phase() const { return m_is_checking_phase; }

  calling_context_table_t &get_calling_context_table() { return m_cc_table; }

  const calling_context_table_t &get_calling_context_table() const {
    return m_cc_table;
  }

  func_fixpoint_map_t &get_func_fixpoint_table() {
    return m_func_fixpoint_table;
  }

  bool keep_cc_invariants() const { return m_keep_cc_invariants; }

  bool keep_invariants() const { return m_keep_invariants; }

  unsigned int get_max_call_contexts() const { return m_max_call_contexts; }

  bool analyze_recursive_functions() const {
    return m_analyze_recursive_functions;
  }

  bool exact_summary_reuse() const { return m_exact_summary_reuse; }

  bool only_main_as_entry() const { return m_only_main_as_entry; }

  unsigned int get_widening_delay() const { return m_widening_delay; }

  unsigned int get_descending_iters() const { return m_descending_iters; }

  unsigned int get_thresholds_size() const { return m_thresholds_size; }

  // context-insensitive invariants for each function (if
  // m_keep_invariants enabled)

  global_invariant_map_t &get_global_pre_invariants() {
    return m_pre_invariants;
  }

  const global_invariant_map_t &get_global_pre_invariants() const {
    return m_pre_invariants;
  }

  global_invariant_map_t &get_global_post_invariants() {
    return m_post_invariants;
  }

  const global_invariant_map_t &get_global_post_invariants() const {
    return m_post_invariants;
  }

  void join_invariants_with(callgraph_node_t cg_node,
                            invariant_map_t &pre_invariants,
                            invariant_map_t &post_invariants) {
    if (m_keep_invariants) {
      join_with(get_global_pre_invariants(), cg_node.get_cfg(), pre_invariants);
      join_with(get_global_post_invariants(), cg_node.get_cfg(),
                post_invariants);
    }
  }
};

template <typename CallGraphNode, typename IntraCallSemAnalyzer,
          typename GlobalCtx>
void check_function(CallGraphNode cg_node, IntraCallSemAnalyzer &analyzer,
                    GlobalCtx &ctx) {
  if (ctx.run_checker()) {
    using intra_checker_t = checker::intra_checker<IntraCallSemAnalyzer>;
    using assertion_checker_t =
        checker::assert_property_checker<IntraCallSemAnalyzer>;
    /**
     * Make sure that the checker is run **after** the whole WTO component
     * has been stabilized.
     **/
    bool &is_checking_phase = ctx.get_is_checking_phase();
    assert(!is_checking_phase);
    is_checking_phase = true;
    typename intra_checker_t::prop_checker_ptr p(
        new assertion_checker_t(ctx.get_checker_verbosity()));
    intra_checker_t checker(analyzer, {p});
    checker.run();
    ctx.get_checks_db() += checker.get_all_checks();
    assert(is_checking_phase);
    is_checking_phase = false;
  }
}

// Wrapper to call the intra-procedural analysis with inter-procedural
// semantics for call and return statements.
template <typename CallGraphNode, typename IntraCallSemAnalyzer>
IntraCallSemAnalyzer *
analyze_function(CallGraphNode cg_node,
                 typename IntraCallSemAnalyzer::abs_tr_t &abs_tr,
                 unsigned iteration) {
  using cfg_t = typename CallGraphNode::cfg_t;
  using abs_dom_t = typename IntraCallSemAnalyzer::abs_dom_t;

  TD_INTER_TIMER_STATS(TimerInterAnalyzeFunc);
  
  auto &ctx = abs_tr.get_context();
  cfg_t cfg = cg_node.get_cfg();
  assert(cfg.has_func_decl());
  auto &func_fixpoint_table = ctx.get_func_fixpoint_table();
  auto &entry = abs_tr.get_abs_value();

  CRAB_VERBOSE_IF(1, get_msg_stream()
                         << "++ Analyzing function  "
                         << cfg.get_func_decl().get_func_name() << "\n";);

  CRAB_LOG("inter-callstack", ctx.print_call_stack());
	   
  if (ctx.analyze_recursive_functions() &&
      ctx.get_widening_set().count(cg_node) > 0) {
    // ### Recursive function ###
    //
    // Create the initial value for the function for each fixpoint
    // iteration.
    auto it = func_fixpoint_table.find(cg_node);
    if (it == func_fixpoint_table.end()) {
      auto exit = abs_tr.get_abs_value().make_bottom();
      func_fixpoint_map_entry<typename IntraCallSemAnalyzer::abs_dom_t>
          fixpoint_start(entry, exit);
      func_fixpoint_table.insert({cg_node, fixpoint_start});
      // increment the counter only the first time
      TD_INTER_COUNT_STATS(CounterRecursiveCalls);	  
    } else {
      entry = it->second.get_entry();
    }
  }
  CRAB_LOG("inter", crab::outs() << "[INTER] Initial abstract state for "
                                 << cfg.get_func_decl().get_func_name() << ": "
                                 << entry << "\n";);

  IntraCallSemAnalyzer *analyzer = nullptr;
  if (abs_tr.has_analyzer(cg_node)) {
    /// --- 1. reuse the intra-analyzer if already exists
    // don't create another analyzer if we already created one
    analyzer = &(abs_tr.get_analyzer(cg_node));
    // make sure no results from previous run
    /// analyzer->clear();
  } else {
    // get liveness symbols for cfg if available
    const live_and_dead_analysis<cfg_t> *live = nullptr;
    if (ctx.get_live_map()) {
      auto it = ctx.get_live_map()->find(cfg);
      if (it != ctx.get_live_map()->end()) {
        live = it->second;
      }
    }
    // get wto for cfg if available
    const ikos::wto<cfg_t> *wto = nullptr;
    if (ctx.get_wto_cfg_map()) {
      auto it = ctx.get_wto_cfg_map()->find(cfg);
      if (it != ctx.get_wto_cfg_map()->end()) {
        wto = it->second;
      }
    }
    /// -- 2. Create intra analyzer (with inter-procedural semantics for
    /// call/return)
    std::unique_ptr<IntraCallSemAnalyzer> new_analyzer(new IntraCallSemAnalyzer(
        cfg, wto, &abs_tr, live, ctx.get_widening_delay(),
        ctx.get_descending_iters(), ctx.get_thresholds_size()));
    analyzer = &(abs_tr.add_analyzer(cg_node, std::move(new_analyzer)));
  }

  // entry is a reference that will be modified by the analysis
  // unused if function not recursive
  abs_dom_t old_entry(entry);

  /// -- 3. Run the analyzer
  analyzer->run_forward();

  // ### Recursive function ###
  // Re-analyze the function until fixpoint
  auto it = func_fixpoint_table.find(cg_node);
  if (it != func_fixpoint_table.end()) {
    auto &fixpo_val = it->second;
    abs_dom_t new_entry = fixpo_val.get_entry();
    abs_dom_t old_exit = fixpo_val.get_exit();
    abs_dom_t new_exit = old_exit.make_bottom();
    if (cfg.has_exit()) {
      new_exit = analyzer->get_post(cfg.exit());
    }
    bool fixpoint_reached = (new_entry <= old_entry && new_exit <= old_exit);
    CRAB_LOG(
        "inter", crab::outs()
                     << "[INTER] Checking fixpoint termination\n\tOLD entry="
                     << old_entry << "\n\tNEW entry=" << new_entry
                     << "\n\tOLD exit=" << old_exit << "\n\tNEW exit="
                     << new_exit << "\nRES=" << fixpoint_reached << "\n";);

    if (fixpoint_reached) {
      func_fixpoint_table.erase(cg_node);
      CRAB_VERBOSE_IF(1, get_msg_stream()
                             << "++ Fixpoint reached for recursive function "
                             << cfg.get_func_decl().get_func_name() << "!\n";);
      // Don't check invariants with the last iteration
      return nullptr;
    } else {
      CRAB_VERBOSE_IF(1, get_msg_stream()
                             << "++ Widening " << iteration
                             << " for recursive function "
                             << cfg.get_func_decl().get_func_name() << "\n";);

      ///// widen entry and exit
      if (iteration >= ctx.get_widening_delay()) {
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Initial abstract state for "
                              << cfg.get_func_decl().get_func_name() << ":\n"
                              << "\tWIDEN\n\t  OLD=" << old_entry
                              << "\n\t  NEW=" << new_entry << "\n";);
        new_entry = old_entry || new_entry;
        CRAB_LOG("inter", crab::outs() << "\t  RES=" << new_entry << "\n";);
        new_exit = old_exit || new_exit;
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Exit abstract state for "
                              << cfg.get_func_decl().get_func_name() << ":\n"
                              << "\tWIDEN\n\t  OLD=" << old_exit
                              << "\n\t  NEW=" << new_exit << "\n";);
        CRAB_LOG("inter", crab::outs() << "\t  RES=" << new_exit << "\n";);
      } else {
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Initial abstract state for "
                              << cfg.get_func_decl().get_func_name() << ":\n"
                              << "\tJOIN\n\t  OLD=" << old_entry
                              << "\n\t  NEW=" << new_entry << "\n";);
        new_entry |= old_entry;
        CRAB_LOG("inter", crab::outs() << "\t  RES=" << new_entry << "\n";);
        new_exit |= old_exit;
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Exit abstract state for "
                              << cfg.get_func_decl().get_func_name() << ":\n"
                              << "\tJOIN\n\t  OLD=" << old_exit
                              << "\n\t  NEW=" << new_exit << "\n";);
        CRAB_LOG("inter", crab::outs() << "\t  RES=" << new_exit << "\n";);
      }

      fixpo_val.set_entry(std::move(new_entry));
      fixpo_val.set_exit(std::move(new_exit));
      // continue next fixpoint iteration
      if (auto fixpo_analyzer =
              analyze_function<CallGraphNode, IntraCallSemAnalyzer>(
                  cg_node, abs_tr, ++iteration)) {
        return fixpo_analyzer;
      }
    }
  }

  /// If we are here is because either the function was not recursive
  /// or we are in the penultimate fixpoint iteration (i.e,
  /// post-fixpoint).

  /// -- 4. Store the invariants (merging with other contexts)
  ctx.join_invariants_with(cg_node, analyzer->get_pre_invariants(),
                           analyzer->get_post_invariants());

  CRAB_VERBOSE_IF(1, get_msg_stream()
		  << "++ Finished analysis of function  "
		  << cfg.get_func_decl().get_func_name() << "\n";);
  CRAB_LOG("inter", 
	   crab::outs() << "Entry=" << analyzer->get_pre(cfg.entry()) << "\n"
	                << "Exit=" << analyzer->get_post(cfg.entry())
	                << "\n";);

  CRAB_LOG("fixpo-trace",
           /// useful for debugging: reset the counters so that they are not
           /// accumulated from run to run.
           reset_wto_cycle_counter_visitor<cfg_t> vis;
           auto &wto = analyzer->get_wto(); wto.accept(&vis););

  // Note that we cannot free analyzer's state here because it is used
  // after.
  return analyzer;
}

/// @brief wrapper to group regions variables in parameter lists based on
/// target's dsa intrinsics
/// @param tgt_cfg the cfg for target function
/// @param src_params the paramters appeared in the source function
/// @param tgt_params the paramters appeared in the target function
/// @param src_cls the equivalence classes for grouping source regions
/// @param tgt_cls the equivalence classes for grouping target regions
template <typename CFG, typename Variable>
void group_regions_by_dsa_intrinsics(
    const CFG &tgt_cfg, const std::vector<Variable> &src_params,
    const std::vector<Variable> &tgt_params,
    std::vector<std::vector<Variable>> &src_cls,
    std::vector<std::vector<Variable>> &tgt_cls) {

  using basic_block_t = typename CFG::basic_block_t;
  using intrinsic_statement_t = typename basic_block_t::intrinsic_t;
  using variable_or_constant_t =
      typename intrinsic_statement_t::variable_or_constant_t;
  // Helper lambda functions for logging
  auto print_var_vector = [](crab::crab_os &o,
                             const std::vector<Variable> &vec) {
    typename std::vector<Variable>::const_iterator it;
    o << "[";
    for (it = vec.begin(); it != vec.end(); it++) {
      if (it != vec.begin())
        o << ",";
      o << (*it);
    }
    o << "]";
  };
  auto print_var_or_const_vector =
      [](crab::crab_os &o, const std::vector<variable_or_constant_t> &vec) {
        typename std::vector<variable_or_constant_t>::const_iterator it;
        o << "[";
        for (it = vec.begin(); it != vec.end(); it++) {
          if (it != vec.begin())
            o << ",";
          o << (*it);
        }
        o << "]";
      };

  CRAB_LOG("inter-group", errs() << "Ins: ";
           print_var_vector(errs(), src_params); errs() << "\n";);
  CRAB_LOG("inter-group", errs() << "Outs: ";
           print_var_vector(errs(), tgt_params); errs() << "\n";);

  // Step 1. get all dsa intrinsic statements
  const basic_block_t &entry_bb = tgt_cfg.get_node(tgt_cfg.entry());
  std::vector<intrinsic_statement_t> intrinsic_stmts;
  intrinsic_stmts.reserve(entry_bb.size());
  for (auto const &stmt : entry_bb) {
    if (stmt.is_intrinsic()) {
      auto intrinsic_stmt = dynamic_cast<const intrinsic_statement_t &>(stmt);
      if (intrinsic_stmt.get_intrinsic_name() == "regions_from_memory_object") {
        intrinsic_stmts.push_back(intrinsic_stmt);
      }
    }
  }

  // Step 2. construct the equivalence classes for grouping input/output
  // regions Note that, the order of the output regions in the parameter list
  // is the same as the output list in the dsa intrinsic
  // The worst case of the following loop is O(m * n) where m is the number of
  // parameters and n is the number of dsa intrinsic statements
  unsigned int index = 0, param_lst_sz = src_params.size();
  while (index < param_lst_sz) {
    if (!tgt_params[index].get_type().is_region()) {
      // skip for reference and scalar variables
      index++;
      continue;
    }
    for (auto const &intrinsic_stmt : intrinsic_stmts) {
      // each intrinsic contains only region variables
      std::vector<Variable> srcs, tgts;
      const unsigned intrinsic_arg_sz = intrinsic_stmt.get_num_args();
      auto const &args = intrinsic_stmt.get_args();
      srcs.reserve(intrinsic_arg_sz);
      tgts.reserve(intrinsic_arg_sz);
      CRAB_LOG("inter-group", errs() << "Intrinsic: ";
                 print_var_or_const_vector(errs(), args); errs() << "\n";);
      for (unsigned j = 0; j < intrinsic_arg_sz; j++) {
        // If the region whose current index refers belongs to some dsa
        // intrinsic, add this group into the equivalence sets
        if (index < param_lst_sz && args[j].get_variable() == tgt_params[index]) {
          srcs.push_back(src_params[index]);
          tgts.push_back(tgt_params[index]);
          index++;
        }
      }
      if (srcs.size() > 0) {
        src_cls.push_back(srcs);
        tgt_cls.push_back(tgts);
      }
    }
    index++;
  }
}

// Wrapper to retrieve all callsite info from intrisinc statements.
// The callsite info will be stored inside an object of class
// domains::callsite_info for each callsite statement. This analysis is
// performed before the inter-procedural analysis. input: global_context a.
// iterate all cfg inside wto of cfgs b. for each cfg, visit callsite statement
// c. for each callsite statement, get callee'cfg
// d. store callsite info
template <typename CallGraph, typename GlobalCtx>
void extract_callsite_info(CallGraph &cg, GlobalCtx &ctx) {
  using cfg_t = typename CallGraph::cfg_t;
  using callsite_t = typename cfg_t::basic_block_t::callsite_t;
  using fdecl_t = typename cfg_t::fdecl_t;
  using callsite_info_t = typename GlobalCtx::callsite_info_t;
  using region_equiv_class_t = typename GlobalCtx::region_equiv_class_t;

  // Iterate all cfgs through the call graph
  for (auto const &node : boost::make_iterator_range(cg.nodes())) {
    const cfg_t &caller_cfg = node.get_cfg();
    for (auto const &bb : caller_cfg) {
      for (auto const &stmt : bb) {
        if (stmt.is_callsite()) {
          const callsite_t &cs = static_cast<const callsite_t &>(stmt);
          // get callee's cfg
          if (cg.has_callee(cs)) {
            const cfg_t &callee_cfg = cg.get_callee(cs).get_cfg();
            const fdecl_t &fdecl = callee_cfg.get_func_decl();
            region_equiv_class_t formal_cls, actual_cls, ret_cls, lhs_cls;
            group_regions_by_dsa_intrinsics(callee_cfg, cs.get_args(),
                                            fdecl.get_inputs(), actual_cls,
                                            formal_cls);
            group_regions_by_dsa_intrinsics(caller_cfg, fdecl.get_outputs(),
                                            cs.get_lhs(), ret_cls, lhs_cls);
            callsite_info_t callsite_info{cs.get_func_name(),
                                          cs.get_args(),
                                          cs.get_lhs(),
                                          fdecl.get_inputs(),
                                          fdecl.get_outputs(),
                                          std::move(actual_cls),
                                          std::move(lhs_cls),
                                          std::move(formal_cls),
                                          std::move(ret_cls)};
            auto &callsite_info_map = ctx.get_callsite_info_map();
            callsite_info_map.insert({&cs, callsite_info});
          }
        }
      }
    }
  }
}

/*
   All the state of the inter-procedural analysis is kept in this
   transformer. This allows us to use the intra-procedural analysis
   fully as a black box.
*/
template <typename CallGraph, typename AbsDom>
class top_down_inter_transformer final
    : public intra_abs_transformer<
          typename CallGraph::node_t::cfg_t::basic_block_t, AbsDom> {
  
  using cg_node_t = typename CallGraph::node_t;
  using cfg_t = typename cg_node_t::cfg_t;
  using wto_cg_t = ikos::wto<CallGraph>;
  using basic_block_t = typename cfg_t::basic_block_t;
  using fdecl_t = typename cfg_t::fdecl_t;
  using variable_t = typename cfg_t::variable_t;
  using intra_abs_transformer_t = intra_abs_transformer<basic_block_t, AbsDom>;
  using abs_transform_api_t =
      typename intra_abs_transformer_t::abs_transform_api_t;

  using this_type = top_down_inter_transformer<CallGraph, AbsDom>;
  // -- intra-procedural analysis with inter-procedural semantics for
  //    call/return statements
  using intra_analyzer_with_call_semantics_t =
      analyzer_internal_impl::fwd_analyzer<cfg_t, this_type>;

public:
  using abs_dom_t = AbsDom;
  using callsite_t = typename abs_transform_api_t::callsite_t;
  // -- global context
  // XXX: we cannot use
  // intra_analyzer_with_call_semantics_t::invariant_map_t because cyclic
  // dependencies
  using invariant_map_t =
      std::unordered_map<typename cfg_t::basic_block_label_t, abs_dom_t>;
  using global_context_t =
      top_down_inter_impl::global_context<CallGraph, abs_dom_t,
                                          invariant_map_t>;
  // -- calling context stuff
  using calling_context_t = typename global_context_t::calling_context_t;
  using calling_context_ptr = typename global_context_t::calling_context_ptr;
  using calling_context_collection_t =
      typename global_context_t::calling_context_collection_t;
  using calling_context_table_t =
      typename global_context_t::calling_context_table_t;
  using func_fixpoint_map_entry_t =
      typename global_context_t::func_fixpoint_map_entry_t;

private:
  // -- the callgraph
  CallGraph &m_cg;
  // -- global parameters of the analysis
  //    the caller owns the pointer
  global_context_t &m_ctx;
  // -- to avoid starting from scratch an analyzer:
  //    saving wto computation, type checking, etc.
  // XXX: it cannot be in m_ctx because of cyclic dependencies.
  std::unordered_map<cg_node_t,
                     std::unique_ptr<intra_analyzer_with_call_semantics_t>>
      m_intra_analyzer;
  // recursion depth (for debugging)
  unsigned m_depth;
  // Produce warning if m_depth >= m_max_depth.
  // If m_max_depth is large enough then it's a strong indication that
  // we failed in detecting a cycle in the call graph and the analysis
  // will not terminate.
  const unsigned m_max_depth = 500;

  inline abs_dom_t make_top() {
    auto const &dom = this->get_abs_value();
    return dom.make_top();
  }

  inline abs_dom_t make_bottom() {
    auto const &dom = this->get_abs_value();
    return dom.make_bottom();
  }

  calling_context_collection_t &get_calling_contexts(cfg_t fun) {
    return m_ctx.get_calling_context_table()[fun];
  }

  void add_calling_context(cfg_t fun, calling_context_ptr cc) {
    auto it = m_ctx.get_calling_context_table().find(fun);
    if (it == m_ctx.get_calling_context_table().end()) {
      calling_context_collection_t ccs;
      ccs.push_back(std::move(cc));
      // XXX: don't use braced-init-list. For some compiler versions,
      // initializer_list only allows const access to its elements so
      // it will try to copy the unique ptr.
      m_ctx.get_calling_context_table().insert(
          std::make_pair(fun, std::move(ccs)));
    } else {
      // the policy decides how to add the new context.
      auto &cs_policy = m_ctx.get_context_sensitivity_policy();
      calling_context_collection_t &ccs = it->second;
      cs_policy.add(ccs, std::move(cc));
    }
  }

  static void get_fdecl_parameters(const fdecl_t &fdecl,
                                   std::vector<variable_t> &out) {
    out.reserve(fdecl.get_num_inputs() + fdecl.get_num_outputs());
    out.insert(out.end(), fdecl.get_inputs().begin(), fdecl.get_inputs().end());
    out.insert(out.end(), fdecl.get_outputs().begin(),
               fdecl.get_outputs().end());
  }

  // Return true if the same function has been analyzed with with an
  // abstract state more general than callee_at_exit.
  bool check_cache(const callsite_t &cs,
		   cfg_t callee_cfg,
		   const bool inside_recursive_call,
		   const AbsDom &callee_at_entry,
		   AbsDom &callee_at_exit) const {
    TD_INTER_TIMER_STATS(TimerInterCheckCache);    
    
    
    auto it = m_ctx.get_calling_context_table().find(callee_cfg);
    if (it != m_ctx.get_calling_context_table().end()) {
      auto &call_contexts = it->second;
      CRAB_LOG("inter-subsume",
	  crab::outs() << "[INTER] Subsumption check at " << cs << "\n";);      
      CRAB_LOG("inter-subsume", if (call_contexts.empty()) {
	  crab::outs() << "[INTER] There is no call contexts stored.\n";
	});
      for (unsigned i = 0, e = call_contexts.size(); i < e; ++i) {
        // If the call is recursive then we cannot use exact
        // subsumption. Otherwise, it's very likely that subsumption
        // never succeeds. Apart from not having reusing, it will
        // create problems during the checking phase which assumes
        // that all function calls are always cached.
        const bool use_exact_subsumption =
            (!inside_recursive_call && m_ctx.exact_summary_reuse());

        CRAB_LOG("inter-subsume",
		 if (use_exact_subsumption) {
		   crab::outs() << "Exact ";
		 } else {
		   crab::outs() << "Approximate ";		   
		 }
                 crab::outs() << "checking if\n"
		              << callee_at_entry << "\nis subsumed by summary "
		              << i << "\n";
                 call_contexts[i]->write(crab::outs()); crab::outs() << "\n";);
	
        if (call_contexts[i]->is_subsumed(callee_at_entry,
                                          use_exact_subsumption)) {
          CRAB_LOG("inter-subsume", crab::outs() << "succeed!\n";);
          CRAB_VERBOSE_IF(1, get_msg_stream()
                                 << "++ Skip redundant analysis of function  "
                                 << cs.get_func_name() << "\n";);

          callee_at_exit = call_contexts[i]->get_post_summary();
	  return true;
          break;
        } else {
          CRAB_LOG("inter-subsume", crab::outs() << "failed!\n";);
        }
      }
    } else {
      CRAB_LOG("inter-subsume", 
	       crab::outs() << "[INTER] There is no call contexts stored.\n";);
    }
    return false;
  }

  
  /* Analysis of a callsite */
  void analyze_callee(const callsite_t &cs, cg_node_t callee_cg_node) {
    // 1. Get CFG from the callee
    cfg_t callee_cfg = callee_cg_node.get_cfg();

    assert(callee_cfg.has_func_decl());
    const fdecl_t &fdecl = callee_cfg.get_func_decl();

    // 2. Generate initial abstract state for the callee
    AbsDom caller(this->get_abs_value());

    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = caller.is_bottom();
    }

    AbsDom callee_at_entry = make_top();
    if (m_ctx.analyze_recursive_functions() ||
	m_ctx.get_widening_set().count(callee_cg_node) <= 0) {
      // If we do not analyze precisely recursive functions then we
      // must start the analysis of a recursive procedure without
      // propagating from caller to callee (i.e., top).
      const domains::callsite_info<variable_t> &callsite =
          m_ctx.get_call_site_info(cs);
      callee_at_entry.callee_entry(callsite, caller);
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = callee_at_entry.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR(
            "Invariant became bottom after propagating from actual to formals ",
            cs);
      }
    }

    AbsDom callee_at_exit = make_top();
    std::vector<variable_t> callee_at_exit_vars;
    get_fdecl_parameters(fdecl, callee_at_exit_vars);

    const bool inside_recursive_call =
        (m_ctx.analyze_recursive_functions() &&
         m_ctx.get_widening_set().count(callee_cg_node) > 0);
    
    // 3. Check if the same call context has been seen already
    bool cache_hit = check_cache(cs, callee_cfg,
				 inside_recursive_call,
				 callee_at_entry, callee_at_exit);    
    if (cache_hit) {
      if (!m_ctx.get_is_checking_phase()) {
        TD_INTER_COUNT_STATS(CounterReusedCalls);
      }
    } else {
      // if (m_ctx.get_is_checking_phase()) {
      //   CRAB_ERROR("in checking phase we should not analyze the callsite ", cs);
      // }

      intra_analyzer_with_call_semantics_t *callee_analysis = nullptr;
      auto &func_fixpoint_table = m_ctx.get_func_fixpoint_table();
      auto it = func_fixpoint_table.find(callee_cg_node);
      if (it != func_fixpoint_table.end()) {
        // ### Precise analysis of recursive call ###
        // 4.a Replace recursive call with a pre-fixpoint.
        callee_at_exit = it->second.get_exit();
        abs_dom_t callee_at_entry_copy(callee_at_entry);
	callee_at_entry_copy |= it->second.get_entry();
	// Super important for termination
        it->second.set_entry(std::move(callee_at_entry_copy)); 
        CRAB_LOG("inter", callee_at_entry = it->second.get_entry();
                 crab::outs()
                 << "[INTER] Replaced recursive call  \"" << cs
                 << "\" with pre-fixpoint=" << callee_at_exit << "\n"
                 << "Next analysis with entry=" << callee_at_entry << "\n");
	assert(!callee_analysis);
      } else {
	if (m_ctx.find_call_stack(callee_cg_node)) {
	  // ### Imprecise analysis of recursive call ###
	  // 4.b Replace recursive call with top.
	  //
	  // This is reachable also during checking phase. E.g. with
	  // this code:
	  //   foo { ... bar(); ... }
	  //   bar { ... foo(); ... }
	  // 
	  // When the checker runs on bar we won't have a summary for
	  // foo since its analysis is not completed yet.
	  callee_at_exit.set_to_top();
          TD_INTER_COUNT_STATS(CounterRecursiveCalls);
	  CRAB_VERBOSE_IF(1, get_msg_stream()
			  << "++ Skipped analysis of recursive callee \""
			  << cs.get_func_name() << "\"\n";);
	  assert(!callee_analysis);	  
	} else {
	  // ### Non-recursive call ###
	  
	  if (m_ctx.get_is_checking_phase()) {
	    m_ctx.print_call_stack();
	    CRAB_ERROR("in checking phase we should not analyze the callsite ", cs);
	  }
	  TD_INTER_COUNT_STATS(CounterAnalyzedCalls);
	  
	  // 4.c Run intra analyzer on the callee
	  // 
	  // The callee can be a recursive function but this call does
	  // not produce a cycle yet so we can analyze it.
	  
	  abs_dom_t callee_at_entry_copy(callee_at_entry);
	  this->set_abs_value(std::move(callee_at_entry_copy));

	  CRAB_LOG("inter", crab::outs() << "[INTER] Started \"";
		   crab::outs()
		   << cs << "\" with entry=" << callee_at_entry << "\n";);

	  m_ctx.get_call_stack().push_back(callee_cg_node);
	  callee_analysis = top_down_inter_impl::analyze_function<
	    typename CallGraph::node_t, intra_analyzer_with_call_semantics_t>(
		   callee_cg_node, *this, 0);
	  m_ctx.get_call_stack().pop_back();

	  // callee_analysis should be only null if the function is
	  // recursive and its fixpoint converges in one iteration which
	  // should only happen if the invariants are bottom.
	  if (callee_analysis && callee_cfg.has_exit()) {
	    callee_at_exit = callee_analysis->get_post(callee_cfg.exit());
	    CRAB_LOG("inter", crab::outs()
		     << "[INTER] Finished \"" << cs
		     << "\" with exit=" << callee_at_exit << "\n";);
	    
	  } else {
	    // if the callee has not exit is because it's a noreturn function.
	    callee_at_exit.set_to_bottom();
	    CRAB_LOG("inter", crab::outs()
		     << "[INTER] Finished \"" << cs
		     << "\" with exit=" << callee_at_exit
		     << ": the callee has no exit block.\n";);
	  }
	  
	  if (callee_analysis && inside_recursive_call) {
	    // After the analysis of the recursive function converges,
	    // We cannot store the summary with the initial abstract
	    // state before the fixpoint started. Instead, we need to
	    // update "callee_at_entry" by taking the invariant at the
	    // entry of the function after the fixpoint converged.
	    callee_at_entry = callee_analysis->get_pre(callee_cfg.entry());
	  }
	  
	  CRAB_LOG("inter2", if (callee_analysis) {
	      callee_analysis->write(crab::outs());
	      crab::outs() << "\n";
	    });
	}
      }  // end analysis of callee
      
      callee_at_exit.project(callee_at_exit_vars);


      // Check if the fixpoint of node has been stabilized.
      auto has_been_stabilized = [this](const typename CallGraph::node_t &node) {
	  if (!m_ctx.analyze_recursive_functions()) {
	    return true;
	  }
	  // node is being currently processed.
	  auto &func_fixpoint_table = m_ctx.get_func_fixpoint_table();
	  if (func_fixpoint_table.find(node) != func_fixpoint_table.end()) {
	    return false;
	  }

	  // If node is part of a WTO nesting then we need to check
	  // all the nesting's elements have been stabilized.
	  if (boost::optional<typename global_context_t::wto_cg_nesting_t> nesting_opt =
	      m_ctx.get_wto_cg_map()[m_ctx.get_current_entry()]->nesting(node)) {
	    for (auto it=(*nesting_opt).begin(), et=(*nesting_opt).end() ; it!=et; ++it) {
	      if (func_fixpoint_table.find(*it) != func_fixpoint_table.end()) {
		return false;
	      }
	    }
	  }
	  return true;
      };

      if (callee_analysis && has_been_stabilized(callee_cg_node)) {
	// 5. Add the new calling context.
	TD_INTER_TIMER_STATS(TimerInterStoreSum);	
	invariant_map_t pre_invariants, post_invariants;
	calling_context_ptr cc(new calling_context_t(
            fdecl, callee_at_entry, callee_at_exit,
            false /*ignore pre_invariants, post_invariants*/,
            std::move(pre_invariants), std::move(post_invariants)));

	add_calling_context(callee_cfg, std::move(cc));
	CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Stored summary for "
			<< fdecl.get_func_name() << "\n";);
	CRAB_LOG("inter",
		 crab::outs() << fdecl << "\n"
		 << "=== Preconditions  ===\n" << callee_at_entry << "\n"
		 << "=== Postconditions ===\n" << callee_at_exit << "\n"
		 << "=====================\n";);
      }

      if (callee_analysis && !m_ctx.included_nested_wto_component(callee_cg_node)) {
	/***
	 *** We delay running the checker until the node does not
	 *** belong to any nested WTO component.
	 ***/
      
	// 6. Check assertions within the whole WTO component.
	auto check = [this](cg_node_t node) {		       
		       if (has_analyzer(node)) {
			 // make sure that the analysis of the
			 // function has been completed
			 if (m_ctx.find_call_stack(node)){
			   return;
			 }
			 CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Running checker on "
					 << node.name() << "\n";);
			 auto &analysis = get_analyzer(node);
			 check_function(node, analysis, m_ctx);
			 CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Finished checker on "
					 << node.name() << "\n";);
		       }
		     };
	
	/// Run the checker on the callee
	check(callee_cg_node);
	
	if (m_ctx.get_widening_set().count(callee_cg_node) > 0) {
	  /// Run the checker recursively on all the callee's nested
	  /// components.
	  CRAB_VERBOSE_IF(1, get_msg_stream()
			  << "++ Running RECURSIVELY the checker starting from "
			  << callee_cg_node << " on all WTO nested components.\n";);
	  m_ctx.apply_fn_to_nested_wto_component(callee_cg_node, check);
	}

	/// Free the callee
	
	// We should be able to free the analysis of the callee if it
	// is not part of any wto nested component.
	
	// CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Free analysis of "
	//		  << callee_cg_node.name() << "\n";);

	// FIXME(07/23/21): free the analysis causes problems with curl program.
	// callee_analysis->clear();
      }
    } // end cache_hit

    pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = callee_at_exit.is_bottom();
    }


    CRAB_LOG("inter-extend",
             crab::outs() << "[INTER] Computing continuation for " << cs << "\n"
                          << "\tCaller before the call=" << caller << "\n"
                          << "\tCallee exit=" << callee_at_exit << "\n";);
    
    // 6. Generate abstract state for the continuation at the caller
    const domains::callsite_info<variable_t> &callsite =
        m_ctx.get_call_site_info(cs);
    caller.caller_continuation(callsite, callee_at_exit);
						     
    if (::crab::CrabSanityCheckFlag && !inside_recursive_call) {
      // If the function is recursive, then caller_cont_dom can be
      // bottom, even after the first fixpoint iteration.
      bool post_bot = caller.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after propagating from formals to "
                   "actuals ",
                   cs);
      }
    }

    CRAB_LOG("inter", crab::outs() << "[INTER] Continuation for \"" << cs
                                   << "\"=" << caller << "\n";);
    this->set_abs_value(std::move(caller));
  }

public:
  top_down_inter_transformer(CallGraph &cg, global_context_t &ctx, AbsDom init)
    : intra_abs_transformer_t(init), m_cg(cg), m_ctx(ctx), m_depth(0) {}

  top_down_inter_transformer(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  const global_context_t &get_context() const { return m_ctx; }

  global_context_t &get_context() { return m_ctx; }

  virtual void exec(callsite_t &cs) override {
    if (!m_cg.has_callee(cs)) {
      CRAB_ERROR("Cannot find callee CFG for ", cs);
    }

    if (this->get_abs_value().is_bottom()) {
      CRAB_LOG("inter", crab::outs() << "[INTER] skip callsite " << cs
                                     << " because abstract state is bottom\n";);
      return;
    }

    if (!m_ctx.get_is_checking_phase()) {
      TD_INTER_COUNT_STATS(CounterCalls);
    }

    cg_node_t callee_cg_node = m_cg.get_callee(cs);
    auto analyze_callee_with_depth_check = [this, &cs](cg_node_t callee) {
      ++m_depth;
      if (m_depth >= m_max_depth) {
        CRAB_ERROR("The recursion depth exceeded ", m_max_depth, ". ",
                   "This might be legitimate if the program is too deep, ",
                   "but it can also indicate a non-termination problem so we "
                   "prefer to abort the analysis.");
      }
      CRAB_LOG("inter-depth", for (unsigned i = 0; i < m_depth;
                                   ++i) { crab::outs() << "--"; } crab::outs()
                                  << "| depth=" << m_depth
                                  << " call=" << cs.get_func_name() << "\n";);
      analyze_callee(cs, callee);
      --m_depth;
    };
    
    if (!m_ctx.analyze_recursive_functions() ||
	m_ctx.get_widening_set().count(callee_cg_node) <= 0) {
      analyze_callee_with_depth_check(callee_cg_node);
    } else {
      // This is recursive call.
      // 
      // TODO: use m_depth to break unexpected infinite loops inside
      // recursive functions
      analyze_callee(cs, callee_cg_node);
    }
  }

  intra_analyzer_with_call_semantics_t &
  add_analyzer(cg_node_t cg_node,
               std::unique_ptr<intra_analyzer_with_call_semantics_t> analysis) {
    auto res =
        m_intra_analyzer.insert(std::make_pair(cg_node, std::move(analysis)));
    return *((res.first)->second);
  }

  bool has_analyzer(cg_node_t cg_node) const {
    return m_intra_analyzer.find(cg_node) != m_intra_analyzer.end();
  }

  intra_analyzer_with_call_semantics_t &get_analyzer(cg_node_t cg_node) {
    auto it = m_intra_analyzer.find(cg_node);
    if (it == m_intra_analyzer.end()) {
      CRAB_ERROR("cannot find analysis for ",
                 cg_node.get_cfg().get_func_decl().get_func_name());
    }
    return *(it->second);
  }
};
} // end namespace top_down_inter_impl
} // end namespace analyzer
} // end namespace crab

namespace crab {
namespace analyzer {

/* The top-down inter-procedural analysis */
template <typename CallGraph, typename AbsDom>
class top_down_inter_analyzer:
    public inter_analyzer_api<CallGraph, AbsDom, AbsDom> {
  using cg_node_t = typename CallGraph::node_t;
  using this_type = top_down_inter_analyzer<CallGraph, AbsDom>;

public:
  using abs_dom_t = AbsDom;
  using cfg_t = typename cg_node_t::cfg_t;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;
  using varname_t = typename cfg_t::varname_t;
  using number_t = typename cfg_t::number_t;
  using variable_t = typename cfg_t::variable_t;
  using liveness_t = live_and_dead_analysis<cfg_t>;
  using liveness_map_t = std::unordered_map<cfg_t, const liveness_t *>;
  using wto_t = ikos::wto<cfg_t>;
  using wto_map_t = std::unordered_map<cfg_t, const wto_t *>;
  using params_t = inter_analyzer_parameters<CallGraph>;
  using checks_db_t = checker::checks_db;
  using summary_t = typename inter_analyzer_api<CallGraph,AbsDom,AbsDom>::summary_t;
  
private:
  // abstract transformer and analysis
  using td_inter_abs_tr_t =
      top_down_inter_impl::top_down_inter_transformer<CallGraph, abs_dom_t>;
  using intra_analyzer_with_call_semantics_t =
      analyzer_internal_impl::fwd_analyzer<cfg_t, td_inter_abs_tr_t>;
  using invariant_map_t =
      typename intra_analyzer_with_call_semantics_t::invariant_map_t;

  // global context stuff
  using global_context_t =
      top_down_inter_impl::global_context<CallGraph, abs_dom_t,
                                          invariant_map_t>;
  using calling_context_t = typename global_context_t::calling_context_t;
  using calling_context_ptr = typename global_context_t::calling_context_ptr;
  using calling_context_table_t =
      typename global_context_t::calling_context_table_t;
  using calling_context_collection_t =
      typename global_context_t::calling_context_collection_t;
  using global_invariant_map_t =
      typename global_context_t::global_invariant_map_t;

  // detect widening points in the call graph
  using wto_cg_t = typename global_context_t::wto_cg_t;
  struct widening_set_builder : public wto_cg_t::wto_component_visitor_t {
    using wto_vertex_t = typename wto_cg_t::wto_vertex_t;
    using wto_cycle_t = typename wto_cg_t::wto_cycle_t;
    using widening_set_t = std::set<typename CallGraph::node_t>;
    widening_set_t &m_widening_set;
    widening_set_builder(widening_set_t &widening_set)
        : m_widening_set(widening_set) {}
    virtual void visit(wto_cycle_t &cycle) override {
      m_widening_set.insert(cycle.head());
      for (auto &wto_component : cycle) {
        wto_component.accept(this);
      }
    }
    virtual void visit(wto_vertex_t &vertex) override {
    }
  };

  inline AbsDom make_bottom() const {
    auto const &dom = m_abs_tr->get_abs_value();
    return dom.make_bottom();
  }

  AbsDom get_invariant(const global_invariant_map_t &global_map,
                       const cfg_t &cfg, basic_block_label_t bb) const {
    auto g_it = global_map.find(cfg);
    if (g_it == global_map.end()) {
      return make_bottom(); // dead function
    }
    const invariant_map_t &m = g_it->second;
    auto it = m.find(bb);
    if (it == m.end()) {
      return make_bottom(); // dead block
    }
    return it->second;
  }

  CallGraph &m_cg;
  global_context_t m_ctx;
  std::unique_ptr<td_inter_abs_tr_t> m_abs_tr;

public:
  top_down_inter_analyzer(CallGraph &cg, abs_dom_t init,
                          const params_t &params = params_t())
      : m_cg(cg),
        m_ctx(params.live_map, params.wto_map, params.run_checker,
              params.checker_verbosity, params.keep_cc_invariants,
              params.keep_invariants, params.max_call_contexts,
              params.analyze_recursive_functions, params.exact_summary_reuse,
              params.only_main_as_entry, params.widening_delay,
              params.descending_iters, params.thresholds_size),
        m_abs_tr(new td_inter_abs_tr_t(m_cg, m_ctx, std::move(init))) {
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking call graph ... ";);
    cg.type_check();
    CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);
  }

  top_down_inter_analyzer(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  /** ===== Run the inter-procedural analysis
   *
   * The top-down analysis runs multiple analyses, one per callgraph
   * entry, starting with init.
   **/
  void run(abs_dom_t init) override {
    CRAB_SCOPED_TIMER_STATS(TimerInter, 1)
    CRAB_VERBOSE_IF(
        1, get_msg_stream()
               << "Computing weak topological ordering of the callgraph ... \n";);

    // Compute all the widening points in the callgraph
    auto &widening_set = m_ctx.get_widening_set();

    std::vector<typename CallGraph::node_t> entries = m_cg.entries();
    if (entries.empty()) {
      // HACK for callgraph limitation: only nodes without
      // predecessors are considered entries.
      auto p = m_cg.nodes();
      entries.insert(entries.end(), p.first, p.second);
    }

    auto &wto_cg_map = m_ctx.get_wto_cg_map();
    for (auto entry : entries) {
      std::unique_ptr<wto_cg_t> wto_cg(new wto_cg_t(m_cg, entry));
      widening_set_builder widen_builder(widening_set);
      wto_cg->accept(&widen_builder);
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Call graph WTO for entry "
                                          << entry << "=" << *wto_cg << "\n";);
      wto_cg_map[entry] = std::move(wto_cg);
    }

    CRAB_LOG(
        "inter", crab::outs() << "Widening points={"; for (auto &cg_node
                                                           : widening_set) {
          crab::outs() << cg_node.name() << ";";
        } crab::outs() << "}\n";);
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Done.\n";);

    top_down_inter_impl::extract_callsite_info<CallGraph, global_context_t>(
        m_cg, m_ctx);

    if (entries.empty()) {
      CRAB_WARN("Found no entry points in the call graph.");
    } else {
      CRAB_VERBOSE_IF(1,
		      if (!m_ctx.only_main_as_entry()) {
			get_msg_stream()
			  << "Started inter-procedural analysis considering all entry points.\n";
		      } else {
			get_msg_stream()
			  << "Started inter-procedural analysis *only* from main.\n";
		      });
    
      for (auto cg_node : entries) {
        if (m_ctx.only_main_as_entry()) {
          if (cg_node.name() != "main") {
	    CRAB_VERBOSE_IF(1, get_msg_stream()
			    << "Skipped analysis of call graph WTO starting from "
			    << cg_node.name() << "\n";);
            continue;
          }
        }

        abs_dom_t dom(init);
        m_abs_tr->set_abs_value(std::move(dom));
        m_ctx.get_call_stack().push_back(cg_node);
        intra_analyzer_with_call_semantics_t *entry_analysis =
            top_down_inter_impl::analyze_function<
                cg_node_t, intra_analyzer_with_call_semantics_t>(cg_node,
                                                                 *m_abs_tr, 0);
        assert(entry_analysis);
        top_down_inter_impl::check_function(cg_node, *entry_analysis, m_ctx);
        m_ctx.get_call_stack().pop_back();
        entry_analysis->clear();
      }
      CRAB_VERBOSE_IF(1, get_msg_stream()
                             << "Finished inter-procedural analysis\n";);

      if (::crab::CrabSanityCheckFlag) {
        if (!m_ctx.get_call_stack().empty()) {
          CRAB_ERROR("Something is wrong with the call stack");
        }
      }
    }
  }

  // return the analyzed call graph
  CallGraph &get_call_graph() override { return m_cg; }

  // Return the *context-insensitive* invariants that hold at the entry of b in
  // cfg.
  AbsDom get_pre(const cfg_t &cfg, const basic_block_label_t &bb) const override {
    return get_invariant(m_ctx.get_global_pre_invariants(), cfg, bb);
  }

  // Return the *context-insensitive* invariants that hold at the exit of b in
  // cfg.
  AbsDom get_post(const cfg_t &cfg, const basic_block_label_t &bb) const override {
    return get_invariant(m_ctx.get_global_post_invariants(), cfg, bb);
  }

  void clear() override {
    CRAB_WARN(
        "clear operation of the inter-procedural analysis not implemented yet");
  }

  // TODO: caching
  summary_t get_summary(const cfg_t &cfg) const override {
    summary_t summary(cfg.get_func_decl());
    auto it = m_ctx.get_calling_context_table().find(cfg);
    if (it != m_ctx.get_calling_context_table().end()) {
      auto &ccs = it->second;
      for (auto &cc:  ccs) {
	summary.add(cc->get_pre_summary(), cc->get_post_summary());
      }
    }
    return summary;
  }
  
  const checks_db_t get_all_checks() const { return m_ctx.get_checks_db(); }

  checks_db_t get_all_checks() { return m_ctx.get_checks_db(); }

  void print_checks(crab::crab_os &o) const {
    if (m_ctx.run_checker()) {
      m_ctx.get_checks_db().write(o);
    }
  }
};
} // namespace analyzer
} // namespace crab
