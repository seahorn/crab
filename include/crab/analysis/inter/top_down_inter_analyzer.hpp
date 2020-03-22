#pragma once

/*
   Standard top-down interprocedural analysis with memoization.
   TODO: precise analysis of recursive functions.
*/

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/cg/cg_bgl.hpp> // for wto of callgraphs
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/iterators/wto.hpp>

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
namespace top_down_inter_impl {

/**
 *  This class represents the calling context of a function F. A
 *  calling context C consists of a summary (input-output pair of
 *  abstract states) that summarizes the analysis of F for a
 *  particular callsite. The calling context can also store optionally
 *  all the invariants inferred during the analysis of F.
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

  // A summary is a pair (I,O) of abstract states. I's variables
  // should only contain the input parameters of the function. I is
  // typically used to check if a summary can be reused. O's variables
  // should contain both input and output parameters of a function. We
  // also include the inputs in O's variables so that the underlying
  // abstract domain can keep track of relationships between input and
  // outputs.

  fdecl_t m_fdecl;
  abs_dom_t m_input;  // must be projected on m_fdecl inputs
  abs_dom_t m_output; // must be projected on m_fdecl inputs and outputs

  // invariants that hold at the entry of each function's block
  invariant_map_t m_pre_invariants;
  // invariants that hold at the exit of each function's block
  invariant_map_t m_post_invariants;

  // if true then the calling context has not been joined yet with
  // other contexts.
  bool m_exact;

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
    return std::move(out);
  }

  // Private constructor used to join calling contexts
  calling_context(fdecl_t fdecl, abs_dom_t input, abs_dom_t output,
                  invariant_map_t &&pre_invariants,
                  invariant_map_t &&post_invariants)
      : m_fdecl(fdecl), m_input(input), m_output(output),
        m_pre_invariants(std::move(pre_invariants)),
        m_post_invariants(std::move(post_invariants)), m_exact(false) {}

public:
  calling_context(fdecl_t fdecl, abs_dom_t input, abs_dom_t output,
                  // if true then the calling context maintains only
                  // invariants for the entry and exit blocks.
                  bool minimize_invariants, invariant_map_t &&pre_invariants,
                  invariant_map_t &&post_invariants)
      : m_fdecl(fdecl), m_input(input), m_output(output),
        m_pre_invariants(std::move(pre_invariants)),
        m_post_invariants(std::move(post_invariants)), m_exact(true) {

    if (minimize_invariants) {
      // XXX: don't remove the invariants, just make them top.
      for (auto &kv : m_pre_invariants) {
        kv.second = abs_dom_t::top();
      }
      for (auto &kv : m_post_invariants) {
        kv.second = abs_dom_t::top();
      }
    }
  }

  calling_context(const calling_context_t &o) = delete;

  calling_context_t &operator=(const calling_context_t &o) = delete;

  // A summary is a pair (I,O) of abstract states. By construction, I
  // must be projected only on function input parameters and O must be
  // projected on function's input and output parameters.

  fdecl_t get_fdecl() const { return m_fdecl; }

  // return I from the pair (I,O)
  abs_dom_t get_input() const { return m_input; }

  // return O from the pair (I,O)
  abs_dom_t get_output() const { return m_output; }

  // Check if d entails the summary input
  bool is_subsumed(abs_dom_t d) const {
    // XXX: we choose for now to use equality for reusing the
    // strongest possible summaries. But we might want to use
    // inclusion by default so that we increase the chance of reusing
    // summaries.
    abs_dom_t input = get_input();
    if (m_exact) {
      return (d <= input && input <= d);
    } else {
      return (d <= input);
    }
  }

  // Join this and other
  std::unique_ptr<calling_context_t> join_with(calling_context_t &other) {
    if (!(get_fdecl() == other.get_fdecl())) {
      CRAB_ERROR("Cannot join calling contexts because they are from different "
                 "functions");
    }

    return std::unique_ptr<calling_context_t>(new calling_context_t(
        m_fdecl, m_input | other.get_input(), m_output | other.get_output(),
        std::move(join(m_pre_invariants, other.m_pre_invariants)),
        std::move(join(m_post_invariants, other.m_post_invariants))));
  }

  // invariants that hold at the entry of basic block bb
  abs_dom_t get_pre(basic_block_label_t bb) const {
    auto it = m_pre_invariants.find(bb);
    if (it == m_pre_invariants.end()) {
      return abs_dom_t::bottom(); // dead block under particular context
    } else {
      return it->second;
    }
  }

  // invariants that hold at the exit of basic block bb
  abs_dom_t get_post(basic_block_label_t bb) const {
    auto it = m_post_invariants.find(bb);
    if (it == m_post_invariants.end()) {
      return abs_dom_t::bottom(); // dead block under particular context
    } else {
      return it->second;
    }
  }

  void write(crab_os &o) const {
    abs_dom_t input = get_input();
    abs_dom_t output = get_output();
    o << "SUMMARY"
      << "\n\tI=" << input << "\n\tO=" << output << "\n";
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
    crab::ScopedCrabStats __st__("Interprocedural.add_context");

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
      crab::CrabStats::count(
          "Interprocedural.num_max_calling_contexts_reached");
      CRAB_LOG("inter",
               crab::outs()
                   << "[INTER] joining two oldest calling contexts\n";);
    }
    ccs.push_back(std::move(cc));

    if (compress) {
      // -- remove redundant contexts
      assert(!ccs.empty());
      calling_context_ptr_deque new_ccs;
      abs_dom_t joined_input = ccs.front()->get_input();
      abs_dom_t joined_output = ccs.front()->get_output();
      auto it = ccs.begin();
      // the first one is the joined calling context so we keep it.
      new_ccs.push_back(std::move(*it));
      ++it;
      for (auto et = ccs.end(); it != et; ++it) {
        // discard any input-output pair that is subsumed by the joined
        // calling context
        if (!((*it)->get_input() <= joined_input &&
              (*it)->get_output() <= joined_output)) {
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

template <typename CallGraphNode, typename AbsDom, typename InvariantMap>
class global_context {
  using this_type = global_context<CallGraphNode, AbsDom, InvariantMap>;

public:
  using cfg_t = typename CallGraphNode::cfg_t;
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
  using liveness_map_t = std::unordered_map<cfg_t, const liveness<cfg_t> *>;
  using wto_cfg_map_t = std::unordered_map<cfg_t, const ikos::wto<cfg_t> *>;
  using checks_db_t = checker::checks_db;
  using global_invariant_map_t = std::unordered_map<cfg_t, invariant_map_t>;

private:
  // -- recursive functions
  std::set<CallGraphNode> m_widening_set;
  // -- to break cycles
  std::vector<CallGraphNode> m_call_stack;
  // -- liveness symbols for each function
  const liveness_map_t *m_live_map;
  // -- wto for each function
  const wto_cfg_map_t *m_wto_cfg_map; // to avoid recomputing wto of cfgs
  // -- context-insensitive invariants (for external queries)
  global_invariant_map_t m_pre_invariants;
  global_invariant_map_t m_post_invariants;
  // -- the policy for making tractable the number of calling contexts
  std::unique_ptr<context_sensitivity_policy_t> m_cs_policy;
  // -- enable checking interleaved with analysis
  bool m_enable_checker;
  checks_db_t m_checks_db;
  unsigned m_checker_verbosity;
  bool m_is_checking_phase;
  // -- all calling contexts
  calling_context_table_t m_cc_table;
  // -- minimize the number of invariants per calling context
  bool m_minimize_invariants;
  // -- max number of calling contexts
  unsigned int m_max_call_contexts;
  // -- fixpoint parameters
  unsigned int m_widening_delay;
  unsigned int m_descending_iters;
  unsigned int m_thresholds_size;

  void join_with(global_invariant_map_t &global_table, cfg_t cfg,
                 invariant_map_t &other) {
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
      global_table.insert({cfg, std::move(invariants)});
    }
  }

public:
  global_context(const liveness_map_t *live_map,
                 const wto_cfg_map_t *wto_cfg_map, bool enable_checker,
                 unsigned checker_verbosity, bool minimize_invariants,
                 unsigned int max_call_contexts, unsigned int widening_delay,
                 unsigned int descending_iters, unsigned int thresholds_size)
      : m_live_map(live_map), m_wto_cfg_map(wto_cfg_map),
        m_cs_policy(
            new default_context_sensitivity_policy_t(max_call_contexts)),
        m_enable_checker(enable_checker),
        m_checker_verbosity(checker_verbosity), m_is_checking_phase(false),
        m_minimize_invariants(minimize_invariants),
        m_max_call_contexts(max_call_contexts),
        m_widening_delay(widening_delay), m_descending_iters(descending_iters),
        m_thresholds_size(thresholds_size) {}

  global_context(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  std::set<CallGraphNode> &get_widening_set() { return m_widening_set; }

  const std::set<CallGraphNode> &get_widening_set() const {
    return m_widening_set;
  }

  std::vector<CallGraphNode> &get_call_stack() { return m_call_stack; }

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

  bool minimize_invariants() const { return m_minimize_invariants; }

  unsigned int get_max_call_contexts() const { return m_max_call_contexts; }

  unsigned int get_widening_delay() const { return m_widening_delay; }

  unsigned int get_descending_iters() const { return m_descending_iters; }

  unsigned int get_thresholds_size() const { return m_thresholds_size; }

  // context-insensitive invariants for each function

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

  void join_invariants_with(CallGraphNode cg_node,
                            invariant_map_t &pre_invariants,
                            invariant_map_t &post_invariants) {
    join_with(get_global_pre_invariants(), cg_node.get_cfg(), pre_invariants);
    join_with(get_global_post_invariants(), cg_node.get_cfg(), post_invariants);
  }
};

// Wrapper to call the intra-procedural analysis with inter-procedural
// semantics for call and return statements.
template <typename CallGraphNode, typename IntraCallSemAnalyzer>
std::shared_ptr<IntraCallSemAnalyzer> get_inter_analysis(
    CallGraphNode cg_node,
    std::shared_ptr<typename IntraCallSemAnalyzer::abs_tr_t> abs_tr) {

  using cfg_t = typename CallGraphNode::cfg_t;
  // -- intra-procedural checker stuff
  using intra_checker_t = checker::intra_checker<IntraCallSemAnalyzer>;
  using assertion_checker_t =
      checker::assert_property_checker<IntraCallSemAnalyzer>;
  using checks_db_t = checker::checks_db;

  cfg_t cfg = cg_node.get_cfg();
  assert(abs_tr);
  assert(cfg.has_func_decl());

  CRAB_VERBOSE_IF(1, get_msg_stream()
                         << "++ Analyzing function  "
                         << cfg.get_func_decl().get_func_name() << "\n";);

  auto &ctx = abs_tr->get_context();
  std::shared_ptr<IntraCallSemAnalyzer> analyzer = nullptr;

  if (abs_tr->has_analyzer(cg_node)) {
    /// --- 1. reuse the intra-analyzer if already exists
    // don't create another analyzer if we already created one
    analyzer = abs_tr->get_analyzer(cg_node);
    // make sure no results from previous run
    analyzer->reset();
  } else {
    // get liveness symbols for cfg if available
    const liveness<cfg_t> *live = nullptr;
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
    analyzer = std::make_shared<IntraCallSemAnalyzer>(
        cfg, wto, abs_tr->get_shared_ptr(), live, ctx.get_widening_delay(),
        ctx.get_descending_iters(), ctx.get_thresholds_size());
    abs_tr->add_analyzer(cg_node, analyzer);
  }

  /// -- 3. Run the analyzer
  analyzer->run_forward();

  /// -- 4. Store the invariants (merging with other contexts)
  ctx.join_invariants_with(cg_node, analyzer->get_pre_invariants(),
                           analyzer->get_post_invariants());

  /// -- 5. Optionally, run the checker after the analysis. Here we
  /// run the checker with the specific invariants for the particular
  /// context. We can also do the checking once the analysis phase
  /// finishes but we would use the context-insensitive invariants
  /// stored in step 4.

  if (ctx.run_checker()) {
    bool &is_checking_phase = ctx.get_is_checking_phase();
    assert(!is_checking_phase);
    is_checking_phase = true;
    typename intra_checker_t::prop_checker_ptr p(
        new assertion_checker_t(ctx.get_checker_verbosity()));
    intra_checker_t checker(*analyzer, {p});
    checker.run();
    ctx.get_checks_db() += checker.get_all_checks();
    assert(is_checking_phase);
    is_checking_phase = false;
  }

  return analyzer;
}

/*
   All the state of the inter-procedural analysis is kept in this
   transformer. This allows us to use the intra-procedural analysis
   fully as a black box.
*/
template <typename CallGraph, typename AbsDom>
class top_down_inter_transformer final
    : public intra_abs_transformer<AbsDom>,
      public std::enable_shared_from_this<
          top_down_inter_transformer<CallGraph, AbsDom>> {

  using this_type = top_down_inter_transformer<CallGraph, AbsDom>;

  using cg_node_t = typename CallGraph::node_t;
  using cfg_t = typename cg_node_t::cfg_t;
  using fdecl_t = typename cfg_t::fdecl_t;
  using variable_t = typename cfg_t::variable_t;
  using intra_abs_transformer_t = intra_abs_transformer<AbsDom>;
  using abs_transform_api_t =
      typename intra_abs_transformer_t::abs_transform_api_t;
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
      top_down_inter_impl::global_context<cg_node_t, abs_dom_t,
                                          invariant_map_t>;
  // -- calling context stuff
  using calling_context_t = typename global_context_t::calling_context_t;
  using calling_context_ptr = typename global_context_t::calling_context_ptr;
  using calling_context_collection_t =
      typename global_context_t::calling_context_collection_t;
  using calling_context_table_t =
      typename global_context_t::calling_context_table_t;

private:
  // -- the callgraph
  CallGraph m_cg;
  // -- global parameters of the analysis
  global_context_t &m_ctx;
  // -- to avoid starting from scratch an analyzer:
  //    saving wto computation, type checking, etc.
  // XXX: it cannot be in m_ctx because of cyclic dependencies.
  std::unordered_map<cg_node_t,
                     std::shared_ptr<intra_analyzer_with_call_semantics_t>>
      m_intra_analyzer;

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
    // crab::CrabStats::count("Interprocedural.num_calling_contexts");
  }

  static void get_fdecl_parameters(fdecl_t &fdecl,
                                   std::vector<variable_t> &out) {
    out.reserve(fdecl.get_num_inputs() + fdecl.get_num_outputs());
    out.insert(out.end(), fdecl.get_inputs().begin(), fdecl.get_inputs().end());
    out.insert(out.end(), fdecl.get_outputs().begin(),
               fdecl.get_outputs().end());
  }

  /**
   *  Restrict operation.
   * 
   * Produce a new abstract state at the callee for the entry block.
   *  - caller_dom: abstract state at the caller before the call.
   *  - callee_dom: initial state at the callee (by default top)
   **/
  static AbsDom get_callee_entry(callsite_t &cs, fdecl_t &fdecl,
                                 AbsDom caller_dom,
                                 AbsDom callee_dom = AbsDom::top()) {
    crab::ScopedCrabStats __st__("Interprocedural.get_callee_entry");

    // caller_dom.normalize();
    if (caller_dom.is_bottom()) {
      return caller_dom;
    }
    CRAB_LOG("inter-restrict",
             errs() << "Inv at the caller: " << caller_dom << "\n");
    // propagate from actual to formal parameters
    CRAB_LOG("inter-restrict", errs()
                                 << "Unifying formal and actual parameters\n";);
    for (unsigned i = 0, e = fdecl.get_inputs().size(); i < e; ++i) {
      variable_t formal = fdecl.get_inputs()[i];
      variable_t actual = cs.get_args()[i];
      if (!(formal == actual)) {
        CRAB_LOG("inter-restrict",
                 errs() << "\t" << formal << " and " << actual << "\n";);
        inter_transformer_helpers<AbsDom>::unify(caller_dom, formal, actual);
      }
    }
    CRAB_LOG("inter-restrict", errs() << "Inv after formal/actual unification: "
                                    << caller_dom << "\n";);
    // Meet
    callee_dom = caller_dom & callee_dom;
    CRAB_LOG("inter-restrict",
             errs() << "Inv after meet with callee  " << callee_dom << "\n";);
    // project onto **input** formal parameters
    callee_dom.project(fdecl.get_inputs());
    CRAB_LOG("inter-restrict",
             errs() << "Inv at the callee after projecting onto formals: ";
             for (auto &v
                  : fdecl.get_inputs()) { errs() << v << ";"; } errs()
             << "\n"
             << callee_dom << "\n";);

    return callee_dom;
  }

  static std::vector<variable_t> set_difference(std::vector<variable_t> v1,
                                                std::vector<variable_t> v2) {
    std::vector<variable_t> out;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),
                        std::inserter(out, out.end()));
    return out;
  }

  /**
   * Extend operation
   *
   * Produce a new abstract state at the caller for the call
   * continuation.
   * - caller_dom: abstract state at the caller before the call
   * - sum_out_dom: abstract state at the exit block of the callee
   *   but already projected onto formal parameters of the function.
   * - sum_out_variables is fdecl.get_inputs() U
   *   fdecl.get_outputs(). We pass them explicitly to avoid to
   *   concatenate them again here.
   * 
   * This code should work even if callsite lhs variables and actual
   * parameters are not disjoint.
   **/
  static AbsDom
  get_caller_continuation(callsite_t &cs, fdecl_t &fdecl, AbsDom caller_dom,
                          const std::vector<variable_t> &sum_out_variables,
                          AbsDom sum_out_dom) {
    crab::ScopedCrabStats __st__("Interprocedural.get_caller_continuation");

    if (caller_dom.is_bottom()) {
      return caller_dom;
    }

    CRAB_LOG("inter-extend",
	     crab::outs() << "Caller before " << cs << "=" << caller_dom << "\n";);
    
    // make sure **output** actual parameters are unconstrained
    caller_dom.forget(cs.get_lhs());

    CRAB_LOG("inter-extend",
	     crab::outs() << "Caller after forgetting lhs variables=" << caller_dom << "\n";);
    
    // Meet
    caller_dom = caller_dom & sum_out_dom;

    CRAB_LOG("inter-extend",
	     crab::outs() << "Caller after meet with callee's summary=" << caller_dom << "\n";);
    
    // propagate from callee's outputs to caller's lhs of the callsite
    for (unsigned i = 0, e = fdecl.get_outputs().size(); i < e; ++i) {
      variable_t formal = fdecl.get_outputs()[i];
      variable_t actual = cs.get_lhs()[i];
      if (!(formal == actual)) {
	CRAB_LOG("inter-extend",
		 crab::outs() << "Unifying " << formal << " and " << actual << "\n";);
        inter_transformer_helpers<AbsDom>::unify(caller_dom, actual, formal);
      }
    }

    CRAB_LOG("inter-extend",
	     crab::outs() << "Caller after unifying formal/actual paramters=" << caller_dom << "\n";);
    
    // Remove the callee variables from the caller continuation
    // XXX: don't forget a variable if it appears on cs.get_lhs()
    caller_dom.forget(set_difference(sum_out_variables, cs.get_lhs()));

    return caller_dom;
  }

  /* Analysis of a callsite */
  void analyze_callee(callsite_t &cs, cg_node_t callee_cg_node) {

    // 1. Get CFG from the callee
    cfg_t callee_cfg = callee_cg_node.get_cfg();

    assert(callee_cfg.has_func_decl());
    fdecl_t fdecl = callee_cfg.get_func_decl();

    // 2. Generate initial abstract state for the callee
    AbsDom caller_dom(this->get_abs_value());

    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = caller_dom.is_bottom();
    }

    AbsDom callee_entry;
    // If we start the analysis of a recursive procedure we must do it
    // without propagating from caller to callee
    if (m_ctx.get_widening_set().count(callee_cg_node) <= 0) {
      callee_entry = get_callee_entry(cs, fdecl, caller_dom);
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = callee_entry.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR(
            "Invariant became bottom after propagating from actual to formals ",
            cs);
      }
    }

    AbsDom callee_exit;
    std::vector<variable_t> callee_exit_vars;
    get_fdecl_parameters(fdecl, callee_exit_vars);

    // 3. Check if the same call context has been seen already
    bool call_context_already_seen = false;
    auto it = m_ctx.get_calling_context_table().find(callee_cfg);
    if (it != m_ctx.get_calling_context_table().end()) {
      auto &call_contexts = it->second;
      CRAB_LOG("inter", if (call_contexts.empty()) {
        crab::outs() << "There is no call contexts stored for " << cs << "\n";
      });
      for (unsigned i = 0, e = call_contexts.size(); i < e; ++i) {
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Checking at " << cs << " if\n"
                              << callee_entry << "\nis subsumed by\n";
                 call_contexts[i]->write(crab::outs()); crab::outs() << "\n";);
        if (call_contexts[i]->is_subsumed(callee_entry)) {
          CRAB_LOG("inter", crab::outs()
                                << "[INTER] Reusing call context at \"" << cs
                                << "\" with entry=" << callee_entry << "\n";);
          callee_exit = call_contexts[i]->get_output();
          call_context_already_seen = true;
          break;
        }
      }
    }

    if (call_context_already_seen) {
      if (!m_ctx.get_is_checking_phase()) {
        crab::CrabStats::count("Interprocedural.num_reused_callsites");
      }
    } else {
      if (m_ctx.get_is_checking_phase()) {
        CRAB_ERROR("in checking phase we should not analyze the callsite ", cs);
      }
      crab::CrabStats::count("Interprocedural.num_analyzed_callsites");
      // 4. Run intra analyzer
      CRAB_LOG("inter", crab::outs()
                            << "[INTER] Started \"" << cs
                            << "\" with entry=" << callee_entry << "\n";);

      abs_dom_t callee_init(callee_entry);
      this->set_abs_value(std::move(callee_init));
      std::shared_ptr<intra_analyzer_with_call_semantics_t> callee_analysis =
          top_down_inter_impl::get_inter_analysis<
              typename CallGraph::node_t, intra_analyzer_with_call_semantics_t>(
              m_cg.get_callee(cs), get_shared_ptr());

      if (callee_cfg.has_exit()) {
        callee_exit = callee_analysis->get_post(callee_cfg.exit());
        CRAB_LOG("inter", crab::outs()
                              << "[INTER] Finished \"" << cs
                              << "\" with exit=" << callee_exit << "\n";);

      } else {
        // if the callee has not exit is because it's a noreturn function.
        callee_exit = AbsDom::bottom();
        CRAB_LOG("inter", crab::outs() << "[INTER] Finished \"" << cs
                                       << "\" with exit=" << callee_exit
                                       << ": the callee has no exit block.\n";);
      }

      CRAB_LOG("inter2", callee_analysis->write(crab::outs());
               crab::outs() << "\n";);

      // project callee_exit onto function input-output parameters
      callee_exit.project(callee_exit_vars);

      // 5. Add the new calling context
      auto &pre_invariants = callee_analysis->get_pre_invariants();
      auto &post_invariants = callee_analysis->get_post_invariants();
      calling_context_ptr cc(new calling_context_t(
          fdecl, callee_entry, callee_exit, m_ctx.minimize_invariants(),
          std::move(pre_invariants), std::move(post_invariants)));
      add_calling_context(callee_cfg, std::move(cc));
      callee_analysis->reset();
    }

    pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = callee_exit.is_bottom();
    }

    // 6. Generate abstract state for the continuation at the caller
    AbsDom caller_cont_dom = get_caller_continuation(
        cs, fdecl, caller_dom, callee_exit_vars, callee_exit);

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = caller_cont_dom.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after propagating from actuals to "
                   "formals ",
                   cs);
      }
    }

    CRAB_LOG("inter", crab::outs() << "[INTER] Continuation for \"" << cs
                                   << "\"=" << caller_cont_dom << "\n";);
    this->set_abs_value(std::move(caller_cont_dom));
  }

public:
  top_down_inter_transformer(CallGraph cg, global_context_t &ctx, AbsDom init)
      : intra_abs_transformer_t(init), m_cg(cg), m_ctx(ctx) {}

  top_down_inter_transformer(const this_type &&o)
      : intra_abs_transformer_t(std::move(o.get_abs_value())),
        m_cg(std::move(o.m_cg)), m_ctx(std::move(o.m_ctx)) {}

  top_down_inter_transformer(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  std::shared_ptr<this_type> get_shared_ptr() {
    return this->shared_from_this();
  }

  const global_context_t &get_context() const { return m_ctx; }

  global_context_t &get_context() { return m_ctx; }

  virtual void exec(callsite_t &cs) override {
    if (!m_cg.has_callee(cs)) {
      CRAB_ERROR("Cannot find callee CFG for ", cs);
    }

    if (this->get_abs_value().is_bottom()) {
      return;
    }

    if (!m_ctx.get_is_checking_phase()) {
      crab::CrabStats::count("Interprocedural.num_callsites");
    }

    cg_node_t callee_cg_node = m_cg.get_callee(cs);
    auto &call_stack = m_ctx.get_call_stack();
    if (std::find(call_stack.begin(), call_stack.end(), callee_cg_node) ==
        call_stack.end()) {
      if (m_ctx.get_widening_set().count(callee_cg_node) > 0) {
        call_stack.push_back(callee_cg_node);
      }
      analyze_callee(cs, callee_cg_node);
      if (m_ctx.get_widening_set().count(callee_cg_node) > 0) {
        call_stack.pop_back();
      }
    } else {
      // we break the cycle by ignoring the effect of the recursive
      // call.
      cfg_t callee_cfg = callee_cg_node.get_cfg();
      assert(callee_cfg.has_func_decl());
      fdecl_t fdecl = callee_cfg.get_func_decl();
      AbsDom caller_dom(this->get_abs_value());
      if (!m_ctx.get_is_checking_phase()) {
        crab::CrabStats::count("Interprocedural.num_recursive_callsites");
      }
      CRAB_VERBOSE_IF(1, get_msg_stream()
                             << "++ Skipped analysis of recursive callee \""
                             << cs << "\"\n";);
      AbsDom callee_exit = AbsDom::top();
      std::vector<variable_t> callee_exit_vars;
      get_fdecl_parameters(fdecl, callee_exit_vars);
      AbsDom caller_cont_dom = get_caller_continuation(
          cs, fdecl, caller_dom, callee_exit_vars, callee_exit);
      CRAB_LOG("inter", crab::outs() << "[INTER] Continuation for \"" << cs
                                     << "\"=" << caller_cont_dom << "\n";);
      this->set_abs_value(std::move(caller_cont_dom));
    }
  }

  void
  add_analyzer(cg_node_t cg_node,
               std::shared_ptr<intra_analyzer_with_call_semantics_t> analysis) {
    m_intra_analyzer.insert({cg_node, analysis});
  }

  bool has_analyzer(cg_node_t cg_node) const {
    return m_intra_analyzer.find(cg_node) != m_intra_analyzer.end();
  }

  std::shared_ptr<intra_analyzer_with_call_semantics_t>
  get_analyzer(cg_node_t cg_node) {
    auto it = m_intra_analyzer.find(cg_node);
    if (it == m_intra_analyzer.end()) {
      CRAB_ERROR("cannot find analysis for ",
                 cg_node.get_cfg().get_func_decl().get_func_name());
    }
    return it->second;
  }
};
} // end namespace top_down_inter_impl
} // end namespace analyzer
} // end namespace crab

namespace crab {
namespace analyzer {

/* Customize parameters for the top-down inter-procedural analysis */
template <typename CallGraph> struct top_down_inter_analyzer_parameters {
  using cg_node_t = typename CallGraph::node_t;
  using cfg_t = typename cg_node_t::cfg_t;
  using liveness_t = liveness<cfg_t>;
  using liveness_map_t = std::unordered_map<cfg_t, const liveness_t *>;
  using wto_t = ikos::wto<cfg_t>;
  using wto_map_t = std::unordered_map<cfg_t, const wto_t *>;

  top_down_inter_analyzer_parameters()
      : run_checker(true), checker_verbosity(0), minimize_invariants(true),
        max_call_contexts(UINT_MAX), widening_delay(2),
        descending_iters(2 /*UINT_MAX*/), thresholds_size(20 /*0*/),
        live_map(nullptr), wto_map(nullptr) {}

  // whether interleave analysis with checking
  bool run_checker;
  unsigned checker_verbosity;
  // whether minimize the number of invariants stored per analyzed
  // function.
  bool minimize_invariants;
  // maximum number of calling contexts per callsite
  unsigned int max_call_contexts;
  // fixpoint parameters
  unsigned int widening_delay;
  unsigned int descending_iters;
  unsigned int thresholds_size;
  // take liveness symbols if already available
  const liveness_map_t *live_map;
  // take wto's if already available
  const wto_map_t *wto_map;
};

/* The top-down inter-procedural analysis */
template <typename CallGraph, typename AbsDom> class top_down_inter_analyzer {
  using cg_node_t = typename CallGraph::node_t;
  using this_type = top_down_inter_analyzer<CallGraph, AbsDom>;

public:
  using abs_dom_t = AbsDom;

  using cfg_t = typename cg_node_t::cfg_t;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;
  using varname_t = typename cfg_t::varname_t;
  using number_t = typename cfg_t::number_t;
  using variable_t = typename cfg_t::variable_t;

  using liveness_t = liveness<cfg_t>;
  using liveness_map_t = std::unordered_map<cfg_t, const liveness_t *>;

  using wto_t = ikos::wto<cfg_t>;
  using wto_map_t = std::unordered_map<cfg_t, const wto_t *>;

  using params_t = top_down_inter_analyzer_parameters<CallGraph>;

  using checks_db_t = checker::checks_db;

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
      top_down_inter_impl::global_context<typename CallGraph::node_t, abs_dom_t,
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
  using wto_cg_t = ikos::wto<CallGraph>;
  struct find_widening_set_visitor
      : public ikos::wto_component_visitor<CallGraph> {
    using wto_vertex_t = ikos::wto_vertex<CallGraph>;
    using wto_cycle_t = ikos::wto_cycle<CallGraph>;
    using widening_set_t = std::set<typename CallGraph::node_t>;
    widening_set_t &m_widening_set;
    find_widening_set_visitor(widening_set_t &widening_set)
        : m_widening_set(widening_set) {}
    void visit(wto_cycle_t &cycle) {
      m_widening_set.insert(cycle.head());
      for (auto &wto_component : cycle) {
        wto_component.accept(this);
      }
    }
    void visit(wto_vertex_t &vertex) {}
  };

  AbsDom get_invariant(const global_invariant_map_t &global_map,
                       const cfg_t &cfg, basic_block_label_t bb) const {
    auto g_it = global_map.find(cfg);
    if (g_it == global_map.end()) {
      return AbsDom::bottom(); // dead function
    }
    const invariant_map_t &m = g_it->second;
    auto it = m.find(bb);
    if (it == m.end()) {
      return AbsDom::bottom(); // dead block
    }
    return it->second;
  }

  CallGraph m_cg;
  global_context_t m_ctx;

public:
  top_down_inter_analyzer(CallGraph cg, const params_t &params = params_t())
      : m_cg(cg), m_ctx(params.live_map, params.wto_map, params.run_checker,
                        params.checker_verbosity, params.minimize_invariants,
                        params.max_call_contexts, params.widening_delay,
                        params.descending_iters, params.thresholds_size) {
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking call graph ... ";);
    crab::CrabStats::resume("CallGraph type checking");
    cg.type_check();
    crab::CrabStats::stop("CallGraph type checking");
    CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);
  }

  top_down_inter_analyzer(const this_type &o) = delete;

  this_type &operator=(const this_type &o) = delete;

  void run(abs_dom_t init = abs_dom_t::top()) {
    crab::ScopedCrabStats __st__("Inter");

    CRAB_VERBOSE_IF(
        1, get_msg_stream()
               << "Computing weak topological ordering of the callgraph ... ";);
    wto_cg_t wto_cg(m_cg);
    auto &widening_set = m_ctx.get_widening_set();
    find_widening_set_visitor widening_vis(widening_set);
    wto_cg.accept(&widening_vis);
    CRAB_LOG("inter", crab::outs() << "Call graph WTO=" << wto_cg << "\n";);
    CRAB_LOG(
        "inter", crab::outs() << "Widening points={"; for (auto &cg_node
                                                           : widening_set) {
          crab::outs() << cg_node.name() << ";";
        } crab::outs() << "}\n";);
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Done.";);

    std::vector<typename CallGraph::node_t> entries = m_cg.entries();
    if (entries.empty()) {
      CRAB_WARN("Found no entry points in the call graph.");
    } else {
      CRAB_VERBOSE_IF(1, get_msg_stream()
                             << "Started inter-procedural analysis\n";);
      auto inter_abs_tr =
          std::make_shared<td_inter_abs_tr_t>(m_cg, m_ctx, init);
      for (auto cg_node : entries) {
        if (widening_set.count(cg_node) > 0) {
          CRAB_ERROR("Entry point cannot be recursive");
        }
        std::shared_ptr<intra_analyzer_with_call_semantics_t> entry_analysis =
            top_down_inter_impl::get_inter_analysis<
                cg_node_t, intra_analyzer_with_call_semantics_t>(cg_node,
                                                                 inter_abs_tr);
        entry_analysis->reset();
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
  CallGraph get_call_graph() { return m_cg; }

  // Return the *context-insensitive* invariants that hold at the entry of b in
  // cfg.
  AbsDom get_pre(const cfg_t &cfg, basic_block_label_t bb) const {
    return get_invariant(m_ctx.get_global_pre_invariants(), cfg, bb);
  }

  // Return the *context-insensitive* invariants that hold at the exit of b in
  // cfg.
  AbsDom get_post(const cfg_t &cfg, basic_block_label_t bb) const {
    return get_invariant(m_ctx.get_global_post_invariants(), cfg, bb);
  }

  // TODO: return context-sensitive invariants if some client requires them.

  void clear() {
    CRAB_WARN(
        "clear operation of the inter-procedural analysis not implemented yet");
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
