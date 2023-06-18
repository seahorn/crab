#pragma once

/*
   A standard two-phase approach for summary-based,
   _context-insensitive_ inter-procedural analysis.
*/

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>
#include <crab/analysis/inter/inter_analyzer_api.hpp>
#include <crab/analysis/inter/inter_params.hpp>
#include <crab/cfg/cfg.hpp>   // callsite_or_fdecl wrapper
#include <crab/cg/cg_bgl.hpp> // for sccg.hpp
#include <crab/domains/generic_abstract_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace crab {
namespace analyzer {
namespace inter_analyzer_impl {

// The analysis supports two different domains for the bottom-up and
// top-down passes, respectively. At some point, we need to convert
// from one to the other.
template <typename Domain1, typename Domain2>
void convert_domains(Domain1 from, Domain2 &to);

/* Store the calling contexts of each function */
template <typename CFG, typename AbsDomain> class call_ctx_table {
public:
  using basic_block_t = typename CFG::basic_block_t;
  using callsite_t = typename basic_block_t::callsite_t;
  using fdecl_t = typename CFG::fdecl_t;
  using callsite_or_fdecl_t = crab::cfg::callsite_or_fdecl<CFG>;
  using varname_t = typename CFG::varname_t;
  using variable_t = typename CFG::variable_t;
  using abs_domain_t = AbsDomain;

private:
  using call_table_t = crab::cfg::callsite_or_fdecl_map<CFG, abs_domain_t>;
  call_table_t m_call_table;
  AbsDomain m_top;

  // XXX: assume context-insensitive analysis so it will merge all
  // calling contexts using abstract domain's join keeping a
  // single calling context per function.
  void insert_helper(callsite_or_fdecl_t key, AbsDomain inv) {
    auto it = m_call_table.find(key);
    if (it != m_call_table.end()) {
      it->second = it->second | inv;
    } else {
      m_call_table.insert(std::make_pair(key, inv));
    }
  }

public:
  call_ctx_table(AbsDomain top) : m_top(top) { m_top.set_to_top(); }

  call_ctx_table(const call_ctx_table<CFG, AbsDomain> &o) = delete;

  call_ctx_table<CFG, AbsDomain> &
  operator=(const call_ctx_table<CFG, AbsDomain> &o) = delete;

  void insert(const callsite_t &cs, AbsDomain inv) {
    insert_helper(&cs, inv);
  }

  void insert(const fdecl_t &d, AbsDomain inv) {
    insert_helper(&d, inv);
  }

  AbsDomain get_call_ctx(const fdecl_t &d) const {
    auto it = m_call_table.find(&d);
    if (it != m_call_table.end()) {
      return it->second;
    } else { 
      return m_top;
    }
  }
  
  void clear() {
    m_call_table.clear();
  }
};

// A summary is an input-output relationship between function
// parameters. The relationship can be as expressive as
// AbsDomain.
template <typename CFG, typename AbsDomain> class summary {

  using fdecl_t = typename CFG::fdecl_t;
  using abs_domain_t = AbsDomain;
  using variable_t = typename CFG::variable_t;
  using varname_t = typename variable_t::varname_t;

  // --- function info
  const fdecl_t &m_fdecl;
  // --- Summary involving only m_params + m_ret_vals variables
  abs_domain_t m_sum;
  // --- Keep all the input original parameters of the function
  std::vector<variable_t> m_inputs;
  // --- Keep a copy of all input original parameters of the function
  std::vector<variable_t> m_internal_inputs;
  // --- Keep all the output original parameters of the function
  std::vector<variable_t> m_outputs;
  // --- Keep a copy of all output original parameters of the function
  std::vector<variable_t> m_internal_outputs;

  // - m_sum is defined in terms of m_internal_inputs  and m_internal_outputs.
  // - m_fdecl is defined in terms of m_inputs and m_outputs.

  // helper to rename summaries
  void rename(abs_domain_t &abs, const std::vector<variable_t> &from_inputs,
              const std::vector<variable_t> &from_outputs,
              const std::vector<variable_t> &to_inputs,
              const std::vector<variable_t> &to_outputs) const {

    assert(from_inputs.size() == to_inputs.size());
    assert(from_outputs.size() == to_outputs.size());

    // append inputs and outputs
    std::vector<variable_t> from_vars, to_vars;
    from_vars.reserve(from_inputs.size() + from_outputs.size());
    from_vars.insert(from_vars.end(), from_inputs.begin(), from_inputs.end());
    from_vars.insert(from_vars.end(), from_outputs.begin(), from_outputs.end());
    to_vars.reserve(to_inputs.size() + to_outputs.size());
    to_vars.insert(to_vars.end(), to_inputs.begin(), to_inputs.end());
    to_vars.insert(to_vars.end(), to_outputs.begin(), to_outputs.end());
    abs.rename(from_vars, to_vars);
  }

public:
  summary(const fdecl_t &fdecl,
          // summary defined only in terms of inputs and outputs
          abs_domain_t sum, std::vector<variable_t> inputs,
          std::vector<variable_t> outputs)
      : m_fdecl(fdecl), m_sum(sum), m_inputs(inputs), m_outputs(outputs) {

    m_internal_inputs.reserve(m_inputs.size());
    m_internal_outputs.reserve(m_outputs.size());

    std::string prefix = "$";
    unsigned i = 0;
    for (auto v : m_inputs) {
      auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
      variable_t fresh_v(vfac.get(prefix + std::to_string(i)), v.get_type());
      m_internal_inputs.push_back(fresh_v);
      i++;
    }
    for (auto v : m_outputs) {
      auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
      variable_t fresh_v(vfac.get(prefix + std::to_string(i)), v.get_type());
      m_internal_outputs.push_back(fresh_v);
      i++;
    }

    if (m_fdecl.get_num_inputs() != m_inputs.size())
      CRAB_ERROR(
          "mismatch between function declaration and summary parameters");
    if (m_fdecl.get_num_outputs() != m_outputs.size())
      CRAB_ERROR(
          "mismatch between function declaration and summary return vals");
    if (m_inputs.size() != m_internal_inputs.size())
      CRAB_ERROR("internal error in summary class");
    if (m_outputs.size() != m_internal_outputs.size())
      CRAB_ERROR("internal error in summary class");
  }

  // -- The summary, input, and output variables contain the
  // -- original variable names.
  abs_domain_t get_sum() const { return m_sum; }

  abs_domain_t make_top() const { return m_sum.make_top(); }

  abs_domain_t make_bottom() const { return m_sum.make_bottom(); }

  // return the input variables of the summary
  const std::vector<variable_t> &get_inputs() const { return m_inputs; }

  // return the output variables of the summary
  const std::vector<variable_t> &get_outputs() const { return m_outputs; }

  // -- The summary, input, and output variables are renamed so
  //    that they have unique variable names. This avoids naming
  //    clashes when summaries are used in the interprocedural
  //    analysis.
  //
  abs_domain_t get_renamed_sum() const {
    abs_domain_t res(m_sum);
    rename(res, m_inputs, m_outputs, m_internal_inputs, m_internal_outputs);
    return res;
  }

  const std::vector<variable_t> &get_renamed_inputs() const {
    return m_internal_inputs;
  }

  const std::vector<variable_t> &get_renamed_outputs() const {
    return m_internal_outputs;
  }

  void write(crab_os &o) const {
    o << m_fdecl.get_func_name() << "(IN:{";
    for (unsigned i = 0; i < m_fdecl.get_num_inputs(); i++) {
      o << m_inputs[i] << ":" << m_fdecl.get_input_type(i);
      if (i != m_fdecl.get_num_inputs() - 1)
        o << ",";
    }
    o << "},OUT:{";
    for (unsigned i = 0; i < m_fdecl.get_num_outputs(); i++) {
      o << m_outputs[i] << ":" << m_fdecl.get_output_type(i);
      if (i != m_fdecl.get_num_outputs() - 1)
        o << ",";
    }
    o << "}) ==>\n";
    abs_domain_t tmp(m_sum);
    o << tmp;
  }

  friend crab_os &operator<<(crab_os &o, const summary<CFG, AbsDomain> &sum) {
    sum.write(o);
    return o;
  }
};

/* Store the summaries for each function*/
template <typename CFG, typename AbsDomain> class summary_table {

public:
  using summary_t = summary<CFG, AbsDomain>;
  using basic_block_t = typename CFG::basic_block_t;
  using callsite_t = typename basic_block_t::callsite_t;
  using fdecl_t = typename CFG::fdecl_t;
  using callsite_or_fdecl_t = crab::cfg::callsite_or_fdecl<CFG>;  
  using abs_domain_t = AbsDomain;
  using variable_t = typename CFG::variable_t;

private:
  using summary_table_t = crab::cfg::callsite_or_fdecl_map<CFG, summary_t>;
  summary_table_t m_sum_table;

public:
  summary_table() {}

  summary_table(const summary_table<CFG, AbsDomain> &o) = delete;

  summary_table<CFG, AbsDomain> &
  operator=(const summary_table<CFG, AbsDomain> &o) = delete;

  // insert summary information
  void insert(const fdecl_t &d, AbsDomain sum, std::vector<variable_t> inputs,
              std::vector<variable_t> outputs) {

    std::vector<variable_t> ins(inputs.begin(), inputs.end());
    std::vector<variable_t> outs(outputs.begin(), outputs.end());
    summary_t sum_tuple(d, sum, ins, outs);
    m_sum_table.insert(std::make_pair(callsite_or_fdecl_t(&d),std::move(sum_tuple)));
  }

  // return true if there is a summary
  bool has_summary(const callsite_t &cs) const {
    auto it = m_sum_table.find(&cs);
    return (it != m_sum_table.end());
  }

  bool has_summary(const fdecl_t &d) const {
    auto it = m_sum_table.find(&d);
    return (it != m_sum_table.end());
  }

  // get the summary
  summary_t get(const callsite_t &cs) const {
    auto it = m_sum_table.find(&cs);
    assert(it != m_sum_table.end());
    return (it->second);
  }

  summary_t get(const fdecl_t &d) const {
    auto it = m_sum_table.find(&d);
    assert(it != m_sum_table.end());
    return (it->second);
  }

  void clear() {
    m_sum_table.clear();
  }
  
  void write(crab_os &o) const {
    o << "--- Begin summary table ---\n";
    for (auto const &p : m_sum_table) {
      p.second.write(o);
      o << "\n";
    }
    o << "--- End summary table ---\n";
  }

  friend crab_os &operator<<(crab_os &o,
                             const summary_table<CFG, AbsDomain> &t) {
    t.write(o);
    return o;
  }
};

/**
 * Abstract transformer specialized for computing summaries.
 * class SumTable stores the summaries.
 **/
template <class SumTable>
class bu_summ_abs_transformer final
    : public intra_abs_transformer<typename SumTable::basic_block_t,
                                   typename SumTable::abs_domain_t> {

public:
  using summ_abs_domain_t = typename SumTable::abs_domain_t;
  using summary_t = typename SumTable::summary_t;
  using abs_dom_t = summ_abs_domain_t;
  using number_t = typename abs_dom_t::number_t;

private:
  using basic_block_t = typename SumTable::basic_block_t;
  using intra_abs_transform_t = intra_abs_transformer<basic_block_t, abs_dom_t>;
  typedef
      typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
  using varname_t = typename abs_dom_t::varname_t;
  using variable_t = typename abs_dom_t::variable_t;
  using linear_expression_t = typename abs_dom_t::linear_expression_t;

  SumTable *m_sum_tbl;

public:
  using typename abs_transform_api_t::callsite_t;

  static void reuse_summary(abs_dom_t &caller, const callsite_t &cs,
                            const summary_t &summ) {
    // Before a summary is plug-in we rename it with unique variable
    // names so we avoid naming clashes in cases like for instance
    // summary variables have same names as lhs of callsites.

    CRAB_LOG("inter", crab::outs() << "    Reuse summary at " << cs << "\n";
             crab::outs() << "    Summary:" << summ << "\n";);

    std::set<variable_t> actuals, formals;
    // --- matching formal and actual parameters
    auto inputs = summ.get_renamed_inputs();
    unsigned i = 0;
    // XXX: propagating down
    for (auto const &p : inputs) {
      auto const &a = cs.get_arg_name(i);
      if (!(a == p)) {
        CRAB_LOG("inter", crab::outs() << "\t\tPropagate from caller to callee "
                                       << p << ":=" << a << "\n");
        inter_transformer_helpers<abs_dom_t>::unify(caller, p, a);
      }
      ++i;
      actuals.insert(a);
      formals.insert(p);
    }

    // --- meet caller's inv with summ
    auto sum_inv = summ.get_renamed_sum();
    caller = caller & sum_inv;
    CRAB_LOG("inter", crab::outs() << "\t\tAfter meet caller and summary: "
                                   << caller << "\n");

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
    for (; caller_it != caller_et; ++caller_it, ++callee_it) {
      auto const &vt = *caller_it;
      auto const &r = *callee_it;
      CRAB_LOG("inter", crab::outs() << "\t\tPropagate from callee to caller "
                                     << vt << ":=" << r << "\n");
      inter_transformer_helpers<abs_dom_t>::unify(caller, vt, r);
      actuals.insert(vt);
      formals.insert(r);
    }

    // --- remove from caller only formal parameters so we can keep
    //     as much context from the caller as possible
    std::vector<variable_t> vs;
    std::set_difference(formals.begin(), formals.end(), actuals.begin(),
                        actuals.end(), std::back_inserter(vs));
    caller.forget(vs);

    CRAB_LOG("inter", crab::outs()
                          << "\t\tAfter forgetting formal parameters {";
             for (auto const &v
                  : vs) crab::outs()
             << v << ";";
             crab::outs() << "}=" << caller << "\n");
  }

public:
  bu_summ_abs_transformer(abs_dom_t inv, SumTable *sum_tbl)
      : intra_abs_transform_t(inv), m_sum_tbl(sum_tbl) {}

  ~bu_summ_abs_transformer() = default;

  virtual void exec(callsite_t &cs) override {
    if (!m_sum_tbl) {
      CRAB_WARN("The summary table is empty: ignored analysis of callsite");
    } else {
      if (m_sum_tbl->has_summary(cs)) {
        auto summ = m_sum_tbl->get(cs);
        reuse_summary(this->m_inv, cs, summ);
      } else {
        CRAB_WARN("summary not found for ", cs);
        for (auto vt : cs.get_lhs()) {
          this->m_inv -= vt; // havoc
        }
      }
    }
  }
};

//////////////////////////////////////////////////////////////////////
/// Generic conversion between top-down and bottom-up domains.
//////////////////////////////////////////////////////////////////////
/// FIXME: the conversion is by converting to linear constraints.
/// This will lose any array, disjunctive, or non-linear information
/// that "from" carries on.
//////////////////////////////////////////////////////////////////////
template <typename Domain1, typename Domain2>
inline void convert_domains_impl(Domain1 from, Domain2 &to) {
  if (from.is_bottom()) {
    to.set_to_bottom();
    return;
  }

  CRAB_LOG("inter", crab::outs()
                        << "Converting from " << from.domain_name() << " to "
                        << to.domain_name() << " might lose precision if "
                        << from.domain_name() << " is more precise than "
                        << to.domain_name() << "\n");

  bool pre_bot = false;
  if (::crab::CrabSanityCheckFlag) {
    pre_bot = from.is_bottom();
  }

  for (auto const &cst : from.to_linear_constraint_system()) {
    to += cst;
  }

  if (::crab::CrabSanityCheckFlag) {
    bool post_bot = to.is_bottom();
    if (!(pre_bot || !post_bot)) {
      CRAB_ERROR("Invariant became bottom after conversion between domains");
    }
  }
}

template <typename Variable>
inline void convert_domains(crab::domains::abstract_domain<Variable> from,
                            crab::domains::abstract_domain<Variable> &to) {
  // We don't know what is inside "from" or "to" so we do
  // conservatively the conversion. But note that we might lose
  // unnecessarily precision if from and to are actually the same
  // domain.
  convert_domains_impl(from, to);
}

template <typename Variable>
inline void convert_domains(crab::domains::abstract_domain_ref<Variable> from,
                            crab::domains::abstract_domain_ref<Variable> &to) {
  // We don't know what is inside "from" or "to" so we do
  // conservatively the conversion. But note that we might lose
  // unnecessarily precision if from and to are actually the same
  // domain.
  convert_domains_impl(from, to);
}
  
template <typename Domain>
inline void convert_domains(Domain from, Domain &to) {
  // do nothing if they are the same domain but not
  // abstract_domain.
  to = from;
}

template <typename Domain1, typename Domain2>
inline void convert_domains(Domain1 from, Domain2 &to) {
  convert_domains_impl(from, to);
}

/**
 * Abstract transformer specialized for performing top-down forward
 * traversal while reusing summaries at the callsites.
 * class CallCtxTable stores the calling context.
 **/
template <class SumTable, class CallCtxTable>
class td_summ_abs_transformer final
    : public intra_abs_transformer<typename CallCtxTable::basic_block_t,
                                   typename CallCtxTable::abs_domain_t> {

public:
  using summ_abs_domain_t = typename SumTable::abs_domain_t;
  using call_abs_domain_t = typename CallCtxTable::abs_domain_t;
  using basic_block_t = typename CallCtxTable::basic_block_t;
  using abs_dom_t = call_abs_domain_t;
  using number_t = typename abs_dom_t::number_t;

private:
  static_assert(std::is_same<typename CallCtxTable::basic_block_t,
                             typename SumTable::basic_block_t>::value,
                "Summary table and calling context table must have same basic "
                "block type");

  static_assert(
      std::is_same<typename summ_abs_domain_t::number_t,
                   typename call_abs_domain_t::number_t>::value,
      "Bottom-up and top-down abstract domains must have same number type");

  static_assert(std::is_same<typename summ_abs_domain_t::varname_t,
                             typename call_abs_domain_t::varname_t>::value,
                "Bottom-up and top-down abstract domains must have same "
                "variable name type");

  using intra_abs_transform_t = intra_abs_transformer<basic_block_t, abs_dom_t>;

public:
  typedef
      typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
  using typename abs_transform_api_t::callsite_t;

private:
  SumTable *m_sum_tbl;
  CallCtxTable *m_call_tbl;

public:
  td_summ_abs_transformer(abs_dom_t inv, SumTable *sum_tbl,
                          CallCtxTable *call_tbl)
      : intra_abs_transform_t(inv), m_sum_tbl(sum_tbl), m_call_tbl(call_tbl) {}

  ~td_summ_abs_transformer() = default;

  virtual void exec(callsite_t &cs) override {
    if (!m_sum_tbl) {
      CRAB_WARN("The summary table is empty");
    } else if (m_sum_tbl->has_summary(cs)) {
      auto summ = m_sum_tbl->get(cs);

      CRAB_LOG("inter", crab::outs()
                            << "    Pluging caller context into callee\n"
                            << "    Summary: " << summ << "\n");

      ///////
      /// Generate the callee context and store it.
      ///////
      using variable_t = typename abs_dom_t::variable_t;

      abs_dom_t callee_ctx_inv(this->m_inv);
      // --- matching formal and actual parameters
      // XXX: propagating down
      unsigned i = 0;
      const std::vector<variable_t> &inputs = summ.get_inputs();
      for (const variable_t &p : inputs) {
        const variable_t &a = cs.get_arg_name(i);
        if (!(a == p)) {
          inter_transformer_helpers<abs_dom_t>::unify(callee_ctx_inv, p, a);
        }
        ++i;
      }

      // --- project only onto formal parameters
      callee_ctx_inv.project(inputs);
      // --- store the callee context
      CRAB_LOG("inter", crab::outs() << "\t\tCallee context stored: "
                                     << callee_ctx_inv << "\n");

      m_call_tbl->insert(cs, callee_ctx_inv);

      /////
      // Generate the continuation at the caller
      /////

      // --- convert this->m_inv to the language of summ_abs_dom_t(summ)
      summ_abs_domain_t caller_ctx_inv = summ.make_top();
      convert_domains(this->m_inv, caller_ctx_inv);
      // forget the variables of the lhs of the callsite, otherwise
      // caller_ctx_inv and m_inv may be inconsistent if the lhs
      // variables are constrained before the callsite (e.g., x=-5;
      // x := abs(x);)
      for (auto const &vt : cs.get_lhs()) {
        this->m_inv -= vt;
      }

      CRAB_LOG("inter", crab::outs() << "\t\tCaller context: " << caller_ctx_inv
                                     << "\n");

      // --- reuse summary to do the continuation
      bu_summ_abs_transformer<SumTable>::reuse_summary(caller_ctx_inv, cs,
                                                       summ);
      CRAB_LOG("inter", crab::outs()
                            << "\t\tCaller context after plugin summary: "
                            << caller_ctx_inv << "\n");

      // --- convert back inv to the language of abs_dom_t
      convert_domains(caller_ctx_inv, this->m_inv);
      CRAB_LOG("inter", crab::outs() << "\t\tCaller continuation after " << cs
                                     << "=" << this->m_inv << "\n");
      return;
    } else {
      // --- no summary found: havoc output variables
      CRAB_WARN("No summary found during top-down pass for ", cs);
    }

    // We could not reuse a summary or summary not found. We havoc
    // output variables of the callsite.
    for (auto const &vt : cs.get_lhs()) {
      this->m_inv -= vt;
    }
  }
};
} // namespace inter_analyzer_impl

template <typename CallGraph,
          // abstract domain used for the bottom-up phase
          typename BU_Dom,
          // abstract domain used for the top-down phase
          typename TD_Dom>
class bottom_up_inter_analyzer:
    public inter_analyzer_api<CallGraph, TD_Dom, BU_Dom> {

  using cg_node_t = typename CallGraph::node_t;
  using cg_edge_t = typename CallGraph::edge_t;

  using this_type = bottom_up_inter_analyzer<CallGraph, BU_Dom, TD_Dom>;

public:
  using cg_t = CallGraph;
  using cfg_t = typename cg_node_t::cfg_t;
  using varname_t = typename cfg_t::varname_t;
  using number_t = typename cfg_t::number_t;
  using variable_t = typename cfg_t::variable_t;
  using liveness_t = live_and_dead_analysis<cfg_t>;
  using liveness_map_t = std::unordered_map<cfg_t, const liveness_t *>;
  using params_t = inter_analyzer_parameters<CallGraph>;
  
private:
  using summ_tbl_t = inter_analyzer_impl::summary_table<cfg_t, BU_Dom>;
  using call_tbl_t = inter_analyzer_impl::call_ctx_table<cfg_t, TD_Dom>;
  using cg_ref_t = crab::cg::call_graph_ref<cg_t>;
  using callsite_or_fdecl_t = crab::cfg::callsite_or_fdecl<cfg_t>;  
public:
  using bu_abs_tr = inter_analyzer_impl::bu_summ_abs_transformer<summ_tbl_t>;
  using td_abs_tr =
      inter_analyzer_impl::td_summ_abs_transformer<summ_tbl_t, call_tbl_t>;
  using bu_analyzer = analyzer_internal_impl::fwd_analyzer<cfg_t, bu_abs_tr>;
  using td_analyzer = analyzer_internal_impl::fwd_analyzer<cfg_t, td_abs_tr>;
  using summary_t = typename inter_analyzer_api<CallGraph,TD_Dom,BU_Dom>::summary_t;
  
  // for checkers
  using abs_dom_t = TD_Dom;
  using abs_tr_t = td_abs_tr;

private:
  using td_analyzer_ptr = std::unique_ptr<td_analyzer>;
  using invariant_map_t = crab::cfg::callsite_or_fdecl_map<cfg_t, td_analyzer_ptr>;

  CallGraph &m_cg;
  // These two used to call make_top() and make_bottom().
  TD_Dom m_td_absval_fac;
  BU_Dom m_bu_absval_fac;
  const liveness_map_t *m_live;
  invariant_map_t m_inv_map;
  summ_tbl_t m_summ_tbl;
  call_tbl_t m_call_tbl;
  fixpoint_parameters m_fixpo_params;
  std::unique_ptr<abs_tr_t> m_abs_tr;

  const liveness_t *get_live(const cfg_t &c) {
    if (m_live) {
      auto it = m_live->find(c);
      if (it != m_live->end())
        return it->second;
    }
    return nullptr;
  }

  inline TD_Dom make_td_bottom() const { return m_td_absval_fac.make_bottom(); }
  inline TD_Dom make_td_top() const { return m_td_absval_fac.make_top(); }
  inline BU_Dom make_bu_bottom() const { return m_bu_absval_fac.make_bottom(); }
  inline BU_Dom make_bu_top() const { return m_bu_absval_fac.make_top(); }

public:
  bottom_up_inter_analyzer(CallGraph &cg,
                           TD_Dom td_absval_fac, BU_Dom bu_absval_fac,
			   const params_t &params = params_t())
    : m_cg(cg), m_td_absval_fac(td_absval_fac), m_bu_absval_fac(bu_absval_fac),
      m_live(params.live_map),
      m_call_tbl(make_td_top()),
      m_abs_tr(new abs_tr_t(make_td_top(), &m_summ_tbl, &m_call_tbl)) {
    
    m_fixpo_params.get_widening_delay() = params.widening_delay;
    m_fixpo_params.get_descending_iterations() = params.descending_iters;
    m_fixpo_params.get_max_thresholds() = params.thresholds_size;
      
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking call graph ... ";);
    crab::CrabStats::resume("CallGraph type checking");
    cg.type_check();
    crab::CrabStats::stop("CallGraph type checking");
    CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);
  }

  bottom_up_inter_analyzer(const this_type &other) = delete;

  this_type &operator=(const this_type &other) = delete;

  /** =====  Run the inter-procedural analysis
   *
   * init is the initial abstract state used for the top-down analysis
   * The bottom-up analysis will start with top (i.e., no
   * preconditions).
   *
   * - During the bottom-up analysis, it computes a summary for each
   *   function except for main.
   * - The top-down analysis reuses the callsite's summary each time
   *   it visits the callsite.
   *
   *   FIXME: the top-down analysis runs only once starting from the
   *   _first_ SCC encountered after topologically ordering all the callgraph
   *   SCCs. It's also assuming that first SCC has only one
   *   component. In general, the callgraph might have more than one
   *   entry point and an entry SCC can have multiple components.
   **/
  void run(TD_Dom init) override {

    CRAB_VERBOSE_IF(1, get_msg_stream()
                           << "Started inter-procedural analysis\n";);
    CRAB_LOG("inter", m_cg.write(crab::outs()); crab::outs() << "\n");

    crab::ScopedCrabStats __st__("Inter");

    bool has_noedges = true;
    for (auto const &n : boost::make_iterator_range(m_cg.nodes())) {
      if (m_cg.num_succs(n) > 0) {
        has_noedges = false;
        break;
      }
    }

    if (has_noedges) {
      // -- Special case when the program has been inlined.
      CRAB_LOG("inter",
               crab::outs()
                   << "Callgraph has no edges so no summaries are computed.\n");

      CRAB_LOG("inter", m_cg.write(crab::outs()); crab::outs() << "\n";);

      for (auto &v : boost::make_iterator_range(m_cg.nodes())) {
        crab::ScopedCrabStats __st__("Inter.TopDown");

        auto cfg = v.get_cfg();
        assert(cfg.has_func_decl());
        auto const &fdecl = cfg.get_func_decl();
        const std::string &fun_name = fdecl.get_func_name();
        if (fun_name != "main")
          continue;

        CRAB_LOG("inter", crab::outs() << "++ Analyzing function "
                                       << fdecl.get_func_name() << "\n");

        td_analyzer_ptr a(new td_analyzer(cfg, &*m_abs_tr, m_td_absval_fac, get_live(cfg), m_fixpo_params)); 
        a->run_forward(init);
        m_inv_map.insert(std::make_pair(callsite_or_fdecl_t(&fdecl), std::move(a)));
      }
      return;
    }

    // -- General case
    std::vector<cg_node_t> rev_order;
    graph_algo::scc_graph<cg_ref_t> Scc_g(m_cg);
    graph_algo::rev_topo_sort(Scc_g, rev_order);

    CRAB_VERBOSE_IF(1, get_msg_stream() << "== Bottom-up phase ...\n";);
    for (auto const &n : rev_order) {
      crab::ScopedCrabStats __st__("Inter.BottomUp");
      std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members(n);
      for (auto m : scc_mems) {

        auto cfg = m.get_cfg();
        assert(cfg.has_func_decl());
        auto const &fdecl = cfg.get_func_decl();
        const std::string &fun_name = fdecl.get_func_name();
        if (fun_name != "main") {
          CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Computing summary for "
                                              << fun_name << "...\n";);

          // --- collect inputs and outputs
          std::vector<variable_t> formals, inputs, outputs;
          formals.reserve(fdecl.get_num_inputs() + fdecl.get_num_outputs());
          inputs.reserve(fdecl.get_num_inputs());
          outputs.reserve(fdecl.get_num_outputs());

          for (unsigned i = 0; i < fdecl.get_num_inputs(); i++) {
            inputs.push_back(fdecl.get_input_name(i));
            formals.push_back(fdecl.get_input_name(i));
          }
          for (unsigned i = 0; i < fdecl.get_num_outputs(); i++) {
            outputs.push_back(fdecl.get_output_name(i));
            formals.push_back(fdecl.get_output_name(i));
          }

          if (outputs.empty()) {
            BU_Dom summary = make_bu_top();
            m_summ_tbl.insert(fdecl, summary, inputs, outputs);
            CRAB_WARN("Skipped summary because function ", fun_name,
                      " has no output parameters");
          } else if (!cfg.has_exit()) {
            CRAB_WARN("Skipped summary because function ", fun_name,
                      " has no exit block");
          } else {
            // --- run the analysis
            bu_abs_tr abs_tr(std::move(make_bu_top()), &m_summ_tbl);
            bu_analyzer a(cfg, &abs_tr, m_bu_absval_fac, get_live(cfg), m_fixpo_params);
            a.run_forward(make_bu_top());

            // --- project onto formal parameters and return values
            BU_Dom summary = a.get_post(cfg.exit());
            // crab::CrabStats::count(BU_Dom::getDomainName() +
            // ".count.project");
            summary.project(formals);
            m_summ_tbl.insert(fdecl, summary, inputs, outputs);
          }
        }
      }
    }

    CRAB_VERBOSE_IF(1, get_msg_stream() << "== Top-down phase ...\n";);
    bool is_root = true;
    for (auto n :
         boost::make_iterator_range(rev_order.rbegin(), rev_order.rend())) {
      crab::ScopedCrabStats __st__("Inter.TopDown");
      std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members(n);

      // The SCC is recursive if it has more than one element or
      // there is only one that calls directly to itself.
      bool is_recursive =
          (scc_mems.size() > 1) ||
          std::any_of(m_cg.succs(n).first, m_cg.succs(n).second,
                      [n](const cg_edge_t &e) { return (n == e.dest()); });

      for (auto m : scc_mems) {
        auto cfg = m.get_cfg();
        assert(cfg.has_func_decl());
        auto const &fdecl = cfg.get_func_decl();
        CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Analyzing function "
                                            << fdecl.get_func_name() << "\n";);
        if (is_recursive) {
          // If the SCC is recursive then what we have in
          // m_call_tbl is incomplete and therefore it is unsound
          // to use it. To remedy it, we insert another calling
          // context with top value that approximates all the
          // possible calling contexts during the recursive calls.
          m_call_tbl.insert(fdecl, make_td_top());
        }

        auto init_inv = init;
        if (is_root) {
          is_root = false;
        } else {
          init_inv = m_call_tbl.get_call_ctx(fdecl);
        }

        CRAB_LOG("inter", crab::outs() << "    Starting analysis of " << fdecl
                                       << " with " << init_inv << "\n");

        if (init_inv.is_bottom()) {
          crab::outs() << "Top-down analysis for " << fdecl.get_func_name()
                       << " started with bottom (i.e., dead function).\n";
        }
        td_analyzer_ptr a(new td_analyzer(cfg, &*m_abs_tr, m_td_absval_fac, get_live(cfg), 
                                          m_fixpo_params));
        a->run_forward(init_inv);
        m_inv_map.insert(std::make_pair(callsite_or_fdecl_t(&fdecl), std::move(a))); 
      }
    }
    CRAB_VERBOSE_IF(1, get_msg_stream()
                           << "Finished inter-procedural analysis\n";);
  }

  //! return the analyzed call graph
  CallGraph &get_call_graph() override { return m_cg; }

  //! Return the invariants that hold at the entry of b in cfg
  TD_Dom get_pre(const cfg_t &cfg,
                 const typename cfg_t::basic_block_label_t &b) const override {
    assert(cfg.has_func_decl());
    auto const &fdecl = cfg.get_func_decl();
    auto const it = m_inv_map.find(&fdecl);
    if (it != m_inv_map.end()) {
      return it->second->get_pre(b);
    } else {
      return make_td_top();
    }
  }

  //! Return the invariants that hold at the exit of b in cfg
  TD_Dom get_post(const cfg_t &cfg,
                  const typename cfg_t::basic_block_label_t &b) const override {
    assert(cfg.has_func_decl());
    auto const &fdecl = cfg.get_func_decl();
    auto const it = m_inv_map.find(&fdecl);
    if (it != m_inv_map.end()) {
      return it->second->get_post(b);
    } else {
      return make_td_top();
    }
  }

  // clear all the analysis' state
  void clear() override {
    m_inv_map.clear();
    m_summ_tbl.clear();
    m_call_tbl.clear();
    m_abs_tr->get_abs_value().set_to_top();
  }

  summary_t get_summary(const cfg_t &cfg) const override {
    // TODO: caching

    // We convert from our internal representation of a summary to the
    // summary representation exposed by inter_analysis_api
    summary_t summary(cfg.get_func_decl());
    auto const& fdecl = cfg.get_func_decl();
    if (m_summ_tbl.has_summary(fdecl)) {
      auto _summary =  m_summ_tbl.get(fdecl);
      // The summary is context-insensitive. This means that the
      // precondition is top.
      summary.add(_summary.make_top(), _summary.get_sum());
    }
    return summary;
  }
  
  // /*DEPRECATED*/ Propagate inv through statements: 
  abs_tr_t &get_abs_transformer() {
    assert(m_abs_tr);
    return *m_abs_tr;
  } 
};

} // namespace analyzer
} // namespace crab
