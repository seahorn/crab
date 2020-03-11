#pragma once

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/graphs/dominance.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/iterators/interleaved_fixpoint_iterator.hpp>

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

namespace crab {

namespace analyzer {

/**
 * Compute necessary preconditions by computing the abstract least
 * fixpoint of the reversed CFG. Optionally, at each step the
 * precondition can be refined with an invariant computed by a
 * forward analysis.
 **/
template <typename CFG, typename AbsDom>
class necessary_preconditions_fixpoint_iterator
    : private ikos::interleaved_fwd_fixpoint_iterator<crab::cfg::cfg_rev<CFG>,
                                                      AbsDom> {

  typedef ikos::interleaved_fwd_fixpoint_iterator<crab::cfg::cfg_rev<CFG>,
                                                  AbsDom>
      fixpoint_iterator_t;
  typedef typename CFG::basic_block_label_t bb_label_t;
  typedef typename CFG::statement_t stmt_t;
  typedef std::unordered_map<bb_label_t, AbsDom> bb_abstract_map_t;
  typedef std::unordered_map<stmt_t *, AbsDom> pp_abstract_map_t;

  /// forward and backward transformers
  typedef intra_abs_transformer<AbsDom> abs_fwd_tr_t;
  typedef intra_necessary_preconditions_abs_transformer<AbsDom,
                                                        pp_abstract_map_t>
      abs_bwd_tr_t;

  CFG m_cfg;
  // postcondition (i.e, final states) that we want to propagate backwards
  AbsDom m_postcond;
  // necessary preconditions
  bb_abstract_map_t m_preconditions;
  // forward invariants
  bb_abstract_map_t m_invariants;
  // preconditions from good states, otherwise from bad states
  bool m_good_states;

  /**
   * Compute necessary preconditions for a basic block
   **/
  virtual AbsDom analyze(bb_label_t node, AbsDom &&precond) override {
    auto &bb = m_cfg.get_node(node);

    CRAB_LOG("backward-fixpoint", crab::outs() << "Post at "
                                               << cfg_impl::get_label_str(node)
                                               << ": " << precond << "\n");

    // invariants that hold at the entry of the block
    AbsDom invariant = m_invariants[node];
    // rebuild local invariants that hold at each program point.
    abs_fwd_tr_t F(invariant);
    pp_abstract_map_t pp_invariants;
    for (auto &s : boost::make_iterator_range(bb.begin(), bb.end())) {
      auto inv = F.get_abs_value();
      pp_invariants.insert(std::make_pair(&s, inv));
      CRAB_LOG("backward-fixpoint", crab::outs()
                                        << "\tRebuilding at statement " << s
                                        << " inv=" << inv << "\n");
      s.accept(&F);
    }

    CRAB_LOG("backward-fixpoint",
             crab::outs() << "Done forward propagation at each program point \n"
                          << "Starting backward propagation ... \n");

    // compute precondition at the entry of the block
    abs_bwd_tr_t B(std::move(precond), &pp_invariants, m_good_states);
    for (auto &s : boost::make_iterator_range(bb.rbegin(), bb.rend())) {
      s.accept(&B);
    }
    precond = std::move(B.preconditions());
    CRAB_LOG("backward-fixpoint", crab::outs() << "Pre at "
                                               << cfg_impl::get_label_str(node)
                                               << ": " << precond << "\n");
    return std::move(precond);
  }

  virtual void process_pre(bb_label_t /*node*/, AbsDom /*postcond*/) override {}

  /**
   *  Store necessary preconditions
   **/
  virtual void process_post(bb_label_t node, AbsDom precond) override {
    m_preconditions.insert(std::make_pair(node, precond));
  }

public:
  typedef typename fixpoint_iterator_t::wto_t wto_t;
  typedef bb_abstract_map_t precond_map_t;
  typedef typename precond_map_t::iterator iterator;
  typedef typename precond_map_t::const_iterator const_iterator;

  // This constructor computes necessary preconditions from error
  // states.
  necessary_preconditions_fixpoint_iterator(
      CFG cfg, AbsDom postcond,
      /* fixpoint parameters */
      unsigned int widening_delay = 1,
      unsigned int descending_iterations = UINT_MAX, size_t jump_set_size = 0)
      : fixpoint_iterator_t(crab::cfg::cfg_rev<CFG>(cfg), nullptr,
                            widening_delay, descending_iterations,
                            jump_set_size),
        m_cfg(cfg), m_postcond(postcond), m_good_states(false) {}

  // This constructor computes necessary preconditions from
  // safe/good (error) states if good_states is true (false).
  necessary_preconditions_fixpoint_iterator(
      CFG cfg, const wto_t *wto, AbsDom postcond, bool good_states,
      /* fixpoint parameters */
      unsigned int widening_delay = 1,
      unsigned int descending_iterations = UINT_MAX, size_t jump_set_size = 0)
      : fixpoint_iterator_t(crab::cfg::cfg_rev<CFG>(cfg), wto, widening_delay,
                            descending_iterations, jump_set_size),
        m_cfg(cfg), m_postcond(postcond), m_good_states(good_states) {}

  void run_backward() { this->run(m_postcond); }

  void run_backward(const std::unordered_map<typename CFG::basic_block_label_t,
                                             AbsDom> &fwd_invariants) {
    m_invariants.insert(fwd_invariants.begin(), fwd_invariants.end());
    run_backward();
  }

  iterator begin() { return m_preconditions.begin(); }

  iterator end() { return m_preconditions.end(); }

  const_iterator begin() const { return m_preconditions.begin(); }

  const_iterator end() const { return m_preconditions.end(); }

  // return the preconditions at basic block node
  AbsDom operator[](bb_label_t node) const {
    auto it = m_preconditions.find(node);
    if (it != m_preconditions.end())
      return it->second;
    else
      return AbsDom::top();
  }

  // clear preconditions and forward invariants (if any)
  void clear() {
    m_preconditions.clear();
    m_invariants.clear();
  }

  const wto_t &get_WTO() const { return this->get_wto(); }
};

/**
 * A refining forward-backward analyzer based on Cousot&Cousot
 * (JLP'92 and ASE'99, section 4).
 *
 * Do first a forward pass computing superset of reachable states
 * from initial states. Then do a backward pass, computing
 * superset of reachable meet with co-reachable from error
 * states. Repeat this until no more refinement or maximum of
 * iterations reached.
 *
 * The implementation follows mc_5 algorithm described in page 82,
 * ASE'99. It computes an over-approximation of the intersection
 * between the set of reachable states starting from the initial
 * states and the set of co-reachable states from error
 * states.
 *
 * This forward-backward refinement implementation produces the
 * same invariants than a classical forward analyzer. The
 * difference is that it can prove more assertions.
 **/
template <typename CFG, typename AbsDom> class intra_forward_backward_analyzer {
public:
  typedef CFG cfg_t;
  typedef typename CFG::basic_block_t bb_t;
  typedef typename CFG::basic_block_label_t bb_label_t;
  typedef typename CFG::statement_t stmt_t;
  typedef typename CFG::varname_t varname_t;
  typedef typename CFG::number_t number_t;
  // used for checkers
  typedef AbsDom abs_dom_t;
  typedef intra_abs_transformer<AbsDom> abs_tr_t;

private:
  typedef intra_fwd_analyzer<CFG, AbsDom> fwd_analyzer_t;
  typedef necessary_preconditions_fixpoint_iterator<CFG, AbsDom>
      bwd_fixpoint_iterator_t;
  typedef typename bwd_fixpoint_iterator_t::precond_map_t precond_map_t;
  typedef typename bwd_fixpoint_iterator_t::wto_t bwd_wto_t;
  typedef liveness<CFG> liveness_t;
  typedef std::unordered_map<bb_label_t, std::set<bb_label_t>> idom_tree_t;
  typedef typename bb_t::assert_t assert_t;
  typedef typename bb_t::bool_assert_t bool_assert_t;
  typedef typename bb_t::ptr_assert_t ptr_assert_t;

public:
  typedef typename fwd_analyzer_t::assumption_map_t assumption_map_t;
  typedef typename fwd_analyzer_t::invariant_map_t invariant_map_t;
  // bwd_wto_t and wto_t are different types because bwd_wto_t is
  // over the reversed CFG.
  typedef typename fwd_analyzer_t::wto_t wto_t;

private:
  // -- the cfg
  CFG m_cfg;
  // we keep the two wto's (from forward and reversed CFGs) to
  // avoid recompute them during the below iterative process.
  // Only the forward wto is exposed to outside clients.
  const wto_t *m_wto;
  const bwd_wto_t *m_b_wto;
  // keep track of which assertions cannot be proven.
  std::vector<std::pair<bb_label_t, stmt_t *>> m_unproven_assertions;
  // keep track of which assertions have been proven
  std::set<stmt_t *> m_proved_assertions;
  // keep the results of the first forward iteration.
  invariant_map_t m_pre_invariants;
  invariant_map_t m_post_invariants;

  void store_analysis_results(fwd_analyzer_t &f) {
    for (auto &kv : f.get_pre_invariants()) {
      m_pre_invariants.insert({kv.first, std::move(kv.second)});
    }

    for (auto &kv : f.get_post_invariants()) {
      m_post_invariants.insert({kv.first, std::move(kv.second)});
    }
  }

  void gather_assertions() {
    for (auto it = m_cfg.begin(), et = m_cfg.end(); it != et; ++it) {
      for (auto &s : *it) {
        if (s.is_assert() || s.is_bool_assert() || s.is_ptr_assert()) {
          m_unproven_assertions.push_back({it->label(), &s});
        }
      }
    }
  }

  // idom induces a tree each key-value pair means that key
  // (block) is an immediate dominator of each element in value
  // (set of basic blocks). TODO: caching
  bool dominates(bb_label_t u, bb_label_t v, const idom_tree_t &idom) {
    auto it = idom.find(u);
    if (it == idom.end()) {
      // u does not dominate any block
      return false;
    }

    // immediate dominated by u
    const std::set<bb_label_t> &idom_u = it->second;
    if (idom_u.count(v) > 0) {
      // u dominates v
      return true;
    }

    for (const bb_label_t w : idom_u) {
      // a successor of u dominates v
      if (dominates(w, v, idom)) {
        return true;
      }
    }

    return false;
  }

  /**
   * For each pair (bb, dom) in refined_assumptions we know that
   * dom is the meet of the over-approximated set of reachable
   * states from bb and the over-approximated co-reachable states
   * from the error states. Therefore, if dom is bottom then any
   * assertion dominated (in the sense of CFG dominance) by bb
   * must hold and hence, it is removed from
   * m_unproven_assertions and added to m_proved_assertions.
   **/
  void discharge_assertions(assumption_map_t &refined_assumptions,
                            const idom_tree_t &idom) {
    if (m_unproven_assertions.empty()) {
      return;
    }

    if (!idom.empty()) {
      // If it's not possible starting from the initial states to
      // reach a block B that leads to an assertion violation then
      // we can discharge all the assertions dominated by B.
      //
      // Once we know that an assertion must hold, we could add in
      // m_pre_invariants an extra constraint from the assertion
      // condition. The problem is that the condition might not be
      // expressed in the abstract domain (e.g., a disequality). So
      // that wouldn't be enough for the checker to discharge the
      // assertion, even if we know already the assertion must
      // hold. Instead, we tell explicitly the checker which
      // assertions were proved by calling safe_assertions method.

      for (auto it = m_cfg.begin(), et = m_cfg.end(); it != et; ++it) {
        bb_label_t n = it->label();
        auto jt = refined_assumptions.find(n);
        if (jt == refined_assumptions.end()) {
          continue;
        }
        if (jt->second.is_bottom()) {
          m_unproven_assertions.erase(
              std::remove_if(
                  m_unproven_assertions.begin(), m_unproven_assertions.end(),
                  [&n, &idom, this](std::pair<bb_label_t, stmt_t *> &kv) {
                    bb_label_t m = kv.first;
                    stmt_t *s = kv.second;
                    bool res = this->dominates(n, m, idom);
                    if (res) {
                      this->m_proved_assertions.insert(s);
                    }
                    return res;
                  }),
              m_unproven_assertions.end());
        }
      }
    } else {
      // if dominance information is not available then we only
      // discharge all assertions only if starting from the entry
      // block we know that we cannot violate them.
      auto it = refined_assumptions.find(m_cfg.entry());
      if (it->second.is_bottom()) {
        for (auto &kv : m_unproven_assertions) {
          m_proved_assertions.insert(kv.second);
        }
        m_unproven_assertions.clear();
      }
    }
  }

public:
  typedef typename bwd_fixpoint_iterator_t::iterator iterator;
  typedef typename bwd_fixpoint_iterator_t::const_iterator const_iterator;

  intra_forward_backward_analyzer(CFG cfg)
      : m_cfg(cfg), m_wto(nullptr), m_b_wto(nullptr) {}

  ~intra_forward_backward_analyzer() {
    if (m_wto)
      delete m_wto;
    if (m_b_wto)
      delete m_b_wto;
  }

  intra_forward_backward_analyzer(
      const intra_forward_backward_analyzer<CFG, AbsDom> &o) = delete;
  intra_forward_backward_analyzer<CFG, AbsDom> &
  operator=(const intra_forward_backward_analyzer<CFG, AbsDom> &o) = delete;

  /**
   * Perform the refining forward-backward loop.
   **/
  void run(AbsDom init_states,
           // behaves as a standard forward analysis
           bool only_forward,
           // assumptions
           const assumption_map_t &assumptions,
           // liveness information
           const liveness_t *live,
           // parameters for each forward or backward analysis
           unsigned int widening_delay = 1,
           unsigned int descending_iters = UINT_MAX, size_t jump_set_size = 0) {
    run(m_cfg.entry(), init_states, only_forward, assumptions, live,
        widening_delay, descending_iters, jump_set_size);
  }

  void run(bb_label_t entry, // only used for the forward pass.
           AbsDom init_states,
           // behaves as a standard forward analysis
           bool only_forward,
           // assumptions
           const assumption_map_t &assumptions,
           // liveness information
           const liveness_t *live,
           // parameters for each forward or backward analysis
           unsigned int widening_delay = 1,
           unsigned int descending_iters = UINT_MAX, size_t jump_set_size = 0) {

    CRAB_LOG("backward", crab::outs()
                             << "Initial states=" << init_states << "\n");

    if (!only_forward && !m_cfg.has_exit()) {
      CRAB_WARN("cannot run backward analysis because CFG has no exit block");
      only_forward = true;
    }

    crab::CrabStats::count("CombinedForwardBackward.invocations");

    // maximum number of refinement iterations
    const unsigned max_num_iters = 5;
    // number of refinement iterations
    unsigned iters = 0;
    // CFG assertions
    std::vector<std::pair<bb_label_t, stmt_t *>> assertions;
    // immediate dominance tree
    idom_tree_t idom_tree;

    crab::CrabStats::resume("CombinedForwardBackward.GatherAssertions");
    gather_assertions();
    CRAB_LOG("backward", crab::outs()
                             << "Found " << m_unproven_assertions.size()
                             << " assertions.\n";);
    crab::CrabStats::stop("CombinedForwardBackward.GatherAssertions");

    if (!m_unproven_assertions.empty() && !only_forward) {
      crab::CrabStats::resume("CombinedForwardBackward.DominatorTree");
      std::unordered_map<bb_label_t, bb_label_t> idom_map;
      crab::analyzer::graph_algo::dominator_tree(m_cfg, m_cfg.entry(),
                                                 idom_map);
      // build idom_tree
      for (auto &kv : idom_map) {
        if (kv.second == boost::graph_traits<CFG>::null_vertex()) {
          continue;
        }
        auto it = idom_tree.find(kv.second);
        if (it == idom_tree.end()) {
          std::set<bb_label_t> reachable({kv.first});
          idom_tree.insert({kv.second, reachable});
        } else {
          it->second.insert(kv.first);
        }
      }

      CRAB_LOG("backward", crab::outs() << "Computed dominance tree:\n";
               for (auto &kv
                    : idom_tree) {
                 crab::outs() << "\t" << cfg_impl::get_label_str(kv.first)
                              << " dominates={";
                 for (auto d : kv.second) {
                   crab::outs() << d << ";";
                 }
                 crab::outs() << "}\n";
               });
      crab::CrabStats::stop("CombinedForwardBackward.DominatorTree");
    }

    assumption_map_t refined_assumptions(assumptions.begin(),
                                         assumptions.end());
    while (true) {
      iters++;
      crab::CrabStats::count("CombinedForwardBackward.iterations");
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Iteration " << iters << "\n"
                                          << "Started forward analysis.\n";);

      crab::CrabStats::resume("CombinedForwardBackward.ForwardPass");
      // run forward analysis refined with preconditions from error states
      fwd_analyzer_t F(m_cfg, init_states, live, m_wto, widening_delay,
                       descending_iters, jump_set_size);
      F.run(entry, refined_assumptions);

      if (iters == 1) {
        store_analysis_results(F);
      }
      crab::CrabStats::stop("CombinedForwardBackward.ForwardPass");

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished forward analysis.\n";);

      CRAB_LOG(
          "backward", crab::outs() << "Forward analysis: \n";
          for (auto &kv
               : boost::make_iterator_range(F.pre_begin(), F.pre_end())) {
            crab::outs() << cfg_impl::get_label_str(kv.first) << ":\n"
                         << kv.second << "\n";
          } crab::outs()
          << "\n";);

      if (only_forward) {
        CRAB_VERBOSE_IF(1, get_msg_stream()
                               << "\nSkipped backward refinement.\n";);
        break;
      } else if (m_unproven_assertions.empty()) {
        CRAB_VERBOSE_IF(
            1, get_msg_stream()
                   << "\nNo assertions found: skipped backward refinement.\n";);
        break;
      }

      // reuse wto for next iteration
      if (iters == 1) {
        m_wto = new wto_t(F.get_wto());
      }

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Started backward analysis.\n";);

      crab::CrabStats::resume("CombinedForwardBackward.BackwardPass");
      // run backward analysis computing necessary preconditions
      // refined with results from the forward analysis.
      AbsDom final_states = AbsDom::bottom();
      bwd_fixpoint_iterator_t B(m_cfg, m_b_wto,
                                // A final state is safe so here means bottom
                                final_states,
                                // negate assertions:
                                // preconditions from error states
                                false, widening_delay, descending_iters,
                                jump_set_size);

      const invariant_map_t &fwd_invariants = F.get_pre_invariants();
      crab::CrabStats::resume("CombinedForwardBackward.MinimizeInvariants");
      // Important for apron and elina domains.
      invariant_map_t minimized_fwd_invariants;
      for (auto &kv : fwd_invariants) {
        AbsDom dom(kv.second);
        dom.minimize();
        minimized_fwd_invariants.insert({kv.first, dom});
      }
      crab::CrabStats::stop("CombinedForwardBackward.MinimizeInvariants");
      B.run_backward(minimized_fwd_invariants);
      crab::CrabStats::stop("CombinedForwardBackward.BackwardPass");

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished backward analysis.\n";);

      CRAB_LOG(
          "backward", crab::outs() << "Backward analysis:\n";
          for (auto &kv
               : boost::make_iterator_range(B.begin(), B.end())) {
            crab::outs() << cfg_impl::get_label_str(kv.first) << ":\n"
                         << kv.second << "\n";
          } crab::outs()
          << "\n");

      crab::CrabStats::resume("CombinedForwardBackward.CheckRefinement");
      assumption_map_t new_refined_assumptions;
      bool more_refinement = false;
      for (auto it = m_cfg.begin(), et = m_cfg.end(); it != et; ++it) {
        AbsDom x = refined_assumptions[it->label()];
        AbsDom y = B[it->label()];
        AbsDom x_narrowing_y = x && y;
        more_refinement |= (!(x <= x_narrowing_y));
        new_refined_assumptions.insert({it->label(), x_narrowing_y});
      }
      if (more_refinement) {
        refined_assumptions.clear();
        refined_assumptions.insert(new_refined_assumptions.begin(),
                                   new_refined_assumptions.end());
        CRAB_LOG("backward",
                 crab::outs()
                     << "Backward analysis refined forward analysis.\n");
      } else {
        CRAB_LOG(
            "backward",
            crab::outs()
                << "Backward analysis cannot refine more forward analysis.\n");
      }
      crab::CrabStats::stop("CombinedForwardBackward.CheckRefinement");

      if (!more_refinement || iters > max_num_iters) {
        if (more_refinement) {
          CRAB_LOG("backward",
                   crab::outs()
                       << "Limit of backward refinements already reached.\n");
        }

        // If refined_assumptions can prove that s is safe then
        // remove assertion s from m_unproven_assertions
        discharge_assertions(refined_assumptions, idom_tree);
        break;
      }

      // reuse wto for next iteration
      if (iters == 1) {
        m_b_wto = new bwd_wto_t(B.get_WTO());
      }

    } // end while true

    CRAB_VERBOSE_IF(1, get_msg_stream()
                           << "Combined forward+backward analysis done after "
                           << iters << " iterations.\n";);
  }

  // Return the invariants that hold at the entry of b
  AbsDom operator[](bb_label_t b) const { return get_pre(b); }

  // Return the invariants that hold at the entry of b
  AbsDom get_pre(bb_label_t b) const {
    auto it = m_pre_invariants.find(b);
    if (it == m_pre_invariants.end())
      return abs_dom_t::top();
    else
      return it->second;
  }

  // Return the invariants that hold at the exit of b
  AbsDom get_post(bb_label_t b) const {
    auto it = m_post_invariants.find(b);
    if (it == m_post_invariants.end())
      return abs_dom_t::top();
    else
      return it->second;
  }

  // Return the wto of the cfg
  const wto_t &get_wto() const {
    assert(m_wto);
    return *m_wto;
  }

  // clear internal state
  void clear() {
    m_pre_invariants.clear();
    m_post_invariants.clear();
  }

  /** Extra API for checkers **/

  CFG get_cfg(void) { return m_cfg; }

  std::shared_ptr<abs_tr_t> get_abs_transformer(AbsDom &&inv) {
    return std::make_shared<abs_tr_t>(std::move(inv));
  }

  void get_safe_assertions(std::set<const stmt_t *> &out) const {
    out.insert(m_proved_assertions.begin(), m_proved_assertions.end());
  }
};

} // namespace analyzer
} // namespace crab
