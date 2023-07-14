#pragma once

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/bwd_np_analyzer.hpp>
#include <crab/analysis/fwd_bwd_params.hpp>
#include <crab/analysis/graphs/dominance.hpp>

#include <boost/range/iterator_range.hpp>
#include <algorithm>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

namespace crab {

namespace analyzer {

#define FWD_BACK_FWD_PHASE "Fwd+Back.FwdPhase"
#define FWD_BACK_BACK_PHASE "Fwd+Back.BackPhase"
#define FWD_BACK_MIN_INVARIANTS "Fwd+Back.MinInv"
#define FWD_BACK_CHECK_REFINE "Fwd+Back.CheckRefine"
#define FWD_BACK_DOM_TREE "Fwd+Back.DomTree"  
const std::string FWD_BACK_TOTAL_RUNS("Fwd+Back.totalruns");  
const std::string FWD_BACK_ITERATIONS("Fwd+Back.iterations");

#define FWD_BACK_TIMER_STATS(NAME) CRAB_SCOPED_TIMER_STATS(NAME, 0)
#define FWD_BACK_COUNT_STATS(NAME) CRAB_COUNT_STATS(NAME, 0)

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
  using cfg_t = CFG;
  using basic_block_t = typename CFG::basic_block_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using statement_t = typename CFG::statement_t;
  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;
  // used for checkers
  using abs_dom_t = AbsDom;
  using abs_tr_t = intra_abs_transformer<basic_block_t, AbsDom>;

private:
  using fwd_analyzer_t = intra_fwd_analyzer<CFG, AbsDom>;
  using bwd_analyzer_t =
      necessary_preconditions_fixpoint_iterator<CFG, AbsDom>;
  using precond_map_t = typename bwd_analyzer_t::precond_map_t;
  using liveness_t = live_and_dead_analysis<CFG>;
  using idom_tree_t =
      std::unordered_map<basic_block_label_t, std::set<basic_block_label_t>>;
  using assert_t = typename basic_block_t::assert_t;
  using bool_assert_t = typename basic_block_t::bool_assert_t;
  using assert_ref_t = typename basic_block_t::assert_ref_t;

public:
  using assumption_map_t = typename fwd_analyzer_t::assumption_map_t;
  using invariant_map_t = typename fwd_analyzer_t::invariant_map_t;
  using wto_t = typename fwd_analyzer_t::wto_t;

private:
  // the cfg
  CFG m_cfg;
  // create bottom/top abstract values
  AbsDom m_absval_fac;
  // keep track of which assertions cannot be proven.
  std::vector<std::pair<basic_block_label_t, statement_t *>>
      m_unproven_assertions;
  // keep track of which assertions have been proven
  std::set<statement_t *> m_proved_assertions;
  // keep the results of the first forward iteration.
  mutable invariant_map_t m_pre_invariants;
  mutable invariant_map_t m_post_invariants;
  // The forward wto is kept to be exposed to outside clients.
  std::unique_ptr<wto_t> m_wto;
  // to be used by checker
  std::unique_ptr<abs_tr_t> m_abs_tr;


  void store_results(fwd_analyzer_t &f) {
    for (auto &kv : f.get_pre_invariants()) {
      m_pre_invariants.insert({kv.first,kv.second});
    }

    for (auto &kv : f.get_post_invariants()) {
      m_post_invariants.insert({kv.first, kv.second});
    }
    m_wto = std::unique_ptr<wto_t>(new wto_t(f.get_wto().clone()));
  }

  void gather_assertions() {
    for (auto it = m_cfg.begin(), et = m_cfg.end(); it != et; ++it) {
      for (auto &s : *it) {
        if (s.is_assert() || s.is_bool_assert() || s.is_ref_assert()) {
          m_unproven_assertions.push_back({it->label(), &s});
        }
      }
    }
  }

  // idom induces a tree each key-value pair means that key
  // (block) is an immediate dominator of each element in value
  // (set of basic blocks). TODO: caching
  bool dominates(const basic_block_label_t &u, const basic_block_label_t &v,
                 const idom_tree_t &idom) {
    auto it = idom.find(u);
    if (it == idom.end()) {
      // u does not dominate any block
      return false;
    }

    // immediate dominated by u
    const std::set<basic_block_label_t> &idom_u = it->second;
    if (idom_u.count(v) > 0) {
      // u dominates v
      return true;
    }

    for (const basic_block_label_t &w : idom_u) {
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
        const basic_block_label_t &n = it->label();
        auto jt = refined_assumptions.find(n);
        if (jt == refined_assumptions.end()) {
          continue;
        }
        if (jt->second.is_bottom()) {
          m_unproven_assertions.erase(
              std::remove_if(
                  m_unproven_assertions.begin(), m_unproven_assertions.end(),
                  [&n, &idom,
                   this](std::pair<basic_block_label_t, statement_t *> &kv) {
                    const basic_block_label_t &m = kv.first;
                    statement_t *s = kv.second;
                    bool res = this->dominates(n, m, idom);
                    if (res) {
		      CRAB_LOG("backward",
			       crab::outs() << "Backward analysis proved " << *s
			       << " because " << n << " dominates " << m << "\n";);
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
  using iterator = typename bwd_analyzer_t::iterator;
  using const_iterator = typename bwd_analyzer_t::const_iterator;

  intra_forward_backward_analyzer(CFG cfg, AbsDom absval_fac)
    : m_cfg(cfg), m_absval_fac(absval_fac), m_wto(nullptr),
      m_abs_tr(new abs_tr_t(absval_fac.make_top())) {} 
        
  ~intra_forward_backward_analyzer() = default;
  intra_forward_backward_analyzer(const intra_forward_backward_analyzer &o) = delete;
  intra_forward_backward_analyzer &
  operator=(const intra_forward_backward_analyzer &o) = delete;

  /**
   * Perform the refining forward-backward loop.
   **/
  void run(AbsDom init_states,
           // assumptions
           const assumption_map_t &assumptions,
           // liveness information
           const liveness_t *live,
           // users parameters 
	   const fixpoint_parameters &fixpo_params,
	   const fwd_bwd_parameters &params) {
    run(m_cfg.entry(), init_states, assumptions, live, fixpo_params, params);
  }

  void run(const basic_block_label_t &entry, // only used for the forward pass.
           AbsDom init_states,
           // assumptions
           const assumption_map_t &assumptions,
           // liveness information
           const liveness_t *live,
           // user parameters 
	   const fixpoint_parameters &fixpo_params,
	   const fwd_bwd_parameters &params) {

    CRAB_LOG("backward", crab::outs() << "Running forward+backward analysis ";
	     if (m_cfg.has_func_decl()) {
	       auto const& fdecl = m_cfg.get_func_decl();
	       crab::outs() << " on " << fdecl << " ";
	     }
	     crab::outs() << "with initial states=" << init_states << "\n");

    // return true if fixpo.get_pre(node) is strictly more precise than
    // old fixpo.get_pre(node)
    auto refine =
        [](const basic_block_label_t &node, const assumption_map_t &old_table,
           const bwd_analyzer_t &fixpo, assumption_map_t &new_table) {
          AbsDom new_val = fixpo.get_pre(node);
          auto it = old_table.find(node);
          if (it == old_table.end()) {
	    CRAB_LOG("backward-refinement",
		     crab::outs() << "New assumptions for fwd analysis at " <<
		     basic_block_traits<basic_block_t>::to_string(node) << "\n";
		     crab::outs() << "New value(from backward)=" << new_val << "\n";
		     );
            new_table.insert({node, std::move(new_val)});	    
            return true;
          } else {
            AbsDom old_val = it->second;
            AbsDom refined_val = old_val && new_val;
	    auto res = (!(old_val <= refined_val));
	    CRAB_LOG("backward-refinement",
		     if (res) {
		       crab::outs() << "Refined assumptions for fwd analysis at " <<
		       basic_block_traits<basic_block_t>::to_string(node) << "\n";
		       crab::outs() << "Old value=" << old_val << "\n";
		       crab::outs() << "New value(from backward)=" << new_val << "\n";
		       crab::outs() << "Refined value=" << refined_val << "\n";
		     });
            new_table.insert({node, std::move(refined_val)});
	    return res;
	  }
        };

    bool only_forward = !params.is_enabled_backward();
    if (!only_forward && !m_cfg.has_exit()) {
      CRAB_WARN("cannot run backward analysis because CFG has no exit block");
      only_forward = true;
    }

    FWD_BACK_COUNT_STATS(FWD_BACK_TOTAL_RUNS);

    // number of refinement iterations
    unsigned iters = 0;
    // CFG assertions
    std::vector<std::pair<basic_block_label_t, statement_t *>> assertions;
    // immediate dominance tree
    idom_tree_t idom_tree;

    gather_assertions();
    CRAB_LOG("backward", crab::outs()
                             << "Found " << m_unproven_assertions.size()
                             << " assertions.\n";);

    if (!m_unproven_assertions.empty() && !only_forward) {
      //crab::ScopedCrabStats __st__(FWD_BACK_DOM_TREE, false);
      std::unordered_map<basic_block_label_t, basic_block_label_t> idom_map;
      crab::analyzer::graph_algo::dominator_tree(m_cfg, m_cfg.entry(),
                                                 idom_map);
      // build idom_tree
      for (auto &kv : idom_map) {
        if (kv.second == boost::graph_traits<CFG>::null_vertex()) {
          continue;
        }
        auto it = idom_tree.find(kv.second);
        if (it == idom_tree.end()) {
          std::set<basic_block_label_t> reachable({kv.first});
          idom_tree.insert({kv.second, reachable});
        } else {
          it->second.insert(kv.first);
        }
      }

      CRAB_LOG("backward-dom-tree", crab::outs() << "Computed dominance tree:\n";
               for (auto &kv
                    : idom_tree) {
                 crab::outs()
                     << "\t"
                     << basic_block_traits<basic_block_t>::to_string(kv.first)
                     << " dominates={";
                 for (auto d : kv.second) {
                   crab::outs() << d << ";";
                 }
                 crab::outs() << "}\n";
               });
    }
    assumption_map_t refined_assumptions(assumptions.begin(), assumptions.end());

    fwd_analyzer_t F(m_cfg, m_absval_fac, live, fixpo_params);
    std::unique_ptr<bwd_analyzer_t> B = nullptr;
    while (true) {
      iters++;
      FWD_BACK_COUNT_STATS(FWD_BACK_ITERATIONS);
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Iteration " << iters << "\n"
                                          << "Started forward analysis.\n";);

      F.clear();
      if(B) {
	B->clear();
      }
      
      // run forward analysis refined with preconditions from error
      // states
      {
	FWD_BACK_TIMER_STATS(FWD_BACK_FWD_PHASE);
	F.run(entry, init_states, refined_assumptions);
      }
      if (!params.get_use_refined_invariants() && iters == 1) {
        store_results(F);
      }

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished forward analysis.\n";);

      CRAB_LOG(
          "backward-fwd-results", crab::outs() << "Forward analysis: \n";
          for (auto &kv
               : boost::make_iterator_range(F.pre_begin(), F.pre_end())) {
            crab::outs() << basic_block_traits<basic_block_t>::to_string(
                                kv.first)
                         << ":\n"
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

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Started backward analysis.\n";);
      if (B == nullptr) {
	B = std::unique_ptr<bwd_analyzer_t>(new bwd_analyzer_t
			    (m_cfg, m_absval_fac, 
			     // negate assertions: preconditions from error states
			     false, fixpo_params));
      }

      // run backward analysis computing necessary preconditions
      // refined with results from the forward analysis.
      const invariant_map_t &fwd_invariants = F.get_pre_invariants();
      // Important for apron and elina domains.
      invariant_map_t minimized_fwd_invariants;
      for (auto &kv : fwd_invariants) {
	FWD_BACK_TIMER_STATS(FWD_BACK_MIN_INVARIANTS);
        AbsDom dom(kv.second);
        dom.minimize();
        minimized_fwd_invariants.insert({kv.first, dom});
      }
      // A final state is safe so we start from bottom
      {
	FWD_BACK_TIMER_STATS(FWD_BACK_BACK_PHASE);      
	B->run_backward(m_absval_fac.make_bottom(), minimized_fwd_invariants);
      }

      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished backward analysis.\n";);

      CRAB_LOG(
          "backward-bwd-results", crab::outs() << "Backward analysis:\n";
          for (auto &kv
               : boost::make_iterator_range(B->begin(), B->end())) {
            crab::outs() << basic_block_traits<basic_block_t>::to_string(
                                kv.first)
                         << ":\n"
                         << kv.second << "\n";
          } crab::outs()
          << "\n");

      assumption_map_t new_refined_assumptions;
      bool more_refinement = false;
      for (auto it = m_cfg.begin(), et = m_cfg.end(); it != et; ++it) {
	FWD_BACK_TIMER_STATS(FWD_BACK_CHECK_REFINE);
        more_refinement |= refine(it->label(), refined_assumptions, *B,
                                  new_refined_assumptions);
      }
      if (more_refinement) {
        refined_assumptions.clear();
        refined_assumptions.insert(new_refined_assumptions.begin(),
                                   new_refined_assumptions.end());
        CRAB_LOG("backward",
		 if (iters > 1) {
		   crab::outs() << "Backward analysis refined forward analysis.\n";
		 });
      } else {
        CRAB_LOG(
            "backward",
            crab::outs()
                << "Backward analysis cannot refine more forward analysis.\n");
      }

      if (!more_refinement || iters > params.get_max_refine_iterations()) {
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
    } // end while true


    if (params.get_use_refined_invariants()) {
      // We store the refined (last iteration) forward invariants 
      store_results(F);
    }
    
    CRAB_VERBOSE_IF(1, get_msg_stream()
                           << "Combined forward+backward analysis done after "
                           << iters << " iterations.\n";);
  }

  // Return the invariants that hold at the entry of b
  const AbsDom& operator[](const basic_block_label_t &b) const { return get_pre(b); }

  // Return the invariants that hold at the entry of b
  const AbsDom& get_pre(const basic_block_label_t &b) const {
    auto it = m_pre_invariants.find(b);
    if (it == m_pre_invariants.end()) {
      it->second = m_absval_fac.make_top();
    }
    return it->second;
  }

  // Return the invariants that hold at the exit of b
  const AbsDom& get_post(const basic_block_label_t &b) const {
    auto it = m_post_invariants.find(b);
    if (it == m_post_invariants.end()) {
      it->second = m_absval_fac.make_top();
    }
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
    m_proved_assertions().clear();
    m_unproven_assertions().clear();
    m_abs_tr.get_abs_value().set_to_top();
  }

  /** Extra API for checkers **/

  CFG get_cfg(void) { return m_cfg; }

  abs_tr_t &get_abs_transformer() {
    assert(m_abs_tr);
    return *m_abs_tr;
  }

  void get_safe_assertions(std::set<const statement_t *> &out) const {
    out.insert(m_proved_assertions.begin(), m_proved_assertions.end());
  }
};

} // namespace analyzer
} // namespace crab
