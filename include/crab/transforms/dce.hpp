#pragma once

#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/transforms/transform.hpp>

namespace crab {
namespace transforms {

/**
 * Naive implementation of dead code elimination.
 *
 * The current version is just a proof of concept and it won't
 * scale. We need an incremental liveness analysis so that we can
 * update liveness while we remove statements and avoid running
 * liveness analysis from scratch each time.
 **/
template <typename CFG> class dead_code_elimination : public transform<CFG> {
  // some typedefs here
  using variable_t = typename CFG::variable_t;
  using statement_t = typename CFG::statement_t;
  using basic_block_t = typename CFG::basic_block_t;
  using live_t = typename statement_t::live_t;
  using bb_reverse_iterator_t = typename CFG::basic_block_t::reverse_iterator;
  using cfg_iterator_t = typename CFG::iterator;
  using liveness_t = typename analyzer::liveness_analysis<CFG>;
  using varset_domain_t = typename liveness_t::varset_domain_t;

  int m_max_iterations;

  bool keep_conservatively(const statement_t &s) const {
    if (s.is_callsite()) {
      // If s is a callsite then the callee can still have
      // side-effects (assertions).
      //
      // TODO: use the assertion crawler analysis to remove the callee
      // if it doesn't have any assertion.
      return true;
    } else if (s.is_arr_init() || s.is_arr_write()) {
      return true;
    } else if (s.is_ref_store() || 
	       s.is_region_init() ||
	       s.is_ref_remove()) {
      return true;
    } else if (s.is_assert() || s.is_ref_assert() || s.is_bool_assert()) {
      return true;
    } else {
      return false;
    }
  }

  bool included_defs(const live_t &live, const varset_domain_t &vars) const {
    for (auto it=live.defs_begin(), et=live.defs_end(); it!=et; ++it) {
      if (!(varset_domain_t(*it) <= vars)) {
	return false;
      }
    }
    return true;
  }
  
public:
  dead_code_elimination(int max_iterations = 10)
      : transform<CFG>(), m_max_iterations(max_iterations) {}

  virtual bool run(CFG &cfg) override {
    bool apply_dce = false;
    bool change;
    int n = m_max_iterations;
    do {
      change = false;
      liveness_t live(cfg);
      live.exec();
      for (cfg_iterator_t b_it = cfg.begin(), b_et = cfg.end(); b_it != b_et;
           ++b_it) {
        basic_block_t &bb = *b_it;
        varset_domain_t out_live = live.get(bb.label());
        std::vector<const statement_t *> to_remove;
        for (bb_reverse_iterator_t s_it = bb.rbegin(), s_et = bb.rend();
             s_it != s_et; ++s_it) {
          statement_t &s = *s_it;
          auto const& s_live_vars = s.get_live();
          // if DEF(s) \not\subseteq out_live then remove s
          if (!keep_conservatively(s) && !(included_defs(s_live_vars, out_live))) {
	    // mark s to be removed
	    change = true;
	    apply_dce = true;
	    to_remove.push_back(&s);
	    CRAB_LOG("dce", crab::outs() << "DEAD: " << s << "\n";);
          }
          // update out_live for the next statement:
          //   out_live = (out_live \ DEFS(s)) U USES(s)
          for (auto it = s_live_vars.defs_begin(), et = s_live_vars.defs_end();
               it != et; ++it) {
            out_live -= *it;
          }
          for (auto it = s_live_vars.uses_begin(), et = s_live_vars.uses_end();
               it != et; ++it) {
            out_live += *it;
          }
        } // end inner for

        // remove the dead statements from the block
        while (!to_remove.empty()) {
          const statement_t *s = to_remove.back();
          to_remove.pop_back();
          bb.remove(s, false /* do not update liveness yet*/);
        }
        // we update uses and defs sets for the basic block once we
        // have removed all dead statements.
        bb.update_uses_and_defs();
      }
      --n;

    } while (change && n > 0);

    return apply_dce;
  }

  virtual std::string get_name() const override {
    return "Dead Code Elimination";
  }
};

} // namespace transforms
} // namespace crab
