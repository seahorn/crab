#pragma once

#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
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
  using bb_reverse_iterator_t = typename CFG::basic_block_t::reverse_iterator;
  using cfg_iterator_t = typename CFG::iterator;
  using liveness_t = typename analyzer::liveness_analysis<CFG>;
  using varset_domain_t = typename liveness_t::varset_domain_t;

  int m_max_iterations;

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
          auto const &s_live_vars = s.get_live();

          // if DEF(s) = {d} and d \not\in out_live then remove s
	  if (!s.is_callsite()) {
	    // XXX: if s is a callsite the callee can still have
	    // side-effects (assertions).
	    // 
	    // TODO: use the assertion crawler analysis to remove the
	    // callee if it doesn't have any assertion.
	    if (s_live_vars.num_defs() == 1) {
	      variable_t def = *(s_live_vars.defs_begin());
	      if (!def.is_array_type()) {
		// XXX: An array variable contain actually multiple
		// definitions so we conservatively skip it.
		if (!(varset_domain_t(def) <= out_live)) {
		  // mark s to be removed
		  change = true;
		  apply_dce = true;
		  to_remove.push_back(&s);
		  CRAB_LOG("dce", crab::outs() << "DEAD: " << s << "\n";);
		}
	      }
	    }
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
