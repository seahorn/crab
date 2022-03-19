#pragma once

/**
 * Specialized fixpoint iterators for kill-gen problems.
 **/

//#include <crab/cfg/basic_block_traits.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>
#include <crab/cfg/cfg_bgl.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
// API for a kill-gen analysis operations
template <class CFG, class Dom> class killgen_operations_api {

public:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using killgen_domain_t = Dom;

protected:
  CFG m_cfg;

public:
  killgen_operations_api(CFG cfg) : m_cfg(cfg) {}

  virtual ~killgen_operations_api() {}

  // whether forward or backward analysis
  virtual bool is_forward() = 0;

  // initial state
  virtual Dom entry() = 0;

  // (optional) initialization for the fixpoint
  virtual void init_fixpoint() = 0;

  // confluence operator
  virtual Dom merge(Dom, Dom) = 0;

  // analyze a basic block
  virtual Dom analyze(const basic_block_label_t &, Dom) = 0;

  // analysis name
  virtual std::string name() = 0;
};

// A simple fixpoint for a killgen analysis
template <class CFG, class AnalysisOps> class killgen_fixpoint_iterator {

public:
  using basic_block_t = typename CFG::basic_block_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using killgen_domain_t = typename AnalysisOps::killgen_domain_t;
  using inv_map_t = std::unordered_map<basic_block_label_t, killgen_domain_t>;
  using iterator = typename inv_map_t::iterator;
  using const_iterator = typename inv_map_t::const_iterator;

protected:
  CFG m_cfg;
  inv_map_t m_in_map;
  inv_map_t m_out_map;

private:
   // to be constructed by the client so that m_analysis can have
   // arbitrary internal state.
  AnalysisOps &m_analysis;

  /**
   * run_bwd_fixpo(G) is in theory equivalent to
   * run_fwd_fixpo(reverse(G)).
   *  However, weak_rev_topo_sort(G) != weak_topo_sort(reverse(G))
   *  For instance, for a G=(V,E) where
   *   V= {v1,v2, v3, v4, v5},
   *   E= {(v1,v2), (v1,v3), (v2,v4), (v4,v1), (v3,v5)}
   *  (1) weak_rev_topo_sort(cfg)=[v5,v3,v4,v2,v1]
   *  (2) weak_topo_sort(reverse(cfg))=[v5,v3,v2,v4,v1] or even
   *      worse [v5,v3,v1,v4,v2] if vertices in the same scc are
   *      traversed in preorder.
   *  For a backward analysis, (1) will converge faster.
   *  For all of this, we decide not to reverse graphs and have
   *  two dual versions for the forward and backward analyses.
   **/

  void run_fwd_fixpo(std::vector<typename CFG::node_t> &order,
                     unsigned &iterations) {

    order = crab::analyzer::graph_algo::weak_topo_sort(m_cfg);
    assert((int)order.size() == std::distance(m_cfg.begin(), m_cfg.end()));
    bool change = true;
    iterations = 0;
    while (change) {
      change = false;
      ++iterations;
      for (unsigned i = 0, e = order.size(); i < e; ++i) {
        auto const &n = order[i];
        auto in = (i == 0 ? m_analysis.entry() : killgen_domain_t::bottom());
        for (auto const &p : m_cfg.prev_nodes(n))
          in = m_analysis.merge(in, m_out_map[p]);
        auto old_out = m_out_map[n];
        auto out = m_analysis.analyze(n, in);
        if (!(out <= old_out)) {
          m_out_map[n] = m_analysis.merge(out, old_out);
          change = true;
        } else
          m_in_map[n] = in;
      }
    }
  }

  void run_bwd_fixpo(std::vector<typename CFG::node_t> &order,
                     unsigned &iterations) {

    order = crab::analyzer::graph_algo::weak_rev_topo_sort(m_cfg);
    assert((int)order.size() == std::distance(m_cfg.begin(), m_cfg.end()));
    bool change = true;
    iterations = 0;
    while (change) {
      change = false;
      ++iterations;
      for (unsigned i = 0, e = order.size(); i < e; ++i) {
        auto const &n = order[i];
        auto out = (i == 0 ? m_analysis.entry() : killgen_domain_t::bottom());
        for (auto const &p : m_cfg.next_nodes(n))
          out = m_analysis.merge(out, m_in_map[p]);
        auto old_in = m_in_map[n];
        auto in = m_analysis.analyze(n, out);
        if (!(in <= old_in)) {
          m_in_map[n] = m_analysis.merge(in, old_in);
          change = true;
        } else
          m_out_map[n] = out;
      }
    }
  }

public:
  killgen_fixpoint_iterator(CFG cfg, AnalysisOps &analysis)
    : m_cfg(cfg), m_analysis(analysis) {
    // don't do anything with m_analysis in the contructor except
    // storing the reference.
  }

  void release_memory() {
    m_in_map.clear();
    m_out_map.clear();
  }

  void run() {
    crab::ScopedCrabStats __st__(m_analysis.name());

    m_analysis.init_fixpoint();

    std::vector<typename CFG::node_t> order;
    unsigned iterations = 0;

    if (m_analysis.is_forward()) {
      run_fwd_fixpo(order, iterations);
    } else {
      run_bwd_fixpo(order, iterations);
    }

    CRAB_LOG(m_analysis.name(),
	     crab::outs() << m_analysis.name();
	     if (m_cfg.has_func_decl()) {
	       crab::outs() << " for " << m_cfg.get_func_decl();
	     }
	     crab::outs() << "\n";
	     crab::outs() << "fixpoint ordering={";
             bool first = true; for (auto &v
                                     : order) {
               if (!first)
                 crab::outs() << ",";
               first = false;
               crab::outs() << basic_block_traits<basic_block_t>::to_string(v);
             } crab::outs() << "}\n";);

    CRAB_LOG(m_analysis.name(), crab::outs() << m_analysis.name() << ": "
                                             << "fixpoint reached in "
                                             << iterations << " iterations.\n");

    CRAB_LOG(
        m_analysis.name(), crab::outs() << m_analysis.name() << " sets:\n";
        for (auto it = m_cfg.label_begin(), et = m_cfg.label_end(); it != et;
             ++it) {
          crab::outs() << basic_block_traits<basic_block_t>::to_string(*it)
                       << " "
                       << "IN=" << m_in_map[*it] << " "
                       << "OUT=" << m_out_map[*it] << "\n";
        } crab::outs()
        << "\n";);
  }

  iterator in_begin() { return m_in_map.begin(); }
  iterator in_end() { return m_in_map.end(); }
  const_iterator in_begin() const { return m_in_map.begin(); }
  const_iterator in_end() const { return m_in_map.end(); }

  iterator out_begin() { return m_out_map.begin(); }
  iterator out_end() { return m_out_map.end(); }
  const_iterator out_begin() const { return m_out_map.begin(); }
  const_iterator out_end() const { return m_out_map.end(); }

  // return null if not found
  const killgen_domain_t *get_in(const basic_block_label_t &bb) const {
    auto it = m_in_map.find(bb);
    if (it != m_in_map.end()) {
      return &(it->second);
    } else {
      return nullptr;
    }
  }

  // return null if not found
  const killgen_domain_t *get_out(const basic_block_label_t &bb) const {
    auto it = m_out_map.find(bb);
    if (it != m_out_map.end()) {
      return &(it->second);
    } else {
      return nullptr;
    }
  }
};
} // end namespace crab
