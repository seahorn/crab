#pragma once

#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/cfg/type_checker.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/fixpoint/fixpoint_params.hpp>
#include <crab/fixpoint/interleaved_fixpoint_iterator.hpp>

#include <algorithm>
#include <memory>

namespace crab {
namespace analyzer {
namespace analyzer_internal_impl {
/**
 * Only for crab internal use.
 * Implementation of an intra-procedural forward analysis.
 *
 * Perform a standard forward flow-sensitive analysis. AbsTr
 * defines the abstract transfer functions as well as which
 * operations are modeled.
 **/
template <typename CFG, typename AbsTr>
class fwd_analyzer : public ikos::interleaved_fwd_fixpoint_iterator<
                         CFG, typename AbsTr::abs_dom_t> {
public:
  using cfg_t = CFG;
  using basic_block_t = typename CFG::basic_block_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using varname_t = typename CFG::varname_t;
  using variable_t = typename CFG::variable_t;
  using number_t = typename CFG::number_t;
  using stmt_t = typename CFG::statement_t;
  using abs_dom_t = typename AbsTr::abs_dom_t;
  using abs_tr_t = AbsTr;

private:
  using fixpo_iterator_t =
      ikos::interleaved_fwd_fixpoint_iterator<CFG, abs_dom_t>;

public:
  using invariant_map_t = typename fixpo_iterator_t::invariant_table_t;
  using assumption_map_t = typename fixpo_iterator_t::assumption_map_t;
  using liveness_t = live_and_dead_analysis<CFG>;
  using wto_t = typename fixpo_iterator_t::wto_t;
  using iterator = typename fixpo_iterator_t::iterator;
  using const_iterator = typename fixpo_iterator_t::const_iterator;

private:
  using live_set_t = typename liveness_t::set_t;

  abs_tr_t *m_abs_tr; // the abstract transformer owned by the caller
  const liveness_t *m_live;
  live_set_t m_formals;

  inline abs_dom_t make_top() const {
    auto const &top_dom = m_abs_tr->get_abs_value();
    return top_dom.make_top();
  }

  void prune_dead_variables(const basic_block_label_t &node, abs_dom_t &inv) {
    if (!m_live) {
      return;
    }
    crab::ScopedCrabStats __st__("Pruning dead variables");
    if (inv.is_bottom() || inv.is_top()) {
      return;
    }
    auto dead = m_live->dead_exit(node);
    dead -= m_formals;
    std::vector<variable_t> dead_vec(dead.begin(), dead.end());
    inv.forget(dead_vec);
  }

  //! Given a basic block and the invariant at the entry it produces
  //! the invariant at the exit of the block.
  abs_dom_t analyze(const basic_block_label_t &node, abs_dom_t &&inv) override {
    auto &b = get_cfg().get_node(node);
    m_abs_tr->set_abs_value(std::move(inv));
    for (auto &s : b) {
      s.accept(&*m_abs_tr);
    }
    abs_dom_t &res = m_abs_tr->get_abs_value();
    prune_dead_variables(node, res);
    return res;
  }

  void process_pre(const basic_block_label_t &node, abs_dom_t inv) override {}

  void process_post(const basic_block_label_t &node, abs_dom_t inv) override {}

  void init_fwd_analyzer() {
    assert(m_abs_tr);
    CRAB_VERBOSE_IF(1, crab::outs() << "CFG with " << get_cfg().size()
                                    << " basic blocks\n";);
    CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking CFG ... ";);
    crab::CrabStats::resume("CFG type checking");
    crab::cfg::type_checker<CFG> tc(get_cfg());
    tc.run();
    crab::CrabStats::stop("CFG type checking");
    CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);

    if (::crab::CrabSanityCheckFlag) {
      // -- This sanity check typically flags whether a variable is
      //    used without proper definition.
      if (get_cfg().has_func_decl()) {
        crab::CrabStats::resume("Live symbols sanity check");
        CRAB_VERBOSE_IF(1, get_msg_stream()
                               << "Live symbols sanity check ... ";);
        liveness_analysis<CFG> live_symbols(get_cfg(), false /* keep IN sets*/);
        live_symbols.exec();
        if (auto const *entry_ls = live_symbols.get_in(get_cfg().entry())) {
          auto const &fdecl = get_cfg().get_func_decl();
          typename liveness_analysis<CFG>::varset_domain_t suspicious_vars(
              *entry_ls);
          for (unsigned i = 0; i < fdecl.get_num_inputs(); i++) {
            suspicious_vars -= fdecl.get_input_name(i);
          }
          if (!suspicious_vars.is_bottom()) {
            crab::outs() << "\n*** Sanity check failed: " << suspicious_vars
                         << " might not be initialized  in "
                         << get_cfg().get_func_decl().get_func_name() 
			 << "\n    This lack of initialization can be legitimate with region statements "
			 << "such as make_ref, gep_ref, int_to_ref "
			 << "if the used regions are not read or modified in the code under analysis.\n";
          } else {
            CRAB_VERBOSE_IF(1, crab::outs() << "OK";);
          }
        }
        CRAB_VERBOSE_IF(1, crab::outs() << "\n";);
        crab::CrabStats::stop("Live symbols sanity check");
      }
    }

    if (m_live) {
      // --- collect input and output parameters for later use
      if (get_cfg().has_func_decl()) {
        auto const &fdecl = get_cfg().get_func_decl();
        for (unsigned i = 0; i < fdecl.get_num_inputs(); i++)
          m_formals += fdecl.get_input_name(i);
        for (unsigned i = 0; i < fdecl.get_num_outputs(); i++)
          m_formals += fdecl.get_output_name(i);
      }
    }    
  }
  
public:
  fwd_analyzer(CFG cfg, abs_tr_t *abs_tr,
	       // to create top and bottom abstract values
	       abs_dom_t absval_fac,
               // live can be nullptr if no live info is available
               const liveness_t *live_and_dead_symbols,
	       const fixpoint_parameters &params)
    : fixpo_iterator_t(cfg, absval_fac, params,
		       false /*disable processor*/),
      m_abs_tr(abs_tr), m_live(live_and_dead_symbols) {
    init_fwd_analyzer();
  }

  fwd_analyzer(const fwd_analyzer &o) = delete;

  fwd_analyzer &operator=(const fwd_analyzer &o) = delete;

  //! Trigger the fixpoint computation
  void run_forward(abs_dom_t init) {
    this->run(init);
  }

  void run_forward(const basic_block_label_t &entry,
		   abs_dom_t init,
                   const assumption_map_t &assumptions) {
    this->run(entry, init, assumptions);
  }

  //! Return the invariants that hold at the entry of b
  inline abs_dom_t operator[](const basic_block_label_t &b) const {
    return get_pre(b);
  }

  //! Return the invariants that hold at the entry of b
  abs_dom_t get_pre(const basic_block_label_t &b) const {
    return fixpo_iterator_t::get_pre(b);
  }

  //! Return the invariants that hold at the exit of b
  abs_dom_t get_post(const basic_block_label_t &b) const {
    return fixpo_iterator_t::get_post(b);
  }

  //! Return the WTO of the CFG. The WTO contains also how many
  //! times each head was visited by the fixpoint iterator.
  wto_t &get_wto() { return fixpo_iterator_t::get_wto(); }
  const wto_t &get_wto() const { return fixpo_iterator_t::get_wto(); }

  void clear() {
    fixpo_iterator_t::clear_pre();
    fixpo_iterator_t::clear_post();
  }

  CFG get_cfg() const { return this->m_cfg; }

  abs_tr_t &get_abs_transformer() { return *m_abs_tr; }

  void get_safe_assertions(std::set<const stmt_t *> &out) const {}

  void write(crab::crab_os &o) {
    for (auto &kv : this->get_pre_invariants()) {
      auto &pre = kv.second;
      auto &b = get_cfg().get_node(kv.first);
      o << basic_block_traits<basic_block_t>::to_string(kv.first) << ":\n";
      if (b.size() == 0) {
        o << "/**\n"
          << "* " << pre << "\n"
          << "**/\n";
      } else {
        auto post = get_post(kv.first);
        o << "**\n"
          << "* PRE: " << pre << "\n"
          << "**/\n";
        for (auto &s : b) {
          o << "\t" << s << "\n";
        }
        o << "/**\n"
          << "* POST: " << post << "\n"
          << "**/\n";
      }
    }
  }
};

/**
 * Wrapper for fwd_analyzer_class.
 *
 * The main difference with fwd_analyzer class is that here we
 * create an abstract transformer instance while fwd_analyzer does
 * not.
 **/
template <typename CFG, typename AbsDomain, typename AbsTr>
class intra_fwd_analyzer_wrapper {
  using fwd_analyzer_t = fwd_analyzer<CFG, AbsTr>;

public:
  using abs_dom_t = AbsDomain;
  using liveness_t = live_and_dead_analysis<CFG>;
  using cfg_t = CFG;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;
  using stmt_t = typename CFG::statement_t;
  using abs_tr_t = typename fwd_analyzer_t::abs_tr_t;
  using wto_t = typename fwd_analyzer_t::wto_t;
  using assumption_map_t = typename fwd_analyzer_t::assumption_map_t;
  using invariant_map_t = typename fwd_analyzer_t::invariant_map_t;
  using iterator = typename fwd_analyzer_t::iterator;
  using const_iterator = typename fwd_analyzer_t::const_iterator;

private:
  std::unique_ptr<abs_tr_t> m_abs_tr;
  fwd_analyzer_t m_analyzer;

public:
  intra_fwd_analyzer_wrapper(CFG cfg,
			     // create bottom and top instances
			     AbsDomain absval_fac,
                             // live variables
                             const liveness_t *live,
			     const fixpoint_parameters &fixpo_params)
    : m_abs_tr(new abs_tr_t(absval_fac.make_top())),
      m_analyzer(cfg, &*m_abs_tr, absval_fac, live,  fixpo_params) {}
  
  intra_fwd_analyzer_wrapper(const intra_fwd_analyzer_wrapper &o) = delete;

  intra_fwd_analyzer_wrapper &
  operator=(const intra_fwd_analyzer_wrapper &o) = delete;

  void run(abs_dom_t init) {
    m_analyzer.run_forward(init);
  }

  void run(const basic_block_label_t &entry,
	   abs_dom_t init,
           const assumption_map_t &assumptions) {
    m_analyzer.run_forward(entry, init, assumptions);
  }

  iterator pre_begin() { return m_analyzer.pre_begin(); }
  iterator pre_end() { return m_analyzer.pre_end(); }
  const_iterator pre_begin() const { return m_analyzer.pre_begin(); }
  const_iterator pre_end() const { return m_analyzer.pre_end(); }

  iterator post_begin() { return m_analyzer.post_begin(); }
  iterator post_end() { return m_analyzer.post_end(); }
  const_iterator post_begin() const { return m_analyzer.post_begin(); }
  const_iterator post_end() const { return m_analyzer.post_end(); }

  const invariant_map_t &get_pre_invariants() const {
    return m_analyzer.get_pre_invariants();
  }

  const invariant_map_t &get_post_invariants() const {
    return m_analyzer.get_post_invariants();
  }

  abs_dom_t operator[](const basic_block_label_t &b) const {
    return m_analyzer.get_pre(b);
  }

  abs_dom_t get_pre(const basic_block_label_t &b) const {
    return m_analyzer.get_pre(b);
  }

  abs_dom_t get_post(const basic_block_label_t &b) const {
    return m_analyzer.get_post(b);
  }

  wto_t &get_wto() { return m_analyzer.get_wto(); }
  const wto_t &get_wto() const { return m_analyzer.get_wto(); }

  void clear() { m_analyzer.clear(); }

  CFG get_cfg() { return m_analyzer.get_cfg(); }

  abs_tr_t &get_abs_transformer() {
    assert(m_abs_tr);
    return *m_abs_tr;
  }

  void get_safe_assertions(std::set<const stmt_t *> &out) const {}
};
} // namespace analyzer_internal_impl

/**
 * External api
 **/
template <typename CFG, typename AbsDomain>
using intra_fwd_analyzer = analyzer_internal_impl::intra_fwd_analyzer_wrapper<
    CFG, AbsDomain,
    intra_abs_transformer<typename CFG::basic_block_t, AbsDomain>>;

} // namespace analyzer
} // namespace crab
