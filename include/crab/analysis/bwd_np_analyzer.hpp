#pragma once

#include <crab/analysis/abs_transformer.hpp>
#include <crab/fixpoint/interleaved_fixpoint_iterator.hpp>
#include <crab/fixpoint/fixpoint_params.hpp>

#include <boost/range/iterator_range.hpp>
#include <unordered_map>

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
    : private ikos::interleaved_fwd_fixpoint_iterator<crab::cfg::cfg_rev<CFG>, AbsDom> {

  using fixpoint_iterator_t =
      ikos::interleaved_fwd_fixpoint_iterator<crab::cfg::cfg_rev<CFG>, AbsDom>;
  using basic_block_t = typename CFG::basic_block_t;
  using bb_label_t = typename CFG::basic_block_label_t;
  using stmt_t = typename CFG::statement_t;
  using bb_abstract_map_t = std::unordered_map<bb_label_t, AbsDom>;
  using pp_abstract_map_t = std::unordered_map<const stmt_t *, AbsDom>;

  /// forward and backward transformers
  using abs_fwd_tr_t = intra_abs_transformer<basic_block_t, AbsDom>;
  using abs_bwd_tr_t =
      intra_necessary_preconditions_abs_transformer<basic_block_t, AbsDom,
                                                    pp_abstract_map_t>;

  // the original CFG (i.e., not reversed)
  CFG m_cfg;
  // to create bottom/top abstract values
  AbsDom m_absval_fac;
  // necessary preconditions
  bb_abstract_map_t m_preconditions;
  // forward invariants
  bb_abstract_map_t m_invariants;
  // preconditions from good states, otherwise from bad states
  bool m_good_states;

  /**
   * Compute necessary preconditions for a basic block
   **/
  virtual AbsDom analyze(const bb_label_t &node, AbsDom &&precond) override {
    auto &bb = m_cfg.get_node(node);

    CRAB_LOG("backward-fixpoint",
             crab::outs() << "Post at "
                          << basic_block_traits<basic_block_t>::to_string(node)
                          << ": " << precond << "\n");

    // invariants that hold at the entry of the block
    AbsDom invariant = m_absval_fac.make_top();
    auto it = m_invariants.find(node);
    if (it != m_invariants.end()) {
      invariant = it->second;
    }
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
    CRAB_LOG("backward-fixpoint",
             crab::outs() << "Pre at "
                          << basic_block_traits<basic_block_t>::to_string(node)
                          << ": " << precond << "\n");
    return std::move(precond);
  }

  virtual void process_pre(const bb_label_t & /*node*/,
                           const AbsDom &/*postcond*/) override {}

  /**
   *  Store necessary preconditions
   **/
  virtual void process_post(const bb_label_t &node, const AbsDom &precond) override {
    AbsDom copy_precond(precond);
    m_preconditions.insert({node, std::move(copy_precond)});
  }

public:
  using wto_t = typename fixpoint_iterator_t::wto_t;
  using precond_map_t = bb_abstract_map_t;
  using iterator = typename precond_map_t::iterator;
  using const_iterator = typename precond_map_t::const_iterator;

  // This constructor computes necessary preconditions from error
  // states.
  necessary_preconditions_fixpoint_iterator(
      CFG cfg, AbsDom absval_fac,
      const fixpoint_parameters &fixpo_params)
    : fixpoint_iterator_t(crab::cfg::cfg_rev<CFG>(cfg), absval_fac, 
			  fixpo_params),
      m_cfg(cfg), m_absval_fac(absval_fac), m_good_states(false) {}

  // This constructor computes necessary preconditions from
  // safe/good (error) states if good_states is true (false).
  necessary_preconditions_fixpoint_iterator(
      CFG cfg, AbsDom absval_fac, bool good_states,
      const fixpoint_parameters &fixpo_params)
    : fixpoint_iterator_t(crab::cfg::cfg_rev<CFG>(cfg), absval_fac, 
			  fixpo_params),
      m_cfg(cfg), m_absval_fac(absval_fac), m_good_states(good_states) {}
  
  
  // postcond: final states that we want to propagate backwards  
  void run_backward(AbsDom postcond) { 
    this->run(postcond);
  }

  // postcond: final states that we want to propagate backwards  
  void run_backward(AbsDom postcond,
		    const std::unordered_map<typename CFG::basic_block_label_t,
		    AbsDom> &fwd_invariants) {
    m_invariants.insert(fwd_invariants.begin(), fwd_invariants.end());
    run_backward(postcond);
  }

  iterator begin() { return m_preconditions.begin(); }

  iterator end() { return m_preconditions.end(); }

  const_iterator begin() const { return m_preconditions.begin(); }

  const_iterator end() const { return m_preconditions.end(); }

  // return the preconditions at basic block node
  AbsDom operator[](const bb_label_t &node) const {
    auto it = m_preconditions.find(node);
    if (it != m_preconditions.end())
      return it->second;
    else
      return m_absval_fac.make_top();
  }

  // clear preconditions and forward invariants (if any)
  void clear() {
    m_preconditions.clear();
    m_invariants.clear();
  }

  const wto_t &get_wto() const { return this->get_wto(); }
};


} // namespace analyzer
} // namespace crab
