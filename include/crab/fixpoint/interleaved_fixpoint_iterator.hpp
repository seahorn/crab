/*******************************************************************************
 *
 * The interleaved fixpoint iterator is described in G. Amato and F. Scozzari's
 * paper: Localizing widening and narrowing. In Proceedings of SAS 2013,
 * pages 25-42. LNCS 7935, 2013.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 * Contributors: Jorge A. Navas (jorge.navas@sri.com)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#pragma once

#include <crab/cfg/cfg_bgl.hpp> // needed by wto
#include <crab/fixpoint/fixpoint_iterators_api.hpp>
#include <crab/fixpoint/fixpoint_params.hpp>
#include <crab/fixpoint/thresholds.hpp>
#include <crab/fixpoint/wto.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <unordered_map>

namespace ikos {

namespace interleaved_fwd_fixpoint_iterator_impl {

template <typename CFG, typename AbstractValue> class wto_iterator;

template <typename CFG, typename AbstractValue> class wto_processor;

} // namespace interleaved_fwd_fixpoint_iterator_impl

// for debugging only
template <typename CFG> static std::string func_name(CFG cfg) {
  if (cfg.has_func_decl()) {
    return cfg.get_func_decl().get_func_name();
  } else {
    return "";
  }
}

template <typename CFG, typename AbstractValue>
class interleaved_fwd_fixpoint_iterator
    : public fixpoint_iterator<CFG, AbstractValue> {

  friend class interleaved_fwd_fixpoint_iterator_impl::wto_iterator<
      CFG, AbstractValue>;

public:
  using basic_block_t = typename CFG::basic_block_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using wto_t = wto<CFG>;
  using assumption_map_t =
      std::unordered_map<basic_block_label_t, AbstractValue>;
  using invariant_table_t =
      std::unordered_map<basic_block_label_t, AbstractValue>;

private:
  using wto_iterator_t =
      interleaved_fwd_fixpoint_iterator_impl::wto_iterator<CFG, AbstractValue>;
  using wto_processor_t =
      interleaved_fwd_fixpoint_iterator_impl::wto_processor<CFG, AbstractValue>;
  using thresholds_t = crab::thresholds<typename CFG::number_t>;
  using wto_thresholds_t = crab::wto_thresholds<CFG>;
  using thresholds_map_t = typename wto_thresholds_t::thresholds_map_t;
  
protected:
  using iterator = typename invariant_table_t::iterator;
  using const_iterator = typename invariant_table_t::const_iterator;

  CFG m_cfg;
  wto_t m_wto;

  // We need to keep it so that we can call make_top(), make_bottom()
  AbstractValue m_absval_fac;
  // user-defined fixpoint parameters
  const crab::fixpoint_parameters &m_params;
  // set of thresholds to jump during widening
  thresholds_map_t m_thresholds_per_cycle;
  // enable post-processing of the invariants
  bool m_enable_processor;
  
private:
  // We don't want derived classes to access directly to m_pre and
  // m_post in case we make internal changes
  invariant_table_t m_pre, m_post;

  inline void set_pre(basic_block_label_t node, const AbstractValue &v) {
    crab::CrabStats::count("Fixpo.invariant_table.update");
    crab::ScopedCrabStats __st__("Fixpo.invariant_table.update");
    // To avoid calling the default constructor
    // m_pre[node] = v;

    auto res = m_pre.insert({node, v});
    if (!res.second) {
      // if the insertion failed then we just update value
      res.first->second = v;
    }
  }

  inline void set_post(basic_block_label_t node, AbstractValue &&v) {
    crab::CrabStats::count("Fixpo.invariant_table.update");
    crab::ScopedCrabStats __st__("Fixpo.invariant_table.update");
    // To avoid calling the default constructor
    // m_post[node] = std::move(v);

    auto it = m_post.find(node);
    if (it == m_post.end()) {
      m_post.insert({node, std::move(v)});
    } else {
      it->second = std::move(v);
    }
  }

  inline AbstractValue get(const invariant_table_t &table,
                           basic_block_label_t node) const {
    crab::CrabStats::count("Fixpo.invariant_table.lookup");
    crab::ScopedCrabStats __st__("Fixpo.invariant_table.lookup");
    return table.at(node);
  }

  inline AbstractValue extrapolate(basic_block_label_t node,
                                   unsigned int iteration,
                                   AbstractValue &before,
                                   AbstractValue &after) {
    crab::CrabStats::count("Fixpo.extrapolate");
    crab::ScopedCrabStats __st__("Fixpo.extrapolate");

    CRAB_VERBOSE_IF(
        1, crab::get_msg_stream()
               << "Widening " << iteration << " at " << func_name(m_cfg) << "::"
               << crab::basic_block_traits<basic_block_t>::to_string(node)
               << "\n";);

    if (iteration <= m_params.get_widening_delay()) {
      auto widen_res = before | after;
      CRAB_VERBOSE_IF(3, crab::outs() << "Prev   : " << before << "\n"
                                      << "Current: " << after << "\n"
                                      << "Res    : " << widen_res << "\n");
      return widen_res;
    } else {
      CRAB_VERBOSE_IF(3,
                      // To avoid closure on the left operand
                      AbstractValue before_copy(before);
                      crab::outs() << "Prev   : " << before_copy << "\n"
                                   << "Current: " << after << "\n");

      if (m_params.get_max_thresholds() > 0) {
        auto it = m_thresholds_per_cycle.find(node);
        if (it == m_thresholds_per_cycle.end()) {
          CRAB_ERROR("no thresholds found for ",
                     crab::basic_block_traits<basic_block_t>::to_string(node));
        }
        thresholds_t thresholds = it->second;
        auto widen_res = before.widening_thresholds(after, thresholds);
        CRAB_VERBOSE_IF(3,
                        // To avoid closure on the result
                        AbstractValue widen_res_copy(widen_res);
                        crab::outs() << "Res    : " << widen_res_copy << "\n");
        return widen_res;
      } else {
        auto widen_res = before || after;
        CRAB_VERBOSE_IF(3,
                        // To avoid closure on the result
                        AbstractValue widen_res_copy(widen_res);
                        crab::outs() << "Res    : " << widen_res_copy << "\n");
        return widen_res;
      }
    }
  }

  inline AbstractValue refine(basic_block_label_t node, unsigned int iteration,
                              AbstractValue &before, AbstractValue &after) {
    crab::CrabStats::count("Fixpo.refine");
    crab::ScopedCrabStats __st__("Fixpo.refine");

    CRAB_VERBOSE_IF(
        1, crab::get_msg_stream()
               << "Decreasing iteration=" << iteration << "\n"
               << "Narrowing at " << func_name(m_cfg) << "::"
               << crab::basic_block_traits<basic_block_t>::to_string(node)
               << "\n";);

    if (iteration == 1) {
      auto narrow_res = before & after;
      CRAB_VERBOSE_IF(3, crab::outs() << "Prev   : " << before << "\n"
                                      << "Current: " << after << "\n"
                                      << "Res    : " << narrow_res << "\n");
      return narrow_res;
    } else {
      auto narrow_res = before && after;
      CRAB_VERBOSE_IF(3, crab::outs() << "Prev   : " << before << "\n"
                                      << "Current: " << after << "\n"
                                      << "Res    : " << narrow_res << "\n");
      return narrow_res;
    }
  }

  void initialize_thresholds(size_t jump_set_size) {
    if (m_params.get_max_thresholds() > 0) {
      crab::CrabStats::resume("Fixpo");
      // select statically some widening points to jump to.
      wto_thresholds_t wto_thresholds(m_cfg, jump_set_size);
      m_wto.accept(&wto_thresholds);
      m_thresholds_per_cycle = wto_thresholds.get_thresholds_map();
      CRAB_VERBOSE_IF(2, crab::outs() << "Thresholds\n"
                                      << wto_thresholds << "\n");
      crab::CrabStats::stop("Fixpo");
    }
  }

  void initialize_invariant_tables() {
    clear();
    for (auto it = m_cfg.label_begin(), et = m_cfg.label_end(); it != et; ++it) {
      auto const &label = *it;
      m_pre.emplace(label, std::move(m_absval_fac.make_bottom()));
      m_post.emplace(label, std::move(m_absval_fac.make_bottom()));
    }
  }
  
public:
  interleaved_fwd_fixpoint_iterator(CFG cfg, AbstractValue absval_fac,
				    const crab::fixpoint_parameters &params,
                                    bool enable_processor = true)
      : m_cfg(cfg), m_wto(cfg), m_absval_fac(absval_fac),
        m_params(params),
        m_enable_processor(enable_processor) {
    initialize_thresholds(m_params.get_max_thresholds());
  }

  virtual ~interleaved_fwd_fixpoint_iterator() {}

  CFG get_cfg() const { return m_cfg; }

  wto_t &get_wto() { return m_wto; }
  const wto_t &get_wto() const { return m_wto; }

  /* Begin access methods for getting invariants */
  AbstractValue get_pre(basic_block_label_t node) const {
    return get(m_pre, node);
  }
  AbstractValue get_post(basic_block_label_t node) const {
    return get(m_post, node);
  }
  const invariant_table_t &get_pre_invariants() const { return m_pre; }
  const invariant_table_t &get_post_invariants() const { return m_post; }
  invariant_table_t &get_pre_invariants() { return m_pre; }
  invariant_table_t &get_post_invariants() { return m_post; }
  iterator pre_begin() { return m_pre.begin(); }
  iterator pre_end() { return m_pre.end(); }
  const_iterator pre_begin() const { return m_pre.begin(); }
  const_iterator pre_end() const { return m_pre.end(); }
  iterator post_begin() { return m_post.begin(); }
  iterator post_end() { return m_post.end(); }
  const_iterator post_begin() const { return m_post.begin(); }
  const_iterator post_end() const { return m_post.end(); }
  /* End access methods for getting invariants */

  void run(AbstractValue init) {
    crab::ScopedCrabStats __st__("Fixpo");

    initialize_invariant_tables();
    
    CRAB_VERBOSE_IF(1, crab::get_msg_stream() << "== Started analysis of "
                                              << func_name(m_cfg) << "\n");
    set_pre(m_cfg.entry(), init);
    wto_iterator_t iterator(this, m_absval_fac);
    m_wto.accept(&iterator);
    if (m_enable_processor) {
      wto_processor_t processor(this);
      m_wto.accept(&processor);
    }
    CRAB_VERBOSE_IF(1, crab::get_msg_stream() << "== Finished analysis of "
                                              << func_name(m_cfg) << "\n");
    CRAB_LOG("fixpo-trace", crab::get_msg_stream()
                                << "Fixpoint trace " << func_name(m_cfg) << ":\n"
                                << m_wto << "\n";);
  }

  void run(basic_block_label_t entry, AbstractValue init,
           const assumption_map_t &assumptions) {
    crab::ScopedCrabStats __st__("Fixpo");

    initialize_invariant_tables();
    
    CRAB_VERBOSE_IF(
        1, crab::get_msg_stream()
               << "== Started fixpoint at block " << func_name(m_cfg) << "::"
               << crab::basic_block_traits<basic_block_t>::to_string(entry)
               << " with initial value=" << init << "\n";);
    set_pre(entry, init);
    wto_iterator_t iterator(this, entry, m_absval_fac, &assumptions);
    m_wto.accept(&iterator);
    if (m_enable_processor) {
      wto_processor_t processor(this);
      m_wto.accept(&processor);
    }
    CRAB_VERBOSE_IF(1, crab::get_msg_stream() << "== Fixpoint reached for "
                                              << func_name(m_cfg) << "\n");
    CRAB_LOG("fixpo-trace", crab::get_msg_stream()
                                << "Fixpoint trace " << func_name(m_cfg) << ":\n"
                                << m_wto << "\n";);
  }

  void clear_pre() { m_pre.clear(); }

  void clear_post() { m_post.clear(); }

  void clear() {
    clear_pre();
    clear_post();
  }

}; // class interleaved_fwd_fixpoint_iterator

namespace interleaved_fwd_fixpoint_iterator_impl {

template <typename CFG, typename AbstractValue>
class wto_iterator : public wto_component_visitor<CFG> {

public:
  using interleaved_iterator_t =
      interleaved_fwd_fixpoint_iterator<CFG, AbstractValue>;
  using wto_vertex_t = wto_vertex<CFG>;
  using wto_cycle_t = wto_cycle<CFG>;
  using wto_t = wto<CFG>;
  using basic_block_t = typename CFG::basic_block_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using wto_nesting_t = typename wto_t::wto_nesting_t;
  using assumption_map_t = typename interleaved_iterator_t::assumption_map_t;

private:
  interleaved_iterator_t *m_iterator;
  // Initial entry point of the analysis
  basic_block_label_t m_entry;
  // To be able to create bottom and top abstract values
  const AbstractValue &m_absval_fac;
  const assumption_map_t *m_assumptions;
  // Used to skip the analysis until m_entry is found
  bool m_skip;

  inline AbstractValue make_top() const { return m_absval_fac.make_top(); }

  inline AbstractValue make_bottom() const { return m_absval_fac.make_bottom(); }

  inline AbstractValue strengthen(basic_block_label_t n, AbstractValue inv) {
    crab::CrabStats::count("Fixpo.strengthen");
    crab::ScopedCrabStats __st__("Fixpo.strengthen");

    if (m_assumptions) {
      auto it = m_assumptions->find(n);
      if (it != m_assumptions->end()) {
        CRAB_VERBOSE_IF(3, crab::outs() << "Before assumption at " << n << ":"
                                        << inv << "\n");

        inv = inv & it->second;
        CRAB_VERBOSE_IF(3, crab::outs() << "After assumption at " << n << ":"
                                        << inv << "\n");
      }
    }
    return inv;
  }

  inline void compute_post(basic_block_label_t node, AbstractValue inv) {
    crab::CrabStats::resume("Fixpo.analyze_block");
    CRAB_VERBOSE_IF(
        2, crab::get_msg_stream()
               << "Analyzing node " << func_name(m_iterator->m_cfg) << "::"
               << crab::basic_block_traits<basic_block_t>::to_string(node);
        auto &n = m_iterator->m_cfg.get_node(node);
        crab::outs() << " size=" << n.size() << "\n";);

    CRAB_VERBOSE_IF(4, crab::outs() << "PRE Invariants:\n" << inv << "\n");
    inv = m_iterator->analyze(node, std::move(inv));
    CRAB_VERBOSE_IF(3, crab::outs() << "POST Invariants:\n" << inv << "\n");
    crab::CrabStats::stop("Fixpo.analyze_block");

    m_iterator->set_post(node, std::move(inv));
  }

  // Simple visitor to check if node is a member of the wto component.
  class member_component_visitor : public wto_component_visitor<CFG> {
    basic_block_label_t _node;
    bool _found;

  public:
    member_component_visitor(basic_block_label_t node)
        : _node(node), _found(false) {}

    virtual void visit(wto_vertex_t &c) override {
      if (!_found) {
        _found = (c.node() == _node);
      }
    }

    virtual void visit(wto_cycle_t &c) override {
      if (!_found) {
        _found = (c.head() == _node);
        if (!_found) {
          for (typename wto_cycle_t::iterator it = c.begin(), et = c.end();
               it != et; ++it) {
            if (_found)
              break;
            it->accept(this);
          }
        }
      }
    }

    bool is_member() const { return _found; }
  };

public:
  wto_iterator(interleaved_iterator_t *iterator, const AbstractValue &absval_fac)
      : m_iterator(iterator), m_entry(m_iterator->get_cfg().entry()),
        m_absval_fac(absval_fac), m_assumptions(nullptr), m_skip(true) {}

  wto_iterator(interleaved_iterator_t *iterator, basic_block_label_t entry,
               const AbstractValue &absval_fac, const assumption_map_t *assumptions)
      : m_iterator(iterator), m_entry(entry), m_absval_fac(absval_fac),
        m_assumptions(assumptions), m_skip(true) {}

  virtual void visit(wto_vertex_t &vertex) override {
    basic_block_label_t node = vertex.node();

    /** decide whether skip vertex or not **/
    if (m_skip && (node == m_entry)) {
      m_skip = false;
    }
    if (m_skip) {
      CRAB_VERBOSE_IF(
          2, crab::outs() << "** Skipped analysis of  "
                          << func_name(m_iterator->m_cfg) << "::"
                          << crab::basic_block_traits<basic_block_t>::to_string(
                                 node)
                          << "\n");
      return;
    }

    AbstractValue pre = std::move(make_top());
    if (node == m_entry) {
      pre = m_iterator->get_pre(node);
      if (m_assumptions && !m_assumptions->empty()) {
        // no necessary but it might avoid copies
        pre = strengthen(node, pre);
        m_iterator->set_pre(node, pre);
      }
    } else {
      auto prev_nodes = m_iterator->m_cfg.prev_nodes(node);
      crab::CrabStats::resume("Fixpo.join_predecessors");
      pre = std::move(make_bottom());
      for (basic_block_label_t prev : prev_nodes) {
        pre |= m_iterator->get_post(prev);
      }
      crab::CrabStats::stop("Fixpo.join_predecessors");
      if (m_assumptions && !m_assumptions->empty()) {
        // no necessary but it might avoid copies
        pre = strengthen(node, pre);
      }
      m_iterator->set_pre(node, pre);
    }

    compute_post(node, pre);
  }

  virtual void visit(wto_cycle_t &cycle) override {
    basic_block_label_t head = cycle.head();

    auto get_nesting = [this](basic_block_label_t n) {
      boost::optional<wto_nesting_t> nesting = m_iterator->m_wto.nesting(n);
      if (nesting) {
        return *nesting;
      } else {
        CRAB_ERROR("WTO nesting: node ", n, " not found");
      }
    };

    /** decide whether skip cycle or not **/
    bool entry_in_this_cycle = false;
    if (m_skip) {
      // We only skip the analysis of cycle is m_entry is not a
      // component of it, included nested components.
      member_component_visitor vis(m_entry);
      cycle.accept(&vis);
      entry_in_this_cycle = vis.is_member();
      m_skip = !entry_in_this_cycle;
      if (m_skip) {
        CRAB_VERBOSE_IF(
            2, crab::outs()
                   << "** Skipped analysis of WTO cycle rooted at  "
                   << crab::basic_block_traits<basic_block_t>::to_string(head)
                   << "\n");
        return;
      }
    }

    CRAB_VERBOSE_IF(
        1, crab::get_msg_stream()
               << "** Analyzing loop with head "
               << func_name(m_iterator->m_cfg) << "::"
               << crab::basic_block_traits<basic_block_t>::to_string(head);
        auto &n = m_iterator->m_cfg.get_node(head);
        crab::outs() << " size=" << n.size() << "\n";);

    auto prev_nodes = m_iterator->m_cfg.prev_nodes(head);
    AbstractValue pre = std::move(make_bottom());
    wto_nesting_t cycle_nesting = get_nesting(head);

    if (entry_in_this_cycle) {
      CRAB_VERBOSE_IF(
          2, crab::outs() << "Skipped predecessors of "
                          << crab::basic_block_traits<basic_block_t>::to_string(
                                 head)
                          << "\n");
      pre = m_iterator->get_pre(m_entry);
    } else {
      crab::CrabStats::count("Fixpo.join_predecessors");
      crab::ScopedCrabStats __st__("Fixpo.join_predecessors");
      for (basic_block_label_t prev : prev_nodes) {
        if (!(get_nesting(prev) > cycle_nesting)) {
          pre |= m_iterator->get_post(prev);
        }
      }
    }
    if (m_assumptions && !m_assumptions->empty()) {
      // no necessary but it might avoid copies
      pre = strengthen(head, pre);
    }

    for (unsigned int iteration = 1;; ++iteration) {
      // keep track of how many times the cycle is visited by the fixpoint
      cycle.increment_fixpo_visits();

      // Increasing iteration sequence with widening
      m_iterator->set_pre(head, pre);
      compute_post(head, pre);
      for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end();
           ++it) {
        it->accept(this);
      }
      crab::CrabStats::resume("Fixpo.join_predecessors");
      AbstractValue new_pre = std::move(make_bottom());
      for (basic_block_label_t prev : prev_nodes) {
        new_pre |= m_iterator->get_post(prev);
      }
      crab::CrabStats::stop("Fixpo.join_predecessors");
      crab::CrabStats::resume("Fixpo.check_fixpoint");
      bool fixpoint_reached = new_pre <= pre;
      crab::CrabStats::stop("Fixpo.check_fixpoint");
      if (fixpoint_reached) {
        // Post-fixpoint reached
        CRAB_VERBOSE_IF(1, crab::get_msg_stream() << "post-fixpoint reached\n");
        m_iterator->set_pre(head, new_pre);
        pre = std::move(new_pre);
        break;
      } else {
        pre = m_iterator->extrapolate(head, iteration, pre, new_pre);
      }
    }

    if (m_iterator->m_params.get_descending_iterations() == 0) {
      // no narrowing
      return;
    }

    CRAB_VERBOSE_IF(1, crab::get_msg_stream() << "Started narrowing phase\n";);

    for (unsigned int iteration = 1;; ++iteration) {
      // Decreasing iteration sequence with narrowing
      compute_post(head, pre);
      for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end();
           ++it) {
        it->accept(this);
      }
      crab::CrabStats::resume("Fixpo.join_predecessors");
      AbstractValue new_pre = std::move(make_bottom());
      for (basic_block_label_t prev : prev_nodes) {
        new_pre |= m_iterator->get_post(prev);
      }
      crab::CrabStats::stop("Fixpo.join_predecessors");
      crab::CrabStats::resume("Fixpo.check_fixpoint");
      bool no_more_refinement = pre <= new_pre;
      crab::CrabStats::stop("Fixpo.check_fixpoint");
      if (no_more_refinement) {
        CRAB_VERBOSE_IF(1, crab::get_msg_stream()
                               << "No more refinement possible.\n");
        // No more refinement possible(pre == new_pre)
        break;
      } else {
        if (iteration > m_iterator->m_params.get_descending_iterations())
          break;
        pre = m_iterator->refine(head, iteration, pre, new_pre);
        m_iterator->set_pre(head, pre);
      }
    }
    CRAB_VERBOSE_IF(
        1,
        crab::get_msg_stream()
            << "** Finished loop with head " << func_name(m_iterator->m_cfg)
            << "::" << crab::basic_block_traits<basic_block_t>::to_string(head)
            << "\n");
  }

}; // class wto_iterator

template <typename CFG, typename AbstractValue>
class wto_processor : public wto_component_visitor<CFG> {

public:
  using interleaved_iterator_t =
      interleaved_fwd_fixpoint_iterator<CFG, AbstractValue>;
  using wto_vertex_t = wto_vertex<CFG>;
  using wto_cycle_t = wto_cycle<CFG>;
  using basic_block_label_t = typename CFG::basic_block_label_t;

private:
  interleaved_iterator_t *m_iterator;

public:
  wto_processor(interleaved_iterator_t *iterator) : m_iterator(iterator) {}

  virtual void visit(wto_vertex_t &vertex) override {
    crab::CrabStats::count("Fixpo.process_invariants");
    crab::ScopedCrabStats __st__("Fixpo.process_invariants");

    basic_block_label_t node = vertex.node();
    m_iterator->process_pre(node, m_iterator->get_pre(node));
    m_iterator->process_post(node, m_iterator->get_post(node));
  }

  virtual void visit(wto_cycle_t &cycle) override {
    crab::CrabStats::count("Fixpo.process_invariants");
    crab::ScopedCrabStats __st__("Fixpo.process_invariants");

    basic_block_label_t head = cycle.head();
    m_iterator->process_pre(head, m_iterator->get_pre(head));
    m_iterator->process_post(head, m_iterator->get_post(head));
    for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end();
         ++it) {
      it->accept(this);
    }
  }

}; // class wto_processor

} // namespace interleaved_fwd_fixpoint_iterator_impl
} // namespace ikos
