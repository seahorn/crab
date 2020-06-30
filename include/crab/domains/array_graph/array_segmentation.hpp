#pragma once

/*
   Quick analysis to collect the set of variables and constants that might
   define array segment boundaries.
*/

#include <crab/cfg/cfg.hpp>
#include <crab/domains/killgen_domain.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>

namespace crab {

namespace analyzer {

template <class CFG>
class array_segment_ops
    : public crab::iterators::killgen_operations_api<
          CFG, domains::flat_killgen_domain<typename CFG::variable_t>> {

  using killgen_operations_api_t = crab::iterators::killgen_operations_api<
      CFG, domains::flat_killgen_domain<typename CFG::variable_t>>;

public:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using N = typename CFG::number_t;
  using V = typename CFG::varname_t;
  using variable_t = typename CFG::variable_t;
  using array_segment_domain_t =
      typename killgen_operations_api_t::killgen_domain_t;

private:
  class array_segment_visitor
      : public crab::cfg::statement_visitor<basic_block_label_t, N, V> {
    using bin_op_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                           N, V>::bin_op_t;
    using assign_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                           N, V>::assign_t;
    using assume_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                           N, V>::assume_t;
    using select_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                           N, V>::select_t;
    using assert_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                           N, V>::assert_t;

    using havoc_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                          N, V>::havoc_t;
    using unreach_t = typename crab::cfg::statement_visitor<basic_block_label_t,
                                                            N, V>::unreach_t;
    using callsite_t =
        typename crab::cfg::statement_visitor<basic_block_label_t, N,
                                              V>::callsite_t;

    using arr_init_t =
        typename crab::cfg::statement_visitor<basic_block_label_t, N,
                                              V>::arr_init_t;
    using arr_load_t =
        typename crab::cfg::statement_visitor<basic_block_label_t, N,
                                              V>::arr_load_t;
    typedef typename crab::cfg::statement_visitor<basic_block_label_t, N,
                                                  V>::arr_store_t arr_store_t;

    // assume all statements have the same type expression_t;
    using linear_expression_t = typename bin_op_t::linear_expression_t;
    using linear_constraint_t = typename assume_t::linear_constraint_t;

    array_segment_domain_t get_variables(const linear_expression_t &e) const {
      array_segment_domain_t res;
      for (auto const &v : e.variables())
        res += v;
      return res;
    }

    array_segment_domain_t get_variables(const linear_constraint_t &c) const {
      array_segment_domain_t res;
      for (auto const &v : c.variables())
        res += v;
      return res;
    }

  public:
    array_segment_domain_t _indexes;

    array_segment_visitor(array_segment_domain_t inv) : _indexes(inv) {}

    void visit(bin_op_t &s) {
      if (!(_indexes & s.lhs()).is_bottom()) {
        _indexes += get_variables(s.left());
        _indexes += get_variables(s.right());
      }
    }

    void visit(assign_t &s) {
      if (!(_indexes & s.lhs()).is_bottom()) {
        _indexes += get_variables(s.rhs());
      }
    }

    void visit(assume_t &s) {
      auto vars = get_variables(s.constraint());
      if (!(_indexes & vars).is_bottom())
        _indexes += vars;
    }

    void visit(arr_init_t &s) {
      _indexes += get_variables(s.lb_index());
      _indexes += get_variables(s.ub_index());
    }

    void visit(arr_load_t &s) { _indexes += get_variables(s.index()); }

    void visit(arr_store_t &s) {
      _indexes += get_variables(s.lb_index());
      _indexes += get_variables(s.ub_index());
    }

    void visit(unreach_t &) {}

    // FIXME: implement these
    void visit(select_t &) {}
    void visit(callsite_t &) {}
    void visit(havoc_t &) {}
    void visit(assert_t &) {}
  };

public:
  array_segment_ops(CFG cfg) : killgen_operations_api_t(cfg) {}

  virtual bool is_forward() override { return false; }

  virtual std::string name() override { return "array_segment"; }

  virtual void init_fixpoint() override {}

  virtual array_segment_domain_t entry() override {
    return array_segment_domain_t::bottom();
  }

  virtual array_segment_domain_t merge(array_segment_domain_t d1,
                                       array_segment_domain_t d2) override {
    return d1 | d2;
  }

  virtual array_segment_domain_t analyze(const basic_block_label_t &bb_id,
                                         array_segment_domain_t in) override {
    auto &bb = this->m_cfg.get_node(bb_id);
    array_segment_visitor vis(in);
    for (auto &s : bb) {
      s.accept(&vis);
    }
    return vis._indexes;
  }
};

template <class CFG>
class array_segmentation : public crab::iterators::killgen_fixpoint_iterator<
                               CFG, array_segment_ops<CFG>> {
public:
  using array_segment_ops_t = array_segment_ops<CFG>;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using statement_t = typename CFG::statement_t;
  using varname_t = typename CFG::varname_t;
  using array_segment_domain_t =
      typename array_segment_ops_t::array_segment_domain_t;

private:
  using killgen_fixpoint_iterator_t =
      crab::iterators::killgen_fixpoint_iterator<CFG, array_segment_ops_t>;
  using segment_map_t =
      boost::unordered_map<basic_block_label_t, array_segment_domain_t>;

  segment_map_t _segment_map;

public:
  array_segmentation(CFG cfg) : killgen_fixpoint_iterator_t(cfg) {}

  array_segmentation(const array_segmentation<CFG> &o) = delete;

  array_segmentation<CFG> &operator=(const array_segmentation<CFG> &o) = delete;

  void exec() {
    this->run();

    CRAB_LOG("array-segment",
             crab::outs()
                 << "Array segment variables alive at the block entries\n";);

    for (auto p :
         boost::make_iterator_range(this->in_begin(), this->in_end())) {
      CRAB_LOG("array-segment", crab::outs()
                                    << p.first << ":" << p.second << "\n";);
      _segment_map.insert(std::make_pair(p.first, p.second));
    }
    this->release_memory();
  }

  array_segment_domain_t get_variables(const basic_block_label_t &bb) const {
    auto it = _segment_map.find(bb);
    if (it == _segment_map.end())
      return array_segment_domain_t();
    else
      return it->second;
  }

  void write(crab_os &o) const {}
};

// Visitor for finding constants that might appear as array
// segment boundaries.
template <class CFG, class ArraySegmentDom>
class array_constant_segment_visitor
    : public crab::cfg::statement_visitor<typename CFG::basic_block_label_t,
                                          typename CFG::number_t,
                                          typename CFG::varname_t> {

  using BBL = typename CFG::basic_block_label_t;
  using N = typename CFG::number_t;
  using V = typename CFG::varname_t;
  using bin_op_t = typename crab::cfg::statement_visitor<BBL, N, V>::bin_op_t;
  using assign_t = typename crab::cfg::statement_visitor<BBL, N, V>::assign_t;
  using assume_t = typename crab::cfg::statement_visitor<BBL, N, V>::assume_t;
  using havoc_t = typename crab::cfg::statement_visitor<BBL, N, V>::havoc_t;
  using unreach_t = typename crab::cfg::statement_visitor<BBL, N, V>::unreach_t;
  using select_t = typename crab::cfg::statement_visitor<BBL, N, V>::select_t;
  using callsite_t =
      typename crab::cfg::statement_visitor<BBL, N, V>::callsite_t;
  using assert_t = typename crab::cfg::statement_visitor<BBL, N, V>::assert_t;

  using arr_init_t =
      typename crab::cfg::statement_visitor<BBL, N, V>::arr_init_t;
  using arr_load_t =
      typename crab::cfg::statement_visitor<BBL, N, V>::arr_load_t;
  using arr_store_t =
      typename crab::cfg::statement_visitor<BBL, N, V>::arr_store_t;

  using linear_expression_t = typename bin_op_t::linear_expression_t;
  using number_t = typename linear_expression_t::number_t;

public:
  using constant_set_t = std::vector<number_t>;

private:
  ArraySegmentDom _dom;
  constant_set_t _csts;

public:
  array_constant_segment_visitor(ArraySegmentDom dom) : _dom(dom) {}

  constant_set_t get_constants() { return _csts; }

  /// XXX: we focus for now only on assignments
  void visit(assign_t &s) {
    if (!(_dom & s.lhs()).is_bottom()) {
      if (s.rhs().is_constant() && s.rhs().constant() >= 0 &&
          std::find(_csts.begin(), _csts.end(), s.rhs().constant()) ==
              _csts.end()) {
        _csts.push_back(s.rhs().constant());
        // for decrementing loops we need to include the upper bound index i and
        // i+1
        _csts.push_back(s.rhs().constant() + 1);
      }
    }
  }

  void visit(bin_op_t &s) {}
  void visit(assume_t &s) {}
  void visit(arr_init_t &s) {}
  void visit(arr_load_t &s) {}
  void visit(arr_store_t &s) {}
  void visit(select_t &) {}
  void visit(callsite_t &) {}
  void visit(havoc_t &) {}
  void visit(assert_t &) {}
  void visit(unreach_t &) {}
};

} // end namespace analyzer
} // end namespace crab
