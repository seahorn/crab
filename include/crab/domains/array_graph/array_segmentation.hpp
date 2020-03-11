#pragma once

/*
   Quick analysis to collect the set of variables and constants that might
   define array segment boundaries.
*/

#include <crab/cfg/cfg.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/killgen_domain.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>

#include <boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>

namespace crab {

namespace analyzer {

template <class CFG>
class array_segment_ops
    : public crab::iterators::killgen_operations_api<
          CFG, domains::flat_killgen_domain<typename CFG::variable_t>> {

  typedef crab::iterators::killgen_operations_api<
      CFG, domains::flat_killgen_domain<typename CFG::variable_t>>
      killgen_operations_api_t;

public:
  typedef typename CFG::basic_block_label_t basic_block_label_t;
  typedef typename CFG::number_t N;
  typedef typename CFG::varname_t V;
  typedef typename CFG::variable_t variable_t;
  typedef typename killgen_operations_api_t::killgen_domain_t
      array_segment_domain_t;

private:
  class array_segment_visitor : public crab::cfg::statement_visitor<N, V> {
    typedef typename crab::cfg::statement_visitor<N, V>::bin_op_t bin_op_t;
    typedef typename crab::cfg::statement_visitor<N, V>::assign_t assign_t;
    typedef typename crab::cfg::statement_visitor<N, V>::assume_t assume_t;
    typedef typename crab::cfg::statement_visitor<N, V>::select_t select_t;
    typedef typename crab::cfg::statement_visitor<N, V>::assert_t assert_t;

    typedef typename crab::cfg::statement_visitor<N, V>::havoc_t havoc_t;
    typedef typename crab::cfg::statement_visitor<N, V>::unreach_t unreach_t;
    typedef typename crab::cfg::statement_visitor<N, V>::callsite_t callsite_t;

    typedef typename crab::cfg::statement_visitor<N, V>::arr_init_t arr_init_t;
    typedef typename crab::cfg::statement_visitor<N, V>::arr_load_t arr_load_t;
    typedef
        typename crab::cfg::statement_visitor<N, V>::arr_store_t arr_store_t;

    // assume all statements have the same type expression_t;
    typedef typename bin_op_t::linear_expression_t linear_expression_t;
    typedef typename assume_t::linear_constraint_t linear_constraint_t;

    array_segment_domain_t get_variables(linear_expression_t e) {
      array_segment_domain_t res;
      for (auto v : e.variables())
        res += v;
      return res;
    }

    array_segment_domain_t get_variables(linear_constraint_t c) {
      array_segment_domain_t res;
      for (auto v : c.variables())
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

  virtual array_segment_domain_t analyze(basic_block_label_t bb_id,
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
  typedef array_segment_ops<CFG> array_segment_ops_t;
  typedef typename CFG::basic_block_label_t basic_block_label_t;
  typedef typename CFG::statement_t statement_t;
  typedef typename CFG::varname_t varname_t;
  typedef typename array_segment_ops_t::array_segment_domain_t
      array_segment_domain_t;

private:
  typedef crab::iterators::killgen_fixpoint_iterator<CFG, array_segment_ops_t>
      killgen_fixpoint_iterator_t;
  typedef boost::unordered_map<basic_block_label_t, array_segment_domain_t>
      segment_map_t;

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

  array_segment_domain_t get_variables(basic_block_label_t bb) const {
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
    : public crab::cfg::statement_visitor<typename CFG::number_t,
                                          typename CFG::varname_t> {

  typedef typename CFG::number_t N;
  typedef typename CFG::varname_t V;
  typedef typename crab::cfg::statement_visitor<N, V>::bin_op_t bin_op_t;
  typedef typename crab::cfg::statement_visitor<N, V>::assign_t assign_t;
  typedef typename crab::cfg::statement_visitor<N, V>::assume_t assume_t;
  typedef typename crab::cfg::statement_visitor<N, V>::havoc_t havoc_t;
  typedef typename crab::cfg::statement_visitor<N, V>::unreach_t unreach_t;
  typedef typename crab::cfg::statement_visitor<N, V>::select_t select_t;
  typedef typename crab::cfg::statement_visitor<N, V>::callsite_t callsite_t;
  typedef typename crab::cfg::statement_visitor<N, V>::assert_t assert_t;

  typedef typename crab::cfg::statement_visitor<N, V>::arr_init_t arr_init_t;
  typedef typename crab::cfg::statement_visitor<N, V>::arr_load_t arr_load_t;
  typedef typename crab::cfg::statement_visitor<N, V>::arr_store_t arr_store_t;

  typedef typename bin_op_t::linear_expression_t linear_expression_t;
  typedef typename linear_expression_t::number_t number_t;

public:
  typedef std::vector<number_t> constant_set_t;

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
