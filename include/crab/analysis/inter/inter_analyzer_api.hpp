#pragma once

#include <crab/support/debug.hpp>
#include <vector>

namespace crab{
namespace analyzer {

template<class CFG, class Dom>  
class inter_analyzer_summary;

/*** Any inter-procedural analysis must implement this API ***/  
template<class CallGraph, class InvariantDom, class SummaryDom>
class inter_analyzer_api {
public:
  /* API type definitions */  
  using cg_t = CallGraph;
  using cg_node_t = typename CallGraph::node_t;  
  using cfg_t = typename CallGraph::cfg_t;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;  
  using basic_block_t = typename cfg_t::basic_block_t;
  using invariant_abs_dom_t = InvariantDom;
  using summary_abs_dom_t = SummaryDom;
  using summary_t = inter_analyzer_summary<cfg_t, summary_abs_dom_t>;
  
  /* API methods */
  virtual ~inter_analyzer_api() {}
  virtual void run(invariant_abs_dom_t init) = 0;
  virtual cg_t& get_call_graph() = 0;
  virtual invariant_abs_dom_t
  get_pre(const cfg_t &cfg, const basic_block_label_t &bb) const = 0;
  virtual invariant_abs_dom_t
  get_post(const cfg_t &cfg, const basic_block_label_t &bb) const = 0;
  virtual summary_t get_summary(const cfg_t &cfg) const = 0;
  virtual void clear() = 0;
};

  
/*** Definition of function summaries (aka contracts) ****/
/* A summary is a collection of pre- and post-condition pairs and a
 *  set of input and output variables.  Given the collection
 *  {(pre1,post1),...,(preN,postN)} and inputs I and outputs O then a
 *  summary is interpreted as
 *
 *   (pre1(I) => post1(I,O)) or ... or (preN(I) => postN(I,O))
 */   
template<class CFG, class AbsDom>  
class inter_analyzer_summary {
public:
  class pre_post_pair {
  public:
    // we intentionally pass them by value
    pre_post_pair(AbsDom pre, AbsDom post): m_pre(pre), m_post(post) {}
    const AbsDom& get_pre() const { return m_pre; }
    AbsDom& get_pre() { return m_pre; }
    const AbsDom& get_post() const { return m_post; }
    AbsDom& get_post() { return m_post; }
  private:
    AbsDom m_pre;
    AbsDom m_post;
  };

private:  
  using pre_post_pair_vector = std::vector<pre_post_pair>;

public:  
  using function_decl_t = typename CFG::fdecl_t;
  using variable_t = typename CFG::variable_t;
  using const_iterator = typename pre_post_pair_vector::const_iterator;
  using iterator = typename pre_post_pair_vector::iterator;
  
  inter_analyzer_summary(const function_decl_t &fdecl)
    : m_fdecl(fdecl) {
    // Initially empty summary
  }

  // Add a new pair pre- and post-conditions.
  // A __precondition__ is a condition (represented as an abstract
  // state) that must be satisfied so that the function can be
  // invoked.  A precondition should be expressed only on the
  // function's input parameters.
  //
  // A __postcondition__ states what is known after the function has
  // been invoked.  A postcondition is an abstract state that relates
  // input with output parameters.
  void add(const AbsDom &pre, const AbsDom &post) {
    m_pre_post_pairs.emplace_back(pre_post_pair(pre, post));
  }

  iterator begin() { return m_pre_post_pairs.begin();}
  iterator end() { return m_pre_post_pairs.end();}  
  const_iterator begin() const { return m_pre_post_pairs.begin();}
  const_iterator end() const { return m_pre_post_pairs.end();}  
  
  const function_decl_t &get_function_declaration() const
  { return m_fdecl; }

  void write(crab_os &o) const {
    o << m_fdecl << "\n";
    if (m_pre_post_pairs.empty()) {
      o << "empty\n";
    } else {
      for (unsigned i=0,sz=m_pre_post_pairs.size();i<sz;) {
	o << "Pre: " <<  m_pre_post_pairs[i].get_pre() << " -- Post: "
	  <<  m_pre_post_pairs[i].get_post();
	++i;
	if (i<sz) {
	o << "\n";
	}
      }
    }
  }

  friend crab_os &operator<<(crab_os&o,  const inter_analyzer_summary &sum) {
    sum.write(o);
    return o;
  }
  
private:
  const function_decl_t &m_fdecl;  
  std::vector<pre_post_pair> m_pre_post_pairs;
};

} // end namespace analyzer
} // end namespace crab
