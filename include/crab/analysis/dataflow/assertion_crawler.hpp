#pragma once

#include <crab/analysis/graphs/cdg.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cg/cg_bgl.hpp> // for sccg.hpp
#include <crab/domains/discrete_domains.hpp>
#include <crab/fixpoint/killgen_fixpoint_iterator.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/indexable.hpp>

#include <boost/range/iterator_range.hpp>
#include <unordered_map>
#include <unordered_set>

/**
 *  Dataflow analysis that for each block b it computes facts <i,V>
 *  meaning that there exists a path emanating from b that will reach
 *  assertion with id=i and its evaluation depends on the set of
 *  variables V.
 **/
namespace crab {
namespace analyzer {
  
using namespace crab::domains;
using namespace crab::cfg;

namespace assertion_crawler_impl{
// Convenient wrapper for using an assert statement as a key in
// patricia tries.
template <typename CFG>
class assert_wrapper : public indexable {
public:
  using statement_t = statement<typename CFG::basic_block_label_t,
				typename CFG::number_t,
				typename CFG::varname_t>;
private:  
  using this_type = assert_wrapper<CFG>;
  // unique identifier for the assert needed for being used as key
  ikos::index_t m_id;
  // the assert statement
  statement_t *m_assert;
public:
  assert_wrapper(ikos::index_t id, statement_t *assert)
    : m_id(id), m_assert(assert) {
    if (!m_assert) {
      CRAB_ERROR("assert_wrapper with null assert statement");
    }
    if (!(m_assert->is_assert() || m_assert->is_ref_assert() || m_assert->is_bool_assert())) {
      CRAB_ERROR("assert_wrapper must be assert, ref_assert, or bool_assert");
    }
  }
  
  statement_t &get() {
    return *m_assert;
  }
  const statement_t &get() const {
    return *m_assert;
  }
  virtual ikos::index_t index() const override {
    return m_id;
  }
  bool operator==(this_type o) const {
    return m_id == o.m_id;
  }
  bool operator<(this_type o) const {
    return m_id < o.m_id;
  }
  virtual void write(crab::crab_os &o) const override {
    o << "\"" << *m_assert << "\"";
  }
};

template<typename CFG>
class assertion_crawler_domain {
 public:
  using variable_t = typename CFG::variable_t;
  using number_t = typename CFG::number_t;
  using assert_wrapper_t = typename assertion_crawler_impl::assert_wrapper<CFG>;  
  using var_dom_t = ikos::discrete_domain<variable_t>;
  
  using first_domain_t = discrete_pair_domain<assert_wrapper_t, var_dom_t>;
  using second_domain_t = discrete_pair_domain<variable_t, var_dom_t>;
  using this_type = assertion_crawler_domain<CFG>;

 private:
  // the dataflow solution: map assertions to set of variables.
  first_domain_t m_first;
  // for inter-procedural analysis: map output variables to the set of
  // variables whose data may flow to.
  second_domain_t m_second;
  
 public:

  assertion_crawler_domain(const std::vector<variable_t> &vars)
    : m_first(first_domain_t::bottom()), m_second(second_domain_t::bottom()) {
    
    for (variable_t v: vars) {
      m_second.set(v, var_dom_t(v));
    }
  }
  
  assertion_crawler_domain()
    : m_first(first_domain_t::bottom()), m_second(second_domain_t::bottom()) {}

  void set_to_bottom() {
    m_first = first_domain_t::bottom();
    m_second = second_domain_t::bottom();
  }

  void set_to_top() {
    m_first = first_domain_t::top();
    m_second = second_domain_t::top();
  }
  
  static this_type top() {
    this_type res;
    res.get_first() = first_domain_t::top();
    res.get_second() = second_domain_t::top();
    return res;
  }

  static this_type bottom() {
    this_type res;
    res.get_first() = first_domain_t::bottom();
    res.get_second() = second_domain_t::bottom();
    return res;
  }

  first_domain_t &get_first() {
    return m_first;
  }

  const first_domain_t &get_first() const {
    return m_first;
  }
  
  second_domain_t &get_second() {
    return m_second;
  }
  
  const second_domain_t &get_second() const {
    return m_second;
  }

  bool operator<=(const this_type& other) const {
    return (get_first() <= other.get_first() && get_second() <= other.get_second());
  }
  
  this_type merge(const this_type& other) const {
    this_type res(*this);
    res.get_first() = res.get_first() | other.get_first();
    res.get_second() = res.get_second() | other.get_second();
    return res;
  }

  void write(crab::crab_os &o) const {
    o << "(";
    m_first.write(o);
    o << ",";
    m_second.write(o);
    o << ")";
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, const assertion_crawler_domain<CFG> &dom) {
    dom.write(o);
    return o;
  }
};

// Helper that applies a function to each pair of discrete_pair_domain.
template<typename Key, typename Value>
struct transform_discrete_pair_domain
  : public std::unary_function<discrete_pair_domain<Key,Value>,
			       discrete_pair_domain<Key,Value>> {
  using discrete_pair_domain_t = discrete_pair_domain<Key, Value>;
  
  template<typename Op>
  discrete_pair_domain_t apply_on_value(discrete_pair_domain_t sol, Op &op) {
    // All of this is needed because discrete_pair_domain_t cannot
    // be modified in-place
    if (sol.is_bottom() || sol.is_top()) {
      return sol;
    }
    for (auto kv : sol) {
      std::pair<Value,bool> p = op(kv.second);
      if (p.second) {
	sol.set(kv.first, p.first);
      }
    }
    return sol;
  }
  
  template<typename Op>      
  discrete_pair_domain_t apply_on_key_and_value(discrete_pair_domain_t sol, Op &op) {
    // All of this is needed because discrete_pair_domain_t cannot
    // be modified in-place
    if (sol.is_bottom() || sol.is_top()) {
      return sol;
    }
    for (auto kv : sol) {
      std::pair<Value,bool> p = op(kv.first, kv.second);
      if (p.second) {
	sol.set(kv.first, p.first);
      }
    }
    return sol;
  }  
};

// Wrapper to store a function summary that is represented by
// assertion_crawler_domain.
template<class CFG>  
class summary {
public:
  using variable_t = typename CFG::variable_t;
  using abs_dom_t = assertion_crawler_domain<CFG>;
  
private:  
  std::vector<variable_t> m_inputs;
  std::vector<variable_t> m_outputs;
  abs_dom_t m_sum;
public:
  summary(const std::vector<variable_t> &inputs,
	  const std::vector<variable_t> &outputs,
	  abs_dom_t sum)
    : m_sum(sum) {
    std::copy(inputs.begin(), inputs.end(), std::back_inserter(m_inputs));
    std::copy(outputs.begin(), outputs.end(), std::back_inserter(m_outputs));    
  }
  const std::vector<variable_t> &get_inputs() const {
    return m_inputs;
  }
  const std::vector<variable_t> &get_outputs() const {
    return m_outputs;
  }
  const abs_dom_t &get_summary() const {
    return m_sum;
  }
  abs_dom_t &get_summary() {
    return m_sum;
  }

  void write(crab_os &o) const {
    o << "Inputs:{";
    for(unsigned i=0,sz=m_inputs.size();i<sz;) {
      o << m_inputs[i];
      ++i;
      if (i<sz) {
	o << ",";
      }
    }
    o << "} Outputs:{";
    for(unsigned i=0,sz=m_outputs.size();i<sz;) {
      o << m_outputs[i];
      ++i;
      if (i<sz) {
	o << ",";
      }
    }
    o << "} " << m_sum;
  }
};
} // end namespace assertion_crawler_impl

template<class CFG>
class assertion_crawler;

// Define the operations of the dataflow analysis
template <class CFG>
class assertion_crawler_operations
  : public killgen_operations_api<CFG, assertion_crawler_impl::assertion_crawler_domain<CFG>> {

  friend class assertion_crawler<CFG>;
  
private:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using basic_block_t = typename CFG::basic_block_t;
  using statement_t = typename basic_block_t::statement_t;
  using variable_t = typename CFG::variable_t;    
  using varname_t = typename CFG::varname_t;  
  using number_t = typename CFG::number_t;
  using fdecl_t = typename CFG::fdecl_t;  
  using live_t = typename statement_t::live_t;
public:
  using assertion_crawler_domain_t = assertion_crawler_impl::assertion_crawler_domain<CFG>;
  using var_dom_t = typename assertion_crawler_domain_t::var_dom_t;  
  using assert_wrapper_t = assertion_crawler_impl::assert_wrapper<CFG>;  
  using assert_map_t = typename std::unordered_map<statement_t *, assert_wrapper_t>;
  // control-dependency graph: map a CFG block to the set of blocks
  // which control-dependent on it.
  using cdg_t = std::unordered_map<basic_block_label_t, std::vector<basic_block_label_t>>;
  // inter-procedural analysis
  using summary_map_t = callsite_or_fdecl_map<CFG, assertion_crawler_impl::summary<CFG>>;
private:  
  using killgen_operations_api_t = killgen_operations_api<CFG, assertion_crawler_domain_t>;
  using assert_map_domain_t =  typename assertion_crawler_domain_t::first_domain_t;
  using summary_dependencies_domain_t =  typename assertion_crawler_domain_t::second_domain_t;
public:
  ////============================================================////
  //            Propagate data/control dependencies
  ////============================================================////  
  class transfer_function : public statement_visitor<basic_block_label_t, number_t, varname_t> {
    using visitor_t = statement_visitor<basic_block_label_t, number_t, varname_t>; 
    using bin_op_t = typename visitor_t::bin_op_t;
    using assign_t = typename visitor_t::assign_t;
    using assume_t = typename visitor_t::assume_t;
    using select_t = typename visitor_t::select_t;
    using assert_t = typename visitor_t::assert_t;
    using int_cast_t = typename visitor_t::int_cast_t;
    using havoc_t = typename visitor_t::havoc_t;
    using unreach_t = typename visitor_t::unreach_t;
    using callsite_t = typename visitor_t::callsite_t;
    using intrinsic_t = typename visitor_t::intrinsic_t;
    using make_ref_t = typename visitor_t::make_ref_t;
    using remove_ref_t = typename visitor_t::remove_ref_t;
    using region_init_t = typename visitor_t::region_init_t;
    using region_copy_t = typename visitor_t::region_copy_t;
    using region_cast_t = typename visitor_t::region_cast_t;
    using load_from_ref_t = typename visitor_t::load_from_ref_t;
    using store_to_ref_t = typename visitor_t::store_to_ref_t;
    using gep_ref_t = typename visitor_t::gep_ref_t;
    using assume_ref_t = typename visitor_t::assume_ref_t;
    using assert_ref_t = typename visitor_t::assert_ref_t;
    using select_ref_t = typename visitor_t::select_ref_t;
    using bool_bin_op_t = typename visitor_t::bool_bin_op_t;
    using bool_assign_cst_t = typename visitor_t::bool_assign_cst_t;
    using bool_assign_var_t = typename visitor_t::bool_assign_var_t;
    using bool_assume_t = typename visitor_t::bool_assume_t;
    using bool_select_t = typename visitor_t::bool_select_t;
    using bool_assert_t = typename visitor_t::bool_assert_t;
    using variable_t = typename CFG::variable_t;

    /** Add data-dependencies **/
    class add_data_deps:
      public std::unary_function<var_dom_t,  std::pair<var_dom_t, bool>> {
						    
      var_dom_t m_uses;
      var_dom_t m_defs;
    public:
      add_data_deps(const live_t &l): m_uses(var_dom_t::bottom()), m_defs(var_dom_t::bottom()) {
        for (auto v : boost::make_iterator_range(l.uses_begin(), l.uses_end())) {
          m_uses += v;
	}
        for (auto v : boost::make_iterator_range(l.defs_begin(), l.defs_end())) {
          m_defs += v;
	}
      }
      
      add_data_deps(const var_dom_t &uses, const var_dom_t &defs) {
	for (auto v: uses) { m_uses += v;}
	for (auto v: defs) { m_defs += v;}	
      }

      add_data_deps(const add_data_deps &o) = delete;

      std::pair<var_dom_t, bool> operator()(var_dom_t d) {
        bool change = false;
        if (m_defs.is_bottom() && !m_uses.is_bottom() && !(m_uses & d).is_bottom()) {
          d += m_uses;
          change = true;
        }
        if (!(d & m_defs).is_bottom()) {
          d -= m_defs;
          d += m_uses;
          change = true;
        }
        return std::make_pair(d, change);
      }
    };

    /** Add control-dependencies **/
    class add_control_deps :
      public std::binary_function<assert_wrapper_t, var_dom_t, std::pair<var_dom_t, bool>> {
				  
      const cdg_t &cdg;
      const std::vector<basic_block_label_t> &roots;
      var_dom_t uses;

      // return true if we find a path in cdg from root to target
      // FIXME: do caching for the queries
      bool reach(const basic_block_label_t &root,
                 const basic_block_label_t &target,
                 std::unordered_set<basic_block_label_t> &visited) {
        if (root == target)
          return true;

        // break cycles
        if (visited.find(root) != visited.end())
          return false;

        visited.insert(root);
        auto it = cdg.find(root);
        if (it == cdg.end())
          return false;

        for (auto child : it->second) {
          if (reach(child, target, visited))
            return true;
        }
        return false;
      }

      bool reach(const basic_block_label_t &target) {
        std::unordered_set<basic_block_label_t> visited;
        for (auto r : roots)
          if (reach(r, target, visited))
            return true;
        return false;
      }

    public:
      add_control_deps(const cdg_t &_cdg,
                       const std::vector<basic_block_label_t> &_roots,
                       const live_t &l)
          : cdg(_cdg), roots(_roots), uses(var_dom_t::bottom()) {
        for (auto v : boost::make_iterator_range(l.uses_begin(), l.uses_end())) {
          uses += v;
	}
      }

      add_control_deps(const add_data_deps &o) = delete;

      std::pair<var_dom_t, bool> operator()(assert_wrapper_t w, var_dom_t d) {
        bool change = false;
        if (reach(w.get().get_parent()->label())) {
          d += uses;
          change = true;
        }
        return std::make_pair(d, change);
      }
    };

    /** Remove data-dependencies **/
    class remove_deps
      : public std::unary_function<var_dom_t, std::pair<var_dom_t, bool>> {
      
      var_dom_t vars;
    public:
      remove_deps(const variable_t &v) : vars(var_dom_t::bottom()) {
        vars += v;
      }

      remove_deps(const std::vector<variable_t> &vs)
          : vars(var_dom_t::bottom()) {
        for (auto v : vs) {
          vars += v;
        }
      }

      remove_deps(const remove_deps &o) = delete;

      std::pair<var_dom_t, bool> operator()(var_dom_t d) {
        bool change = false;
        if (!(d & vars).is_bottom()) {
          d -= vars;
          change = true;
        }
        return std::make_pair(d, change);
      }
    };

    using transform_sol_t =
	 assertion_crawler_impl::transform_discrete_pair_domain<assert_wrapper_t, var_dom_t>;
    using transform_sum_deps_t =
	 assertion_crawler_impl::transform_discrete_pair_domain<variable_t, var_dom_t>;
    
    assertion_crawler_domain_t m_sol;
    // map each assertion to a unique identifier
    assert_map_t &m_assert_map;
    // control-dependence graph (it can be empty)
    const cdg_t &m_cdg;
    // summaries from other functions
    summary_map_t &m_summaries;
    bool m_opt_ignore_region_offset;
    
    inline void process_assertion(statement_t &s) {
      assert(s.is_assert() || s.is_ref_assert() || s.is_bool_assert());
      
      auto it = m_assert_map.find(&s);
      if (it != m_assert_map.end())
        return;

      var_dom_t vdom = var_dom_t::bottom();
      auto const &l = s.get_live();
      for (auto v : boost::make_iterator_range(l.uses_begin(), l.uses_end())) {
        vdom += v;
      }

      unsigned id = m_assert_map.size();
      assert_wrapper_t val(id, &s);
      m_assert_map.insert(typename assert_map_t::value_type(&s, val));
      m_sol.get_first().set(val, vdom);
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "*** " << s << "\n"
                                             << "\tAdded " << vdom << "\n";);
    }
    
    inline void propagate_data(statement_t &s) {
      CRAB_LOG("assertion-crawler-step",
	       crab::outs() << "*** " << s << "\n" << "\tBEFORE: " << m_sol << "\n");
      
      if (m_opt_ignore_region_offset) {
	std::unique_ptr<add_data_deps> op = nullptr;

	if (s.is_ref_make()) {
	  auto make_ref = static_cast<make_ref_t*>(&s);
	  var_dom_t defs = var_dom_t::bottom();
	  var_dom_t uses = var_dom_t::bottom();
	  defs += make_ref->lhs();
	  uses += make_ref->region();
	  op.reset(new add_data_deps(uses, defs));
	} else if (s.is_ref_remove()) {
	  auto remove_ref = static_cast<remove_ref_t*>(&s);
	  var_dom_t defs = var_dom_t::bottom();
	  var_dom_t uses = var_dom_t::bottom();
	  uses += remove_ref->region();
	  op.reset(new add_data_deps(uses, defs));
	} else if (s.is_ref_load()) { 
	  auto load_ref = static_cast<load_from_ref_t*>(&s);
	  var_dom_t defs = var_dom_t::bottom();
	  var_dom_t uses = var_dom_t::bottom();
	  defs += load_ref->lhs();
	  uses += load_ref->region();
	  op.reset(new add_data_deps(uses, defs));
	} else if (s.is_ref_store()) {
	  auto store_ref = static_cast<store_to_ref_t*>(&s);
	  var_dom_t defs = var_dom_t::bottom();
	  var_dom_t uses = var_dom_t::bottom();
	  defs += store_ref->region();
	  uses += store_ref->region();
	  if (store_ref->val().is_variable()) {
	    uses += store_ref->val().get_variable();
	  }
	  op.reset(new add_data_deps(uses, defs));
	} else if (s.is_ref_gep()) {
	  auto gep_ref = static_cast<gep_ref_t*>(&s);
	  if (gep_ref->lhs() != gep_ref->rhs()) {
	    var_dom_t defs = var_dom_t::bottom();
	    var_dom_t uses = var_dom_t::bottom();
	    defs += gep_ref->lhs();
	    uses += gep_ref->rhs();
	    op.reset(new add_data_deps(uses, defs));	  
	  }
	} else {
	  op.reset(new add_data_deps(s.get_live()));	  	
	}
	if (op) {
	  transform_sol_t f1; 
	  m_sol.get_first() = std::move(f1.apply_on_value(m_sol.get_first(), *op));
	  transform_sum_deps_t f2;
	  m_sol.get_second() = std::move(f2.apply_on_value(m_sol.get_second(), *op));	
	}
      }  else {
	add_data_deps op(s.get_live());
	transform_sol_t f1;
	m_sol.get_first() = std::move(f1.apply_on_value(m_sol.get_first(), op));
	transform_sum_deps_t f2;
	m_sol.get_second() = std::move(f2.apply_on_value(m_sol.get_second(), op));
      }
      
      CRAB_LOG("assertion-crawler-step", crab::outs() << "\tAFTER " << m_sol << "\n";);
    }

    inline void propagate_data_and_control(statement_t &s) {
      CRAB_LOG("assertion-crawler-step", crab::outs()
	       << "*** " << s << "\n"
	       << "\tBEFORE: " << m_sol << "\n");
      // -- add data dependencies
      add_data_deps op(s.get_live());
      transform_sol_t df1;
      m_sol.get_first() = std::move(df1.apply_on_value(m_sol.get_first(), op));
      transform_sum_deps_t df2;
      m_sol.get_second() = std::move(df2.apply_on_value(m_sol.get_second(), op));      
      CRAB_LOG("assertion-crawler-step",
               crab::outs() << "\tAFTER data-dep " << m_sol << "\n";);

      if (!m_cdg.empty()) {
	// -- add control dependencies
	for (auto const &pred : boost::make_iterator_range(s.get_parent()->prev_blocks())) {
	  // it->second is the set of basic blocks that control
	  // dependent on s' block
	  auto it = m_cdg.find(pred);
	  if (it != m_cdg.end()) {
	    auto const &children = it->second;
	    CRAB_LOG("assertion-crawler-step-control", crab::outs() << "{";
		     for (auto &c : children) {
		       crab::outs() << c << ";";
		     }
		     crab::outs() << "} control-dependent on " << pred << "\n";);
	    add_control_deps op(m_cdg, children, s.get_live());
	    transform_sol_t cf;
	    m_sol.get_first() = std::move(cf.apply_on_key_and_value(m_sol.get_first(), op));
	    CRAB_LOG("assertion-crawler-step-control",
		     crab::outs() << "\tAFTER control-dep " << m_sol << "\n";);
	  }
	}
      }
    }
    
  public:
    transfer_function(assertion_crawler_domain_t init,
		      const cdg_t &g,
		      assert_map_t &assert_map,
		      summary_map_t &summaries,
		      bool ignore_region_offset)
      : m_sol(init), m_assert_map(assert_map), m_cdg(g), m_summaries(summaries),
	m_opt_ignore_region_offset(ignore_region_offset) {}

    assertion_crawler_domain_t get_solution() const {
      return m_sol;
    }

    virtual void visit(bin_op_t &s) override {
      propagate_data(s);
    }

    virtual void visit(assign_t &s) override {
      propagate_data(s);
    }

    virtual void visit(assume_t &s) override {
      propagate_data_and_control(s);
    }

    virtual void visit(select_t &s) override {
      propagate_data(s);
    }

    virtual void visit(assert_t &s) override {
      process_assertion(s);
    }

    virtual void visit(int_cast_t &s) override {
      propagate_data(s);
    }
    
    virtual void visit(unreach_t &) override {
      m_sol.set_to_bottom();
    }

    virtual void visit(havoc_t &s) override {
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "*** " << s << "\n"
                                             << "\tBEFORE: " << m_sol << "\n");
      remove_deps op(s.get_variable());
      transform_sol_t f1;
      m_sol.get_first() = std::move(f1.apply_on_value(m_sol.get_first(), op));
      transform_sum_deps_t f2;      
      m_sol.get_second() = std::move(f2.apply_on_value(m_sol.get_second(), op));      
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "\tAFTER " << m_sol << "\n";);
    }

    // For each key-value pair rename from variables with to variables
    // *only* on the value. This is an expensive operation because it
    // creates a new discrete_pair_domain from scratch.
    template<typename Key>
    static void rename(discrete_pair_domain<Key, var_dom_t> &dpd,
		       const std::vector<variable_t> &from, const std::vector<variable_t> &to) {
      using discrete_pair_domain_t = discrete_pair_domain<Key, var_dom_t>;
      if (from.empty() || dpd.is_top() || dpd.is_bottom()) {
	return;
      }
      if (from.size() != to.size()) {
	CRAB_ERROR("discrete_pair_domain::rename with input vectors of different sizes");
      }
      discrete_pair_domain_t res = discrete_pair_domain_t::bottom();
      for (auto kv: dpd) {
	auto key = kv.first;
	auto dom = kv.second;
	dom.rename(from, to);
	res.set(key, dom);
      }
      std::swap(dpd, res);
    }

    // dpd is of the form [key1 -> {o1,x1,...}, key2 -> {o2,x2,...}, ...]
    //     
    // sdd is a map that relates outputs with inputs.
    //     It is of the form [o1 -> {i1,i2}, o2  -> {i3,i4}]
    // This function replaces the occurrences of each sdd's key in dpd with its value.
    // That is, [key1 -> {i1,i2,x1,...}, key2 -> {i3,i4,x2,...}, ...]
    template<typename Key>
    static void apply_summary(discrete_pair_domain<Key, var_dom_t> &dpd,
			      summary_dependencies_domain_t &sdd) {
      auto outputs_to_inputs = [&sdd](const var_dom_t &var_dom) {
				 if (var_dom.is_top()) {
				   return var_dom_t::top();
				 } else if (var_dom.is_bottom()) {
				   return var_dom_t::bottom();
				 } else {
				   var_dom_t res = var_dom_t::bottom();
				   for (auto it=var_dom.begin(), et=var_dom.end(); it!=et; ++it) {
				     variable_t v = *it;
				     var_dom_t inputs(sdd[v]);
				     if (inputs.is_bottom()) {
				       // the callee doesn't know
				       // about v so we just propagate
				       // it as it's.
				       var_dom_t vs(v);
				       res += vs; // not in sdd
				     } else {
				       res += inputs;
				     }
				   }
				   return res;
				 }};

      using discrete_pair_domain_t = discrete_pair_domain<Key, var_dom_t>;      
      discrete_pair_domain_t out = discrete_pair_domain_t::bottom();
      for (auto kv: dpd) {
	var_dom_t res = outputs_to_inputs(kv.second);
	out.set(kv.first, res);
      }
      std::swap(dpd, out);
    }
				      
    
    virtual void visit(callsite_t &s) override {

      CRAB_LOG("assertion-crawler-step-cs", crab::outs()
	       << "*** " << s << "\n"
	       << "\tBEFORE: " << m_sol << "\n");
      
      auto it = m_summaries.find(&s);
      if (it != m_summaries.end()) {
	auto callee_summary_info = it->second;
	const std::vector<variable_t> &callsite_inputs = s.get_args();
	const std::vector<variable_t> &callsite_outputs = s.get_lhs();
	const std::vector<variable_t> &callee_inputs = callee_summary_info.get_inputs();
	const std::vector<variable_t> &callee_outputs = callee_summary_info.get_outputs();	
	auto &callee_summary = callee_summary_info.get_summary();

	CRAB_LOG("assertion-crawler-step-cs",
		 crab::outs() << "\tSummary at the callee: "
		 << callee_summary.get_second() << "\n";
		 crab::outs() << "\tCallee input variables: {";
		 for (auto const& v: callee_inputs) {
		   crab::outs() << v << ";";
		 }
		 crab::outs() << "}\n\tCallee output variables: {";
		 for (auto const& v: callee_outputs) {
		   crab::outs() << v << ";";
		 }
		 crab::outs() << "}\n";);

	// -- Update assertion map domain at the caller
	assert_map_domain_t amd(m_sol.get_first());
	rename(amd, callsite_outputs, callee_outputs);
	apply_summary(amd, callee_summary.get_second());
	rename(amd, callee_inputs, callsite_inputs);
	// Propagate the assertion map domain from the callee to the caller
	assert_map_domain_t callee_amd(callee_summary.get_first());
	rename(callee_amd, callee_inputs, callsite_inputs);
	m_sol.get_first() = amd | callee_amd;
	
	// -- Update the summary dependencies at the caller
	summary_dependencies_domain_t sdm(m_sol.get_second());
	rename(sdm, callsite_outputs, callee_outputs);
	apply_summary(sdm, callee_summary.get_second());
	rename(sdm, callee_inputs, callsite_inputs);
	m_sol.get_second() = sdm;
      } else {
	CRAB_LOG("assertion-crawler-step-cs",
		 CRAB_WARN("assertion-crawler did not find summary for callee at ", s));

	// crab::outs() << "Summary table:\n";
	// for (auto &kv: m_summaries) {
	//   crab::outs() << "\t";
	//   crab::outs() << &(kv.first) << " ";	  
	//   kv.first.write(crab::outs());
	//   crab::outs() << " -> ";
	//   crab::outs() << &(kv.second) << " ";	  
	//   kv.second.write(crab::outs());
	//   crab::outs() << "\n";
	// }
	propagate_data(s);
      }
      
      CRAB_LOG("assertion-crawler-step-cs", crab::outs()
	       << "\tAFTER " << m_sol << "\n";);
    }

    virtual void visit(intrinsic_t &s) override {
      propagate_data(s);      
    }

    virtual void visit(make_ref_t &s) override {
      propagate_data(s);
    }

    virtual void visit(remove_ref_t &s) override {
      propagate_data(s);
    }

    virtual void visit(region_init_t &s) override {
      propagate_data(s);
    }

    virtual void visit(region_copy_t &s) override {
      propagate_data(s);
    }

    virtual void visit(region_cast_t &s) override {
      propagate_data(s);
    }

    virtual void visit(load_from_ref_t &s) override {
      propagate_data(s);             
    }

    virtual void visit(store_to_ref_t &s) override {
      propagate_data(s);       
    }

    virtual void visit(gep_ref_t &s) override {
      propagate_data(s); 
    }

    virtual void visit(assume_ref_t &s) override {
      propagate_data_and_control(s);
    }

    virtual void visit(assert_ref_t &s) override {
      process_assertion(s);
    }
    
    virtual void visit(select_ref_t &s) override {
      propagate_data(s);
    }    

    virtual void visit(bool_bin_op_t &s) override {
      propagate_data(s);
    }

    virtual void visit(bool_assign_var_t &s) override {
      propagate_data(s);
    }

    virtual void visit(bool_assign_cst_t &s) override {
      propagate_data(s);
    }

    virtual void visit(bool_select_t &s) override {
      propagate_data(s);
    }
    
    virtual void visit(bool_assume_t &s) override {
      propagate_data_and_control(s);
    }
    
    virtual void visit(bool_assert_t &s) override {
      process_assertion(s);
    }
  };

private:
  assert_map_t &m_assert_map;
  summary_map_t &m_summaries;
  cdg_t m_cdg; // control dependencies

  // only data-dependencies (no control)  
  bool m_opt_only_data; 
  // ignore dependencies from references used to access regionsa  
  bool m_opt_ignore_region_offset; 
  
public:
  assertion_crawler_operations
  (CFG cfg, assert_map_t &assert_map, summary_map_t &summaries,
   bool only_data, bool ignore_region_offset)
    : killgen_operations_api_t(cfg),
      m_assert_map(assert_map),
      m_summaries(summaries),
      m_opt_only_data(only_data),
      m_opt_ignore_region_offset(ignore_region_offset) {}

  /* === Begin killgen_operations_api === */
  virtual bool is_forward() override {
    // This is a backward analysis
    return false;
  }

  virtual std::string name() override {
    return "assertion-crawler";
  }

  virtual void init_fixpoint() override {
    if (!m_opt_only_data) {
      crab::ScopedCrabStats __st__("Control-Dependency Graph");
      crab::analyzer::graph_algo::control_dep_graph(this->m_cfg, m_cdg);
    }
  }

  virtual assertion_crawler_domain_t entry() override {
    if (this->m_cfg.has_func_decl()) {
      return assertion_crawler_domain_t(this->m_cfg.get_func_decl().get_outputs());
    } else {
      return assertion_crawler_domain_t();
    }
  }

  virtual assertion_crawler_domain_t
  merge(assertion_crawler_domain_t d1, assertion_crawler_domain_t d2) override {
    return d1.merge(d2);
  }

  virtual assertion_crawler_domain_t analyze(const basic_block_label_t &bb_id,
					     assertion_crawler_domain_t out) override {
    auto &bb = this->m_cfg.get_node(bb_id);
    transfer_function vis(out, m_cdg, m_assert_map, m_summaries,
			  m_opt_ignore_region_offset);
    for (auto &s : boost::make_iterator_range(bb.rbegin(), bb.rend())) {
      s.accept(&vis);
    }
    return vis.get_solution();
  }
  /* ===  End killgen_operations_api === */
  
  const cdg_t &get_cdg() const {
    return m_cdg;
  }
  assert_map_t &get_assert_map() {
    return m_assert_map;
  }
  const assert_map_t &get_assert_map() const {
    return m_assert_map;
  }
  summary_map_t &get_summary_map() {
    return m_summaries;
  }
  const summary_map_t &get_summary_map() const {
    return m_summaries;
  }
};

template<typename CallGraph>
class inter_assertion_crawler;

/**
 * Intra-procedural assertion crawler dataflow analysis
 *
 * Compute for each basic block b a set of facts (i,V) such that
 * there exists a path from b that will check assertion i and its
 * evaluation depends on the set of variables V.
 **/
template <class CFG>
class assertion_crawler : public killgen_fixpoint_iterator<
                              CFG, assertion_crawler_operations<CFG>> {

  template<typename T> friend class inter_assertion_crawler;
public:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using varname_t = typename CFG::varname_t;

private:
  using assertion_crawler_op_t = assertion_crawler_operations<CFG>;
  using fixpo_t = killgen_fixpoint_iterator<CFG, assertion_crawler_op_t>; 
  using basic_block_t = typename CFG::basic_block_t;

public:
  using assertion_crawler_domain_t = typename assertion_crawler_op_t::assertion_crawler_domain_t;
  using assert_map_domain_t = typename assertion_crawler_domain_t::first_domain_t;
  using assert_map_t = typename assertion_crawler_op_t::assert_map_t;
  using summary_map_t = typename assertion_crawler_op_t::summary_map_t;    
  using results_map_t = std::unordered_map<basic_block_label_t, assertion_crawler_domain_t>;
private:

  results_map_t m_results;  
  assertion_crawler_op_t m_assert_crawler_op;

  // used by inter_assertion_crawler
  assertion_crawler_domain_t &lookup_results(const basic_block_label_t &bb) {
    auto it = m_results.find(bb);
    if (it == m_results.end()) {
      auto res = m_results.insert({bb, assertion_crawler_domain_t::top()});
      return res.first->second;
    } else {
      return it->second;
    }
  }
  
public:
  assertion_crawler(CFG cfg, assert_map_t &assert_map, summary_map_t &summaries,
		    bool only_data = false,
		    bool ignore_region_offset = false) 
    : fixpo_t(cfg, m_assert_crawler_op),
      m_assert_crawler_op(cfg, assert_map, summaries, only_data, ignore_region_offset) {
  }

  assertion_crawler(const assertion_crawler<CFG> &o) = delete;
  assertion_crawler<CFG> &operator=(const assertion_crawler<CFG> &o) = delete;

  void exec(void) {
    this->run();
    for (auto p : boost::make_iterator_range(this->in_begin(), this->in_end())) {
      m_results.insert(std::make_pair(p.first, p.second));
    }
    this->release_memory();
  }

  // return the dataflow facts that hold at the entry of the block bb
  assert_map_domain_t get_results(const basic_block_label_t &bb) const {
    auto it = m_results.find(bb);
    if (it == m_results.end()) {
      return assert_map_domain_t::top();
    } else {
      return it->second.get_first();
    }
  }

  // return the dataflow facts of the pre-state at each program point in bb
  void get_results(
      const basic_block_label_t &b,
      std::map<typename CFG::statement_t *, assert_map_domain_t> &res) const {
    auto it = m_results.find(b);
    if (it != m_results.end()) {
      if (!it->second.get_first().is_bottom()) {
        auto &bb = this->m_cfg.get_node(b);
        typename assertion_crawler_op_t::transfer_function vis
	  (it->second /* OUT dataflow facts */,
	   m_assert_crawler_op.m_cdg,
	   m_assert_crawler_op.m_assert_map,
	   m_assert_crawler_op.m_summaries,
	   m_assert_crawler_op.m_opt_ignore_region_offset);
        for (auto &s : boost::make_iterator_range(bb.rbegin(), bb.rend())) {
          s.accept(&vis);
          auto in = vis.get_solution().get_first();
          res.insert(std::make_pair(&s, in));
        }
      }
    }
  }

  void write(crab_os &o) const {
    o << "Assertion Crawler Analysis ";    
    if (this->m_cfg.has_func_decl()) {
      auto const &fdecl = this->m_cfg.get_func_decl();
      o << "for " << fdecl.get_func_name() <<"\n";
    } else {
      o << "\n";
    }
    
    // Print invariants in DFS to enforce a fixed order
    std::set<basic_block_label_t> visited;
    std::vector<basic_block_label_t> worklist;
    worklist.push_back(this->m_cfg.entry());
    visited.insert(this->m_cfg.entry());
    while (!worklist.empty()) {
      auto cur_label = worklist.back();
      worklist.pop_back();
      auto it = m_results.find(cur_label);
      assert(it != m_results.end());
      auto inv = it->second.get_first();
      crab::outs() << basic_block_traits<basic_block_t>::to_string(cur_label)
                   << "=" << inv << "\n";
      auto const &cur_node = this->m_cfg.get_node(cur_label);
      for (auto const &kid_label :
           boost::make_iterator_range(cur_node.next_blocks())) {
        if (visited.insert(kid_label).second) {
          worklist.push_back(kid_label);
        }
      }
    }
  }
};

template <typename CFG>
inline crab_os &operator<<(crab_os &o, const assertion_crawler<CFG> &ac) {
  ac.write(o);
  return o;
}


/** Context-insensitive interprocedural analysis **/
template <typename CallGraph>
class inter_assertion_crawler {
  using cg_node_t = typename CallGraph::node_t;
  using cg_edge_t = typename CallGraph::edge_t;
  using this_type = inter_assertion_crawler<CallGraph>;

public:
  using cg_t = CallGraph;
  using cg_ref_t = crab::cg::call_graph_ref<cg_t>;
  using cfg_t = typename cg_node_t::cfg_t;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;  
  using varname_t = typename cfg_t::varname_t;
  using number_t = typename cfg_t::number_t;
  using variable_t = typename cfg_t::variable_t;
  using assertion_crawler_t = assertion_crawler<cfg_t>;
  using assertion_crawler_domain_t = typename assertion_crawler_t::assertion_crawler_domain_t;   
  using assert_map_domain_t = typename assertion_crawler_t::assert_map_domain_t;
private:
  using summary_map_t = typename assertion_crawler_t::summary_map_t;
  using summary_t = typename summary_map_t::mapped_type;
  using assert_map_t = typename assertion_crawler_t::assert_map_t;      
  using inter_results_map_t = std::unordered_map<cfg_t, std::unique_ptr<assertion_crawler_t>>;
  
  cg_t &m_cg; // the call graph
  assert_map_t m_assert_map; // to assign id's to assert statements
  summary_map_t m_summaries; // summaries
  inter_results_map_t m_inter_results; // the results of the analysis
  bool m_opt_only_data; // whether consider only data dependencies
  // ignore dependencies from references used to access regions
  bool m_opt_ignore_region_offset;
public:
  inter_assertion_crawler(cg_t &cg, bool only_data = false, bool ignore_region_offset = false)
    : m_cg(cg), m_opt_only_data(only_data), m_opt_ignore_region_offset(ignore_region_offset) {
  }

  inter_assertion_crawler(const this_type &other) = delete;
  this_type &operator=(const this_type &other) = delete;
  
  void run(bool do_type_checking = true) {
    if (do_type_checking) {
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Type checking call graph ... ";);
      crab::CrabStats::resume("CallGraph type checking");
      m_cg.type_check();
      crab::CrabStats::stop("CallGraph type checking");
      CRAB_VERBOSE_IF(1, get_msg_stream() << "OK\n";);
    }
    
    CRAB_VERBOSE_IF(1, get_msg_stream()
                           << "Started inter-procedural assertion crawler analysis\n";);
    CRAB_LOG("assertion-crawler", m_cg.write(crab::outs()); crab::outs() << "\n");
    crab::ScopedCrabStats __st__("InterAssertionCrawler");

    auto analyze = [this](cfg_t cfg) {
		     assert(cfg.has_func_decl());
		     auto const &fdecl = cfg.get_func_decl();
		     const std::string &fun_name = fdecl.get_func_name();
		     CRAB_VERBOSE_IF(1, get_msg_stream() << "++ Analyzing " << fun_name << "...\n";);
		     // --- Run the intra-procedural analysis to compute summaries.
		     std::unique_ptr<assertion_crawler_t> intra_analysis
		       (new assertion_crawler_t(cfg, m_assert_map, m_summaries,
						m_opt_only_data, m_opt_ignore_region_offset));
		     intra_analysis->exec();
		     // --- Store the cfg's summary 
		     assertion_crawler_domain_t &results = intra_analysis->lookup_results(cfg.entry());
		     summary_t summary(fdecl.get_inputs(), fdecl.get_outputs(), results);
		     auto it = m_summaries.find(&fdecl);
		     if (it == m_summaries.end()) {
		       m_summaries.insert(std::make_pair(&fdecl, std::move(summary)));		       
		     } else {
		       it->second = std::move(summary);
		     }
		     
		     return std::move(intra_analysis);	   
		   };

    auto hasStabilized = [](assertion_crawler_t &old, assertion_crawler_t &current, cfg_t cfg) {
			     for (auto it = cfg.label_begin(), et=cfg.label_end(); it!=et; ++it) {
			       const basic_block_label_t &bb = *it;
			       const assertion_crawler_domain_t &old_d = old.lookup_results(bb);
			       const assertion_crawler_domain_t &current_d = current.lookup_results(bb);
			       if (!(current_d <= old_d)) {
				 return false;
			       }
			     }
			     return true;
			   };

    // If the fixpoint converges very slowly we give up after
    // max_fixpo_iters iterations
    const unsigned max_fixpo_iters = 10;
    
    std::vector<cg_node_t> rev_order;
    graph_algo::scc_graph<cg_ref_t> Scc_g(m_cg);
    graph_algo::rev_topo_sort(Scc_g, rev_order);
    for (auto const &n : rev_order) {
      std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members(n);

      // The SCC is recursive if it has more than one element or
      // there is only one that calls directly to itself.
      bool isRecursive =
	(scc_mems.size() > 1) ||
	std::any_of(m_cg.succs(n).first, m_cg.succs(n).second,
		    [n](const cg_edge_t &e) { return (n == e.dest()); });
      
      if (isRecursive) {
	// standard fixpoint among all the SCC members
	unsigned iter = 0;
	bool change = true;
	while (change || iter > max_fixpo_iters) {
	  change = false;
	  for (unsigned i=0, num_sccs=scc_mems.size();i<num_sccs;++i) {
	    cfg_t cfg = scc_mems[i].get_cfg();
	    std::unique_ptr<assertion_crawler_t> current = analyze(cfg);
	    if (iter == 0) {
	      m_inter_results.insert(std::make_pair(cfg, std::move(current)));
	      change = true;
	    } else {
	      auto it = m_inter_results.find(cfg);
	      assert(it != m_inter_results.end());
	      assertion_crawler_t &last = *(it->second);
	      change |= !hasStabilized(last, *current, cfg);
	      it->second = std::move(current);
	    }
	  }
	  ++iter;
	}
	if (iter > max_fixpo_iters) {
	  CRAB_WARN("Inter-procedural Assertion Crawler analysis was abruptly stopped on a ",
		    "recursive SCC after ", max_fixpo_iters, " iterations");
	}
      } else {
	for (auto m : scc_mems) {
	  cfg_t cfg = m.get_cfg();
	  std::unique_ptr<assertion_crawler_t> intra_analysis = analyze(cfg);
	  m_inter_results.insert(std::make_pair(cfg, std::move(intra_analysis)));
	}
      } 
    }

    CRAB_VERBOSE_IF(1, get_msg_stream()
		    << "Finished inter-procedural assertion crawler analysis\n";);
  }

  // Return the dataflow facts that hold at the entry of the basic block b in cfg.
  assert_map_domain_t get_results(const cfg_t &cfg, const basic_block_label_t &bb) const {
    auto it = m_inter_results.find(cfg);
    if (it != m_inter_results.end()) {
      return it->second->get_results(bb);
    } else {
      return assert_map_domain_t::top();
    }
  }

  void write(crab_os &o) const {
    for (auto &kv: m_inter_results) {
      kv.second->write(o);
    }
  }
  
};

} // namespace analyzer
} // namespace crab
