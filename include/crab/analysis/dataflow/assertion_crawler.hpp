#pragma once

#include <crab/analysis/graphs/cdg.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>
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

// Useful if the assertion crawler is used for warning analysis (diagnosis)  
// #define UNDERAPPROXIMATE_DATA_DEPS
  
using namespace crab::iterators;
using namespace crab::domains;
using namespace crab::cfg;

namespace assertion_crawler_details{
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
  
} // end namespace assertion_crawler_details

// forward declaration
template <typename CFG> class assertion_crawler;

// Define the operations of the dataflow analysis
template <class CFG>
class assertion_crawler_operations
    : public killgen_operations_api<
       CFG, discrete_pair_domain<assertion_crawler_details::assert_wrapper<CFG>,
				 ikos::discrete_domain<typename CFG::variable_t>>> {

private:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using basic_block_t = typename CFG::basic_block_t;
  using statement_t = typename basic_block_t::statement_t;  
  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;
  using fdecl_t = typename CFG::fdecl_t;  
  // set of uses and definitions of an instruction
  using live_t = crab::cfg::live<number_t, varname_t>;
public:
  using assert_wrapper_t = assertion_crawler_details::assert_wrapper<CFG>;
  using var_dom_t = ikos::discrete_domain<typename CFG::variable_t>;  
  using discrete_pair_domain_t = discrete_pair_domain<assert_wrapper_t, var_dom_t>;
  using assert_map_t = typename std::unordered_map<statement_t *, assert_wrapper_t>;
  // control-dependency graph: map a CFG block to the set of blocks
  // which control-dependent on it.
  using cdg_t = std::unordered_map<basic_block_label_t, std::vector<basic_block_label_t>>;
private:  
  using killgen_operations_api_t = killgen_operations_api<CFG, discrete_pair_domain_t>;  
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

    // Helper that applies a function to each pair of discrete_pair_domain.
    template <typename BinaryFunction>
    struct apply_discrete_pair
        : public std::unary_function<discrete_pair_domain_t, discrete_pair_domain_t> {
      using this_type = apply_discrete_pair<BinaryFunction>;
      using key_value_pair = std::pair<typename discrete_pair_domain_t::key_type,
				       typename discrete_pair_domain_t::value_type>;
      using function_type = std::binary_function<assert_wrapper_t, var_dom_t,
                                                 std::pair<var_dom_t, bool>>;
      static_assert(std::is_base_of<function_type, BinaryFunction>::value,
                    "Function must be subclass of type F");
      BinaryFunction &m_f;

    public:
      apply_discrete_pair(BinaryFunction &f) : m_f(f) {}
      apply_discrete_pair(const this_type &o) = delete;
      
      discrete_pair_domain_t operator()(discrete_pair_domain_t sol) {
	// All of this is needed because discrete_pair_domain_t cannot
	// be modified in-place
        if (sol.is_bottom() || sol.is_top()) {
          return sol;
	}
        for (auto kv : sol) {
          auto p = m_f(kv.first, kv.second);
          if (p.second) {
            sol.set(kv.first, p.first);
	  }
        }
        return sol;
      }
    };

    /** Add data-dependencies **/
    class add_data_deps: public std::binary_function<assert_wrapper_t, var_dom_t, 
						     std::pair<var_dom_t, bool>> {
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

      std::pair<var_dom_t, bool> operator()(assert_wrapper_t /*w*/, var_dom_t d) {
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
    class add_control_deps : public std::binary_function<assert_wrapper_t, var_dom_t,
							 std::pair<var_dom_t, bool>> {
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
        : public std::binary_function<assert_wrapper_t, var_dom_t,
                                      std::pair<var_dom_t, bool>> {
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
      
      std::pair<var_dom_t, bool> operator()(assert_wrapper_t /*w*/,
                                            var_dom_t d) {
        bool change = false;
        if (!(d & vars).is_bottom()) {
          d -= vars;
          change = true;
        }
        return std::make_pair(d, change);
      }
    };

    using apply_add_data_t = apply_discrete_pair<add_data_deps>;
    using apply_add_control_t = apply_discrete_pair<add_control_deps>;
    using apply_remove_t = apply_discrete_pair<remove_deps>;

    // solution: map blocks to pairs of assertion id and set of variables.
    discrete_pair_domain_t m_sol;
    // map each assertion to a unique identifier
    assert_map_t &m_assert_map;
    // control-dependence graph (it can be empty)
    const cdg_t &m_cdg;

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
      m_sol.set(val, vdom);
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "*** " << s << "\n"
                                             << "\tAdded " << vdom << "\n";);
      
    }
    
    inline void propagate_data(statement_t &s) {
      CRAB_LOG("assertion-crawler-step",
	       crab::outs() << "*** " << s << "\n" << "\tBEFORE: " << m_sol << "\n");

#ifdef UNDERAPPROXIMATE_DATA_DEPS
      std::unique_ptr<add_data_deps> op = nullptr;
      /// Under-approximation of region statements 
      if (s.is_ref_make()) {
	auto make_ref = static_cast<make_ref_t*>(&s);
	var_dom_t defs = var_dom_t::bottom();
	var_dom_t uses = var_dom_t::bottom();
	defs += make_ref->lhs();
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
	apply_add_data_t f(*op);
	m_sol = std::move(f(m_sol));
      }
#else
      add_data_deps op(s.get_live());
      apply_add_data_t f(op);
      m_sol = std::move(f(m_sol));
#endif      

      CRAB_LOG("assertion-crawler-step", crab::outs() << "\tAFTER " << m_sol << "\n";);
    }

    inline void propagate_data_and_control(statement_t &s) {
      CRAB_LOG("assertion-crawler-step", crab::outs()
	       << "*** " << s << "\n"
	       << "\tBEFORE: " << m_sol << "\n");
      // -- add data dependencies
      add_data_deps op(s.get_live());
      apply_add_data_t df(op);
      m_sol = std::move(df(m_sol));
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
	    apply_add_control_t cf(op);
	    m_sol = std::move(cf(m_sol));
	    CRAB_LOG("assertion-crawler-step-control",
		     crab::outs() << "\tAFTER control-dep " << m_sol << "\n";);
	  }
	}
      }
    }
    
  public:
    transfer_function(discrete_pair_domain_t init,
		      const cdg_t &g,
		      assert_map_t &assert_map)
      : m_sol(init), m_assert_map(assert_map), m_cdg(g) {}

    discrete_pair_domain_t get_solution() const {
      return m_sol;
    }

    void visit(bin_op_t &s) {
      propagate_data(s);
    }

    void visit(assign_t &s) {
      propagate_data(s);
    }

    void visit(assume_t &s) {
      propagate_data_and_control(s);
    }

    void visit(select_t &s) {
      propagate_data(s);
    }

    void visit(assert_t &s) {
      process_assertion(s);
    }

    void visit(int_cast_t &s) {
      propagate_data(s);
    }
    
    void visit(unreach_t &) {
      m_sol = discrete_pair_domain_t::bottom();
    }

    void visit(havoc_t &s) {
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "*** " << s << "\n"
                                             << "\tBEFORE: " << m_sol << "\n");
      remove_deps op(s.get_variable());
      apply_remove_t f(op);
      m_sol = std::move(f(m_sol));
      CRAB_LOG("assertion-crawler-step", crab::outs()
                                             << "\tAFTER " << m_sol << "\n";);
    }

    void visit(callsite_t &s) {
      // Intra-procedural
      
      CRAB_LOG("assertion-crawler-step", crab::outs()
	       << "*** " << s << "\n"
	       << "\tBEFORE: " << m_sol << "\n");
      
      propagate_data(s);
      
      CRAB_LOG("assertion-crawler-step", crab::outs()
	       << "\tAFTER " << m_sol << "\n";);
    }

    void visit(intrinsic_t &s) {
      propagate_data(s);      
    }

    void visit(make_ref_t &s) {
      propagate_data(s);
    }

    void visit(remove_ref_t &s) {
      propagate_data(s);
    }

    void visit(region_init_t &s) {
      propagate_data(s);
    }

    void visit(region_copy_t &s) {
      propagate_data(s);
    }

    void visit(region_cast_t &s) {
      propagate_data(s);
    }

    void visit(load_from_ref_t &s) {
      propagate_data(s);             
    }

    void visit(store_to_ref_t &s) {
      propagate_data(s);       
    }

    void visit(gep_ref_t &s) {
      propagate_data(s); 
    }

    void visit(assume_ref_t &s) {
      propagate_data_and_control(s);
    }

    void visit(assert_ref_t &s) {
      process_assertion(s);
    }
    
    void visit(select_ref_t &s) {
      propagate_data(s);
    }    

    void visit(bool_bin_op_t &s) {
      propagate_data(s);
    }

    void visit(bool_assign_var_t &s) {
      propagate_data(s);
    }

    void visit(bool_assign_cst_t &s) {
      propagate_data(s);
    }

    void visit(bool_select_t &s) {
      propagate_data(s);
    }
    
    void visit(bool_assume_t &s) {
      propagate_data_and_control(s);
    }
    
    void visit(bool_assert_t &s) {
      process_assertion(s);
    }
  };

private:
  assert_map_t &m_assert_map;
  cdg_t m_cdg; // control dependencies 
  bool m_only_data; // only data-dependencies (no control)
public:
  assertion_crawler_operations(CFG cfg, assert_map_t &assert_map, bool only_data)
    : killgen_operations_api_t(cfg),
      m_assert_map(assert_map), 
      m_only_data(only_data) {}

  /* === Begin killgen_operations_api === */
  virtual bool is_forward() override {
    // This is a backward analysis
    return false;
  }

  virtual std::string name() override {
    return "assertion-crawler";
  }

  virtual void init_fixpoint() override {
    if (!m_only_data) {
      crab::ScopedCrabStats __st__("Control-Dependency Graph");
      crab::analyzer::graph_algo::control_dep_graph(this->m_cfg, m_cdg);
    }
  }

  virtual discrete_pair_domain_t entry() override {
    return discrete_pair_domain_t::bottom();
  }

  virtual discrete_pair_domain_t
  merge(discrete_pair_domain_t d1, discrete_pair_domain_t d2) override {
    return d1 | d2;
  }

  virtual discrete_pair_domain_t analyze(const basic_block_label_t &bb_id,
					 discrete_pair_domain_t out) override {
    auto &bb = this->m_cfg.get_node(bb_id);
    transfer_function vis(out, m_cdg, m_assert_map);
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
};

/**
 * Intra-procedural assertion crawler dataflow analysis
 *
 * Compute for each basic block b a set of facts (i,V) such that
 * there exists a path from b that will check assertion i and its
 * evaluation depends on the set of variables V.
 **/
template <class CFG>
class assertion_crawler : public crab::iterators::killgen_fixpoint_iterator<
                              CFG, assertion_crawler_operations<CFG>> {

public:
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using varname_t = typename CFG::varname_t;

private:
  using assertion_crawler_op_t = assertion_crawler_operations<CFG>;
  using fixpo_t = crab::iterators::killgen_fixpoint_iterator<CFG, assertion_crawler_op_t>; 
  using basic_block_t = typename CFG::basic_block_t;

public:
  using discrete_pair_domain_t = typename assertion_crawler_op_t::discrete_pair_domain_t;
  using assert_map_t = typename assertion_crawler_op_t::assert_map_t;  
  using results_map_t = std::unordered_map<basic_block_label_t, discrete_pair_domain_t>;
private:

  results_map_t m_results;  
  assertion_crawler_op_t m_assert_crawler_op;

public:
  assertion_crawler(CFG cfg, assert_map_t &assert_map, bool only_data = false) 
    : fixpo_t(cfg, m_assert_crawler_op),
      m_assert_crawler_op(cfg, assert_map, only_data) {}

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
  discrete_pair_domain_t get_results(const basic_block_label_t &bb) {
    auto it = m_results.find(bb);
    if (it == m_results.end()) {
      return discrete_pair_domain_t::top();
    } else {
      return it->second;
    }
  }

  // return the dataflow facts of the pre-state at each program point in bb
  void get_results(
      const basic_block_label_t &b,
      std::map<typename CFG::statement_t *, discrete_pair_domain_t> &res) {
    auto it = m_results.find(b);
    if (it != m_results.end()) {
      if (!it->second.is_bottom()) {
        auto &bb = this->m_cfg.get_node(b);
        typename assertion_crawler_op_t::transfer_function vis
	  (it->second /* OUT dataflow facts */,
	   m_assert_crawler_op.get_cdg(),
	   m_assert_crawler_op.get_assert_map());
        for (auto &s : boost::make_iterator_range(bb.rbegin(), bb.rend())) {
          s.accept(&vis);
          auto in = vis.get_solution();
          res.insert(std::make_pair(&s, in));
        }
      }
    }
  }

  void write(crab_os &o) const {
    o << "Assertion Crawler Analysis\n";    
    if (this->m_cfg.has_func_decl()) {
      auto const &fdecl = this->m_cfg.get_func_decl();
      o << "for " << fdecl.get_func_name();
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
      auto inv = it->second;
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
  
} // namespace analyzer
} // namespace crab
