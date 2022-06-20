#pragma once

/**
   Dataflow analysis that computes the set of unjustified assumptions
   taken by Crab while proving the assertions.
 **/

#include <crab/analysis/dataflow/assertion_crawler.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace crab {

namespace analyzer {

// A class to represent an unjustified assumption
template <class CFG> class assumption {

public:
  using var_dom_t = ikos::discrete_domain<typename CFG::variable_t>;
  using statement_t = typename CFG::basic_block_t::statement_t;

protected:
  // unique id
  unsigned _id;
  // statement where the unsound assumption was assumed.
  const statement_t *_s;
  // a subset of s' variables on which the unsound assumption was assumed.
  var_dom_t _vars;

public:
  assumption(unsigned id, const statement_t *s, var_dom_t vars)
      : _id(id), _s(s), _vars(vars) {}

  virtual ~assumption() {}

  var_dom_t get_vars() const { return _vars; }

  const statement_t *get_statement() const { return _s; }

  std::string get_id_str() const {
    std::string res("a");
    return res + std::to_string(_id);
  }

  virtual void write(crab::crab_os &o) const = 0;
};

template <typename CFG>
inline crab_os &operator<<(crab::crab_os &o, const assumption<CFG> &a) {
  a.write(o);
  return o;
}

// Subclass of assumption for arithmetic overflow
template <class CFG> class overflow_assumption : public assumption<CFG> {

public:
  using base_type = assumption<CFG>;
  using this_type = overflow_assumption<CFG>;

  using var_dom_t = typename base_type::var_dom_t;
  using statement_t = typename base_type::statement_t;

private:
  var_dom_t _get_vars(const statement_t *s) const {
    var_dom_t vars;
    auto lvars = s->get_live();
    for (auto v :
         boost::make_iterator_range(lvars.uses_begin(), lvars.uses_end())) {
      vars += v;
    }
    return vars;
  }

public:
  overflow_assumption(unsigned id, const statement_t *s)
      : base_type(id, s, _get_vars(s)) {}

  virtual void write(crab::crab_os &o) const override {
    o << "assume NotOverflow(" << this->_vars << ")"
      << " as " << this->get_id_str();
  }
};

/** The analysis of unjustified assumptions.
 *
 * Compute for each assertion the set of unjustified assumptions,
 * that is, the unproved assumptions taken by the analyzer while
 * proving assertions.
 *
 * For now, we only take into consideration integer overflows.
 **/
template <class CFG> class assumption_analysis {

public:
  struct stats_t {
    unsigned num_assertions;
    unsigned max_assumptions_per_assertion;
    float avg_assumptions_per_assertion;

    stats_t(unsigned n1, unsigned n2, float n3)
        : num_assertions(n1), max_assumptions_per_assertion(n2),
          avg_assumptions_per_assertion(n3) {}

    void write(crab::crab_os &o) const {
      o << "BRUNCH_STAT Assertions " << num_assertions << "\n";
      o << "BRUNCH_STAT Assumptions "
        << avg_assumptions_per_assertion * num_assertions << "\n";
      o << "BRUNCH_STAT MaxAssumptionsPerAssertion "
        << max_assumptions_per_assertion << "\n";
      o << "BRUNCH_STAT AvgAssumptionsPerAssertion "
        << avg_assumptions_per_assertion << "\n";
    }

    friend crab_os &
    operator<<(crab::crab_os &o,
               const typename assumption_analysis<CFG>::stats_t &stats) {
      stats.write(o);
      return o;
    }
  };

public:
  using statement_t = typename CFG::statement_t;
  using assumption_t = assumption<CFG>;
  using assumption_ptr = std::shared_ptr<assumption_t>;
  using vector_assumption_ptr = std::vector<assumption_ptr>;

protected:
  // map an assert statement (assert, ref_assert, and bool_assert) to its unjustified assumptions
  using assumption_map_t =
      std::unordered_map<const statement_t *, vector_assumption_ptr>;
  // map an arbitrary statement to the assumptions originated in that statement
  using assumption_origin_map_t =
      std::unordered_map<const statement_t *, vector_assumption_ptr>;

  /** the cfg **/
  CFG m_cfg;
  /** map each assertion to its unjustified assumptions **/
  assumption_map_t m_assumption_map;
  /** map a statement to the set of assumptions originated by that statement **/
  assumption_origin_map_t m_origin_map;
  /** generate id's for assumptions **/
  unsigned m_id;

  static bool can_overflow(statement_t &s) {
    // TODO: cover select/assume with conditions that can overflow

    if (s.is_int_cast()) {
      auto cast_stmt =
          static_cast<const typename CFG::basic_block_t::int_cast_t *>(&s);
      return (cast_stmt->op() == crab::cfg::CAST_TRUNC);
    } else if (s.is_bin_op()) {
      auto bin_op_s =
          static_cast<const typename CFG::basic_block_t::bin_op_t *>(&s);
      return (bin_op_s->op() == crab::cfg::BINOP_ADD ||
              bin_op_s->op() == crab::cfg::BINOP_SUB ||
              bin_op_s->op() == crab::cfg::BINOP_MUL);
    } else {
      return false;
    }
  }

protected:
  void insert_assumption(const statement_t *a, assumption_ptr as) {
    auto it = m_assumption_map.find(a);
    if (it != m_assumption_map.end()) {
      it->second.push_back(as);
    } else {
      vector_assumption_ptr assumes = {as};
      m_assumption_map.insert(std::make_pair(a, assumes));
    }
  }

  // create a new overflow assumption
  assumption_ptr new_assumption(const statement_t *s) {
    assumption_ptr assume =
        std::make_shared<overflow_assumption<CFG>>(m_id++, s);
    auto it = m_origin_map.find(s);
    if (it != m_origin_map.end()) {
      it->second.push_back(assume);
    } else {
      vector_assumption_ptr assumes = {assume};
      m_origin_map.insert(std::make_pair(s, assumes));
    }
    return assume;
  }

public:
  assumption_analysis(CFG cfg) : m_cfg(cfg), m_id(0) {}

  virtual ~assumption_analysis() {}

  virtual void exec() = 0;

  // return the assumptions originated in s
  void get_originated_assumptions(const statement_t *s,
                                  std::vector<assumption_ptr> &out) {
    auto it = m_origin_map.find(s);
    if (it != m_origin_map.end()) {
      std::copy(it->second.begin(), it->second.end(), std::back_inserter(out));
    }
  }

  // return the unjustified assumptions for assertion s
  void get_assumptions(const statement_t *s, std::vector<assumption_ptr> &out) {
    auto it = m_assumption_map.find(s);
    if (it != m_assumption_map.end()) {
      std::copy(it->second.begin(), it->second.end(), std::back_inserter(out));
    }
  }

  stats_t get_stats() const {
    unsigned num_assertions = 0;
    unsigned max_assumptions_per_assertion = 0;
    float avg_assumptions_per_assertion = 0.0;
    unsigned num_assumptions = 0;

    for (auto &kv : m_assumption_map) {
      num_assertions++;
      num_assumptions += kv.second.size();
      if (kv.second.size() > max_assumptions_per_assertion)
        max_assumptions_per_assertion = kv.second.size();
    }
    if (num_assertions > 0)
      avg_assumptions_per_assertion =
          (float)num_assumptions / (float)num_assertions;

    return stats_t(num_assertions, max_assumptions_per_assertion,
                   avg_assumptions_per_assertion);
  }

  void write(crab::crab_os &o) const {
    // We order the assertions so that we always print them in the same order.
    std::vector<const statement_t *> assertions;
    assertions.reserve(m_assumption_map.size());
    for (auto &kv : m_assumption_map) {
      assertions.push_back(kv.first);
    }
    std::sort(assertions.begin(), assertions.end(),
              [](const statement_t *a1, const statement_t *a2) {
                crab_string_os oss1, oss2;
                oss1 << *a1;
                oss2 << *a2;
                return oss1.str() < oss2.str();
              });

    o << "******* UNJUSTIFIED ASSUMPTION ANALYSIS ********\n";
    for (const statement_t *assertion : assertions) {
      o << *assertion << " with unjustified assumptions:\n";
      auto it = m_assumption_map.find(assertion);
      assert(it != m_assumption_map.end());
      auto assumptions = it->second;
      for (assumption_ptr a : assumptions) {
        o << "\t" << *a << " in " << *(a->get_statement()) << "\n";
      }
    }
    // for (auto &kv: m_assumption_map) {
    //   o << *(kv.first) << " with unjustified assumptions:\n";
    //   for (auto a: kv.second) {
    //     o << "\t" << *a << " in " << *(a->get_statement()) << "\n";
    //   }
    // }

    o << "******************* STATS **********************\n";
    auto stats = get_stats();
    o << stats << "\n";
  }
};

template <typename CFG>
inline crab_os &operator<<(crab::crab_os &o,
                           const assumption_analysis<CFG> &aa) {
  aa.write(o);
  return o;
}

/**
    Based on "Guiding dynamic symbolic execution toward unverified program
    executions" by Christakis, Muller, and Wustholz (ICSE'16).

    For each assertion it considers all the assumptions
    backward-reachable in the CFG.
**/
template <class CFG>
class assumption_naive_analysis : public assumption_analysis<CFG> {

  using assumption_analysis_t = assumption_analysis<CFG>;
  using typename assumption_analysis_t::statement_t;
  using typename assumption_analysis_t::assumption_ptr;
  using typename assumption_analysis_t::vector_assumption_ptr;
  using bb_label_t = typename CFG::basic_block_label_t;
  using label_set_t = std::unordered_set<bb_label_t>;
  using bb_t = typename CFG::basic_block_t;

  void backward_reachable(const bb_label_t &r, const statement_t *a,
                          label_set_t &visited) {
    auto ret = visited.insert(r);
    if (!ret.second)
      return; // already in visited

    bb_t &bb = this->m_cfg.get_node(r);
    for (auto &s : boost::make_iterator_range(bb.begin(), bb.end())) {
      if (&s == a) {
        break;
      }
      
      if (assumption_analysis_t::can_overflow(s)) {
        assumption_ptr assume = this->new_assumption(&s);
        assumption_analysis_t::insert_assumption(a, assume);
      }
    }

    for (auto child : this->m_cfg.prev_nodes(r)) {
      backward_reachable(child, a, visited);
    }
  }

public:
  assumption_naive_analysis(CFG cfg) : assumption_analysis_t(cfg) {}

  virtual ~assumption_naive_analysis() {}

  virtual void exec() {
    for (auto &bb :
         boost::make_iterator_range(this->m_cfg.begin(), this->m_cfg.end())) {
      for (auto &s : boost::make_iterator_range(bb.begin(), bb.end())) {
        if (s.is_assert()) {
          label_set_t visited;
          backward_reachable(bb.label(), &s, visited);
        } else if (s.is_ref_assert()) {
	  // TODO
	} else if (s.is_bool_assert()) {
	  // TODO
	}
      }
    }
  }
};

/** Uses the assertion crawler analysis to improve precision **/
template <class CFG>
class assumption_dataflow_analysis : public assumption_analysis<CFG> {

  using assumption_analysis_t = assumption_analysis<CFG>;
  using typename assumption_analysis_t::assumption_ptr;
  using typename assumption_analysis_t::statement_t;
  using assertion_crawler_t = crab::analyzer::assertion_crawler<CFG>;
  using assert_map_domain_t = typename assertion_crawler_t::assert_map_domain_t;
  using typename assumption_analysis_t::vector_assumption_ptr;
  using pp_inv_map_t = std::map<statement_t *, assert_map_domain_t>;
  using basic_block_t = typename CFG::basic_block_t;

public:
  assumption_dataflow_analysis(CFG cfg) : assumption_analysis_t(cfg) {}

  virtual ~assumption_dataflow_analysis() override {}

  virtual void exec() override {
    typename assertion_crawler_t::assert_map_t assert_map;
    typename assertion_crawler_t::summary_map_t summaries;
    assertion_crawler_t assert_crawler(this->m_cfg, assert_map, summaries);
    assert_crawler.exec();

    for (auto &bb :
         boost::make_iterator_range(this->m_cfg.begin(), this->m_cfg.end())) {

      // find out if the block is of interest
      bool op_can_overflow = false;
      for (auto &s : bb) {
        if (assumption_analysis_t::can_overflow(s)) {
          op_can_overflow = true;
          break;
        }
      }

      if (!op_can_overflow)
        continue;

      auto bb_out = assert_crawler.get_results(bb.label());

      if (bb_out.is_bottom())
        continue;

      CRAB_LOG("aa",
               crab::outs()
                   << "Basic block    : "
                   << basic_block_traits<basic_block_t>::to_string(bb.label())
                   << "\n";
               crab::outs() << "OUT invariants : " << bb_out << "\n";);

      // get dataflow information at the pre-state at each program point
      pp_inv_map_t pp_in_map;
      assert_crawler.get_results(bb.label(), pp_in_map);

      for (auto &s : boost::make_iterator_range(bb.begin(), bb.end())) {
        if (!assumption_analysis_t::can_overflow(s))
          continue;

        auto in = pp_in_map[&s];
        CRAB_LOG("aa", crab::outs() << "Check s = " << s << "\n";
                 crab::outs() << "IN(s)   = " << in << "\n";);

        // -- Traverse all pairs (i,V) where i is the id for the
        //    assertion and V is the variables computed by the
        //    dataflow analysis. If V intersect with the variables
        //    on the rhs of the arithmetic operation then we
        //    assume an assumption to the assertion.
        assumption_ptr assume = this->new_assumption(&s);
        for (auto kv : boost::make_iterator_range(in.begin(), in.end())) {
          auto vars = kv.second;
          if (!(vars & assume->get_vars()).is_bottom())
            assumption_analysis_t::insert_assumption(&(kv.first.get()), assume);
        }
      }
    }
  }
};

} // namespace analyzer
} // namespace crab
