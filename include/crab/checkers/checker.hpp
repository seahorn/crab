#pragma once 

/* 
   A generic forward checker for properties
 */

#include <crab/common/types.hpp>
#include <crab/common/stats.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/inter_fwd_analyzer.hpp>

#include <memory>

namespace crab {

  namespace checker {

  template<typename Analyzer>
  class checker {
    /*
      The checker propagates the invariants that hold at the entry of
      each block (assume forward analysis) to each program point while
      checking for the properties. The analysis results are shared by
      all the property checkers so the analysis needs to be run only
      once.
     */

    checker(const checker<Analyzer>& other); // non construction copyable
    checker<Analyzer>& operator=(const checker<Analyzer>& other); // non-copyable
    
   public:

    typedef std::shared_ptr<property_checker<Analyzer>> prop_checker_ptr;
    typedef std::vector<prop_checker_ptr> prop_checker_vector;

   protected:

    prop_checker_vector m_checkers;
    
   public:

    checker(prop_checker_vector checkers): m_checkers(checkers) {}

    virtual ~checker() {}
    
    virtual void run() = 0;
    
    virtual void show(crab_os& o) {
      for (auto prop_checker: m_checkers) {
        prop_checker->write(o);
      }
    }

    // merge all the databases in one: useful for crab clients
    virtual checks_db get_all_checks() const {
      checks_db res;
      for (auto prop_checker: m_checkers) {
        res += prop_checker->get_db();
      }
      return res;
    }

  };

  template<typename Analyzer>  
  class intra_checker: public checker<Analyzer> {
   public:

    typedef checker<Analyzer> base_checker_t;
    using typename base_checker_t::prop_checker_ptr;
    using typename base_checker_t::prop_checker_vector;

   private:

    typedef typename Analyzer::cfg_t cfg_t;
    typedef typename cfg_t::statement_t statement_t;
    typedef typename Analyzer::abs_dom_t abs_dom_t;
    typedef typename Analyzer::abs_tr_t abs_tr_t;

    Analyzer& m_analyzer;

   public:

    intra_checker(Analyzer& analyzer, prop_checker_vector checkers)
      : base_checker_t(checkers)
      , m_analyzer(analyzer) { }
      
    
    virtual void run() override {
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Started property checker.\n";);      
      crab::ScopedCrabStats __st__("Checker");
      cfg_t cfg = m_analyzer.get_cfg();

      // In some cases, the analyzer might know that an assertion is
      // safe but it cannot be proven by propagating only
      // invariants. This is possible if the analyzer is based on a
      // forward/backward refinement loop.
      std::set<const statement_t*> safe_assertions;
      m_analyzer.get_safe_assertions(safe_assertions);
      
      for (auto &bb: cfg) {
        for (auto checker: this->m_checkers) {
	  if (checker->is_interesting(bb)) {
	    crab::ScopedCrabStats __st__("Checker." + checker->get_property_name());
	    abs_dom_t inv = m_analyzer[bb.label()];
	    std::shared_ptr<abs_tr_t> abs_tr = m_analyzer.get_abs_transformer(&inv);
	    // propagate forward the invariants from the block entry 
	    // while checking the property
	    checker->set(&*abs_tr, safe_assertions);
	    for (auto &stmt: bb) {
	      stmt.accept(&*checker);
	    }
	  }
        }
      }
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished property checker.\n";);
    }
  };


  template< typename Analyzer>
  class inter_checker: public checker<Analyzer> {
   public:

    typedef checker<Analyzer> base_checker_t;
    using typename base_checker_t::prop_checker_ptr;
    using typename base_checker_t::prop_checker_vector;

   private:

    typedef typename Analyzer::cg_t cg_t;
    typedef typename Analyzer::cfg_t cfg_t;
    typedef typename cfg_t::statement_t statement_t;    
    typedef typename Analyzer::abs_dom_t abs_dom_t;
    typedef typename Analyzer::abs_tr_t abs_tr_t;
    
    Analyzer& m_analyzer;
    
   public:

    inter_checker(Analyzer& analyzer, prop_checker_vector checkers)
        : base_checker_t(checkers), m_analyzer(analyzer) { }
    
    virtual void run() override {
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Started property checker.\n";);
      crab::ScopedCrabStats __st__("Checker");
      cg_t& cg = m_analyzer.get_call_graph(); 
      for (auto &v: boost::make_iterator_range(vertices(cg))) {
        cfg_t cfg = v.get_cfg();
	
        for (auto &bb: cfg) {
	  // the forward+backward analyzer only works for
	  // intra-procedural analysis.
	  std::set<const statement_t*> safe_assertions;
          for (auto checker: this->m_checkers) {
            crab::ScopedCrabStats __st__("Checker." + checker->get_property_name());
            abs_dom_t inv = m_analyzer.get_pre(cfg, bb.label());
	    std::shared_ptr<abs_tr_t> abs_tr = m_analyzer.get_abs_transformer(&inv);
            // propagate forward the invariants from the block entry 
            // while checking the property
            checker->set(&*abs_tr, safe_assertions);
            for (auto &stmt: bb) {
              stmt.accept(&*checker);
	    }
          }
        }
      }
      CRAB_VERBOSE_IF(1, get_msg_stream() << "Finished property checker.";);

      
    }
  };

  } // end namespace
} // end namespace
