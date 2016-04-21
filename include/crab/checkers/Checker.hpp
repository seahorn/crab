#ifndef CRAB_CHECKER_HPP
#define CRAB_CHECKER_HPP

/* 
   A generic forward checker for properties
 */

#include <crab/common/stats.hpp>
#include <crab/checkers/BaseProperty.hpp>
#include <crab/analysis/FwdAnalyzer.hpp>
#include <crab/analysis/InterFwdAnalyzer.hpp>

namespace crab {

  namespace checker {

  template<typename Analyzer>
  class Checker {
    /*
      The checker propagates the invariants that hold at the entry of
      each block (assume forward analysis) to each program point while
      checking for the properties. The analysis results are shared by
      all the property checkers so the analysis needs to be run only
      once.
     */
   public:

    typedef boost::shared_ptr<PropertyChecker<Analyzer> > prop_checker_ptr;
    typedef std::vector<prop_checker_ptr> prop_checker_vector;

   protected:

    prop_checker_vector m_checkers;
    
   public:

    Checker (prop_checker_vector checkers): m_checkers (checkers) { }
    
    virtual void Run () = 0;
    
    virtual void Show (std::ostream& o) {
      for (auto checker: this->m_checkers) {
        checker->write (o);
      }
    }
  };

  template<typename Analyzer>  
  class IntraChecker: public Checker <Analyzer> {
   public:

    typedef Checker<Analyzer> base_checker_t;
    using typename base_checker_t::prop_checker_ptr;
    using typename base_checker_t::prop_checker_vector;

   private:

    typedef typename Analyzer::cfg_t cfg_t;
    typedef typename Analyzer::abs_dom_t abs_dom_t;
    typedef typename Analyzer::abs_tr_ptr abs_tr_ptr;

    Analyzer& m_analyzer;

   public:

    IntraChecker (Analyzer& analyzer, prop_checker_vector checkers): 
        base_checker_t (checkers), m_analyzer (analyzer) { }
    
    virtual void Run () override {
      crab::ScopedCrabStats __st__("Checker");
      cfg_t cfg = m_analyzer.get_cfg ();
      for (auto &bb: cfg) {
        for (auto checker: this->m_checkers) {
          crab::ScopedCrabStats __st__("Checker." + checker->get_property_name());
          abs_dom_t inv = m_analyzer [bb.label ()];
          abs_tr_ptr abs_tr = m_analyzer.get_abs_transformer (inv);
          // propagate forward the invariants from the block entry 
          // while checking the property
          checker->set (abs_tr);
          for (auto &stmt: bb)
            stmt.accept (&*checker);
        }
      }
    }
  };


  template< typename Analyzer>
  class InterChecker: public Checker<Analyzer> {
   public:

    typedef Checker<Analyzer> base_checker_t;
    using typename base_checker_t::prop_checker_ptr;
    using typename base_checker_t::prop_checker_vector;

   private:

    typedef typename Analyzer::cg_t cg_t;
    typedef typename Analyzer::cfg_t cfg_t;
    typedef typename Analyzer::abs_tr_ptr abs_tr_ptr;
    typedef typename Analyzer::abs_dom_t abs_dom_t;

    Analyzer& m_analyzer;
    
   public:

    InterChecker (Analyzer& analyzer, prop_checker_vector checkers): 
        base_checker_t (checkers), m_analyzer (analyzer) { }
    
    virtual void Run () override {
      crab::ScopedCrabStats __st__("Checker");
      cg_t& cg = m_analyzer.get_call_graph (); 
      for (auto &v: boost::make_iterator_range (vertices (cg))) {
        cfg_t& cfg = v.getCfg ();
        for (auto &bb: cfg) {
          for (auto checker: this->m_checkers) {
            crab::ScopedCrabStats __st__("Checker." + checker->get_property_name());
            abs_dom_t inv = m_analyzer.get_pre (cfg, bb.label ());
            abs_tr_ptr abs_tr = m_analyzer.get_abs_transformer (inv);
            // propagate forward the invariants from the block entry 
            // while checking the property
            checker->set (abs_tr);
            for (auto &stmt: bb)
              stmt.accept (&*checker);
          }
        }
      }
    }
  };

  } // end namespace
} // end namespace
#endif 
