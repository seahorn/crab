#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/Liveness.hpp>
#include <ikos/analysis/AbsTransformer.hpp>
#include <ikos/domains/domain_traits.hpp>

namespace analyzer
{

  using namespace cfg;
  using namespace std;

  //! Perform a forward flow-sensitive analysis.
  template< typename CFG, typename AbsTr>
  class FwdAnalyzer: 
      public interleaved_fwd_fixpoint_iterator< typename CFG::basic_block_label_t, 
                                                CFG, 
                                                typename AbsTr::abs_dom_t > 
  {
    typedef typename CFG::basic_block_label_t basic_block_label_t;
    typedef typename CFG::varname_t varname_t;
    typedef typename AbsTr::abs_dom_t abs_dom_t;

    typedef interleaved_fwd_fixpoint_iterator<basic_block_label_t, CFG, abs_dom_t> fwd_iterator_t;

    typedef boost::unordered_map<basic_block_label_t, abs_dom_t> invariant_map_t;    
    typedef Liveness<CFG> liveness_t;     
    typedef typename liveness_t::live_set_t live_set_t;     
    
    CFG              m_cfg;
    liveness_t       m_live;
    invariant_map_t  m_pre_map;
    
    //! Given a basic block and the invariant at the entry it produces
    //! the invariant at the exit of the block.
    abs_dom_t analyze (basic_block_label_t node, abs_dom_t pre) 
    { 
      auto &b = m_cfg.get_node (node);
      AbsTr vis (pre);
      for (auto &s : b) { s.accept (&vis); }
      abs_dom_t post = vis.inv ();
      
      // prune dead variables 
      if (post.is_bottom() || post.is_top()) return post;
      
      auto dead = m_live.dead_exit (node);
      domain_traits::forget (post, dead.begin (), dead.end ());
      return post;
    } 
    
    void process_pre (basic_block_label_t node, abs_dom_t inv) 
    {//cout << "Pre at " << node << ": " << inv << endl; 
      auto it = m_pre_map.find (node);
       if (it == m_pre_map.end())
         m_pre_map.insert(typename invariant_map_t::value_type (node, inv));
    }
    
    void process_post (basic_block_label_t node, abs_dom_t inv) 
    {//cout << "Post at " << node << ": " << inv << endl; 
    }
    
   public:
    
    FwdAnalyzer (CFG cfg, bool runLive=false): 
        fwd_iterator_t (cfg), m_cfg (cfg), m_live (m_cfg)
    { 
      if (runLive)
         m_live.exec ();
    }
    
    //! Trigger the fixpoint computation 
    void Run (abs_dom_t inv) 
    { 
      this->run (inv); 
    }      
     
    
    //! Propagate invariants at the statement level
    template < typename Statement >
    abs_dom_t AnalyzeStmt (Statement s, abs_dom_t pre) 
    {
      AbsTr vis (pre);
      vis.visit (s);
      return vis.inv ();
    }
    
    //! Return the invariants that hold at the entry of b
    abs_dom_t operator[] (basic_block_label_t b) const
    {
      auto it = m_pre_map.find (b);
      if (it == m_pre_map.end ())
        return abs_dom_t::top ();
      else
        return it->second;
    }
  }; 

  //! Specialized type for a numerical forward analyzer
  template<typename CFG, typename AbsNumDomain>  
  struct NumFwdAnalyzer {
    typedef NumAbsTransformer <typename CFG::varname_t, AbsNumDomain> num_abs_tr_t;
    typedef FwdAnalyzer <CFG, num_abs_tr_t> type;
  };

} // end namespace

#endif /* FWD_ANALYZER_HPP*/
