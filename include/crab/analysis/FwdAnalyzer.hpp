#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <crab/cfg/Cfg.hpp>
#include <crab/cfg/VarFactory.hpp>
#include <crab/analysis/Liveness.hpp>
#include <crab/analysis/AbsTransformer.hpp>
#include <crab/analysis/InvTable_traits.hpp>
#include <crab/domains/domain_traits.hpp>

namespace crab {

  namespace analyzer {

    using namespace cfg;

    //! Perform a forward flow-sensitive analysis.
    template< typename CFG, typename AbsTr, typename VarFactory,
              typename InvTblValTy>
    class FwdAnalyzer: 
        public interleaved_fwd_fixpoint_iterator< typename CFG::basic_block_label_t, 
                                                  CFG, 
                                                  typename AbsTr::abs_dom_t > 
    {
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename AbsTr::abs_dom_t abs_dom_t;
      
      typedef interleaved_fwd_fixpoint_iterator<basic_block_label_t, CFG, abs_dom_t> fwd_iterator_t;
      typedef typename AbsTr::z_lin_cst_t z_lin_cst_t;
      typedef boost::unordered_map<basic_block_label_t, InvTblValTy> invariant_map_t;    
      typedef Liveness<CFG> liveness_t;     
      typedef typename liveness_t::live_set_t live_set_t;     
      
      CFG m_cfg;
      VarFactory&  m_vfac;
      liveness_t m_live;
      bool m_keep_shadows;
      // Preserve invariants at the entry and exit. This might be
      // expensive in terms of memory. As an alternative, we could
      // compute the invariants at the exit by propagating locally from
      // the invariants at the entry.
      invariant_map_t  m_pre_map;
      invariant_map_t  m_post_map;
      
      
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
      {
        auto it = m_pre_map.find (node);
        if (it == m_pre_map.end())
        {
          if (!m_keep_shadows) {
            auto shadows = m_vfac.get_shadow_vars ();
            domain_traits::forget (inv, shadows.begin (), shadows.end ());
          }
          
          typedef inv_tbl_traits<abs_dom_t,InvTblValTy> inv_tbl_t;
          if (inv.is_bottom ())
            m_pre_map.insert(make_pair (node, inv_tbl_t::bot()));
          else if (inv.is_top ())
            m_pre_map.insert(make_pair (node, inv_tbl_t::top()));
          else
            m_pre_map.insert(make_pair (node, inv_tbl_t::marshall(inv)));
        }
      }
      
      void process_post (basic_block_label_t node, abs_dom_t inv) 
      {
        auto it = m_post_map.find (node);
        if (it == m_post_map.end())
        {
          if (!m_keep_shadows) {
            auto shadows = m_vfac.get_shadow_vars ();
            domain_traits::forget (inv, shadows.begin (), shadows.end ());
          }
          
          typedef inv_tbl_traits<abs_dom_t,InvTblValTy> inv_tbl_t;
          if (inv.is_bottom ())
            m_post_map.insert(make_pair (node, inv_tbl_t::bot()));
          else if (inv.is_top ())
            m_post_map.insert(make_pair (node, inv_tbl_t::top()));
          else
            m_post_map.insert(make_pair (node, inv_tbl_t::marshall(inv)));
        }
      }
      
     public:
      
      typedef typename invariant_map_t::iterator iterator;        
      typedef typename invariant_map_t::const_iterator const_iterator;        
      
      FwdAnalyzer (CFG cfg, VarFactory& vfac, 
                   bool runLive, bool keep_shadows=false): 
          fwd_iterator_t (cfg), m_cfg (cfg), 
          m_vfac (vfac), m_live (m_cfg), 
          m_keep_shadows (keep_shadows) {
        if (runLive)
          m_live.exec ();
      }
      
      iterator       pre_begin ()       { return m_pre_map.begin(); } 
      iterator       pre_end ()         { return m_pre_map.end();   }
      const_iterator pre_begin () const { return m_pre_map.begin(); }
      const_iterator pre_end ()   const { return m_pre_map.end();   }
      
      iterator       post_begin ()       { return m_post_map.begin(); } 
      iterator       post_end ()         { return m_post_map.end();   }
      const_iterator post_begin () const { return m_post_map.begin(); }
      const_iterator post_end ()   const { return m_post_map.end();   }
      
      //! Trigger the fixpoint computation 
      void Run (abs_dom_t inv)  {
        this->run (inv); 
      }      
      
      //! Propagate invariants at the statement level
      template < typename Statement >
      abs_dom_t AnalyzeStmt (Statement s, abs_dom_t pre) {
        AbsTr vis (pre);
        vis.visit (s);
        return vis.inv ();
      }
      
      //! Return the invariants that hold at the entry of b
      InvTblValTy operator[] (basic_block_label_t b) const {
        auto it = m_pre_map.find (b);
        if (it == m_pre_map.end ())
          return inv_tbl_traits<abs_dom_t,InvTblValTy>::top();
        else
          return it->second;
      }
      
      //! Return the invariants that hold at the exit of b
      InvTblValTy get_pre (basic_block_label_t b) const { 
        return this->operator[](b);
      }
      
      //! Return the invariants that hold at the exit of b
      InvTblValTy get_post (basic_block_label_t b) const {
        auto it = m_post_map.find (b);
        if (it == m_post_map.end ())
          return inv_tbl_traits<abs_dom_t,InvTblValTy>::top();
        else 
          return it->second;      
      }
    }; 

    //! Specialized type for a numerical forward analyzer
    template<typename CFG, typename AbsNumDomain, typename VarFactory, 
             typename InvTblValTy = AbsNumDomain>  
    struct NumFwdAnalyzer {
      typedef NumAbsTransformer <AbsNumDomain> num_abs_tr_t;
      typedef FwdAnalyzer <CFG, num_abs_tr_t, VarFactory, InvTblValTy> type;
    };
  
  } // end namespace
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
