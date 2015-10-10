#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include "boost/range/algorithm/set_algorithm.hpp"

#include <crab/cfg/Cfg.hpp>
#include <crab/cfg/VarFactory.hpp>
#include <crab/analysis/Liveness.hpp>
#include <crab/analysis/AbsTransformer.hpp>
#include <crab/analysis/InvTable_traits.hpp>
#include <crab/analysis/InterDS.hpp>
#include <crab/domains/domain_traits.hpp>

namespace crab {

  namespace analyzer {

    using namespace cfg;
  
    template<typename CFG>
    inline boost::optional<typename CFG::varname_t> findReturn (CFG cfg)
    {
      typedef typename CFG::varname_t varname_t;

      if (cfg.has_exit ()) {
        auto const & bb = cfg.get_node (cfg.exit ());
        for (auto const& s : boost::make_iterator_range (bb.begin (),
                                                         bb.end ())) {
          if (s.isReturn ()) {                            
            const Return <varname_t>* ret_stmt = 
                static_cast< const Return <varname_t> *> (&s);
            return  ret_stmt->get_ret_var ();
          }
        }
      }
      return boost::optional<varname_t> ();
    }
        
    //! Perform a forward flow-sensitive analysis.
    //  AbsTr defines the abstract transfer functions as well as which
    //  operations are modelled.
    template< typename CFG, typename AbsTr, typename VarFactory, typename InvTblValTy>
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
      typedef typename liveness_t::set_t live_set_t;     


     public:

      // datastructure types in case of inter-procedural analysis
      typedef SummaryTable <CFG, typename AbsTr::summ_abs_domain_t> summ_tbl_t;
      typedef CallCtxTable <CFG, typename AbsTr::call_abs_domain_t> call_tbl_t;

     private:

      VarFactory&  m_vfac;
      liveness_t m_live;
      // datastructures needed to perform interprocedural analysis
      summ_tbl_t* m_summ_tbl;
      call_tbl_t* m_call_tbl;
      live_set_t m_formals;
      // Preserve invariants at the entry and exit. This might be
      // expensive in terms of memory. As an alternative, we could
      // compute the invariants at the exit by propagating locally from
      // the invariants at the entry.
      invariant_map_t  m_pre_map;
      invariant_map_t  m_post_map;

      void prune_dead_variables (abs_dom_t &inv, basic_block_label_t node) {
        // prune dead variables 
        if (inv.is_bottom() || inv.is_top()) return;
        auto dead = m_live.dead_exit (node);       

        dead -= m_formals;
        domain_traits::forget (inv, dead.begin (), dead.end ());
      }

      //! Given a basic block and the invariant at the entry it produces
      //! the invariant at the exit of the block.
      abs_dom_t analyze (basic_block_label_t node, abs_dom_t pre) 
      { 
        auto &b = this->get_cfg().get_node (node);
        AbsTr vis (pre, m_summ_tbl, m_call_tbl);
        for (auto &s : b) { s.accept (&vis); }
        abs_dom_t post = vis.inv ();
        prune_dead_variables (post, node);
        return post;
      } 
      
      void process_pre (basic_block_label_t node, abs_dom_t inv) 
      {
        auto it = m_pre_map.find (node);
        if (it == m_pre_map.end())
        {          
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

     public:

      // --- intra-procedural version
      FwdAnalyzer (CFG cfg, VarFactory& vfac, bool runLive):
          fwd_iterator_t (cfg), 
          m_vfac (vfac), m_live (cfg), 
          m_summ_tbl (nullptr), m_call_tbl (nullptr) {

        if (runLive) {
          m_live.exec ();
        }
      }

      // --- inter-procedural version
      FwdAnalyzer (CFG cfg, VarFactory& vfac, bool runLive, 
                   summ_tbl_t* sum_tbl, call_tbl_t* call_tbl):
          fwd_iterator_t (cfg), 
          m_vfac (vfac), m_live (cfg), 
          m_summ_tbl (sum_tbl), m_call_tbl (call_tbl) {
        if (runLive) {
          m_live.exec ();

          // --- collect formal parameters and return value (if any)
          auto fdecl = this->get_cfg ().get_func_decl ();
          assert (fdecl);
          for (unsigned i=0; i < (*fdecl).get_num_params();i++)
            m_formals += (*fdecl).get_param_name (i); 

          if (auto ret_val = findReturn (this->get_cfg ()))
            m_formals += *ret_val; 
        }
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

      //! Return the invariants that hold at the entry of b
      InvTblValTy operator[] (basic_block_label_t b) const {
        return get_pre (b);
      }
      
      //! Return the invariants that hold at the entry of b
      InvTblValTy get_pre (basic_block_label_t b) const { 
        auto it = m_pre_map.find (b);
        if (it == m_pre_map.end ())
          return inv_tbl_traits<abs_dom_t,InvTblValTy>::top();
        else
          return it->second;
      }
      
      //! Return the invariants that hold at the exit of b
      InvTblValTy get_post (basic_block_label_t b) const {
        auto it = m_post_map.find (b);
        if (it == m_post_map.end ())
          return inv_tbl_traits<abs_dom_t,InvTblValTy>::top();
        else 
          return it->second;      
      }

      //! Propagate invariants at the statement level
      template < typename Statement >
      abs_dom_t AnalyzeStmt (Statement s, abs_dom_t pre) {
        AbsTr vis (pre, m_summ_tbl, m_call_tbl); 
        vis.visit (s);
        return vis.inv ();            
      }

    }; 

    //! Specialized type for a numerical forward analyzer
    template<typename CFG, typename AbsNumDomain, typename VarFactory, 
             typename InvTblValTy = AbsNumDomain>  
    class NumFwdAnalyzer {
     private:

      typedef NumAbsTransformer <AbsNumDomain,
                                 SummaryTable<CFG,AbsNumDomain>,
                                 CallCtxTable<CFG,AbsNumDomain> > num_abs_tr_t; 
     public:

      typedef FwdAnalyzer <CFG, num_abs_tr_t, VarFactory, InvTblValTy> type;

    };
  
  } // end namespace
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
