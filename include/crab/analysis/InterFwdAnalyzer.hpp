#ifndef INTER_FWD_ANALYZER_HPP
#define INTER_FWD_ANALYZER_HPP

/* 
   A standard two-phase approach for context-insensitive
   inter-procedural analysis. 
*/

#include "boost/noncopyable.hpp"

///Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/cfg/Cfg.hpp>
#include <crab/cfg/VarFactory.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/analysis/graphs/Sccg.hpp>
#include <crab/analysis/graphs/TopoOrder.hpp>
#include <crab/analysis/FwdAnalyzer.hpp>
#include <crab/analysis/Liveness.hpp>
#include <crab/analysis/InterDS.hpp>

namespace crab {

  namespace analyzer {

    using namespace cfg;
    using namespace cg;
    using namespace graph_algo;

    template< typename CG, typename VarFactory,
              // abstract domain used for the bottom-up phase
              typename BUAbsDomain, 
              // abstract domain used for the top-down phase
              typename TDAbsDomain>
    class InterFwdAnalyzer: public boost::noncopyable {

      typedef typename CG::node_t cg_node_t;
      typedef typename cg_node_t::cfg_t cfg_t;
      typedef typename cfg_t::varname_t varname_t;
      typedef Liveness<cfg_t> liveness_t;     

     public:
      typedef boost::unordered_map <cfg_t, const liveness_t*> liveness_map_t;

     private:
      typedef SummaryTable <cfg_t, BUAbsDomain> summ_tbl_t;
      typedef CallCtxTable <cfg_t, TDAbsDomain> call_tbl_t;
      typedef BottomUpSummNumAbsTransformer<summ_tbl_t, call_tbl_t> bu_abs_tr;
      typedef TopDownSummNumAbsTransformer<summ_tbl_t, call_tbl_t> td_abs_tr;
      typedef FwdAnalyzer <cfg_t, bu_abs_tr, VarFactory> bu_analyzer;
      typedef FwdAnalyzer <cfg_t, td_abs_tr, VarFactory> td_analyzer;

      typedef boost::shared_ptr <td_analyzer> td_analyzer_ptr;
      typedef boost::unordered_map <std::size_t, td_analyzer_ptr> invariant_map_t;

      CG m_cg;
      VarFactory&  m_vfac;
      const liveness_map_t* m_live;
      invariant_map_t m_inv_map;
      summ_tbl_t m_summ_tbl;
      unsigned int m_widening_thres;
      unsigned int m_narrowing_iters;

      const liveness_t* get_live (cfg_t c) {
        if (m_live) {
          auto it = m_live->find (c);
          if (it != m_live->end ())
            return it->second;
        }
        return nullptr;
      }
      
     public:
      
      InterFwdAnalyzer (CG cg, VarFactory& vfac, const liveness_map_t* live,
                        unsigned int widening_thres=1,
                        unsigned int narrowing_iters=UINT_MAX): 
          m_cg (cg), m_vfac (vfac), m_live (live),
          m_widening_thres (widening_thres), 
          m_narrowing_iters (narrowing_iters) {
      }
      
      //! Trigger the whole analysis
      void Run (TDAbsDomain init = TDAbsDomain::top ())  {

        std::vector<cg_node_t> rev_order;
        SccGraph <CG> Scc_g (m_cg);
        rev_topo_sort < SccGraph <CG> > (Scc_g, rev_order);
        call_tbl_t call_tbl;
       
        CRAB_DEBUG ("Bottom-up phase ...");
        
        for (auto n: rev_order) {
          vector<cg_node_t> &scc_mems = Scc_g.getComponentMembers (n);
          for (auto m: scc_mems) {

            cfg_t& cfg = m.getCfg ();
            auto fdecl = cfg.get_func_decl ();            
            assert (fdecl);

            std::string fun_name = boost::lexical_cast<string> ((*fdecl).get_func_name ());
            if (fun_name != "main" && cfg.has_exit ()) {
              CRAB_DEBUG ("--- Analyzing ", (*fdecl).get_func_name ());
              // --- run the analysis
              bu_analyzer a (cfg, m_vfac, get_live (cfg), 
                             &m_summ_tbl, &call_tbl,
                             m_widening_thres, m_narrowing_iters) ; 
              a.Run (BUAbsDomain::top ());
              // --- build the summary
              std::vector<varname_t> formals;
              formals.reserve ((*fdecl).get_num_params());
              for (unsigned i=0; i < (*fdecl).get_num_params();i++)
                formals.push_back ((*fdecl).get_param_name (i));
              auto ret_val_opt = findReturn (cfg);
              if (ret_val_opt) 
                formals.push_back (*ret_val_opt);
              // --- project onto formal parameters and return 
              auto inv = a.get_post (cfg.exit ());
              domain_traits::project (inv, formals.begin (), formals.end ());            
              if (ret_val_opt) 
                formals.pop_back ();
              m_summ_tbl.insert (*fdecl, inv, ret_val_opt, formals);
            }
          }
        } 

        CRAB_DEBUG ("Top-down phase ...");
        bool is_root = true;
        for (auto n: boost::make_iterator_range (rev_order.rbegin(),
                                                 rev_order.rend ())) {

          vector<cg_node_t> &scc_mems = Scc_g.getComponentMembers (n);
          for (auto m: scc_mems) {
            cfg_t& cfg = m.getCfg ();
            auto fdecl = cfg.get_func_decl ();
            assert (fdecl);
            CRAB_DEBUG ("--- Analyzing ", (*fdecl).get_func_name ());
            
            if (scc_mems.size () > 1) {
              // If the node is recursive then what we have in call_tbl
              // is incomplete and therefore it is unsound to use it. To
              // remedy it, we insert another calling context with top
              // value that approximates all the possible calling
              // contexts during the recursive calls.
              call_tbl.insert (*fdecl, TDAbsDomain::top ());
            }
            
            td_analyzer_ptr a (new td_analyzer (cfg, m_vfac, get_live (cfg), 
                                                &m_summ_tbl, &call_tbl,
                                                m_widening_thres,
                                                m_narrowing_iters));           
            if (is_root) {
              a->Run (init);
              is_root = false;
            }
            else {
              CRAB_DEBUG("Starting analysis of ", 
                         *fdecl, " with ", 
                         call_tbl.get_call_ctx (*fdecl));
              a->Run (call_tbl.get_call_ctx (*fdecl));
            }
            
            m_inv_map.insert (make_pair (CfgHasher<cfg_t>::hash(*fdecl), a));
          }
        }
      }

      //! Return the invariants that hold at the entry of b in cfg
      TDAbsDomain get_pre (const cfg_t &cfg, 
                           typename cfg_t::basic_block_label_t b) const { 

        if (auto fdecl = cfg.get_func_decl ()) {
          auto const it = m_inv_map.find (CfgHasher<cfg_t>::hash(*fdecl));
          if (it != m_inv_map.end ())
            return it->second->get_pre (b);
        }
        return TDAbsDomain::top ();
      }
      
      //! Return the invariants that hold at the exit of b in cfg
      TDAbsDomain get_post (const cfg_t &cfg, 
                            typename cfg_t::basic_block_label_t b) const {
        
        if (auto fdecl = cfg.get_func_decl ()) {
          auto const it = m_inv_map.find (CfgHasher<cfg_t>::hash(*fdecl));
          if (it != m_inv_map.end ())
            return it->second->get_post (b);
        }
        return TDAbsDomain::top ();
      }

      //! Return true if there is a summary for cfg
      bool has_summary (const cfg_t &cfg) const {
        if (auto fdecl = cfg.get_func_decl ())
          return m_summ_tbl.hasSummary (*fdecl);
        return false;
      }

      //! Return the summary for cfg
      BUAbsDomain get_summary(const cfg_t &cfg) const {
        if (auto fdecl = cfg.get_func_decl ()) {
          if (m_summ_tbl.hasSummary (*fdecl)) {
            auto summ = m_summ_tbl.get (*fdecl);
            return summ.get_sum ();
          }
        }
        
        CRAB_WARN("Summary not found");
        return BUAbsDomain::top ();
      }
        
    }; 

  } // end namespace
} // end namespace

#endif /* INTER_FWD_ANALYZER_HPP*/
