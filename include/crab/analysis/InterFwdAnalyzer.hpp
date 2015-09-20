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
#include <crab/cg/TopoOrder.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/analysis/FwdAnalyzer.hpp>
#include <crab/analysis/InvTable_traits.hpp>
#include <crab/analysis/InterDS.hpp>

namespace crab {

  namespace analyzer {

    using namespace cfg;
    using namespace cg;

    template< typename CG, typename VarFactory,
              // abstract domain used for the bottom-up phase
              typename BUAbsDomain, 
              // abstract domain used for the top-down phase
              typename TDAbsDomain,
              typename InvTblValTy = TDAbsDomain>
    class InterFwdAnalyzer: public boost::noncopyable {

      typedef typename CG::node_t cg_node_t;
      typedef typename cg_node_t::cfg_t cfg_t;
      typedef typename cfg_t::varname_t varname_t;

      typedef SummaryTable <cfg_t, BUAbsDomain> summ_tbl_t;
      typedef CallCtxTable <cfg_t, TDAbsDomain> call_tbl_t;
      typedef BottomUpSummNumAbsTransformer<summ_tbl_t, call_tbl_t> bu_abs_tr;
      typedef TopDownSummNumAbsTransformer<summ_tbl_t, call_tbl_t> td_abs_tr;
      typedef FwdAnalyzer <cfg_t, bu_abs_tr, VarFactory, BUAbsDomain> bu_analyzer;
      typedef FwdAnalyzer <cfg_t, td_abs_tr, VarFactory, InvTblValTy> td_analyzer;

      typedef boost::shared_ptr <td_analyzer> td_analyzer_ptr;
      typedef boost::unordered_map <std::size_t, td_analyzer_ptr> invariant_map_t;
 
      CG m_cg;
      VarFactory&  m_vfac;
      bool m_live;
      bool m_keep_shadows;
      invariant_map_t m_inv_map;

     public:
      
      InterFwdAnalyzer (CG cg, VarFactory& vfac, 
                        bool runLive, bool keep_shadows=false): 
          m_cg (cg), m_vfac (vfac),
          m_live (runLive), m_keep_shadows (keep_shadows) {
      }
      
      //! Trigger the whole analysis
      void Run (TDAbsDomain init = TDAbsDomain::top ())  {

        std::vector<cg_node_t> rev_order;
        SccGraph <CG> Scc_g (m_cg);
        rev_topo_sort (Scc_g, rev_order);
        summ_tbl_t summ_tbl;
        call_tbl_t call_tbl;
       
        CRAB_DEBUG ("Bottom-up phase ...");
        
        for (auto n: rev_order) {
          vector<cg_node_t> &scc_mems = Scc_g.getComponentMembers (n);
          for (auto m: scc_mems) {

            cfg_t& cfg = m.getCfg ();
            auto fdecl = cfg.get_func_decl ();
            assert (fdecl);

            std::string fun_name = boost::lexical_cast<string> ((*fdecl).get_func_name ());
            if (fun_name != "main") {
              CRAB_DEBUG ("--- Analyzing ", (*fdecl).get_func_name ());
              // --- run the analysis
              bu_analyzer a (cfg, m_vfac, m_live, &summ_tbl, &call_tbl) ; 
              a.Run (BUAbsDomain::top ());
              // --- build the summary
              std::set<varname_t> formals;
              for (unsigned i=0; i < (*fdecl).get_num_params();i++)
                formals.insert ((*fdecl).get_param_name (i));
              auto ret_val_opt = findReturn (cfg);
              if (ret_val_opt) 
                formals.insert (*ret_val_opt);
              // --- project onto formal parameters and return 
              auto inv = a.get_post (cfg.exit ());
              domain_traits::project (inv, formals.begin (), formals.end ());            
              if (ret_val_opt) 
                formals.erase (*ret_val_opt);
              summ_tbl.insert (*fdecl, inv, ret_val_opt, formals);
            }
          }
        } 

        //cout << summ_tbl << endl;

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
            
            td_analyzer_ptr a (new td_analyzer (cfg, m_vfac, m_live, 
                                                &summ_tbl, &call_tbl)) ;           
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

      //! Return the invariants that hold at the exit of b in cfg
      InvTblValTy get_pre (const cfg_t &cfg, 
                           typename cfg_t::basic_block_label_t b) const { 
        auto fdecl = cfg.get_func_decl ();
        assert (fdecl);
        auto const it = m_inv_map.find (CfgHasher<cfg_t>::hash(*fdecl));
        if (it != m_inv_map.end ())
          return it->second->get_pre (b);
        else
          CRAB_ERROR("InterFwdAnalyzer: Cfg not found");
      }
      
      //! Return the invariants that hold at the exit of b in cfg
      InvTblValTy get_post (const cfg_t &cfg, 
                            typename cfg_t::basic_block_label_t b) const {
        auto fdecl = cfg.get_func_decl ();
        assert (fdecl);
        auto const it = m_inv_map.find (CfgHasher<cfg_t>::hash(*fdecl));
        if (it != m_inv_map.end ())
          return it->second->get_post (b);
        else
          CRAB_ERROR("InterFwdAnalyzer: Cfg not found");
      }
        
    }; 

  } // end namespace
} // end namespace

#endif /* INTER_FWD_ANALYZER_HPP*/
