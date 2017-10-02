#ifndef INTER_FWD_ANALYZER_HPP
#define INTER_FWD_ANALYZER_HPP

/* 
   A standard two-phase approach for _context-insensitive_
   inter-procedural analysis.
*/

#include "boost/noncopyable.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cg/cg.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/inter_fwd_analyzer_ds.hpp>
#include <crab/analysis/dataflow/liveness.hpp>

namespace crab {

  namespace analyzer {

    template< typename CG, 
              // abstract domain used for the bottom-up phase
              typename BU_Dom, 
              // abstract domain used for the top-down phase
              typename TD_Dom>
    class inter_fwd_analyzer: public boost::noncopyable {

      typedef typename CG::node_t cg_node_t;

     public:

      typedef CG cg_t;
      typedef typename cg_node_t::cfg_t cfg_t; // cfg_t is actually a wrap to Cfg&
      typedef typename cfg_t::varname_t varname_t;
      typedef typename cfg_t::number_t number_t;
      typedef liveness<cfg_t> liveness_t;     
      typedef boost::unordered_map <cfg_t, const liveness_t*> liveness_map_t;

     private:

      typedef summary_table <cfg_t, BU_Dom> summ_tbl_t;
      typedef call_ctx_table <cfg_t, TD_Dom> call_tbl_t;

     public:

      typedef bu_summ_abs_transformer<summ_tbl_t> bu_abs_tr;
      typedef td_summ_abs_transformer<summ_tbl_t, call_tbl_t> td_abs_tr;
      typedef boost::shared_ptr<td_abs_tr> td_abs_tr_ptr;
      typedef fwd_analyzer<cfg_t, bu_abs_tr> bu_analyzer;
      typedef fwd_analyzer<cfg_t, td_abs_tr> td_analyzer;
      typedef typename summ_tbl_t::Summary summary_t;
      typedef boost::shared_ptr<summary_t> summary_ptr;
      
     public:

      // for communication with checkers
      typedef TD_Dom abs_dom_t;
      typedef td_abs_tr_ptr abs_tr_ptr;

     private:

      typedef boost::shared_ptr <td_analyzer> td_analyzer_ptr;
      typedef boost::unordered_map <std::size_t, td_analyzer_ptr> invariant_map_t;

      CG m_cg;
      const liveness_map_t* m_live;
      invariant_map_t m_inv_map;
      summ_tbl_t m_summ_tbl;
      call_tbl_t m_call_tbl;
      unsigned int m_widening_delay;
      unsigned int m_descending_iters;
      size_t m_jump_set_size; // max size of the jump set (=0 if jump set disabled)
      
      const liveness_t* get_live (const cfg_t& c) {
        if (m_live) {
          auto it = m_live->find (c);
          if (it != m_live->end ())
            return it->second;
        }
        return nullptr;
      }
      
     public:
      
      inter_fwd_analyzer (CG cg, const liveness_map_t* live,
                          unsigned int widening_delay=1,
                          unsigned int descending_iters=UINT_MAX,
                          size_t jump_set_size=0)
          : m_cg (cg), m_live (live),
            m_widening_delay (widening_delay), 
            m_descending_iters (descending_iters),
            m_jump_set_size (jump_set_size) { }
      
      //! Trigger the whole analysis
      void Run (TD_Dom init = TD_Dom::top ())  {

	CRAB_VERBOSE_IF(1, crab::outs() << "Started inter-procedural analysis\n";);
        CRAB_LOG("inter", 
                 m_cg.write (crab::outs()); crab::outs () << "\n");
                 
        crab::ScopedCrabStats __st__("Inter");

        bool has_noedges = true;
        for (auto const &v: boost::make_iterator_range (vertices (m_cg))) {
          if (out_degree (v, m_cg) > 0) {
            has_noedges = false;
            break;
          }
        } 
        
        if (has_noedges) {
          // -- Special case when the program has been inlined.
          CRAB_LOG("inter",
                   crab::outs() << "Callgraph has no edges so no summaries are computed.\n");

          CRAB_LOG("inter", 
                   m_cg.write(crab::outs());
                   crab::outs () << "\n";);
                   

          for (auto &v: boost::make_iterator_range (vertices (m_cg))) {
            crab::ScopedCrabStats __st__("Inter.TopDown");

            auto cfg = v.get_cfg ();
            auto fdecl = cfg.get_func_decl ();
            assert (fdecl);
            std::string fun_name = (*fdecl).get_func_name ().str();
            if (fun_name != "main") continue;
            
            CRAB_LOG ("inter",
                   crab::outs() << "++ Analyzing function " << (*fdecl).get_func_name() << "\n");

	    auto abs_tr = boost::make_shared<td_abs_tr> (&init, &m_summ_tbl, &m_call_tbl);
            auto a = boost::make_shared<td_analyzer> (cfg, nullptr, &*abs_tr,
                                                      m_widening_delay,
                                                      m_descending_iters,
                                                      m_jump_set_size,
						      get_live (cfg));
						      
            a->Run ();
            m_inv_map.insert (std::make_pair (cfg::cfg_hasher<cfg_t>::hash(*fdecl), a));
          }
          return;
        }

        // -- General case 
        std::vector<cg_node_t> rev_order;
	graph_algo::scc_graph<CG> Scc_g (m_cg);
        graph_algo::rev_topo_sort<graph_algo::scc_graph<CG> > (Scc_g, rev_order);

        CRAB_VERBOSE_IF(1,crab::outs() << "== Bottom-up phase ...\n";);	
        for (auto n: rev_order) {
          crab::ScopedCrabStats __st__("Inter.BottomUp");
          std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members (n);
          for (auto m: scc_mems) {

            auto cfg = m.get_cfg ();
            auto fdecl = cfg.get_func_decl ();            
            assert (fdecl);

            std::string fun_name = (*fdecl).get_func_name ().str();
            if (fun_name != "main" && cfg.has_exit ()) {
	      CRAB_VERBOSE_IF(1, crab::outs() << "++ Analyzing function "
			                      << (*fdecl).get_func_name () << "\n";);
              // --- run the analysis
	      auto init_inv = BU_Dom::top ();
	      bu_abs_tr abs_tr (&init_inv, &m_summ_tbl);
              bu_analyzer a (cfg, nullptr, &abs_tr, 
                             m_widening_delay, m_descending_iters, m_jump_set_size,
			     get_live (cfg)) ; 
              a.Run ();
	      
              // --- build the summary
              std::vector<varname_t> formals, inputs, outputs;
              formals.reserve ((*fdecl).get_num_inputs() + (*fdecl).get_num_outputs());
              inputs.reserve ((*fdecl).get_num_inputs());
              outputs.reserve ((*fdecl).get_num_outputs());

              for (unsigned i=0; i < (*fdecl).get_num_inputs();i++) {
                inputs.push_back ((*fdecl).get_input_name (i));
                formals.push_back ((*fdecl).get_input_name (i));
              }
              for (unsigned i=0; i < (*fdecl).get_num_outputs();i++) {
                outputs.push_back ((*fdecl).get_output_name (i));
                formals.push_back ((*fdecl).get_output_name (i));
              }
	      
              // --- project onto formal parameters and return values
              auto inv = a.get_post (cfg.exit ());
              //crab::CrabStats::count (BU_Dom::getDomainName() + ".count.project");
              domains::domain_traits<BU_Dom>::project (inv,
                                                       formals.begin (), 
                                                       formals.end ());            

              m_summ_tbl.insert (*fdecl, inv, inputs, outputs);
            }
          }
        } 

        CRAB_VERBOSE_IF(1, crab::outs() << "== Top-down phase ...\n";);
        bool is_root = true;
        for (auto n: boost::make_iterator_range (rev_order.rbegin(),
                                                 rev_order.rend ())) {
          crab::ScopedCrabStats __st__("Inter.TopDown");
          std::vector<cg_node_t> &scc_mems = Scc_g.get_component_members (n);
          for (auto m: scc_mems) {
            auto cfg = m.get_cfg ();
            auto fdecl = cfg.get_func_decl ();
            assert (fdecl);
	    CRAB_VERBOSE_IF(1, crab::outs() << "++ Analyzing function " 
			                    << (*fdecl).get_func_name () << "\n";);
            if (scc_mems.size () > 1) {
              // If the node is recursive then what we have in m_call_tbl
              // is incomplete and therefore it is unsound to use it. To
              // remedy it, we insert another calling context with top
              // value that approximates all the possible calling
              // contexts during the recursive calls.
              m_call_tbl.insert (*fdecl, TD_Dom::top ());
            }
	   
	    auto init_inv = init;	    
            if (is_root)
              is_root = false;
            else
	    {
	      init_inv = m_call_tbl.get_call_ctx (*fdecl);	      
              CRAB_LOG("inter",
                       crab::outs() << "    Starting analysis of "
  		                    << *fdecl <<  " with " << init_inv << "\n");
            }

	    auto abs_tr = boost::make_shared<td_abs_tr>(&init_inv, &m_summ_tbl, &m_call_tbl);
            auto a = boost::make_shared<td_analyzer> (cfg, nullptr, &*abs_tr, 
                                                      m_widening_delay,
						      m_descending_iters,
						      m_jump_set_size,
						      get_live (cfg));
	    a->Run();
            m_inv_map.insert (std::make_pair (cfg::cfg_hasher<cfg_t>::hash(*fdecl), a));
          }
        }
	CRAB_VERBOSE_IF(1,crab::outs() << "Finished inter-procedural analysis\n";);	
      }

      //! return the analyzed call graph
      CG& get_call_graph () {
        return m_cg;
      }

      //! Return the invariants that hold at the entry of b in cfg
      TD_Dom get_pre (const cfg_t &cfg, 
                      typename cfg_t::basic_block_label_t b) const { 

        if (auto fdecl = cfg.get_func_decl ()) {
          auto const it = m_inv_map.find (cfg::cfg_hasher<cfg_t>::hash(*fdecl));
          if (it != m_inv_map.end ())
            return it->second->get_pre (b);
        }
        return TD_Dom::top ();
      }
      
      //! Return the invariants that hold at the exit of b in cfg
      TD_Dom get_post (const cfg_t &cfg, 
                       typename cfg_t::basic_block_label_t b) const {
        
        if (auto fdecl = cfg.get_func_decl ()) {
          auto const it = m_inv_map.find (cfg::cfg_hasher<cfg_t>::hash(*fdecl));
          if (it != m_inv_map.end ())
            return it->second->get_post (b);
        }
        return TD_Dom::top ();
      }

      //! Propagate inv through statements
      td_abs_tr_ptr get_abs_transformer (TD_Dom &inv) {
        return boost::make_shared<td_abs_tr>(&inv, &m_summ_tbl, &m_call_tbl);        
      }

      //! Return true if there is a summary for a function
      bool has_summary (const cfg_t &cfg) const {
        if (auto fdecl = cfg.get_func_decl ())
          return m_summ_tbl.hasSummary (*fdecl);
        return false;
      }

      //! Return the summary for a function
      summary_ptr get_summary(const cfg_t &cfg) const {
        if (auto fdecl = cfg.get_func_decl ()) {
          if (m_summ_tbl.hasSummary (*fdecl)) {
            summary_t summ = m_summ_tbl.get (*fdecl);
            return boost::make_shared<summary_t>(summ);
          }
        }
        CRAB_WARN("Summary not found");
        return nullptr;
      }
        
    }; 

  } // end namespace
} // end namespace

#endif /* INTER_FWD_ANALYZER_HPP*/
