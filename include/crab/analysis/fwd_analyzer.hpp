#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/iterators/fwd_fixpoint_iterators.hpp>
#include <crab/analysis/liveness.hpp>
#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/inter_ds.hpp>
#include <crab/domains/domain_traits.hpp>

#include "boost/range/algorithm/set_algorithm.hpp"

namespace crab {

  namespace analyzer {

    using namespace cfg;
  
    template<typename CFG>
    inline boost::optional<typename CFG::varname_t> find_return_var (const CFG& cfg)
    {
      typedef typename CFG::varname_t varname_t;

      if (cfg.has_exit ()) {
        auto const &bb = cfg.get_node (cfg.exit ());
        for (auto const &s : boost::make_iterator_range (bb.begin(), bb.end())) {
          if (s.is_return ()) {                            
            const return_stmt<varname_t>* ret_stmt = static_cast<const return_stmt<varname_t> *> (&s);
            return ret_stmt->get_ret_var ();
          }
        }
      }
      return boost::optional<varname_t> ();
    }
        
    //! Perform a forward flow-sensitive analysis.
    //  AbsTr defines the abstract transfer functions as well as which
    //  operations are modelled.
    template< typename CFG, typename AbsTr, typename VarFactory>
    class fwd_analyzer: 
        public ikos::interleaved_fwd_fixpoint_iterator< typename CFG::basic_block_label_t, 
                                                        CFG, 
                                                        typename AbsTr::abs_dom_t >,
        public boost::noncopyable
    {

     public:

      typedef CFG cfg_t;
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename AbsTr::abs_dom_t abs_dom_t;
      typedef boost::shared_ptr<AbsTr> abs_tr_ptr;
      
     private:

      typedef interleaved_fwd_fixpoint_iterator<basic_block_label_t, CFG, abs_dom_t> fwd_iterator_t;
      typedef boost::unordered_map<basic_block_label_t, abs_dom_t> invariant_map_t;    

     public:

      // liveness info
      typedef liveness<CFG> liveness_t;     

     private:

      typedef typename liveness_t::set_t live_set_t;     

     public:

      // datastructure types in case of inter-procedural analysis
      typedef summary_table<CFG, typename AbsTr::summ_abs_domain_t> summ_tbl_t;
      typedef call_ctx_table<CFG, typename AbsTr::call_abs_domain_t> call_tbl_t;

      typedef typename invariant_map_t::iterator iterator;        
      typedef typename invariant_map_t::const_iterator const_iterator;        

     private:

      VarFactory&  m_vfac;
      const liveness_t* m_live;
      // Datastructures needed to perform interprocedural analysis
      // m_summ_tbl and m_call_tbl are preserved in memory so it could
      // be expensive.
      summ_tbl_t* m_summ_tbl;
      call_tbl_t* m_call_tbl;
      live_set_t m_formals;
      // Preserve invariants at the entry and exit. This might be
      // expensive in terms of memory. To mitigate this, we could
      // compute the invariants at the exit by propagating locally
      // from the invariants at the entry.
      invariant_map_t  m_pre_map;
      invariant_map_t  m_post_map;

      void prune_dead_variables (abs_dom_t &inv, basic_block_label_t node) {
        if (!m_live) return;

        if (inv.is_bottom() || inv.is_top()) return;
        auto dead = m_live->dead_exit (node);       

        dead -= m_formals;
        domains::domain_traits<abs_dom_t>::forget (inv, 
                                                   dead.begin (), 
                                                   dead.end ());
      }

      //! Given a basic block and the invariant at the entry it produces
      //! the invariant at the exit of the block.
      void analyze (basic_block_label_t node, abs_dom_t &inv) 
      { 
        auto &b = this->get_cfg().get_node (node);
        AbsTr vis (inv, m_summ_tbl, m_call_tbl);
        for (auto &s : b) { s.accept (&vis); }
        prune_dead_variables (inv, node);
      } 
      
      void process_pre (basic_block_label_t node, abs_dom_t inv) 
      {
        auto it = m_pre_map.find (node);
        if (it == m_pre_map.end())
        {          
          if (inv.is_bottom ())
            m_pre_map.insert(make_pair (node, abs_dom_t::bottom ()));
          else if (inv.is_top ())
            m_pre_map.insert(make_pair (node, abs_dom_t::top ()));
          else
            m_pre_map.insert(make_pair (node, inv));
        }
      }
      
      void process_post (basic_block_label_t node, abs_dom_t inv) 
      {
        auto it = m_post_map.find (node);
        if (it == m_post_map.end())
        {
          if (inv.is_bottom ())
            m_post_map.insert(make_pair (node, abs_dom_t::bottom ()));
          else if (inv.is_top ())
            m_post_map.insert(make_pair (node, abs_dom_t::top ()));
          else
            m_post_map.insert(make_pair (node, inv));
        }
      }
      
     public:

      // --- intra-procedural version
      // live can be nullptr if no live information is available
      fwd_analyzer(CFG cfg, VarFactory& vfac, 
                   const liveness_t* live,
                   unsigned int widening_delay=1,
                   unsigned int descending_iters=UINT_MAX,
                   size_t jump_set_size=0)
          : fwd_iterator_t (cfg, widening_delay, descending_iters, jump_set_size), 
            m_vfac (vfac), m_live (live), 
            m_summ_tbl (nullptr), m_call_tbl (nullptr) { }
      
      // --- inter-procedural version
      // live can be nullptr if no live information is available
      fwd_analyzer (CFG cfg, VarFactory& vfac, 
                    const liveness_t* live, 
                    summ_tbl_t* sum_tbl, call_tbl_t* call_tbl,
                    unsigned int widening_delay=1,
                    unsigned int descending_iters=UINT_MAX,
                    size_t jump_set_size=0)
          : fwd_iterator_t (cfg, widening_delay, descending_iters, jump_set_size), 
            m_vfac (vfac), m_live (live), 
            m_summ_tbl (sum_tbl), m_call_tbl (call_tbl) {
        
        if (live) {
          // --- collect formal parameters and return value (if any)
          auto fdecl = this->get_cfg ().get_func_decl ();
          assert (fdecl);
          for (unsigned i=0; i < (*fdecl).get_num_params();i++)
            m_formals += (*fdecl).get_param_name (i); 
        
          if (auto ret_val = find_return_var (this->get_cfg ()))
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

      //! Propagate inv through statements
      abs_tr_ptr get_abs_transformer (abs_dom_t &inv) {
        // pass inv by ref to avoid copies
        auto vis = boost::make_shared<AbsTr>(inv, m_summ_tbl, m_call_tbl);  
        return vis;
      }

      //! Return the invariants that hold at the entry of b
      abs_dom_t operator[] (basic_block_label_t b) const {
        return get_pre (b);
      }
      
      //! Return the invariants that hold at the entry of b
      abs_dom_t get_pre (basic_block_label_t b) const { 
        auto it = m_pre_map.find (b);
        if (it == m_pre_map.end ())
          return abs_dom_t::top ();
        else
          return it->second;
      }
      
      //! Return the invariants that hold at the exit of b
      abs_dom_t get_post (basic_block_label_t b) const {
        auto it = m_post_map.find (b);
        if (it == m_post_map.end ())
          return abs_dom_t::top ();
        else 
          return it->second;      
      }
    }; 

    //! Specialized type for a numerical forward analyzer with array
    //! and pointers
    // XXX: num_fwd_analyzer should be renamed to something like
    //      fwd_analyzer_impl. The prefix num is for historical
    //      reasons (originally supported only numerical operations)
    template<typename CFG, typename AbsNumDomain, typename VarFactory>
    class num_fwd_analyzer {
     private:

      typedef num_abs_transformer<AbsNumDomain,
                                  summary_table<CFG,AbsNumDomain>,
                                  call_ctx_table<CFG,AbsNumDomain> > num_abs_tr_t; 
     public:

      typedef fwd_analyzer<CFG, num_abs_tr_t, VarFactory> type;

    };
  
  } // end namespace
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
