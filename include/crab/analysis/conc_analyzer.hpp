#ifndef CONCURRENT_ANALYZER_HPP
#define CONCURRENT_ANALYZER_HPP

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

#include <crab/cfg/conc_sys.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/common/debug.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/analysis/fwd_analyzer.hpp>

#include <boost/noncopyable.hpp>

namespace crab {

  namespace conc {

    using namespace cfg;
    using namespace analyzer;

    //! Global fixpoint among the threads.
    //  Currently very simple flow-insensitive abstraction `a la` Mine
    //  (VMCAI'14).
    template< typename ThreadId, typename CFG, 
              typename AbsDomain, typename VarFactory>
    class conc_analyzer: public boost::noncopyable {
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename CFG::basic_block_t basic_block_t;
      typedef typename fwd_analyzer_impl<CFG, AbsDomain, VarFactory>::type fwd_analyzer_t;
      typedef typename fwd_analyzer_t::liveness_t liveness_t;
      typedef conc_sys<ThreadId, CFG> conc_sys_t;
      
     public:
      typedef boost::unordered_map<basic_block_label_t, AbsDomain> inv_map_t;

     private:
      
      typedef boost::shared_ptr<inv_map_t> inv_map_ptr;
      typedef boost::unordered_map<ThreadId, inv_map_ptr > global_inv_map_t;
      
      conc_sys_t& m_conc_sys;
      VarFactory& m_vfac;
      bool m_run_live;
      global_inv_map_t m_global_inv;
      
      unsigned int m_widening_delay;
      unsigned int m_descending_iters;
      size_t m_jump_set_size; 

     public:
      
      conc_analyzer (conc_sys_t& conc_sys, VarFactory& vfac, bool runLive,
                     unsigned int widening_delay=1,
                     unsigned int descending_iters=UINT_MAX,
                     size_t jump_set_size=0)
          :  m_conc_sys (conc_sys), m_vfac (vfac), m_run_live (runLive),
             m_widening_delay (widening_delay),
             m_descending_iters (descending_iters),
             m_jump_set_size (jump_set_size)
      { }
      
      //! Return analysis results
      inv_map_t& getInvariants (ThreadId t) 
      {
        return *(m_global_inv [t]);
      }
      
      //! Trigger the global fixpoint computation 
      void Run (AbsDomain inv) 
      { 
        /// --- inv contains the interferences 
        bool change = true;
        unsigned num_iter = 0;
        while (change)
        {
          num_iter++;
          change=false;
          for (auto const p: m_conc_sys)
          {
            /// --- run the thread separately
            liveness_t* live = nullptr;

            if (m_run_live) {
              liveness_t ls (p.second);
              ls.exec ();
              live = &ls;
            }
              
            fwd_analyzer_t thread_analyzer (p.second, m_vfac, live,
                                            m_widening_delay,
                                            m_descending_iters,
                                            m_jump_set_size);
            thread_analyzer.Run (inv);
            
            auto locals = boost::make_iterator_range (m_conc_sys.get_locals (p.first));
            
            /// --- check if there is a change in the thread invariants
            auto it = m_global_inv.find (p.first);
            if (it == m_global_inv.end ())
            {
              inv_map_ptr inv_map (new inv_map_t ());             
              for (auto bb: boost::make_iterator_range (p.second.label_begin (), 
                                                        p.second.label_end ()))
              {
                AbsDomain bb_inv = thread_analyzer [bb];
                inv_map->insert (make_pair (bb, bb_inv));
                /// -- flow insensitive abstraction 
                inv = inv | bb_inv;
                domains::domain_traits<AbsDomain>::forget (inv, 
                                                           locals.begin (), 
                                                           locals.end ());
              }
              m_global_inv [p.first] = inv_map;
              change=true;
            }
            else
            {
              inv_map_ptr inv_map = it->second;
              for (auto bb: boost::make_iterator_range (p.second.label_begin (), 
                                                        p.second.label_end ()))
             {
               AbsDomain& old_inv = (*inv_map) [bb];           
               AbsDomain new_inv = thread_analyzer [bb];
               if (!(new_inv <= old_inv && old_inv <= new_inv))
               {
                 change=true;
                 old_inv = new_inv;
               }
               /// -- flow insensitive abstraction
               inv = inv | new_inv;
               domains::domain_traits<AbsDomain>::forget (inv, 
                                                          locals.begin (), 
                                                          locals.end ());
             }
            }
          } // end analysis of all threads
        }
        crab::outs() << "Global fixpoint reached in " << num_iter << " iterations\n";
      }      
    }; 
  } // end namespace
} // end namespace

#endif /* CONC_ANALYZER_HPP*/
