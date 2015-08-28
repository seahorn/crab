#ifndef CONCURRENT_ANALYZER_HPP
#define CONCURRENT_ANALYZER_HPP

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

#include <ikos/cfg/ConcSys.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

namespace conc
{
   using namespace cfg;
   using namespace analyzer;
   using namespace std;

   //! Global fixpoint among the threads.
   //  Currently very simple flow-insensitive abstraction `a la` Mine
   //  (VMCAI'14).
   template< typename ThreadId, typename CFG, 
             typename AbsDomain, typename VarFactory,
             typename InvTableTy>
   class ConcAnalyzer {
     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::varname_t varname_t;
     typedef typename CFG::basic_block_t basic_block_t;
     typedef typename NumFwdAnalyzer<CFG, AbsDomain, VarFactory, InvTableTy>::type fwd_analyzer_t;
     typedef ConcSys <ThreadId, CFG> conc_sys_t;

     typedef inv_tbl_traits <AbsDomain, InvTableTy> inv_tbl_t;

    public:
     typedef boost::unordered_map<basic_block_label_t, AbsDomain> inv_map_t;

    private:

     typedef boost::shared_ptr<inv_map_t> inv_map_ptr;
     typedef boost::unordered_map<ThreadId, inv_map_ptr > global_inv_map_t;

     conc_sys_t& m_conc_sys;
     VarFactory& m_vfac;
     bool m_run_live;
     global_inv_map_t m_global_inv;
     bool m_keep_shadows;

    public:

     ConcAnalyzer (conc_sys_t& conc_sys, VarFactory& vfac,  
                   bool runLive, bool keep_shadows=false): 
         m_conc_sys (conc_sys), m_vfac (vfac), m_run_live (runLive), 
         m_keep_shadows (keep_shadows) 
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
           fwd_analyzer_t thread_analyzer (p.second, m_vfac, 
                                           m_run_live, m_keep_shadows);
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
               AbsDomain bb_inv = inv_tbl_t::unmarshall (thread_analyzer [bb]);
               inv_map->insert (make_pair (bb, bb_inv));
               /// -- flow insensitive abstraction 
               inv = inv | bb_inv;
               domain_traits::forget (inv, locals.begin (), locals.end ());
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
               AbsDomain new_inv = inv_tbl_t::unmarshall (thread_analyzer [bb]);
               if (!(new_inv <= old_inv && old_inv <= new_inv))
               {
                 change=true;
                 old_inv = new_inv;
               }
               /// -- flow insensitive abstraction
               inv = inv | new_inv;
               domain_traits::forget (inv, locals.begin (), locals.end ());
             }
           }
         } // end analysis of all threads
       }
       cout << "Global fixpoint reached in " << num_iter << " iterations\n";
     }      

   }; 
} // end namespace

#endif /* CONC_ANALYZER_HPP*/
