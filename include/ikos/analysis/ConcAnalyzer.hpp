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

   // Global fixpoint among the threads
   template< typename BasicBlockLabel, typename VariableName,
             typename VariableFactory, typename AbsDomain>
   class ConcAnalyzer
   {
     
     typedef Cfg <BasicBlockLabel,VariableName> cfg_t;
     typedef FwdAnalyzer<BasicBlockLabel, VariableName, cfg_t, VariableFactory, AbsDomain> fwd_analyzer_t;
     typedef ConcSys <BasicBlockLabel, VariableName> conc_sys_t;

    public:
     typedef typename conc_sys_t::thread_t thread_t;
     typedef boost::unordered_map<BasicBlockLabel, AbsDomain> inv_map_t;

    private:
     typedef boost::shared_ptr<inv_map_t> inv_map_ptr;
     typedef boost::unordered_map<const thread_t*, inv_map_ptr > global_inv_map_t;

     conc_sys_t& m_conc_sys;
     VariableFactory&  m_vfac;
     bool m_run_live;
     global_inv_map_t m_global_inv;


    public:

     ConcAnalyzer (conc_sys_t& conc_sys,  VariableFactory &vfac, bool runLive=false): 
         m_conc_sys (conc_sys), m_vfac (vfac), m_run_live (runLive) 
     { }

     //! Trigger the global fixpoint computation 
     void Run (AbsDomain inv) 
     { 
       bool change = true;
       unsigned num_iter = 0;
       while (change)
       {
         num_iter++;
         change=false;
         for (const thread_t* t: m_conc_sys)
         {
           /// --- run the thread separately
           fwd_analyzer_t thread_analyzer (*t, m_vfac, m_run_live);
           thread_analyzer.Run (inv);

           /// --- check if there is a change in the thread invariants
           auto it = m_global_inv.find (t);
           if (it == m_global_inv.end ())
           {
             inv_map_ptr inv_map (new inv_map_t ());             
             for (auto bb: boost::make_iterator_range (t->label_begin (), 
                                                       t->label_end ()))
               inv_map->insert (make_pair (bb, thread_analyzer[bb]));
             m_global_inv [t] = inv_map;
             change=true;
           }
           else
           {
             inv_map_ptr inv_map = it->second;
             for (auto bb: boost::make_iterator_range (t->label_begin (), 
                                                       t->label_end ()))
             {
               AbsDomain& old_inv = (*inv_map) [bb];           
               AbsDomain new_inv = thread_analyzer [bb];
               if (!(new_inv <= old_inv && old_inv <= new_inv))
               {
                 change=true;
                 old_inv = new_inv;
               }
             }
           }
         } // end analysis of all threads
       }
     }      

     inv_map_t& getInvariants (const thread_t* t) 
     {
       return *(m_global_inv [t]);
     }

   }; 
} // end namespace

#endif /* CONC_ANALYZER_HPP*/
