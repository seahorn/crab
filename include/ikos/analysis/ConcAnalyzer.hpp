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
   //  It performs a very simple flow-insensitive abstraction `a la`
   //  Mine (VMCAI'14).
   template< typename CFG, typename AbsDomain>
   class ConcAnalyzer
   {     
     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::varname_t varname_t;
     typedef typename CFG::BasicBlock_t basic_block_t;

     typedef typename NumFwdAnalyzer<CFG, AbsDomain>::type fwd_analyzer_t;
     typedef ConcSys <basic_block_label_t, varname_t> conc_sys_t;

    public:

     typedef typename conc_sys_t::thread_t thread_t;
     typedef boost::unordered_map<basic_block_label_t, AbsDomain> inv_map_t;

    private:

     typedef boost::shared_ptr<inv_map_t> inv_map_ptr;
     typedef boost::unordered_map<const thread_t*, inv_map_ptr > global_inv_map_t;

     conc_sys_t& m_conc_sys;
     bool m_run_live;
     global_inv_map_t m_global_inv;
     
    public:

     ConcAnalyzer (conc_sys_t& conc_sys,  bool runLive=false): 
         m_conc_sys (conc_sys), m_run_live (runLive) 
     { }

     //! Return analysis results
     inv_map_t& getInvariants (const thread_t* t) 
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
         for (const thread_t* t: m_conc_sys)
         {
           /// --- run the thread separately
           fwd_analyzer_t thread_analyzer (*t, m_run_live);
           thread_analyzer.Run (inv);

           auto locals = boost::make_iterator_range (m_conc_sys.get_locals (t));

           /// --- check if there is a change in the thread invariants
           auto it = m_global_inv.find (t);
           if (it == m_global_inv.end ())
           {
             inv_map_ptr inv_map (new inv_map_t ());             
             for (auto bb: boost::make_iterator_range (t->label_begin (), 
                                                       t->label_end ()))
             {
               inv_map->insert (make_pair (bb, thread_analyzer[bb]));
               /// -- flow insensitive abstraction 
               inv = inv | thread_analyzer[bb];
               domain_traits::forget (inv, locals.begin (), locals.end ());
             }
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
