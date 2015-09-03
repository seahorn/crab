#ifndef CONCURRENT_SYSTEM_HPP
#define CONCURRENT_SYSTEM_HPP

#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include "boost/range/algorithm/set_algorithm.hpp"

namespace crab {

   namespace conc_impl 
   {
     // To convert a thread id to a string
     template< typename T >
     inline std::string get_thread_id_str (T t);
   } 

  namespace conc {
     using namespace cfg;

     template <typename Set>
     void set_difference (Set &s1, Set &s2)
     {
       Set s3;
       boost::set_difference (s1, s2, std::inserter (s3, s3.end ()));
       std::swap (s3, s1);
     }
  
  
     //! A concurrent system consists of a set of CFGs and a set of
     //  shared global variables. Each CFG corresponds to a thread. The
     //  communication between threads is via shared global variables.
     template<class ThreadId, class CFG>
     class ConcSys : public boost::noncopyable
     {
       
       typedef typename CFG::varname_t varname_t;
       typedef boost::unordered_map <ThreadId, CFG> cfg_map_t;
       typedef boost::unordered_map <ThreadId, set<varname_t> > var_map_t;
       
       cfg_map_t m_cfg_map;
       var_map_t m_gv_map;
       var_map_t m_lv_map;
       
      public:
       
       typedef typename cfg_map_t::iterator iterator;
       typedef typename cfg_map_t::const_iterator const_iterator;
       typedef typename set<varname_t>::iterator var_iterator;
       typedef typename set<varname_t>::const_iterator const_var_iterator;
       
       ConcSys () { }
       
       // [begin, end) contains all the global shared variables 
       // read/written by t
       template<typename VarNameIt>
       void add_thread (ThreadId t, CFG cfg, 
                        VarNameIt sh_begin, VarNameIt sh_end)
           
       {
         auto it = m_cfg_map.find (t);
         if (it != m_cfg_map.end ())
           return;
         
         set<varname_t> shared_vs;
         shared_vs.insert (sh_begin, sh_end);
         
         set<varname_t> local_vs;
         for (auto &b: boost::make_iterator_range (cfg.begin (), cfg.end ()))
         {
           for (auto &s: b)
           {
             auto ls = s.getLive ();
             if (ls.defs_begin () != ls.defs_end ())
               local_vs.insert (ls.defs_begin (), ls.defs_end ());
             if (ls.uses_begin () != ls.uses_end ())
               local_vs.insert (ls.uses_begin (), ls.uses_end ());
           }
         }
         set_difference (local_vs, shared_vs);
         
         m_cfg_map.insert (make_pair (t, cfg));
         m_gv_map.insert (make_pair (t, shared_vs));
         m_lv_map.insert (make_pair (t, local_vs));
       }
       
       iterator begin () { return m_cfg_map.begin (); }
       iterator end () { return m_cfg_map.end (); }
       const_iterator begin () const { return m_cfg_map.begin (); }
       const_iterator end () const { return m_cfg_map.end (); }
       
       pair<const_var_iterator, const_var_iterator> 
       get_globals (ThreadId t) const
       {
         auto const it = m_gv_map.find (t);
         if (it == m_gv_map.end ())
           CRAB_ERROR ("Thread not found");
         return make_pair (it->second.begin (), it->second.end ());
       }
       
       pair<const_var_iterator, const_var_iterator> 
       get_locals (ThreadId t) const
       {
         auto const it = m_lv_map.find (t);
         if (it == m_lv_map.end ())
           CRAB_ERROR ("Thread not found");
         return make_pair (it->second.begin (), it->second.end ());
       }
       
       void write (ostream& o) const
       {
         for (auto const p: m_cfg_map)
         {
           o << crab::conc_impl::get_thread_id_str (p.first) << "\n";
           
           typename var_map_t::const_iterator it = m_gv_map.find (p.first);
           assert (it != m_gv_map.end () && "thread not found");
           
           auto const &shared_vs = it->second;
           o << "Shared global vars {";
           for (auto const v: shared_vs)
             o << v << " " ;
           o << "}\n";
           o << "CFG: \n" << p.second;
         }
       }
       
       friend ostream& operator<<(ostream& o, 
                                  const ConcSys<ThreadId, CFG>& sys) {
         sys.write (o);
         return o;
       }
     };
  } // end namespace conc
}  // end namespace crab

#endif 
