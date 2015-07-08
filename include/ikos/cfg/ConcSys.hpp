#ifndef CONCURRENT_SYSTEM_HPP
#define CONCURRENT_SYSTEM_HPP

#include <boost/functional/hash.hpp>
#include <boost/noncopyable.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>
#include "boost/range/algorithm/set_algorithm.hpp"

namespace conc
{
  using namespace std;
  using namespace cfg;

  template <typename Set>
  void set_difference (Set &s1, Set &s2)
  {
    Set s3;
    boost::set_difference (s1, s2, std::inserter (s3, s3.end ()));
    std::swap (s3, s1);
  }
  

  //! A concurrent system consists of a set of threads and a set of
  //  shared global variables. The communication between threads is
  //  only through the shared global variables.
  template<class BasicBlockLabel, class VariableName>
  class ConcSys : public boost::noncopyable
  {

   public:
    
    typedef Cfg<BasicBlockLabel, VariableName> thread_t;
    
   private:
    
    typedef boost::unordered_map <const thread_t *, 
                                  pair <std::set<VariableName>,
                                        std::set<VariableName> > > thread_map_t;

    typedef typename thread_map_t::value_type Binding;
    
    struct getFirst : public std::unary_function<Binding, const thread_t*>
    {
      getFirst () { }
      const thread_t* operator () (const Binding &p) const { return p.first; }
    }; 

    thread_map_t m_thread_map;

   public:

    //! iterators for threads
    typedef boost::transform_iterator<getFirst, 
                                      typename thread_map_t::iterator> iterator;

    typedef boost::transform_iterator<getFirst, 
                                      typename thread_map_t::const_iterator> const_iterator;

    //! iterators for thread's variables
    typedef typename std::set <VariableName>::iterator var_iterator;

    typedef typename std::set <VariableName>::const_iterator const_var_iterator;

   public:

    ConcSys () { }

    // [it, end) contains all the global shared variables 
    // read/written by t
    template<typename VariableNameIt>
    void addThread (thread_t* t, 
                    VariableNameIt shared_begin, 
                    VariableNameIt shared_end)
    {
      auto it = m_thread_map.find (t);
      if (it != m_thread_map.end ())
        return;
      
      std::set<VariableName> shared_vs;
      shared_vs.insert (shared_begin, shared_end);
      
      std::set<VariableName> local_vs;
      for (auto &b: boost::make_iterator_range (t->begin (), t->end ()))
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
      m_thread_map.insert (make_pair (t, make_pair (shared_vs, local_vs)));
    }

    iterator begin () 
    {
      return boost::make_transform_iterator (m_thread_map.begin (), 
                                             getFirst ());
    }

    iterator end () 
    {
      return boost::make_transform_iterator (m_thread_map.end (), 
                                             getFirst ());
    }
      
    const_iterator begin () const 
    {
      return boost::make_transform_iterator (m_thread_map.begin (), 
                                             getFirst ());
    }

    const_iterator end () const 
    {
      return boost::make_transform_iterator (m_thread_map.end (), 
                                             getFirst ());
    }

    pair<const_var_iterator, const_var_iterator> 
    get_globals (const thread_t *t) const
    {
      auto const it = m_thread_map.find (t);
      if (it == m_thread_map.end ())
        IKOS_ERROR ("Thread not found");

      auto const &p = it->second;
      return make_pair (p.first.begin (), p.first.end ());
    }

    pair<const_var_iterator, const_var_iterator> 
    get_locals (const thread_t *t) const
    {
      auto const it = m_thread_map.find (t);
      if (it == m_thread_map.end ())
        IKOS_ERROR ("Thread not found");

      auto const &p = it->second; 
      return make_pair (p.second.begin (), p.second.end ());
    }

    void write (ostream& o) const
    {
      for (auto const p: m_thread_map)
      {
        o << "Shared global vars {";
        for (auto const v: p.second.first)
          o << v << " " ;
        o << "}\n";
        o << "Thread: \n" << *(p.first);
      }
    }

    friend ostream& operator<<(ostream& o, 
                               const ConcSys<BasicBlockLabel, VariableName>& sys)
    {
      sys.write (o);
      return o;
    }

  };

}

#endif 
