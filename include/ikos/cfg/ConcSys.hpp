#ifndef CONCURRENT_SYSTEM_HPP
#define CONCURRENT_SYSTEM_HPP

#include <boost/functional/hash.hpp>
#include <boost/noncopyable.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/unordered_map.hpp>

namespace conc
{
  using namespace std;
  using namespace cfg;

  //! A concurrent system consists of a set of threads and a set of
  //  shared global variables. The communication between threads is
  //  only through the shared variables.
  template<class BasicBlockLabel, class VariableName>
  class ConcSys : public boost::noncopyable
  {

   public:

    typedef Cfg<BasicBlockLabel, VariableName> thread_t;

   private:
    
    typedef boost::unordered_map <const thread_t *, vector<VariableName> > thread_map_t;

    typedef typename thread_map_t::value_type Binding;
    
    struct getFirst : public std::unary_function<Binding, const thread_t*>
    {
      getFirst () { }
      const thread_t* operator () (const Binding &p) const { return p.first; }
    }; 

    // struct getRefFirst : public std::unary_function<Binding, thread_t>
    // {
    //   getRefFirst () { }
    //   thread_t& operator () (const Binding &p) const { return *(p.first); }
    // }; 

    // size_t hash_value (FunctionDecl<VariableName> decl)
    // {
    //   boost::hash<const std::string> hasher; 
    //   size_t res = hasher (decl.get_func_name().name ());
    //   boost::hash_combine(res, decl.get_num_params ())
    //   return res;
    // }

    // unsigned getThreadId (const thread_t *t)
    // {
    //   auto it = m_thread_ids.find (t);
    //   if (it != m_thread_ids.end ()) return it->second;
    //   unsigned id = m_thread_ids.size ();
    //   m_thread_ids[t] = id;
    //   return id;
    // }

    thread_map_t m_thread_map;

   public:

    typedef boost::transform_iterator<getFirst, 
                                      typename thread_map_t::iterator> iterator;

    typedef boost::transform_iterator<getFirst, 
                                      typename thread_map_t::const_iterator> const_iterator;

   public:

    ConcSys () { }

    // [it, end) contains all the global shared variables 
    // read/written by t
    template<typename VariableNameIt>
    void addThread (thread_t* t, VariableNameIt begin, VariableNameIt end)
    {
      auto it = m_thread_map.find (t);
      if (it != m_thread_map.end ())
        return;
      
      vector<VariableName> shared_vs;
      std::copy (begin, end, std::back_inserter (shared_vs));

      m_thread_map.insert (make_pair (t, shared_vs));
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
   
    void write (ostream& o) const
    {
      for (auto const p: m_thread_map)
      {
        o << "Shared global vars {";
        for (auto const v: p.second)
          o << v << " " ;
        o << "}\n";
        o << "Thread: \n" << *(p.first);
      }
    }

    friend ostream& operator<<(ostream& o, const ConcSys<BasicBlockLabel, VariableName>& sys)
    {
      sys.write (o);
      return o;
    }

  };

}

#endif 
