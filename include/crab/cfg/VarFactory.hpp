#ifndef VARIABLE_FACTORY_HPP
#define VARIABLE_FACTORY_HPP

/*
 * Factories for variable names
 */

#include <ikos/common/types.hpp>
//#include <ikos/common/bignums.hpp>

#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

using namespace std;

namespace cfg 
{

  namespace var_factory_impl
  {

    namespace indexed_string_impl 
    {
      template< typename T >
      inline string get_str(T e);

      template<> inline string get_str(string e) { return e; }
    } 

    // This variable factory creates a new variable associated to an
    // element of type T. It can also create variables that are not
    // associated to an element of type T. We call them shadow
    // variables.
    // 
    // The factory uses a counter of type index_t to generate variable
    // id's that always increases.
    template< class T>
    class VariableFactory : public boost::noncopyable
    {
       typedef VariableFactory< T > VariableFactory_t;

     public:
      
      class IndexedString 
      {
        
        template< typename Any>
        friend class VariableFactory;

       public:
        typedef ikos::index_t index_t; 
        // FIXME: we should use some unlimited precision type to avoid
        // overflow. However, this change is a bit involving since we
        // need to change the algorithm api's in patricia_trees.hpp because
        // they assume ikos::index_t.
        // typedef ikos::z_number index_t;
        
       private:
        boost::shared_ptr< T > _s;
        index_t _id;
        VariableFactory* _vfac;

        IndexedString();

        IndexedString(index_t id, VariableFactory *vfac): 
            _s(0), _id(id), _vfac(vfac) { }

        IndexedString(boost::shared_ptr< T > s, index_t id, VariableFactory *vfac): 
            _s(s), _id(id), _vfac(vfac) { }
        
       public:

        IndexedString(const IndexedString& is): _s(is._s), _id(is._id), _vfac(is._vfac) { }
        
        IndexedString& operator=(IndexedString is) {
          _s = is._s;
          _id = is._id;
          _vfac = is._vfac;
          return *this;
        }
        
        index_t index() const { return this->_id; }

        string str() const 
        { 
          if (_s)
          {  return indexed_string_impl::get_str< T >(*_s);  }
          else
          { // unlikely prefix
            return "@shadow.var._" + boost::lexical_cast<string> (_id);
          }
        }

        boost::optional<T> get(){ 
          if (_s) return *_s; 
          else return boost::optional<T> ();
        }

        VariableFactory& getVarFactory () { return *_vfac; }

        bool operator<(IndexedString s)  const 
        { return (_id < s._id); }

        bool operator==(IndexedString s) const 
        { return (_id == s._id);}

        ostream& write(ostream& o) 
        {
          o << str();
          return o;
        }
        
        friend ostream& operator<<(ostream& o, IndexedString s) 
        {
          o << s.str ();
          return o;
        }

        friend size_t hash_value (IndexedString  s)
        {
          boost::hash<index_t> hasher;
          return hasher(s.index());
        }
        
      }; 

     public:

      typedef typename IndexedString::index_t index_t;

     private:

      typedef boost::unordered_map< T, IndexedString >   t_map_t;      
      typedef boost::unordered_map< index_t, IndexedString > shadow_map_t;      

      index_t _next_id;
      t_map_t _map;
      shadow_map_t _shadow_map;
      vector<IndexedString> _shadow_vars;

     public:
      typedef IndexedString variable_t;
      typedef boost::iterator_range<typename vector<IndexedString>::iterator> var_range;
      typedef boost::iterator_range<typename vector<IndexedString>::const_iterator> const_var_range;

     public:
      VariableFactory (): _next_id (1) { }
      
      VariableFactory (index_t start_id): _next_id (start_id) { }

      // return all the shadow variables created by the factory.
      const_var_range get_shadow_vars () const 
      {
        return boost::make_iterator_range (_shadow_vars.begin (),
                                           _shadow_vars.end ());
      }

      // special purpose: for generating IndexedString's without being
      // associated with a particular T (w/o caching).
      IndexedString get ()
      {
        IndexedString is (_next_id++, this);
        _shadow_vars.push_back (is);
        return is;
      }

      // special purpose: for generating IndexedString's without being
      // associated with a particular T (w/ caching).
      IndexedString get (index_t key)
      {
        auto it = _shadow_map.find (key);
        if (it == _shadow_map.end()) 
        {
          IndexedString is (_next_id++, this);
          _shadow_map.insert (typename shadow_map_t::value_type (key, is));
          _shadow_vars.push_back (is);
          return is;
        }
        else 
        return it->second;
      }

      IndexedString operator[](T s) 
      {
        auto it = _map.find (s);
        if (it == _map.end()) 
        {
          IndexedString is (boost::make_shared<T>(s), _next_id++, this);
          _map.insert (typename t_map_t::value_type (s, is));
          return is;
        }
        else 
        return it->second;
      }
    }; 

    //! Specialized factory for strings
    class StrVariableFactory : public boost::noncopyable  
    {
      typedef VariableFactory< std::string > StrVariableFactory_t;
      std::unique_ptr< StrVariableFactory_t > m_factory; 
      
     public: 
      
      typedef StrVariableFactory_t::variable_t varname_t;
      typedef StrVariableFactory_t::const_var_range const_var_range;
    
      StrVariableFactory(): m_factory (new StrVariableFactory_t()){ }
      
      varname_t operator[](std::string v) {
        return (*m_factory)[v];
      }
    
      const_var_range get_shadow_vars () const  {
        return m_factory->get_shadow_vars ();
      }
    }; 

    //! Specialized factory for integers
    class IntVariableFactory : public boost::noncopyable  
    {
     public: 
      typedef int varname_t;
      
      IntVariableFactory() { }
      
      varname_t operator[](int v) { return v; }
    }; 
  
    inline int fresh_colour(int col_x, int col_y) {
      switch(col_x)
      {
        case 0:
          {
            return col_y == 1 ? 2 : 1;
          }
        case 1:
          {
          return col_y == 0 ? 2 : 0;
          }
        case 2:
          {
            return col_y == 0 ? 1 : 0;
          }
        default:
          assert (false && "Unreachable");
          return 0;
      }
    }
  
    //! Three-coloured variable allocation. So the number of variables
    //  is bounded by 3|Tbl|, rather than always increasing.
    class StrVarAlloc_col {
      static const char** col_prefix;
     public:
      
      typedef StrVariableFactory::varname_t varname_t;
      static StrVariableFactory vfac;
      
      StrVarAlloc_col()
          : colour(0), next_id(0)
      { }
      
      StrVarAlloc_col(const StrVarAlloc_col& o)
          : colour(o.colour), next_id(o.next_id)
      { }
      
      StrVarAlloc_col(const StrVarAlloc_col& x, const StrVarAlloc_col& y)
          : colour(fresh_colour(x.colour, y.colour)),
            next_id(0)
      {
        assert(colour != x.colour);
        assert(colour != y.colour);
      }
      
      StrVarAlloc_col& operator=(const StrVarAlloc_col& x)
      {
        colour = x.colour;
        next_id = x.next_id;
        return *this;
      }
      
      StrVariableFactory::varname_t next() {
        std::stringstream ss;
        ss << col_prefix[colour] << next_id++;
        return vfac[ss.str()];
      }
      
     protected:
      int colour;
      int next_id;
    };
  

    class IntVarAlloc_col {
     public:
      typedef int varname_t;
      static IntVariableFactory vfac;

      
      IntVarAlloc_col()
          : colour(0), next_id(0)
      { }
      
      IntVarAlloc_col(const IntVarAlloc_col& o)
          : colour(o.colour), next_id(o.next_id)
      { }
      
      IntVarAlloc_col(const IntVarAlloc_col& x, const IntVarAlloc_col& y)
          : colour(fresh_colour(x.colour, y.colour)),
            next_id(0)
      {
        assert(colour != x.colour);
        assert(colour != y.colour);
      }
      
      IntVarAlloc_col& operator=(const IntVarAlloc_col& x)
      {
        colour = x.colour;
        next_id = x.next_id;
        return *this;
      }
      
      IntVariableFactory::varname_t next() {
        int id = next_id++;
        return 3*id + colour;
      }
      
     protected:
      int colour;
      int next_id;
    };

  } // end namespace var_factory_impl

} // end namespace cfg

#endif 
