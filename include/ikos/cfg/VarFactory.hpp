#ifndef VARIABLE_FACTORY_HPP
#define VARIABLE_FACTORY_HPP

/*
 * A factory for variable names
 */

#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/iterator_range.hpp>

namespace cfg 
{

  namespace var_factory_impl
  {

    typedef uint64_t index_t;

    namespace indexed_string_impl 
    {
      template< typename T >
      inline std::string get_str(T e);

      template<> inline std::string get_str(std::string e) { return e; }
    } 

    // This variable factory creates a new variable associated to an
    // element of type T. It can also create variables that are not
    // associated to an element of type T. We call them shadow
    // variables.
    // 
    // The factory uses an integer counter to generate variable id's
    // that always increases.
    template< class T>
    class VariableFactory 
    {
       typedef VariableFactory< T > VariableFactory_t;

     public:
      
      class IndexedString 
      {
        
        template< typename Any>
        friend class VariableFactory;
        
       private:
        boost::shared_ptr< T > _s;
        index_t _id;
        VariableFactory* _vfac;

       private:
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

        std::string str() const 
        { 
          if (_s)
          {  return indexed_string_impl::get_str< T >(*_s);  }
          else
          { // unlikely prefix
            return "@shadow.var._" + boost::lexical_cast<string> (_id);
          }
        }

        boost::optional<T> get(){ 
          if (_s) return boost::optional<T>(*_s); 
          else return boost::optional<T> ();
        }

        VariableFactory& getVarFactory () { return *_vfac; }

        bool operator<(IndexedString s)  const 
        { return (_id < s._id); }

        bool operator==(IndexedString s) const 
        { return (_id == s._id);}

        std::ostream& write(std::ostream& o) 
        {
          o << str();
          return o;
        }
        
        friend std::ostream& operator<<(std::ostream& o, IndexedString s) 
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

     private:
      typedef boost::unordered_map< T, IndexedString >   t_map_t;      
      typedef boost::unordered_map< int, IndexedString > shadow_map_t;      

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
      IndexedString get (int key)
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
  } // end namespace var_factory_impl

} // end namespace

#endif 
