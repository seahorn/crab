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
            return "__tmpXYZ__ " + boost::lexical_cast<string> (_id);
          }
        }

        //T get(){ return *this->_s; }

        boost::optional<T> get(){ return *_s; }

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

     public:
      typedef IndexedString variable_t;

     private:
      typedef boost::unordered_map< T, IndexedString > map_t;      
      index_t _next_id;
      map_t   _map;
      
     public:
      VariableFactory(): _next_id (1) { }
      
      VariableFactory(index_t start_id): _next_id (start_id) { }

      // for generating temporary variables without being associated
      // with a particular T: we do not cache tem
      IndexedString operator()()
      {
        IndexedString is (_next_id++, this);
        return is;
      }
      
      IndexedString operator[](T s) 
      {
        typename map_t::iterator it = _map.find (s);
        if (it == _map.end()) 
        {
          IndexedString is (boost::make_shared<T>(s), _next_id++, this);
          _map.insert (typename map_t::value_type (s, is));
          return is;
        }
        else 
        return it->second;
      }
    }; 
  } // end namespace var_factory_impl

} // end namespace

#endif 
