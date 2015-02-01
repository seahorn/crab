/*******************************************************************************
 * Common types and classes used in IKOS.
 ******************************************************************************/


#ifndef IKOS_COMMON_HPP
#define IKOS_COMMON_HPP

#include <stdint.h>
#include <string>
#include <iostream>
#include <boost/noncopyable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>

namespace ikos 
{
  
  // Numerical type for indexed objects
  typedef uint64_t index_t;

  
  // Interface for writeable objects
  class writeable {

  public:
    virtual std::ostream& write(std::ostream& o) = 0;

    virtual ~writeable() { }

  }; // class writeable

  inline std::ostream& operator<<(std::ostream& o, writeable& x) {
    return (x.write(o));
  }


  // Exception for IKOS internal errors
  class error: public writeable {
    
  protected:
    std::string msg;
    
  private:
    error();
    
  public:
    error(std::string msg_): msg(msg_) { }
    
    std::string message() {
      return msg;
    }

    std::ostream& write(std::ostream& o) {
      o << message();
      return o;
    }
    
  }; // class error


  // Container data structure for typed variables
  template< typename Type, typename VariableName >
  class variable: public writeable {
    
  public:
    typedef variable< Type, VariableName > variable_t;

  public:
    VariableName _n;
    
  public:
    variable(VariableName n): _n(n) { }

    variable(const variable_t& v): writeable(), _n(v._n) { }
    
    variable_t& operator=(variable_t v) {
      this->_n = v._n;
      return *this;
    }

    VariableName name() {
      return _n;
    }
    
    index_t index() {
      return _n.index();
    }

    bool operator<(const variable_t& v) const {
      variable_t v1 = const_cast< variable_t& >(*this);
      variable_t v2 = const_cast< variable_t& >(v);
      return v1._n.index() < v2._n.index();
    }

    std::ostream& write(std::ostream& o) {
      o << _n;
      return o;
    }
    
  }; // class variable


  //! Simple management for variable names
  class VariableFactory : public boost::noncopyable   
  {
    
   public:
    
    class indexed_string 
    {
      friend class VariableFactory;
      
      boost::shared_ptr< std::string > _s;
      index_t _id;
      VariableFactory* _vfac;
      
      indexed_string (boost::shared_ptr< std::string > s, index_t id, VariableFactory* vfac): 
          _s(s), _id(id), _vfac (vfac) { }
      
     public:
      
      indexed_string(const indexed_string& is): 
      _s(is._s), _id(is._id), _vfac (is._vfac) { }
      
      indexed_string& operator= (indexed_string is) 
      {
        _s = is._s;
        _id = is._id;
        _vfac = is._vfac;
        
        return *this;
      }
      
      index_t index()     const { return this->_id; }
      std::string name () const { return *this->_s; }
      
      bool operator<(indexed_string s)  const { return (this->_id < s._id); }
      bool operator==(indexed_string s) const { return (this->_id == s._id); }
      
      VariableFactory& getVarFactory () { return *_vfac; }
      
      std::ostream& write(std::ostream& o) 
      {
        o << *_s;
        return o;
      }
      
      friend std::ostream& operator<<(std::ostream& o, indexed_string s) 
      {
        return s.write(o);
      }
      
    }; // class indexed_string
    
   private:
    
    typedef boost::unordered_map< std::string, indexed_string > map_t;
    
    index_t _next_id;
    map_t   _map;
    
   public: 
    
    typedef indexed_string varname_t; 
    
    VariableFactory (): _next_id(1) { }
    
    VariableFactory (index_t start_id): _next_id(start_id) { }
    
    indexed_string operator[] (std::string s) 
    {
      map_t::iterator it = _map.find(s);
      if (it == _map.end()) 
      {
        indexed_string is (boost::shared_ptr< std::string > (new std::string(s)), 
                           _next_id++, this);
        _map.insert(std::pair< std::string, indexed_string >(s, is));
        return is;
      }
      else 
        return it->second;
    }
  }; 

  typedef typename VariableFactory::varname_t varname_t;

  inline std::size_t hash_value (varname_t v)
  {
    boost::hash<index_t> hasher;
    return hasher(v.index());
  }

  // Enumeration type for basic arithmetic operations
  typedef enum {
    OP_ADDITION,
    OP_SUBTRACTION,
    OP_MULTIPLICATION,
    OP_DIVISION
  } operation_t;
  
} // namespace ikos

#endif // IKOS_COMMON_HPP
