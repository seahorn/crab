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

    VariableName name() const {
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

  // Enumeration type for basic arithmetic operations
  typedef enum {
    OP_ADDITION,
    OP_SUBTRACTION,
    OP_MULTIPLICATION,
    OP_DIVISION
  } operation_t;
  
} // namespace ikos

#endif // IKOS_COMMON_HPP
