#ifndef __CFG_IMPL__
#define __CFG_IMPL__

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>

#include <ikos/common/types.hpp>

/* 
  A simple variable factory based on strings and specialized
  definition of a CFG. 
*/

namespace cfg_impl
{
  using namespace cfg;
  using namespace std;

  template<> inline std::string get_label_str(std::string e) 
  { return e; }

  class StrVariableFactory : public boost::noncopyable  
  {
    typedef var_factory_impl::VariableFactory< std::string > StrVariableFactory_t;
    std::unique_ptr< StrVariableFactory_t > m_factory; 
    
   public: 

    typedef StrVariableFactory_t::variable_t varname_t;
    typedef StrVariableFactory_t::const_var_range const_var_range;

    StrVariableFactory(): m_factory (new StrVariableFactory_t()){ }

    varname_t operator[](std::string v)
    { 
      return (*m_factory)[v];
    }

    const_var_range get_shadow_vars () const 
    {
      return m_factory->get_shadow_vars ();
    }
  }; 

  // A variable factory based on strings
  typedef StrVariableFactory VariableFactory;
  typedef typename VariableFactory::varname_t varname_t;

  // CFG
  typedef variable< z_number, varname_t >      z_var;
  typedef std::string                          basic_block_label_t;
  typedef Cfg< basic_block_label_t, varname_t> cfg_t;
  typedef cfg_t::basic_block_t                 basic_block_t;

}

#endif 
