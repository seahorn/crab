#ifndef __CFG_IMPL__
#define __CFG_IMPL__

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>

#include <ikos/common/types.hpp>

/* 
  Used for tests
*/

namespace cfg_impl
{
  using namespace cfg;
  using namespace std;

  template<> inline std::string get_label_str(std::string e) 
  { return e; }

  // A variable factory based on strings
  typedef cfg::var_factory_impl::StrVariableFactory VariableFactory;
  typedef typename VariableFactory::varname_t varname_t;

  // CFG
  typedef variable< z_number, varname_t >      z_var;
  typedef std::string                          basic_block_label_t;
  typedef Cfg< basic_block_label_t, varname_t> cfg_t;
  typedef cfg_t::basic_block_t                 basic_block_t;

}

#endif 
