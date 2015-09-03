#include <crab/cfg/VarFactory.hpp>

namespace crab {
  namespace cfg {
     namespace var_factory_impl {
       StrVariableFactory StrVarAlloc_col::vfac;
       static const char* col_prefix_data[] = { "_x", "_y", "_z" };
       const char** StrVarAlloc_col::col_prefix = col_prefix_data;
     }
  }
}
