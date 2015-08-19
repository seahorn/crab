#include <ikos/domains/term/term_util.hpp>

namespace ikos {
  namespace term {

  StrVariableFactory StrVarAlloc_col::vfac;
  static const char* col_prefix_data[] = { "_x", "_y", "_z" };
  const char** StrVarAlloc_col::col_prefix = col_prefix_data;

  }
}
