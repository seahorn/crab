#include <crab/types/varname_factory.hpp>

namespace crab {
namespace var_factory_impl {
static const char *col_prefix_data[] = {"_x", "_y", "_z"};
const char **str_var_alloc_col::col_prefix = col_prefix_data;
} // namespace var_factory_impl
} // namespace crab
