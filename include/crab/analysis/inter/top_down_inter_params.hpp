#pragma once
#include "inter_params.hpp"
/* DEPRECATED: use inter_analyzer_parameters instead */

namespace crab {
namespace analyzer {
template<typename CallGraph>
using top_down_inter_analyzer_parameters = inter_analyzer_parameters<CallGraph>;
} // namespace analyzer
} // namespace crab
