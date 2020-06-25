#ifndef __CRAB_LANGUAGE__
#define __CRAB_LANGUAGE__

#include <crab/cfg/cfg.hpp>
#include <crab/cg/cg.hpp>
#include <crab/config.h>
#include <crab/support/debug.hpp>
#include <crab/types/varname_factory.hpp>

namespace crab {

namespace cfg_impl {

/// BEGIN MUST BE DEFINED BY CRAB CLIENT
// A variable factory based on strings
using variable_factory_t = cfg::var_factory_impl::str_variable_factory;
using varname_t = typename variable_factory_t::varname_t;

// CFG basic block labels
using basic_block_label_t = std::string;
template <>
inline const std::string &get_label_str(const basic_block_label_t &bb) {
  return bb;
}
/// END MUST BE DEFINED BY CRAB CLIENT

/// To define CFG over integers
using z_cfg_t = cfg::cfg<basic_block_label_t, varname_t, ikos::z_number>;
using z_cfg_ref_t = cfg::cfg_ref<z_cfg_t>;
using z_cfg_rev_t = cfg::cfg_rev<z_cfg_ref_t>;
using z_basic_block_t = z_cfg_t::basic_block_t;
using z_var = variable<ikos::z_number, varname_t>;
using z_lin_t = ikos::linear_expression<ikos::z_number, varname_t>;
using z_lin_cst_t = ikos::linear_constraint<ikos::z_number, varname_t>;
using z_ref_cst_t = reference_constraint<ikos::z_number, varname_t>;
/// To define CFG over rationals
using q_cfg_t = cfg::cfg<basic_block_label_t, varname_t, ikos::q_number>;
using q_cfg_ref_t = cfg::cfg_ref<q_cfg_t>;
using q_cfg_rev_t = cfg::cfg_rev<q_cfg_ref_t>;
using q_basic_block_t = q_cfg_t::basic_block_t;
using q_var = variable<ikos::q_number, varname_t>;
using q_lin_t = ikos::linear_expression<ikos::q_number, varname_t>;
using q_lin_cst_t = ikos::linear_constraint<ikos::q_number, varname_t>;
using q_ref_cst_t = reference_constraint<ikos::q_number, varname_t>;
} // namespace cfg_impl

namespace cg_impl {
/// To define CG over integers
using z_cg_t = cg::call_graph<cfg_impl::z_cfg_ref_t>;
using z_cg_ref_t = cg::call_graph_ref<z_cg_t>;
} // namespace cg_impl

} // namespace crab
#endif /*__CRAB_LANGUAGE__*/
