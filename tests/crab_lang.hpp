#pragma once

#include <crab/cfg/basic_block_traits.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/cg/cg.hpp>
#include <crab/config.h>
#include <crab/support/debug.hpp>
#include <crab/types/varname_factory.hpp>
#include <crab/types/tag.hpp>

/** 
 * Here we define control flow graphs, call graphs, basic blocks,
 * variable, linear expresssions, etc. Crab depends on three basic
 * parametric types that any client must be instantiated:
 *
 * - A variable name: varname_t. From a variable name, Crab builds
 *   variables which consist of a variable name and a type.
 * 
 * - A basic block label: basic_block_label_t. From a basic block
 *   label, Crab builds a basic block and from there a CFG and a call
 *   graph.
 *
 * - A number representation: number_t. This type defines the
 *   representation for numbers (gmp, machine arithmetic, etc).
 *
 **/

namespace crab {
namespace cfg_impl {

/* ===== BEGIN TO BE DEFINED BY CRAB CLIENT ===== */
// A variable factory based on strings
using variable_factory_t = var_factory_impl::str_variable_factory;
using varname_t = typename variable_factory_t::varname_t;
// CFG basic block labels
using basic_block_label_t = std::string;
/* ===== END TO BE DEFINED BY CRAB CLIENT ===== */

/// To define CFG over integers
using z_cfg_t = cfg::cfg<basic_block_label_t, varname_t, ikos::z_number>;
using z_cfg_ref_t = cfg::cfg_ref<z_cfg_t>;
using z_cfg_rev_t = cfg::cfg_rev<z_cfg_ref_t>;
using z_basic_block_t = z_cfg_t::basic_block_t;
using z_var = variable<ikos::z_number, varname_t>;
using z_var_or_cst_t = variable_or_constant<ikos::z_number, varname_t>;
using z_lin_exp_t = ikos::linear_expression<ikos::z_number, varname_t>;
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


/* ===== BEGIN TO BE DEFINED BY CRAB CLIENT ===== */

template<>
class variable_name_traits<std::string> {
public:
  static std::string to_string(std::string varname) {
    return varname;
  }
};
  
template<>
class basic_block_traits<cfg_impl::z_basic_block_t> {
public:
  using bb_label_t = typename cfg_impl::z_basic_block_t::basic_block_label_t;  
  static std::string to_string(const bb_label_t &bbl) {
    return bbl;
  }
};

template<>
class basic_block_traits<cfg_impl::q_basic_block_t> {
public:
  using bb_label_t = typename cfg_impl::q_basic_block_t::basic_block_label_t;
  static std::string to_string(const bb_label_t &bbl) {
    return bbl;
  }
};
/* ===== END TO BE DEFINED BY CRAB CLIENT ===== */
  
} // namespace crab
