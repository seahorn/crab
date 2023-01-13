#pragma once

#include "./crab_lang.hpp"

#include <crab/domains/apron_domains.hpp>
#include <crab/domains/array_adaptive.hpp>
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/boxes.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/combined_congruences.hpp>
#include <crab/domains/constant_domain.hpp>
#include <crab/domains/dis_intervals.hpp>
#include <crab/domains/elina_domains.hpp>
#include <crab/domains/fixed_tvpi_domain.hpp>
#include <crab/domains/flat_boolean_domain.hpp>
#include <crab/domains/generic_abstract_domain.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/lookahead_widening_domain.hpp>
#include <crab/domains/powerset_domain.hpp>
#include <crab/domains/region_domain.hpp>
#include <crab/domains/sign_domain.hpp>
#include <crab/domains/sign_constant_domain.hpp>
#include <crab/domains/sparse_dbm.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/domains/split_oct.hpp>
#include <crab/domains/term_equiv.hpp>
#include <crab/domains/wrapped_interval_domain.hpp>

namespace crab {

namespace domain_impl {

using namespace crab::cfg_impl;
using namespace crab::domains;
using namespace ikos;

using z_ref_cst_t = reference_constraint<z_number, varname_t>;
using z_lin_cst_sys_t = linear_constraint_system<z_number, varname_t>;
using q_lin_cst_sys_t = linear_constraint_system<q_number, varname_t>;
using z_interval_t = interval<z_number>;
using q_interval_t = interval<q_number>;

/*===================================================================*/
// Numerical domains over integers
/*===================================================================*/
using z_interval_domain_t = interval_domain<z_number, varname_t>;
using z_constant_domain_t = constant_domain<z_number, varname_t>;  
using z_ric_domain_t = numerical_congruence_domain<z_interval_domain_t>;
using z_dbm_graph_t = DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::adapt_ss>;
class SparseDBMParams {
public:
  enum { implement_inter_transformers = 1 };
};
class SplitDBMParams {
public:
  enum { implement_inter_transformers = 1 };
};
class SplitOctParams {
public:
  enum { implement_inter_transformers = 1 };
};
using z_dbm_domain_t = sparse_dbm_domain<z_number, varname_t, z_dbm_graph_t, SparseDBMParams>;
using z_sdbm_domain_t = split_dbm_domain<z_number, varname_t, z_dbm_graph_t, SplitDBMParams>;
using z_soct_domain_t = split_oct_domain<z_number, varname_t, z_dbm_graph_t, SplitOctParams>;  
using z_boxes_domain_t = boxes_domain<z_number, varname_t>;
using z_dis_interval_domain_t = dis_interval_domain<z_number, varname_t>;
using z_box_apron_domain_t = apron_domain<z_number, varname_t, APRON_INT>;
using z_oct_apron_domain_t = apron_domain<z_number, varname_t, APRON_OCT>;
using z_pk_apron_domain_t = apron_domain<z_number, varname_t, APRON_PK>;
using z_poly_pplite_domain_t = apron_domain<z_number, varname_t, APRON_PPLITE_POLY>;
using z_fpoly_pplite_domain_t = apron_domain<z_number, varname_t, APRON_PPLITE_FPOLY>;
using z_pset_pplite_domain_t = apron_domain<z_number, varname_t, APRON_PPLITE_PSET>;
using z_zones_elina_domain_t = elina_domain<z_number, varname_t, ELINA_ZONES>;
using z_oct_elina_domain_t = elina_domain<z_number, varname_t, ELINA_OCT>;
using z_pk_elina_domain_t = elina_domain<z_number, varname_t, ELINA_PK>;
using z_term_domain_t =
    term_domain<term::TDomInfo<z_number, varname_t, z_interval_domain_t>>;
using z_term_dbm_t =
    term_domain<term::TDomInfo<z_number, varname_t, z_sdbm_domain_t>>;
using z_term_dis_int_t =
    term_domain<term::TDomInfo<z_number, varname_t, z_dis_interval_domain_t>>;
using z_num_domain_t =
    reduced_numerical_domain_product2<z_term_dis_int_t, z_sdbm_domain_t>;
//using z_fixed_tvpi_domain_t = fixed_tvpi_domain<z_soct_domain_t>;
using z_fixed_tvpi_domain_t = fixed_tvpi_domain<z_sdbm_domain_t>;  
// Boolean-numerical domain over integers
using z_bool_num_domain_t = flat_boolean_numerical_domain<z_dbm_domain_t>;
using z_bool_interval_domain_t =
    flat_boolean_numerical_domain<z_interval_domain_t>;
// lookahead widening
using z_soct_domain_lw_t = lookahead_widening_domain<z_soct_domain_t>;    
/*===================================================================*/
// Arrays domains
/*===================================================================*/
using z_aa_int_t = array_adaptive_domain<z_interval_domain_t>;
using z_aa_term_int_t = array_adaptive_domain<z_term_domain_t>;
using z_aa_bool_int_t = array_adaptive_domain<z_bool_interval_domain_t>;
using z_aa_sdbm_t = array_adaptive_domain<z_sdbm_domain_t>;
using z_aa_box_apron_t = array_adaptive_domain<z_box_apron_domain_t>;
using z_aa_zones_elina_t = array_adaptive_domain<z_zones_elina_domain_t>;
using z_as_dis_int_t = array_smashing<z_dis_interval_domain_t>;
using z_as_sdbm_t = array_smashing<z_sdbm_domain_t>;
using z_as_bool_num_t = array_smashing<z_bool_num_domain_t>;
// completion disjunctive domains
using z_pow_aa_int_t = powerset_domain<z_aa_int_t>;
// Machine integer arithmetic domains
using z_wrapped_interval_domain_t =
    wrapped_interval_domain<z_number, varname_t>;
/*===================================================================*/
// Region domain
/*===================================================================*/
using var_allocator = crab::var_factory_impl::str_var_alloc_col;
template<class BaseAbsDom>
struct TestRegionParams {
  using number_t = z_number;
  using varname_t = crab::cfg_impl::varname_t;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;  
  using base_abstract_domain_t = BaseAbsDom;
  using base_varname_t = typename BaseAbsDom::varname_t;
};
using z_rgn_aa_int_params_t = TestRegionParams<
  array_adaptive_domain<
    interval_domain<z_number, typename var_allocator::varname_t>>>;
using z_rgn_int_params_t = TestRegionParams<
  interval_domain<z_number, typename var_allocator::varname_t>>;
using z_rgn_bool_int_params_t = TestRegionParams<
  flat_boolean_numerical_domain<
    interval_domain<z_number, typename var_allocator::varname_t>>>;
using z_rgn_sdbm_params_t = TestRegionParams<
  split_dbm_domain<z_number, typename var_allocator::varname_t, z_dbm_graph_t>>;
using z_rgn_constant_params_t = TestRegionParams<
  constant_domain<z_number, typename var_allocator::varname_t>>;
using z_rgn_sign_params_t = TestRegionParams<
  sign_domain<z_number, typename var_allocator::varname_t>>;
using z_rgn_sign_cst_params_t = TestRegionParams<
  sign_constant_domain<z_number, typename var_allocator::varname_t>>;  
using z_rgn_int_t = region_domain<z_rgn_int_params_t>;
using z_rgn_bool_int_t = region_domain<z_rgn_bool_int_params_t>;
using z_rgn_sdbm_t = region_domain<z_rgn_sdbm_params_t>;
using z_rgn_aa_int_t = region_domain<z_rgn_aa_int_params_t>;
using z_rgn_constant_t = region_domain<z_rgn_constant_params_t>;
using z_rgn_sign_t = region_domain<z_rgn_sign_params_t>;
using z_rgn_sign_constant_t = region_domain<z_rgn_sign_cst_params_t>;    
/*===================================================================*/
/// Numerical domains over real
/*===================================================================*/
using q_interval_domain_t = interval_domain<q_number, varname_t>;
using q_boxes_domain_t = boxes_domain<q_number, varname_t>;
using q_pk_apron_domain_t = apron_domain<q_number, varname_t, APRON_PK>;
using q_poly_pplite_domain_t = apron_domain<q_number, varname_t, APRON_PPLITE_POLY>;
using q_fpoly_pplite_domain_t = apron_domain<q_number, varname_t, APRON_PPLITE_FPOLY>;
using q_pset_pplite_domain_t = apron_domain<q_number, varname_t, APRON_PPLITE_PSET>;
using q_oct_apron_domain_t = apron_domain<q_number, varname_t, APRON_OCT>;
using q_pk_elina_domain_t = elina_domain<q_number, varname_t, ELINA_PK>;
using q_oct_elina_domain_t = elina_domain<q_number, varname_t, ELINA_OCT>;

/*===================================================================*/
// Wrapper for an arbitrary abstract domain
/*===================================================================*/
using z_abs_domain_t = abstract_domain_ref<z_var>;
using q_abs_domain_t = abstract_domain_ref<q_var>;
} // namespace domain_impl
} // namespace crab
