#pragma once

#include "./crab_lang.hpp"

#include <crab/domains/apron_domains.hpp>
#include <crab/domains/array_adaptive.hpp>
#include <crab/domains/array_expansion.hpp>
#include <crab/domains/array_graph.hpp>
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/boxes.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/dis_intervals.hpp>
#include <crab/domains/elina_domains.hpp>
#include <crab/domains/flat_boolean_domain.hpp>
#include <crab/domains/generic_abstract_domain.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/powerset_domain.hpp>
#include <crab/domains/region_domain.hpp>
#include <crab/domains/sparse_dbm.hpp>
#include <crab/domains/split_dbm.hpp>
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
using z_ric_domain_t = numerical_congruence_domain<z_interval_domain_t>;
using z_SparseGraph =
    DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::adapt_ss>;
using z_dbm_domain_t = sparse_dbm_domain<z_number, varname_t, z_SparseGraph>;
using z_SplitGraph =
    DBM_impl::DefaultParams<z_number, DBM_impl::GraphRep::adapt_ss>;
using z_sdbm_domain_t = split_dbm_domain<z_number, varname_t, z_SplitGraph>;
using z_boxes_domain_t = boxes_domain<z_number, varname_t>;
using z_dis_interval_domain_t = dis_interval_domain<z_number, varname_t>;
using z_box_apron_domain_t = apron_domain<z_number, varname_t, APRON_INT>;
using z_oct_apron_domain_t = apron_domain<z_number, varname_t, APRON_OCT>;
using z_pk_apron_domain_t = apron_domain<z_number, varname_t, APRON_PK>;
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
// Boolean-numerical domain over integers
using z_bool_num_domain_t = flat_boolean_numerical_domain<z_dbm_domain_t>;
using z_bool_interval_domain_t =
    flat_boolean_numerical_domain<z_interval_domain_t>;
/*===================================================================*/
// Arrays domains
/*===================================================================*/
class ArrayAdaptParams {
public:
  enum { is_smashable = 1 };
  enum { smash_at_nonzero_offset = 0 };
  enum { max_smashable_cells = 64 };
  enum { max_array_size = 512 };
};
using z_aa_int_t = array_adaptive_domain<z_interval_domain_t, ArrayAdaptParams>;
using z_aa_term_int_t =
    array_adaptive_domain<z_term_domain_t, ArrayAdaptParams>;
using z_aa_bool_int_t =
    array_adaptive_domain<z_bool_interval_domain_t, ArrayAdaptParams>;
using z_ae_int_t = array_expansion_domain<z_interval_domain_t>;
using z_ae_term_int_t = array_expansion_domain<z_term_domain_t>;
using z_ae_sdbm_t = array_expansion_domain<z_sdbm_domain_t>;
using z_ae_box_apron_t = array_expansion_domain<z_box_apron_domain_t>;
using z_ae_zones_elina_t = array_expansion_domain<z_zones_elina_domain_t>;
using z_ag_sdbm_intv_t =
    array_graph_domain<z_sdbm_domain_t, z_interval_domain_t>;
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
using z_rgn_aa_int_params_t = region_domain_impl::Params<
    z_number, varname_t,
    array_adaptive_domain<
        interval_domain<z_number, typename var_allocator::varname_t>,
        ArrayAdaptParams>>;
using z_rgn_int_params_t = region_domain_impl::Params<
    z_number, varname_t,
    interval_domain<z_number, typename var_allocator::varname_t>>;
using z_rgn_sdbm_params_t = region_domain_impl::Params<
    z_number, varname_t,
    split_dbm_domain<z_number, typename var_allocator::varname_t,
                     z_SplitGraph>>;
using z_rgn_int_t = region_domain<z_rgn_int_params_t>;
using z_rgn_sdbm_t = region_domain<z_rgn_sdbm_params_t>;
using z_rgn_aa_int_t = region_domain<z_rgn_aa_int_params_t>;
/*===================================================================*/
/// Numerical domains over real
/*===================================================================*/
using q_interval_domain_t = interval_domain<q_number, varname_t>;
using q_boxes_domain_t = boxes_domain<q_number, varname_t>;
using q_pk_apron_domain_t = apron_domain<q_number, varname_t, APRON_PK>;
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
