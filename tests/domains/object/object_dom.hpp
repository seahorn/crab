#pragma once

#include "../../common.hpp"
#include <crab/domains/object_domain.hpp>

namespace crab {

namespace object_domain_impl {
using namespace crab::cfg_impl;
using namespace crab::domains;
using namespace ikos;
template<class BaseAbsDom>
struct TestObjectParams {
  using number_t = z_number;
  using varname_t = crab::cfg_impl::varname_t;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;  
  using base_abstract_domain_t = BaseAbsDom;
  using field_abstract_domain_t = BaseAbsDom;
};

using z_obj_sdbm_params_t = TestObjectParams<
  split_dbm_domain<z_number, typename domain_impl::var_allocator::varname_t, domain_impl::z_dbm_graph_t>>;
using z_obj_sdbm_t = object_domain<z_obj_sdbm_params_t>;
using z_obj_zones_params_t = TestObjectParams<domain_impl::z_soct_domain_t>;
using z_obj_zones_t = object_domain<z_obj_zones_params_t>;
} // namespace object_domain_impl
} // namespace crab