#include "../../common.hpp"
#include "../../program_options.hpp"
#include <crab/domains/uf_domain.hpp>

using namespace std;
using namespace crab::domain_impl;

namespace {
using z_uf_domain_t = uf_domain<z_number, varname_t>;
using z_rgn_zones_params_t = TestRegionParams<z_soct_domain_t>;
using z_rgn_zones_t = region_domain<z_rgn_zones_params_t>;
using z_mru_rgn_zones_t = mru_region_domain<z_rgn_zones_params_t>;
using variable_vector_t = std::vector<z_var>;
using variable_t = z_var;
using variable_map_t = ikos::patricia_tree<z_var, z_var>;
} // namespace

void simulate_obj_map_update_per_object(z_mru_rgn_zones_t &mru_mem, 
                                        variable_vector_t rgn_vars) {
  variable_t representative = rgn_vars[0];
  for (int i = 0, sz = rgn_vars.size(); i < sz; ++i) {
    // update object region map
    if (i == 0) {
      mru_mem.m_obj_rgn_map.insert({rgn_vars[i], rgn_vars});
    }
    mru_mem.m_rev_obj_rgn_map.insert({rgn_vars[i], representative});
  }
}

int main(int argc, char **argv) {
  region_domain_params p(true /*allocation_sites*/, true /*deallocation*/,
                         true /*refine_uninitialized_regions*/,
                         true /*tag_analysis*/, false /*is_dereferenceable*/,
                         true /*skip_unknown_regions*/);
  crab_domain_params_man::get().update_params(p);

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  z_mru_rgn_zones_t mru_mem;
  variable_vector_t obj1;
  variable_vector_t obj2;

  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  { // expand op
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);
    z_var ref_addr(vfac["ref.base"], crab::INT_TYPE, 32);

    // obj1
    z_var rgn1(vfac["V_1"], crab::REG_INT_TYPE, 32);
    obj1.push_back(rgn1);
    z_var rgn2(vfac["V_2"], crab::REG_INT_TYPE, 32);
    obj1.push_back(rgn2);

    // obj2
    z_var rgn3(vfac["U_1"], crab::REG_INT_TYPE, 32);
    obj2.push_back(rgn3);
    z_var rgn4(vfac["U_2"], crab::REG_INT_TYPE, 32);
    obj2.push_back(rgn4);

    z_var cused(vfac["C_used"], crab::BOOL_TYPE, 1);
    z_var cdirty(vfac["C_dirty"], crab::BOOL_TYPE, 1);
    z_var cache_addr(vfac["C_base"], crab::INT_TYPE, 32);

    // mem init
    mru_mem += (x >= z_number(0));
    mru_mem += (x <= z_number(1));
    mru_mem += (rgn1 >= z_number(1));
    mru_mem += (rgn1 <= z_number(2));
    mru_mem += (rgn1 - rgn2 < z_number(0));
    mru_mem += (rgn3 - rgn4 < z_number(0));
    mru_mem += (x - y < z_number(0));
    // mem: { 1 <= V_1 <= 2; V_1 < V_2; 0 <= x <= 1; x < y }

    // cache init
    mru_mem.m_regs.set_to_top();
    mru_mem.m_addrs.set_to_top();
    mru_mem.m_mem.assign_bool_cst(cused, z_lin_cst_t::get_false());
    mru_mem.m_cache_lines.set_to_bottom();
    // cache: ( bot, top, top, { C_used = false } )

    simulate_obj_map_update_per_object(mru_mem, obj1);
    simulate_obj_map_update_per_object(mru_mem, obj2);

    crab::outs() << "Pre State: \nMem";
    mru_mem.write(crab::outs());
    mru_mem.expand(rgn2);

    crab::outs() << "Post State: \nMem";
    mru_mem.write(crab::outs());

  }

  return 0;
}