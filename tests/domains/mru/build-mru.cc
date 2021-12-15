#include "../../common.hpp"
#include "../../program_options.hpp"
#include <crab/domains/uf_domain.hpp>

using namespace std;
using namespace crab::domain_impl;

namespace {
using z_rgn_zones_params_t = TestRegionParams<z_soct_domain_t>;
using z_rgn_zones_t = region_domain<z_rgn_zones_params_t>;
using z_uf_domain_t = uf_domain<z_number, varname_t>;
using z_cache_zones_params_t = region_domain<z_rgn_zones_params_t>;
// using z_cache_reg_domain_t =
//   reduced_numerical_domain_product2<z_cache_zones_params_t, z_uf_domain_t>;
} // namespace

typedef struct cache {
  z_cache_zones_params_t cache_lines;
  z_uf_domain_t regs;
  z_uf_domain_t addrs;
  z_bool_interval_domain_t flags;
} z_cache_zones_domain_t;

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

  z_rgn_zones_t mem;
  z_cache_zones_domain_t cache;

  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  { // cache update
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);
    z_var ref_addr(vfac["ref.base"], crab::INT_TYPE, 32);
    z_var rgn1(vfac["V_1"], crab::REG_INT_TYPE, 32);
    z_var rgn2(vfac["V_2"], crab::REG_INT_TYPE, 32);
    z_var_or_cst_t n5_32(z_number(5), crab::variable_type(crab::INT_TYPE, 32));

    z_var cused(vfac["C_used"], crab::BOOL_TYPE, 1);
    z_var cdirty(vfac["C_dirty"], crab::BOOL_TYPE, 1);
    z_var cache_line1(vfac["C_1"], crab::REG_INT_TYPE, 32);
    z_var cache_line2(vfac["C_2"], crab::REG_INT_TYPE, 32);
    z_var cache_addr(vfac["C_base"], crab::INT_TYPE, 32);

    // mem init
    mem += (x >= z_number(0));
    mem += (x <= z_number(1));
    mem.region_init(rgn1);
    mem.region_init(rgn2);
    // Simulate cache_lines
    cache.cache_lines.region_init(cache_line1);
    cache.cache_lines.region_init(cache_line2);
    mem += (rgn1 >= z_number(1));
    mem += (rgn1 <= z_number(2));
    mem += (rgn1 - rgn2 < z_number(0));
    // mem: { 1 <= V_1 <= 2; V_1 < V_2; 0 <= x <= 1 }

    // cache init
    cache.cache_lines.make_bottom();
    cache.regs.make_top();
    cache.addrs.make_top();
    cache.flags.make_top();
    cache.flags.assign_bool_cst(cused, z_lin_cst_t::get_false());
    // cache: ( bot, top, top, { C_used = false } )

    // make_ref(@V_1:region(int),8,as_0)
    mem.ref_make(ref, rgn1, size8, as_man.mk_tag());
    // Simulate ref_store(ref, rgn2, n5_32)
    // use cache
    cache.cache_lines += (cache_line1 >= z_number(1));
    cache.cache_lines += (cache_line1 <= z_number(2));
    cache.cache_lines += (cache_line1 - cache_line2 < z_number(0));
    cache.cache_lines.ref_store(ref, cache_line2, n5_32);

    cache.addrs.assign(cache_addr, ref_addr);
    cache.flags.assign_bool_cst(cdirty, z_lin_cst_t::get_true());
    cache.flags.assign_bool_cst(cused, z_lin_cst_t::get_true());
    /*
      cache: ( { 1 <= C_1 <= 2; C_1 < C_2; C_2 = 5 },
      top,
      { C_base == ref.base },
      { C_dirty = true; C_used = true } )
    */

    crab::outs() << "Mem:\n\t" << mem << "\n";
    crab::outs() << "Cache:\n\t (" << cache.cache_lines << ",\t\n";
    crab::outs() << "\t" << cache.regs << ", \t\n";
    crab::outs() << "\t" << cache.addrs << ", \t\n";
    crab::outs() << "\t" << cache.flags << ")\n";
  }

  return 0;
}