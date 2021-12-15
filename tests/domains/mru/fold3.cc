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
}

typedef struct cache {
  z_soct_domain_t cache_lines;
  z_uf_domain_t regs;
  z_uf_domain_t addrs;
  z_bool_interval_domain_t flags;
} z_cache_zones_domain_t;

using variable_vector_t = std::vector<z_var>;
using variable_map_t = ikos::patricia_tree<z_var, z_var>;

variable_vector_t Vs;
variable_vector_t Cs;

template <class Dom>
void print_dom(std::string name, Dom &dom) {
  crab::outs() << name << ":\n\t" << dom << "\n";
}

void print_cache(z_cache_zones_domain_t &cache) {
  crab::outs() << "Cache:\n\t (" << cache.cache_lines << ",\t\n";
  crab::outs() << "\t" << cache.regs << ", \t\n";
  crab::outs() << "\t" << cache.addrs << ", \t\n";
  crab::outs() << "\t" << cache.flags << ")\n";
}

template <class Dom>
variable_vector_t vars(Dom const &dom) {
  variable_map_t vars;
  auto csts = dom.to_linear_constraint_system();
  for (auto it = csts.begin(); it != csts.end(); ++it) {
    auto var_range = it->variables();
    for (auto itv = var_range.begin(); itv != var_range.end(); ++itv) {
      if (vars.find(*itv) != nullptr)
        continue;
      vars.insert(*itv, *itv);
    }
  }
  variable_vector_t res;
  for (auto it = vars.begin(); it != vars.end(); ++it) {
    res.push_back(it->first);
  }
  return res;
}

template <class Dom1, class Dom2>
Dom1 meet(Dom1 dom1, Dom2 const &dom2) {
  auto csts = dom2.to_linear_constraint_system();
  dom1 += csts;
  return dom1;
}

variable_vector_t lcs_vars(variable_vector_t left, variable_vector_t right) {
  variable_map_t vars;
  for (auto it = left.begin(); it != left.end(); ++it) {
    vars.insert(*it, *it);
  }
  for (auto it = right.begin(); it != right.end(); ++it) {
    if (vars.find(*it) != nullptr)
        continue;
    vars.insert(*it, *it);
  }
  variable_vector_t res;
  for (auto it = vars.begin(); it != vars.end(); ++it) {
    res.push_back(it->first);
  }
  return res;
}

variable_vector_t vec_minus(variable_vector_t left, variable_vector_t right) {
  variable_map_t vars;
  for (auto it = left.begin(); it != left.end(); ++it) {
    vars.insert(*it, *it);
  }
  for (auto it = right.begin(); it != right.end(); ++it) {
    if (vars.find(*it) != nullptr)
        vars.remove(*it);
  }
  variable_vector_t res;
  for (auto it = vars.begin(); it != vars.end(); ++it) {
    res.push_back(it->first);
  }
  return res;
}


void fold (z_soct_domain_t &mem, z_cache_zones_domain_t &cache) {
  auto cach_lines = cache.cache_lines;
  cach_lines.rename(Cs, Vs);
  // print_dom("cache_lines", cach_lines);
  variable_vector_t vars_cache_lines = vars(cach_lines);
  cache.regs.rename(Cs, Vs);
  auto C = meet(cach_lines, cache.regs);
  // print_dom("C", C);
  variable_vector_t vars_cache = vars(C);
  z_soct_domain_t MC = mem;
  MC.project(vars_cache);
  // print_dom("MC", MC);
  auto A = MC | C;
  // print_dom("A", A);
  auto B = mem;
  variable_vector_t mem_vars = vars(mem);
  B.project(vec_minus(mem_vars, vars_cache_lines));
  // print_dom("B", B);
  mem = A & B;
}

int main(int argc, char **argv) {
  region_domain_params p(true/*allocation_sites*/,
             true/*deallocation*/,
             true/*refine_uninitialized_regions*/,
             true/*tag_analysis*/,
             false/*is_dereferenceable*/,
             true/*skip_unknown_regions*/);
  crab_domain_params_man::get().update_params(p);

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  z_soct_domain_t mem;
  z_cache_zones_domain_t cache;

  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  
  { // fold op
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var c(vfac["c"], crab::INT_TYPE, 32);
    z_var ref(vfac["ref"], crab::REF_TYPE, 32);
    z_var ref_addr(vfac["ref.base"], crab::INT_TYPE, 32);
    z_var rgn1(vfac["V_1"], crab::INT_TYPE, 32);
    Vs.push_back(rgn1);
    z_var rgn2(vfac["V_2"], crab::INT_TYPE, 32);
    Vs.push_back(rgn2);
    z_var rgn3(vfac["V_3"], crab::INT_TYPE, 32);
    Vs.push_back(rgn3);

    z_var cused(vfac["C_used"], crab::BOOL_TYPE, 1);
    z_var cdirty(vfac["C_dirty"], crab::BOOL_TYPE, 1);
    z_var cache_line1(vfac["C_1"], crab::REG_INT_TYPE, 32);
    Cs.push_back(cache_line1);
    z_var cache_line2(vfac["C_2"], crab::REG_INT_TYPE, 32);
    Cs.push_back(cache_line2);
    z_var cache_line3(vfac["C_3"], crab::REG_INT_TYPE, 32);
    Cs.push_back(cache_line3);
    z_var cache_addr(vfac["C_base"], crab::INT_TYPE, 32);

    // mem init
    mem += (y >= z_number(6));
    mem += (y - x > z_number(0));
    mem += (c <= z_number(5));
    mem += (rgn3 >= z_number(0));
    mem += (rgn3 <= z_number(1));
    mem += (rgn1 >= z_number(0));
    mem += (rgn2 <= c);
    mem += (rgn1 - rgn2 <= z_number(0));
    // mem: { 0 <= V_1; V_1 <= V_2; V_2 <= c; 0 <= V_3 <= 1; c <= 5; y > x; y >= 6 }

    // cache init
    cache.regs.set_to_top();
    cache.addrs.set_to_top();
    cache.flags.set_to_top();
    cache.flags.assign_bool_cst(cused, z_lin_cst_t::get_false());
    // cache: ( bot, top, top, { C_used = false } )

    // Simulate ref_store(ref + 8, rgn3, 2)
    // update cache
    cache.cache_lines += (cache_line1 >= z_number(0));
    cache.cache_lines += (cache_line1 <= cache_line2);
    cache.cache_lines += (cache_line2 <= z_number(5));
    cache.cache_lines += (cache_line3 >= z_number(0));
    cache.cache_lines += (cache_line3 <= z_number(1));
    cache.addrs.assign(cache_addr, ref_addr);
    cache.flags.assign_bool_cst(cdirty, z_lin_cst_t::get_true());
    cache.flags.assign_bool_cst(cused, z_lin_cst_t::get_true());
    // strong update
    cache.cache_lines.assign(cache_line3, z_number(2));

    /*
      cache: ( { 0 <= C_1; C_1 <= C_2; C_2 <= 5; C_3 = 2 }, 
      top, 
      { C_base == ref.base },
      { C_dirty = true; C_used = true } )
    */
    crab::outs() << "Pre State: \n";
    print_dom("Mem", mem);
    print_cache(cache);
    fold(mem, cache);

    crab::outs() << "Post State: \n";
    print_dom("Mem", mem);

  }
  
  return 0;
}