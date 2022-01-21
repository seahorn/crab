#pragma once

#include <crab/domains/flat_boolean_domain.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/types.hpp>
#include <crab/domains/uf_domain.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/types/tag.hpp>
#include <crab/types/varname_factory.hpp>

#include "region/ghost_variables.hpp"
#include "region/tags.hpp"

#include <unordered_map>

namespace crab {
namespace domains {
namespace mru_region_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;
  using base_abstract_domain_t = BaseAbsDom;
  using base_varname_t = typename BaseAbsDom::varname_t;

  static_assert(std::is_same<Number, typename BaseAbsDom::number_t>::value,
                "Number type and BaseAbsDom::number_t must be the same");
  // This is a strong requirement
  static_assert(
      std::is_same<base_varname_t,
                   typename varname_allocator_t::varname_t>::value,
      "BaseAbsDom::varname_t and allocator_varname_t must be the same");
};
} // namespace mru_region_domain_impl

template <typename Params>
class mru_region_domain
    : public abstract_domain_api<mru_region_domain<Params>> {
  using mru_region_domain_t = mru_region_domain<Params>;
  using abstract_domain_t = abstract_domain_api<mru_region_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
  using base_varname_t = typename Params::base_varname_t;
  using base_uf_domain_t = uf_domain<number_t, base_varname_t>;
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_variable_vector_t =
      typename base_abstract_domain_t::variable_vector_t;
  using base_variable_t = typename base_abstract_domain_t::variable_t;
  using base_variable_or_constant_t =
      typename base_abstract_domain_t::variable_or_constant_t;
  using base_linear_expression_t =
      typename base_abstract_domain_t::linear_expression_t;
  using base_linear_constraint_t =
      typename base_abstract_domain_t::linear_constraint_t;
  using base_linear_constraint_system_t =
      typename base_abstract_domain_t::linear_constraint_system_t;
  using base_varname_allocator_t = typename Params::varname_allocator_t;

  /**------------------ Begin type definitions ------------------**/
  using ghost_variables_t = region_domain_impl::ghost_variables<
      abstract_domain_t, base_abstract_domain_t, base_varname_allocator_t>;
  using tag_t = region_domain_impl::tag<number_t>;
  // Environment domains: map regions to finite domain
  using rgn_counting_env_t = ikos::separate_domain<variable_t, small_range>;
  using rgn_bool_env_t = ikos::separate_domain<variable_t, boolean_value>;
  using rgn_type_env_t = ikos::separate_domain<variable_t, type_value>;
  // Union-find where equivalence classes are attached to boolean values
  using rgn_dealloc_t = union_find_domain<variable_t, boolean_value>;
  // Map from variable to its ghost variables
  using var_map_t = std::unordered_map<variable_t, ghost_variables_t>;
  // Reverse map used only for pretty printing: a CrabIR variable has
  // associated a set of ghost variables. This reverse map maps back
  // each ghost variable to its CrabIR variable.
  using ghost_var_id = unsigned;
  using rev_var_map_t =
      std::unordered_map<base_variable_t, std::pair<variable_t, ghost_var_id>>;
  // Map from a representative region variable to its equivalent class
  using obj_rgn_map_t = std::unordered_map<variable_t, variable_vector_t>;
  // Reverse map used for rgn look up, by given a region variable,
  // we return its representative
  // FIXME: this could be done by using above map only.
  //      Q: does values (vectors) in map share the same value?
  using rev_obj_rgn_map_t = std::unordered_map<variable_t, variable_t>;
  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;
  using base_uf_dom_binop_t =
      std::function<base_uf_domain_t(base_uf_domain_t, base_uf_domain_t)>;
  // Map each reference or region variable to a set of allocation sites
  using alloc_site_env_t =
      separate_discrete_domain<variable_t, allocation_site>;
  using allocation_sites = typename alloc_site_env_t::value_type;
  // Map variables to sets of tags
  using tag_env_t = separate_discrete_domain<variable_t, tag_t>;
  using tag_set = typename tag_env_t::value_type;
  /**------------------ End type definitions ------------------**/

  /**------------------ Begin field definitions ------------------**/
  /*
      abs_domain mem;
      abs_domain cache_lines;
      uf_domain regs;
      uf_domain addrs;
    The implementation follows the most recently used memory model including
    cache and memory. The memory in our abstraction is represented as a domain
    and cache is a structure includeing cache lines, registers, and addresses.
    See the following comments for details:
  */

  // a special symbol to represent bottom state of the mru region domain
  // mru region domain is bottom iff mem is bottom.
  bool m_is_bottom;

  // To create ghost variables for each abstract subdomain.
  base_varname_allocator_t m_alloc;
  // Map a variable_t to its ghost variables:
  var_map_t m_var_map;
  // Reverse map from ghost variables to variable_t
  rev_var_map_t m_rev_var_map;
  // Map a region variable to its ghost variables
  // This map is used to track regions within the mru object
  var_map_t m_cache_var_map;

public:
  // The base abstract domains (over ghost variables): the domains only
  // contain variables, there is no variable of region / reference type.
  // NOTE: For each subdomain, all base variables remian as the same.
  // Any other operations do not require renaming, such as:
  //  - cache lines and regs reduction;
  //  - join / meet over cache lines and memory

  // an abstract domain that contains all the memory and register state.
  base_abstract_domain_t m_mem;
  // an abstract domain that contains all properties for the cache lines.
  // limitation:
  //    - cache is used only for one object
  base_abstract_domain_t m_cache_lines;
  // an abstract domain that models equalities between scalar variables
  // in m_mem and variables in m_cache_lines
  base_uf_domain_t m_regs;
  // an abstract domain that keeps equalities of base address between references
  // we assume each object is allocated by some address called base address.
  // each reference refers to an object so each reference has additional info
  // about base address. This domain tracks all base address alloacted so far.
  // Thus, two variables are mapped to the same uninterpreted symbol
  // if they belong to the same memory object.
  base_uf_domain_t m_addrs;

  variable_t *m_base;
  variable_t *m_used;
  variable_t *m_dirty;

private:
  // Abstract domain to count how many addresses are in a region.
  // This allows us to decide when strong update is sound: only if
  // one address per region (i.e., singleton).
  rgn_counting_env_t m_rgn_counting_dom;
  // Whether some data might have been written to any address within
  // the region.
  rgn_bool_env_t m_rgn_init_dom;
  // For each unknown region we keep track of the (dynamic) type of
  // its last written value.
  rgn_type_env_t m_rgn_type_dom;
  // Map each reference to its set of possible allocation sites. It
  // also maps each region variable to the union of all possible
  // allocation sites from all possible references stored in that
  // region. Regions need to be tracked because references are stored
  // in regions.
  alloc_site_env_t m_alloc_site_dom;
  // Keep track of whether some memory within a region has been
  // deallocated.
  //
  // To reason about deallocation, we need to know which regions might
  // belong to the same allocated memory object.  Each call to
  // ref_make models a new allocation returning a reference to the
  // base address of the allocated block. Then, ref_gep is used to
  // perform pointer arithmetic within an allocated block. Very
  // importantly, ref_gep can switch between regions although we
  // assume that those regions always belong to the same allocated
  // block.
  //
  // We partition region variables into equivalence classes attached
  // to a boolean value. Two region variables are in the same
  // equivalence class if they might belong to the same allocated
  // block. The boolean flag indicates whether the allocated block
  // might have been deallocated.
  //
  // The partitioning is done as follows:
  //
  //   region_init(rgn): create a singleton equivalence class with rgn.
  //   ref_gep(reg1, rgn1, ref2, rgn2, o): join together the
  //                                       equivalence classes of rgn1
  //                                       and rgn2.
  rgn_dealloc_t m_rgn_dealloc_dom;
  // Tag analysis: map each (any type) variable to a set of tags
  tag_env_t m_tag_env;

public: // FIXME: do not expose after unit testing
  obj_rgn_map_t m_obj_rgn_map;
  rev_obj_rgn_map_t m_rev_obj_rgn_map;

private:
  /**------------------ End field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  mru_region_domain(
      base_varname_allocator_t &&alloc, var_map_t &&var_map,
      rev_var_map_t &&rev_var_map, var_map_t &&cache_var_map,
      base_abstract_domain_t &&mem, base_abstract_domain_t &&cache_lines,
      base_uf_domain_t &&regs, base_uf_domain_t &&addrs,
      rgn_counting_env_t &&rgn_counting_dom, rgn_bool_env_t &&rgn_init_dom,
      rgn_type_env_t &&rgn_type_dom, alloc_site_env_t &&alloc_site_dom,
      rgn_dealloc_t &&rgn_dealloc_dom, tag_env_t &&tag_env,
      obj_rgn_map_t &&obj_rgn_map, rev_obj_rgn_map_t &&rev_obj_rgn_map)
      : m_is_bottom(mem.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_cache_var_map(std::move(cache_var_map)), m_mem(std::move(mem)),
        m_cache_lines(std::move(cache_lines)), m_regs(std::move(regs)),
        m_addrs(std::move(addrs)),
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
        m_rgn_type_dom(std::move(rgn_type_dom)),
        m_alloc_site_dom(std::move(alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
        m_tag_env(std::move(tag_env)), m_obj_rgn_map(std::move(obj_rgn_map)),
        m_rev_obj_rgn_map(std::move(rev_obj_rgn_map)) {}

  bool can_propagate_initialized_regions(
      const rgn_bool_env_t &left_init_dom, const var_map_t &right_varmap,
      variable_vector_t &regions,
      std::vector<ghost_variables_t> &right_base_vars) const {
    bool propagate = false;
    for (auto kv : left_init_dom) {
      if (kv.second.is_false()) {
        auto it = right_varmap.find(kv.first);
        if (it != right_varmap.end()) {
          // uninitialized on the left but initilized on the right
          regions.push_back(it->first);
          right_base_vars.push_back(it->second);
          propagate = true;
        }
      }
    }
    return propagate;
  }

  // Special step: if one region is uninitialized on left operand but
  // initialized on the right then we extract the information from the
  // right and add it to the left. This can happen if the first store
  // to a region happens inside a loop. At the entry of the loop
  // nothing is known about the region so the join there would lose
  // all the information about the region.
  //
  // This step is not optimal because it keeps only non-relational
  // information about regions but no relational information between
  // regions and other variables.
  void refine_regions(const variable_vector_t &left_regions,
                      /* the ghost variables for each regions[i] on right_dom*/
                      const std::vector<ghost_variables_t> &old_right_gvars,
                      /* left operand */
                      var_map_t &left_varmap, rev_var_map_t &left_rev_varmap,
                      base_varname_allocator_t &left_alloc,
                      base_abstract_domain_t &left_dom,
                      /* right operand */
                      const base_abstract_domain_t &right_dom,
                      base_varname_allocator_t old_right_alloc) const {
    assert(left_regions.size() == old_right_gvars.size());
    if (left_regions.empty()) {
      return;
    }

    // This is needed to avoid variable clashing with the left operand
    base_varname_allocator_t right_alloc(left_alloc, old_right_alloc);

    /* Modify the left by creating a variable in the base domain for
       the region*/
    base_variable_vector_t new_left_vars, new_right_vars;
    base_variable_vector_t old_right_vars;

    new_left_vars.reserve(left_regions.size());
    new_right_vars.reserve(left_regions.size());
    for (unsigned i = 0, sz = left_regions.size(); i < sz; ++i) {
      ghost_variables_t new_left_gvars = get_or_insert_gvars(
          left_regions[i], left_varmap, left_rev_varmap, left_alloc);
      // FIXME/TODO[if skip_unknown_regions = 0]: must pass the
      // dynamic type instead of left_regions[i].get_type(), otherwise
      // it will crash.
      ghost_variables_t new_right_gvars = ghost_variables_t::create(
          right_alloc, left_regions[i].get_type(), __LINE__);

      new_left_gvars.add(new_left_vars);
      new_right_gvars.add(new_right_vars);
      old_right_gvars[i].add(old_right_vars);
    }

    assert(old_right_vars.size() == new_right_vars.size());
    assert(new_left_vars.size() == new_right_vars.size());

    /* Propagate invariants on the region from the right to the left */
    base_abstract_domain_t dom(right_dom);
    dom.project(old_right_vars);

    // Renaming in two steps to avoid variable clashing with the left
    // operand. The assumption is that old_right_vars and
    // new_left_vars might have common names. (it might be not
    // necessary anymore).
    dom.rename(old_right_vars, new_right_vars);
    dom.rename(new_right_vars, new_left_vars);
    left_dom = left_dom & dom;
  }

  // Perform *this = join(*this, right)
  void do_join(const mru_region_domain_t &right) {

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(m_rgn_type_dom | right.m_rgn_type_dom);
    alloc_site_env_t out_alloc_site_dom(m_alloc_site_dom |
                                        right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);

    // Merge allocator
    base_varname_allocator_t out_alloc(m_alloc, right.m_alloc);

    // Perform renaming before join
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t right_addrs(right.m_addrs);

    var_map_t out_var_map; // TODO: if we agree on cache_var_map, need update
    rev_var_map_t out_rev_var_map;
    // -- Compute common renamings
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    // perform the common renamings
    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() + ".count.join.skipped.inconsistent_dynamic_type");
          continue;
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        // Store all ghost variables into each vector
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);
      }
    }

    // might need project the variables used in common because common variables
    // are subset of the variables used in domain
    m_mem.project(left_vars);
    m_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    m_cache_lines.project(left_vars);
    m_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    m_regs.project(left_vars);
    m_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    m_addrs.project(left_vars);
    m_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    m_mem |= right_mem;
    m_cache_lines |= right_cache_lines;
    m_regs |= right_regs;
    m_addrs |= right_addrs;
    m_is_bottom = m_mem.is_bottom();
    std::swap(m_alloc, out_alloc);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_rgn_counting_dom, out_rgn_counting_dom);
    std::swap(m_rgn_init_dom, out_rgn_init_dom);
    std::swap(m_rgn_type_dom, out_rgn_type_dom);
    std::swap(m_alloc_site_dom, out_alloc_site_dom);
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_tag_env, out_tag_env);

    CRAB_LOG("mru-region", crab::outs() << *this << "\n");
  }

  mru_region_domain_t do_join_or_widening(
      const mru_region_domain_t &left, const mru_region_domain_t &right,
      const bool is_join /*unused*/, base_dom_binop_t base_dom_op,
      base_uf_dom_binop_t base_uf_dom_op) const {

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(m_rgn_type_dom | right.m_rgn_type_dom);
    alloc_site_env_t out_alloc_site_dom(m_alloc_site_dom |
                                        right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);

    // Merge allocator
    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);

    // Perform renaming before join
    base_abstract_domain_t left_mem(left.m_mem);
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t left_cache_lines(left.m_cache_lines);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t left_regs(left.m_regs);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t left_addrs(left.m_addrs);
    base_uf_domain_t right_addrs(right.m_addrs);
    var_map_t out_var_map, out_cache_var_map;
    rev_var_map_t out_rev_var_map;
    obj_rgn_map_t out_obj_rgn_map;
    rev_obj_rgn_map_t out_rev_obj_rgn_map;
    // -- Compute common renamings
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = left.m_var_map.size();
    out_cache_var_map = std::move(left.m_cache_var_map);
    out_obj_rgn_map = std::move(left.m_obj_rgn_map);
    out_rev_obj_rgn_map = std::move(left.m_rev_obj_rgn_map);
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : left.m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() +
              ".count.join_or_widen.skipped.inconsistent_dynamic_type");
          continue;
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);
      }
    }

    // JN: project might be necessary to avoid keeping variables that
    // exist in the base domain but they doen't exist on m_var_map.
    //
    // If such a variable exists only on either left_dom or right_dom
    // the join removes it.  However, if we have the same variable in
    // both left_dom and righ_dom then the join will preserve it.
    //
    left_mem.project(left_vars);
    left_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    left_cache_lines.project(left_vars);
    left_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    left_regs.project(left_vars);
    left_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    left_addrs.project(left_vars);
    left_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    // Final join or widening
    base_abstract_domain_t out_mem(base_dom_op(left_mem, right_mem));
    base_abstract_domain_t out_cache_lines(
        base_dom_op(left_cache_lines, right_cache_lines));
    base_uf_domain_t out_regs(base_uf_dom_op(left_regs, right_regs));
    base_uf_domain_t out_addrs(base_uf_dom_op(left_addrs, right_addrs));

    mru_region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_cache_var_map),
        std::move(out_mem), std::move(out_cache_lines), std::move(out_regs),
        std::move(out_addrs), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
        std::move(out_tag_env), std::move(out_obj_rgn_map),
        std::move(out_rev_obj_rgn_map));
    res.m_is_bottom = res.m_mem.is_bottom();
    return res;
  }

  mru_region_domain_t do_meet_or_narrowing(
      const mru_region_domain_t &left, const mru_region_domain_t &right,
      const bool is_meet /*unused*/, base_dom_binop_t base_dom_op,
      base_uf_dom_binop_t base_uf_dom_op) const {

    // The following domains do not require common renaming
    // (i.e. no ghost variables). The domains are finite.
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(m_rgn_type_dom | right.m_rgn_type_dom);
    alloc_site_env_t out_alloc_site_dom(m_alloc_site_dom |
                                        right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);

    // Merge allocator
    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);

    // Perform renaming before meet
    base_abstract_domain_t left_mem(left.m_mem);
    base_abstract_domain_t right_mem(right.m_mem);
    base_abstract_domain_t left_cache_lines(left.m_cache_lines);
    base_abstract_domain_t right_cache_lines(right.m_cache_lines);
    base_uf_domain_t left_regs(left.m_regs);
    base_uf_domain_t right_regs(right.m_regs);
    base_uf_domain_t left_addrs(left.m_addrs);
    base_uf_domain_t right_addrs(right.m_addrs);
    var_map_t out_var_map, out_cache_var_map;
    rev_var_map_t out_rev_var_map;
    obj_rgn_map_t out_obj_rgn_map;
    rev_obj_rgn_map_t out_rev_obj_rgn_map;

    // -- Compute common renamings
    // out vars are used for output domain
    base_variable_vector_t left_vars, right_vars, out_vars;
    // Only left (right) track variables used in left (right) dom only (not in
    // common)
    base_variable_vector_t only_left_vars, only_left_out_vars;
    base_variable_vector_t only_right_vars, only_right_out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t left_renamings = left.m_var_map.size();
    size_t right_renamings = right.m_var_map.size();
    size_t num_renamings = left_renamings + right_renamings;
    out_cache_var_map = std::move(left.m_cache_var_map);
    out_obj_rgn_map = std::move(left.m_obj_rgn_map);
    out_rev_obj_rgn_map = std::move(left.m_rev_obj_rgn_map);
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    only_left_vars.reserve(left_renamings);
    only_left_out_vars.reserve(left_renamings);
    only_right_vars.reserve(right_renamings);
    only_right_out_vars.reserve(right_renamings);

    for (auto &kv : left.m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() +
              ".count.meet_or_narrowing.skipped.inconsistent_dynamic_type");
          continue;
        }
        // common renaming
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_var_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, v);
      } else {
        // only on left
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(only_left_vars);
        out_gvars.add(only_left_out_vars);
        out_var_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, kv.first);
      }
    }

    // add missing maps from right operand
    for (auto &kv : right.m_var_map) {
      const variable_t &v = kv.first;
      auto it = left.m_var_map.find(v);
      if (it == left.m_var_map.end()) {
        // only on right
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(only_right_vars);
        out_gvars.add(only_right_out_vars);
        out_var_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_var_map, kv.first);
      }
    }

    // append common into only left
    only_left_vars.insert(only_left_vars.end(), left_vars.begin(),
                          left_vars.end());
    // append common into only left out
    only_left_out_vars.insert(only_left_out_vars.end(), out_vars.begin(),
                              out_vars.end());
    // need to project before renaming
    left_mem.project(only_left_vars);
    left_mem.rename(only_left_vars, only_left_out_vars);
    left_cache_lines.project(only_left_vars);
    left_cache_lines.rename(only_left_vars, only_left_out_vars);
    left_regs.project(only_left_vars);
    left_regs.rename(only_left_vars, only_left_out_vars);
    left_addrs.project(only_left_vars);
    left_addrs.rename(only_left_vars, only_left_out_vars);

    // append common into only right
    only_right_vars.insert(only_right_vars.end(), right_vars.begin(),
                           right_vars.end());
    // append common into only right out
    only_right_out_vars.insert(only_right_out_vars.end(), out_vars.begin(),
                               out_vars.end());
    // need to project before renaming
    right_mem.project(only_left_vars);
    right_mem.rename(only_left_vars, only_left_out_vars);
    right_cache_lines.project(only_right_vars);
    right_cache_lines.rename(only_right_vars, only_right_out_vars);
    right_regs.project(only_right_vars);
    right_regs.rename(only_right_vars, only_right_out_vars);
    right_addrs.project(only_left_vars);
    right_addrs.rename(only_left_vars, only_left_out_vars);

    // Final meet or narrowing
    base_abstract_domain_t out_mem(base_dom_op(left_mem, right_mem));
    base_abstract_domain_t out_cache_lines(
        base_dom_op(left_cache_lines, right_cache_lines));
    base_uf_domain_t out_regs(base_uf_dom_op(left_regs, right_regs));
    base_uf_domain_t out_addrs(base_uf_dom_op(left_addrs, right_addrs));

    mru_region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_cache_var_map),
        std::move(out_mem), std::move(out_cache_lines), std::move(out_regs),
        std::move(out_addrs), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
        std::move(out_tag_env), std::move(out_obj_rgn_map),
        std::move(out_rev_obj_rgn_map));
    res.m_is_bottom = res.m_mem.is_bottom();
    return res;
  }

  const ghost_variables_t &get_or_insert_gvars(const variable_t &v) {
    return get_or_insert_gvars(v, m_var_map, m_rev_var_map, m_alloc);
  }

  const ghost_variables_t &
  get_or_insert_gvars(const variable_t &v, var_map_t &varmap,
                      rev_var_map_t &rev_varmap,
                      base_varname_allocator_t &alloc) const {
    auto it = varmap.find(v);
    if (it != varmap.end()) {
      return it->second;
    } else {
      variable_type vty = get_dynamic_type_or_fail(v);
      auto gvars = ghost_variables_t::create(alloc, vty, __LINE__);
      auto res = varmap.insert({v, gvars});
      if (res.second) {
        gvars.update_rev_varmap(rev_varmap, v);
      }
      return res.first->second;
    }
  }

  ghost_variables_t get_gvars_or_fail(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      CRAB_ERROR("get_gvars_or_fail failed");
    }
  }

  boost::optional<ghost_variables_t> get_gvars(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      return boost::none;
    }
  }

  base_variable_or_constant_t
  rename_variable_or_constant(const variable_or_constant_t &v) {
    if (v.is_constant()) {
      return base_variable_or_constant_t(v.get_constant(), v.get_type());
    } else {
      base_variable_t bv = get_or_insert_gvars(v.get_variable()).get_var();
      return base_variable_or_constant_t(bv);
    }
  }

  // Rename linear expression
  base_linear_expression_t rename_linear_expr(const linear_expression_t &e,
                                              bool is_cache_rename = false) {
    base_linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      const ghost_variables_t &new_gvars = get_or_insert_gvars(v);
      base_variable_t nv = new_gvars.get_var();
      m_cache_var_map.insert({v, new_gvars});
      out = out + (coef * nv);
    }
    return out;
  }

  // Rename linear constraint
  base_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst,
                                             bool is_cache_rename = false) {
    if (cst.is_inequality() || cst.is_strict_inequality()) {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression(), is_cache_rename),
          (typename base_linear_constraint_t::kind_t)cst.kind(),
          cst.is_signed());

    } else {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression(), is_cache_rename),
          (typename base_linear_constraint_t::kind_t)cst.kind());
    }
  }

  // used only for to_linear_constraint_system()
  boost::optional<variable_t> rev_rename_var(const base_variable_t &v,
                                             bool ignore_references) const {
    auto it = m_rev_var_map.find(v);
    if (it != m_rev_var_map.end()) {
      variable_t v = it->second.first;
      ghost_var_id ghost_id = it->second.second;
      if (!ignore_references || !v.get_type().is_reference()) {
        if (ghost_id == 1) { // we only get the first ghost variable
                             // which represents the address of the
                             // reference. Skipping other ghost
                             // variables such as offsets and sizes.
          return v;
        }
      }
    }
    return boost::none;
  }

  // used only for to_linear_constraint_system()
  boost::optional<linear_expression_t>
  rev_rename_linear_expr(const base_linear_expression_t &e,
                         bool ignore_references) const {
    linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const base_variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      if (boost::optional<variable_t> v_opt =
              rev_rename_var(v, ignore_references)) {
        out = out + (coef * (*v_opt));
      } else {
        return boost::optional<linear_expression_t>();
      }
    }
    return out;
  }

  // used only for to_linear_constraint_system()
  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const base_linear_constraint_t &cst,
                        bool ignore_references) const {
    if (boost::optional<linear_expression_t> e =
            rev_rename_linear_expr(cst.expression(), ignore_references)) {

      if (cst.is_inequality() || cst.is_strict_inequality()) {
        return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind(),
            cst.is_signed());

      } else {
        return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind());
      }
    } else {
      return boost::optional<linear_constraint_t>();
    }
  }

  // Rename linear constraints
  base_linear_constraint_system_t
  rename_linear_cst_sys(const linear_constraint_system_t &csts,
                        bool is_cache_rename = false) {
    base_linear_constraint_system_t out;
    for (auto const &cst : csts) {
      out += rename_linear_cst(cst, is_cache_rename);
    }
    return out;
  }

  base_linear_constraint_t
  convert_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst) {
    if (ref_cst.is_tautology()) {
      return base_linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return base_linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        base_variable_t x = get_or_insert_gvars(ref_cst.lhs()).get_var();
        if (ref_cst.is_equality()) {
          return base_linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return base_linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return base_linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return base_linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return base_linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return base_linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
        base_variable_t x = get_or_insert_gvars(ref_cst.lhs()).get_var();
        base_variable_t y = get_or_insert_gvars(ref_cst.rhs()).get_var();
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return base_linear_constraint_t(x == y + offset);
        } else if (ref_cst.is_disequality()) {
          return base_linear_constraint_t(x != y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return base_linear_constraint_t(x <= y + offset);
        } else if (ref_cst.is_less_than()) {
          return base_linear_constraint_t(x < y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return base_linear_constraint_t(x >= y + offset);
        } else if (ref_cst.is_greater_than()) {
          return base_linear_constraint_t(x > y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
  }

  template <class RangeVars>
  void merge_tags(const variable_t &x, RangeVars vars) {
    tag_set tags = tag_set::bottom();
    for (auto const &v : vars) {
      tags = tags | m_tag_env[v];
    }
    m_tag_env.set(x, tags);
  }

  bool is_tracked_region(const variable_t &v) const {
    auto ty = v.get_type();
    if (!ty.is_region()) {
      return false;
    }
    if (!ty.is_unknown_region()) {
      return true;
    }
    return (crab_domain_params_man::get().region_skip_unknown_regions()
                ? false
                : has_dynamic_type(v));
  }

  static bool is_tracked_unknown_region(const variable_t &v) {
    return (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            v.get_type().is_unknown_region());
  }

  // if the variable has a static type (different from unknown region)
  // then returns true.  If the static type is an unknown region then
  // the actual variable's type must be inferred dynamically (as part
  // of the abstract domain). In this case, it will return true if the
  // type inference done by the abstract domain succeeds (type is
  // inferred different from unknown region).
  bool has_dynamic_type(const variable_t &v) const {
    if (!v.get_type().is_unknown_region()) {
      return true;
    }

    if (crab_domain_params_man::get().region_skip_unknown_regions()) {
      return false;
    } else {
      type_value dyn_ty = m_rgn_type_dom[v];
      if (dyn_ty.is_top() || dyn_ty.is_bottom()) {
        return false;
      }
      return (!dyn_ty.get().is_unknown_region());
    }
  }

  // if has_dynamic_type(v) returns true then it returns the type,
  // otherwise it will fail.
  variable_type get_dynamic_type_or_fail(const variable_t &v) const {
    assert(has_dynamic_type(v));

    if (!v.get_type().is_unknown_region()) {
      return v.get_type();
    }

    type_value dyn_ty = m_rgn_type_dom[v];
    if (dyn_ty.is_top() || dyn_ty.is_bottom()) {
      CRAB_ERROR("get_dynamic_type_or_fail cannot be called on top or bottom");
    }
    if (dyn_ty.get().is_unknown_region()) {
      CRAB_ERROR("get_dynamic_type_or_fail should not return unknown region");
    }
    return dyn_ty.get();
  }

  static void ERROR_IF_NOT_REGION(const variable_t &v, unsigned line) {
    if (!v.get_type().is_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a region at line ", line);
    }
  }

  static void ERROR_IF_ARRAY_REGION(const variable_t &v, unsigned line) {
    if (v.get_type().is_array_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " cannot contain an array at line ",
                 line);
    }
  }

  static void ERROR_IF_NOT_REF(const variable_t &v, unsigned line) {
    if (!v.get_type().is_reference()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a reference at line ", line);
    }
  }

  static void ERROR_IF_NOT_INT(const variable_t &v, unsigned line) {
    if (!v.get_type().is_integer()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not an integer at line ", line);
    }
  }

  /**------------------ End helper method definitions ------------------**/

  /**------------- Begin MRU helper method definitions -------------------**/
  // perform a similiar work for cache.lines.constraint(cache.regs)
  // This one requires to be renaming because mem base dom used different base
  // vars
  base_abstract_domain_t cache_lines_regs_reduce() {
    base_abstract_domain_t res(m_cache_lines);
    auto csts = m_regs.to_linear_constraint_system();
    res += csts;
    return res;
  }

  // get all used Crab IR vars by region dom
  void vars(variable_vector_t &used_vars) const {
    used_vars.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      used_vars.push_back(kv.first);
    }
  }

  void region_vars_in_mem(variable_vector_t &used_vars) const {
    used_vars.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      if (kv.first.get_type().is_region()) {
        used_vars.push_back(kv.first);
      }
    }
  }

  void region_vars_in_cache(variable_vector_t &used_vars) const {
    used_vars.reserve(m_cache_var_map.size());
    for (auto &kv : m_cache_var_map) {
      if (kv.first.get_type().is_region()) {
        used_vars.push_back(kv.first);
      }
    }
  }

  void vars_in_cache(variable_vector_t &used_vars) const {
    used_vars.reserve(m_cache_var_map.size());
    for (auto &kv : m_cache_var_map) {
      used_vars.push_back(kv.first);
    }
  }

  variable_vector_t vec_minus(variable_vector_t left, variable_vector_t right) {
    using variable_map_t = ikos::patricia_tree<variable_t, variable_t>;
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

  // NOTE: For debug use
  void print_base_vec(const base_variable_vector_t &vec) {
    for (auto elem : vec)
      crab::outs() << elem << " ";
    crab::outs() << "\n";
  }

  // NOTE: For debug use
  void print_var_map() {
    crab::outs() << "Cache varmap\t\n";
    for (auto &kv : m_var_map)
      crab::outs() << kv.first << " ";
    crab::outs() << "\n";
  }

  // NOTE: For debug use
  void print_vec(const variable_vector_t &vec) {
    for (auto elem : vec)
      crab::outs() << elem << " ";
    crab::outs() << "\n";
  }

  // TODO: need to implement this on next iteration
  void invalidate_cache() {}

  // fresh cache as: { cache lines: bot, regs: top, addrs: addrs -= C_base }
  // for addrs, we keep all equalities of base address for each reference
  // also need to assign cache used flag C_used as false
  void make_fresh_cache() {
    m_cache_lines.set_to_bottom();
    m_regs.set_to_top();
    m_addrs.set_to_top();
    if (m_base) {
      forget(*m_base);
    }
    if (m_used) {
      auto false_cst = base_linear_constraint_t::get_false();
      base_variable_t c_used = get_or_insert_gvars(*m_used).get_var();
      m_mem.assign_bool_cst(c_used, false_cst);
    }
  }

  // commit cache contents into memory
  // require to fresh cache after commit.
  void commit_cache() {
    fold();
    make_fresh_cache();
  }

  void update_cache(variable_t rgn, variable_t base) {
    // clear cache var map
    m_cache_var_map.clear();
    expand(rgn);  // perform cache_line expand
    if (m_base) { // C_base == base
      base_variable_t c_base = get_or_insert_gvars(*m_base).get_var();
      base_variable_t base_in_base = get_or_insert_gvars(base).get_var();
      m_addrs.assign(c_base, base_in_base);
    }
    if (m_used) { // C_used = true
      auto true_cst = base_linear_constraint_t::get_true();
      base_variable_t c_used = get_or_insert_gvars(*m_used).get_var();
      m_mem.assign_bool_cst(c_used, true_cst);
    }
    if (m_dirty) { // C_dirty = false
      auto false_cst = base_linear_constraint_t::get_false();
      base_variable_t c_dirty = get_or_insert_gvars(*m_dirty).get_var();
      m_mem.assign_bool_cst(c_dirty, false_cst);
    }
  }
  /**------------- End MRU helper method definitions -------------------**/

public:
  /**------------------ Begin domain API definitions ------------------**/
  mru_region_domain_t make_top() const override {
    return mru_region_domain_t(true);
  }

  mru_region_domain_t make_bottom() const override {
    return mru_region_domain_t(false);
  }

  void set_to_top() override {
    mru_region_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    mru_region_domain_t abs(false);
    std::swap(*this, abs);
  }

  mru_region_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  mru_region_domain(const mru_region_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_alloc(o.m_alloc), m_var_map(o.m_var_map),
        m_rev_var_map(o.m_rev_var_map), m_cache_var_map(o.m_cache_var_map),
        m_mem(o.m_mem), m_cache_lines(o.m_cache_lines), m_regs(o.m_regs),
        m_addrs(o.m_addrs), m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom), m_rgn_type_dom(o.m_rgn_type_dom),
        m_alloc_site_dom(o.m_alloc_site_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom), m_tag_env(o.m_tag_env),
        m_obj_rgn_map(o.m_obj_rgn_map), m_rev_obj_rgn_map(o.m_rev_obj_rgn_map) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  mru_region_domain(mru_region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_cache_var_map(std::move(o.m_cache_var_map)),
        m_mem(std::move(o.m_mem)), m_cache_lines(std::move(o.m_cache_lines)),
        m_regs(std::move(o.m_regs)), m_addrs(std::move(o.m_addrs)),
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
        m_rgn_type_dom(std::move(o.m_rgn_type_dom)),
        m_alloc_site_dom(std::move(o.m_alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(o.m_rgn_dealloc_dom)),
        m_tag_env(std::move(o.m_tag_env)),
        m_obj_rgn_map(std::move(o.m_obj_rgn_map)),
        m_rev_obj_rgn_map(std::move(o.m_rev_obj_rgn_map)) {}

  mru_region_domain_t &operator=(const mru_region_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_cache_var_map = o.m_cache_var_map;
      m_mem = o.m_mem;
      m_cache_lines = o.m_cache_lines;
      m_regs = o.m_regs;
      m_addrs = o.m_addrs;
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_rgn_init_dom = o.m_rgn_init_dom;
      m_rgn_type_dom = o.m_rgn_type_dom;
      m_alloc_site_dom = o.m_alloc_site_dom;
      m_rgn_dealloc_dom = o.m_rgn_dealloc_dom;
      m_tag_env = o.m_tag_env;
      m_obj_rgn_map = o.m_obj_rgn_map;
      m_rev_obj_rgn_map = o.m_rev_obj_rgn_map;
    }
    return *this;
  }

  mru_region_domain_t &operator=(mru_region_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = std::move(o.m_alloc);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_cache_var_map = std::move(o.m_cache_var_map);
      m_mem = std::move(o.m_mem);
      m_cache_lines = std::move(o.m_cache_lines);
      m_regs = std::move(o.m_regs);
      m_addrs = std::move(o.m_addrs);
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_rgn_init_dom = std::move(o.m_rgn_init_dom);
      m_rgn_type_dom = std::move(o.m_rgn_type_dom);
      m_alloc_site_dom = std::move(o.m_alloc_site_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_tag_env = std::move(o.m_tag_env);
      m_obj_rgn_map = std::move(o.m_obj_rgn_map);
      m_rev_obj_rgn_map = std::move(o.m_rev_obj_rgn_map);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    bool res = !is_bottom() && m_mem.is_top() && m_cache_lines.is_top() &&
               m_regs.is_top() && m_addrs.is_top() &&
               m_rgn_counting_dom.is_top();
    if (crab_domain_params_man::get().region_allocation_sites()) {
      res = res && m_alloc_site_dom.is_top();
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      res = res && m_rgn_dealloc_dom.is_top();
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      res = res && m_tag_env.is_top();
    }
    return res;
  }

  bool operator<=(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    CRAB_LOG("mru-region-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                            << "\n\t" << o << "\n";);
    if (is_bottom() || o.is_top()) { // this indeed <= o
      CRAB_LOG("mru-region-leq", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) { // this indeed > o
      CRAB_LOG("mru-region-leq", crab::outs() << "Result2=0\n";);
      return false;
    }

    if (!(m_rgn_counting_dom <= o.m_rgn_counting_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result3=0\n";);
      return false;
    }

    if (!(m_rgn_init_dom <= o.m_rgn_init_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result4=0\n";);
      return false;
    }

    if (!(m_rgn_type_dom <= o.m_rgn_type_dom)) {
      CRAB_LOG("mru-region-leq", crab::outs() << "Result5=0\n";);
      return false;
    }

    if (crab_domain_params_man::get().region_allocation_sites()) {
      if (!(m_alloc_site_dom <= o.m_alloc_site_dom)) {
        CRAB_LOG("mru-region-leq", crab::outs() << "Result6=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_deallocation()) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("mru-region-leq", crab::outs() << "Result7=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_tag_analysis()) {
      if (!(m_tag_env <= o.m_tag_env)) {
        CRAB_LOG("mru-region-leq", crab::outs() << "Result8=0\n";);
        return false;
      }
    }

    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);

    base_abstract_domain_t left_mem(m_mem);
    base_abstract_domain_t right_mem(o.m_mem);
    base_abstract_domain_t left_cache_lines(m_cache_lines);
    base_abstract_domain_t right_cache_lines(o.m_cache_lines);
    base_uf_domain_t left_regs(m_regs);
    base_uf_domain_t right_regs(o.m_regs);
    base_uf_domain_t left_addrs(m_addrs);
    base_uf_domain_t right_addrs(o.m_addrs);

    // perform the common renaming
    base_variable_vector_t left_vars, right_vars, out_vars;
    // reserve vec space by upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = o.m_var_map.find(v);
      if (it != o.m_var_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          CRAB_ERROR(domain_name(), "::operator<= ", " dynamic type of ", v,
                     " must be the same in both left and right operands");
        }
        ghost_variables_t out_gvars(
            ghost_variables_t::create(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
      }
    }

    left_mem.project(left_vars);
    left_mem.rename(left_vars, out_vars);
    right_mem.project(right_vars);
    right_mem.rename(right_vars, out_vars);

    left_cache_lines.project(left_vars);
    left_cache_lines.rename(left_vars, out_vars);
    right_cache_lines.project(right_vars);
    right_cache_lines.rename(right_vars, out_vars);

    left_regs.project(left_vars);
    left_regs.rename(left_vars, out_vars);
    right_regs.project(right_vars);
    right_regs.rename(right_vars, out_vars);

    left_addrs.project(left_vars);
    left_addrs.rename(left_vars, out_vars);
    right_addrs.project(right_vars);
    right_addrs.rename(right_vars, out_vars);

    CRAB_LOG("mru-region-leq", crab::outs()
                                   << "Inclusion test (after renaming):\n\t"
                                   << "cache line:" << left_cache_lines
                                   << "\n\t" << right_cache_lines << "\n"
                                   << "reg" << left_regs << "\n\t" << right_regs
                                   << "\n";);

    bool res = left_mem <= right_mem && left_cache_lines <= right_cache_lines &&
               left_regs <= right_regs && left_addrs <= right_addrs;
    CRAB_LOG("mru-region-leq", crab::outs() << "Result9=" << res << "\n";);
    return res;
  }

  void operator|=(const mru_region_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) { // bot | o = o
      *this = o;
      return;
    } else if (o.is_bottom()) { // this | bot = this
      return;
    } else if (is_top() || o.is_top()) { // top | o or this | top, return top
      set_to_top();
      return;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Join " << *this << " and " << o << "=");

    variable_vector_t left_regions, right_regions;
    std::vector<ghost_variables_t> regions_left_base_vars,
        regions_right_base_vars;

    bool refine_left = false;
    bool refine_right = false;
    if (crab_domain_params_man::get().region_refine_uninitialized_regions()) {
      refine_left = can_propagate_initialized_regions(
          m_rgn_init_dom, o.m_var_map, left_regions, regions_right_base_vars);
      refine_right = can_propagate_initialized_regions(
          o.m_rgn_init_dom, m_var_map, right_regions, regions_left_base_vars);
    }

    // The code is complicated to achieve a zero-cost abstraction. If
    // we cannot improve invariants on the right operand we avoid
    // making a copy of it.
    if (refine_left && !refine_right) {
      refine_regions(left_regions, regions_right_base_vars, m_var_map,
                     m_rev_var_map, m_alloc, m_mem, o.m_mem, o.m_alloc);
      do_join(o);
    } else if (!refine_left && refine_right) {
      mru_region_domain_t right(o);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_mem, m_mem,
                     m_alloc);
      do_join(right);
    } else if (refine_left && refine_right) {
      mru_region_domain_t right(o);
      refine_regions(left_regions, regions_right_base_vars, m_var_map,
                     m_rev_var_map, m_alloc, m_mem, right.m_mem, right.m_alloc);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_mem, m_mem,
                     m_alloc);
      do_join(right);
    } else {
      do_join(o);
    }
  }

  mru_region_domain_t operator|(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) { // bot | o = o
      return o;
    } else if (o.is_bottom()) { // this | bot = this
      return *this;
    } else if (is_top() || o.is_top()) { // top | o or this | top, return top
      mru_region_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Join " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 | v2; };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 | v2; };

    variable_vector_t left_regions, right_regions;
    std::vector<ghost_variables_t> regions_left_base_vars,
        regions_right_base_vars;

    bool refine_left = false;
    bool refine_right = false;
    if (crab_domain_params_man::get().region_refine_uninitialized_regions()) {
      refine_left = can_propagate_initialized_regions(
          m_rgn_init_dom, o.m_var_map, left_regions, regions_right_base_vars);
      refine_right = can_propagate_initialized_regions(
          o.m_rgn_init_dom, m_var_map, right_regions, regions_left_base_vars);
    }

    // The code is complicated to achieve a zero-cost abstraction. If
    // we cannot improve invariants then we try to avoid making copies
    // of left and/or right operands.

    if (refine_left && !refine_right) {
      // Refine left by propagating information from right's regions
      mru_region_domain_t left(*this);
      refine_regions(left_regions, regions_right_base_vars, left.m_var_map,
                     left.m_rev_var_map, left.m_alloc, left.m_mem, o.m_mem,
                     o.m_alloc);
      mru_region_domain_t res(std::move(do_join_or_widening(
          left, o, true /*is join*/, base_dom_op, base_uf_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (!refine_left && refine_right) {
      // Refine right by propagating information from left's regions
      mru_region_domain_t right(o);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_mem, m_mem,
                     m_alloc);
      mru_region_domain_t res(std::move(do_join_or_widening(
          *this, right, true /*is join*/, base_dom_op, base_uf_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (refine_left && refine_right) {
      // Refine both left and right
      mru_region_domain_t left(*this);
      mru_region_domain_t right(o);
      refine_regions(left_regions, regions_right_base_vars, left.m_var_map,
                     left.m_rev_var_map, left.m_alloc, left.m_mem, right.m_mem,
                     right.m_alloc);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_mem,
                     left.m_mem, left.m_alloc);
      mru_region_domain_t res(std::move(do_join_or_widening(
          left, right, true /*is join*/, base_dom_op, base_uf_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else {
      mru_region_domain_t res(std::move(do_join_or_widening(
          *this, o, true /*is join*/, base_dom_op, base_uf_dom_op)));

      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;
    }
  }

  mru_region_domain_t operator&(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 & v2; };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 & v2; };

    mru_region_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, true /*is meet*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t operator||(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 || v2;
    };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 || v2; };

    mru_region_domain_t res(std::move(do_join_or_widening(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t widening_thresholds(
      const mru_region_domain_t &o,
      const iterators::thresholds<number_t> &thresholds) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("mru-region", crab::outs() << "Widening with threshold " << *this
                                        << " and " << o << "=");

    auto base_dom_op = [&thresholds](const base_abstract_domain_t &v1,
                                     const base_abstract_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };
    auto base_uf_dom_op = [&thresholds](const base_uf_domain_t &v1,
                                        const base_uf_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    mru_region_domain_t res(std::move(do_join_or_widening(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  mru_region_domain_t operator&&(const mru_region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("mru-region", crab::outs()
                               << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 && v2;
    };
    auto base_uf_dom_op = [](const base_uf_domain_t &v1,
                             const base_uf_domain_t &v2) { return v1 && v2; };

    mru_region_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, false /*is widen*/, base_dom_op, base_uf_dom_op)));

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    return res;
  }

  // x := e operates on m_mem
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_mem.assign(get_or_insert_gvars(x).get_var(), b_e);
    }
  }

  // x := e operates on m_cache_lines
  void cache_lines_assign(const variable_t &x, const linear_expression_t &e) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    // TODO: must check x and e are region-based assignment
    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      const ghost_variables_t &new_gvars = get_or_insert_gvars(x);
      m_cache_var_map.insert({x, new_gvars});
      m_cache_lines.assign(get_or_insert_gvars(x).get_var(), b_e);
    }
  }

  // x := e operates on m_regs
  void regs_assign(const variable_t &x, const linear_expression_t &e) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    // TODO: must check x is a region variable and e is a scalar variable
    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_regs.assign(get_or_insert_gvars(x).get_var(), b_e);
    }
  }

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      auto b_csts = rename_linear_cst_sys(csts);
      m_mem += b_csts;
      m_is_bottom = m_mem.is_bottom();
    }
  }

  void add_cons_into_cache_lines(const linear_constraint_system_t &csts) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      auto b_csts = rename_linear_cst_sys(csts);
      m_cache_lines += b_csts;
    }
  }

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = rename_linear_cst(rhs);
      m_mem.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);

      if (crab_domain_params_man::get().region_tag_analysis()) {
        merge_tags(lhs, rhs.variables());
      }
    }
  }

  // FIXME: The followings are UNDEFINED METHODS
  boolean_value is_null_ref(const variable_t &ref) override {
    return boolean_value::bottom();
  }

  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {}
  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {}

  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {}
  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {}
  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {}

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {}

  /********************** Boolean operations **********************/
  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {}
  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {}
  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {}
  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {}

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {}

  /********************** Array operations **********************/
  // make a fresh array with contents a[j] initialized to val such that
  // j \in [lb_idx,ub_idx] and j % elem_size == val.
  // elem_size is in bytes.
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  // lhs := a[i] where elem_size is in bytes
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {}
  // a[i] := val where elem_size is in bytes
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &val,
                   bool is_strong_update) override {}
  // forall i<=k<j and k % elem_size == 0 :: a[k] := val.
  // elem_size is in bytes
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {}
  // forall i :: a[i] := b[i]
  void array_assign(const variable_t &a, const variable_t &b) override {}

  /********************* Regions and reference operations *********************/
  // Initialize a region
  void region_init(const variable_t &reg) override {}
  // Make a copy of a region
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {}
  // Cast between regions of different types
  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {}
  // Create a new reference ref associated with as within region reg
  void ref_make(const variable_t &ref, const variable_t &reg,
                /* size of the allocation in bytes */
                const variable_or_constant_t &size,
                /* identifier for the allocation site */
                const allocation_site &as) override {}
  // Remove a reference ref within region reg
  void ref_free(const variable_t &reg, const variable_t &ref) override {}
  // Read the content of reference ref within reg. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {}
  // Write the content of val to the address pointed by ref in region
  // reg.
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {}
  // Create a new reference ref2 to region reg2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {}
  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const variable_t &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {}
  // Treat region as an array and perform an array store.
  void ref_store_to_array(const variable_t &ref, const variable_t &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {}
  // Add constraints between references
  void ref_assume(const reference_constraint_t &cst) override {}
  // Convert a reference to an integer variable
  void ref_to_int(const variable_t &reg, const variable_t &ref,
                  const variable_t &int_var) override {}
  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref) override {}
  // if (cond) ref_gep(ref1, rgn1, lhs_ref, lhs_rgn, 0) else
  //           ref_gep(ref2, rgn2, lhs_ref, lhs_rgn, 0)
  // cond is a boolean variable
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {}

  /********************** Backward numerical operations **********************/
  // x = y op z
  // Substitute x with y op z in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const mru_region_domain_t &invariant) override {}
  // x = y op k
  // Substitute x with y op k in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const mru_region_domain_t &invariant) override {}
  // x = e
  // Substitute x with e in the abstract value
  // The result is meet with invariant.
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const mru_region_domain_t &invariant) override {}

  /********************** Backward boolean operations **********************/
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const mru_region_domain_t &invariant) override {
  }
  void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const mru_region_domain_t &invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const mru_region_domain_t &invariant) override {
  }
  void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const mru_region_domain_t &invariant) override {}

  /********************** Backward array operations **********************/
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const mru_region_domain_t &invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const mru_region_domain_t &invariant) override {}
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const mru_region_domain_t &invariant) override {}
  void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const mru_region_domain_t &invariant) override {}
  void backward_array_assign(const variable_t &a, const variable_t &b,
                             const mru_region_domain_t &invariant) override {}

  /********************** Miscellaneous operations **********************/
  // Forget v
  void operator-=(const variable_t &v) override {}

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {
    return interval_t::bottom();
  }

  // Normalize the abstract domain if such notion exists.
  void normalize() override {}

  // Reduce the size of the abstract domain representation.
  void minimize() override {}

  // Make a new copy of var without relating var with new_var
  void expand(const variable_t &var, const variable_t &new_var) override {}

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const mru_region_domain_t &invariant) override {}

  // FIXME: The above methods are UNDEFINED METHODS

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    for (auto &v : variables) {
      if (crab_domain_params_man::get().region_allocation_sites()) {
        if (v.get_type().is_reference() || v.get_type().is_region()) {
          m_alloc_site_dom -= v;
        }
      }
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;
        if (is_tracked_unknown_region(v)) {
          m_rgn_type_dom -= v;
        }
        if (crab_domain_params_man::get().region_deallocation()) {
          m_rgn_dealloc_dom.forget(v);
        }
      }
      if (crab_domain_params_man::get().region_tag_analysis()) {
        m_tag_env -= v;
      }
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        it->second.add(base_vars);
        it->second.remove_rev_varmap(m_rev_var_map);
        m_var_map.erase(it);
      }
    }
    // forget variable into each subdomain
    m_mem.forget(base_vars);
    m_cache_lines.forget(base_vars);
    m_regs.forget(base_vars);
    m_addrs.forget(base_vars);
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    variable_vector_t sorted_variables(variables.begin(), variables.end());
    std::sort(sorted_variables.begin(), sorted_variables.end());

    // -- project in the base domain
    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    for (auto const &v : variables) {
      if (v.get_type().is_region() && !is_tracked_region(v)) {
        continue;
      }
      get_or_insert_gvars(v).add(base_vars);
    }

    // project mem domain
    m_mem.project(base_vars);
    // project cache_lines domain
    m_cache_lines.project(base_vars);
    // project regs domain
    m_regs.project(base_vars);
    // project addrs domain
    m_addrs.project(base_vars);

    // -- update m_var_map and m_rev_var_map
    std::vector<variable_t> var_map_to_remove;
    var_map_to_remove.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      if (!std::binary_search(sorted_variables.begin(), sorted_variables.end(),
                              kv.first)) {
        var_map_to_remove.push_back(kv.first);
        kv.second.remove_rev_varmap(m_rev_var_map);
      }
    }
    for (auto &v : var_map_to_remove) {
      m_var_map.erase(v);
      m_cache_var_map.erase(v);
    }

    m_rgn_counting_dom.project(sorted_variables);
    m_rgn_init_dom.project(sorted_variables);
    if (!crab_domain_params_man::get().region_skip_unknown_regions()) {
      m_rgn_type_dom.project(sorted_variables);
    }
    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom.project(sorted_variables);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.project(sorted_variables);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.project(sorted_variables);
    }
  }

  // NOTE: region might init without sorting but rename method required
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_bottom() || is_top()) {
      return;
    }

    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename different lengths");
    }

    // update var_map and rev_var_map
    // extract base_variable_vector for from and to
    base_variable_vector_t old_base_vars, new_base_vars;
    for (unsigned i = 0; i < from.size(); ++i) {
      variable_t old_var = from[i];
      auto search = m_var_map.find(old_var);
      if (search == m_var_map.end()) { // NOT FOUND
        continue;
      }
      ghost_variables_t old_ghost_vars = search->second;
      // append base_var, [offset, size] into old_base_vars
      old_ghost_vars.add(old_base_vars);
      // remove from rev_var_map
      old_ghost_vars.remove_rev_varmap(m_rev_var_map);
      // remove from var_map
      m_var_map.erase(old_var);

      variable_t new_var = to[i];
      // create/get new ghost variable for new_var
      const ghost_variables_t &new_gvars = get_or_insert_gvars(new_var);
      // append base_var, [offset, size] into new_base_vars
      new_gvars.add(new_base_vars);
    }

    // rename mem domain
    m_mem.rename(old_base_vars, new_base_vars);
    // rename cache_lines domain
    m_cache_lines.rename(old_base_vars, new_base_vars);
    // rename regs domain
    m_regs.rename(old_base_vars, new_base_vars);
    // rename addrs domain
    m_addrs.rename(old_base_vars, new_base_vars);

    // rename the rest of environments
    m_rgn_counting_dom.rename(from, to);
    m_rgn_init_dom.rename(from, to);
    if (!crab_domain_params_man::get().region_skip_unknown_regions()) {
      m_rgn_type_dom.rename(from, to);
    }
    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom.rename(from, to);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.rename(from, to);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.rename(from, to);
    }
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    //=================================================================//
    //       Special intrinsics supported by the region domain
    //=================================================================//
    // ---DSA region analysis---
    //      This analysis indicates which regions might belong to the
    //      same memory object. The intrinstic is added only if object
    //      has more than one field.
    // TODO: need to support other analysis
    auto error_if_not_rgn = [&name](const variable_t &var) {
      if (!var.get_type().is_region()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var,
                   " should be a region");
      }
    };

    if (is_bottom()) {
      return;
    }

    if (name == "regions_from_memory_object") {
      // why pass to the outputs vector?
      // https://github.com/seahorn/clam/blob/ee2f50723a40377864fee7618b90df92d1413a58/lib/Clam/CfgBuilder.cc#L4501
      // pass region varaibles into cache_map
      variable_t representative = outputs[0];
      for (int i = 0, sz = outputs.size(); i < sz; ++i) {
        error_if_not_rgn(outputs[i]);
        if (outputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        // update object region map
        if (i == 0) {
          m_obj_rgn_map.insert({outputs[i], outputs});
        }
        m_rev_obj_rgn_map.insert({outputs[i], representative});
      }
    }
  }

  // Convert the abstract state into a conjunction of linear constraints.
  // Only integer or boolean variables in memory are exposed.
  //
  // Variables that shadow memory regions and reference variables are
  // ignored.
  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      return linear_constraint_t::get_false();
    } else if (is_top()) {
      return linear_constraint_t::get_true();
    } else {
      linear_constraint_system_t out_csts;
      const bool ignore_references = true;
      for (base_linear_constraint_t cst : m_mem.to_linear_constraint_system()) {
        if (boost::optional<linear_constraint_t> out_cst =
                rev_rename_linear_cst(cst, ignore_references)) {
          out_csts += *(out_cst);
        }
      }
      return out_csts;
    }
  }

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR("mru_region_domain::to_disjunctive_linear_constraint_system not "
               "implemented");
  }

  std::string domain_name() const override {
    return "MRURegionDomain(" + m_cache_lines.domain_name() + ", " +
           m_regs.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("mru-region-print", o << "Mem=\t" << m_mem << ",\n"
                                     << "( CacheLine=\t" << m_cache_lines
                                     << ",\n"
                                     << "Regs=\t" << m_regs << ",\n"
                                     << "Addrs=\t" << m_addrs << " )\n";);

      std::unordered_map<std::string, std::string> renaming_map;
      auto get_type_fn = [](const variable_t &v) -> variable_type {
        return variable_type(variable_type_kind::UNK_TYPE);
      };
      // Assigns a string name to each ghost variable
      ghost_variables_t::mk_renaming_map(m_rev_var_map, get_type_fn,
                                         renaming_map);
      m_alloc.add_renaming_map(renaming_map);
      o << "Mem=\t" << m_mem << ",\n"
        << "( CacheLine=\t" << m_cache_lines << ",\n"
        << "Regs=\t" << m_regs << ",\n"
        << "Addrs=\t" << m_addrs << " )\n";
      m_alloc.clear_renaming_map();
    }
  }
  /**------------- End domain API definitions -------------------**/

  /**------------- Begin MRU API definitions -------------------**/

  // The expand performs a copy of region vars from mem into cache_lines
  //  i.e. by given a most recently used object o, we copy all regions with o
  //  in mem into the cache:
  // For instance, o used V1, V2:
  //    mem: { 1 <= V1 <= 2; V1 < V2; ... }
  // The cache line will be:
  //    cache_lines: { 1 <= V1 <= 2; V1 < V2; }
  void expand(const variable_t &rgn) {
    variable_vector_t Vs_for_new_obj;
    auto it_var = m_rev_obj_rgn_map.find(rgn);
    if (it_var == m_rev_obj_rgn_map.end()) {
      // this should not happen
      CRAB_ERROR("mru_region_domain::expand the mru object should be found");
    }
    variable_t representative = it_var->second;
    auto it_rev_vars = m_obj_rgn_map.find(representative);
    if (it_rev_vars == m_obj_rgn_map.end()) {
      // this should not happen
      CRAB_ERROR("mru_region_domain::expand the mru object should be found");
    }
    Vs_for_new_obj = it_rev_vars->second;

    // expand does not need to perform base variable renaming
    mru_region_domain_t tmp = mru_region_domain_t(*this);
    tmp.project(Vs_for_new_obj);

    // perform expand
    m_cache_lines |= tmp.m_mem;

    CRAB_LOG("mru-region", crab::outs() << *this << "\n");
  }

  // The cache is used by given an allocated object.
  //  The fold operation performs smash all regions (fields) into the summary
  //  object in mem. That is, for each region used in cache line, join with
  //  region in mem.
  //  In addition, we also update regs used in cache via reduction.
  void fold() {
    variable_vector_t Vs;
    region_vars_in_mem(Vs);

    // Put together cache lines and regs
    // Step 1.
    //  Perform regs: { r == Cj ... }
    //  constrain with cache_lines: { Ci <= Cj <= Ck ... }, i != j != k
    // the reduced cache_lines will be:
    //    { Ci <= Cj <= Ck; Ci <= r <= Ck; ...  }
    base_abstract_domain_t reduced = cache_lines_regs_reduce();

    // Step 2.
    //  Our design is joining cache line with memory.
    //  We need the projection before join because we only fold regions & regs
    //  We do not affect other properties in mem.
    variable_vector_t vars_cache;
    vars_in_cache(vars_cache); // vars for cache lines + registers
    mru_region_domain_t A = mru_region_domain_t(*this);
    A.project(vars_cache); // memory object represented in the cache
    A.m_mem |= reduced;    // join cache with mem,
                           // this step might lost other variables

    // Step 3.
    // So far we performed the fold operation, but the rest properties need to
    // put it back. This required we project all variables used in cache
    // Finally, we perform meet to recover.
    mru_region_domain_t B = mru_region_domain_t(*this);
    variable_vector_t mem_vars;
    vars(mem_vars); // All variables in mem
    variable_vector_t vars_cache_lines;
    region_vars_in_cache(vars_cache_lines); // All region variables in cache
    B.project(vec_minus(mem_vars,
                        vars_cache_lines)); // get all properties except those
                                            // relevant to the cache lines
    mru_region_domain_t res = A & B; // put back the rest of the properties

    CRAB_LOG("mru-region", crab::outs() << res << "\n");
    std::swap(*this, res);
  }

  /**------------- End MRU API definitions -------------------**/

}; // end class mru_region_domain

template <typename Params>
struct abstract_domain_traits<mru_region_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab