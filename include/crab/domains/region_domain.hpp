#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/tag.hpp>
#include <crab/types/varname_factory.hpp>

#include <algorithm>
#include <functional>
#include <unordered_map>

namespace crab {
namespace domains {
namespace region_domain_impl {
template <class Number, class VariableName, class BaseAbsDom>
class Params {
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
} // end namespace region_domain_impl

////////////////////////////////////////////////////////////////////////
// Abstract domain for regions and references.
////////////////////////////////////////////////////////////////////////
template <typename Params>
class region_domain final : public abstract_domain_api<region_domain<Params>> {
  using region_domain_t = region_domain<Params>;
  using abstract_domain_t = abstract_domain_api<region_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
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
  
  // Environment domains: map regions to finite domain
  using rgn_counting_env_t = ikos::separate_domain<variable_t, small_range>;
  using rgn_bool_env_t = ikos::separate_domain<variable_t, boolean_value>;
  // Union-find where equivalence classes are attached to boolean values
  using rgn_dealloc_t = union_find_domain<variable_t, boolean_value>;
  // Map from variable to base variable
  using var_map_t = std::unordered_map<variable_t, base_variable_t>;
  using rev_var_map_t = std::unordered_map<base_variable_t, variable_t>;
  // Map each reference or region variable to a set of allocation sites
  using alloc_site_env_t = separate_discrete_domain<variable_t, allocation_site>;
  using allocation_sites = typename alloc_site_env_t::value_type;
  // Tags
  class tag_t : public indexable {
    ikos::index_t m_id;
  public:
    tag_t(number_t n): m_id(0) {
      if (n < 0)
      { CRAB_ERROR("Cannot use negative numbers for tags");}
      if (!n.fits_int64())
      { CRAB_ERROR("Too large value for a tag"); }
      m_id = (int64_t) n;
    }
    bool operator<(const tag_t &as) const { return m_id < as.m_id; }
    bool operator==(const tag_t &as) const { return m_id == as.m_id; }
    virtual ikos::index_t index() const override { return m_id; }
    void write(crab_os &o) const override { o << "TAG_" << m_id; }
    friend crab_os &operator<<(crab_os &o, const tag_t &as) {
      as.write(o);
      return o;
    }
  };
  // Map variables to sets of tags
  using tag_env_t = separate_discrete_domain<variable_t, tag_t>;
  using tag_set = typename tag_env_t::value_type;
  
  bool m_is_bottom; // special symbol for bottom

  /** Begin base domain **/
  // To create ghost variables for the base domain.
  base_varname_allocator_t m_alloc;
  // Map a variable_t to base domain's variable:
  //  - int/real/bool variable is mapped while preserving its type
  //  - ref variables are mapped to **integer** variables to ensure
  //    that the base domain is completely unaware of references.
  //  - regions of int/real/bool are mapped to int/real/bool
  //  - regions of ref are mapped to int
  //  - regions of array of int/bool are mapped to arrays of int/bool
  var_map_t m_var_map;
  rev_var_map_t m_rev_var_map;
  // The base abstract domain: all the heavy lifting is done here.
  // m_base_dom does not have any variable of reference type.
  base_abstract_domain_t m_base_dom;
  /** End base domain **/
  
  // Abstract domain to count how many addresses are in a region.
  // This allows us to decide when strong update is sound: only if
  // one address per region (i.e., singleton).
  rgn_counting_env_t m_rgn_counting_dom;
  // Whether some data might have been written to any address within
  // the region.
  rgn_bool_env_t m_rgn_init_dom;
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
  
  region_domain(base_varname_allocator_t &&alloc, var_map_t &&var_map,
                rev_var_map_t &&rev_var_map,
                rgn_counting_env_t &&rgn_counting_dom,
                rgn_bool_env_t &&rgn_init_dom, 
                base_abstract_domain_t &&base_dom,
		alloc_site_env_t &&alloc_site_dom,
		rgn_dealloc_t &&rgn_dealloc_dom,
		tag_env_t &&tag_env)
      : m_is_bottom(base_dom.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
	m_base_dom(std::move(base_dom)),	
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
	m_alloc_site_dom(std::move(alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
	m_tag_env(std::move(tag_env)) {}

  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;

  bool can_propagate_initialized_regions(
      const rgn_bool_env_t &left_init_dom, const var_map_t &right_varmap,
      variable_vector_t &regions,
      base_variable_vector_t &right_base_vars) const {
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
  // This step is not perfect because it keeps only non-relational
  // information about regions but no relational information between
  // regions and other variables.
  void refine_regions(const variable_vector_t &left_regions,
                      /* the base variable for each regions[i] on right_dom*/
                      const base_variable_vector_t &old_right_base_vars,
                      /* left operand */
                      var_map_t &left_varmap, rev_var_map_t &left_rev_varmap,
                      base_varname_allocator_t &left_alloc,
                      base_abstract_domain_t &left_dom,
                      /* right operand */
                      const base_abstract_domain_t &right_dom,
                      base_varname_allocator_t old_right_alloc) const {
    assert(left_regions.size() == old_right_base_vars.size());
    if (left_regions.empty())
      return;

    // This is needed to avoid variable clashing with the left operand
    base_varname_allocator_t right_alloc(left_alloc, old_right_alloc);

    /* Modify the left by creating a variable in the base domain for
       the region*/
    base_variable_vector_t left_base_vars, right_base_vars;
    left_base_vars.reserve(left_regions.size());
    for (unsigned i = 0, sz = left_regions.size(); i < sz; ++i) {
      left_base_vars.push_back(rename_var(left_regions[i], left_varmap,
                                          left_rev_varmap, left_alloc));
      right_base_vars.push_back(base_variable_t(
          right_alloc.next(), old_right_base_vars[i].get_type()));
    }

    /* Propagate invariants on the region from the right to the
       left */
    base_abstract_domain_t dom(right_dom);
    dom.project(old_right_base_vars);

    // Renaming in two steps to avoid variable clashing with the left
    // operand. The assumption is that old_right_base_vars and
    // left_base_vars might have common names.
    dom.rename(old_right_base_vars, right_base_vars);
    dom.rename(right_base_vars, left_base_vars);
    left_dom = left_dom & dom;
  }

  // Perform *this = join(*this, right)
  void do_join(const region_domain_t &right) {
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    alloc_site_env_t out_alloc_site_dom(m_alloc_site_dom |
					right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);

    base_varname_allocator_t out_alloc(m_alloc, right.m_alloc);
    base_abstract_domain_t right_dom(right.m_base_dom);
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;
    base_variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    // perform the common renamings
    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_var_map.insert({v, out_v});
        out_rev_var_map.insert({out_v, v});
      }
    }

    // See comments in do_join_or_widening
    m_base_dom.project(left_vars);
    m_base_dom.rename(left_vars, out_vars);
    right_dom.project(right_vars);
    right_dom.rename(right_vars, out_vars);

    m_base_dom |= right_dom;
    m_is_bottom = m_base_dom.is_bottom();
    std::swap(m_alloc, out_alloc);
    std::swap(m_rgn_counting_dom, out_rgn_counting_dom);
    std::swap(m_rgn_init_dom, out_rgn_init_dom);
    std::swap(m_alloc_site_dom, out_alloc_site_dom);    
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_tag_env, out_tag_env);

    CRAB_LOG("region", crab::outs() << *this << "\n");
  }

  region_domain_t do_join_or_widening(const region_domain_t &left,
                                      const region_domain_t &right,
                                      base_dom_binop_t base_dom_op) const {

    // rgn_counting_dom does not require common renaming (i.e., no
    // ghost variables).  The domain is finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    // rgn_init_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom | right.m_rgn_init_dom);

    // alloc_site_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom | right.m_alloc_site_dom);
    
    // rgn_dealloc_dom does not require common renaming (i.e., no
    // ghost variables).  The domain is finite
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    // tag_env does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    tag_env_t out_tag_env(left.m_tag_env | right.m_tag_env);
    
    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);
    base_abstract_domain_t left_dom(left.m_base_dom);
    base_abstract_domain_t right_dom(right.m_base_dom);
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;

    // -- Compute common renamings
    base_variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings = left.m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : left.m_var_map) {
      const variable_t &v = kv.first;
      auto it = right.m_var_map.find(v);
      if (it != right.m_var_map.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_var_map.insert({v, out_v});
        out_rev_var_map.insert({out_v, v});
      }
    }

    // JN: project might be necessary to avoid keeping variables that
    // exist in the base domain but they doen't exist on m_var_map.
    //
    // If such a variable exists only on either left_dom or right_dom
    // the join removes it.  However, if we have the same variable in
    // both left_dom and righ_dom then the join will preserve it.
    //
    left_dom.project(left_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.project(right_vars);
    right_dom.rename(right_vars, out_vars);

    // Final join or widening
    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_base_dom),
	std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
	std::move(out_tag_env));
    return res;
  }

  region_domain_t do_meet_or_narrowing(const region_domain_t &left,
                                       const region_domain_t &right,
                                       base_dom_binop_t base_dom_op) const {

    // these domains are finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom &
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom & right.m_rgn_init_dom);
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom &
					right.m_alloc_site_dom);    
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom &
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(left.m_tag_env & right.m_tag_env);
    
    // This shouldn't happen but just in case ...
    if (out_rgn_counting_dom.is_bottom() || out_rgn_init_dom.is_bottom() ||
        out_rgn_dealloc_dom.is_bottom() || out_alloc_site_dom.is_bottom()) {
      return make_bottom();
    }

    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);
    base_abstract_domain_t left_dom(left.m_base_dom);
    base_abstract_domain_t right_dom(right.m_base_dom);
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;

    base_variable_vector_t left_vars, right_vars, out_vars;
    base_variable_vector_t only_left_vars, only_left_out_vars;
    base_variable_vector_t only_right_vars, only_right_out_vars;
    // upper bound to avoid reallocations
    size_t left_renamings = left.m_var_map.size();
    size_t right_renamings = right.m_var_map.size();
    size_t num_renamings = left_renamings + right_renamings;
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
        // common renaming
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_var_map.insert({v, out_v});
        out_rev_var_map.insert({out_v, v});
      } else {
        // only on left
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        only_left_vars.push_back(kv.second);
        only_left_out_vars.push_back(out_v);
        out_var_map.insert({kv.first, out_v});
        out_rev_var_map.insert({out_v, kv.first});
      }
    }

    // add missing maps from right operand
    for (auto &kv : right.m_var_map) {
      const variable_t &v = kv.first;
      auto it = left.m_var_map.find(v);
      if (it == left.m_var_map.end()) {
        // only on right
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        only_right_vars.push_back(kv.second);
        only_right_out_vars.push_back(out_v);
        out_var_map.insert({kv.first, out_v});
        out_rev_var_map.insert({out_v, kv.first});
      }
    }

    // append only_left_vars and left_vars
    only_left_vars.insert(only_left_vars.end(), left_vars.begin(),
                          left_vars.end());
    // append only_left_out_vars end out_vars
    only_left_out_vars.insert(only_left_out_vars.end(), out_vars.begin(),
                              out_vars.end());
    // need to project before renaming
    left_dom.project(only_left_vars);
    left_dom.rename(only_left_vars, only_left_out_vars);

    // append only_right_vars and right_vars
    only_right_vars.insert(only_right_vars.end(), right_vars.begin(),
                           right_vars.end());
    // append only_right_out_vars end out_vars
    only_right_out_vars.insert(only_right_out_vars.end(), out_vars.begin(),
                               out_vars.end());
    // need to project before renaming
    right_dom.project(only_right_vars);
    right_dom.rename(only_right_vars, only_right_out_vars);

    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_base_dom),
	std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
	std::move(out_tag_env));

    return res;
  }

  // Create a fresh ghost variable in the base domain to shadow v
  base_variable_t make_base_variable(base_varname_allocator_t &var_allocator,
                                     const variable_t &v) const {
    if (v.get_type().is_reference()) {
      base_variable_t bv(var_allocator.next(), crab::INT_TYPE, 32);
      return bv;
    } else if (v.get_type().is_region()) {
      variable_type_kind ty;
      unsigned bitwidth = 0;
      auto vty = v.get_type();
      if (vty.is_bool_region()) {
        ty = BOOL_TYPE;
        bitwidth = 1;
      } else if (vty.is_integer_region()) {
        ty = INT_TYPE;
        bitwidth = vty.get_integer_region_bitwidth();
      } else if (vty.is_reference_region()) {
        ty = INT_TYPE;
        bitwidth = 32;
      } else if (vty.is_real_region()) {
        ty = REAL_TYPE;
      } else if (vty.is_bool_array_region()) {
        ty = ARR_BOOL_TYPE;
      } else if (vty.is_int_array_region()) {
        ty = ARR_INT_TYPE;
      } else if (vty.is_real_array_region()) {
        ty = ARR_REAL_TYPE;
      } else {
	assert(false);
        CRAB_ERROR(domain_name() + "::make_base_variable: unreachable");
      }
      base_variable_t bv(var_allocator.next(), ty, bitwidth);
      return bv;
    } else {
      base_variable_t bv(var_allocator.next(), v.get_type());
      return bv;
    }
  }

  // Get the ghost variable in the base domain that represents v
  boost::optional<base_variable_t> get_var(const variable_t &v) const {
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      return boost::optional<base_variable_t>();
    }
  }

  /*
   * To rename variables so that they have the same variable's type
   * used by the base abstract domain.
   */

  // Rename variable
  const base_variable_t &rename_var(const variable_t &v) {
    return rename_var(v, m_var_map, m_rev_var_map, m_alloc);
  }

  const base_variable_t &rename_var(const variable_t &v, var_map_t &varmap,
                                    rev_var_map_t &rev_varmap,
                                    base_varname_allocator_t &alloc) const {
    auto it = varmap.find(v);
    if (it != varmap.end()) {
      return it->second;
    } else {
      base_variable_t bv = make_base_variable(alloc, v);
      auto res = varmap.insert({v, bv});
      if (res.second) {
        rev_varmap.insert({bv, v});
      }
      return res.first->second;
    }
  }

  // used only for pretty-printing and to_linear_constraint_system()
  boost::optional<variable_t> rev_rename_var(const base_variable_t &v,
                                             bool ignore_references) const {
    auto it = m_rev_var_map.find(v);
    if (it != m_rev_var_map.end()) {
      if (!ignore_references || !it->second.get_type().is_reference()) {
        return it->second;
      }
    }
    return boost::optional<variable_t>();
  }

  // Rename linear expression
  base_linear_expression_t rename_linear_expr(const linear_expression_t &e) {
    // TODO caching
    base_linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      const base_variable_t &nv = rename_var(v);
      out = out + (coef * nv);
    }
    return out;
  }

  // used only for pretty-printing and to_linear_constraint_system()
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

  // Rename linear constraint
  base_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst) {
    // TODO caching
    if (cst.is_inequality() || cst.is_strict_inequality()) {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression()),
          (typename base_linear_constraint_t::kind_t)cst.kind(),
          cst.is_signed());

    } else {
      return base_linear_constraint_t(
          rename_linear_expr(cst.expression()),
          (typename base_linear_constraint_t::kind_t)cst.kind());
    }
  }

  // used only for pretty-printing and to_linear_constraint_system()
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
  rename_linear_cst_sys(const linear_constraint_system_t &csts) {
    // TODO caching
    base_linear_constraint_system_t out;
    for (auto const &cst : csts) {
      out += rename_linear_cst(cst);
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
        base_variable_t x = rename_var(ref_cst.lhs());
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
        base_variable_t x = rename_var(ref_cst.lhs());
        base_variable_t y = rename_var(ref_cst.rhs());
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

  template<class RangeVars>
  void merge_tags(const variable_t &x, RangeVars vars) {
    tag_set tags = tag_set::bottom();
    for (auto const &v: vars) {
      tags = tags | m_tag_env[v];
    }
    m_tag_env.set(x, tags);
  }

  static void ERROR_IF_NOT_REGION(const variable_t& v, unsigned line) {
    if (!v.get_type().is_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a region at line ", line);
    }
  }

  static void ERROR_IF_ARRAY_REGION(const variable_t& v, unsigned line) {
    if (v.get_type().is_array_region()) {
      CRAB_ERROR(v, ":", v.get_type(), " cannot contain an array at line ", line);
    }
  }
  
  static void ERROR_IF_NOT_REF(const variable_t& v, unsigned line) {
    if (!v.get_type().is_reference()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not a reference at line ", line);
    }
  }

  static void ERROR_IF_NOT_INT(const variable_t& v, unsigned line) {
    if (!v.get_type().is_integer()) {
      CRAB_ERROR(v, ":", v.get_type(), " is not an integer at line ", line);
    }
  }
  
public:
  region_domain_t make_top() const override { return region_domain_t(true); }

  region_domain_t make_bottom() const override {
    return region_domain_t(false);
  }

  void set_to_top() override {
    region_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    region_domain_t abs(false);
    std::swap(*this, abs);
  }

  region_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  region_domain(const region_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_alloc(o.m_alloc), m_var_map(o.m_var_map),
        m_rev_var_map(o.m_rev_var_map),
	m_base_dom(o.m_base_dom),	
        m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom),
	m_alloc_site_dom(o.m_alloc_site_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom),
	m_tag_env(o.m_tag_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  region_domain(region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
	m_base_dom(std::move(o.m_base_dom)),	
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
	m_alloc_site_dom(std::move(o.m_alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(o.m_rgn_dealloc_dom)),
	m_tag_env(std::move(o.m_tag_env)) {}

  region_domain_t &operator=(const region_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_base_dom = o.m_base_dom;      
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_rgn_init_dom = o.m_rgn_init_dom;
      m_alloc_site_dom = o.m_alloc_site_dom;
      m_rgn_dealloc_dom = o.m_rgn_dealloc_dom;
      m_tag_env = o.m_tag_env;
    }
    return *this;
  }

  region_domain_t &operator=(region_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_alloc = std::move(o.m_alloc);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_base_dom = std::move(o.m_base_dom);      
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_rgn_init_dom = std::move(o.m_rgn_init_dom);
      m_alloc_site_dom = std::move(o.m_alloc_site_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_tag_env = std::move(o.m_tag_env);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    bool res = (!is_bottom() && m_base_dom.is_top() && m_rgn_counting_dom.is_top());
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

  bool operator<=(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    CRAB_LOG("region-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                        << "\n\t" << o << "\n";);
    if (is_bottom() || o.is_top()) {
      CRAB_LOG("region-leq", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) {
      CRAB_LOG("region-leq", crab::outs() << "Result2=0\n";);
      return false;
    }

    if (!(m_rgn_counting_dom <= o.m_rgn_counting_dom)) {
      CRAB_LOG("region-leq", crab::outs() << "Result3=0\n";);
      return false;
    }

    if (!(m_rgn_init_dom <= o.m_rgn_init_dom)) {
      CRAB_LOG("region-leq", crab::outs() << "Result4=0\n";);
      return false;
    }

    if (crab_domain_params_man::get().region_allocation_sites()) {
      if (!(m_alloc_site_dom <= o.m_alloc_site_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result5=0\n";);
        return false;
      }
    }
    
    if (crab_domain_params_man::get().region_deallocation()) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result6=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_tag_analysis()) {
      if (!(m_tag_env <= o.m_tag_env)) {
        CRAB_LOG("region-leq", crab::outs() << "Result7=0\n";);
        return false;
      }
    }

    // perform the common renaming
    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);
    base_abstract_domain_t left_dom(m_base_dom);
    base_abstract_domain_t right_dom(o.m_base_dom);
    base_variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings = m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = o.m_var_map.find(v);
      if (it != o.m_var_map.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
      }
    }

    left_dom.forget(out_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.forget(out_vars);
    right_dom.rename(right_vars, out_vars);

    CRAB_LOG("region-leq", crab::outs()
                               << "Inclusion test (after renaming):\n\t"
                               << left_dom << "\n\t" << right_dom << "\n";);

    bool res = left_dom <= right_dom;
    CRAB_LOG("region-leq", crab::outs() << "Result8=" << res << "\n";);
    return res;
  }

  void operator|=(const region_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) {
      *this = o;
      return;
    } else if (o.is_bottom()) {
      return;
    } else if (is_top() | o.is_top()) {
      set_to_top();
      return;
    }

    CRAB_LOG("region", crab::outs() << "Join " << *this << " and " << o << "=");

    variable_vector_t left_regions, right_regions;
    base_variable_vector_t regions_left_base_vars, regions_right_base_vars;

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
                     m_rev_var_map, m_alloc, m_base_dom, o.m_base_dom,
                     o.m_alloc);
      do_join(o);
    } else if (!refine_left && refine_right) {
      region_domain_t right(o);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_base_dom,
                     m_base_dom, m_alloc);
      do_join(right);
    } else if (refine_left && refine_right) {
      region_domain_t right(o);
      refine_regions(left_regions, regions_right_base_vars, m_var_map,
                     m_rev_var_map, m_alloc, m_base_dom, right.m_base_dom,
                     right.m_alloc);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_base_dom,
                     m_base_dom, m_alloc);
      do_join(right);
    } else {
      do_join(o);
    }
  }

  region_domain_t operator|(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() | o.is_top()) {
      region_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("region", crab::outs() << "Join " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 | v2; };

    variable_vector_t left_regions, right_regions;
    base_variable_vector_t regions_left_base_vars, regions_right_base_vars;

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
      region_domain_t left(*this);
      refine_regions(left_regions, regions_right_base_vars, left.m_var_map,
                     left.m_rev_var_map, left.m_alloc, left.m_base_dom,
                     o.m_base_dom, o.m_alloc);
      region_domain_t res(std::move(do_join_or_widening(left, o, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (!refine_left && refine_right) {
      // Refine right by propagating information from left's regions
      region_domain_t right(o);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_base_dom,
                     m_base_dom, m_alloc);
      region_domain_t res(
          std::move(do_join_or_widening(*this, right, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (refine_left && refine_right) {
      // Refine both left and right
      region_domain_t left(*this);
      region_domain_t right(o);
      refine_regions(left_regions, regions_right_base_vars, left.m_var_map,
                     left.m_rev_var_map, left.m_alloc, left.m_base_dom,
                     right.m_base_dom, right.m_alloc);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_base_dom,
                     left.m_base_dom, left.m_alloc);
      region_domain_t res(
          std::move(do_join_or_widening(left, right, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else {
      region_domain_t res(
          std::move(do_join_or_widening(*this, o, base_dom_op)));

      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;
    }
  }

  region_domain_t operator&(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("region", crab::outs() << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return (v1 & v2);
    };
    region_domain_t res(std::move(do_meet_or_narrowing(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t operator||(const region_domain_t &o) const override {
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

    CRAB_LOG("region", crab::outs()
                           << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 || v2;
    };
    region_domain_t res(std::move(do_join_or_widening(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t widening_thresholds(
      const region_domain_t &o,
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

    CRAB_LOG("region", crab::outs()
                           << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [&thresholds](const base_abstract_domain_t &v1,
                                     const base_abstract_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    region_domain_t res(std::move(do_join_or_widening(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t operator&&(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("region", crab::outs() << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return (v1 && v2); };
    region_domain_t res(std::move(do_meet_or_narrowing(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (!is_bottom()) {
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;	
        if (crab_domain_params_man::get().region_deallocation()) {
          m_rgn_dealloc_dom.forget(v);
        }
	
      }
      if (crab_domain_params_man::get().region_allocation_sites()) {
	if (v.get_type().is_reference() || v.get_type().is_region()) {
	  m_alloc_site_dom -= v;
	}
      }
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env -= v;
      } 
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        m_base_dom -= it->second;
        m_rev_var_map.erase(it->second);
        m_var_map.erase(it);
      }
    }
  }

  // Initialize region rgn.
  void region_init(const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_init");
    crab::ScopedCrabStats __st__(domain_name() + ".region_init");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    
    if (is_bottom()) {
      return;
    }

    small_range count_num = m_rgn_counting_dom[rgn];
    if (count_num <= small_range::oneOrMore()) {
      CRAB_ERROR("region_domain::init_region: ", rgn,
                 " cannot be initialized twice");
    }

    // Set to zero the number of references
    m_rgn_counting_dom.set(rgn, small_range::zero());

    // No stores in the region
    m_rgn_init_dom.set(rgn, boolean_value::get_false());

    if (crab_domain_params_man::get().region_deallocation()) {
      // No deallocated objects in the region
      m_rgn_dealloc_dom.set(rgn, boolean_value::get_false());
    }
    if (crab_domain_params_man::get().region_allocation_sites()) {
      // Region does not contain any allocation site
      m_alloc_site_dom.set(rgn, alloc_site_env_t::value_type::bottom());
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      // Region does not contain any tag
      m_tag_env.set(rgn, tag_env_t::value_type::bottom());
    }
    
    if (!rgn.get_type().is_unknown_region()) {
      // Assign a ghost variable to rgn for modeling its content.
      make_base_variable(m_alloc, rgn);
    }
    CRAB_LOG("region", crab::outs() << "After region_init(" << rgn
                                    << ")=" << *this << "\n";);
  }

  void region_copy(const variable_t &lhs_rgn,
                   const variable_t &rhs_rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_copy");
    crab::ScopedCrabStats __st__(domain_name() + ".region_copy");

    ERROR_IF_NOT_REGION(lhs_rgn, __LINE__);
    ERROR_IF_NOT_REGION(rhs_rgn, __LINE__);
    if (lhs_rgn.get_type() != rhs_rgn.get_type()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn, ":=", rhs_rgn,
                 " with different types");
    }
    
    if (is_bottom()) {
      return;
    }
    
    m_rgn_counting_dom.set(lhs_rgn, m_rgn_counting_dom[rhs_rgn]);
    m_rgn_init_dom.set(lhs_rgn, m_rgn_init_dom[rhs_rgn]);
    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom.set(lhs_rgn, m_alloc_site_dom[rhs_rgn]);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.add(rhs_rgn, lhs_rgn);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.set(lhs_rgn, m_tag_env[rhs_rgn]);
    }

    if (lhs_rgn.get_type().is_unknown_region()) {
      return;
    }
    
    const base_variable_t &base_lhs = rename_var(lhs_rgn);
    const base_variable_t &base_rhs = rename_var(rhs_rgn);
    auto num_refs = m_rgn_counting_dom[rhs_rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      // the rhs region is a singleton then we use assign 
      auto ty = lhs_rgn.get_type();
      if (ty.is_bool_region()) {
	m_base_dom.assign_bool_var(base_lhs, base_rhs, false);
      } else if (ty.is_integer_region() || ty.is_real_region() ||
		 ty.is_reference_region()) {
	m_base_dom.assign(base_lhs, base_rhs);
      } else if (ty.is_bool_array_region() ||
		 ty.is_int_array_region() ||
		 ty.is_real_array_region()) {
	m_base_dom.array_assign(base_lhs, base_rhs);
      } else {
	CRAB_ERROR(domain_name() + "::region_copy with unexpected region type");
      }
    } else {
      // the rhs region might not be a singleton so we use expand.
      m_base_dom.forget({base_lhs});
      m_base_dom.expand(base_rhs, base_lhs);
    }
  }


  void region_cast(const variable_t &src_rgn,
		   const variable_t &dst_rgn) override {
    // A region_cast is used to cast unknown regions to typed regions,
    // or viceversa.
                   
    crab::CrabStats::count(domain_name() + ".count.region_cast");
    crab::ScopedCrabStats __st__(domain_name() + ".region_cast");
    
    ERROR_IF_NOT_REGION(src_rgn, __LINE__);
    ERROR_IF_NOT_REGION(dst_rgn, __LINE__);

    if (is_bottom()) {
      return;
    }

    m_rgn_counting_dom.set(dst_rgn, m_rgn_counting_dom[src_rgn]);
    m_rgn_init_dom.set(dst_rgn, m_rgn_init_dom[src_rgn]);
    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom.set(dst_rgn, m_alloc_site_dom[src_rgn]);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.add(src_rgn, dst_rgn);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.set(dst_rgn, m_tag_env[src_rgn]);
    }

    // Since we currently don't model unknown regions we bail out
    // here.

    CRAB_LOG("region",
	     CRAB_WARN("TODO region_cast ", src_rgn, ":", src_rgn.get_type(), " to ",
		       dst_rgn, ":", dst_rgn.get_type(), " in base domain"););    
  }
  
  // Create a new reference ref to region rgn.
  void ref_make(const variable_t &ref, const variable_t &rgn,
		const variable_or_constant_t &size /*unused*/,
		const allocation_site &as) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);
    
    if (is_bottom()) {
      return;
    }

    // Update region counting
    auto num_refs = m_rgn_counting_dom[rgn];
    m_rgn_counting_dom.set(rgn, num_refs.increment());

    if (crab_domain_params_man::get().region_allocation_sites()) {
      // Associate allocation site as to ref
      m_alloc_site_dom.set(ref, as);
    }
    
    // Assign a base domain variable to ref
    rename_var(ref);

    CRAB_LOG("region", crab::outs() << "After ref_make(" << ref << "," << rgn << ","
	                            << size << "," << as << ")=" << *this << "\n";);
  }

  void ref_free(const variable_t &rgn, const variable_t &ref) override {
    crab::CrabStats::count(domain_name() + ".count.ref_free");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_free");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);
    
    if (is_bottom()) {
      return;
    }

    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom -= ref;
    }
    
    // We conservatively mark the region's equivalence class as
    // possibly deallocated
    if (crab_domain_params_man::get().region_deallocation()) {
      if (!m_rgn_dealloc_dom.contains(rgn)) {
        CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn,
                                     " in ", m_rgn_dealloc_dom));
      } else {
        m_rgn_dealloc_dom.set(rgn, boolean_value::top());
      }
    }
    
    CRAB_LOG("region", crab::outs() << "After ref_free(" << rgn << "," << ref
                                    << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within rgn. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);    
    ERROR_IF_NOT_REF(ref, __LINE__);
    if ((rgn.get_type().is_bool_region() && !res.get_type().is_bool()) ||
	(rgn.get_type().is_integer_region() &&
	 !res.get_type().is_integer()) ||
	(rgn.get_type().is_real_region() && !res.get_type().is_real()) ||
	(rgn.get_type().is_reference_region() &&
           !res.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::ref_load: type of lhs ", res, " (",
		 res.get_type(), ") is not compatible with region ", rgn,
		 " (", rgn.get_type(), ")");
    }
    
    if (is_bottom()) {
      return;
    }

    const base_variable_t &base_res = rename_var(res);

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-load",
               CRAB_WARN("region_domain::ref_load: reference ", ref, " is null."););
      // set_to_bottom();
      m_base_dom -= base_res;
      return;
    }
    
    if (crab_domain_params_man::get().region_allocation_sites()) {
      if (res.get_type().is_reference()) {
	m_alloc_site_dom.set(res, m_alloc_site_dom[rgn]);
      }
    }

    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.set(res, m_tag_env[rgn]);
    }
    
    if (rgn.get_type().is_unknown_region()) {
      CRAB_LOG(
          "region-load",
          CRAB_WARN("region_domain::ref_load: skip unknown region ", rgn););
      m_base_dom -= base_res;
      return;
    }
        
    auto ref_load = [&rgn, &base_res, this](const base_variable_t &rgn_var) {
      // At this point references and region of references have been
      // translated to integers and region of integers, respectively.
      if (rgn_var.get_type() != base_res.get_type()) {
        CRAB_ERROR("region_domain::ref_load: ", "Type of region ", rgn,
                   " does not match with ", base_res);
      }
      if (base_res.get_type().is_bool()) {
        this->m_base_dom.assign_bool_var(base_res, rgn_var, false);
      } else if (base_res.get_type().is_integer() ||
                 base_res.get_type().is_real()) {
        this->m_base_dom.assign(base_res, rgn_var);
      } else {
        // variables of type base_variable_t cannot be REF_TYPE or ARR_TYPE
        CRAB_ERROR("region_domain::ref_load: unsupported type in ", base_res);
      }
    };
    
    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      // strong read
      CRAB_LOG("region-load", crab::outs() << "Reading from singleton\n";);
      if (auto region_var_opt = get_var(rgn)) {
        ref_load(*region_var_opt);
      } else {
        m_base_dom -= base_res;
      }
    } else {
      // weak read
      CRAB_LOG("region-load", crab::outs() << "Reading from non-singleton\n";);
      if (auto region_var_opt = get_var(rgn)) {
        base_variable_t fresh_region_var = make_base_variable(m_alloc, rgn);
        m_base_dom.expand(*region_var_opt, fresh_region_var);
        ref_load(fresh_region_var);
      } else {
        m_base_dom -= base_res;
      }
    }

    CRAB_LOG("region-load", crab::outs() << "After " << res << ":="
                                         << "ref_load(" << ref << "," << rgn
                                         << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in rgn.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);    
    ERROR_IF_NOT_REF(ref, __LINE__);
    if ((rgn.get_type().is_bool_region() && !val.get_type().is_bool()) ||
	(rgn.get_type().is_integer_region() &&
	 !val.get_type().is_integer()) ||
	(rgn.get_type().is_real_region() && !val.get_type().is_real()) ||
	(rgn.get_type().is_reference_region() &&
           !val.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::ref_store: type of value ", val, " (",
		 val.get_type(), ") is not compatible with region ", rgn,
		 " (", rgn.get_type(), ")");
    }
    
    if (is_bottom()) {
      return;
    }
    
    auto forget_region_if_tracked = [this](const variable_t &v) {
      if (v.get_type().is_region() && !v.get_type().is_unknown_region()) {
	base_variable_t rgn_var = rename_var(v);
	m_base_dom -= rgn_var;
      }
    };
    
    bool is_uninitialized_rgn = m_rgn_init_dom[rgn].is_false();
    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-store",
               CRAB_WARN("region_domain::ref_store: reference ", ref,
                         " is null. "););
      // set_to_bottom();
      forget_region_if_tracked(rgn);      
      return;
    }
    
    auto ref_store = [&val, this](base_abstract_domain_t &base_dom,
				  const base_variable_t &base_rgn_var) {
      // At this point region of references has been translated to region of
      // integers. val can be any scalar including a reference.
      if (val.get_type().is_bool()) {
        if (val.is_constant()) {
          base_dom.assign_bool_cst(
              base_rgn_var,
              (val.is_bool_true() ? base_linear_constraint_t::get_true()
                                  : base_linear_constraint_t::get_false()));
        } else {
          base_dom.assign_bool_var(base_rgn_var, rename_var(val.get_variable()),
                                   false);
        }
      } else if (val.get_type().is_integer() || val.get_type().is_real() ||
                 val.get_type().is_reference()) {
        if (val.is_constant()) {
          base_dom.assign(base_rgn_var, val.get_constant());
        } else {
          base_dom.assign(base_rgn_var, rename_var(val.get_variable()));
        }
      } else {
        CRAB_ERROR(domain_name(), "::ref_store: unsupported type ", val.get_type());
      }
    };

    auto num_refs = m_rgn_counting_dom[rgn];
    // A region can have more than one reference but we can still
    // perform a strong update as long as nobody wrote yet in the
    // region. The use of "is_uninitialized_rgn" avoids fooling the
    // abstract domain with code like this:
    //
    //   r1 := make_ref(R);
    //   r2 := make_ref(R);
    //   store_ref(R, r1, v1);
    //   store_ref(R, r2, v2);
    //
    if (is_uninitialized_rgn || num_refs.is_zero() || num_refs.is_one()) {
      /* strong update */
      CRAB_LOG("region-store", crab::outs() << "Performing strong update\n";);

      if (!rgn.get_type().is_unknown_region()) {
	base_variable_t rgn_var = rename_var(rgn);
	ref_store(m_base_dom, rgn_var);
      } else {
	CRAB_LOG(
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      } 
      
      if (crab_domain_params_man::get().region_allocation_sites()) {
	if (val.get_type().is_reference()) {
	  if (val.is_reference_null()) {
	    // reset to empty set because of strong update
	    m_alloc_site_dom.set(rgn, allocation_sites::bottom());
	  } else {
	    assert(val.is_variable());
	    m_alloc_site_dom.set(rgn, m_alloc_site_dom[val.get_variable()]);
	  }
	}
      }

      if (crab_domain_params_man::get().region_tag_analysis()) {
	if (val.is_variable()) {
	  m_tag_env.set(rgn, m_tag_env[val.get_variable()]);
	}
      }
	
    } else {
      /* weak update */
      CRAB_LOG("region-store", crab::outs() << "Performing weak update\n";);

      if (!rgn.get_type().is_unknown_region()) {
	base_variable_t rgn_var = rename_var(rgn);	
	base_abstract_domain_t tmp(m_base_dom);
	ref_store(tmp, rgn_var);
	m_base_dom |= tmp;
      } else {
	CRAB_LOG( 
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      }
      
      if (crab_domain_params_man::get().region_allocation_sites()) {
	if (val.get_type().is_reference()) {
	  if (val.is_variable()) {
	    m_alloc_site_dom.set(rgn, m_alloc_site_dom[rgn] |
				 m_alloc_site_dom[val.get_variable()]);
	  }
	}
      }
      if (crab_domain_params_man::get().region_tag_analysis()) {
	if (val.is_variable()) {
	  m_tag_env.set(rgn, m_tag_env[rgn] | m_tag_env[val.get_variable()]);
	}
      }
      
    }
    CRAB_LOG("region-store", crab::outs()
                                 << "After ref_store(" << ref << "," << rgn
                                 << "," << val << ")=" << *this << "\n";);
  }

  // Create a new reference ref2 to region rgn2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &rgn1,
               const variable_t &ref2, const variable_t &rgn2,
               const linear_expression_t &offset) override {
    crab::CrabStats::count(domain_name() + ".count.ref_gep");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_gep");

    ERROR_IF_NOT_REGION(rgn1, __LINE__);
    ERROR_IF_NOT_REGION(rgn2, __LINE__);
    ERROR_IF_NOT_REF(ref1, __LINE__);
    ERROR_IF_NOT_REF(ref2, __LINE__);
    
    auto eval = [this](const linear_expression_t &e) {
      interval_t r = e.constant();
      for (auto p : e) {
        r += p.first * m_base_dom.operator[](rename_var(p.second));
      }
      return r;
    };

    if (!is_bottom()) {
      auto boffset = rename_linear_expr(offset);
      m_base_dom.assign(rename_var(ref2), rename_var(ref1) + boffset);
      if (!(rgn1 == rgn2 && (eval(offset) == (number_t(0))))) {
        // Update region counting
        auto num_refs = m_rgn_counting_dom[rgn2];
        m_rgn_counting_dom.set(rgn2, num_refs.increment());
      }

      if (crab_domain_params_man::get().region_allocation_sites()) {
	m_alloc_site_dom.set(ref2, m_alloc_site_dom[ref1]);
      }

      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(ref2, m_tag_env[ref1]);
      }
      
      if (crab_domain_params_man::get().region_deallocation()) {
        bool found1 = m_rgn_dealloc_dom.contains(rgn1);
        bool found2 = m_rgn_dealloc_dom.contains(rgn2);
        if (!found1) {
          CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn1,
                                       " in ", m_rgn_dealloc_dom));
        }
        if (!found2) {
          CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn2,
                                       " in ", m_rgn_dealloc_dom));
        }
	
        if (found1 && found2) {
          m_rgn_dealloc_dom.join(rgn1, rgn2);
        } else if (!found1 && found2) {
	  //m_rgn_dealloc_dom.remove_equiv_class(rgn2);
	} else if (found1 && !found2) {
	  //m_rgn_dealloc_dom.remove_equiv_class(rgn1);	  
	} else {
	  // do nothing since rgn1 and rgn2 are already untracked.
	}
      }
    }

    CRAB_LOG("region", crab::outs()
                           << "After (" << rgn2 << "," << ref2
                           << ") := ref_gep(" << rgn1 << "," << ref1 << " + "
                           << offset << ")=" << *this << "\n";);
  }

  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(const variable_t &lhs,
                           const variable_t &ref /*unused*/,
                           const variable_t &rgn,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load_from_array");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load_from_array");

    if (is_bottom()) {
      return;
    }

    const base_variable_t &base_lhs = rename_var(lhs);

    /***
     * These syntactic checks are only needed if the abstract domain
     * API is called directly.
     **/
    if (!lhs.get_type().is_integer() && !lhs.get_type().is_bool() &&
        !lhs.get_type().is_real()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_load_from_array: ", lhs,
                                   " must be bool/int/real type"););
      m_base_dom -= base_lhs;
    }

    if (!rgn.get_type().is_array_region()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_load_from_array: ", rgn,
                                   " must be an array region"););
      m_base_dom -= base_lhs;
      return;
    }

    // TODO: check that the array element matches the type of lhs
    // (bool or int)

    if (boost::optional<base_variable_t> arr_var_opt = get_var(rgn)) {
      auto num_refs = m_rgn_counting_dom[rgn];
      auto b_elem_size = rename_linear_expr(elem_size);
      if (num_refs.is_zero() || num_refs.is_one()) {
        CRAB_LOG("region", crab::outs() << "Reading from singleton\n";);
        auto b_index = rename_linear_expr(index);
        m_base_dom.array_load(base_lhs, *arr_var_opt, b_elem_size, b_index);
      } else {
        // TODO: get the bitwidth from index
        base_variable_t unknown_index(m_alloc.next(), crab::INT_TYPE, 32);
        m_base_dom.array_load(base_lhs, *arr_var_opt, b_elem_size,
                              unknown_index);
      }
    } else {
      m_base_dom -= base_lhs;
    }

    // TODO: tag analysis
  }

  // Treat region as an array and perform an array store.
  void ref_store_to_array(const variable_t &ref /*unused*/,
                          const variable_t &rgn,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store_to_array");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store_to_array");

    if (is_bottom()) {
      return;
    }

    /***
     * These syntactic checks are only needed if the abstract domain
     * API is called directly.
     **/
    if (!rgn.get_type().is_array_region()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_store_to_array: ", rgn,
                                   " must be an array region"););
      return;
    }

    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    // TODO: check that the array element matches the type of val
    // (bool or integer)

    auto num_refs = m_rgn_counting_dom[rgn];
    base_variable_t arr_var = rename_var(rgn);
    auto b_elem_size = rename_linear_expr(elem_size);
    auto b_val = rename_linear_expr(val);
    if (num_refs.is_zero() || num_refs.is_one()) {
      CRAB_LOG("region", crab::outs() << "Reading from singleton\n";);
      auto b_index = rename_linear_expr(index);
      m_base_dom.array_store(arr_var, b_elem_size, b_index, b_val, false);

    } else {
      // TODO: get the bitwidth from index
      base_variable_t unknown_index(m_alloc.next(), crab::INT_TYPE, 32);
      m_base_dom.array_store(arr_var, b_elem_size, unknown_index, b_val, false);
    }

    // TODO: tag analysis    
  }

  void ref_assume(const reference_constraint_t &ref_cst) override {
    crab::CrabStats::count(domain_name() + ".count.ref_assume");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_assume");

    if (!is_bottom()) {
      if (ref_cst.is_tautology()) {
        return;
      }
      if (ref_cst.is_contradiction()) {
        set_to_bottom();
        return;
      }

      if (crab_domain_params_man::get().region_allocation_sites()) {
	if (ref_cst.is_equality() && ref_cst.is_binary()) {
	  auto lhs_as = m_alloc_site_dom[ref_cst.lhs()];
	  auto rhs_as = m_alloc_site_dom[ref_cst.rhs()];
	  // **Important note about soundness**: if the program
	  // environment is not modeled precisely enough we can
	  // conclude *incorrectly* that p == q is definitely not
	  // true. At least here, we skip the cases where either p or
	  // q does not have at least one allocation site.
	  if (!lhs_as.is_bottom() && !rhs_as.is_bottom()) {	  
	    allocation_sites inter = lhs_as & rhs_as;	    
	    if (inter.is_bottom()) {
	      // if they do not have any common allocation site then
	      // they cannot be the same address.	      
	      set_to_bottom();
	      return;
	    }
	    // If ref_cst is a disequality we cannot make the state
	    // bottom by looking only at allocation sites. Even if
	    // both references point to the same allocation site they
	    // can still be different.
	  }
	}
      }
      
      auto b_lin_csts = convert_ref_cst_to_linear_cst(ref_cst);
      m_base_dom += b_lin_csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
    CRAB_LOG("region",
             crab::outs() << "ref_assume(" << ref_cst << ")" << *this << "\n";);
  }

  void ref_to_int(const variable_t &rgn, const variable_t &ref_var,
                  const variable_t &int_var) override {
    crab::CrabStats::count(domain_name() + ".count.ref_to_int");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_to_int");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);
    
    if (!is_bottom()) {
      base_variable_t src_var = rename_var(ref_var);
      base_variable_t dst_var = rename_var(int_var);
      m_base_dom.assign(dst_var, src_var);

      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(int_var, m_tag_env[ref_var]);
      }
    }
  }

  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    crab::CrabStats::count(domain_name() + ".count.int_to_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".int_to_ref");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      base_variable_t src_var = rename_var(int_var);
      base_variable_t dst_var = rename_var(ref_var);

      m_base_dom.assign(dst_var, src_var);

      if (crab_domain_params_man::get().region_allocation_sites()) {
	m_alloc_site_dom -= ref_var;
      }

      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(ref_var, m_tag_env[int_var]);
      }
      
      // Update region counting
      auto num_refs = m_rgn_counting_dom[rgn];
      m_rgn_counting_dom.set(rgn, num_refs.increment());
    }
  }

  // This default implementation is expensive because it will call the
  // join.
  DEFAULT_SELECT_REF(region_domain_t)
  
  // arithmetic operations
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }
    }
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), k);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, m_tag_env[y]);
      }      
    }
  }
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_base_dom.assign(rename_var(x), b_e);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	merge_tags(x, e.variables());
      }
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    if (!is_bottom()) {
      auto b_e1 = rename_linear_expr(e1);
      auto b_e2 = rename_linear_expr(e2);
      auto b_cond = rename_linear_cst(cond);
      m_base_dom.select(rename_var(lhs), b_cond, b_e1, b_e2);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	tag_set tags = tag_set::bottom();
	for (auto const &v: e1.variables()) {
	  tags = tags | m_tag_env[v];
	}
	for (auto const &v: e2.variables()) {
	  tags = tags | m_tag_env[v];
	}
	m_tag_env.set(lhs, tags);
      }      
    }
  }
  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_base_dom.backward_assign(rename_var(x), b_e, invariant.m_base_dom);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }      
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, rename_var(x), rename_var(y), z,
                                invariant.m_base_dom);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, rename_var(x), rename_var(y), rename_var(z),
                                invariant.m_base_dom);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      auto b_csts = rename_linear_cst_sys(csts);
      m_base_dom += b_csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  // not part of the numerical_domains api but it should be
  void set(const variable_t &x, interval_t intv) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      if (intv.is_bottom()) {
        set_to_bottom();
      } else {
        m_base_dom.set(rename_var(x), intv);
      }
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, tag_env_t::value_type::bottom()); 
      }      
    }
  }

  interval_t operator[](const variable_t &x) override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom[rename_var(x)];
    }
  }

  // int cast operations and bitwise operations
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(dst), rename_var(src));
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(dst, m_tag_env[src]);
      }            
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }            
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), z);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, m_tag_env[y]);
      }            
    }
  }

  // boolean operations
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = rename_linear_cst(rhs);
      m_base_dom.assign_bool_cst(rename_var(lhs), b_rhs);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	merge_tags(lhs, rhs.variables());
      }            
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {

      if (crab_domain_params_man::get().region_allocation_sites()) {
	if ((rhs.is_equality() || rhs.is_disequality()) && rhs.is_binary()) {
	  auto op1_as = m_alloc_site_dom[rhs.lhs()];
	  auto op2_as = m_alloc_site_dom[rhs.rhs()];
	  // -- See note about soundness in ref_asume.
	  if (!op1_as.is_bottom() && !op2_as.is_bottom()) {	  
	    allocation_sites inter = op1_as & op2_as;
	    if (inter.is_bottom()) {
	      // if they do not have any common allocation site then
	      // they cannot be the same address.
	      if (rhs.is_equality()) {
		auto false_cst =  base_linear_constraint_t::get_false();
		m_base_dom.assign_bool_cst(rename_var(lhs), false_cst);
	      } else {
		assert(rhs.is_disequality());
		auto true_cst =  base_linear_constraint_t::get_true();
		m_base_dom.assign_bool_cst(rename_var(lhs), true_cst);
	      }
	      return;
	    }
	  }
	}
      }
      
      auto b_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.assign_bool_cst(rename_var(lhs), b_rhs);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	merge_tags(lhs, rhs.variables());
      }            
    }
  }

  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.assign_bool_var(rename_var(lhs), rename_var(rhs), is_not_rhs);
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(lhs, m_tag_env[rhs]);
      }            
    }
  }
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.apply_binary_bool(op, rename_var(x), rename_var(y),
                                   rename_var(z));
      
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }            
    }
  }
  void assume_bool(const variable_t &v, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_bottom()) {
      m_base_dom.assume_bool(rename_var(v), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }
  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      m_base_dom.select_bool(rename_var(lhs), rename_var(cond),
			     rename_var(b1), rename_var(b2));

      if (crab_domain_params_man::get().region_tag_analysis()) {
	// it can be improve if we know if cond is true or false
	m_tag_env.set(lhs, m_tag_env[b1] | m_tag_env[b2]);
      }                  
    }
  }  
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = rename_linear_cst(rhs);
      m_base_dom.backward_assign_bool_cst(rename_var(lhs), b_rhs,
                                          invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }      
    }
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.backward_assign_bool_cst(rename_var(lhs), b_rhs,
                                          invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_var(rename_var(lhs), rename_var(rhs),
                                          is_not_rhs, invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.backward_apply_binary_bool(op, rename_var(x), rename_var(y),
                                            rename_var(z),
                                            invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }

  // array operations
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".array_init");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_lb_idx = rename_linear_expr(lb_idx);
      auto b_ub_idx = rename_linear_expr(ub_idx);
      auto b_val = rename_linear_expr(val);
      m_base_dom.array_init(rename_var(a), b_elem_size, b_lb_idx, b_ub_idx,
                            b_val);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }      
      
    }
  }
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".array_load");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      m_base_dom.array_load(rename_var(lhs), rename_var(a), b_elem_size, b_i);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      auto b_v = rename_linear_expr(v);
      m_base_dom.array_store(rename_var(a), b_elem_size, b_i, b_v,
                             is_strong_update);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store_range");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      auto b_j = rename_linear_expr(j);
      auto b_v = rename_linear_expr(v);
      m_base_dom.array_store_range(rename_var(a), b_elem_size, b_i, b_j, b_v);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".array_assign");

    if (!is_bottom()) {
      m_base_dom.array_assign(rename_var(lhs), rename_var(rhs));
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_init");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_lb_idx = rename_linear_expr(lb_idx);
      auto b_ub_idx = rename_linear_expr(ub_idx);
      auto b_val = rename_linear_expr(val);
      m_base_dom.backward_array_init(rename_var(a), b_elem_size, b_lb_idx,
                                     b_ub_idx, b_val, invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_load");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      m_base_dom.backward_array_load(rename_var(lhs), rename_var(a),
                                     b_elem_size, b_i, invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      auto b_v = rename_linear_expr(v);
      m_base_dom.backward_array_store(rename_var(a), b_elem_size, b_i, b_v,
                                      is_strong_update, invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store_range");

    if (!is_bottom()) {
      auto b_elem_size = rename_linear_expr(elem_size);
      auto b_i = rename_linear_expr(i);
      auto b_j = rename_linear_expr(j);
      auto b_v = rename_linear_expr(v);
      m_base_dom.backward_array_store_range(rename_var(a), b_elem_size, b_i,
                                            b_j, b_v, invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_assign");

    if (!is_bottom()) {
      m_base_dom.backward_array_assign(rename_var(lhs), rename_var(rhs),
                                       invariant.m_base_dom);
      if (crab_domain_params_man::get().region_tag_analysis()) {
	// TODO tag analysis
      }            
    }
  }

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
        if (crab_domain_params_man::get().region_deallocation()) {
          m_rgn_dealloc_dom.forget(v);
        }	
      }
      if (crab_domain_params_man::get().region_tag_analysis()) {
	m_tag_env -= v;
      }
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        base_vars.push_back(it->second);
        m_rev_var_map.erase(it->second);
        m_var_map.erase(it);
      }
    }
    m_base_dom.forget(base_vars);
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
      if (v.get_type().is_unknown_region()) {
        // skip those for now
        continue;
      }
      base_vars.push_back(rename_var(v));
    }
    m_base_dom.project(base_vars);

    // -- update m_var_map and m_rev_var_map
    std::vector<variable_t> var_map_to_remove;
    var_map_to_remove.reserve(m_var_map.size());
    for (auto &kv : m_var_map) {
      if (!std::binary_search(sorted_variables.begin(), sorted_variables.end(),
                              kv.first)) {
        var_map_to_remove.push_back(kv.first);
        m_rev_var_map.erase(kv.second);
      }
    }
    for (auto &v : var_map_to_remove) {
      m_var_map.erase(v);
    }

    m_rgn_counting_dom.project(sorted_variables);
    m_rgn_init_dom.project(sorted_variables);
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

  void normalize() override {}
  void minimize() override {}

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_ERROR("region_domain::expand not implemented");
  }

  // The code might be simpler and more efficient if m_var_map is a
  // std::map.  However, this would slow down lookups.  The
  // expectation is that this operation should not been called often
  // since the region domain is usually at the root of the hierarchy
  // of domains.
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

    //
    // rename m_var_map and m_rev_var_map
    // 

    std::vector<std::pair<variable_t, variable_t>> zipped_vars;
    for (unsigned i=0,sz=from.size();i<sz;++i) {
      zipped_vars.push_back({from[i],to[i]});
    }

    auto compare_zipped_pair =
      [](const std::pair<variable_t, variable_t> &p1,
	 const std::pair<variable_t, variable_t> &p2) {
	return p1.first < p2.first;
      };
    

    std::sort(zipped_vars.begin(), zipped_vars.end(), compare_zipped_pair);

    // Needed because we cannot insert/delete while iterating over
    // m_var_map
    std::vector<variable_t> var_map_to_remove;
    std::vector<std::pair<variable_t, base_variable_t>> var_map_to_insert;
    
    for (auto &kv: m_var_map) {
      const variable_t &old_var = kv.first;
      const base_variable_t &old_base_var = kv.second;

      std::pair<variable_t,variable_t> p(old_var, old_var);
      auto lower = std::lower_bound(zipped_vars.begin(), zipped_vars.end(), p,
				    compare_zipped_pair);
      
      if (!(lower == zipped_vars.end() || old_var < (*lower).first)) {
	// found
	variable_t new_var = (*lower).second;
	base_variable_t new_base_var = rename_var(new_var);
	
	var_map_to_remove.push_back(old_var);
	var_map_to_insert.push_back({new_var, new_base_var});
	
	m_rev_var_map.erase(old_base_var);
	m_rev_var_map.insert({new_base_var, new_var});
      }
    }

    for (auto &v: var_map_to_remove) {
      m_var_map.erase(v);
    }
    for (auto &p: var_map_to_insert) {
      m_var_map.insert(p);
    }


    // Rename the rest 
    m_rgn_counting_dom.rename(from, to);
    m_rgn_init_dom.rename(from, to);
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
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    //=================================================================//
    //       Special intrinsics supported by the region domain
    //=================================================================//
    // ---Use-after-free analysis----
    // 
    // b := is_unfreed_or_null(rgn, ref)
    //       b is true if the reference does not point to a freed
    //       memory object or it is null.
    // unfreed_or_null(rgn, ref)
    //       ensures that the reference does not point to a freed memory object
    // 
    // ---Nullity analysis---
    // 
    // nonnull(ref)
    //       ensures that the reference is not null
    //
    // ---Array bounds analysis----
    //
    // b := is_dereferenceable(rgn, ref, sz)
    //       b is true if the memory pointed by [ref, ref+sz] is
    //       dereferenceable
    // 
    // ---Tag analysis---
    /// 
    // add_tag(rgn, ref, TAG)
    //       Add TAG to the set of tags already associated with the
    //       data pointed by ref within rgn.
    // b := does_not_have_tag(rgn, ref, TAG)
    //       b is true if the data pointed by ref within rgn does
    //       *definitely* not have tag TAG.
    //=================================================================//
    auto error_if_not_arity = [&name, &inputs, &outputs](unsigned num_inputs,
                                                         unsigned num_outputs) {
      if (inputs.size() != num_inputs || outputs.size() != num_outputs) {
        CRAB_ERROR("Intrinsics ", name, " unexpected number of parameters");
      }
    };
    auto error_if_not_variable = [&name](const variable_or_constant_t &vc) {
      if (!vc.is_variable()) {
	CRAB_ERROR("Intrinsics ", name, " expected a variable input");
      }
    };      
    auto error_if_not_constant = [&name](const variable_or_constant_t &vc) {
      if (!vc.is_constant()) {
	CRAB_ERROR("Intrinsics ", name, " expected a constant input");
      }
    };      
    
    auto error_if_not_bool = [&name](const variable_t &var) {
      if (!var.get_type().is_bool()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var, " should be Bool");
      }
    };
    auto error_if_not_rgn = [&name](const variable_t &var) {
      if (!var.get_type().is_region()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var,
                   " should be a region");
      }
    };
    auto error_if_not_ref = [&name](const variable_t &var) {
      if (!var.get_type().is_reference()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var,
                   " should be a reference");
      }
    };

    auto set_bool_var_to_true = [this](const variable_t &bool_var) {
      /// Require that the base domain can reason about booleans
      operator-=(bool_var);
      assume_bool(bool_var, false /*not negated*/);
    };
    
    if (is_bottom()) {
      return;
    }
    
    if (name == "is_unfreed_or_null") {
      if (crab_domain_params_man::get().region_deallocation()) {
	error_if_not_arity(2, 1);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);      
	variable_t rgn = inputs[0].get_variable();
	variable_t ref = inputs[1].get_variable();
	variable_t bv = outputs[0];
	error_if_not_rgn(rgn);
	error_if_not_ref(ref);
	error_if_not_bool(bv);
	// the reference is definitely null
	if (is_null_ref(ref).is_true()) {
	  set_bool_var_to_true(bv);
	} else {
	  bool is_allocated = (m_rgn_dealloc_dom.contains(rgn) &&
			       m_rgn_dealloc_dom[rgn].is_false());
	  /// the reference belongs to a memory region that doesn't have
	  /// any deallocated memory object.
	  if (is_allocated) {
	    set_bool_var_to_true(bv);
	  }
	}
      }
    } else if (name == "unfreed_or_null") {
      if (crab_domain_params_man::get().region_deallocation()) {
	error_if_not_arity(2, 0);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	variable_t rgn = inputs[0].get_variable();
	variable_t ref = inputs[1].get_variable();
	error_if_not_rgn(rgn);
	error_if_not_ref(ref);
	if (!is_null_ref(ref).is_true()) {
	  // We can only mark the region as "deallocated" if it doesn't
	  // have any reference yet. Even if it has only one we cannot
	  // tell whether it's the same ref than ref.
	  auto num_refs = m_rgn_counting_dom[rgn];
	  if (num_refs.is_zero()) {
	    m_rgn_dealloc_dom.set(rgn, boolean_value::get_false());
	  }
	}
      }
    } else if (name == "nonnull") {
      error_if_not_arity(1, 0);
      error_if_not_variable(inputs[0]);      
      variable_t ref = inputs[0].get_variable();
      error_if_not_ref(ref);
      auto nonnull_cst = reference_constraint_t::mk_gt_null(ref);
      ref_assume(nonnull_cst);
    } else if (name == "add_tag") {
      if (crab_domain_params_man::get().region_tag_analysis()) {
	error_if_not_arity(3, 0);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	error_if_not_constant(inputs[2]);
	variable_t rgn = inputs[0].get_variable();
	error_if_not_rgn(rgn);
	variable_t ref = inputs[1].get_variable();
	error_if_not_ref(ref);
	tag_t tag(inputs[2].get_constant());

	// We ignore ref and merge the tag with the tags already in
	// the region
	m_tag_env.set(rgn, m_tag_env[rgn] | tag_set(tag));
      }
    } else if (name == "does_not_have_tag") {
      if (crab_domain_params_man::get().region_tag_analysis()) {
	error_if_not_arity(3, 1);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	error_if_not_constant(inputs[2]);
	error_if_not_bool(outputs[0]);
	variable_t rgn = inputs[0].get_variable();
	error_if_not_rgn(rgn);	
	variable_t ref = inputs[1].get_variable();
	error_if_not_ref(ref);	
	tag_t tag(inputs[2].get_constant());
	variable_t bv = outputs[0];
	tag_set tags = m_tag_env[rgn];
	if (!(tag_set(tag) <= tags)) {
	  set_bool_var_to_true(bv);
	} else {
	  operator-=(bv);			  
	} 
      }
    } else if (name == "is_dereferenceable") {
      error_if_not_arity(3, 1);
      CRAB_WARN("This region domain ignores the is_dereferenceable intrinsics");
      variable_t v = outputs[0];
      operator-=(v);
    } else {
      // pass the intrinsics to the base domain
      //=== base domain ===/
      std::vector<base_variable_or_constant_t> base_inputs;
      std::vector<base_variable_t> base_outputs;
      base_inputs.reserve(inputs.size());
      base_outputs.reserve(outputs.size());
      for (unsigned i = 0, sz = inputs.size(); i < sz; ++i) {
	if (inputs[i].get_type().is_unknown_region()) {
	  // rename_var does not support unknown regions so we bail
	  // out. The base domain shouldn't care about regions anyway.
	  return;
	}
        base_inputs.push_back(inputs[i].is_variable() ?
	     base_variable_or_constant_t(rename_var(inputs[i].get_variable())) :
	     base_variable_or_constant_t(inputs[i].get_constant(),
					 inputs[i].get_type()));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
	if (outputs[i].get_type().is_unknown_region()) {
	  // rename_var does not support unknown regions so we bail
	  // out. The base domain shouldn't care about regions anyway.
	  return;
	}	
        base_outputs.push_back(rename_var(outputs[i]));
      }
      m_base_dom.intrinsic(name, base_inputs, base_outputs);
    }
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const region_domain_t &invariant) override {
    if (!is_bottom()) {
      //=== base domain ===/
      std::vector<base_variable_or_constant_t> base_inputs;
      std::vector<base_variable_t> base_outputs;      
      base_inputs.reserve(inputs.size());
      base_outputs.reserve(outputs.size());
      for (unsigned i = 0, sz = inputs.size(); i < sz; ++i) {
        base_inputs.push_back(inputs[i].is_variable() ?
	     base_variable_or_constant_t(rename_var(inputs[i].get_variable())):
	     base_variable_or_constant_t(inputs[i].get_constant(),
					 inputs[i].get_type()));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
        base_outputs.push_back(rename_var(outputs[i]));
      }
      m_base_dom.backward_intrinsic(name, base_inputs, base_outputs,
                                    invariant.m_base_dom);
    }
  }
  /* end intrinsics operations */

  /* Begin abstract_domain_results_api */

  // Return a 3-valued boolean. If true then ref is definitely
  // null. If false then ref is definitely non-null. Otherwise, we do
  // not know.
  boolean_value is_null_ref(const variable_t &ref) override {
    if (is_bottom()) {
      return boolean_value::bottom();
    }
    
    if (!ref.get_type().is_reference()) {
      return boolean_value::get_false();
    }
    
    if (auto var_opt = get_var(ref)) {
      interval_t ival = m_base_dom[*var_opt];
      number_t zero(0);
      
      if (!(interval_t(zero) <= ival)) {
	return boolean_value::get_false();
      }
      
      boost::optional<number_t> x = ival.lb().number();
      boost::optional<number_t> y = ival.ub().number();
      if (x && y && *x == zero && *y == zero) {
	return boolean_value::get_true();
      }
    }
    
    return boolean_value::top();
  }

  virtual bool get_allocation_sites(const variable_t &ref,
				    std::vector<allocation_site> &alloc_sites) override {

    if (!crab_domain_params_man::get().region_allocation_sites() || !ref.get_type().is_reference()) {
      return false;
    }
    
    allocation_sites out = m_alloc_site_dom[ref];
    
    if (out.is_top()) {
      return false;
    }
    
    if (!out.is_bottom()) {
      for (auto it = out.begin(), et = out.end(); it!=et; ++it) {
	alloc_sites.push_back(*it);
      }
    }
    return true;
  }

  virtual bool get_tags(const variable_t &rgn,
			const variable_t &ref /*unused*/,
			std::vector<uint64_t> &tags) override {

    
    if (!crab_domain_params_man::get().region_tag_analysis() ||
	(!ref.get_type().is_reference() || !rgn.get_type().is_region())) {
      return false;
    }

    tag_set tag_set = m_tag_env[rgn];    
    
    if (tag_set.is_top()) {
      return false;
    }
    
    if (!tag_set.is_bottom()) {
      for (auto it = tag_set.begin(), et = tag_set.end(); it!=et; ++it) {
	tags.push_back((*it).index());
      }
    }
    return true;
  }
  
  /* End abstract_domain_results_api */
  
  // Only integer or boolean variables are exposed.
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
      for (base_linear_constraint_t cst :
           m_base_dom.to_linear_constraint_system()) {
        if (boost::optional<linear_constraint_t> out_cst =
                rev_rename_linear_cst(cst, ignore_references)) {
          out_csts += *(out_cst);
        }
      }
      return out_csts;
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR("region_domain::to_disjunctive_linear_constraint_system not "
               "implemented");
  }

  std::string domain_name() const override {
    return "RegionDomain(" + m_base_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("region-print",
	       o << "(RgnCounter=" << m_rgn_counting_dom << ","
	       << "RgnInit=" << m_rgn_init_dom;
	       if (crab_domain_params_man::get().region_allocation_sites()) {
		 o << "," << "AllocSites=" << m_alloc_site_dom;
	       }
	       if (crab_domain_params_man::get().region_deallocation()) {
		 o << "," << "RgnDealloc=" << m_rgn_dealloc_dom;
	       }
	       if (crab_domain_params_man::get().region_tag_analysis()) {
		 o << "," << "Tags=" << m_tag_env;		 
	       }
	       o << ")\n";);
      
      CRAB_LOG(
          "region-print2", o << "MapVar={"; for (auto &kv
                                                 : m_var_map) {
            o << kv.first << "->" << kv.second << ";";
          } o << "},"
              << "BaseDom=" << m_base_dom << ")\n";);
      std::unordered_map<std::string, std::string> renaming_map;
      for (auto &kv : m_rev_var_map) {
        renaming_map[kv.first.name().str()] = kv.second.name().str();
      }
      m_alloc.add_renaming_map(renaming_map);
      o << m_base_dom;
      m_alloc.clear_renaming_map();
    }
  }
}; // class region_domain

template <typename Params>
struct abstract_domain_traits<region_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab
