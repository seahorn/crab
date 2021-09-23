#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/types.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/tag.hpp>
#include <crab/types/varname_factory.hpp>

#include <algorithm>
#include <functional>
#include <unordered_map>

////////////////////////////////////////////////////////////////////////
// Abstract domain for regions and references.
//
// This version improves region_domain.hpp in two ways:
//
// - it allows to reason about the is_dereferenceable intrinsics.
//
// - it can also be more precise with unknown regions by performing
//   type inference during the analysis. The type of an unknown region
//   can be narrowed to a more precise type by inspecting the region
//   writes.
////////////////////////////////////////////////////////////////////////

namespace crab {
namespace domains {
/**
 * Params must be defined similarly to region_domain_impl::Params
 * class in region_domain.hpp
 **/
template <typename Params>
class region_with_is_deref_domain final
    : public abstract_domain_api<region_with_is_deref_domain<Params>> {
  using region_domain_t = region_with_is_deref_domain<Params>;
  using abstract_domain_t = abstract_domain_api<region_domain_t>;

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
  using bound_t = typename interval_t::bound_t;
  using base_varname_allocator_t = typename Params::varname_allocator_t;

  // References can be modeled by more than one ghost
  // variable. Currently, we model three things about references: (1)
  // the address of the reference, (2) its offset (distance to the
  // base address of the memory object to which the reference points
  // to), and (3) the size of the memory object. The ghost variables
  // (2)-(3) are optional.
  //
  // To make the code of the region domain simpler and extensible, we
  // group together all the ghost variables in a class and provide
  // abstract operations so that the rest of the code is mostly
  // unaware of how many ghost variables the domain is keeping track
  // and how they should be updated.
  class ghost_variables_t {
  public:
    class object_offset_size_t {
      friend class ghost_variables_t;

      base_variable_t m_offset;
      base_variable_t m_size;

      object_offset_size_t(base_variable_t &&offset, base_variable_t &&size)
          : m_offset(std::move(offset)), m_size(std::move(size)) {
        if (!m_offset.get_type().is_integer() ||
            !m_size.get_type().is_integer()) {
          CRAB_ERROR("ghost offset and size variables must be integers");
        }
      }
      object_offset_size_t(base_variable_t offset, base_variable_t size)
          : m_offset(offset), m_size(size) {
        if (!m_offset.get_type().is_integer() ||
            !m_size.get_type().is_integer()) {
          CRAB_ERROR("ghost offset and size variables must be integers");
        }
      }

    public:
      base_variable_t get_offset() const { return m_offset; }

      base_variable_t get_size() const { return m_size; }

      void init(base_abstract_domain_t &dom,
                const base_variable_or_constant_t &size) const {
        dom.assign(m_offset, number_t(0));
        if (size.is_constant()) {
          dom.assign(m_size, size.get_constant());
        } else {
          dom.assign(m_size, size.get_variable());
        }
      }

      void assign(base_abstract_domain_t &dom,
                  const object_offset_size_t &rhs) const {
        if (m_offset.get_type() != rhs.m_offset.get_type()) {
          CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                     "variables");
        }
        if (m_size.get_type() != rhs.m_size.get_type()) {
          CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                     "variables");
        }

        dom.assign(m_offset, rhs.m_offset);
        dom.assign(m_size, rhs.m_size);
      }

      void assign(base_abstract_domain_t &dom, const object_offset_size_t &rhs,
                  const base_linear_expression_t &offset) const {
        if (m_offset.get_type() != rhs.m_offset.get_type()) {
          CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                     "variables");
        }
        if (m_size.get_type() != rhs.m_size.get_type()) {
          CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                     "variables");
        }

        dom.assign(m_offset, rhs.m_offset + offset);
        dom.assign(m_size, rhs.m_size);
      }

      bool is_deref(base_abstract_domain_t &dom,
                    const base_variable_or_constant_t &byte_sz) const {
        // TODO: we don't check that m_offset >= 0
        base_abstract_domain_t tmp(dom);
        if (byte_sz.is_constant()) {
          tmp += base_linear_constraint_t(m_offset + byte_sz.get_constant() >
                                          m_size);
        } else {
          tmp += base_linear_constraint_t(m_offset + byte_sz.get_variable() >
                                          m_size);
        }
        return tmp.is_bottom();
      }

      void expand(base_abstract_domain_t &dom,
                  const object_offset_size_t &rhs) const {
        if (m_offset.get_type() != rhs.m_offset.get_type()) {
          CRAB_ERROR("Type inconsistency while expand ghost object_offset_size "
                     "variables");
        }
        if (m_size.get_type() != rhs.m_size.get_type()) {
          CRAB_ERROR("Type inconsistency while expand ghost object_offset_size "
                     "variables");
        }

        dom.expand(m_offset, rhs.m_offset);
        dom.expand(m_size, rhs.m_size);
      }

      void forget(base_abstract_domain_t &dom) const {
        dom -= m_offset;
        dom -= m_size;
      }
    }; /* end class object_offset_size_t */

    // - if the variable is a reference then `m_var` is its symbolic
    //   address and its type is an integer.  note that the base domain
    //   is completely unaware of regions and references.
    // - if the variable is T=integer|real|bool or arrays of T then
    //   `m_var` is its value and the type is preserved.
    // - if the variable is a region of type T=integer|real|bool then
    //   `m_var` represents its contents and its type is T.
    // - if the variable is a region of type reference then `m_var`
    //   represents its content and its type is an integer.
    base_variable_t m_var;
    // if the variable is a reference or a region of references then
    // m_object_offset_size.first is the `offset`, the distance between its
    // address and the base address of the memory object to which
    // address belongs to, and m_object_offset_size.second is the `size` of the
    // memory object.
    boost::optional<object_offset_size_t> m_object_offset_size;
    // The CrabIR type of the variable which is being shadowed
    variable_type m_vty;

    void check_types() {
      if (m_var.get_type().is_region() || m_var.get_type().is_reference()) {
        CRAB_ERROR("Ghost variables cannot be regions or references");
      }

      if (m_object_offset_size) {
        if (!(*m_object_offset_size).get_offset().get_type().is_integer()) {
          CRAB_ERROR("ghost offset can only be integer");
        }
        if (!(*m_object_offset_size).get_size().get_type().is_integer()) {
          CRAB_ERROR("ghost size can only be integer");
        }
      }
    }

    ghost_variables_t(base_variable_t bv, variable_type vty)
        : m_var(bv), m_object_offset_size(boost::none), m_vty(vty) {
      check_types();
    }

    ghost_variables_t(base_variable_t bv, base_variable_t offset,
                      base_variable_t size, variable_type vty)
        : m_var(bv), m_object_offset_size(object_offset_size_t(offset, size)),
          m_vty(vty) {
      check_types();
    }

    // Create a fresh ghost variable in the base domain to shadow v
    static base_variable_t
    make_base_variable(base_varname_allocator_t &var_allocator,
                       const variable_type &vty) {
      if (vty.is_reference()) {
        base_variable_t bv(var_allocator.next(), crab::INT_TYPE,
                           32 /*should be defined in Params*/);
        return bv;
      } else if (vty.is_region()) {
        variable_type_kind ty;
        unsigned bitwidth = 0;
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
          CRAB_ERROR("make_base_variable: unreachable");
        }
        base_variable_t bv(var_allocator.next(), ty, bitwidth);
        return bv;
      } else {
        base_variable_t bv(var_allocator.next(), vty);
        return bv;
      }
    }

  public:
    static ghost_variables_t create(base_varname_allocator_t &alloc,
                                    const variable_type &vty, unsigned line) {
      if (vty.is_unknown_region()) {
        CRAB_ERROR("create should not be called at line ", line,
                   " with a type ", vty);
      }

      if (crab_domain_params_man::get().region_is_dereferenceable() &&
          (vty.is_reference() || vty.is_reference_region())) {
        return ghost_variables_t(make_base_variable(alloc, vty),
                                 make_base_variable(alloc, vty),
                                 make_base_variable(alloc, vty), vty);
      } else {
        return ghost_variables_t(make_base_variable(alloc, vty), vty);
      }
    }

    // Create a fresh set of ghost variables from x
    static ghost_variables_t create(base_varname_allocator_t &alloc,
                                    const ghost_variables_t &x) {
      if (x.has_offset_and_size()) {
        object_offset_size_t offset_size = x.get_offset_and_size();
        return ghost_variables_t(
            make_base_variable(alloc, x.get_var().get_type()),
            make_base_variable(alloc, offset_size.get_offset().get_type()),
            make_base_variable(alloc, offset_size.get_size().get_type()),
            x.m_vty);
      } else {
        return ghost_variables_t(
            make_base_variable(alloc, x.get_var().get_type()), x.m_vty);
      }
    }

    // Return true if the types of the ghost variables are
    // syntactically equal.
    bool same_type(const ghost_variables_t &o) const {
      // Check that the types of CrabIR variables are consistent
      if (!(m_vty.is_region() && o.m_vty.is_unknown_region() ||
            m_vty.is_unknown_region() && o.m_vty.is_region() ||
            m_vty == o.m_vty)) {
        return false;
      }

      // Check that the types of the ghost variables are syntactically
      // equal.
      if (get_var().get_type() != o.get_var().get_type()) {
        return false;
      }

      bool x_has_offset_and_size = has_offset_and_size();
      bool y_has_offset_and_size = o.has_offset_and_size();

      if (x_has_offset_and_size && y_has_offset_and_size) {
        object_offset_size_t x_offset_size = get_offset_and_size();
        object_offset_size_t y_offset_size = o.get_offset_and_size();
        return (x_offset_size.get_offset().get_type() ==
                    y_offset_size.get_offset().get_type() &&
                x_offset_size.get_size().get_type() ==
                    y_offset_size.get_size().get_type());
      }

      return (!x_has_offset_and_size && !y_has_offset_and_size);
    }

    // Store all ghost variables into a vector
    void add(base_variable_vector_t &vec) const {
      vec.push_back(m_var);
      if (m_object_offset_size) {
        vec.push_back((*m_object_offset_size).get_offset());
        vec.push_back((*m_object_offset_size).get_size());
      }
    }

    // Abstract domain operations
    void assign(base_abstract_domain_t &dom,
                const ghost_variables_t &rhs) const {
      if (m_var.get_type() != rhs.m_var.get_type()) {
        CRAB_ERROR("Type inconsistency while assign in ghost variables ",
                   m_var.get_type(), " and ", rhs.m_var.get_type());
      }

      if (m_var.get_type().is_bool()) {
        dom.assign_bool_var(m_var, rhs.m_var, false);
      } else {
        dom.assign(m_var, rhs.m_var);
      }

      if (m_object_offset_size && rhs.m_object_offset_size) {
        (*m_object_offset_size).assign(dom, *(rhs.m_object_offset_size));
      }
    }

    void expand(base_abstract_domain_t &dom,
                const ghost_variables_t &new_gvars) const {
      if (m_var.get_type() != new_gvars.m_var.get_type()) {
        CRAB_ERROR("Type inconsistency while expand in ghost variables ",
                   m_var.get_type(), " and ", new_gvars.m_var.get_type());
      }
      dom.expand(m_var, new_gvars.m_var);
      if (m_object_offset_size && new_gvars.m_object_offset_size) {
        (*m_object_offset_size).expand(dom, *(new_gvars.m_object_offset_size));
      }
    }

    void forget(base_abstract_domain_t &dom) const {
      dom -= m_var;
      if (m_object_offset_size) {
        (*m_object_offset_size).forget(dom);
      }
    }

    base_variable_t get_var() const { return m_var; }

    bool has_offset_and_size() const { return (bool)m_object_offset_size; }

    // address: 1, offset: 2, size: 3
    using ghost_var_id = unsigned;

    // Add ghost variables in the reverse map. A reverse map (rev_map)
    // maps ghost variables to CrabIR variables.
    void update_rev_varmap(
        std::unordered_map<base_variable_t, std::pair<variable_t, ghost_var_id>>
            &rev_map,
        variable_t v) const {
      rev_map.insert({m_var, {v, 1}});
      if (m_object_offset_size) {
        rev_map.insert({(*m_object_offset_size).get_offset(), {v, 2}});
        rev_map.insert({(*m_object_offset_size).get_size(), {v, 3}});
      }
    }

    // Remove ghost variables from the reverse map. A reverse map
    // (rev_map) maps ghost variables to CrabIR variables.
    void remove_rev_varmap(
        std::unordered_map<base_variable_t, std::pair<variable_t, ghost_var_id>>
            &rev_map) {
      rev_map.erase(m_var);
      if (m_object_offset_size) {
        rev_map.erase((*m_object_offset_size).get_offset());
        rev_map.erase((*m_object_offset_size).get_size());
      }
    }

    // Assign a string name to each ghost variable
    static void mk_renaming_map(
        const std::unordered_map<base_variable_t,
                                 std::pair<variable_t, ghost_var_id>> &rev_map,
        std::function<variable_type(const variable_t &)> get_type_fn,
        std::unordered_map<std::string, std::string> &out_str_map) {
      for (auto &kv : rev_map) {
        base_variable_t bv = kv.first;
        variable_t v = kv.second.first;
        ghost_var_id ghost_id = kv.second.second;

        variable_type vty = get_type_fn(v);
        if (crab_domain_params_man::get().region_is_dereferenceable() &&
	    (vty.is_reference() || vty.is_reference_region())) {
          if (ghost_id == 1) {
            out_str_map[bv.name().str()] = v.name().str() + ".address";
          } else if (ghost_id == 2) {
            out_str_map[bv.name().str()] = v.name().str() + ".offset";
          } else if (ghost_id == 3) {
            out_str_map[bv.name().str()] = v.name().str() + ".size";
          }
        } else {
          out_str_map[bv.name().str()] = v.name().str();
        }
      }
    }

    object_offset_size_t get_offset_and_size() const {
      if (!m_object_offset_size) {
        CRAB_ERROR(
            "ghost variables for offset and size have not been allocated");
      }
      return *m_object_offset_size;
    }

    void write(crab_os &o) const {
      if (!m_object_offset_size) {
        o << m_var;
      } else {
        o << "{reference:" << m_var
          << ",offset:" << (*m_object_offset_size).get_offset()
          << ",size:" << (*m_object_offset_size).get_size() << "}";
      }
    }

    friend crab_os &operator<<(crab_os &o, const ghost_variables_t &gvars) {
      gvars.write(o);
      return o;
    }
  }; /* end class ghost_variables_t */

  // A simple class for tags (i.e., numerical identifiers). We don't
  // use crab::tag because we want to have the flexibility of creating
  // tags without a tag manager.
  class tag_t : public indexable {
    ikos::index_t m_id;

  public:
    tag_t(number_t n) : m_id(0) {
      if (n < 0) {
        CRAB_ERROR("Cannot use negative numbers for tags");
      }
      if (!n.fits_int64()) {
        CRAB_ERROR("Too large value for a tag");
      }
      m_id = (int64_t)n;
    }
    bool operator<(const tag_t &as) const { return m_id < as.m_id; }
    bool operator==(const tag_t &as) const { return m_id == as.m_id; }
    virtual ikos::index_t index() const override { return m_id; }
    void write(crab_os &o) const override { o << "TAG_" << m_id; }
    friend crab_os &operator<<(crab_os &o, const tag_t &as) {
      as.write(o);
      return o;
    }
  }; /* end class tag_t */

  /**------------- Begin type definitions -------------------**/

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
  // Map each reference or region variable to a set of allocation sites
  using alloc_site_env_t =
      separate_discrete_domain<variable_t, allocation_site>;
  using allocation_sites = typename alloc_site_env_t::value_type;
  // Map variables to sets of tags
  using tag_env_t = separate_discrete_domain<variable_t, tag_t>;
  using tag_set = typename tag_env_t::value_type;
  /**------------- End type definitions -------------------**/

  bool m_is_bottom; // special symbol for bottom
  /** Begin base domain **/
  // To create ghost variables for the base domain.
  base_varname_allocator_t m_alloc;
  // Map a variable_t to its ghost variables:
  var_map_t m_var_map;
  // Reverse map from ghost variables to variable_t
  rev_var_map_t m_rev_var_map;
  // The base abstract domain (over ghost variables): all the heavy
  // lifting is done here.  m_base_dom does not have any variable of
  // region or reference type.
  base_abstract_domain_t m_base_dom;
  /** End base domain **/

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

  region_with_is_deref_domain(
      base_varname_allocator_t &&alloc, var_map_t &&var_map,
      rev_var_map_t &&rev_var_map, rgn_counting_env_t &&rgn_counting_dom,
      rgn_bool_env_t &&rgn_init_dom, rgn_type_env_t &&rgn_type_dom,
      base_abstract_domain_t &&base_dom, alloc_site_env_t &&alloc_site_dom,
      rgn_dealloc_t &&rgn_dealloc_dom, tag_env_t &&tag_env)
      : m_is_bottom(base_dom.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_base_dom(std::move(base_dom)),
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
        m_rgn_type_dom(std::move(rgn_type_dom)),
        m_alloc_site_dom(std::move(alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
        m_tag_env(std::move(tag_env)) {}

  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;

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
  void do_join(const region_domain_t &right) {
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(m_rgn_type_dom | right.m_rgn_type_dom);
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
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
	    !kv.second.same_type(it->second)) {
          crab::CrabStats::count(
              domain_name() + ".count.join.skipped.inconsistent_dynamic_type");
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
    std::swap(m_rgn_type_dom, out_rgn_type_dom);
    std::swap(m_alloc_site_dom, out_alloc_site_dom);
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_tag_env, out_tag_env);

    CRAB_LOG("region", crab::outs() << *this << "\n");
  }

  region_domain_t do_join_or_widening(const region_domain_t &left,
                                      const region_domain_t &right,
                                      const bool is_join /*unused*/,
                                      base_dom_binop_t base_dom_op) const {

    // rgn_counting_dom does not require common renaming (i.e., no
    // ghost variables).  The domain is finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    // rgn_init_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom | right.m_rgn_init_dom);

    // The domain is finite
    rgn_type_env_t out_rgn_type_dom(left.m_rgn_type_dom | right.m_rgn_type_dom);

    // alloc_site_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom |
                                        right.m_alloc_site_dom);

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
    left_dom.project(left_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.project(right_vars);
    right_dom.rename(right_vars, out_vars);

    // Final join or widening
    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    region_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_base_dom), std::move(out_alloc_site_dom),
        std::move(out_rgn_dealloc_dom), std::move(out_tag_env));
    return res;
  }

  region_domain_t do_meet_or_narrowing(const region_domain_t &left,
                                       const region_domain_t &right,
                                       const bool is_meet /*unused*/,
                                       base_dom_binop_t base_dom_op) const {

    // these domains are finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom &
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom & right.m_rgn_init_dom);
    rgn_type_env_t out_rgn_type_dom(left.m_rgn_type_dom & right.m_rgn_type_dom);
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom &
                                        right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom &
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(left.m_tag_env & right.m_tag_env);

    // This shouldn't happen but just in case ...
    if (out_rgn_counting_dom.is_bottom() || out_rgn_init_dom.is_bottom() ||
        out_rgn_type_dom.is_bottom() || out_rgn_dealloc_dom.is_bottom() ||
        out_alloc_site_dom.is_bottom()) {
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
        std::move(out_rgn_init_dom), std::move(out_rgn_type_dom),
        std::move(out_base_dom), std::move(out_alloc_site_dom),
        std::move(out_rgn_dealloc_dom), std::move(out_tag_env));
    return res;
  }

  void ref_store(base_abstract_domain_t &base_dom, const variable_t &rgn_var,
                 const variable_or_constant_t &val) {
    assert(is_tracked_region(rgn_var));
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn_var);

    // At this point region of references has been translated to region of
    // integers. val can be any scalar including a reference.
    if (val.get_type().is_bool()) {
      if (val.is_constant()) {
        base_dom.assign_bool_cst(rgn_gvars.get_var(),
                                 (val.is_bool_true()
                                      ? base_linear_constraint_t::get_true()
                                      : base_linear_constraint_t::get_false()));
      } else {
        base_dom.assign_bool_var(
            rgn_gvars.get_var(),
            get_or_insert_gvars(val.get_variable()).get_var(), false);
      }
    } else if (val.get_type().is_integer() || val.get_type().is_real() ||
               val.get_type().is_reference()) {

      if (val.is_constant()) {
        base_dom.assign(rgn_gvars.get_var(), val.get_constant());

        if (val.get_type().is_reference()) {
          if (rgn_gvars.has_offset_and_size()) {
            // Storing NULL
            assert(val.get_constant() == number_t(0));
            rgn_gvars.get_offset_and_size().forget(base_dom);
          }
        }
      } else {
        ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
        base_dom.assign(rgn_gvars.get_var(), val_gvars.get_var());

        if (val.get_type().is_reference()) {
          if (rgn_gvars.has_offset_and_size() &&
              val_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().assign(
                base_dom, val_gvars.get_offset_and_size());
          } else if (rgn_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().forget(base_dom);
          }
        }
      }
    } else {
      CRAB_ERROR(domain_name(), "::ref_store: unsupported type ",
                 val.get_type());
    }
  }

  /*
   * To rename variables (variable_t) to ghost variables
   * (base_variable_t) so that they can be used in the base abstract
   * domain.
   */

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
  base_linear_expression_t rename_linear_expr(const linear_expression_t &e) {
    base_linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      base_variable_t nv = get_or_insert_gvars(v).get_var();
      out = out + (coef * nv);
    }
    return out;
  }

  // Rename linear constraint
  base_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst) {
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
  rename_linear_cst_sys(const linear_constraint_system_t &csts) {
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
    return (crab_domain_params_man::get().region_skip_unknown_regions() ?
	    false :
	    has_dynamic_type(v));
	    
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

  region_with_is_deref_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  region_with_is_deref_domain(const region_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_alloc(o.m_alloc), m_var_map(o.m_var_map),
        m_rev_var_map(o.m_rev_var_map), m_base_dom(o.m_base_dom),
        m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom), m_rgn_type_dom(o.m_rgn_type_dom),
        m_alloc_site_dom(o.m_alloc_site_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom), m_tag_env(o.m_tag_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  region_with_is_deref_domain(region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_base_dom(std::move(o.m_base_dom)),
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
        m_rgn_type_dom(std::move(o.m_rgn_type_dom)),
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
      m_rgn_type_dom = o.m_rgn_type_dom;
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
      m_rgn_type_dom = std::move(o.m_rgn_type_dom);
      m_alloc_site_dom = std::move(o.m_alloc_site_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_tag_env = std::move(o.m_tag_env);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    bool res =
        (!is_bottom() && m_base_dom.is_top() && m_rgn_counting_dom.is_top());
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

    if (!(m_rgn_type_dom <= o.m_rgn_type_dom)) {
      CRAB_LOG("region-leq", crab::outs() << "Result5=0\n";);
      return false;
    }

    if (crab_domain_params_man::get().region_allocation_sites()) {
      if (!(m_alloc_site_dom <= o.m_alloc_site_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result6=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_deallocation()) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result7=0\n";);
        return false;
      }
    }

    if (crab_domain_params_man::get().region_tag_analysis()) {
      if (!(m_tag_env <= o.m_tag_env)) {
        CRAB_LOG("region-leq", crab::outs() << "Result8=0\n";);
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
    left_dom.forget(out_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.forget(out_vars);
    right_dom.rename(right_vars, out_vars);

    CRAB_LOG("region-leq", crab::outs()
                               << "Inclusion test (after renaming):\n\t"
                               << left_dom << "\n\t" << right_dom << "\n";);

    bool res = left_dom <= right_dom;
    CRAB_LOG("region-leq", crab::outs() << "Result9=" << res << "\n";);
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
      region_domain_t left(*this);
      refine_regions(left_regions, regions_right_base_vars, left.m_var_map,
                     left.m_rev_var_map, left.m_alloc, left.m_base_dom,
                     o.m_base_dom, o.m_alloc);
      region_domain_t res(std::move(
          do_join_or_widening(left, o, true /*is join*/, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (!refine_left && refine_right) {
      // Refine right by propagating information from left's regions
      region_domain_t right(o);
      refine_regions(right_regions, regions_left_base_vars, right.m_var_map,
                     right.m_rev_var_map, right.m_alloc, right.m_base_dom,
                     m_base_dom, m_alloc);
      region_domain_t res(std::move(
          do_join_or_widening(*this, right, true /*is join*/, base_dom_op)));
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
      region_domain_t res(std::move(
          do_join_or_widening(left, right, true /*is join*/, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else {
      region_domain_t res(std::move(
          do_join_or_widening(*this, o, true /*is join*/, base_dom_op)));

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
    region_domain_t res(std::move(
        do_meet_or_narrowing(*this, o, true /*is meet*/, base_dom_op)));

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
    region_domain_t res(std::move(
        do_join_or_widening(*this, o, false /*is widen*/, base_dom_op)));

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

    region_domain_t res(std::move(
        do_join_or_widening(*this, o, false /*is widen*/, base_dom_op)));

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
    region_domain_t res(std::move(
        do_meet_or_narrowing(*this, o, false /*is narrow*/, base_dom_op)));

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
        if (is_tracked_unknown_region(v)) {
          m_rgn_type_dom -= v;
        }
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
        it->second.forget(m_base_dom);
        it->second.remove_rev_varmap(m_rev_var_map);
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

    if (is_tracked_unknown_region(rgn)) {
      m_rgn_type_dom.set(rgn, rgn.get_type());
    }

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

    if (!(rgn.get_type().is_unknown_region())) {
      // Assign ghost variables to rgn for modeling its content
      ghost_variables_t::create(m_alloc, rgn.get_type(), __LINE__);
      crab::CrabStats::count(domain_name() +
                             ".count.region_init.typed_regions");
    } else {
      // We cannot assign ghost variables if the region is unknown
      crab::CrabStats::count(domain_name() +
                             ".count.region_init.unknown_regions");
    }
    CRAB_LOG("region", crab::outs()
                           << "After region_init(" << rgn << ":"
                           << rgn.get_type() << ")=" << *this << "\n";);
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
    if (is_tracked_unknown_region(lhs_rgn)) {
      m_rgn_type_dom.set(lhs_rgn, m_rgn_type_dom[rhs_rgn]);
    }
    if (crab_domain_params_man::get().region_allocation_sites()) {
      m_alloc_site_dom.set(lhs_rgn, m_alloc_site_dom[rhs_rgn]);
    }
    if (crab_domain_params_man::get().region_deallocation()) {
      m_rgn_dealloc_dom.add(rhs_rgn, lhs_rgn);
    }
    if (crab_domain_params_man::get().region_tag_analysis()) {
      m_tag_env.set(lhs_rgn, m_tag_env[rhs_rgn]);
    }

    if (!is_tracked_region(rhs_rgn)) {
      return;
    }

    const ghost_variables_t &base_lhs = get_or_insert_gvars(lhs_rgn);
    const ghost_variables_t &base_rhs = get_or_insert_gvars(rhs_rgn);
    auto num_refs = m_rgn_counting_dom[rhs_rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      // the rhs region is a singleton then we use assign
      base_lhs.assign(m_base_dom, base_rhs);
    } else {
      // the rhs region might not be a singleton so we use expand.
      base_lhs.forget(m_base_dom);
      base_rhs.expand(m_base_dom, base_lhs);
    }
    CRAB_LOG("region", crab::outs()
                           << "After region_copy(" << lhs_rgn << ":"
                           << lhs_rgn.get_type() << "," << rhs_rgn << ":"
                           << rhs_rgn.get_type() << ")=" << *this << "\n";);
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

    auto assign_regions = [this](const variable_t &x, const variable_t &y,
                                 const variable_type &ty) {
      assert(ty.is_region());
      assert(!ty.is_unknown_region());

      if (!is_tracked_region(x) || !is_tracked_region(y)) {
        crab::CrabStats::count(domain_name() + ".count.region_cast.skipped");
        return;
      }

      ghost_variables_t x_gvars = get_or_insert_gvars(x);
      ghost_variables_t y_gvars = get_or_insert_gvars(y);

      if (ty.is_reference_region()) {
        // reference modeled as integer
        m_base_dom.assign(x_gvars.get_var(), y_gvars.get_var());
        if (x_gvars.has_offset_and_size() && y_gvars.has_offset_and_size()) {
          x_gvars.get_offset_and_size().assign(m_base_dom,
                                               y_gvars.get_offset_and_size());
        }
      } else if (ty.is_bool_region()) {
        m_base_dom.assign_bool_var(x_gvars.get_var(), y_gvars.get_var(), false);
      } else if (ty.is_integer_region() || ty.is_real_region()) {
        m_base_dom.assign(x_gvars.get_var(), y_gvars.get_var());
      } else if (ty.is_array_region()) {
        m_base_dom.array_assign(x_gvars.get_var(), y_gvars.get_var());
      }
    };

    bool unknown_dst = is_tracked_unknown_region(dst_rgn);
    bool unknown_src = is_tracked_unknown_region(src_rgn);

    if (!unknown_src && unknown_dst) {
      m_rgn_type_dom.set(dst_rgn, src_rgn.get_type());
      assign_regions(dst_rgn, src_rgn, src_rgn.get_type());
    } else if (unknown_src && !unknown_dst) {
      type_value src_type = m_rgn_type_dom[src_rgn];
      if (type_value(dst_rgn.get_type()) <= src_type) {
        assign_regions(dst_rgn, src_rgn, dst_rgn.get_type());
      } else {
        crab::CrabStats::count(
            domain_name() +
            ".count.region_cast.skipped.inconsistent_dynamic_type");
      }
    } else {
      // We can be here either (1) both regions are typed, or (2)
      // crab_domain_params_man::get().region_skip_unknown_regions() = 1.
      ///
      // We shouldn't use region_cast if both regions are typed so we
      // do nothing although we should print a warning or something.
      crab::CrabStats::count(domain_name() + ".count.region_cast.skipped");
    }

    CRAB_LOG("region", crab::outs()
                           << "After region_cast(" << src_rgn << ":"
                           << src_rgn.get_type() << "," << dst_rgn << ":"
                           << dst_rgn.get_type() << ")=" << *this << "\n";);
  }

  // Create a new reference ref to region rgn.
  void ref_make(const variable_t &ref, const variable_t &rgn,
                const variable_or_constant_t &size,
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

    // Assign ghost variables to ref
    ghost_variables_t ref_gvars = get_or_insert_gvars(ref);

    // initialize ghost variables
    if (ref_gvars.has_offset_and_size()) {
      ref_gvars.get_offset_and_size().init(m_base_dom,
                                           rename_variable_or_constant(size));
    }

    CRAB_LOG("region", crab::outs() << "After ref_make(" << ref << "," << rgn
                                    << ":" << rgn.get_type() << "," << size
                                    << "," << as << ")=" << *this << "\n";);
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
        (rgn.get_type().is_integer_region() && !res.get_type().is_integer()) ||
        (rgn.get_type().is_real_region() && !res.get_type().is_real()) ||
        (rgn.get_type().is_reference_region() &&
         !res.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::ref_load: type of lhs ", res, " (",
                 res.get_type(), ") is not compatible with region ", rgn, " (",
                 rgn.get_type(), ")");
    }

    if (is_bottom()) {
      return;
    }

    const ghost_variables_t &gvars_res = get_or_insert_gvars(res);

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-load", CRAB_WARN("region_domain::ref_load: reference ",
                                        ref, " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      gvars_res.forget(m_base_dom);
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

    if (!is_tracked_region(rgn)) {
      CRAB_LOG(
          "region-load",
          CRAB_WARN("region_domain::ref_load: skip unknown region ", rgn););
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.untracked_region");
      gvars_res.forget(m_base_dom);
      return;
    }

    if (is_tracked_unknown_region(rgn)) {
      type_value rgn_ty = m_rgn_type_dom[rgn];
      if (rgn_ty.is_bottom()) {
        // this shouldn't happen because it can only happen if
        // m_rgn_type_dom is bottom and this is already checked.
        CRAB_ERROR(domain_name(), "::ref_load: the dynamic type of region ",
                   rgn, " is bottom");
      }

      if (rgn_ty.is_top()) {
        gvars_res.forget(m_base_dom);
        crab::CrabStats::count(domain_name() +
                               ".count.ref_load.skipped.dynamic_type_is_top");
        return;
      }

      variable_type rgn_content_ty = rgn_ty.get().get_region_content_type();
      if (res.get_type() != rgn_content_ty) {
        CRAB_LOG("region-load",
                 CRAB_WARN("region_domain::ref_load: skip unknown region ", rgn,
                           " because the region has the dynamic type ", rgn_ty,
                           " which is not compatible with the type of the lhs ",
                           res.get_type()););

        crab::CrabStats::count(
            domain_name() +
            ".count.ref_load.skipped.inconsistent_dynamic_type");
        gvars_res.forget(m_base_dom);
        return;
      }
    }

    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      // strong read
      CRAB_LOG("region-load", crab::outs() << "Reading from singleton\n";);
      if (auto region_gvars_opt = get_gvars(rgn)) {
        gvars_res.assign(m_base_dom, *region_gvars_opt);
      } else {
        gvars_res.forget(m_base_dom);
      }
    } else {
      // weak read
      CRAB_LOG("region-load", crab::outs() << "Reading from non-singleton\n";);
      if (auto region_gvars_opt = get_gvars(rgn)) {
        variable_type rgn_ty = get_dynamic_type_or_fail(rgn);
        ghost_variables_t fresh_region_gvars =
            ghost_variables_t::create(m_alloc, (*region_gvars_opt));
        (*region_gvars_opt).expand(m_base_dom, fresh_region_gvars);
        gvars_res.assign(m_base_dom, fresh_region_gvars);
      } else {
        gvars_res.forget(m_base_dom);
      }
    }

    CRAB_LOG("region-load", crab::outs()
                                << "After " << res << ":="
                                << "ref_load(" << rgn << ":" << rgn.get_type()
                                << "," << ref << ":" << ref.get_type()
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
        (rgn.get_type().is_integer_region() && !val.get_type().is_integer()) ||
        (rgn.get_type().is_real_region() && !val.get_type().is_real()) ||
        (rgn.get_type().is_reference_region() &&
         !val.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::ref_store: type of value ", val, " (",
                 val.get_type(), ") is not compatible with region ", rgn, " (",
                 rgn.get_type(), ")");
    }

    if (is_bottom()) {
      return;
    }

    // forget the ghost variables in the base domain
    auto forget_region = [this](const variable_t &rgn) {
      if (boost::optional<ghost_variables_t> gvars = get_gvars(rgn)) {
        (*gvars).forget(m_base_dom);
      }
    };

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-store", CRAB_WARN("region_domain::ref_store: reference ",
                                         ref, " is null. "););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_store.skipped.null_reference");
      forget_region(rgn);
      return;
    }

    bool is_uninitialized_rgn = m_rgn_init_dom[rgn].is_false();
    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    if (is_tracked_unknown_region(rgn)) {
      type_value rgn_ty = m_rgn_type_dom[rgn];
      if (rgn_ty.is_bottom()) {
        // this shouldn't happen because it can only happen if
        // m_rgn_type_dom is bottom and this is already checked.
        CRAB_ERROR(domain_name(), "::ref_store: the dynamic type of region ",
                   rgn, " is bottom");
      }
      if (rgn_ty.is_top()) {
        crab::CrabStats::count(domain_name() +
                               ".count.ref_store.skipped.dynamic_type_is_top");
        forget_region(rgn);
        return;
      } else if (rgn_ty.get().is_unknown_region()) {
        // 1. Set the dynamic type of the unknown region. From now on,
        // all memory accesses must satisfy this type.
        m_rgn_type_dom.set(rgn, variable_type::mk_region(val.get_type()));
      } else {
        // 2. Check that type of val satisfy the dynamic type of the
        // region.
        variable_type rgn_content_ty = rgn_ty.get().get_region_content_type();
        if (val.get_type() != rgn_content_ty) {
          if (rgn_content_ty.is_reference()) {
            // Skip the store
            //
            // IMPORTANT: we are trying to write some non-reference
            // value into a region where we wrote a reference
            // before. The reason why this is possible is due to
            // pointer analysis imprecision.
            //
            // Skipping the store is sound because any future read of
            // any non-reference value will be soundly ignored (under
            // word-level addressability assumption)
            CRAB_LOG("region-store",
                     crab::outs() << "Skipped the write of a value of type "
                                  << val.get_type() << "\n";);
            crab::CrabStats::count(
                domain_name() +
                ".count.ref_store.skipped.inconsistent_dynamic_type");
            return;
          } else if (val.get_type().is_reference()) {
            // Reinterpret the region
            //
            // IMPORTANT: we throw away previous writes of
            // non-reference values in the region and from now on, we
            // only stores references.
            //
            // The argument to make uninitialized the region is that
            // previous writes must have happened in different
            // locations because those writes wrote values of
            // different types.
            m_rgn_init_dom.set(rgn, boolean_value::get_false());
            m_rgn_type_dom.set(rgn, variable_type::mk_region(val.get_type()));
            if (boost::optional<ghost_variables_t> gvars = get_gvars(rgn)) {
              m_var_map.erase(rgn);
              (*gvars).remove_rev_varmap(m_rev_var_map);
              (*gvars).forget(m_base_dom);
            }
            CRAB_LOG("region-store",
                     crab::outs() << "Reinterpreting the dynamic type of "
                                  << rgn << " to " << val.get_type() << "\n";);
            // We don't return
          } else {
            CRAB_LOG(
                "region-store",
                CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn,
                          " because the region has the dynamic type ", rgn_ty,
                          " and it is not compatible with the type of the lhs ",
                          val.get_type()););
            crab::CrabStats::count(
                domain_name() +
                ".count.ref_store.skipped.inconsistent_dynamic_type");
            forget_region(rgn);
            return;
          }
        }
      }
    }

    bool is_tracked_rgn = is_tracked_region(rgn);
    if (!is_tracked_rgn) {
      CRAB_LOG(
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      crab::CrabStats::count(domain_name() +
                             ".count.ref_store.skipped.untracked_region");
      // don't return yet
    }

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

      if (is_tracked_rgn) {
        ref_store(m_base_dom, rgn, val);
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

      if (is_tracked_rgn) {
        base_abstract_domain_t tmp(m_base_dom);
        ref_store(tmp, rgn, val);
        m_base_dom |= tmp;
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
    CRAB_LOG("region-store",
             crab::outs() << "After ref_store(" << rgn << ":" << rgn.get_type()
                          << "," << ref << ":" << ref.get_type() << "," << val
                          << ":" << val.get_type() << ")=" << *this << "\n";);
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
        r += p.first *
             m_base_dom.operator[](get_or_insert_gvars(p.second).get_var());
      }
      return r;
    };

    if (!is_bottom()) {
      auto boffset = rename_linear_expr(offset);

      ghost_variables_t ref1_gvars = get_or_insert_gvars(ref1);
      ghost_variables_t ref2_gvars = get_or_insert_gvars(ref2);

      m_base_dom.assign(ref2_gvars.get_var(), ref1_gvars.get_var() + boffset);

      if (ref1_gvars.has_offset_and_size() &&
          ref2_gvars.has_offset_and_size()) {
        ref2_gvars.get_offset_and_size().assign(
            m_base_dom, ref1_gvars.get_offset_and_size(), boffset);
      } else if (ref2_gvars.has_offset_and_size()) {
        ref2_gvars.get_offset_and_size().forget(m_base_dom);
      }

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
          // m_rgn_dealloc_dom.remove_equiv_class(rgn2);
        } else if (found1 && !found2) {
          // m_rgn_dealloc_dom.remove_equiv_class(rgn1);
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

    base_variable_t base_lhs = get_or_insert_gvars(lhs).get_var();

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

    if (boost::optional<ghost_variables_t> arr_gvars_opt = get_gvars(rgn)) {
      auto num_refs = m_rgn_counting_dom[rgn];
      auto b_elem_size = rename_linear_expr(elem_size);
      if (num_refs.is_zero() || num_refs.is_one()) {
        CRAB_LOG("region", crab::outs() << "Reading from singleton\n";);
        auto b_index = rename_linear_expr(index);
        m_base_dom.array_load(base_lhs, (*arr_gvars_opt).get_var(), b_elem_size,
                              b_index);
      } else {
        // TODO: get the bitwidth from index
        base_variable_t unknown_index(m_alloc.next(), crab::INT_TYPE, 32);
        m_base_dom.array_load(base_lhs, (*arr_gvars_opt).get_var(), b_elem_size,
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
    base_variable_t arr_var = get_or_insert_gvars(rgn).get_var();
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
      ghost_variables_t src_gvars = get_or_insert_gvars(ref_var);
      ghost_variables_t dst_gvars = get_or_insert_gvars(int_var);
      m_base_dom.assign(dst_gvars.get_var(), src_gvars.get_var());

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
      ghost_variables_t src_gvars = get_or_insert_gvars(int_var);
      ghost_variables_t dst_gvars = get_or_insert_gvars(ref_var);

      m_base_dom.assign(dst_gvars.get_var(), src_gvars.get_var());

      if (dst_gvars.has_offset_and_size()) {
        dst_gvars.get_offset_and_size().forget(m_base_dom);
      }

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
      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());
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
      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), k);

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
      m_base_dom.assign(get_or_insert_gvars(x).get_var(), b_e);

      if (crab_domain_params_man::get().region_tag_analysis()) {
        merge_tags(x, e.variables());
      }
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    if (!is_bottom()) {
      auto b_e1 = rename_linear_expr(e1);
      auto b_e2 = rename_linear_expr(e2);
      auto b_cond = rename_linear_cst(cond);
      m_base_dom.select(get_or_insert_gvars(lhs).get_var(), b_cond, b_e1, b_e2);

      if (crab_domain_params_man::get().region_tag_analysis()) {
        tag_set tags = tag_set::bottom();
        for (auto const &v : e1.variables()) {
          tags = tags | m_tag_env[v];
        }
        for (auto const &v : e2.variables()) {
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
      m_base_dom.backward_assign(get_or_insert_gvars(x).get_var(), b_e,
                                 invariant.m_base_dom);

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
      m_base_dom.backward_apply(op, get_or_insert_gvars(x).get_var(),
                                get_or_insert_gvars(y).get_var(), z,
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
      m_base_dom.backward_apply(op, get_or_insert_gvars(x).get_var(),
                                get_or_insert_gvars(y).get_var(),
                                get_or_insert_gvars(z).get_var(),
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
        m_base_dom.set(get_or_insert_gvars(x).get_var(), intv);
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
      return m_base_dom[get_or_insert_gvars(x).get_var()];
    }
  }

  // int cast operations and bitwise operations
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, get_or_insert_gvars(dst).get_var(),
                       get_or_insert_gvars(src).get_var());

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
      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());

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
      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), z);

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
      m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);

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
                auto false_cst = base_linear_constraint_t::get_false();
                m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                           false_cst);
              } else {
                assert(rhs.is_disequality());
                auto true_cst = base_linear_constraint_t::get_true();
                m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                           true_cst);
              }
              return;
            }
          }
        }
      }

      auto b_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);

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
      m_base_dom.assign_bool_var(get_or_insert_gvars(lhs).get_var(),
                                 get_or_insert_gvars(rhs).get_var(),
                                 is_not_rhs);

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
      m_base_dom.apply_binary_bool(op, get_or_insert_gvars(x).get_var(),
                                   get_or_insert_gvars(y).get_var(),
                                   get_or_insert_gvars(z).get_var());

      if (crab_domain_params_man::get().region_tag_analysis()) {
        m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }
    }
  }
  void assume_bool(const variable_t &v, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_bottom()) {
      m_base_dom.assume_bool(get_or_insert_gvars(v).get_var(), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      m_base_dom.select_bool(get_or_insert_gvars(lhs).get_var(),
                             get_or_insert_gvars(cond).get_var(),
                             get_or_insert_gvars(b1).get_var(),
                             get_or_insert_gvars(b2).get_var());

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
      m_base_dom.backward_assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                          b_rhs, invariant.m_base_dom);
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
      m_base_dom.backward_assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                          b_rhs, invariant.m_base_dom);
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
      m_base_dom.backward_assign_bool_var(get_or_insert_gvars(lhs).get_var(),
                                          get_or_insert_gvars(rhs).get_var(),
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
      m_base_dom.backward_apply_binary_bool(
          op, get_or_insert_gvars(x).get_var(),
          get_or_insert_gvars(y).get_var(), get_or_insert_gvars(z).get_var(),
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
      m_base_dom.array_init(get_or_insert_gvars(a).get_var(), b_elem_size,
                            b_lb_idx, b_ub_idx, b_val);
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
      m_base_dom.array_load(get_or_insert_gvars(lhs).get_var(),
                            get_or_insert_gvars(a).get_var(), b_elem_size, b_i);
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
      m_base_dom.array_store(get_or_insert_gvars(a).get_var(), b_elem_size, b_i,
                             b_v, is_strong_update);
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
      m_base_dom.array_store_range(get_or_insert_gvars(a).get_var(),
                                   b_elem_size, b_i, b_j, b_v);
      if (crab_domain_params_man::get().region_tag_analysis()) {
        // TODO tag analysis
      }
    }
  }
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".array_assign");

    if (!is_bottom()) {
      m_base_dom.array_assign(get_or_insert_gvars(lhs).get_var(),
                              get_or_insert_gvars(rhs).get_var());
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
      m_base_dom.backward_array_init(get_or_insert_gvars(a).get_var(),
                                     b_elem_size, b_lb_idx, b_ub_idx, b_val,
                                     invariant.m_base_dom);
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
      m_base_dom.backward_array_load(get_or_insert_gvars(lhs).get_var(),
                                     get_or_insert_gvars(a).get_var(),
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
      m_base_dom.backward_array_store(get_or_insert_gvars(a).get_var(),
                                      b_elem_size, b_i, b_v, is_strong_update,
                                      invariant.m_base_dom);
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
      m_base_dom.backward_array_store_range(get_or_insert_gvars(a).get_var(),
                                            b_elem_size, b_i, b_j, b_v,
                                            invariant.m_base_dom);
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
      m_base_dom.backward_array_assign(get_or_insert_gvars(lhs).get_var(),
                                       get_or_insert_gvars(rhs).get_var(),
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
      if (v.get_type().is_region() && !is_tracked_region(v)) {
        continue;
      }
      get_or_insert_gvars(v).add(base_vars);
    }
    m_base_dom.project(base_vars);

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
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      zipped_vars.push_back({from[i], to[i]});
    }

    auto compare_zipped_pair = [](const std::pair<variable_t, variable_t> &p1,
                                  const std::pair<variable_t, variable_t> &p2) {
      return p1.first < p2.first;
    };

    std::sort(zipped_vars.begin(), zipped_vars.end(), compare_zipped_pair);

    // Needed because we cannot insert/delete while iterating over
    // m_var_map
    std::vector<variable_t> var_map_to_remove;
    std::vector<std::pair<variable_t, ghost_variables_t>> var_map_to_insert;

    for (auto &kv : m_var_map) {
      variable_t old_var = kv.first;
      ghost_variables_t old_gvars = kv.second;

      std::pair<variable_t, variable_t> p(old_var, old_var);
      auto lower = std::lower_bound(zipped_vars.begin(), zipped_vars.end(), p,
                                    compare_zipped_pair);

      if (!(lower == zipped_vars.end() || old_var < (*lower).first)) {
        // found
        variable_t new_var = (*lower).second;
        ghost_variables_t new_base_gvars = get_or_insert_gvars(new_var);
        ;

        var_map_to_remove.push_back(old_var);
        var_map_to_insert.push_back({new_var, new_base_gvars});

        old_gvars.remove_rev_varmap(m_rev_var_map);
        new_base_gvars.update_rev_varmap(m_rev_var_map, new_var);
      }
    }

    for (auto &v : var_map_to_remove) {
      m_var_map.erase(v);
    }
    for (auto &p : var_map_to_insert) {
      m_var_map.insert(p);
    }

    // Rename the rest
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
      if (crab_domain_params_man::get().region_is_dereferenceable()) {
        error_if_not_arity(3, 1);
        // ignore region variable (inputs[0])
        error_if_not_variable(inputs[1]);
        error_if_not_bool(outputs[0]);
        variable_t bv = outputs[0];
        variable_t ref = inputs[1].get_variable();
        CRAB_LOG("region-domain-is-deref",
                 crab::outs() << bv << ":= is_dereferenceable(" << inputs[0]
                              << "," << inputs[1] << "," << inputs[2] << ")\n"
                              << *this << "\n";);
        if (boost::optional<ghost_variables_t> ref_gvars = get_gvars(ref)) {

          if ((*ref_gvars).has_offset_and_size()) {
            // CRAB_LOG("region-domain-is-deref",
            // 	     crab::outs() << "\t" << "Found ghost variables for " << ref
            // << "\n"
            // 	     << "\toffset=" <<
            // (*ref_gvars).get_offset_and_size().get_offset() << "\n"
            // 	     << "\tsize=" <<
            // (*ref_gvars).get_offset_and_size().get_size() << "\n";);
            if ((*ref_gvars)
                    .get_offset_and_size()
                    .is_deref(m_base_dom,
                              rename_variable_or_constant(inputs[2]))) {
              set_bool_var_to_true(bv);
              CRAB_LOG("region-domain-is-deref", crab::outs()
                                                     << "\tRESULT=TRUE\n");
              return;
            }
          }
        }
        CRAB_LOG("region-domain-is-deref", crab::outs()
                                               << "\tRESULT=UNKNOWN\n");
        operator-=(bv);
      }
    } else {
      // pass the intrinsics to the base domain
      //=== base domain ===/
      std::vector<base_variable_or_constant_t> base_inputs;
      std::vector<base_variable_t> base_outputs;
      base_inputs.reserve(inputs.size());
      base_outputs.reserve(outputs.size());
      for (unsigned i = 0, sz = inputs.size(); i < sz; ++i) {
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        base_inputs.push_back(
            inputs[i].is_variable()
                ? base_variable_or_constant_t(
                      get_or_insert_gvars(inputs[i].get_variable()).get_var())
                : base_variable_or_constant_t(inputs[i].get_constant(),
                                              inputs[i].get_type()));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
        if (outputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        base_outputs.push_back(get_or_insert_gvars(outputs[i]).get_var());
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
        base_inputs.push_back(
            inputs[i].is_variable()
                ? base_variable_or_constant_t(
                      get_or_insert_gvars(inputs[i].get_variable()).get_var())
                : base_variable_or_constant_t(inputs[i].get_constant(),
                                              inputs[i].get_type()));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
        base_outputs.push_back(get_or_insert_gvars(outputs[i]).get_var());
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

    if (auto gvars_opt = get_gvars(ref)) {
      interval_t ival = m_base_dom[(*gvars_opt).get_var()];
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

  virtual bool
  get_allocation_sites(const variable_t &ref,
                       std::vector<allocation_site> &alloc_sites) override {

    if (!crab_domain_params_man::get().region_allocation_sites() || !ref.get_type().is_reference()) {
      return false;
    }

    allocation_sites out = m_alloc_site_dom[ref];

    if (out.is_top()) {
      return false;
    }

    if (!out.is_bottom()) {
      for (auto it = out.begin(), et = out.end(); it != et; ++it) {
        alloc_sites.push_back(*it);
      }
    }
    return true;
  }

  virtual bool get_tags(const variable_t &rgn, const variable_t &ref /*unused*/,
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
      for (auto it = tag_set.begin(), et = tag_set.end(); it != et; ++it) {
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
    return "RegionWithIsDerefDomain(" + m_base_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG(
          "region-print", o << "(RgnCounter=" << m_rgn_counting_dom << ","
                            << "RgnInit=" << m_rgn_init_dom;
          if (!crab_domain_params_man::get().region_skip_unknown_regions()) {
            o << ","
              << "RgnUnkType=" << m_rgn_type_dom;
          } if (crab_domain_params_man::get().region_allocation_sites()) {
            o << ","
              << "AllocSites=" << m_alloc_site_dom;
          } if (crab_domain_params_man::get().region_deallocation()) {
            o << ","
              << "RgnDealloc=" << m_rgn_dealloc_dom;
          } if (crab_domain_params_man::get().region_tag_analysis()) {
            o << ","
              << "Tags=" << m_tag_env;
          } o
          << ")\n";);

      CRAB_LOG(
          "region-print2", o << "MapVar={"; for (auto &kv
                                                 : m_var_map) {
            o << kv.first << "->" << kv.second << ";";
          } o << "},"
              << "BaseDom=" << m_base_dom << ")\n";);
      std::unordered_map<std::string, std::string> renaming_map;
      auto get_type_fn = [this](const variable_t &v) -> variable_type {
        if (has_dynamic_type(v)) {
          return get_dynamic_type_or_fail(v);
        } else {
          return variable_type(variable_type_kind::UNK_TYPE);
        }
      };
      // Assigns a string name to each ghost variable
      ghost_variables_t::mk_renaming_map(m_rev_var_map, get_type_fn,
                                         renaming_map);
      m_alloc.add_renaming_map(renaming_map);
      o << m_base_dom;
      m_alloc.clear_renaming_map();
    }
  }
}; // class region_with_is_deref_domain

template <typename Params>
struct abstract_domain_traits<region_with_is_deref_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab
