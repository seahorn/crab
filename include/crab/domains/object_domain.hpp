#pragma once

#include <boost/optional.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/symbolic_variable_eq_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/varname_factory.hpp>

#include <crab/domains/object/object_info.hpp>
#include <crab/domains/object/odi_map_domain.hpp>
#include <crab/domains/region/ghost_variable_manager.hpp>
#include <crab/domains/region/tags.hpp>

#include <set>
#include <unordered_map>

namespace crab {
namespace domains {
namespace object_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;

  static_assert(std::is_same<Number, typename BaseAbsDom::number_t>::value,
                "Number type and BaseAbsDom::number_t must be the same");
  // This is a strong requirement
  static_assert(
      std::is_same<varname_t, typename varname_allocator_t::varname_t>::value,
      "BaseAbsDom::varname_t and allocator_varname_t must be the same");
};

/// @brief a special log method to print vector
/// @tparam TType
/// @param o crab ostream
/// @param vec the vector for printing
template <typename TType>
void print_vector(crab::crab_os &o, const std::vector<TType> &vec) {
  typename std::vector<TType>::const_iterator it;
  o << "[";
  for (it = vec.begin(); it != vec.end(); it++) {
    if (it != vec.begin())
      o << ",";
    o << (*it);
  }
  o << "]";
}
} // end namespace object_domain_impl

#define OBJECT_DOMAIN_SCOPED_STATS(NAME) CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define OBJECT_COUNT_STATS(NAME) CRAB_DOMAIN_COUNT_STATS(NAME, 0)

/// @brief An abstract domain infers object invariants
/// @tparam Params
// This domain is based on the partitioning memory model,
// an object, more specifically an abstract object, is an object represents
// a number of concrete objects following the same dsa node.
// An object is spawn from a set of fields (or regions in crabIR).
// In addition, we model the memory architecture with a cache.
// The cache tracks the most recently used (MRU) object by memory load / store
// where the MRU object is one concrete object from the summarzied objects.
// The operations such as memory load or store are precisely abstracted
// if the requested properties can be found in the cache.
template <typename Params>
class object_domain final : public abstract_domain_api<object_domain<Params>> {
  using object_domain_t = object_domain<Params>;
  using abstract_domain_t = abstract_domain_api<object_domain_t>;

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
  /**------------------ Begin type definitions ------------------**/
  // type name for base domain
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_dom_varname_t = typename base_abstract_domain_t::varname_t;
  using base_dom_variable_vector_t =
      typename base_abstract_domain_t::variable_vector_t;
  using base_dom_variable_t = typename base_abstract_domain_t::variable_t;
  using base_dom_variable_or_constant_t =
      typename base_abstract_domain_t::variable_or_constant_t;
  using base_dom_linear_expression_t =
      typename base_abstract_domain_t::linear_expression_t;
  using base_dom_linear_constraint_t =
      typename base_abstract_domain_t::linear_constraint_t;
  using base_dom_linear_constraint_system_t =
      typename base_abstract_domain_t::linear_constraint_system_t;
  using base_dom_varname_allocator_t = typename Params::varname_allocator_t;

  // type name for cache & summary domains
  using field_abstract_domain_t = typename Params::field_abstract_domain_t;
  using flds_dom_varname_t = typename field_abstract_domain_t::varname_t;
  using flds_dom_variable_vector_t =
      typename field_abstract_domain_t::variable_vector_t;
  using flds_dom_variable_t = typename field_abstract_domain_t::variable_t;
  using flds_dom_variable_or_constant_t =
      typename field_abstract_domain_t::variable_or_constant_t;
  using flds_dom_linear_expression_t =
      typename field_abstract_domain_t::linear_expression_t;
  using flds_dom_linear_constraint_t =
      typename field_abstract_domain_t::linear_constraint_t;
  using flds_dom_linear_constraint_system_t =
      typename field_abstract_domain_t::linear_constraint_system_t;
  using flds_dom_varname_allocator_t = typename Params::varname_allocator_t;

  // type name for register domain
public:
  using eq_domain_value_t =
      symbolic_variable_equiality_domain<base_abstract_domain_t>;

private:
  using eq_register_domain_t = eq_domain_value_t;
  using usymb_t = typename eq_domain_value_t::domain_t;

  // type name for address domain
  using address_abstract_domain_t = eq_domain_value_t;
  using addr_value_domain_t = typename eq_domain_value_t::domain_t;

  using symb_flds_regs_map_t = std::unordered_map<
      usymb_t,
      std::pair<base_dom_variable_vector_t, base_dom_variable_vector_t>,
      symbolic_variable_equiality_domain_impl::symb_var_hash>;

  // Management of ghost variables
  // We model two types of variables in object domain
  // One is for region variable, another is for reference variable.
  // The former one is used to perform the reduction between region and register
  // , the later one is used to infer sea.is_deref.
  // Note that, we are using fixed naming where
  // variable_t and subdomain::variable_t have the same type
  // WARN: if the type of the domain for object is different from
  // the type for base domain, the ghost variable manager should split into two.
  using ghost_var_num_man_t =
      typename region_domain_impl::ghost_variable_manager_with_fixed_naming<
          object_domain_t, base_abstract_domain_t>;
  using ghost_var_eq_man_t =
      typename region_domain_impl::ghost_variable_manager_with_fixed_naming<
          object_domain_t, eq_domain_value_t>;

public:
  using ghost_variables_eq_t = typename ghost_var_eq_man_t::ghost_variables_t;

private:
  using ghost_variables_t = typename ghost_var_num_man_t::ghost_variables_t;
  using ghost_variable_vector_t = typename std::vector<ghost_variables_t>;
  using ghost_variable_kind = typename ghost_variables_t::ghost_variable_kind;

  // Object identifier / object id
  // ODI map: object id -> ODI product domain
  // Map from an object id to object's product domain:
  // Object infos: <Count, Used, Dirty>
  // Count: count how many allocations are owned by an abstract object.
  // Each call to ref_make models a new allocation.
  // This allows us to decide when only if one address per object
  // (i.e., singleton object).
  // Used: indicate whether the cache domain is used.
  // Dirty: indicate whether the cache domain is dirty.
  //
  // Object vals: <summary dom, cache dom, symb dom>
  // Summary domain: An domain keeps properties for all summarized objects
  // Cache domain: An domain for most recetnly used object
  // symb domain:
  //    keeps what uninterpreted symbols assign to fields
  //    NOTE that: we could keep a map from each field to
  //               an uninterpreted symbol
  // An object id indicates which kinds of object is.
  using obj_id_t = variable_t;
  using odi_map_t =
      object_domain_impl::odi_map_domain<obj_id_t, object_domain_t,
                                         field_abstract_domain_t>;
  using object_info_t = typename odi_map_t::odi_info_t;
  using object_value_t = typename odi_map_t::odi_value_t;
  using odi_domain_product_t = typename odi_map_t::map_raw_value_t;
  using summary_domain_t = typename odi_map_t::summary_domain_t;
  using cache_domain_t = typename odi_map_t::cache_domain_t;
  using eq_fields_domain_t = typename odi_map_t::eq_domain_t;

  // Object fields to id map
  // Map region variables to object's id
  // This map is constructed during intrinsic call and share in common
  using obj_flds_id_map_t =
      std::shared_ptr<std::unordered_map<variable_t, obj_id_t>>;

  // References to corresponding base addresse variable as ghost
  // Map reference variables to variables used in equality domain
  // This map is constructued for address domain, we keep track
  // the variables created for references.
  // The map is constructed dynamically but share accress all states
  using refs_base_addrs_map_t =
      std::shared_ptr<std::unordered_map<variable_t, variable_t>>;
  /**------------------ End type definitions ------------------**/
  // FIXME: the current solution does not support different dom operations
  static_assert(
      std::is_same<base_abstract_domain_t, field_abstract_domain_t>::value,
      "base_abstract_domain_t and field_abstract_domain_t must be the same");

  // static_assert(
  //     std::is_same<eq_register_domain_t, eq_fields_domain_t>::value,
  //     "eq_register_domain_t and eq_fields_domain_t must be the same");

  /**------------------ Begin class field definitions ------------------**/

  // a special symbol to represent bottom state of the object domain
  // object domain is bottom if this flag is set or all subdomains are bottom.
  bool m_is_bottom;

  // The abstract state definition:
  // Base domain:
  // Domain for register, references, booleans
  base_abstract_domain_t m_base_dom;

  // A map from each region variable to corresponding object domain
  odi_map_t m_odi_map;

  // A domain to keep equalities over base addresses of objects.
  address_abstract_domain_t m_addrs_dom;

  // A domain to keep what uninterpreted symbols assign to registers if those
  // are equal to some objects' fields
  eq_register_domain_t m_eq_regs_dom;

  // Map region variables to an object id for an abstract object
  // To determine the object id,
  // we choose the first region variable passed by the intrinsic method.
  // If any regions that not indicating by the intrinsic, we add those regions
  // into this method during the evaluation of make_ref.
  // Invariants: the field -> id map event if we partially rename it.
  // WARN: for now, we did not cover unknown regions.
  obj_flds_id_map_t m_flds_id_map;

  // Map reference variables to base addresses variables
  // The base address is used to represent memory allocation address for a
  // concrete memory object. In addition, each cache that used for an abstract
  // object has a special base address variable for the MRU object
  refs_base_addrs_map_t m_refs_base_addrs_map;

  // A ghost variable manager to create / get ghost variable for
  // each variable
  // This is required if subdomain is any array domain
  // (array_smashing or array_adaptive) or
  // logical-numerical (flat_boolean_domain) or
  // numerical (zones, octagons, pk, etc).
  // They should use only array, integer or boolean variables.
  // So in object domain, we takes typed variable (e.g. region, reference)
  // but in subdomain like base domain, odi map,
  // we only keep array, integer or boolean variables
  // Details: https://github.com/seahorn/crab/wiki/IndexedNamesAndTypedVariables
  ghost_var_num_man_t m_ghost_var_num_man;
  ghost_var_eq_man_t m_ghost_var_eq_man;

/* clang-format off */
  // Domain hierarchy:
  //                      m_base_dom:                       m_addrs_dom:
  //                     ┌─────────────────────────────┐
  // Crab IR variables   │    numerical properties     │
  //                     │       for Regs, Refs        │    ┌─────────────────┐
  //  Regs, Rgns, Refs   │─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─│    │  Ref|-> #symb   │
  //                     │     equalities for Regs     │    └─────────────────┘
  //                     │                             │             ▲
  //                     │Reg|-> #symb ◀ ┬ ─ ─ ─ ─ ─ ─ ┼ ─ ─ ─ ─ ─ ─ | ─ ─ ─ ┐
  //                     └───────────────┼─────────────┘             │       |
  //             ┌─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─│─ ─ ─ ─ ─ ─ ─┬─ ─ ─ ─ ─ ─ ─┘       │
  //  m_odi_map: │                       |             │                     |
  // ┌───────────┼────┬────────────────┐ │┌────────────┼───┬───────────────┐ │
  // │           ▼    │ equalities for │ |│            ▼   │ equalities for│ |
  // │mru_base_addres |      Rgns      │ ││mru_base_addres |      Rgns     │ │
  // │   |-> #symb    │ Rgn|-> #symb  ◀┼-┘│   |-> #symb    │  Rgn|-> #symb ◀─┘
  // │                |         ▲  ▲   │  │                |         ▲  ▲  │
  // ├ ─ ─ ─ ─ ─ ─ ─ ─┴─ ─ ─ ─ ─│─ ┼  ─│  ├ ─ ─ ─ ─ ─ ─ ─ ─┴─ ─ ─ ─ ─│─ ┼  ┤ ...
  // │cache:                    ▼  ▼   │  │cache:                    ▼  ▼  │
  // │  numerical properties for Rgns  │  │  numerical properties for Rgns │
  // ├ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  ─│  ├ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  ┤
  // │summary:                         │  │summary:                        │
  // │  numerical properties for Rgns  │  │  numerical properties for Rgns |
  // └─────────────────────────────────┘  └────────────────────────────────┘
/* clang-format on */
  // In addition, we specify alternative for improvement of precision
  // Enabling a flag object.singletons_in_base=true
  // The domain will keep all singleton objects in base domain,
  // singleton object:
  //    the object info implies the number of reference is exactly one
  // non-singleton object:
  //    the object info shows the number of reference is > 1
  // I.e. for each abstract object,
  //  there is only three states of object reference: 0, 1, and [1, +oo]
  /**------------------ End class field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  // Constructor Definition
  object_domain(base_abstract_domain_t &&base_dom, odi_map_t &&odi_map,
                address_abstract_domain_t &&addrs_dom,
                eq_register_domain_t &&eq_regs_dom,
                obj_flds_id_map_t &&flds_id_map,
                refs_base_addrs_map_t &&refs_base_addrs_map,
                ghost_var_num_man_t &&ghost_var_num_man,
                ghost_var_eq_man_t &&ghost_var_eq_man)
      : m_is_bottom(base_dom.is_bottom()), m_base_dom(base_dom),
        m_odi_map(odi_map), m_addrs_dom(addrs_dom), m_eq_regs_dom(eq_regs_dom),
        m_flds_id_map(flds_id_map), m_refs_base_addrs_map(refs_base_addrs_map),
        m_ghost_var_num_man(ghost_var_num_man),
        m_ghost_var_eq_man(ghost_var_eq_man) {}

  std::function<variable_type(const variable_t &)> get_type() const {
    auto fn = [](const variable_t &v) { return v.get_type(); };
    return fn;
  }

  ghost_variables_t get_or_insert_gvars(const variable_t &v) {
    m_ghost_var_eq_man.get_or_insert(v);
    return m_ghost_var_num_man.get_or_insert(v);
  }

  ghost_variables_eq_t get_or_insert_eq_gvars(const variable_t &v) {
    return m_ghost_var_eq_man.get_or_insert(v);
  }

  boost::optional<ghost_variables_t> get_num_gvars(const variable_t &v) const {
    return m_ghost_var_num_man.get(v);
  }

  boost::optional<ghost_variables_eq_t>
  get_eq_gvars(const variable_t &v) const {
    return m_ghost_var_eq_man.get(v);
  }

  base_dom_variable_or_constant_t
  rename_variable_or_constant(const variable_or_constant_t &v) {
    if (v.is_constant()) {
      return base_dom_variable_or_constant_t(v.get_constant(), v.get_type());
    } else {
      auto v_gvars = get_or_insert_gvars(v.get_variable());
      return base_dom_variable_or_constant_t(v_gvars.get_var());
    }
  }

  boost::optional<base_dom_variable_t>
  rename_variable_optional(boost::optional<variable_t> v_opt) {
    if (v_opt == boost::none) {
      return boost::none;
    } else {
      if (is_unknown_region(*v_opt)) {
        return boost::none;
      }
      ghost_variables_t v_gvars = get_or_insert_gvars(*v_opt);
      return v_gvars.get_var();
    }
  }

  bool is_unknown_region(const variable_t &v) const {
    return v.get_type().is_unknown_region();
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

  /// @brief performs the join operation between current abstract state and the state \p right passed by
  /// @param[in] right an abstract state serves as an input for the join operation
  /// @note Perform *this = join(*this, right)
  /// @pre \p this, \p right are two different abstract states
  /// @post the new \p this satisfies: \gamma(the new \p this) \leq^{concrete} \gamma(join( \p this, \p right ))
  void self_join(const object_domain_t &right) {

    // The join operation performs as follows:
    // 1. if we do not keep singleton objects in the base domain,
    //    the join is a pairwise join on each value stored in the odi map
    // 2. Otherwise, for each abstract object, generally,
    //    we have four cases at the join:

    // Case 1: not singleton, odi map saved properties respectively
    //
    // Join: we join the base dom and odi map
    // E.g.
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 = 3, V_2 <= 5 }
    //  s1 U s2 : Base in s1 U Base in s2,
    //            Object (V_1, V_2) = s1.odi_map U s2.odi_map
    //                              = { V_1 <= V_2, V_2 <= 5 }

    // Case 2:
    // *this (right): singleton object in the base domain,
    //  odi map does not contain properties for that object
    // right (*this): not singleton, odi map stores invariants for that object
    //
    // Join: project singleton object, join with non-singleton in the odi map
    // E.g.
    //  Let s1 be the first state, s2 is second one.
    //  Let Vs be the fields in that abstract object.
    //  For Base in s1,
    //  only_singleton = Base.project(Vs)
    //  no_singleton = Base.forget(Vs)
    //  s1 U s2 :
    //  no_singleton U Base in s2,
    //      Object (Vs) in s2 U only_singleton
    // E.g.
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 <= 3, V_2 <= 5 }
    //  s1 U s2:
    //    Base = s1.Base.forget(V_1, V_2) U s2.Base
    //    Object (V_1, V_2) = { V_1 <= 3, V_2 <= 6 }

    // Other cases (simple):
    // one state: exists such an abstract object
    // another: does not exist that abstract object
    //
    // Join: we join the base dom and odi map
    // Example1:
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 U s2 :
    //    Base = { ... }, Object (V_1, V_2) = empty

    // Example2:
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 U s2 :
    //    Base = Base^s1 U Base^s2,
    //    Object (V_1, V_2) = empty

    // Overall, for both cases, the only second case requires move singleton
    // into odi map. Other cases require pairwise join on each subdomain

    // The DSA node info is passing via a specific intrinsic at the begining statement
    // from each function definition in the crab IR.
    // If there is any missing DSA information, we update it by merging two maps
    // Ideally, DSA information should be static, but here, we gradually obtain
    // information.

    if (!m_flds_id_map) {
      m_flds_id_map = right.m_flds_id_map;
    }
    if (m_flds_id_map && right.m_flds_id_map &&
        m_flds_id_map != right.m_flds_id_map) { // dsa info are inconsistent
      m_flds_id_map->insert(right.m_flds_id_map->begin(),
                            right.m_flds_id_map->end());
    }

    // m_refs_base_addrs_map do not require join,
    // they share accross all states. No need to update
    if (right.m_refs_base_addrs_map && !m_refs_base_addrs_map) {
      m_refs_base_addrs_map = right.m_refs_base_addrs_map;
    }

    /********* pairwise join *********/
    // 1. join odi maps
    // The reason to copy base dom is we need to update base domain
    // if cache is committed. We perform a domain reduction between
    // base domain and cache subdomain in the odi map before cache invalidation.
    // TODO: avoid copying this at first, check cache invalidation, 
    // return new computed states if domain reduction is performed
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    m_odi_map.compound_join(right.m_odi_map, *this, right, m_base_dom,
                            r_base_dom);

    // 2. join base domains
    m_base_dom |= r_base_dom;

    // 3. join address domains
    m_addrs_dom |= right.m_addrs_dom;

    // 4. join regs' equality domains
    m_eq_regs_dom |= right.m_eq_regs_dom;
  }

  /// @brief performs the join / widen operation between two abstract state passed by references
  /// @param[in] left an abstract state serves as an input for the join operation
  /// @param[in] right an abstract state serves as an input for the join operation
  /// @param[in] is_join true for join operation; false for widening
  /// @return the new abstract state computed by \p left and \p right
  /// @pre \p left, \p right are two different abstract states
  /// @post \p \return satisfies: \gamma( \p \return ) \leq^{concrete} \gamma(join( \p left, \p right ))
  object_domain_t join_or_widening(const object_domain_t &left,
                                   const object_domain_t &right,
                                   const bool is_join) const {

    obj_flds_id_map_t out_flds_id_map =
        left.m_flds_id_map ? left.m_flds_id_map : right.m_flds_id_map;
    if (left.m_flds_id_map && right.m_flds_id_map &&
        m_flds_id_map != right.m_flds_id_map) {
      out_flds_id_map->insert(right.m_flds_id_map->begin(),
                              right.m_flds_id_map->end());
    }

    refs_base_addrs_map_t out_refs_base_addrs_map =
        left.m_refs_base_addrs_map ? left.m_refs_base_addrs_map
                                   : right.m_refs_base_addrs_map;

    // 1. join odi maps
    // Copying base domains is required for flushing caches
    base_abstract_domain_t l_base_dom = left.m_base_dom;
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    odi_map_t out_odi_map =
        is_join ? left.m_odi_map.join(right.m_odi_map, left, right, l_base_dom,
                                      r_base_dom)
                : left.m_odi_map.widening(right.m_odi_map, left, right,
                                          l_base_dom, r_base_dom);

    // 2. join base domains
    base_abstract_domain_t out_base_dom(is_join ? l_base_dom | r_base_dom
                                                : l_base_dom || r_base_dom);

    // 3. join address domains
    address_abstract_domain_t out_addrs_dom(
        is_join ? left.m_addrs_dom | right.m_addrs_dom
                : left.m_addrs_dom || right.m_addrs_dom);

    // 4. join regs' equality domains
    eq_register_domain_t out_eq_regs_dom(
        is_join ? left.m_eq_regs_dom | right.m_eq_regs_dom
                : left.m_eq_regs_dom || right.m_eq_regs_dom);

    ghost_var_num_man_t out_ghost_var_num_man(left.m_ghost_var_num_man);
    ghost_var_eq_man_t out_ghost_var_eq_man(left.m_ghost_var_eq_man);

    object_domain_t res(
        std::move(out_base_dom), std::move(out_odi_map),
        std::move(out_addrs_dom), std::move(out_eq_regs_dom),
        std::move(out_flds_id_map), std::move(out_refs_base_addrs_map),
        std::move(out_ghost_var_num_man), std::move(out_ghost_var_eq_man));
    return res;
  }

  /// @brief meet / narrowing two abstract states, return a new state
  /// @param[in] left an abstract state serves as an input for the meet operation
  /// @param[in] right an abstract state serves as an input for the meet operation
  /// @param[in] is_meet true for meet operation; false for narrowing
  /// @return the new abstract state computed by \p left and \p right
  /// @pre \p left, \p right are two different abstract states
  /// @post \p \return satisfies: 
  /// \meet^{concrete}(\gamma( \p left ), \gamma( \p right )) \leq^{concrete} \gamma( \p \return )
  object_domain_t meet_or_narrowing(const object_domain_t &left,
                                    const object_domain_t &right,
                                    const bool is_meet) const {

    // 1. if we do not keep singleton objects in the base domain,
    //    the meet is a pairwise meet on each value stored in the odi map
    // 2. Otherwise, for each abstract object, generally,
    //    we have four cases at the meet:

    // Case 1: not singleton, odi map saved properties respectively
    //
    // meet: we meet the base dom and odi map
    // E.g.
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 = 3, V_2 <= 5 }
    //  s1 meet s2 : s1.Base meet s2.Base,
    //            Object1 (V_1, V_2) = s1.odi_map meet s2.odi_map
    //                               = { V_1 = 3, V_2 <= 5 }

    // Case 2:
    // *this (right): singleton object in the base domain,
    //  odi map does not contain properties for that object
    // right (*this): not singleton, odi map stores invariants for that object
    //
    // meet: meet non-singleton in odi map with singleton in base domain
    //  Let s1 be the first state, s2 is second one.
    //  Let Vs be the fields in that abstract object.
    //  For Base in s1,
    //  only_singleton = object in s2
    //  object in s2 = top (empty)
    //  s1 meet s2 :
    //  only_singleton meet Base in s1,
    //      Object (Vs) in s2 = top
    // E.g.
    //  s1: Base = { V_1 = 3, V_2 = 4 }, Object (V_1, V_2) = empty
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s1 meet s2:
    //    b = s1.Base meet { V_1 <= V_2, V_2 <= 5 }
    //    Base = b meet s2.Base
    //    Object (V_1, V_2) = empty

    // Other cases:
    // one state: exists such an abstract object
    // another: does not exist that abstract object
    // meet
    // Example1:
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 meet s2 :
    //    Base = s1.Base meet s2.Base, Object (V_1, V_2) = empty

    // Example2:
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 meet s2 :
    //   Base = s1.Base meet s2.Base,
    //   Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }

    obj_flds_id_map_t out_flds_id_map =
        left.m_flds_id_map ? left.m_flds_id_map : right.m_flds_id_map;
    if (left.m_flds_id_map && right.m_flds_id_map &&
        m_flds_id_map != right.m_flds_id_map) {
      out_flds_id_map->insert(right.m_flds_id_map->begin(),
                              right.m_flds_id_map->end());
    }

    refs_base_addrs_map_t out_refs_base_addrs_map =
        left.m_refs_base_addrs_map ? left.m_refs_base_addrs_map
                                   : right.m_refs_base_addrs_map;

    // 1. meet odi maps
    // Copying base domains is required for flushing caches
    base_abstract_domain_t l_base_dom = left.m_base_dom;
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    odi_map_t out_odi_map =
        is_meet ? left.m_odi_map.meet(right.m_odi_map, left, right, l_base_dom,
                                      r_base_dom)
                : left.m_odi_map.narrowing(right.m_odi_map, left, right,
                                           l_base_dom, r_base_dom);

    // 2. meet base domains
    base_abstract_domain_t out_base_dom(is_meet ? l_base_dom & r_base_dom
                                                : l_base_dom && r_base_dom);

    // 3. meet address domains
    address_abstract_domain_t out_addrs_dom(
        is_meet ? left.m_addrs_dom & right.m_addrs_dom
                : left.m_addrs_dom && right.m_addrs_dom);

    // 4. meet regs' equality domains
    eq_register_domain_t out_eq_regs_dom(
        is_meet ? left.m_eq_regs_dom & right.m_eq_regs_dom
                : left.m_eq_regs_dom && right.m_eq_regs_dom);

    ghost_var_num_man_t out_ghost_var_num_man(left.m_ghost_var_num_man);
    ghost_var_eq_man_t out_ghost_var_eq_man(left.m_ghost_var_eq_man);

    object_domain_t res(
        std::move(out_base_dom), std::move(out_odi_map),
        std::move(out_addrs_dom), std::move(out_eq_regs_dom),
        std::move(out_flds_id_map), std::move(out_refs_base_addrs_map),
        std::move(out_ghost_var_num_man), std::move(out_ghost_var_eq_man));
    return res;
  }

  // compare two abstract states, return boolean
  bool less_than_eq(const object_domain_t &left,
                    const object_domain_t &right) const {

    bool res = left.m_base_dom <= right.m_base_dom;
    CRAB_LOG("object-leq", crab::outs() << "Result3=" << res << "\n";);
    if (!res) {
      return false;
    }

    res &= (left.m_odi_map <= right.m_odi_map);
    CRAB_LOG("object-leq", crab::outs() << "Result4=" << res << "\n";);
    if (!res) {
      return false;
    }

    res &= (left.m_addrs_dom <= right.m_addrs_dom);
    CRAB_LOG("object-leq", crab::outs() << "Result5=" << res << "\n";);
    if (!res) {
      return false;
    }

    res &= (left.m_eq_regs_dom <= right.m_eq_regs_dom);
    CRAB_LOG("object-leq", crab::outs() << "Result6=" << res << "\n";);
    return res;
  }

  /***************** Fields and Object id operations *****************/
  /// @brief get the corresponding object id
  /// @param[in] rgn a region variable
  /// @return if we find object id, returns it; otherwise, return boost::none
  /// @note object id is a special variable representing an abstract object
  boost::optional<obj_id_t> get_obj_id(const variable_t &rgn) const {
    if (!m_flds_id_map) {
      return boost::none;
    }
    auto it = (*m_flds_id_map).find(rgn);
    if (it == (*m_flds_id_map).end()) {
      return boost::none;
    }
    return it->second;
  }

  /// @brief get the corresponding object id or report error
  /// @param[in] rgn a region variable
  /// @return if we find object id, returns; otherwise, abort execution with an error
  obj_id_t get_obj_id_or_fail(const variable_t &rgn) const {
    if (!m_flds_id_map) {
      CRAB_ERROR(domain_name(), "::", __func__, ", the odi map does not exist");
    }
    auto it = (*m_flds_id_map).find(rgn);
    if (it == (*m_flds_id_map).end()) {
      dump();
      CRAB_ERROR(domain_name(), "::", __func__,
                 ", the odi map does not include the region ", rgn,
                 "(type: ", rgn.get_type(),
                 ") belonging to an abstract object");
    }
    return it->second;
  }

  /// @brief create a new object id to represent an abstract object
  /// @param[in] v a field from a new object
  /// @return either creating a shadow variable to represent id or use passed
  /// field.
  /// @note this function will only be used when interpreter analyzes crab dsa intrinsic
  obj_id_t create_new_obj_id(const variable_t &v) {
#ifdef CREATE_NEW_ID_NAME
    varname_t v_var_name = v.name();
    std::string str_id_name =
        "id_" +
        std::to_string(object_domain_impl::id_val_generator_t::get_next_val());
    varname_t id_name = v_var_name.get_var_factory().get(str_id_name);
    obj_id_t v_id(id_name, crab::REG_UNKNOWN_TYPE, 32);
    return v_id;
#else
    return v;
#endif
  }

  /// @brief get all object fields by giving an id
  /// @param[in] id an object id
  /// @param[in, out] obj_flds the vector stored the result
  /// @note this function should be distinguished with \c get_obj_flds_with_ghost
  void get_obj_flds(const obj_id_t &id, variable_vector_t &obj_flds) const {
    assert(m_flds_id_map);
    for (auto kv : (*m_flds_id_map)) {
      if (kv.second == id) {
        obj_flds.push_back(kv.first);
      }
    }
  }

  /// @brief get all object fields with ghost variables by giving an id
  /// @param[in] id an object id
  /// @param[in, out] obj_ghost_flds the vector stored the result
  /// @note this function has different specs from \c get_obj_flds
  /// The ghost variables here could be the base address, offset, and size
  /// for a region with reference types
  /// e.g. In crab IR, if a region V_3:region(ref), the ghost variables are:
  ///    V_3.address, V_3.offset, V_3.size
  /// \c get_obj_flds_with_ghost will fetch all ghost variables used by object fields
  void get_obj_flds_with_ghost(const obj_id_t &id,
                               base_dom_variable_vector_t &obj_dom_flds) const {

    auto get_obj_ghost_flds = [&](const obj_id_t &id,
                                  ghost_variable_vector_t &obj_ghost_flds) {
      variable_vector_t obj_flds;
      get_obj_flds(id, obj_flds);
      for (auto v : obj_flds) {
        if (is_unknown_region(v)) {
          continue;
        }
        auto gvars_opt = get_num_gvars(v);
        obj_ghost_flds.push_back(*gvars_opt);
      }
    };
    ghost_variable_vector_t obj_ghost_flds;
    get_obj_ghost_flds(id, obj_ghost_flds);
    for (auto v : obj_ghost_flds) {
      if (v.has_offset_and_size()) {
        obj_dom_flds.push_back(v.get_offset_and_size().get_offset());
        obj_dom_flds.push_back(v.get_offset_and_size().get_size());
      }
      obj_dom_flds.push_back(v.get_var());
    }
  }

  /// @brief based on each dsa intrinsic in the crab IR, group region variables as
  /// fields for each abstract object
  /// @note abstract object \equiv dsa node
  ///       region variables used in one dsa node \equiv object fields
  /// @param rgn the region belongs to one abstract object
  /// @param id an object id
  void update_fields_id_map(const variable_t &rgn, const obj_id_t &id) {
    // update m_flds_id_map
    if (!m_flds_id_map) { // create a object fields - id map if it is null
      m_flds_id_map =
          std::make_shared<std::unordered_map<variable_t, obj_id_t>>();
    }
    auto it = (*m_flds_id_map).find(rgn);
    if (it != (*m_flds_id_map).end()) {
      it->second = id;
    } else {
      (*m_flds_id_map).insert({rgn, id});
    }
  }

  /***************** Base address operations *****************/
  /// @brief get or create a variable representing the base address by given a
  /// reference variable
  /// @param[in] v the reference variable or an object id
  /// @return a ghost variable represent its base address
  /// specifically, for abstract object, the base address is used to indicate
  /// the mru object stored on cache subdomain. For mru object, we use a
  /// variable named `<obj_id>_mru_base`.
  /// @pre the input variable \p v must be either a reference or an object id
  /// @note use this function when the base address is required for operation
  /// @warning do not use this function to check whether a based address is
  /// created
  variable_t get_or_insert_base_addr(const variable_t &v) {
    assert(v.get_type().is_reference() || v.get_type().is_region());
    if (!m_refs_base_addrs_map) { // create a ref - base addr map if it is null
      m_refs_base_addrs_map =
          std::make_shared<std::unordered_map<variable_t, variable_t>>();
    }
    auto it = (*m_refs_base_addrs_map).find(v);
    if (it != (*m_refs_base_addrs_map).end()) {
      return it->second;
    }
    varname_t v_var_name = v.name();
    std::string str_base_name =
        v.get_type().is_reference() ? "_base" : "_mru_base";
    varname_t base_name = v_var_name.get_var_factory().get_or_insert_varname(
        v_var_name, str_base_name);
    variable_t v_base(base_name, crab::INT_TYPE, 32);
    (*m_refs_base_addrs_map).insert({v, v_base});
    return v_base;
  }

  /// @brief get the corresponding base address variable
  /// @param[in] v a reference variable or an object id
  /// @return if we find the ghost, returns; otherwise, return boost::none
  /// @note use this function for special handling when we do not find the ghost
  boost::optional<variable_t> get_base_addr(const variable_t &v) const {
    assert(v.get_type().is_reference() || v.get_type().is_region());
    if (!m_refs_base_addrs_map) {
      return boost::none;
    }
    auto it = (*m_refs_base_addrs_map).find(v);
    if (it == (*m_refs_base_addrs_map).end()) {
      return boost::none;
    }
    return it->second;
  }

  // /// @brief get the corresponding base address variable or report error
  // /// @param rgn a reference variable or an object id
  // /// @return if we find the ghost, returns; otherwise, abort with an error
  // variable_t get_base_addr_or_fail(const variable_t &v) const {
  //   assert(v.get_type().is_reference() || v.get_type().is_region());
  //   if (!m_refs_base_addrs_map) {
  //     CRAB_ERROR(domain_name(), "::", __func__,
  //                ": base address map is not created");
  //   }
  //   auto it = (*m_refs_base_addrs_map).find(v);
  //   if (it == (*m_refs_base_addrs_map).end()) {
  //     CRAB_ERROR(domain_name(), "::", __func__, ": ", v, " not found");
  //   }
  //   return it->second;
  // }

  /// @brief check whether two reference variable refer the same memory object
  /// @param[in] x a reference variable
  /// @param[in] y a reference variable
  /// @return true if x_base == y_base satisfies in address domain
  bool test_two_addrs_equality(const variable_t &x, const variable_t &y) {
    return m_addrs_dom.equals(get_or_insert_base_addr(x),
                              get_or_insert_base_addr(y));
  }

  /// @brief check the reference refer the MRU object
  /// @param ref a reference variable
  /// @param id an object id
  /// @return true if ref_base == id_mru_base satisfies in address domain
  bool test_ref_refer_mru_object(const variable_t &ref, const obj_id_t &id) {
    return test_two_addrs_equality(ref, id);
  }

  /***************** ODI map operations *****************/
  /// @brief change singleton object to non-singleton
  /// @pre based on number of references in object_info, the abstract object now becomes a non-singleton object
  /// @post the abstract object referred by \p id is non-singleton
  /// @param[in] id an object id
  /// @par Side Effects:
  ///   This function modifies the abstract state which invokes this function.
  ///   Depending on whether we put properties for singleton object on the base domain,
  ///   changing singleton obejct involves different domain operations 
  void change_object_status(const obj_id_t &id) {

    base_dom_variable_vector_t flds_vec;
    get_obj_flds_with_ghost(id, flds_vec);

    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    if (!prod_ref) {
      CRAB_ERROR(domain_name(), "::", __func__, ": object ", id,
                 " value is not found on the odi map");
    }
    // NOTE: copy the following value is required because
    // updating the odi map is construcing a new tree based on previous tree
    odi_domain_product_t res_prod = *prod_ref;
    object_info_t &obj_info = odi_map_t::object_info_val(res_prod);
    // increase reference counting
    obj_info.refcount_val().increment(id);

    // Two cases here:
    //  a. if we keep singleton in base, we need to copy its properties
    //  from the base to the odi map;
    //     m_odi_map[id].sum = m_base_dom.project(flds)
    //  b. otherwise, although there is no need to update cache,
    //  we still need to copy from the cache to the summary.
    // The copy is required because the summary subdomain represents summary
    // objects which is the abstraction of one or more concrete objects.
    object_value_t &odi_val = odi_map_t::object_odi_val(res_prod);
    if (!crab_domain_params_man::get().singletons_in_base()) {
      // Not necessary to reset cache because there is no ref_store/ref_load
      // occur. However, we still need to commit cache (reduction might need).
      apply_reduction_between_object_and_base(m_base_dom, odi_val, id);
      m_is_bottom = m_base_dom.is_bottom();
      // copy singleton object into summary
      odi_map_t::object_sum_val(odi_val) = odi_map_t::object_cache_val(odi_val);
      odi_map_t::object_eq_val(odi_val).set_to_top();
    } else {
      // TODO: performance may matter if the base dom contains huge dimensions
      // instead of copy the base domain and perform project, it is better to
      // have a special project onto operation returns a new domain after
      // projection.
      base_abstract_domain_t singleton_base(m_base_dom);
      m_ghost_var_num_man.project(flds_vec, singleton_base);

      // copy singleton object into summary
      odi_map_t::object_sum_raw_val(odi_val) = std::move(singleton_base);
      for (auto &fld : flds_vec) {
        m_ghost_var_num_man.forget(fld, m_base_dom);
      }
    }
    // update odi map
    m_odi_map.set(id, std::move(res_prod));

    CRAB_LOG(
        "object-change-status",
        crab::outs() << "singleton to non-sigleton on object: " << id << "\n";
        m_odi_map.odi_write(crab::outs(), res_prod); crab::outs() << "\n";);
  }

  /***************** Cache operations *****************/
  /// @brief invalidate cache for the abstract object \p id if it is missed
  /// @details cache missed when a \p ref refers a different memory object from
  /// the cache since \p ref is going to perform read / write
  /// @param id an object id
  /// @param ref the reference variable
  /// @param rgn the region which reference referred
  /// @param reg_symb optional symbolic variable representing equality between
  ///   register and field
  /// @param offset_size_symb optional symbolic variable for ghost variable
  void invalidate_cache_if_miss(
      obj_id_t &id, const variable_t &ref, const variable_t &rgn,
      boost::optional<usymb_t> &reg_symb,
      boost::optional<std::pair<usymb_t, usymb_t>> &offset_size_symb) {

    // get the variable represented the mru base address
    variable_t mru_obj_base = get_or_insert_base_addr(id);
    CRAB_LOG(
        "object-entailment",
        crab::outs() << mru_obj_base << " == " << get_or_insert_base_addr(ref)
                     << "?\n"
                     << "Addrs = " << m_addrs_dom << "\n"
                     << "Is mru cached? "
                     << (test_two_addrs_equality(id, ref) ? "true" : "false")
                     << "\n";);

    bool update_new_mru = m_odi_map.invalidate_cache_if_miss(
        id, *this, get_or_insert_eq_gvars(rgn), reg_symb, offset_size_symb,
        test_two_addrs_equality(id, ref));

    // update address dom and object info if cache is missed
    if (update_new_mru) {
      m_addrs_dom.add(get_or_insert_base_addr(ref), mru_obj_base);
    }
  }

  /// @brief project \p base_dom onto the dimensionality of \p m_eq_regs_dom
  /// @param[in] base_dom the base domain that used for projection
  /// @return a new domain projected based on the scope of \p m_eq_regs_dom
  base_abstract_domain_t
  project_onto_eq_regs(const base_abstract_domain_t &base_dom) const {
    auto regs_opt = m_eq_regs_dom.get_all_variables();
    if (regs_opt) {
      // To avoid copy, we first perform projection.
      base_abstract_domain_t regs_only_base = base_dom;
      regs_only_base.project(*regs_opt);
      return regs_only_base;
    } else {
      return base_dom.make_top();
    }
  }

  /// @brief project \p cache_dom onto the dimensionality of \p eq_flds_dom
  /// @param[in] cache_dom the numerical domain stores properties of fields
  /// @param[in] eq_flds_dom the equality domain for fields
  /// @return  a new domain projected based on the scope of \p eq_flds_dom
  cache_domain_t
  project_onto_eq_flds(const cache_domain_t &cache_dom,
                       const eq_fields_domain_t &eq_flds_dom) const {
    auto flds_opt = (*eq_flds_dom).get_all_variables();
    if (flds_opt) {
      // To avoid copy, we first perform projection.
      cache_domain_t flds_only_cache = cache_dom;
      flds_only_cache.project(*flds_opt);
      return flds_only_cache;
    } else {
      return cache_dom.make_top();
    }
  }

  /// @brief Applying domain reduction based on compilation flags
  /// @details a method to apply domain reduction before each transfer function by
  /// specifying an domain parameter:
  ///  object.reduction_level=FULL_REDUCTION or
  ///  object.reduction_level=REDUCTION_BEFORE_CHECK && reduce_for_assert = true
  /// @note domain reduction could be disabled if parameters:
  ///  object.reduction_level=NO_REDUCTION
  /// @param[in] reduce_for_assert a special bool indicate reduction is called
  /// before assertion
  void apply_reduction_based_on_flags(bool reduce_for_assert = false) {
    OBJECT_DOMAIN_SCOPED_STATS(".overall_reduction");
    if (crab_domain_params_man::get().reduction_level() ==
            object_domain_params::reduction_level_t::NO_REDUCTION ||
        (crab_domain_params_man::get().reduction_level() ==
             object_domain_params::reduction_level_t::REDUCTION_BEFORE_CHECK &&
         !reduce_for_assert)) {
      // if the flags are set properly, no reduction is needed.
      return;
    }
    CRAB_LOG("object-reduction", crab::outs() << "State Before Reduction:\n"
                                              << *this << "\n";);
    // giving an odi map domain:
    //    {obj_id |-> <eq_fld_dom, eq_reg_dom, cache_dom>}
    // eq_reg_dom and base_dom
    for (auto kv : m_odi_map) {
      // for each abstract object:
      //    there is a variable set of Symb, a set of flds and a set of regs
      //    a domain reduction from cache_dom to base_dom is:
      //  base_dom' = (\exists Symb (\exists flds. eq_fld_dom ^ cache_dom) ^
      //  eq_reg_dom) ^ base_dom
      const obj_id_t &id = kv.first;
      base_dom_variable_vector_t flds_vec;
      get_obj_flds_with_ghost(id, flds_vec);
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (prod_ref) {
        const object_info_t &prod_info = odi_map_t::object_info_val(*prod_ref);
        const small_range &num_refs = prod_info.refcount_val();
        if (num_refs.is_one() &&
            crab_domain_params_man::get()
                .singletons_in_base()) { // singleton object in the base dom
          continue; // no need to perform reduction
        }
        if (prod_info.cache_reg_loaded_val()) {
          // The reduction should only be performed out when some register are
          // loaded values from this object
          object_value_t out_prod_val = odi_map_t::object_odi_val(*prod_ref);
          object_info_t out_prod_info = prod_info;
          out_prod_info.cache_reg_loaded_val() = false;
          apply_reduction_from_object_to_base(m_base_dom, out_prod_val, id);
          m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                               std::move(out_prod_val)));
        }
      }
    }

    base_abstract_domain_t regs_only_base =
        std::move(project_onto_eq_regs(m_base_dom));
    for (auto kv : m_odi_map) {
      // for each abstract object:
      //  a domain reduction from base_dom to cache_dom,
      //  cache_dom' = (\exists Symb (\exists regs. eq_reg_dom ^ base_dom) ^
      //  eq_fld_dom) ^ cache_dom
      const obj_id_t &id = kv.first;
      base_dom_variable_vector_t flds_vec;
      get_obj_flds_with_ghost(id, flds_vec);
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (prod_ref) {
        const object_info_t &prod_info = odi_map_t::object_info_val(*prod_ref);
        const small_range &num_refs = prod_info.refcount_val();
        if (num_refs.is_one() &&
            crab_domain_params_man::get()
                .singletons_in_base()) { // singleton object && object in base
          continue; // no need to perform reduction
        }
        if (prod_info.cache_reg_stored_val()) {
          // if some fields are stored from regs
          // perform reduction from base to obj
          object_value_t out_prod_val = odi_map_t::object_odi_val(*prod_ref);
          object_info_t out_prod_info = prod_info;
          out_prod_info.cache_reg_stored_val() = false;
          apply_reduction_from_base_to_object(regs_only_base, out_prod_val, id);
          m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                                std::move(out_prod_val)));
        }
      }
    }
    CRAB_LOG("object-reduction", crab::outs() << "State After Reduction:\n"
                                              << *this << "\n";);
  }

  /// @brief method for crab intrinsic to perform reduction for a specific
  /// object \p id
  /// @param[in] id an object id
  /// @param[in] is_from_cache_to_base boolean flag indicates reduction direction
  /// @note direction is either from cache to base or from base to cache
  void apply_reduction_per_object(const obj_id_t &id,
                                  bool is_from_cache_to_base) {
    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    if (prod_ref) {
      const object_info_t &prod_info = odi_map_t::object_info_val(*prod_ref);
      const small_range &num_refs = prod_info.refcount_val();
      if (num_refs.is_one() &&
          crab_domain_params_man::get()
              .singletons_in_base()) { // singleton object && object in base
        return; // no need to perform reduction
      }
      if (is_from_cache_to_base) {
        if (prod_info.cache_reg_loaded_val()) {
          // if some registers are loaded values from this object
          // perform reduction from object to base
          object_value_t out_prod_val =
              odi_map_t::object_odi_val(*prod_ref);
          object_info_t out_prod_info = prod_info;
          out_prod_info.cache_reg_loaded_val() = false;
          apply_reduction_from_object_to_base(m_base_dom, out_prod_val, id);
          m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                               std::move(out_prod_val)));
        }
      } else {
        if (prod_info.cache_reg_stored_val()) {
          // if some fields are stored from regs
          // perform reduction from base to obj
          object_value_t out_prod_val = odi_map_t::object_odi_val(*prod_ref);
          object_info_t out_prod_info = prod_info;
          out_prod_info.cache_reg_stored_val() = false;
          base_abstract_domain_t regs_only_base = std::move(project_onto_eq_regs(m_base_dom));
          apply_reduction_from_base_to_object(regs_only_base, out_prod_val, id);
          m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                                std::move(out_prod_val)));
        }
      }
    }
  }

  /// @brief transfer "regs == flds" constriants between MRU cache and base
  /// @param[in,out] base_dom base domain
  /// @param[in,out] prod an object value stores
  ///    a product of <SUM dom, cache dom, EQ dom>
  /// @param[in] id the object id refers \p prod in the \c m_odi_map
  void apply_reduction_between_object_and_base(base_abstract_domain_t &base_dom,
                                               object_value_t &prod,
                                               const obj_id_t &id) const {
    // E.g.
    // pre state:
    //      eq_regs_dom = { x == t1 ; z == t2; k = t3 }
    //      base_dom = { 2 <= z; z <= 3; k = 0; }
    //      Object = (
    //                cache_dom = { 0 <= y; y <= 10 },
    //                eq_flds_dom = { y == t1; w == t2 }
    //               )
    // After reduction:
    //      eq_regs_dom = { x == t1 ; z == t2; k = t3 }
    //      base_dom = { k = 0; 2 <= z; z <= 3; 0 <= x; x <= 10; }
    //      Object = (
    //                cache_dom = { 0 <= y; y <= 10; 2 <= w; w <= 3; },
    //                eq_flds_dom = { y == t1; w == t2 }
    //               )
    apply_reduction_from_object_to_base(base_dom, prod, id);
    base_abstract_domain_t regs_only_base =
        std::move(project_onto_eq_regs(base_dom));
    apply_reduction_from_base_to_object(regs_only_base, prod, id);
  }

  /// @brief construct an equality map based on symbolic variables
  /// @param[in,out] map a map <symb> -> <[regs], [flds]>
  /// @param[in] eq_regs_dom equality domain for registers
  /// @param[in] eq_flds_dom equality domain for fields
  /// @param[in] flds object fields
  /// @pre map is empty
  void build_map(symb_flds_regs_map_t &map,
                 const eq_register_domain_t &eq_regs_dom,
                 const eq_fields_domain_t &eq_flds_dom,
                 const base_dom_variable_vector_t &flds) const {
    OBJECT_DOMAIN_SCOPED_STATS(".build_map");
    // invariants: for any term t that is used for equality,
    //  there \exists a fld \in flds and reg \in regs
    // such that fld |-> t, reg |-> t to represent fld == reg.

    // Worst case about performance is O(m * n) where m and n
    // are number of variables in each domain.
    for (auto &v : flds) {
      // for each term t, fill \p map by {t |-> <[flds], []> }
      // any field variable assigned by t on \p eq_flds_dom
      std::shared_ptr<const usymb_t> t_ptr = (*eq_flds_dom).get(v);
      if (!t_ptr) {
        continue;
      }
      auto it = map.find(*t_ptr);
      if (it != map.end()) {
        std::get<0>(it->second).push_back(v);
      } else {
        map.insert({*t_ptr, {{v}, {}}});
      }
    }

    // for each term t, fill \p map by {t |-> <[flds], [regs]> }
    // any register variable assigned by t on \p eq_regs_dom
    for (auto it = map.begin(); it != map.end(); ++it) {
      const usymb_t &symb = it->first;
      auto regs_opt = eq_regs_dom.get_variables(symb);
      if (regs_opt) {
        std::get<1>(it->second)
            .insert(std::get<1>(it->second).end(), regs_opt->begin(), regs_opt->end());
      }
    }
  }

  /// @brief Reduce equalities between registers and fields from cache to base
  /// @param[in,out] base_dom the base domain at the abstract state invokes this function
  /// @param[in] prod object value stores a product of <SUM dom, cache dom, EQ dom>
  /// @param id the object id refers \p prod in the \c m_odi_map
  void apply_reduction_from_object_to_base(base_abstract_domain_t &base_dom,
                                           const object_value_t &prod,
                                           const obj_id_t &id) const {

    OBJECT_DOMAIN_SCOPED_STATS(".reduction");
    // The following reduction is performed the followings:
    // a. get equalities constraints from 'eq_regs_dom' and 'eq_flds_dom'
    // b. rename a field from 'reduced_cache' to a register 
    //    if there is a register with the same symbol as that field;
    //    if there are multiple registers equal to a field, add equality constriants
    //    into 'reduced_cache'.
    // c. forget all fields in the 'reduced_cache'.
    // d. meet 'reduced_cache' with 'base_dom'
    // E.g.
    // pre state:
    //      eq_regs_dom = { k == t1; w == t2; v == t4; m == t2 }
    //      base_dom = { 3 <= a }
    //      Object = (
    //                summary_dom = {...},
    //                cache_dom = { x <= y; y == z },
    //                eq_flds_dom = { x == t1 ; y == t2; z = t3 }
    //               )
    // After reduction:
    //      eq_regs_dom = { k == t1; w == t2; v == t4; m == t2 }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { x <= y; y == z },
    //                eq_flds_dom = { x == t1 ; y == t2; z = t3 }
    //               )
    CRAB_LOG("object-reduce", crab::outs()
                                  << "Before Reduction from object to base:\n"
                                  << "base = " << base_dom << "\n"
                                  << "eq_regs = " << m_eq_regs_dom << "\n"
                                  << "odi = ";
             m_odi_map.odi_val_write(crab::outs(), prod);
             crab::outs() << "\n";);

    base_dom_variable_vector_t obj_flds;
    get_obj_flds_with_ghost(id, obj_flds);
    const eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(prod);
    const cache_domain_t &cache_dom = odi_map_t::object_cache_val(prod);
    cache_domain_t reduced_cache =
        std::move(project_onto_eq_flds(cache_dom, eq_flds_dom));

    symb_flds_regs_map_t map;
    build_map(map, m_eq_regs_dom, eq_flds_dom, obj_flds);
    CRAB_LOG("object-reduce", symbol_map_write(crab::outs(), map));

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t flds = std::get<0>(it->second);
      const base_dom_variable_vector_t regs = std::get<1>(it->second);
      if (flds.empty() || regs.empty()) {
        // lost any equalities between flds and regs. Do not perform reduction
        continue;
      }
      assert(flds.size() > 0);
      assert(regs.size() > 0);
      const base_dom_variable_t &fld = flds[0];
      auto it_regs = regs.begin();
      const base_dom_variable_t &var = *it_regs;
      while (it_regs != regs.end()) {
        if (it_regs == regs.begin()) {
          // rename fld to the first register equals to
          reduced_cache.rename({fld}, {var});
        } else {
          // for additional registers, add equality
          reduced_cache += ikos::operator==(var, *it_regs);
        }
        ++it_regs;
      }
    }
    map.clear();
    reduced_cache.forget(obj_flds);
    base_dom &= *reduced_cache;
    CRAB_LOG(
        "object-reduce", crab::outs() << "After Reduction:\n"
                                      << "base = " << base_dom << "\n"
                                      << "eq_regs = " << m_eq_regs_dom << "\n"
                                      << "odi = ";
        m_odi_map.odi_val_write(crab::outs(), prod); crab::outs() << "\n";);
  }

  /// @brief Reduce equalities between registers and fields from base to cache
  /// @param[in,out] base_dom the base domain at the abstract state invokes this function
  /// @param[in,out] prod object value stores a product of <SUM dom, cache dom, EQ dom>
  /// @param[in] id the object id refers \p prod in the \c m_odi_map
  void apply_reduction_from_base_to_object(base_abstract_domain_t &base_dom,
                                           object_value_t &prod,
                                           const obj_id_t &id) const {

    OBJECT_DOMAIN_SCOPED_STATS(".reduction");
    // The following reduction is performed the followings:
    // a. get equalities constraints from 'eq_regs_dom' and 'eq_flds_dom'
    // b. rename a register from 'base_dom' to a field
    //    if the field with the same symbol as that field;
    //    if multiple fields equal to the register, add equalities
    // c. project those regs in base domain.
    // d. meet reduced and renamed base domain with cache domain
    // E.g.
    // pre state:
    //      eq_regs_dom = { k == t1; w == t2; v == t4; m == t2; }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { z = 3 },
    //                eq_flds_dom = { x == t1 ; y == t2; }
    //               )
    // After reduction:
    //      eq_regs_dom = { k == t1; w == t2; v == t4; m == t2; }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { z = 3, x <= y },
    //                eq_flds_dom = { x == t1 ; y == t2; }
    //               )
    CRAB_LOG("object-reduce", crab::outs()
                                  << "Before Reduction from base to object:\n"
                                  << "base = " << base_dom << "\n"
                                  << "eq_regs = " << m_eq_regs_dom << "\n"
                                  << "odi = ";
             m_odi_map.odi_val_write(crab::outs(), prod);
             crab::outs() << "\n";);

    base_dom_variable_vector_t obj_flds;
    get_obj_flds_with_ghost(id, obj_flds);
    const eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(prod);
    cache_domain_t &cache_dom = odi_map_t::object_cache_val(prod);
    variable_vector_t remained_regs;
    variable_vector_t renamed_flds;

    symb_flds_regs_map_t map;
    build_map(map, m_eq_regs_dom, eq_flds_dom, obj_flds);
    CRAB_LOG("object-reduce", symbol_map_write(crab::outs(), map));

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t flds = std::get<0>(it->second);
      const base_dom_variable_vector_t regs = std::get<1>(it->second);
      if (flds.empty() || regs.empty()) {
        // lost any equalities between flds and regs. Do not perform reduction
        continue;
      }
      assert(flds.size() > 0);
      assert(regs.size() > 0);
      const base_dom_variable_t &reg = regs[0];
      auto it_flds = flds.begin();
      const base_dom_variable_t &var = *it_flds;
      while (it_flds != flds.end()) {
        if (it_flds == flds.begin()) {
          // base_dom.rename({reg}, {var});
          base_dom += ikos::operator==(var, reg);
        } else {
          base_dom += ikos::operator==(var, *it_flds);
        }
        ++it_flds;
      }
    }
    // FIXME: avoid copying, use efficient project onto
    base_abstract_domain_t reduced_base = base_dom;
    reduced_base.project(obj_flds);
    *cache_dom &= reduced_base;
    CRAB_LOG(
        "object-reduce", crab::outs() << "After Reduction:\n"
                                      << "base = " << base_dom << "\n"
                                      << "eq_regs = " << m_eq_regs_dom << "\n"
                                      << "odi = ";
        m_odi_map.odi_val_write(crab::outs(), prod); crab::outs() << "\n";);
  }

  /***************** Print operations *****************/
  void print_flds_id_map(crab_os &o) const {
    o << "Fields -> id map: ";
    if (!m_flds_id_map) {
      o << "not created";
    } else if ((*m_flds_id_map).empty()) {
      o << "empty";
    } else {
      std::unordered_map<obj_id_t, variable_vector_t> print_map;
      for (auto kv : (*m_flds_id_map)) {
        variable_t field = kv.first;
        obj_id_t id = kv.second;
        auto it = print_map.find(id);
        if (it == print_map.end()) {
          print_map.insert({id, {field}});
        } else {
          auto &s = it->second;
          auto it = std::upper_bound(s.begin(), s.end(), field);
          s.insert(it, field);
        }
      }

      // print map
      for (auto it = print_map.begin(); it != print_map.end();) {
        o << "Object " << it->first;
        object_domain_impl::print_vector(o, it->second);
        ++it;
        if (it != print_map.end()) {
          o << ", ";
        }
      }
    }
  }

  void symbol_map_write(crab_os &o, const symb_flds_regs_map_t &map) const {
    o << "reduce map: ";
    for (const auto &kv : map) {
      o << "[" << kv.first << "] = (";
      object_domain_impl::print_vector(o, kv.second.first);
      o << ", ";
      object_domain_impl::print_vector(o, kv.second.second);
      o << "),";
    }
    o << "\n";
  }

  void object_write(crab_os &o) const { // a special output for object domain
    // not using api from seperate domain
    if (m_odi_map.is_bottom()) {
      o << "Object = _|_";
    } else if (m_odi_map.is_top()) {
      o << "Object = {}";
    } else {
      for (auto it = m_odi_map.begin(); it != m_odi_map.end();) {
        obj_id_t id = it->first;
        variable_vector_t vars;
        get_obj_flds(id, vars);
        std::sort(vars.begin(), vars.end());
        o << "\nObject ";
        object_domain_impl::print_vector(o, vars);
        o << "= ";
        auto prod = it->second;
        m_odi_map.odi_write(o, prod);
        ++it;
        if (it != m_odi_map.end()) {
          o << ",";
        }
      }
    }
  }

  void dump() const { write(crab::outs()); }
  /**------------------ End helper method definitions ------------------**/

public:
  /**------------------ Begin domain API definitions ------------------**/
  object_domain_t make_top() const override { return object_domain_t(true); }

  object_domain_t make_bottom() const override {
    return object_domain_t(false);
  }

  void set_to_top() override {
    object_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    object_domain_t abs(false);
    std::swap(*this, abs);
  }

  object_domain(bool is_top = true)
      : m_is_bottom(!is_top), m_ghost_var_num_man(get_type()),
        m_ghost_var_eq_man(get_type()) {}

  object_domain(const object_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(o.m_base_dom),
        m_odi_map(o.m_odi_map), m_addrs_dom(o.m_addrs_dom),
        m_eq_regs_dom(o.m_eq_regs_dom), m_flds_id_map(o.m_flds_id_map),
        m_refs_base_addrs_map(o.m_refs_base_addrs_map),
        m_ghost_var_num_man(o.m_ghost_var_num_man),
        m_ghost_var_eq_man(o.m_ghost_var_eq_man) {
    OBJECT_DOMAIN_SCOPED_STATS(".copy");
  }

  object_domain(object_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(std::move(o.m_base_dom)),
        m_odi_map(std::move(o.m_odi_map)),
        m_addrs_dom(std::move(o.m_addrs_dom)),
        m_eq_regs_dom(std::move(o.m_eq_regs_dom)),
        m_flds_id_map(std::move(o.m_flds_id_map)),
        m_refs_base_addrs_map(std::move(o.m_refs_base_addrs_map)),
        m_ghost_var_num_man(std::move(o.m_ghost_var_num_man)),
        m_ghost_var_eq_man(std::move(o.m_ghost_var_eq_man)) {
    OBJECT_DOMAIN_SCOPED_STATS(".move");
  }

  object_domain_t &operator=(const object_domain_t &o) {
    OBJECT_DOMAIN_SCOPED_STATS(".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_base_dom = o.m_base_dom;
      m_odi_map = o.m_odi_map;
      m_addrs_dom = o.m_addrs_dom;
      m_eq_regs_dom = o.m_eq_regs_dom;
      m_flds_id_map = o.m_flds_id_map;
      m_refs_base_addrs_map = o.m_refs_base_addrs_map;
      m_ghost_var_num_man = o.m_ghost_var_num_man;
      m_ghost_var_eq_man = o.m_ghost_var_eq_man;
    }
    return *this;
  }

  object_domain_t &operator=(object_domain_t &&o) {
    OBJECT_DOMAIN_SCOPED_STATS(".move");
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_base_dom = std::move(o.m_base_dom);
      m_odi_map = std::move(o.m_odi_map);
      m_addrs_dom = std::move(o.m_addrs_dom);
      m_eq_regs_dom = std::move(o.m_eq_regs_dom);
      m_flds_id_map = std::move(o.m_flds_id_map);
      m_refs_base_addrs_map = std::move(o.m_refs_base_addrs_map);
      m_ghost_var_num_man = std::move(o.m_ghost_var_num_man);
      m_ghost_var_eq_man = std::move(o.m_ghost_var_eq_man);
    };
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    OBJECT_DOMAIN_SCOPED_STATS(".is_top");

    bool res = (!is_bottom() && m_base_dom.is_top() && m_odi_map.is_top());
    return res;
  }

  bool operator<=(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".leq");

    CRAB_LOG("object-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                        << "\n\t" << o << "\n";);
    if (is_bottom() || o.is_top()) {
      CRAB_LOG("object-leq", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) {
      CRAB_LOG("object-leq", crab::outs() << "Result2=0\n";);
      return false;
    }

    return less_than_eq(*this, o);
  }

  void operator|=(const object_domain_t &o) override {
    OBJECT_DOMAIN_SCOPED_STATS(".join");

    // Trivial cases first
    if (is_bottom()) { // this is bot, assign this by o
      if (!o.is_bottom()) {
        *this = o;
      }
      return;
    } else if (o.is_bottom()) { // o is bot, nothing change
      return;
    } else if (is_top() || o.is_top()) { // one is top, set to top
      set_to_top();
      return;
    }

    CRAB_LOG("object-join",
             crab::outs() << "Join " << *this << "\n and " << o << "\n =\n");
    self_join(o);
    CRAB_LOG("object-join", crab::outs() << "Result=" << *this << "\n");
  }

  object_domain_t operator|(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".join");

    // Trivial cases first
    if (is_bottom()) { // this is bot, return o
      return o;
    } else if (o.is_bottom()) { // o is bot, return this
      return *this;
    } else if (is_top() || o.is_top()) { // one is top, set to top
      object_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("object-join",
             crab::outs() << "Join " << *this << "\n and " << o << "\n =\n");
    object_domain_t res(
        std::move(join_or_widening(*this, o, true /*is join*/)));
    CRAB_LOG("object-join", crab::outs() << "Result=" << res << "\n");
    return res;
  }

  object_domain_t operator&(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    }

    CRAB_LOG("object-meet", crab::outs()
                                << "Meet " << *this << " and " << o << " =\n");

    object_domain_t res(
        std::move(meet_or_narrowing(*this, o, true /*is meet*/)));

    CRAB_LOG("object-meet", crab::outs() << res << "\n");
    return res;
  }

  void operator&=(const object_domain_t &o) override {
    OBJECT_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      *this = o;
      return;
    }

    CRAB_LOG("object-meet", crab::outs()
                                << "Meet " << *this << " and " << o << " =\n");

    // TODO: improve this by avoiding the copy of the left operand.
    *this = std::move(meet_or_narrowing(*this, o, true /*is meet*/));
    CRAB_LOG("object-meet", dump(); crab::outs() << "\n");
  }

  object_domain_t operator||(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("object", crab::outs()
                           << "Widening " << *this << " and " << o << " =\n");

    object_domain_t res(
        std::move(join_or_widening(*this, o, false /*is widen*/)));

    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t
  widening_thresholds(const object_domain_t &o,
                      const thresholds<number_t> &thresholds) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("object", crab::outs() << "Widening with threshold " << *this
                                    << " and " << o << " =\n");

    object_domain_t res(
        std::move(join_or_widening(*this, o, false /*is widen*/)));

    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t operator&&(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("object", crab::outs()
                           << "Narrowing " << *this << " and " << o << " =\n");

    object_domain_t res(
        std::move(meet_or_narrowing(*this, o, false /*is narrow*/)));

    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  /***************** Regions and reference operations *****************/

  /// @brief Initialize a region by a region variable \p rgn
  /// @param[in] rgn a region variable in crabIR region_init
  void region_init(const variable_t &rgn) override {
    OBJECT_DOMAIN_SCOPED_STATS(".region_init");

    ERROR_IF_NOT_REGION(rgn, __LINE__);

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    if (get_obj_id(rgn) == boost::none) {
      // if a region does not belong to an object, treat it as an object
      // i.e. treat region as a field, as well as an object id
      obj_id_t id = create_new_obj_id(rgn);
      update_fields_id_map(rgn, id);
    }
    obj_id_t id = get_obj_id_or_fail(rgn);
    // initialize information for an abstract object for object
    object_info_t obj_info = object_info_t(
        // No references owned by the object
        small_range::zero(),
        // Cache is not used
        boolean_value::get_false(),
        // Cache is not dirty
        boolean_value::get_false(), false, false);
    // construct a new odi value, default is top
    odi_domain_product_t obj_odi;
    // set object info
    odi_map_t::object_info_val(obj_odi) = obj_info;
    m_odi_map.set(id, std::move(obj_odi));
    // Assign ghost variables to rgn for modeling its content
    get_or_insert_gvars(rgn);

    CRAB_LOG("object", crab::outs() << "After region_init(" << rgn
                                    << ")=" << *this << "\n";);
  }

  /// @brief Create a new reference ref associated with as within region
  /// @param[in] ref a reference variable
  /// @param[in] rgn a region variable
  /// @param[in] size a symbolic / contant val represented allocation size
  /// @param[in] as allocation site identifier, set tag for allocation (uaf)
  void ref_make(const variable_t &ref, const variable_t &rgn,
                /* size of the allocation in bytes */
                const variable_or_constant_t &size,
                /* identifier for the allocation site */
                const allocation_site &as) override {
    OBJECT_DOMAIN_SCOPED_STATS(".ref_make");
  /*  TODO: determine singleton and non-singleton is based on the number of 
   *  interpreting ref_make for an abstract object
   *  Thus, the following CrabIR pattern cannot be precise:
   *  DSA_INFO(rgn1(int), rgn2(int));
   *  make_ref(ref1, rgn1, 8);
   *  make_ref(ref2, rgn1, 8);
   *  store_ref(ref1, rgn1, 2);
   *  -- Current implementation treat summary subdomain as a top.
   */

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);

    if (is_bottom()) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    // Assign ghost variables to ref
    ghost_variables_t ref_gvars = get_or_insert_gvars(ref);
    // Initialize ghost variables
    // If we model reference's offset and size
    // then, ref.offset = 0 and ref.size = size
    if (ref_gvars.has_offset_and_size()) {
      ref_gvars.get_offset_and_size().init(m_base_dom,
                                           rename_variable_or_constant(size));
    }

    if (auto id_opt = get_obj_id(rgn)) {
      // the object id must exist for each region if the
      // region belongs to an abstract object.

      // retrieve an abstract object info
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object goes to top
        return;
      }
      const object_info_t &obj_prod_info =
          odi_map_t::object_info_val(*prod_ref);
      const small_range &num_refs = obj_prod_info.refcount_val();

      // Check number of references for an abstract object
      if (num_refs.is_zero()) {
        // if object is singleton, add a base address for that object.
        // This is required when join / meet of two states.
        object_info_t out_obj_prod_info = obj_prod_info;
        object_value_t out_obj_prod_val = odi_map_t::object_odi_val(*prod_ref);
        out_obj_prod_info.refcount_val().increment(ref);
        m_odi_map.set(*id_opt,
                      odi_domain_product_t(std::move(out_obj_prod_info),
                                           std::move(out_obj_prod_val)));
      } else if (num_refs.is_one()) {
        // if the abstract object is a singleton object,
        // now since the number of references is increasing,
        // the state need to change the object's status to make it non-singleton
        change_object_status(*id_opt);
      }
    }

    CRAB_LOG("object", crab::outs() << "After ref_make(" << ref << "," << rgn
                                    << ":" << rgn.get_type() << "," << size
                                    << "," << as << ")=" << *this << "\n";);
  }

  /// @brief Read the content of reference ref within rgn. The content is
  ///  stored in res.
  /// @param[in] ref a reference variable
  /// @param[in] rgn a region variable used in object domain
  /// @param[in] res a register variable used in base domain
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    OBJECT_DOMAIN_SCOPED_STATS(".ref_load");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);
    // checks types,
    // the type of region variable should be consistent with the type of
    // register
    // E.g. if region represents integer, the register must be an integer type
    if ((rgn.get_type().is_bool_region() && !res.get_type().is_bool()) ||
        (rgn.get_type().is_integer_region() && !res.get_type().is_integer()) ||
        (rgn.get_type().is_real_region() && !res.get_type().is_real()) ||
        (rgn.get_type().is_reference_region() &&
         !res.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::", __func__, ": type of lhs ", res, " (",
                 res.get_type(), ") is not compatible with region ", rgn, " (",
                 rgn.get_type(), ")");
    }

    if (is_bottom()) {
      return;
    }

    // skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    const ghost_variables_t &res_gvars = get_or_insert_gvars(res);
    // use ghost variable for field
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);

    // The reference ref should not be evaluated as a null pointer.
    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_load: reference ", ref,
                                   " is null."););
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      res_gvars.forget(m_base_dom);
      m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
      return;
    }

    if (auto id_opt = get_obj_id(rgn)) {
      // the object id must exist for each region if the
      // region belongs to an abstract object.

      // retrieve an abstract object info
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object goes to top
        return;
      }
      const object_info_t &obj_info_ref = odi_map_t::object_info_val(*prod_ref);
      const small_range &num_refs = obj_info_ref.refcount_val();

      // In crab IR, the number of references cannot be zero
      //  if ref_load access a not null reference.
      // So zero case should not exist
      assert(!num_refs.is_zero());

      if (num_refs.is_one() && /*a flag to keep singleton in the base domain*/
          crab_domain_params_man::get().singletons_in_base()) {
        // singleton object in the base domain,
        // perform a strong read:
        //   i.e. an assignment on base dom for res := rgn
        res_gvars.assign(m_base_dom, rgn_gvars);
        if (rgn_gvars.has_offset_and_size() &&
            res_gvars.has_offset_and_size()) {
          rgn_gvars.get_offset_and_size().assign(
              m_base_dom, res_gvars.get_offset_and_size());
        } else if (rgn_gvars.has_offset_and_size()) {
          rgn_gvars.get_offset_and_size().forget(m_base_dom);
        }
        // forget register res in the equality register domain
        m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
      } else { // use odi map
        // keep the equality between res == rgn as memory load
        // a. generate or reuse a symbolic variable for rgn
        // b. keep that symbolic variable assigned to res
        boost::optional<usymb_t> reg_symb = boost::none;
        boost::optional<std::pair<usymb_t, usymb_t>> reg_offset_size_symb =
            boost::none;
        invalidate_cache_if_miss((*id_opt), ref, rgn, reg_symb,
                                 reg_offset_size_symb);
        // assigning register with symbolic variable
        boost::optional<ghost_variables_eq_t> res_eq_gvars = get_eq_gvars(res);
        m_eq_regs_dom.set(res_eq_gvars.value().get_var(), *reg_symb);
        if (res_eq_gvars.value().has_offset_and_size()) {
          m_eq_regs_dom.set(
              res_eq_gvars.value().get_offset_and_size().get_offset(),
              std::get<0>(*reg_offset_size_symb));
          m_eq_regs_dom.set(
              res_eq_gvars.value().get_offset_and_size().get_size(),
              std::get<1>(*reg_offset_size_symb));
        }
        // forget old property for register res
        res_gvars.forget(m_base_dom);
      }
    }

    CRAB_LOG("object", crab::outs()
                           << "After " << res << ":="
                           << "ref_load(" << rgn << ":" << rgn.get_type() << ","
                           << ref << ":" << ref.get_type() << ")=" << *this
                           << "\n";);
  }

  /// @brief Write the content of val to the address pointed by ref in region.
  /// @param[in] ref a reference variable
  /// @param[in] rgn a region variable used in object domain
  /// @param[in] val a val, could be a register or constant
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    OBJECT_DOMAIN_SCOPED_STATS(".ref_store");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_ARRAY_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);
    // checks types,
    // the type of region variable should be consistent with the type of
    // register
    if ((rgn.get_type().is_bool_region() && !val.get_type().is_bool()) ||
        (rgn.get_type().is_integer_region() && !val.get_type().is_integer()) ||
        (rgn.get_type().is_real_region() && !val.get_type().is_real()) ||
        (rgn.get_type().is_reference_region() &&
         !val.get_type().is_reference())) {
      CRAB_ERROR(domain_name(), "::", __func__, ": type of value ", val, " (",
                 val.get_type(), ") is not compatible with region ", rgn, " (",
                 rgn.get_type(), ")");
    }

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    // use ghost variable for field
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);

    // The reference ref should not be evaluated as a null pointer.
    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_store: reference ",
                                   ref, " is null."););
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      operator-=(rgn);
      return;
    }

    if (auto id_opt = get_obj_id(rgn)) {
      // the object id must exist for each region if the
      // region belongs to an abstract object.

      // retrieve an abstract object odi
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object does not exist, default is top
        return;
      }
      // retrieve an abstract object info
      const object_info_t &obj_info_ref = odi_map_t::object_info_val(*prod_ref);

      small_range num_refs = obj_info_ref.refcount_val();
      assert(!num_refs.is_zero());

      if (num_refs.is_one() && /*a flag to keep singleton in the base domain*/
          crab_domain_params_man::get().singletons_in_base()) {
        // singleton object in the base domain,
        // perform a strong update
        //   i.e. an assignment on base dom for rgn := val
        if (val.is_constant()) {
          m_base_dom.assign(rgn_gvars.get_var(), val.get_constant());
          if (val.is_reference_null() && rgn_gvars.has_offset_and_size()) {
            // if the val is NULL
            m_base_dom.assign(rgn_gvars.get_offset_and_size().get_offset(),
                              number_t(0));
            m_base_dom.assign(rgn_gvars.get_offset_and_size().get_size(),
                              number_t(0));
          }
        } else {
          ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
          rgn_gvars.assign(m_base_dom, val_gvars);
          if (rgn_gvars.has_offset_and_size() &&
              val_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().assign(
                m_base_dom, val_gvars.get_offset_and_size());
          } else if (rgn_gvars.has_offset_and_size()) {
            rgn_gvars.get_offset_and_size().forget(m_base_dom);
          }
        }
      } else { // use odi map
        // keep the equality between rgn == val as a memory store
        // a. generate or reuse a symbolic variable for val if val is variable
        // b. keep that symbolic variable assigned to rgn
        boost::optional<usymb_t> reg_symb = boost::none;
        boost::optional<std::pair<usymb_t, usymb_t>> reg_offset_size_symb =
            boost::none;
        if (val.is_variable()) { // obtain symbols of the val
          ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
          reg_symb = get_or_insert_symbol(m_eq_regs_dom, val_gvars.get_var());
          if (val_gvars.has_offset_and_size()) {
            reg_offset_size_symb = std::make_pair(
                get_or_insert_symbol(
                    m_eq_regs_dom,
                    val_gvars.get_offset_and_size().get_offset()),
                get_or_insert_symbol(
                    m_eq_regs_dom, val_gvars.get_offset_and_size().get_size()));
          }
        }
        invalidate_cache_if_miss((*id_opt), ref, rgn, reg_symb,
                                 reg_offset_size_symb);
        // retrieve an abstract object odi again since cache may be flushed
        prod_ref = m_odi_map.find((*id_opt));
        object_value_t out_prod = odi_map_t::object_odi_val(*prod_ref);
        object_info_t out_info = odi_map_t::object_info_val(*prod_ref);
        cache_domain_t &flds_dom = odi_map_t::object_cache_val(out_prod);
        eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(out_prod);
        if (val.is_constant()) {
          flds_dom.assign(rgn_gvars.get_var(), val.get_constant());
          if (val.is_reference_null() && rgn_gvars.has_offset_and_size()) {
            // if the val is NULL
            flds_dom.assign(rgn_gvars.get_offset_and_size().get_offset(),
                            number_t(0));
            flds_dom.assign(rgn_gvars.get_offset_and_size().get_size(),
                            number_t(0));
          }
          // forget region rgn in the equality field domain
          m_ghost_var_eq_man.forget(rgn, *eq_flds_dom);
        } else { // val is a variable (i.e. register)
          ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
          boost::optional<ghost_variables_eq_t> val_eq_gvars =
              get_eq_gvars(val.get_variable());
          // assigning register with symbolic variable
          m_eq_regs_dom.set(val_eq_gvars.value().get_var(), *reg_symb);
          if (val_eq_gvars.value().has_offset_and_size()) {
            m_eq_regs_dom.set(
                val_eq_gvars.value().get_offset_and_size().get_offset(),
                std::get<0>(*reg_offset_size_symb));
            m_eq_regs_dom.set(
                val_eq_gvars.value().get_offset_and_size().get_size(),
                std::get<1>(*reg_offset_size_symb));
          }
          // forget old property for region rgn
          rgn_gvars.forget(*flds_dom);
        }
        // update object info
        object_info_t out_obj_info = object_info_t(
            num_refs, boolean_value::get_true() /*Cache is used*/,
            boolean_value::get_true() /*Cache is dirty*/,
            out_info.cache_reg_loaded_val(), out_info.cache_reg_stored_val());
        // update odi map
        m_odi_map.set(*id_opt, odi_domain_product_t(std::move(out_obj_info),
                                                    std::move(out_prod)));
      }
    }

    CRAB_LOG("object", crab::outs()
                           << "After ref_store(" << rgn << ":" << rgn.get_type()
                           << "," << ref << ":" << ref.get_type() << "," << val
                           << ":" << val.get_type() << ")=" << *this << "\n";);
  }

  /// @brief Create a new reference ref2 to region rgn2. The reference ref2 is
  /// created by adding offset to ref1.
  /// @param[in] ref1 a reference variable pointed to the base address
  /// @param[in] rgn1 a region variable for the first field
  /// @param[in] ref2 a new reference variable computed by \p ref1 + \p offset
  /// @param[in] rgn2 the corresponding region referred by \p ref2
  /// @param[in] offset numerical / symbolic offset
  void ref_gep(const variable_t &ref1, const variable_t &rgn1,
               const variable_t &ref2, const variable_t &rgn2,
               const linear_expression_t &offset) override {
    OBJECT_DOMAIN_SCOPED_STATS(".ref_gep");

    ERROR_IF_NOT_REGION(rgn1, __LINE__);
    ERROR_IF_NOT_REGION(rgn2, __LINE__);
    ERROR_IF_NOT_REF(ref1, __LINE__);
    ERROR_IF_NOT_REF(ref2, __LINE__);

    if (is_bottom()) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    ghost_variables_t ref1_gvars = get_or_insert_gvars(ref1);
    ghost_variables_t ref2_gvars = get_or_insert_gvars(ref2);

    if (ref1_gvars.has_offset_and_size() && ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().assign(
          m_base_dom, ref1_gvars.get_offset_and_size(), offset);
    } else if (ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().forget(m_base_dom);
    }

    // assign equality: ref2 == ref1
    // In C memory model,
    // pointer arithmetic cannot be performed from different memory objects.
    // Thus, a precondition is ref2 and ref1 are belongs to a same memory obj.
    m_addrs_dom.add(get_or_insert_base_addr(ref1),
                    get_or_insert_base_addr(ref2));

    m_base_dom.assign(ref2_gvars.get_var(), ref1_gvars.get_var() + offset);

    CRAB_LOG("object", crab::outs()
                           << "After (" << rgn2 << "," << ref2
                           << ") := ref_gep(" << rgn1 << "," << ref1 << " + "
                           << offset << ")=" << *this << "\n";);
  }

  // Add constraints between references
  void ref_assume(const reference_constraint_t &ref_cst) override {
    OBJECT_DOMAIN_SCOPED_STATS(".ref_assume");

    if (!is_bottom()) {
      if (ref_cst.is_tautology()) {
        return;
      }
      if (ref_cst.is_contradiction()) {
        set_to_bottom();
        return;
      }

      // REDUCTION: perform reduction
      apply_reduction_based_on_flags(true);

      auto lin_cst = m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
          ref_cst, ghost_variable_kind::ADDRESS);
      m_base_dom += lin_cst;
      m_is_bottom = m_base_dom.is_bottom();
      /*
       * We cannot maintain any relational information between two pointers 
       * that refer to different memory objects. Even if two pointers point to 
       * the same memory object, we still cannot imply any relationship for 
       * non-equality operations.
       */
      if (!m_is_bottom && ref_cst.is_equality()) {
        auto offset_lin_csts =
            m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
                ref_cst, ghost_variable_kind::OFFSET);
        m_base_dom += offset_lin_csts;
        m_is_bottom = m_base_dom.is_bottom();
      }
      if (!m_is_bottom && ref_cst.is_equality()) {
        auto size_lin_csts =
            m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
                ref_cst, ghost_variable_kind::SIZE);
        m_base_dom += size_lin_csts;
        m_is_bottom = m_base_dom.is_bottom();
      }
      if (!m_is_bottom && ref_cst.is_equality()) {
        if (ref_cst.is_binary()) { // ref_cst is lhs == rhs
          const variable_t &lhs = ref_cst.lhs();
          const variable_t &rhs = ref_cst.rhs();
          m_addrs_dom.add(get_or_insert_base_addr(rhs),
                          get_or_insert_base_addr(lhs));
        }
      }
    }

    CRAB_LOG("object",
             crab::outs() << "ref_assume(" << ref_cst << ")" << *this << "\n";);
  }

  // Convert a reference to an integer variable
  void ref_to_int(const variable_t &rgn, const variable_t &ref_var,
                  const variable_t &int_var) override {

    OBJECT_DOMAIN_SCOPED_STATS(".ref_to_int");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      // We represent reference as numerical in domain
      m_base_dom.assign(int_var, ref_var);
      m_addrs_dom -= get_or_insert_base_addr(ref_var);
    }
  }

  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    OBJECT_DOMAIN_SCOPED_STATS(".int_to_ref");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      ghost_variables_t ref_gvars = get_or_insert_gvars(ref_var);
      ghost_variables_t int_gvars = get_or_insert_gvars(int_var);
      if (ref_gvars.has_offset_and_size()) {
        ref_gvars.get_offset_and_size().forget(m_base_dom);
      }

      ref_gvars.assign(m_base_dom, int_gvars);
    }
  }

  /// @brief Make a copy of a region. i.e., copy from \p rhs_rgn to \p lhs_rgn
  /// @param[in] lhs_rgn The destination region where contents are to be copied.
  /// @param[in] rhs_rgn The source region from which contents are copied.
  void region_copy(const variable_t &lhs_rgn,
                   const variable_t &rhs_rgn) override {
    OBJECT_DOMAIN_SCOPED_STATS(".region_copy");
    // region_copy is limited to object domain, we provide an object_copy
    // function to avoid copy single region each time.
    // This region copy is performed only if lhs_rgn is formed as an object by
    // itself (i.e. object with one field)
    // The rhs_rgn has no such restriction. E.g. if there is a function call
    // passed by fields of a memory object.

    /* These are ensured by well-typed Crab CFGs */
    ERROR_IF_NOT_REGION(lhs_rgn, __LINE__);
    ERROR_IF_NOT_REGION(rhs_rgn, __LINE__);
    if (lhs_rgn.get_type() != rhs_rgn.get_type()) {
      CRAB_ERROR(domain_name() + "::", __func__, ": ", lhs_rgn, ":=", rhs_rgn,
                 " with different types");
    }

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(lhs_rgn) || is_unknown_region(rhs_rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    const ghost_variables_t &base_lhs = get_or_insert_gvars(lhs_rgn);
    const ghost_variables_t &base_rhs = get_or_insert_gvars(rhs_rgn);

    if (auto rhs_id_opt = get_obj_id(rhs_rgn)) {
      // the rhs region belongs to some object
      // the lhs region is the copied one

      auto lhs_id_opt = get_obj_id(lhs_rgn);
      // get obj id for lhs_rgn
      obj_id_t lhs_id = lhs_id_opt ? *lhs_id_opt : create_new_obj_id(lhs_rgn);
      variable_vector_t lhs_flds;
      get_obj_flds(lhs_id, lhs_flds);
      // if the field is one, copy it through this transfer function
      // for other cases, the copy is performed by object_copy
      if (lhs_flds.size() > 1) {
        return;
      }
      update_fields_id_map(lhs_rgn, lhs_id);

      // retrieve rhs abstract object info
      const odi_domain_product_t *prod_ref = m_odi_map.find(*rhs_id_opt);
      if (!prod_ref) { // object goes to top
        return;
      }
      object_info_t obj_info_ref = odi_map_t::object_info_val(*prod_ref);
      object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);

      const small_range &num_refs = obj_info_ref.refcount_val();

      // In crab IR, the number of references cannot be zero
      // So zero case should not exist
      assert(!num_refs.is_zero());

      if (num_refs.is_one() &&
          crab_domain_params_man::get().singletons_in_base()) {
        // singleton object in the base domain
        base_lhs.assign(m_base_dom, base_rhs);
      } else { // num_refs > 1, non-singleton, use odi map
        // copy field
        base_lhs.forget(odi_map_t::object_sum_raw_val(prod_value_ref));
        base_lhs.forget(odi_map_t::object_cache_raw_val(prod_value_ref));
        base_lhs.assign(odi_map_t::object_sum_raw_val(prod_value_ref),
                        base_rhs);
        base_lhs.assign(odi_map_t::object_cache_raw_val(prod_value_ref),
                        base_rhs);
        std::shared_ptr<usymb_t> t_ptr =
            odi_map_t::object_eq_raw_val(prod_value_ref).get(rhs_rgn);
        if (t_ptr) {
          boost::optional<ghost_variables_eq_t> lhs_rgn_eq_gvars =
              get_eq_gvars(lhs_rgn);
          boost::optional<ghost_variables_eq_t> rhs_rgn_eq_gvars =
              get_eq_gvars(rhs_rgn);
          lhs_rgn_eq_gvars.value().assign(
              odi_map_t::object_eq_raw_val(prod_value_ref),
              rhs_rgn_eq_gvars.value());
        }
        m_addrs_dom.expand(get_or_insert_base_addr(*rhs_id_opt),
                           get_or_insert_base_addr(lhs_id));
      }
      m_odi_map.set(lhs_id, odi_domain_product_t(std::move(obj_info_ref),
                                                 std::move(prod_value_ref)));
    }

    CRAB_LOG("object", crab::outs()
                           << "After region_copy(" << lhs_rgn << ":"
                           << lhs_rgn.get_type() << "," << rhs_rgn << ":"
                           << rhs_rgn.get_type() << ")=" << *this << "\n";);
  }

  // Cast between regions of different types
  void region_cast(const variable_t &src_rgn,
                   const variable_t &dst_rgn) override {
    // A region_cast is used to cast unknown regions to typed regions,
    // or viceversa.
    OBJECT_DOMAIN_SCOPED_STATS(".region_cast");

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(dst_rgn)) {
      return;
    }

    if (is_unknown_region(src_rgn)) {
      // create a fresh obj for dst_rgn
      obj_id_t id = create_new_obj_id(dst_rgn);
      // if a region does not belong to an object, treat it as an object
      // treat region as a field, as well as an object id
      update_fields_id_map(dst_rgn, id);
      // initialize information for an abstract object for object
      m_odi_map.set(
          id, odi_domain_product_t(object_info_t(
                                       // No references owned by the object
                                       small_range::zero(),
                                       // Cache is not used
                                       boolean_value::get_false(),
                                       // Cache is not dirty
                                       boolean_value::get_false(),
                                       false, false),
                                   object_value_t()));
    }

    CRAB_LOG("object", crab::outs()
                           << "After region_cast(" << src_rgn << ":"
                           << src_rgn.get_type() << "," << dst_rgn << ":"
                           << dst_rgn.get_type() << ")=" << *this << "\n";);
  }

  // Remove a reference ref within region reg
  // TODO: object domain does not handles reference free. Thus, it cannot
  // check UAF.
  void ref_free(const variable_t &reg, const variable_t &ref) override {}

  // This default implementation is expensive because it will call the
  // join.
  // DEFAULT_SELECT_REF(object_domain_t)
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    OBJECT_DOMAIN_SCOPED_STATS(".select_ref");
    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      ghost_variables_t lhs_ref_gvars = get_or_insert_gvars(lhs_ref);
      base_dom_variable_t lhs_b_rgn =
          is_unknown_region(lhs_rgn) ? lhs_rgn
                                     : get_or_insert_gvars(lhs_rgn).get_var();
      ghost_variables_t cond_gvars = get_or_insert_gvars(cond);
      base_dom_variable_or_constant_t b_ref1 =
          rename_variable_or_constant(ref1);
      base_dom_variable_or_constant_t b_ref2 =
          rename_variable_or_constant(ref2);
      boost::optional<base_dom_variable_t> b_rgn1 =
          rename_variable_optional(rgn1);
      boost::optional<base_dom_variable_t> b_rgn2 =
          rename_variable_optional(rgn2);

      m_base_dom.select_ref(lhs_ref_gvars.get_var(), lhs_b_rgn, cond, b_ref1,
                            b_rgn1, b_ref2, b_rgn2);
    }
  }

  /**************************** Numerical operations *************************/
  // x := y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());
      m_eq_regs_dom -= x;
    }
  }

  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), k);
      m_eq_regs_dom -= x;
    }
  }

  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());
      m_eq_regs_dom -= x;
    }
  }

  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), k);
      m_eq_regs_dom -= x;
    }
  }

  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply(op, get_or_insert_gvars(dst).get_var(),
                       get_or_insert_gvars(src).get_var());
      m_eq_regs_dom -= dst;
    }
  }

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    OBJECT_DOMAIN_SCOPED_STATS(".select");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      auto b_e1 = m_ghost_var_num_man.rename_linear_expr(e1);
      auto b_e2 = m_ghost_var_num_man.rename_linear_expr(e2);
      auto b_cond = m_ghost_var_num_man.rename_linear_cst(cond);
      m_base_dom.select(get_or_insert_gvars(lhs).get_var(), b_cond, b_e1, b_e2);
      m_eq_regs_dom -= lhs;
    }
  }

  // x := e
  void assign(const variable_t &x, const linear_expression_t &e) override {
    OBJECT_DOMAIN_SCOPED_STATS(".assign");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      auto b_e = m_ghost_var_num_man.rename_linear_expr(e);
      m_base_dom.assign(get_or_insert_gvars(x).get_var(), b_e);
      m_eq_regs_dom -= x;
    }
  }

  // join(*this, copy_of_this(x := e))
  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
  }

  bool entails(const linear_constraint_t &rhs) const override { return false; }

  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {}

  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {}

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    OBJECT_DOMAIN_SCOPED_STATS(".add_constraints");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags(true);

      auto b_csts = m_ghost_var_num_man.rename_linear_cst_sys(csts);
      m_base_dom += b_csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  /********************** Boolean operations **********************/

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    OBJECT_DOMAIN_SCOPED_STATS(".assign_bool_cst");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      auto b_rhs = m_ghost_var_num_man.rename_linear_cst(rhs);
      m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);
      m_eq_regs_dom -= lhs;
    }
  }

  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {

    OBJECT_DOMAIN_SCOPED_STATS(".assign_bool_var");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.assign_bool_var(get_or_insert_gvars(lhs).get_var(),
                                 get_or_insert_gvars(rhs).get_var(),
                                 is_not_rhs);
      m_eq_regs_dom -= lhs;
    }
  }

  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    OBJECT_DOMAIN_SCOPED_STATS(".apply_binary_bool");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.apply_binary_bool(op, get_or_insert_gvars(x).get_var(),
                                   get_or_insert_gvars(y).get_var(),
                                   get_or_insert_gvars(z).get_var());
      m_eq_regs_dom -= x;
    }
  }

  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {
    OBJECT_DOMAIN_SCOPED_STATS(".assume_bool");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags(true);

      m_base_dom.assume_bool(get_or_insert_gvars(v).get_var(), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
      m_eq_regs_dom -= v;
    }
  }

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {

    OBJECT_DOMAIN_SCOPED_STATS(".select_bool");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      m_base_dom.select_bool(get_or_insert_gvars(lhs).get_var(),
                             get_or_insert_gvars(cond).get_var(),
                             get_or_insert_gvars(b1).get_var(),
                             get_or_insert_gvars(b2).get_var());
      m_eq_regs_dom -= lhs;
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {

    OBJECT_DOMAIN_SCOPED_STATS(".assign_bool_ref_cst");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      apply_reduction_based_on_flags();

      auto rhs_lin_cst = m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
          rhs, ghost_variable_kind::ADDRESS);
      m_base_dom.assign_bool_cst(lhs, rhs_lin_cst);
      m_eq_regs_dom -= lhs;
    }
  }

  /********************** Array operations **********************/
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(object_domain_t)

  // FIXME: The followings are UNDEFINED METHODS

  /********************** Backward numerical operations **********************/
  // x = y op z
  // Substitute x with y op z in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const object_domain_t &invariant) override {}
  // x = y op k
  // Substitute x with y op k in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const object_domain_t &invariant) override {}
  // x = e
  // Substitute x with e in the abstract value
  // The result is meet with invariant.
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const object_domain_t &invariant) override {}

  /********************** Backward boolean operations **********************/
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const object_domain_t &invariant) override {}
  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const object_domain_t &invariant) override {
  }
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const object_domain_t &invariant) override {}
  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const object_domain_t &invariant) override {}

  /********************** Miscellaneous operations **********************/

  // Normalize the abstract domain if such notion exists.
  void normalize() override {}

  // Reduce the size of the abstract domain representation.
  void minimize() override {}

  // Make a new copy of var without relating var with new_var
  void expand(const variable_t &var, const variable_t &new_var) override {

    if (is_bottom() || is_top()) {
      return;
    }

    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const object_domain_t &invariant) override {}

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::", __func__,
               " not "
               "implemented");
  }

  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {

    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {

    CRAB_ERROR(domain_name(), "::", __func__, " not implemented");
  }

  // FIXME: The above methods are UNDEFINED METHODS

  /// @brief Forget v
  /// @param v the program variable that its property need to be forgot
  /// @note Forgetting a region means forgetting its properties, but not
  /// removing it from the scope.
  void operator-=(const variable_t &v) override {
    OBJECT_DOMAIN_SCOPED_STATS(".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    // To forget a variable, it is easier to forget if variable is not a rgn
    // However, if it is a region, we need to perform the followings:
    // 1. the abstract object that region var belongs to is a singleton,
    //    forget it in base dom if we keep singletons in the base domain.
    // 2. object is not singleton, forget it in odi map
    if (v.get_type().is_region()) {
      // for now skip analysis for unknown region
      if (is_unknown_region(v)) {
        return;
      }
      if (auto id_opt = get_obj_id(v)) {
        // retrieve an abstract object info
        const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
        if (!prod_ref) {
          return;
        }
        object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
        object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);

        const small_range &num_refs = prod_info_ref.refcount_val();
        variable_vector_t obj_flds;
        get_obj_flds(*id_opt, obj_flds);

        if (obj_flds.size() == 1) { // The region is part of an object that has only one region
          m_odi_map -= *id_opt;
        } else {
          // Check number of references for an abstract object
          if (num_refs.is_one() &&
              crab_domain_params_man::get()
                  .singletons_in_base()) { // singleton object in the base
            m_ghost_var_num_man.forget(v, m_base_dom);
          } else { // use odi map
            m_ghost_var_num_man.forget(
                v, odi_map_t::object_sum_raw_val(prod_value_ref));
            m_ghost_var_num_man.forget(
                v, odi_map_t::object_cache_raw_val(prod_value_ref));
            m_ghost_var_eq_man.forget(
                v, odi_map_t::object_eq_raw_val(prod_value_ref));
          }
          m_odi_map.set(*id_opt, odi_domain_product_t(std::move(prod_info_ref),
                                                      std::move(prod_value_ref)));
        }
      }
    } else { // forget a non region variable
      m_ghost_var_num_man.forget(v, m_base_dom);
      if (v.get_type().is_reference()) {
        m_addrs_dom.operator-=(get_or_insert_base_addr(v));
      }
      m_ghost_var_eq_man.forget(v, m_eq_regs_dom);
    }
  }

  /// @brief forget a set of variables
  /// @param variables a set of program variables that their property need to be forgot
  /// @note same functionality as \c operator-=
  void forget(const variable_vector_t &variables) override {
    OBJECT_DOMAIN_SCOPED_STATS(".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    CRAB_LOG("object-forget", crab::outs() << "Forgetting (";
             object_domain_impl::print_vector(crab::outs(), variables);
             crab::outs() << ")=" << *this << "\n";);

    variable_vector_t non_odi_vars;
    variable_vector_t ref_vars;
    non_odi_vars.reserve(variables.size());
    ref_vars.reserve(variables.size());
    // The following map keeps fields that need to be remained
    std::unordered_map<obj_id_t, variable_vector_t> flds_by_id_map;

    // filtering incoming input variable vectors based on their types
    // since we split them into several subdomains
    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        if (is_unknown_region(v)) { // skip unknown regions
          continue;
        }
        if (auto id_opt = get_obj_id(v)) {
          auto it = flds_by_id_map.find(*id_opt);
          if (it != flds_by_id_map.end()) {
            (it->second).push_back(v);
          } else {
            flds_by_id_map.insert({(*id_opt), {v}});
          }
        }
      } else {
        non_odi_vars.push_back(v);
        if (v.get_type().is_reference()) {
          ref_vars.push_back(get_or_insert_base_addr(v));
        }
      }
    }

    for (auto kv : flds_by_id_map) {
      const obj_id_t &id = kv.first;
      const variable_vector_t &flds_to_forget = kv.second;
      variable_vector_t obj_flds;
      get_obj_flds(id, obj_flds);
      if (obj_flds.size() == flds_to_forget.size()) {
        // if forget all fields, forget this object
        m_odi_map -= id;
      } else {
        // retrieve an abstract object info
        const odi_domain_product_t *prod_ref = m_odi_map.find(id);
        if (!prod_ref) {
          continue;
        }
        object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
        const small_range &num_refs = prod_info_ref.refcount_val();
        object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
        if (num_refs.is_one() && crab_domain_params_man::get()
                                    .singletons_in_base()) { // singleton object in the base
          non_odi_vars.insert(non_odi_vars.end(), flds_to_forget.begin(),
                              flds_to_forget.end());
        } else if (obj_flds.size() !=
                  flds_to_forget.size()) { // non-singleton object
          // partially forget field(s)
          for (auto v : flds_to_forget) {
            m_ghost_var_num_man.forget(
                v, odi_map_t::object_sum_raw_val(prod_value_ref));
            m_ghost_var_num_man.forget(
                v, odi_map_t::object_cache_raw_val(prod_value_ref));
            m_ghost_var_eq_man.forget(
                v, odi_map_t::object_eq_raw_val(prod_value_ref));
          }
        }
        m_odi_map.set(id, odi_domain_product_t(std::move(prod_info_ref),
                                               std::move(prod_value_ref)));
      }
    }

    // forget registers, references
    for (auto v : non_odi_vars) {
      m_ghost_var_num_man.forget(v, m_base_dom);
      m_ghost_var_eq_man.forget(v, m_eq_regs_dom);
    }
    // forget references in equality domain
    m_addrs_dom.forget(ref_vars);

    CRAB_LOG("object-forget", crab::outs()
                                  << "After Froget:" << *this << "\n";);
  }

  /// @brief projecting the dimension of the abstract state onto the scope of \p variables
  /// @param variables a set of program variables
  /// @post After projections, the state will only keep properties w.r.t \p variables
  void project(const variable_vector_t &variables) override {
    OBJECT_DOMAIN_SCOPED_STATS(".project");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    CRAB_LOG("object-project", crab::outs() << "Projecting (";
             object_domain_impl::print_vector(crab::outs(), variables);
             crab::outs() << ")=" << *this << "\n";);

    variable_vector_t non_odi_vars;
    variable_vector_t ref_vars;
    non_odi_vars.reserve(variables.size());
    ref_vars.reserve(variables.size());
    // The following map keeps fields that need to be remained
    std::unordered_map<obj_id_t, variable_vector_t> flds_by_id_map;
    // a temporary map to keep results
    odi_map_t out_odi_map;

    // filtering incoming input variable vectors based on their types
    // since we split them into several subdomains
    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        if (is_unknown_region(v)) { // skip unknown regions
          continue;
        }
        if (auto id_opt = get_obj_id(v)) {
          auto it = flds_by_id_map.find(*id_opt);
          if (it != flds_by_id_map.end()) {
            (it->second).push_back(v);
          } else {
            flds_by_id_map.insert({(*id_opt), {v}});
          }
        }
      } else {
        non_odi_vars.push_back(v);
        if (v.get_type().is_reference()) {
          ref_vars.push_back(get_or_insert_base_addr(v));
        }
      }
    }

    // projecting fields need to reconstruct odi map
    // keep odi that only be remained
    for (auto kv : flds_by_id_map) {
      const obj_id_t &id = kv.first;
      const variable_vector_t &flds = kv.second;
      ref_vars.push_back(get_or_insert_base_addr(id));
      // retrieve an abstract object info
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (!prod_ref) {
        continue;
      }
      object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
      object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
      const small_range &num_refs = prod_info_ref.refcount_val();
      if (num_refs.is_one() && crab_domain_params_man::get()
                                   .singletons_in_base()) { // singleton object
        non_odi_vars.insert(non_odi_vars.end(), flds.begin(), flds.end());
      } else { // non-singleton object
        // project based on corresponding field(s)
        m_ghost_var_num_man.project(
            flds, odi_map_t::object_sum_raw_val(prod_value_ref));
        m_ghost_var_num_man.project(
            flds, odi_map_t::object_cache_raw_val(prod_value_ref));
        m_ghost_var_eq_man.project(
            flds, odi_map_t::object_eq_raw_val(prod_value_ref));
      }
      out_odi_map.set(id, odi_domain_product_t(std::move(prod_info_ref),
                                               std::move(prod_value_ref)));
    }
    m_ghost_var_num_man.project(non_odi_vars, m_base_dom);
    m_addrs_dom.project(ref_vars);
    m_ghost_var_eq_man.project(non_odi_vars, m_eq_regs_dom);
    std::swap(m_odi_map, out_odi_map);
    CRAB_LOG("object-project", crab::outs()
                                   << "After Projection:" << *this << "\n";);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    OBJECT_DOMAIN_SCOPED_STATS(".rename");

    if (is_bottom() || is_top()) {
      return;
    }

    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::", __func__, " different lengths");
    }

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    variable_vector_t from_non_odi_vars;
    variable_vector_t to_non_odi_vars;
    variable_vector_t from_ref_vars;
    variable_vector_t to_ref_vars;
    from_non_odi_vars.reserve(from.size());
    to_non_odi_vars.reserve(from.size());
    from_ref_vars.reserve(from.size());
    to_ref_vars.reserve(from.size());
    // The following map keeps fields that need to be renamed
    std::unordered_map<obj_id_t,
                       std::pair<variable_vector_t, variable_vector_t>>
        renamed_flds_by_id_map;

    // filtering incoming input variable vectors based on their types
    // since we split them into several subdomains
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &old_v = from[i];
      const variable_t &new_v = to[i];
      if (old_v.get_type() != new_v.get_type()) {
        CRAB_ERROR(domain_name(), "::", __func__, " ", old_v, " and ", new_v, " must preserve the same type");
      }
      if (old_v.get_type().is_region()) {
        if (auto id_opt = get_obj_id(old_v)) {
          auto it = renamed_flds_by_id_map.find(*id_opt);
          if (it != renamed_flds_by_id_map.end()) {
            std::get<0>(it->second).push_back(old_v);
            std::get<1>(it->second).push_back(new_v);
          } else {
            renamed_flds_by_id_map.insert({(*id_opt), {{old_v}, {new_v}}});
          }
        }
      } else {
        from_non_odi_vars.push_back(from[i]);
        to_non_odi_vars.push_back(to[i]);
        if (old_v.get_type().is_reference()) {
          from_ref_vars.push_back(get_or_insert_base_addr(old_v));
          to_ref_vars.push_back(get_or_insert_base_addr(to[i]));
        }
      }
    }

    m_addrs_dom.rename(from_ref_vars, to_ref_vars);
    m_ghost_var_eq_man.rename(from_non_odi_vars, to_non_odi_vars,
                              m_eq_regs_dom);

    // renaming fields
    for (auto kv : renamed_flds_by_id_map) {
      const obj_id_t &id = kv.first;
      const variable_vector_t &from_flds = std::get<0>(kv.second);
      const variable_vector_t &to_flds = std::get<1>(kv.second);
      // retrieve an abstract object info
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (!prod_ref) {
        continue;
      }
      object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
      object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
      const small_range &num_refs = prod_info_ref.refcount_val();
      // update base_dom or odi_map
      if (num_refs.is_one()) { // singleton object
        m_ghost_var_num_man.rename(from_flds, to_flds, m_base_dom);
      } else { // non-singleton object
        // project based on corresponding field(s)
        m_ghost_var_num_man.rename(
            from_flds, to_flds, odi_map_t::object_sum_raw_val(prod_value_ref));
        m_ghost_var_num_man.rename(
            from_flds, to_flds,
            odi_map_t::object_cache_raw_val(prod_value_ref));
        m_ghost_var_eq_man.rename(from_flds, to_flds,
                                  odi_map_t::object_eq_raw_val(prod_value_ref));
      }
      m_odi_map.set(id, odi_domain_product_t(std::move(prod_info_ref),
                                             std::move(prod_value_ref)));
    }
  }

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {
    OBJECT_DOMAIN_SCOPED_STATS(".to_interval");
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      if (v.get_type().is_region()) {
        CRAB_ERROR(domain_name(), "::", __func__,
            " extracting value of a region is not supported.");
      }
      return m_base_dom[v];
    }
  }

  interval_t at(const variable_t &v) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".to_interval");
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      if (v.get_type().is_region()) {
        CRAB_ERROR(domain_name(), "::", __func__,
            " extracting value of a region is not supported.");
      }
      return m_base_dom.at(v);
    }
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      return linear_constraint_t::get_false();
    } else if (is_top()) {
      return linear_constraint_t::get_true();
    } else {
      linear_constraint_system_t out_csts =
          m_base_dom.to_linear_constraint_system();
      return out_csts;
    }
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    OBJECT_DOMAIN_SCOPED_STATS(".intrinsic");

    //=================================================================//
    //       Special intrinsics supported by the object domain
    //=================================================================//
    // ---DSA region analysis---
    //      This analysis indicates which regions might belong to the
    //      same memory object. The intrinstic is added only if object
    //      has more than one field.
    // TODO: need to support other analysis
    auto error_if_not_variable =
        [&, func = __func__](const variable_or_constant_t &vc) {
          if (!vc.is_variable()) {
            CRAB_ERROR(domain_name(), "::", func, " ", name,
                       " expected a variable input");
          }
        };
    auto error_if_not_constant =
        [&, func = __func__](const variable_or_constant_t &vc) {
          if (!vc.is_constant()) {
            CRAB_ERROR(domain_name(), "::", func, " ", name,
                       " expected a constant input");
          }
        };
    auto error_if_not_bool = [&, func = __func__](const variable_t &var) {
      if (!var.get_type().is_bool()) {
        CRAB_ERROR(domain_name(), "::", func, " ", name, " parameter ", var,
                   " should be Bool");
      }
    };
    auto error_if_not_rgn = [&,
                             func = __func__](const variable_or_constant_t &x) {
      if (!x.is_variable() || !x.get_type().is_region()) {
        // the input vector should only contains region variables
        CRAB_ERROR(domain_name(), "::", func, " ", name, " parameter ", x,
                   " should be a region");
      }
    };
    auto error_if_not_ref = [&, func = __func__](const variable_t &var) {
      if (!var.get_type().is_reference()) {
        CRAB_ERROR(domain_name(), "::", func, " ", name, " parameter ", var,
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

    // REDUCTION: perform reduction
    apply_reduction_based_on_flags();

    if (name == "regions_from_memory_object") {
      // pass region varaibles into object field map
      // the intrinsics is only added in clam if the object has more than one
      // region.
      assert(inputs.size() >= 1);
      error_if_not_rgn(inputs[0]);
      if (inputs[0].get_type().is_unknown_region()) {
        return;
      }
      auto obj_id_opt = get_obj_id(inputs[0].get_variable());
      obj_id_t obj_id = obj_id_opt ? *obj_id_opt : inputs[0].get_variable();
      for (int i = 1, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        // Note that, region initialization is before the intrinsic calls
        // Any obj info set up are not for object id will be removed.
        auto old_id_opt = get_obj_id(inputs[i].get_variable());
        if (old_id_opt && obj_id != *old_id_opt) {
          m_odi_map -= *old_id_opt;
          m_addrs_dom -= get_or_insert_base_addr(*old_id_opt);
          update_fields_id_map(inputs[i].get_variable(), obj_id);
        }
      }
      // No need to update odi map
    } else if (name == "do_reduction") {
      assert(inputs.size() == 3);
      error_if_not_rgn(inputs[0]);
      error_if_not_variable(inputs[1]);
      error_if_not_ref(inputs[1].get_variable());
      error_if_not_constant(inputs[2]);
      variable_t rgn = inputs[0].get_variable();
      variable_t ref = inputs[1].get_variable();
      bool reduce_direction = inputs[2].is_bool_true();
      if (auto id_opt = get_obj_id(rgn)) {
        if (test_ref_refer_mru_object(ref, *id_opt)) {
          // current reference refers the mru object
          apply_reduction_per_object(*id_opt, reduce_direction);
        }
      }
    } else if (name == "is_dereferenceable") {
      if (crab_domain_params_man::get().region_is_dereferenceable()) {
        assert(inputs.size() == 3);
        assert(outputs.size() == 1);
        // ignore region variable (inputs[0])
        error_if_not_variable(inputs[1]);
        error_if_not_bool(outputs[0]);
        variable_t bv = outputs[0];
        variable_t ref = inputs[1].get_variable();
        CRAB_LOG("object-is-deref",
                 crab::outs() << bv << ":= is_dereferenceable(" << inputs[0]
                              << "," << inputs[1] << "," << inputs[2] << ")\n"
                              << *this << "\n";);
        if (auto ref_gvars_opt = get_num_gvars(ref)) {
          if ((*ref_gvars_opt).has_offset_and_size()) {
            CRAB_LOG("object-is-deref",
                     crab::outs()
                         << "\t"
                         << "Reference " << ref << "\n"
                         << "\toffset="
                         << (*ref_gvars_opt).get_offset_and_size().get_offset()
                         << "\n"
                         << "\tsize="
                         << (*ref_gvars_opt).get_offset_and_size().get_size()
                         << "\n";);
            if ((*ref_gvars_opt)
                    .get_offset_and_size()
                    .is_deref(m_base_dom,
                              rename_variable_or_constant(inputs[2]))) {
              set_bool_var_to_true(bv);
              CRAB_LOG("object-is-deref", crab::outs() << "\tRESULT=TRUE\n");
              return;
            }
          }
        }
        CRAB_LOG("object-is-deref", crab::outs() << "\tRESULT=UNKNOWN\n");
        operator-=(bv);
      }
    } else if (name == "commit_cache") {
      assert(inputs.size() >= 1);
      auto obj_id_opt = get_obj_id(inputs[0].get_variable());
      std::unordered_set<obj_id_t> id_set;
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          return;
        }
        // indicate which object needs to commit
        if (auto id_opt = get_obj_id(inputs[i].get_variable())) {
          if (id_set.find(*id_opt) == id_set.end()) {
            id_set.insert(*id_opt);
          }
        }
      }
      // commit the cache for each object
      for (auto &id : id_set) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(id);
        if (!prod_ref) { // object goes to top
          continue;
        }
        object_info_t obj_info = odi_map_t::object_info_val(*prod_ref);
        object_value_t obj_value = odi_map_t::object_odi_val(*prod_ref);
        commit_cache_if_dirty(obj_value, obj_info, id);
        // reset cache
        obj_info.cacheused_val() = boolean_value::get_false();
        obj_info.cachedirty_val() = boolean_value::get_false();
        m_odi_map.set(id, odi_domain_product_t(std::move(obj_info),
                                               std::move(obj_value)));
      }
    } else if (name == "copy_memory_object") {
      assert(inputs.size() >= 1);
      assert(outputs.size() >= 1);
      std::vector<variable_t> input_vars;
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          return;
        }
        // indicate which object needs to commit
        input_vars.push_back(inputs[i].get_variable());
      }
      copy_object(input_vars, outputs, *this, *this);
    }
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    if (is_bottom()) {
      return boolean_value::bottom();
    }

    if (!ref.get_type().is_reference()) {
      return boolean_value::get_false();
    }

    if (auto gvars_opt = get_num_gvars(ref)) {
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

  std::string domain_name() const override {
    field_abstract_domain_t fld_dom;
    const char *base_prefix = "Base:";
    std::string base_name = m_base_dom.domain_name();
    const char *object_prefix = ", Object:";
    std::string odi_name = m_odi_map.domain_name();
    const char *addr_prefix = ", Addrs:";
    std::string addr_name = m_addrs_dom.domain_name();
    const char *prefix = "ObjectDomain";
    std::string name;
    unsigned len = base_name.size() + odi_name.size() + addr_name.size() + 41;
    name.reserve(len);
    name.append(prefix);
    name.append("(");
    name.append(base_prefix);
    name.append(base_name);
    name.append(object_prefix);
    name.append(odi_name);
    name.append(addr_prefix);
    name.append(addr_name);
    name.append(")");
    return name;
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("object-print", o << "("
                                 << "Flds-id map=";
               print_flds_id_map(o); o << ",\n"
                                       << "BaseDom=";
               m_ghost_var_num_man.write(o, m_base_dom);
               o << ","
                 << "Addrs=" << m_addrs_dom;
               o << ","
                 << "Regs=" << m_eq_regs_dom;
               o << ","; object_write(o); o << ")\n"; return;);
      o << "Base = ";
      m_ghost_var_num_man.write(o, m_base_dom);
      o << ",\neq_addrs = ";
      m_addrs_dom.write(o);
      o << ",\neq_regs = ";
      m_ghost_var_eq_man.write(o, m_eq_regs_dom);
      o << ",";
      object_write(o);
    }
  }
  /**------------- End domain API definitions -------------------**/

  /**----------- Begin domain Inter definitions -----------------**/
  using equiv_class_regions_t = std::vector<variable_t>;
  // a vector of classes where each class contains a vector of regions
  using classes_t = std::vector<equiv_class_regions_t>;

  /// @brief copy object based on one dsa node to another
  /// @param[in] src_rgns a set of regions from src dsa node
  /// @param[in] dst_rgns a set of regions from dst dsa node
  /// @param[in] src_dom the source abstract state
  /// @param[in,out] dst_dom the destinated abstract state
  void copy_object(const equiv_class_regions_t &src_rgns,
                   const equiv_class_regions_t &dst_rgns,
                   const object_domain_t &src_dom,
                   object_domain_t &dst_dom) const {
    CRAB_LOG("object-copy", crab::outs() << "Copying ";
             object_domain_impl::print_vector(crab::outs(), src_rgns);
             crab::outs() << " -> ";
             object_domain_impl::print_vector(crab::outs(), dst_rgns);
             crab::outs() << "\n src: " << src_dom << "\ndst: " << dst_dom
                          << "\n";);
    if (src_rgns.size() == 0 || src_rgns.size() != dst_rgns.size()) {
      return;
    }
    // the regions might contain unknown region which are not considered,
    // filter them out
    equiv_class_regions_t src_rgns_no_unknown, dst_rgns_no_unknown;
    src_rgns_no_unknown.reserve(src_rgns.size());
    dst_rgns_no_unknown.reserve(dst_rgns.size());
    for (unsigned i = 0, len = src_rgns.size(); i < len; ++i) {
      if (src_dom.is_unknown_region(src_rgns[i]) &&
          dst_dom.is_unknown_region(dst_rgns[i])) {
        continue;
      }
      src_rgns_no_unknown.push_back(src_rgns[i]);
      dst_rgns_no_unknown.push_back(dst_rgns[i]);
    }

    if (src_rgns_no_unknown.size() == 0 ||
        src_rgns_no_unknown.size() != dst_rgns_no_unknown.size()) {
      return;
    }
    // precondition: the rgns from source is formed as some abstract object
    obj_id_t src_id = src_dom.get_obj_id_or_fail(src_rgns_no_unknown[0]);
    auto dst_id_opt = dst_dom.get_obj_id(dst_rgns_no_unknown[0]);
    if (dst_id_opt == boost::none) {
      dst_id_opt = dst_dom.create_new_obj_id(dst_rgns_no_unknown[0]);
      for (auto fld : dst_rgns) {
        dst_dom.update_fields_id_map(fld, *dst_id_opt);
      }
    }
    const odi_domain_product_t *prod_ref = src_dom.m_odi_map.find(src_id);
    if (!prod_ref) {
      // default is top
      return;
    }
    object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
    object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
    const small_range &num_refs = prod_info_ref.refcount_val();

    assert(!num_refs.is_zero());

    if (num_refs.is_one() &&
        crab_domain_params_man::get().singletons_in_base()) {
      base_abstract_domain_t tmp = src_dom.m_base_dom;
      for (unsigned i = 0, len = src_rgns_no_unknown.size(); i < len; ++i) {
        tmp.expand(src_rgns_no_unknown[i], dst_rgns_no_unknown[i]);
      }
      dst_dom.m_ghost_var_num_man.project(dst_rgns_no_unknown, tmp);
      dst_dom.m_base_dom &= tmp;
    } else { // num_refs > 1
      for (unsigned i = 0, len = src_rgns_no_unknown.size(); i < len; ++i) {
        odi_map_t::object_sum_val(prod_value_ref)
            .expand(src_rgns_no_unknown[i], dst_rgns_no_unknown[i]);
        odi_map_t::object_cache_val(prod_value_ref)
            .expand(src_rgns_no_unknown[i], dst_rgns_no_unknown[i]);
        odi_map_t::object_eq_val(prod_value_ref)
            .expand(src_rgns_no_unknown[i], dst_rgns_no_unknown[i]);
      }
      dst_dom.m_ghost_var_num_man.project(
          dst_rgns_no_unknown, odi_map_t::object_sum_raw_val(prod_value_ref));
      dst_dom.m_ghost_var_num_man.project(
          dst_rgns_no_unknown, odi_map_t::object_cache_raw_val(prod_value_ref));
      dst_dom.m_ghost_var_eq_man.project(
          dst_rgns_no_unknown, odi_map_t::object_eq_raw_val(prod_value_ref));
    }
    if (auto src_base_addr_opt = src_dom.get_base_addr(src_id)) {
      dst_dom.m_addrs_dom.expand(
          *src_base_addr_opt, dst_dom.get_or_insert_base_addr(*dst_id_opt));
    }
    dst_dom.m_odi_map.set(*dst_id_opt,
                          odi_domain_product_t(std::move(prod_info_ref),
                                               std::move(prod_value_ref)));
    CRAB_LOG("object-copy", crab::outs() << "After copying\n";
             crab::outs() << "dst:" << dst_dom << "\n";);
  }

  /// @brief Restrict operation. Return a new abstract state for the entry of
  /// the callee.
  /// @param[in] callsite the call site information, including dsa info
  /// @param[in] caller_dom abstract state at the caller before the call.
  void callee_entry(const callsite_info<variable_t> &callsite,
                    const object_domain_t &caller_dom) override {
    OBJECT_DOMAIN_SCOPED_STATS(".callee_entry");

    auto rgn_in_group = [](const variable_t &rgn,
                           const classes_t &cls) -> bool {
      if (!rgn.get_type().is_region()) {
        return false;
      }
      bool res = false;
      for (auto &c : cls) {
        res = std::find(c.begin(), c.end(), rgn) != c.end();
        if (res) {
          return true;
        }
      }
      return res;
    };

    if (caller_dom.is_bottom()) {
      set_to_bottom();
      return;
    } else if (caller_dom.is_top()) {
      // default: current entry state is top
      return;
    }

    const variable_vector_t &caller_in_params = callsite.get_caller_in_params();
    const variable_vector_t &callee_in_params = callsite.get_callee_in_params();
    const classes_t &caller_actual_cls = callsite.get_caller_in_groups();
    const classes_t &callee_formal_cls = callsite.get_callee_in_groups();

    // Copy caller domain
    object_domain_t caller(caller_dom);
    // Project onto actual parameters
    caller.project(caller_in_params);
    CRAB_LOG("inter-restrict", crab::outs()
                                   << "Inv at the caller: " << caller << "\n");

    // 1. propagate from actual to formal parameters
    //    (excluding regions belonging to dsa)
    for (unsigned i = 0, sz = caller_in_params.size(); i < sz; ++i) {
      const variable_t &formal = callee_in_params[i];
      const variable_t &actual = caller_in_params[i];

      if (!(formal == actual)) {
        // following DSA analysis assumption:
        // the region in the callee cannot be mapped to two distinct regions in the caller.
        // This mean there is no case that callee's regions belong to a dsa node
        // but caller's do not (picture 8 in SEADSA-SAS17 paper).
        // For the case that a src region belonging to a dsa node but dst does not
        // This differs from copying objects, we implemented it in region_copy.
        // Here, we use object_copy for both src regions and dst regions belong
        // to some dsa nodes.
        if (rgn_in_group(formal, callee_formal_cls) && rgn_in_group(actual, caller_actual_cls)) {
          continue;
        }
        CRAB_LOG("inter-restrict", crab::outs()
                                       << "\t" << formal << ":"
                                       << formal.get_type() << " and " << actual
                                       << ":" << actual.get_type() << "\n";);
        inter_transformers_impl::unify(caller, formal, actual);
        if (caller.is_bottom()) {
          CRAB_ERROR("Obtained bottom after unification");
        }
      }
    }

    // 2. propagate abstract objects
    // unify objects
    for (unsigned i = 0, len = callee_formal_cls.size(); i < len; ++i) {
      const equiv_class_regions_t &formal_rgns = callee_formal_cls[i];
      const equiv_class_regions_t &actual_rgns = caller_actual_cls[i];
      copy_object(actual_rgns, formal_rgns, caller_dom, caller);
    }
    CRAB_LOG("inter-restrict",
             crab::outs() << "Inv after formal/actual unification: " << caller
                          << "\n";);

    // 3. Meet with current entry state
    operator&=(caller);
    CRAB_LOG("inter-restrict",
             crab::outs() << "Inv after meet with callee  " << *this << "\n";);

    // 4. Project onto **input** formal parameters
    project(callee_in_params);
    CRAB_LOG(
        "inter-restrict",
        crab::outs() << "Inv at the callee after projecting onto formals: ";
        for (auto &v
             : callee_in_params) { crab::outs() << v << ";"; } crab::outs()
        << "\n"
        << *this << "\n";);
  }

  /// @brief extend operation. Return a new abstract state at the caller after
  /// the call.
  /// @param callsite the call site information, including dsa info
  /// @param callee_dom abstract state at the exit block of the callee but
  /// already projected onto formal parameters of the function.
  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const object_domain_t &callee_dom) override {
    OBJECT_DOMAIN_SCOPED_STATS(".caller_continuation");

    auto rgn_in_group = [](const variable_t &rgn,
                           const classes_t &cls) -> bool {
      if (!rgn.get_type().is_region()) {
        return false;
      }
      bool res = false;
      for (auto &c : cls) {
        res = std::find(c.begin(), c.end(), rgn) != c.end();
        if (res) {
          return true;
        }
      }
      return res;
    };

    if (is_bottom()) {
      return;
    }

    if (callee_dom.is_bottom()) {
      set_to_bottom();
      return;
    }

    object_domain_t callee_at_exit(callee_dom);

    const variable_vector_t &caller_in_params = callsite.get_caller_in_params();
    const variable_vector_t &callee_in_params = callsite.get_callee_in_params();
    const variable_vector_t &caller_out_params =
        callsite.get_caller_out_params();
    const variable_vector_t &callee_out_params =
        callsite.get_callee_out_params();
    const classes_t &caller_actual_cls = callsite.get_caller_in_groups();
    const classes_t &callee_formal_cls = callsite.get_callee_in_groups();
    const classes_t &caller_lhs_cls = callsite.get_caller_out_groups();
    const classes_t &callee_ret_cls = callsite.get_callee_out_groups();
    std::vector<bool> cls_check(caller_actual_cls.size(), false);

    auto rgn_in_group_and_which_subgroup =
        [&cls_check](const variable_t &rgn, const classes_t &cls) -> bool {
      if (!rgn.get_type().is_region()) {
        return false;
      }
      bool res = false;
      for (unsigned i = 0; i < cls.size(); ++i) {
        res = std::find(cls[i].begin(), cls[i].end(), rgn) != cls[i].end();
        if (res) {
          cls_check[i] = true;
          return res;
        }
      }
      return res;
    };

    // 1. forget output parameters at the callsite
    // e.g reassignment in case
    forget(caller_out_params);

    // 2. propagate from callee's outputs to caller's lhs of the callsite
    for (unsigned i = 0, e = callee_out_params.size(); i < e; ++i) {
      const variable_t &out_formal = callee_out_params[i];
      const variable_t &out_actual = caller_out_params[i];
      if (!(out_formal == out_actual)) {
        if (rgn_in_group(out_actual, caller_lhs_cls) && rgn_in_group(out_formal, callee_ret_cls)) {
          continue;
        }
        CRAB_LOG("inter-extend", crab::outs()
                                     << "Unifying output " << out_actual
                                     << ":= " << out_formal << "\n";);
        inter_transformers_impl::unify(callee_at_exit, out_actual, out_formal);
      }
    }

    // 3. propagate from callee's inputs to caller's inputs at callsite
    // This is needed to propagate up new input-output relationships:
    // (1) e.g. V1, V2 = call function_f(V2), where V1 is a new variable,
    // V2 is input and ouput variable;
    // def function_f(V3)
    // In case we have relations between V1 and V2. We need to unify formal
    // parameter V3 with V2.
    // V2 is a variable killed at the callsite
    // (2) is that we pass a variable into callee and this appears
    // both caller and callee's parameters list
    // (3) for any formal parameters (read only)

    // For (1). find callsite input parameters are killed at the callsite
    // Parameters that appear both as callsite argument and callee's
    // formal parameter so they shouldn't be forgotten.
    std::set<variable_t> cs_in_args(caller_in_params.begin(),
                                    caller_in_params.end());
    // killed_cs_in_args stores all variables killed at the callsite
    // those have been unified at step 2.
    auto tmp = inter_transformers_impl::set_intersection(caller_in_params,
                                                         caller_out_params);
    std::set<variable_t> killed_cs_in_args(
        tmp.begin(), tmp.end()); // linear becasue tmp is sorted

    std::vector<variable_t> caller_in_out_params;
    caller_in_out_params.reserve(callee_in_params.size() +
                                 callee_out_params.size());
    for (unsigned i = 0, e = callee_in_params.size(); i < e; ++i) {
      const variable_t &in_formal = callee_in_params[i];
      const variable_t &in_actual = caller_in_params[i];
      if (cs_in_args.find(in_formal) != cs_in_args.end()) {
        // For (2). find all shared inputs and do not forget them
        caller_in_out_params.push_back(in_formal);
      } else if (killed_cs_in_args.find(in_actual) == killed_cs_in_args.end()) {
        // For any variables in case (3), perform unfication to keep new
        // relational properties between outputs in the callee and inputs in the
        // caller
        CRAB_LOG("inter-extend", crab::outs() << "Unifying input " << in_actual
                                              << ":=" << in_formal << "\n";);
        if (rgn_in_group_and_which_subgroup(in_actual, caller_actual_cls)) {
          continue;
        }
        inter_transformers_impl::unify(callee_at_exit, in_actual, in_formal);
      }
    }

    // 4. propagate abstract objects
    // unify out objects
    for (unsigned i = 0, len = caller_lhs_cls.size(); i < len; ++i) {
      const equiv_class_regions_t &lhs_rgns = caller_lhs_cls[i];
      const equiv_class_regions_t &ret_rgns = callee_ret_cls[i];
      copy_object(ret_rgns, lhs_rgns, callee_dom, *this);
    }

    // unify in objects in case (3)
    for (unsigned i = 0, len = caller_actual_cls.size(); i < len; ++i) {
      if (cls_check[i]) {
        const equiv_class_regions_t &formal_rgns = callee_formal_cls[i];
        const equiv_class_regions_t &actual_rgns = caller_actual_cls[i];
        copy_object(formal_rgns, actual_rgns, callee_dom, *this);
      }
    }

    // 5. Forget callee's in and out parameters
    std::vector<variable_t> callee_forget_params;
    callee_forget_params.reserve(callee_in_params.size() +
                                 callee_out_params.size());
    // callee's ins + outs
    callsite.get_all_callee_params(callee_forget_params);
    caller_in_out_params.insert(caller_in_out_params.end(),
                                caller_out_params.begin(),
                                caller_out_params.end());
    callee_forget_params = inter_transformers_impl::set_difference(
        callee_forget_params, caller_in_out_params);
    callee_at_exit.forget(callee_forget_params);
    CRAB_LOG(
        "inter-extend", crab::outs() << "Forgotten all callee parameters {";
        for (auto const &v
             : callee_forget_params) { crab::outs() << v << ";"; } crab::outs()
        << "}\n";);

    // 6. Meet with the callee_at_exit after unification
    operator&=(callee_at_exit);

    CRAB_LOG("inter-extend", crab::outs() << *this << "\n";);
  }

  /**------------- End domain Inter definitions -----------------**/

  // WARN: the following APIs have strong coupling between object domain and odi
  // map domain. These methods are exposed to ODI map domain.

  /***************** Reg-fld domain operations *****************/
  /// @brief get or create a symbolic variable in fld-reg equality domain
  /// @param eq_dom an equality domain value
  /// @param variable a field
  /// @return a symbolic variable such as #var1
  usymb_t get_or_insert_symbol(eq_register_domain_t &eq_dom,
                               const flds_dom_variable_t &variable) {
    std::shared_ptr<usymb_t> symb_ptr = eq_dom.get(variable);
    if (symb_ptr) {
      return *symb_ptr;
    } else {
      return symbolic_variable_equiality_domain_impl::make_fresh_var_symbol<
          symbolic_variable_equiality_domain_impl::symbolic_var>();
    }
  }

  /// @brief a helper method to project singleton object in the base domain
  /// @param id the object id
  /// @return a temporary domain value which only keeps properties for object id
  base_abstract_domain_t project_singleton(const obj_id_t &id) const {
    base_dom_variable_vector_t flds_vec;
    get_obj_flds_with_ghost(id, flds_vec);

    base_abstract_domain_t singleton_base(m_base_dom);
    singleton_base.project(flds_vec);
    return singleton_base;
  }

  /***************** Cache operations *****************/
  /// @brief commit cache contents into summary if it is dirty, perform
  /// reduction to a specific base domain and reset cache
  /// @param prod object value stores a product of <SUM dom, cache dom, EQ dom>
  /// @param obj_info_ref the object info stores reference count & cache status
  /// @param id an object id
  void commit_cache_if_dirty(base_abstract_domain_t &base_dom,
                             object_value_t &prod,
                             const object_info_t &obj_info_ref,
                             const obj_id_t &id) const {
    // (1) perform reduction by transfering regs-flds' constriants between cache
    //  domain and base domain
    apply_reduction_between_object_and_base(base_dom, prod, id);
    m_odi_map.commit_cache_if_dirty(obj_info_ref, prod);
  }

  /// @brief commit cache contents into summary if it is dirty, perform
  /// reduction and reset cache
  /// @param prod object value stores a product of <SUM dom, cache dom, EQ dom>
  /// @param obj_info_ref the object info stores reference count & cache status
  /// @param id an object id to indicate which object
  void commit_cache_if_dirty(object_value_t &prod,
                             const object_info_t &obj_info_ref,
                             const obj_id_t &id) {
    commit_cache_if_dirty(m_base_dom, prod, obj_info_ref, id);
    m_is_bottom = m_base_dom.is_bottom();
  }

  /// @brief a wrapper method for adpaint interface of domain reduction from odi map to base
  /// @param[in] prod the abstract value corresponding to object \p id
  /// @param[in] id the object id
  void apply_reduction_from_object_to_base(const object_value_t &prod, const obj_id_t &id) {
    apply_reduction_from_object_to_base(m_base_dom, prod, id);
  }

  /// @brief a wrapper method for adpaint interface of domain reduction from base to odi map
  /// @param[in,out] prod the abstract value corresponding to object \p id
  /// @param[in] id the object id
  void apply_reduction_from_base_to_object(object_value_t &prod,
                                           const obj_id_t &id) const {
    base_abstract_domain_t regs_only_base =
        std::move(project_onto_eq_regs(m_base_dom));
    apply_reduction_from_base_to_object(regs_only_base, prod, id);
  }

}; // class object_domain

template <typename Params>
struct abstract_domain_traits<object_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab