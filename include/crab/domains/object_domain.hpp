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

#include <crab/domains/array_adaptive.hpp>
#include <crab/domains/flat_boolean_domain.hpp>
#include <crab/domains/object/object_info.hpp>
#include <crab/domains/object/odi_map_domain.hpp>
#include <crab/domains/region/ghost_variable_manager.hpp>
#include <crab/domains/region/tags.hpp>

#include <set>
#include <unordered_map>
#include <unordered_set>

////////////////////////////////////////////////////////////////////////
/// Abstract domain for regions and references based on memory abstraction.
///
/// This domain is based on the memory abstraction described in the paper
/// "Automatic Inference of Relational Object Invariants" (VMCAI'25).
///
/// The gap: when a heap-allocated object is updated field-by-field, its
/// invariant (e.g. len <= cap) is temporarily broken between writes.
/// Standard allocation-site abstraction forces all field updates to be
/// weak -- joining the new value with all previously summarized objects --
/// which causes the inferred invariant to degrade to top. E.g. memory domain in
/// region_domain.hpp is based on smash abstraction, where weak updates occur.
///
/// Our memory abstraction, called recency-use abstraction, addresses this
/// by organizing memory into banks, each with a cache that holds the
/// single most recently used (MRU) object. Field updates on the cached
/// object are handled as strong updates without affecting the summary
/// invariant of the other objects in the same bank. When the analysis
/// moves to a different object, the cache is packed back into the
/// summary. This isolates temporary invariant violations inside the
/// cache and preserves the summary invariant for all other objects.
////////////////////////////////////////////////////////////////////////

namespace crab {
namespace domains {
namespace object_domain_impl {

// FIXME: to get a bool value by a variable should be an API like entails
// the following code is fragile if we change configuration defined in CLAM
template <typename> struct base_is_instance_of_Flat_Boolean {
  static constexpr int value = 0;
};
template <typename Dom>
struct base_is_instance_of_Flat_Boolean<flat_boolean_numerical_domain<Dom>> {
  static constexpr int value = 1;
  using abstract_domain_t = flat_boolean_numerical_domain<Dom>;
  using dom_variable_t = typename abstract_domain_t::variable_t;
  static boolean_value get_bool_val_by_var(abstract_domain_t &abs,
                                           const dom_variable_t &var) {
    return abs.first().get_bool(var);
  }
};
template <typename Dom>
struct base_is_instance_of_Flat_Boolean<
    array_adaptive_domain<flat_boolean_numerical_domain<Dom>>> {
  static constexpr int value = 2;
  using abstract_domain_t =
      array_adaptive_domain<flat_boolean_numerical_domain<Dom>>;
  using dom_variable_t = typename abstract_domain_t::variable_t;
  static boolean_value get_bool_val_by_var(abstract_domain_t &abs,
                                           const dom_variable_t &var) {
    return abs.get_content_domain().get_content_domain().first().get_bool(var);
  }
};
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

/// print \p vec to \p o as [a, b, c]
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
// This domain is based on the partitioning memory model,
// an object, more specifically an abstract object, is an object represents
// a number of concrete objects following the same dsa node.
// An object is spawn from a set of fields (or regions in crabIR).
// In addition, we model the memory architecture with a cache.
// The cache tracks the most recently used (MRU) object by memory load / store
// where the MRU object is one concrete object from the summarized objects.
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
  // The equality domain keeps singleton classes: a field or register may carry
  // a symbol without a partner yet, which is how a reg-fld link is established.
  using eq_domain_value_t = symbolic_variable_equality_domain<
      typename base_abstract_domain_t::number_t,
      typename base_abstract_domain_t::varname_t, SVEQNoNormalizeParams>;

private:
  using eq_register_domain_t = eq_domain_value_t;
  using usymb_t = typename eq_domain_value_t::class_id_t;

  // type name for address domain
  using address_abstract_domain_t = eq_domain_value_t;
  using addr_value_domain_t = typename eq_domain_value_t::class_id_t;

  using symb_flds_regs_map_t =
      std::unordered_map<usymb_t, std::pair<base_dom_variable_vector_t,
                                            base_dom_variable_vector_t>>;

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

  // ODI map: object id -> <object_info, <summary, cache, field-eq>>;
  // see odi_map_domain.hpp for the layout and object_info.hpp for the info
  // tuple. The object id is a designated field variable of the object.
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
  // This map is constructed for address domain, we keep track
  // the variables created for references.
  // The map is constructed dynamically but shared across all states
  using refs_base_addrs_map_t =
      std::shared_ptr<std::unordered_map<variable_t, variable_t>>;
  using ghost_rgn_map_t =
      std::shared_ptr<std::unordered_map<variable_t, variable_t>>;

  // selector for the general binary lattice operation (see combine())
  using combine_kind = object_domain_impl::combine_kind;
  /**------------------ End type definitions ------------------**/
  // FIXME: the current solution does not support different dom operations
  static_assert(
      std::is_same<base_abstract_domain_t, field_abstract_domain_t>::value,
      "base_abstract_domain_t and field_abstract_domain_t must be the same");

  static_assert(
      std::is_same<typename Params::number_t,
                   typename Params::base_abstract_domain_t::number_t>::value,
      "Number type and BaseAbsDom::number_t must be the same");
  // This is a strong requirement
  static_assert(
      std::is_same<typename Params::base_varname_t,
                   typename Params::varname_allocator_t::varname_t>::value,
      "BaseAbsDom::varname_t and allocator_varname_t must be the same");

  static_assert(object_domain_impl::base_is_instance_of_Flat_Boolean<
                    base_abstract_domain_t>::value > 0,
                "base domain should be a flat boolean numerical domain");

  /**------------------ Begin class field definitions ------------------**/

  // The object domain is bottom iff this flag is set. Transfer functions that
  // detect an inconsistent subdomain update this flag; is_bottom() does not
  // consult the subdomains.
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
  // Invariant: old entries are kept on rename -- the map is shared across
  // states and other states may still use the old names.
  // WARN: for now, we did not cover unknown regions.
  obj_flds_id_map_t m_flds_id_map;

  // Map reference variables to base addresses variables
  // The base address is used to represent memory allocation address for a
  // concrete memory object. In addition, each cache that used for an abstract
  // object has a special base address variable for the MRU object
  refs_base_addrs_map_t m_refs_base_addrs_map;
  ghost_rgn_map_t m_ghost_rgn_map;

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
  // For each abstract object the tracked reference count is 0, 1 or
  // [1, +oo] (singleton = exactly one reference).
  // NOTE: unknown regions are ALWAYS skipped (loads through them havoc
  // the destination). Unlike region_domain, this domain has no
  // dynamic-type machinery, so the region.skip_unknown_regions=false
  // setting is not supported here.
  // TODO: re-introduce a flag that keeps a SINGLETON object's fields
  // directly in the base domain for extra precision. The previous
  // object.singletons_in_base option was removed because the odi merge
  // operators never supported it (joins aborted); bringing it back
  // requires implementing the singleton merge cases in
  // odi_map_domain::combine_op.
  /**------------------ End class field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  // Constructor Definition
  object_domain(base_abstract_domain_t &&base_dom, odi_map_t &&odi_map,
                address_abstract_domain_t &&addrs_dom,
                eq_register_domain_t &&eq_regs_dom,
                obj_flds_id_map_t &&flds_id_map,
                refs_base_addrs_map_t &&refs_base_addrs_map,
                ghost_rgn_map_t &&ghost_rgn_map,
                ghost_var_num_man_t &&ghost_var_num_man,
                ghost_var_eq_man_t &&ghost_var_eq_man)
      : m_is_bottom(base_dom.is_bottom()), m_base_dom(std::move(base_dom)),
        m_odi_map(std::move(odi_map)), m_addrs_dom(std::move(addrs_dom)),
        m_eq_regs_dom(std::move(eq_regs_dom)),
        m_flds_id_map(std::move(flds_id_map)),
        m_refs_base_addrs_map(std::move(refs_base_addrs_map)),
        m_ghost_rgn_map(std::move(ghost_rgn_map)),
        m_ghost_var_num_man(std::move(ghost_var_num_man)),
        m_ghost_var_eq_man(std::move(ghost_var_eq_man)) {}

  static std::function<variable_type(const variable_t &)> var_type_fn() {
    return [](const variable_t &v) { return v.get_type(); };
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

  /// @brief merge two shared side maps: keep the non-null side and, when
  /// the DSA info diverged, merge the right entries into it
  /// @note states share these maps; the merge mutates the kept (left) map
  template <typename SharedMapPtr>
  static SharedMapPtr merge_shared_map(const SharedMapPtr &l,
                                       const SharedMapPtr &r) {
    SharedMapPtr out = l ? l : r;
    if (l && r && l != r) {
      out->insert(r->begin(), r->end());
    }
    return out;
  }

  /// @brief performs the join operation between current abstract state and the
  /// state \p right passed by
  /// @param[in] right an abstract state serves as an input for the join
  /// operation
  /// @note Perform *this = join(*this, right)
  /// @pre \p this, \p right are two different abstract states
  /// @post the new \p this satisfies: \gamma(the new \p this) \leq^{concrete}
  /// \gamma(join( \p this, \p right ))
  void self_join(const object_domain_t &right) {

    // The join is pairwise on the subdomains: any pending cache/base
    // reduction is applied on both sides first, then the odi map, base,
    // address and register-equality domains are joined component by
    // component. The per-object case analysis (including the
    // singletons-in-base mode) lives in odi_map_domain's merge operators.

    // The DSA node info arrives via an intrinsic at the beginning of each
    // function in the crab IR. Ideally it would be static, but it is
    // obtained gradually, so missing entries are filled by merging the
    // shared side maps.

    m_flds_id_map = merge_shared_map(m_flds_id_map, right.m_flds_id_map);
    m_refs_base_addrs_map =
        merge_shared_map(m_refs_base_addrs_map, right.m_refs_base_addrs_map);
    m_ghost_rgn_map = merge_shared_map(m_ghost_rgn_map, right.m_ghost_rgn_map);

    apply_reduction_based_on_flags(true);
    auto r_opt = right.cow_apply_reduction();
    const object_domain_t &new_right =
        boost::get_optional_value_or(r_opt, right);
    m_odi_map.compound_join(new_right.m_odi_map);

    m_base_dom |= new_right.m_base_dom;

    m_addrs_dom |= new_right.m_addrs_dom;

    m_eq_regs_dom |= new_right.m_eq_regs_dom;
  }

  /// @brief the general binary lattice operation between two abstract
  /// states: join, widening (optionally with thresholds), meet or narrowing
  /// @param[in] left an abstract state serves as an input for the operation
  /// @param[in] right an abstract state serves as an input for the operation
  /// @param[in] kind selects which operation to perform
  /// @param[in] ts non-null only for widening with thresholds
  /// @return the new abstract state computed by \p left and \p right
  /// @pre \p left, \p right are two different, non-bottom abstract states
  /// @post for join / widening, \p \return satisfies:
  /// \gamma(join( \p left, \p right )) \leq^{concrete} \gamma( \p \return );
  /// for meet / narrowing, \p \return satisfies:
  /// \meet^{concrete}(\gamma( \p left ), \gamma( \p right )) \leq^{concrete}
  /// \gamma( \p \return )
  object_domain_t combine(const object_domain_t &left,
                          const object_domain_t &right, combine_kind kind,
                          const thresholds<number_t> *ts = nullptr) const {

    // The combine is pairwise on the subdomains: any pending cache/base
    // reduction is applied on both sides first, then the odi map, base,
    // address and register-equality domains are combined component by
    // component. The per-object case analysis (including the
    // singletons-in-base mode) lives in odi_map_domain's merge operators.
    const bool is_grow =
        kind == combine_kind::JOIN || kind == combine_kind::WIDENING;

    obj_flds_id_map_t out_flds_id_map =
        merge_shared_map(left.m_flds_id_map, right.m_flds_id_map);
    refs_base_addrs_map_t out_refs_base_addrs_map = merge_shared_map(
        left.m_refs_base_addrs_map, right.m_refs_base_addrs_map);
    ghost_rgn_map_t out_ghost_rgn_map =
        merge_shared_map(left.m_ghost_rgn_map, right.m_ghost_rgn_map);

    auto l_opt = left.cow_apply_reduction();
    auto r_opt = right.cow_apply_reduction();
    const object_domain_t &new_left = boost::get_optional_value_or(l_opt, left);
    const object_domain_t &new_right =
        boost::get_optional_value_or(r_opt, right);

    // combine one subdomain pair; \p dom_ts is null for the subdomains
    // that keep plain widening even when thresholds are given
    auto combine_dom = [kind](const auto &l, const auto &r,
                              const thresholds<number_t> *dom_ts) {
      if (kind == combine_kind::JOIN) {
        return l | r;
      } else if (kind == combine_kind::WIDENING) {
        return dom_ts ? l.widening_thresholds(r, *dom_ts) : l || r;
      } else if (kind == combine_kind::MEET) {
        return l & r;
      } else {
        return l && r;
      }
    };

    odi_map_t out_odi_map =
        new_left.m_odi_map.combine(new_right.m_odi_map, kind, ts);

    base_abstract_domain_t out_base_dom(
        combine_dom(new_left.m_base_dom, new_right.m_base_dom, ts));

    // an inconsistent per-object meet or base meet makes the state bottom
    // (a join / widening of two non-bottom states cannot be bottom)
    if (!is_grow && (out_odi_map.is_bottom() || out_base_dom.is_bottom())) {
      object_domain_t res;
      res.set_to_bottom();
      return res;
    }

    address_abstract_domain_t out_addrs_dom(
        combine_dom(new_left.m_addrs_dom, new_right.m_addrs_dom, nullptr));

    eq_register_domain_t out_eq_regs_dom(
        combine_dom(new_left.m_eq_regs_dom, new_right.m_eq_regs_dom, nullptr));

    ghost_var_num_man_t out_ghost_var_num_man(new_left.m_ghost_var_num_man);
    ghost_var_eq_man_t out_ghost_var_eq_man(new_left.m_ghost_var_eq_man);

    object_domain_t res(
        std::move(out_base_dom), std::move(out_odi_map),
        std::move(out_addrs_dom), std::move(out_eq_regs_dom),
        std::move(out_flds_id_map), std::move(out_refs_base_addrs_map),
        std::move(out_ghost_rgn_map), std::move(out_ghost_var_num_man),
        std::move(out_ghost_var_eq_man));
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
  /// @return if we find object id, returns; otherwise, abort execution with an
  /// error
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

  /// @brief the id representing \p v's abstract object; currently the
  /// field variable itself
  obj_id_t create_new_obj_id(const variable_t &v) { return v; }

  /// @brief get all object fields by giving an id
  /// @param[in] id an object id
  /// @param[in, out] obj_flds output: the fields of \p id
  void get_obj_flds(const obj_id_t &id, variable_vector_t &obj_flds) const {
    assert(m_flds_id_map);
    for (const auto &kv : (*m_flds_id_map)) {
      if (kv.second == id) {
        obj_flds.push_back(kv.first);
      }
    }
  }

  /// @brief get all object fields with ghost variables by giving an id
  /// @param[in] id an object id
  /// @param[in, out] obj_dom_flds output: the base-domain ghosts
  /// @note this function has different specs from \c get_obj_flds
  /// The ghost variables here could be the base address, offset, and size
  /// for a region with reference types
  /// e.g. In crab IR, if a region V_3:region(ref), the ghost variables are:
  ///    V_3.address, V_3.offset, V_3.size
  /// \c get_raw_obj_flds will fetch all ghost variables used by object fields
  void get_raw_obj_flds(const obj_id_t &id,
                        base_dom_variable_vector_t &obj_dom_flds) const {

    auto get_obj_ghost_flds = [&](const obj_id_t &id,
                                  ghost_variable_vector_t &obj_ghost_flds) {
      variable_vector_t obj_flds;
      get_obj_flds(id, obj_flds);
      for (const auto &v : obj_flds) {
        if (is_unknown_region(v)) {
          continue;
        }
        auto gvars_opt = get_num_gvars(v);
        obj_ghost_flds.push_back(*gvars_opt);
      }
    };
    ghost_variable_vector_t obj_ghost_flds;
    get_obj_ghost_flds(id, obj_ghost_flds);
    for (const auto &v : obj_ghost_flds) {
      if (v.has_offset_and_size()) {
        obj_dom_flds.push_back(v.get_offset_and_size().get_offset());
        obj_dom_flds.push_back(v.get_offset_and_size().get_size());
      }
      obj_dom_flds.push_back(v.get_var());
    }
  }

  void get_raw_obj_write_flds_map(
      const obj_id_t &id,
      std::unordered_map<base_dom_variable_t, base_dom_variable_t> &map) const {
    variable_vector_t obj_flds;
    get_obj_flds(id, obj_flds);
    for (const auto &v : obj_flds) {
      if (is_unknown_region(v)) {
        continue;
      }
      auto w_rgn_opt = get_write_region(v);
      if (w_rgn_opt) {
        auto gvars_opt = get_num_gvars(v);
        auto g_w_vars_opt = get_num_gvars(*w_rgn_opt);
        if (gvars_opt->has_offset_and_size()) {
          map.insert({g_w_vars_opt->get_offset_and_size().get_offset(),
                      gvars_opt->get_offset_and_size().get_offset()});
          map.insert({g_w_vars_opt->get_offset_and_size().get_size(),
                      gvars_opt->get_offset_and_size().get_size()});
        }
        map.insert({g_w_vars_opt->get_var(), gvars_opt->get_var()});
      }
    }
  }

  /// @brief based on each dsa intrinsic in the crab IR, group region variables
  /// as fields for each abstract object
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
  /// @brief get or mint the ghost variable `<key>_<suffix>` recorded in a
  /// lazily-created shared map
  /// @param[in,out] map the shared map, created on first use
  /// @param[in] key the variable the ghost is derived from
  /// @param[in] suffix the name suffix of the minted ghost
  /// @param[in] kind the type kind of the minted ghost
  /// @param[in] width the bitwidth of the minted ghost
  /// @return the ghost variable stored in (or added to) \p map
  static variable_t get_or_insert_ghost(
      std::shared_ptr<std::unordered_map<variable_t, variable_t>> &map,
      const variable_t &key, const std::string &suffix, variable_type_kind kind,
      unsigned width) {
    if (!map) {
      map = std::make_shared<std::unordered_map<variable_t, variable_t>>();
    }
    auto it = map->find(key);
    if (it != map->end()) {
      return it->second;
    }
    varname_t key_name = key.name();
    varname_t ghost_name =
        key_name.get_var_factory().get_or_insert_varname(key_name, suffix);
    variable_t ghost(ghost_name, kind, width);
    map->insert({key, ghost});
    return ghost;
  }

  /// @brief look up the ghost variable minted for \p key, if any
  static boost::optional<variable_t> get_ghost(
      const std::shared_ptr<std::unordered_map<variable_t, variable_t>> &map,
      const variable_t &key) {
    if (!map) {
      return boost::none;
    }
    auto it = map->find(key);
    if (it == map->end()) {
      return boost::none;
    }
    return it->second;
  }

  /// @brief get or create a variable representing the base address by given a
  /// reference variable
  /// @param[in] v the reference variable or an object id
  /// @return a ghost variable represent its base address
  /// specifically, for abstract object, the base address is used to indicate
  /// the mru object stored on cache subdomain. For mru object, we use a
  /// variable named `<obj_id>_mru_base`.
  /// @pre the input variable \p v must be either a reference or an object id
  /// @note use this function when the base address is required for operation
  /// @warning do not use this function to check whether a base address is
  /// created
  variable_t get_or_insert_base_addr(const variable_t &v,
                                     bool is_cache = false) {
    assert(v.get_type().is_reference() || v.get_type().is_region());
    if (!is_cache) {
      assert(v.get_type().is_reference() || v.get_type().is_reference_region());
    }
    variable_t tmp = v;
    if (is_cache) {
      varname_t v_var_name = v.name();
      std::string str_mru_name = "_mru";
      varname_t mru_name = v_var_name.get_var_factory().get_or_insert_varname(
          v_var_name, str_mru_name);
      variable_t v_rgn_mru(mru_name, v.get_type());
      tmp = v_rgn_mru;
    }
    return get_or_insert_ghost(m_refs_base_addrs_map, tmp, "_base",
                               crab::INT_TYPE, 32);
  }

  /// @brief get the corresponding base address variable
  /// @param[in] v a reference variable or an object id
  /// @return if we find the ghost, returns; otherwise, return boost::none
  /// @note use this function for special handling when we do not find the ghost
  boost::optional<variable_t> get_base_addr(const variable_t &v) const {
    assert(v.get_type().is_reference() || v.get_type().is_region());
    return get_ghost(m_refs_base_addrs_map, v);
  }

  variable_t get_or_insert_write_region(const variable_t &rgn) {
    assert(rgn.get_type().is_region());
    return get_or_insert_ghost(m_ghost_rgn_map, rgn, "_write",
                               rgn.get_type().get_type_kind(),
                               rgn.get_type().get_bitwidth());
  }

  boost::optional<variable_t> get_write_region(const variable_t &rgn) const {
    assert(rgn.get_type().is_region());
    return get_ghost(m_ghost_rgn_map, rgn);
  }

  /// @brief check the reference refer the MRU object
  /// @param ref a reference variable
  /// @param id an object id
  /// @return true if ref_base == id_mru_base satisfies in address domain
  bool test_ref_refer_mru_object(const variable_t &ref, const obj_id_t &id) {
    auto ref_base = get_or_insert_base_addr(ref);
    auto id_mru_base = get_or_insert_base_addr(id, true);
    return m_addrs_dom.equals(ref_base, id_mru_base);
  }

  /// @brief record base_addr == \p addr in the address domain, dropping any
  /// equalities \p addr held before
  /// @param base_addr an address the domain already relates (or a new one)
  /// @param addr the address being (re)pointed at \p base_addr 's class
  /// @note the equality domain's add() unions the two classes, which keeps
  /// \p addr 's previous partners. That is wrong here: an address is
  /// re-pointed, not aliased with everything it used to equal -- e.g. when a
  /// cache miss moves an object's MRU base address onto another object,
  /// unioning would claim the two objects share a base address. So forget
  /// \p addr first, which leaves it alone in \p base_addr 's class.
  void assign_base_addr(const variable_t &base_addr, const variable_t &addr) {
    m_addrs_dom.operator-=(addr);
    m_addrs_dom.add(base_addr, addr);
  }

  /***************** ODI map operations *****************/
  /// @brief change singleton object to non-singleton
  /// @pre based on number of references in object_info, the abstract object now
  /// becomes a non-singleton object
  /// @post the abstract object referred by \p id is non-singleton
  /// @param[in] id an object id
  /// @par Side Effects:
  ///   This function modifies the abstract state which invokes this function.
  ///   Depending on whether we put properties for singleton object on the base
  ///   domain, changing singleton object involves different domain operations
  void change_object_status(const obj_id_t &id) {

    base_dom_variable_vector_t flds_vec;
    get_raw_obj_flds(id, flds_vec);

    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    if (!prod_ref) {
      CRAB_ERROR(domain_name(), "::", __func__, ": object ", id,
                 " value is not found on the odi map");
    }
    // NOTE: copy the following value is required because
    // updating the odi map is constructing a new tree based on previous tree
    odi_domain_product_t res_prod = *prod_ref;
    object_info_t &obj_info = odi_map_t::object_info_val(res_prod);
    obj_info.refcount_val().increment(id);

    // Two cases here:
    // Although there is no need to update the cache, the singleton's
    // properties must be copied from the cache into the summary: the
    // summary abstracts one or more concrete objects.
    // NOTE: the NEW allocation's (uninitialized) contents are deliberately
    // NOT added to the summary. Analyzed programs are assumed UB-free, so
    // memory is always initialized before it is read; region_domain makes
    // the same assumption (its ref_make only bumps the refcount).
    object_value_t &odi_val = odi_map_t::object_odi_val(res_prod);
    // No cache reset is needed (no ref_store/ref_load occurred), but any
    // pending reduction between the cache and the base must be applied
    // before the singleton is copied into the summary.
    apply_reduction_between_object_and_base(m_base_dom, odi_val, id);
    m_is_bottom = m_base_dom.is_bottom();
    // copy singleton object into summary
    odi_map_t::object_sum_val(odi_val) = odi_map_t::object_cache_val(odi_val);
    obj_info.sumpresence_val() = obj_info.cacheused_val();
    odi_map_t::object_eq_val(odi_val).set_to_top();

    m_odi_map.set(id, std::move(res_prod));

    CRAB_LOG(
        "object-change-status",
        crab::outs() << "singleton to non-singleton on object: " << id << "\n";
        m_odi_map.odi_write(crab::outs(), res_prod); crab::outs() << "\n";);
  }

  /***************** Cache operations *****************/
  /// @brief ensure the MRU cache holds \p id's object for an access through
  /// \p ref: on a miss (\p ref refers a different memory object from the
  /// cached one) delegate to the odi map to commit and reload the cache, then
  /// record the new MRU base-address binding in the address domain
  /// @param id an object id
  /// @param ref the reference variable
  /// @param rgn the region which reference referred
  /// @param reg_symb optional symbolic variable representing equality between
  ///   register and field
  /// @param offset_size_symb optional symbolic variable for ghost variable
  /// @param is_store \c true from transfer function on ref_store, \c false from
  /// ref_load
  /// @return \c true if the cache already held \p id's object (a hit)
  bool ensure_mru_cache(
      obj_id_t &id, const variable_t &ref, const variable_t &rgn,
      boost::optional<usymb_t> &reg_symb,
      boost::optional<std::pair<usymb_t, usymb_t>> &offset_size_symb,
      bool is_store) {

    // get the variable represented the mru base address
    variable_t mru_obj_base = get_or_insert_base_addr(id, true);
    CRAB_LOG(
        "object-entailment",
        crab::outs() << mru_obj_base << " == " << get_or_insert_base_addr(ref)
                     << "?\n"
                     << "Addrs = " << m_addrs_dom << "\n"
                     << "Is mru cached? "
                     << (test_ref_refer_mru_object(ref, id) ? "true" : "false")
                     << "\n";);
    ghost_variables_eq_t rgn_eq_gvars = get_or_insert_eq_gvars(rgn);
    if (is_store && reg_symb != boost::none) {
      rgn_eq_gvars = get_or_insert_eq_gvars(get_or_insert_write_region(rgn));
    }

    bool update_new_mru = m_odi_map.invalidate_cache_if_miss(
        id, *this, std::move(rgn_eq_gvars), reg_symb, offset_size_symb,
        test_ref_refer_mru_object(ref, id), is_store);

    // update address dom and object info if cache is missed
    if (update_new_mru) {
      assign_base_addr(get_or_insert_base_addr(ref), mru_obj_base);
      // remove pointer alias information for fields since cache is updated
      variable_vector_t obj_flds, obj_base_flds;
      get_obj_flds(id, obj_flds);
      for (const auto &v : obj_flds) {
        if (v.get_type().is_reference_region()) {
          obj_base_flds.push_back(get_or_insert_base_addr(v));
        }
      }
      m_addrs_dom.forget(obj_base_flds);
    }
    return !update_new_mru;
  }

  /// @brief a copy of \p dom restricted to the variables tracked by
  /// \p eq_dom (an equality domain); top if it tracks none
  /// @param[in] eq_dom the equality domain whose variables are kept
  template <typename Dom, typename EqDom>
  static Dom project_onto_eq_vars(const Dom &dom, const EqDom &eq_dom) {
    auto vars_opt = eq_dom.get_all_variables();
    if (!vars_opt) {
      return dom.make_top();
    }
    Dom out = dom;
    out.project(*vars_opt);
    return out;
  }

  boost::optional<object_domain_t> cow_apply_reduction() const {
    bool cow = false;
    for (auto kv : m_odi_map) {
      const std::shared_ptr<odi_domain_product_t> &prod = kv.second;
      const object_info_t prod_info = odi_map_t::object_info_val(*prod);
      if (prod_info.cache_reg_stored_val() ||
          prod_info.cache_reg_loaded_val()) {
        cow = true;
        break;
      }
    }
    if (cow) {
      object_domain_t copied = *this;
      copied.apply_reduction_based_on_flags(true);
      return copied;
    } else {
      return boost::none;
    }
  }

  /// @brief Applying domain reduction based on compilation flags
  /// @details a method to apply domain reduction before each transfer function
  /// by specifying an domain parameter:
  ///  object.reduction_level=FULL_REDUCTION or
  ///  object.reduction_level=REDUCTION_BEFORE_CHECK && reduce_for_assert = true
  /// @note domain reduction could be disabled if parameters:
  ///  object.reduction_level=NO_REDUCTION
  /// @param[in] reduce_for_assert a special bool indicate reduction is called
  /// before assertion
  void apply_reduction_based_on_flags(bool reduce_for_assert = false,
                                      bool reset_bool_flag = true) {
    OBJECT_DOMAIN_SCOPED_STATS(".overall_reduction");
    if (crab_domain_params_man::get().reduction_level() ==
            object_domain_params::reduction_level_t::NO_REDUCTION ||
        (crab_domain_params_man::get().reduction_level() ==
             object_domain_params::reduction_level_t::REDUCTION_BEFORE_CHECK &&
         !reduce_for_assert)) {
      // if the flags are set properly, no reduction is needed.
      return;
    }
    bool no_store = true;
    bool no_load = true;
    CRAB_LOG("object-reduction", crab::outs() << "State Before Reduction:\n"
                                              << *this << "\n";);
    // NOTE on `auto kv` over m_odi_map (here and below): the patricia
    // iterator yields a binding_t of two REFERENCES into the current leaf,
    // so `auto kv` copies two references, never the entry, and `const
    // auto &` would save nothing. The iterator itself keeps the visited
    // leaf and its spine alive (tree nodes are shared_ptr-owned), which is
    // also why set() on the map inside the body is safe: the loop keeps
    // walking the original tree snapshot.
    for (auto kv : m_odi_map) {
      // for each abstract object:
      //    there is a variable set of Symb, a set of flds and a set of regs
      //    a domain reduction from cache_dom to base_dom is:
      //  base_dom' = (\exists Symb (\exists flds. eq_fld_dom ^ cache_dom) ^
      //  eq_reg_dom) ^ base_dom
      const obj_id_t &id = kv.first;
      const std::shared_ptr<odi_domain_product_t> &prod = kv.second;
      const object_info_t prod_info = odi_map_t::object_info_val(*prod);
      const small_range &num_refs = prod_info.refcount_val();
      if (prod_info.cache_reg_stored_val()) {
        no_store = false;
      }
      if (prod_info.cache_reg_loaded_val()) {
        // The reduction should only be performed out when some register are
        // loaded values from this object
        no_load = false;
        reduce_object_to_base(id, prod_info, *prod);
      }
    }

    if (no_store && no_load)
      return;

    base_abstract_domain_t regs_only_base =
        project_onto_eq_vars(m_base_dom, m_eq_regs_dom);
    for (auto kv : m_odi_map) {
      // for each abstract object:
      //  a domain reduction from base_dom to cache_dom,
      //  cache_dom' = (\exists Symb (\exists regs. eq_reg_dom ^ base_dom) ^
      //  eq_fld_dom) ^ cache_dom
      const obj_id_t &id = kv.first;
      const std::shared_ptr<odi_domain_product_t> &prod = kv.second;
      const object_info_t prod_info = odi_map_t::object_info_val(*prod);
      const small_range &num_refs = prod_info.refcount_val();
      if (prod_info.cache_reg_stored_val() ||
          prod_info.cache_reg_loaded_val()) {
        // if some fields are stored from regs
        // perform reduction from base to obj
        reduce_base_to_object(id, prod_info, *prod, regs_only_base,
                              reset_bool_flag);
      }
    }
    CRAB_LOG("object-reduction", crab::outs() << "State After Reduction:\n"
                                              << *this << "\n";);
  }

  /// @brief propagate the cached properties of object \p id into the base
  /// domain and store the updated odi back into the map
  /// @param[in] id an object id
  /// @param[in] prod_info the object info of \p id (copied into the map)
  /// @param[in] prod the odi of \p id
  void reduce_object_to_base(const obj_id_t &id, const object_info_t &prod_info,
                             const odi_domain_product_t &prod) {
    object_value_t out_prod_val = odi_map_t::object_odi_val(prod);
    object_info_t out_prod_info = prod_info;
    apply_reduction_from_object_to_base(m_base_dom, out_prod_val, id);
    m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                           std::move(out_prod_val)));
  }

  /// @brief propagate register properties from \p regs_only_base into the
  /// cache of object \p id and store the updated odi back into the map
  /// @param[in] id an object id
  /// @param[in] prod_info the object info of \p id (copied into the map)
  /// @param[in] prod the odi of \p id
  /// @param[in] regs_only_base the base domain already projected onto the
  /// variables tracked by the register-equality domain
  /// @param[in] reset_flags if true, clear the reduction flags on the
  /// stored info
  void reduce_base_to_object(const obj_id_t &id, const object_info_t &prod_info,
                             const odi_domain_product_t &prod,
                             base_abstract_domain_t &regs_only_base,
                             bool reset_flags) {
    object_value_t out_prod_val = odi_map_t::object_odi_val(prod);
    object_info_t out_prod_info = prod_info;
    if (reset_flags) {
      out_prod_info.cache_reg_loaded_val() = false;
      out_prod_info.cache_reg_stored_val() = false;
    }
    apply_reduction_from_base_to_object(regs_only_base, out_prod_val, id);
    m_odi_map.set(id, odi_domain_product_t(std::move(out_prod_info),
                                           std::move(out_prod_val)));
  }

  /// @brief method for crab intrinsic to perform reduction for a specific
  /// object \p id
  /// @param[in] id an object id
  /// @param[in] is_from_cache_to_base boolean flag indicates reduction
  /// direction
  void apply_reduction_per_object(const obj_id_t &id,
                                  bool is_from_cache_to_base) {
    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    if (prod_ref) {
      const object_info_t prod_info = odi_map_t::object_info_val(*prod_ref);
      const small_range &num_refs = prod_info.refcount_val();
      if (is_from_cache_to_base) {
        if (prod_info.cache_reg_loaded_val()) {
          // if some registers are loaded values from this object
          // perform reduction from object to base
          reduce_object_to_base(id, prod_info, *prod_ref);
        }
      } else {
        if (prod_info.cache_reg_stored_val() ||
            prod_info.cache_reg_loaded_val()) {
          // if some fields are stored from regs
          // perform reduction from base to obj
          base_abstract_domain_t regs_only_base =
              project_onto_eq_vars(m_base_dom, m_eq_regs_dom);
          reduce_base_to_object(id, prod_info, *prod_ref, regs_only_base,
                                /*reset_flags=*/false);
        }
      }
    }
  }

  /// @brief transfer "regs == flds" constraints between MRU cache and base
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
        project_onto_eq_vars(base_dom, m_eq_regs_dom);
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
      boost::optional<usymb_t> t_opt = (*eq_flds_dom).get_class_id(v);
      if (!t_opt) {
        continue;
      }
      auto it = map.find(*t_opt);
      if (it != map.end()) {
        std::get<0>(it->second).push_back(v);
      } else {
        map.insert({*t_opt, {{v}, {}}});
      }
    }

    // for each term t, fill \p map by {t |-> <[flds], [regs]> }
    // any register variable assigned by t on \p eq_regs_dom
    for (auto it = map.begin(); it != map.end(); ++it) {
      const usymb_t &symb = it->first;
      auto regs_opt = eq_regs_dom.get_variables(symb);
      if (regs_opt) {
        std::get<1>(it->second)
            .insert(std::get<1>(it->second).end(), regs_opt->begin(),
                    regs_opt->end());
      }
    }
  }

  /// @brief Reduce equalities between registers and fields from cache to base
  /// @param[in,out] base_dom the base domain at the abstract state invokes this
  /// function
  /// @param[in] prod object value stores a product of <SUM dom, cache dom, EQ
  /// dom>
  /// @param id the object id refers \p prod in the \c m_odi_map
  void apply_reduction_from_object_to_base(base_abstract_domain_t &base_dom,
                                           const object_value_t &prod,
                                           const obj_id_t &id) const {

    OBJECT_DOMAIN_SCOPED_STATS(".reduction");
    if (crab_domain_params_man::get().reduction_level() ==
        object_domain_params::reduction_level_t::NO_REDUCTION) {
      // No reduction when object.reduction_level=NONE
      return;
    }
    // The following reduction is performed the followings:
    // a. get equalities constraints from 'eq_regs_dom' and 'eq_flds_dom'
    // b. rename a field from 'reduced_cache' to a register
    //    if there is a register with the same symbol as that field;
    //    if there are multiple registers equal to a field, add equality
    //    constraints into 'reduced_cache'.
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
    get_raw_obj_flds(id, obj_flds);
    const eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(prod);
    const cache_domain_t &cache_dom = odi_map_t::object_cache_val(prod);
    cache_domain_t reduced_cache =
        project_onto_eq_vars(cache_dom, *eq_flds_dom);

    symb_flds_regs_map_t map;
    build_map(map, m_eq_regs_dom, eq_flds_dom, obj_flds);
    CRAB_LOG("object-reduce", symbol_map_write(crab::outs(), map));

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t &flds = std::get<0>(it->second);
      const base_dom_variable_vector_t &regs = std::get<1>(it->second);
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
    if (!(base_dom <= *reduced_cache)) {
      base_dom &= *reduced_cache;
    }
    CRAB_LOG(
        "object-reduce", crab::outs() << "After Reduction:\n"
                                      << "base = " << base_dom << "\n"
                                      << "eq_regs = " << m_eq_regs_dom << "\n"
                                      << "odi = ";
        m_odi_map.odi_val_write(crab::outs(), prod); crab::outs() << "\n";);
  }

  /// @brief Reduce equalities between registers and fields from base to cache
  /// @param[in,out] base_dom the base domain at the abstract state invokes this
  /// function
  /// @param[in,out] prod object value stores a product of <SUM dom, cache dom,
  /// EQ dom>
  /// @param[in] id the object id refers \p prod in the \c m_odi_map
  void apply_reduction_from_base_to_object(base_abstract_domain_t &base_dom,
                                           object_value_t &prod,
                                           const obj_id_t &id) const {

    OBJECT_DOMAIN_SCOPED_STATS(".reduction");
    if (crab_domain_params_man::get().reduction_level() ==
        object_domain_params::reduction_level_t::NO_REDUCTION) {
      // No reduction when object.reduction_level=NONE
      return;
    }
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

    base_dom_variable_vector_t obj_flds, all_flds;
    std::unordered_map<base_dom_variable_t, base_dom_variable_t> rgn_wrt_map;
    get_raw_obj_flds(id, obj_flds);
    get_raw_obj_write_flds_map(id, rgn_wrt_map);
    all_flds = obj_flds;
    std::transform(rgn_wrt_map.begin(), rgn_wrt_map.end(),
                   std::back_inserter(all_flds),
                   [](const std::pair<const base_dom_variable_t,
                                       base_dom_variable_t> &pair) {
                     return pair.first;
                   });
    eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(prod);
    cache_domain_t &cache_dom = odi_map_t::object_cache_val(prod);

    symb_flds_regs_map_t map;
    build_map(map, m_eq_regs_dom, eq_flds_dom, all_flds);
    base_dom_variable_vector_t shadow_flds_used, flds_to_forget;
    CRAB_LOG("object-reduce", symbol_map_write(crab::outs(), map));

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t &flds = std::get<0>(it->second);
      const base_dom_variable_vector_t &regs = std::get<1>(it->second);
      if (flds.empty()) {
        // lost any equalities between flds and regs. Do not perform reduction
        continue;
      }
      assert(flds.size() > 0);
      auto it_flds = flds.begin();
      const base_dom_variable_t &var = *it_flds;
      if (regs.empty()) {
        // a pending store whose register died before this reduction: the
        // store itself still happened, so the field's OLD value must not
        // survive. Route the write ghosts through the rename epilogue with
        // no constraints attached: the field becomes top instead of
        // resurrecting its pre-store value.
        for (const auto &fld : flds) {
          auto it_w = rgn_wrt_map.find(fld);
          if (it_w != rgn_wrt_map.end()) {
            flds_to_forget.push_back(it_w->second);
            shadow_flds_used.push_back(it_w->first);
          }
        }
        continue;
      }
      auto it_map = rgn_wrt_map.find(var);
      if (it_map != rgn_wrt_map.end()) {
        flds_to_forget.push_back(it_map->second);
        shadow_flds_used.push_back(it_map->first);
      }
      assert(regs.size() > 0);
      const base_dom_variable_t &reg = regs[0];
      while (it_flds != flds.end()) {
        if (it_flds == flds.begin()) {
          base_dom += ikos::operator==(var, reg);
        } else {
          base_dom += ikos::operator==(var, *it_flds);
        }
        ++it_flds;
      }
    }
    // FIXME: avoid copying, use efficient project onto
    base_abstract_domain_t reduced_base = base_dom;
    reduced_base.project(all_flds);
    (*cache_dom).forget(flds_to_forget);
    *cache_dom &= reduced_base;
    // renaming: from rgn_w to rgn
    (*cache_dom).forget(flds_to_forget);
    (*eq_flds_dom).forget(flds_to_forget);
    (*cache_dom).rename(shadow_flds_used, flds_to_forget);

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
      for (const auto &kv : (*m_flds_id_map)) {
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
    // not using api from separate domain
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

  /**------------------ End helper method definitions ------------------**/

public:
  void dump() const { write(crab::outs()); }
  /**------------------ Begin domain API definitions ------------------**/
  object_domain_t make_top() const override {
    object_domain_t top = object_domain_t(true);
    top.m_flds_id_map = m_flds_id_map;
    top.m_refs_base_addrs_map = m_refs_base_addrs_map;
    top.m_ghost_rgn_map = m_ghost_rgn_map;
    return top;
  }

  object_domain_t make_bottom() const override {
    return object_domain_t(false);
  }

  void set_to_top() override {
    object_domain_t abs(true);
    abs.m_flds_id_map = m_flds_id_map;
    abs.m_refs_base_addrs_map = m_refs_base_addrs_map;
    abs.m_ghost_rgn_map = m_ghost_rgn_map;
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    object_domain_t abs(false);
    std::swap(*this, abs);
  }

  object_domain(bool is_top = true)
      : m_is_bottom(!is_top), m_ghost_var_num_man(var_type_fn()),
        m_ghost_var_eq_man(var_type_fn()) {}

  object_domain(const object_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(o.m_base_dom),
        m_odi_map(o.m_odi_map), m_addrs_dom(o.m_addrs_dom),
        m_eq_regs_dom(o.m_eq_regs_dom), m_flds_id_map(o.m_flds_id_map),
        m_refs_base_addrs_map(o.m_refs_base_addrs_map),
        m_ghost_rgn_map(o.m_ghost_rgn_map),
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
        m_ghost_rgn_map(std::move(o.m_ghost_rgn_map)),
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
      m_ghost_rgn_map = o.m_ghost_rgn_map;
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
      m_ghost_rgn_map = std::move(o.m_ghost_rgn_map);
      m_ghost_var_num_man = std::move(o.m_ghost_var_num_man);
      m_ghost_var_eq_man = std::move(o.m_ghost_var_eq_man);
    };
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    OBJECT_DOMAIN_SCOPED_STATS(".is_top");

    // m_addrs_dom and m_eq_regs_dom are deliberately not consulted: every
    // transfer function that records an address equality also constrains the
    // base domain (ref_gep/ref_assume write the address ghosts into it), and
    // every register symbol is created behind a successful odi-map lookup.
    // So "base and odi map are top" implies the two equality domains are too.
    bool res = (!is_bottom() && m_base_dom.is_top() && m_odi_map.is_top());
    if (::crab::CrabSanityCheckFlag && res) {
      if (!m_addrs_dom.is_top() || !m_eq_regs_dom.is_top()) {
        CRAB_ERROR(domain_name(), "::is_top: base and odi map are top but an "
                                  "equality domain is not; the invariant "
                                  "documented here no longer holds");
      }
    }
    return res;
  }

  bool operator<=(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".leq");

    CRAB_LOG("object-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                        << "\n<=\n\t" << o << "\n";);
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
      // set_to_top() keeps this side's shared maps; merge the other
      // side's entries too (see operator|)
      set_to_top();
      m_flds_id_map = merge_shared_map(m_flds_id_map, o.m_flds_id_map);
      m_refs_base_addrs_map =
          merge_shared_map(m_refs_base_addrs_map, o.m_refs_base_addrs_map);
      m_ghost_rgn_map = merge_shared_map(m_ghost_rgn_map, o.m_ghost_rgn_map);
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
      // The shared side maps (field->id, base addresses, write regions) are
      // incremental facts about the program, not part of the lattice value:
      // make_top() keeps this side's and BOTH sides' entries are merged so
      // that all states converge on one shared map (operator|= preserves
      // them the same way).
      object_domain_t abs = make_top();
      abs.m_flds_id_map = merge_shared_map(abs.m_flds_id_map, o.m_flds_id_map);
      abs.m_refs_base_addrs_map =
          merge_shared_map(abs.m_refs_base_addrs_map, o.m_refs_base_addrs_map);
      abs.m_ghost_rgn_map =
          merge_shared_map(abs.m_ghost_rgn_map, o.m_ghost_rgn_map);
      return abs;
    }

    CRAB_LOG("object-join",
             crab::outs() << "Join " << *this << "\n and " << o << "\n =\n");
    object_domain_t res(combine(*this, o, combine_kind::JOIN));
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

    object_domain_t res(combine(*this, o, combine_kind::MEET));

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
    *this = combine(*this, o, combine_kind::MEET);
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

    CRAB_LOG("object-widen",
             crab::outs() << "Widening " << *this << " and " << o << " =\n");

    object_domain_t res(combine(*this, o, combine_kind::WIDENING));

    CRAB_LOG("object-widen", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t
  widening_thresholds(const object_domain_t &o,
                      const thresholds<number_t> &thresholds) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".widening.thresholds");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("object-widen", crab::outs() << "Widening with threshold " << *this
                                          << " and " << o << " =\n");

    object_domain_t res(combine(*this, o, combine_kind::WIDENING, &thresholds));

    CRAB_LOG("object-widen", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t operator&&(const object_domain_t &o) const override {
    OBJECT_DOMAIN_SCOPED_STATS(".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("object-narrow",
             crab::outs() << "Narrowing " << *this << " and " << o << " =\n");

    object_domain_t res(combine(*this, o, combine_kind::NARROWING));

    CRAB_LOG("object-narrow", crab::outs() << res << "\n");
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

    apply_reduction_based_on_flags();

    if (get_obj_id(rgn) == boost::none) {
      // if a region does not belong to an object, treat it as an object
      // i.e. treat region as a field, as well as an object id
      obj_id_t id = create_new_obj_id(rgn);
      update_fields_id_map(rgn, id);
    }
    obj_id_t id = get_obj_id_or_fail(rgn);
    // fresh object: no references, not initialized, no summary, unused and
    // clean cache, no pending reg<->fld reduction
    object_info_t obj_info = object_info_t(
        small_range::zero(), /*obj_init=*/boolean_value::get_false(),
        /*sum_presence=*/boolean_value::get_false(),
        /*cache_used=*/boolean_value::get_false(),
        /*cache_dirty=*/boolean_value::get_false(), /*is_loaded=*/false,
        /*is_stored=*/false);
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
  /// @param[in] size a symbolic / constant val represented allocation size
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
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object goes to top
        return;
      }
      const object_info_t obj_prod_info = odi_map_t::object_info_val(*prod_ref);
      const small_range &num_refs = obj_prod_info.refcount_val();
      if (num_refs.is_zero()) {
        // First allocation: the object becomes a singleton (refcount 0 -> 1).
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

    // the loaded content is unknown: havoc the destination
    if (is_unknown_region(rgn)) {
      get_or_insert_gvars(res).forget(m_base_dom);
      m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
      return;
    }

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
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object is top: the load result is unknown
        res_gvars.forget(m_base_dom);
        m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
        return;
      }
      const object_info_t obj_info_ref = odi_map_t::object_info_val(*prod_ref);
      const bool is_loaded = obj_info_ref.cache_reg_loaded_val();
      const small_range &num_refs = obj_info_ref.refcount_val();

      // In crab IR, the number of references cannot be zero
      //  if ref_load access a not null reference.
      // So zero case should not exist
      assert(!num_refs.is_zero());
      if (num_refs.is_zero()) { // FIXME: handle the cast where region cast from
                                // an unknown region to a valid region since it
                                // requires to support unknown region, ignore
                                // load / store for now
        res_gvars.forget(m_base_dom);
        m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
        return;
      }

      // use odi map
      // keep the equality between res == rgn as memory load
      // a. generate or reuse a symbolic variable for rgn
      // b. keep that symbolic variable assigned to res
      boost::optional<usymb_t> reg_symb = boost::none;
      boost::optional<std::pair<usymb_t, usymb_t>> reg_offset_size_symb =
          boost::none;
      bool is_hit = ensure_mru_cache((*id_opt), ref, rgn, reg_symb,
                                     reg_offset_size_symb, false);
      // forget old property for register res
      res_gvars.forget(m_base_dom);
      // Optimize: directly assignment if property is a constant
      prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) {
        // the pre-access reduction inside the cache update discovered
        // infeasibility: the odi map became bottom, and so is the state
        set_to_bottom();
        return;
      }
      unsigned fully_assigned = 0;
      boost::optional<ghost_variables_eq_t> res_eq_gvars = get_eq_gvars(res);
      // own the value (three COW pointer copies, no domain data copied):
      // references into the leaf die when set() replaces it -- see the
      // find() invalidation contract; ref_store already owns its copy
      const object_value_t out_prod = odi_map_t::object_odi_val(*prod_ref);
      const cache_domain_t &flds_dom = odi_map_t::object_cache_val(out_prod);
      const eq_fields_domain_t &eq_flds_dom =
          odi_map_t::object_eq_val(out_prod);
      if (auto rgn_const = flds_dom.at(rgn_gvars.get_var()).singleton()) {
        m_base_dom.assign(res_gvars.get_var(), *rgn_const);
        fully_assigned++;
      }
      // assigning register with symbolic variable
      m_eq_regs_dom.set(res_eq_gvars.value().get_var(), *reg_symb);
      if (res_eq_gvars.value().has_offset_and_size()) {
        auto offset_base_var = rgn_gvars.get_offset_and_size().get_offset();
        auto size_base_var = rgn_gvars.get_offset_and_size().get_size();
        if (auto offset_const = flds_dom.at(offset_base_var).singleton()) {
          m_base_dom.assign(res_gvars.get_offset_and_size().get_offset(),
                            *offset_const);
          fully_assigned++;
        }
        m_eq_regs_dom.set(
            res_eq_gvars.value().get_offset_and_size().get_offset(),
            std::get<0>(*reg_offset_size_symb));
        if (auto sz_const = flds_dom.at(size_base_var).singleton()) {
          m_base_dom.assign(res_gvars.get_offset_and_size().get_size(),
                            *sz_const);
          fully_assigned++;
        }
        m_eq_regs_dom.set(res_eq_gvars.value().get_offset_and_size().get_size(),
                          std::get<1>(*reg_offset_size_symb));
      }
      bool is_res_ref = res_eq_gvars.value().has_offset_and_size();
      bool update_load = (!is_res_ref && fully_assigned == 1) ||
                         (is_res_ref && fully_assigned == 3);
      auto rgn_w = get_or_insert_write_region(rgn);
      auto rgnw_eq_gvars = get_or_insert_eq_gvars(rgn_w);
      bool forget_w = false;
      // equalties for rgn_w should be dropped since if load_ref follows
      // an assume. rgnw == reg will no longer equals.
      // NOTE: offset/size ghosts exist only under
      // region.is_dereferenceable; guard by their presence, not by the
      // region's type
      if (rgnw_eq_gvars.has_offset_and_size()) {
        forget_w |=
            (*eq_flds_dom)
                .find_opt(rgnw_eq_gvars.get_offset_and_size().get_offset()) !=
            boost::none;
        forget_w |=
            (*eq_flds_dom)
                .find_opt(rgnw_eq_gvars.get_offset_and_size().get_size()) !=
            boost::none;
      }
      forget_w |=
          (*eq_flds_dom).find_opt(rgnw_eq_gvars.get_var()) != boost::none;
      if (update_load || forget_w) {
        object_info_t out_info = odi_map_t::object_info_val(*prod_ref);
        if (update_load) {
          out_info.cache_reg_loaded_val() = is_hit ? is_loaded : false;
        }
        if (forget_w) {
          object_value_t out_prod = odi_map_t::object_odi_val(*prod_ref);
          eq_fields_domain_t &eq_flds_dom = odi_map_t::object_eq_val(out_prod);
          m_ghost_var_eq_man.forget(rgn_w, *eq_flds_dom);
          m_odi_map.set(*id_opt, odi_domain_product_t(std::move(out_info),
                                                      std::move(out_prod)));
        } else {
          // update odi map; out_prod is a const view into the map entry,
          // so the copy is required (set() replaces the entry)
          m_odi_map.set(*id_opt,
                        odi_domain_product_t(std::move(out_info),
                                             object_value_t(out_prod)));
        }
      }
      // the loaded reference points wherever the field points: re-point
      // res's base address into rgn's class (assign_base_addr forgets its
      // SECOND argument, so the field's binding must be the kept side)
      if (res.get_type().is_reference()) {
        assign_base_addr(get_or_insert_base_addr(rgn),
                         get_or_insert_base_addr(res));
      }

    } else {
      // the region is not tracked by any abstract object: the load result
      // is unknown
      res_gvars.forget(m_base_dom);
      m_ghost_var_eq_man.forget(res, m_eq_regs_dom);
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

    apply_reduction_based_on_flags();

    // use ghost variable for field
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);

    // The reference ref should not be evaluated as a null pointer.
    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_store: reference ",
                                   ref, " is null."););
      crab::CrabStats::count(domain_name() +
                             ".count.ref_store.skipped.null_reference");
      operator-=(rgn);
      return;
    }

    if (auto id_opt = get_obj_id(rgn)) {

      // retrieve an abstract object odi
      const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
      if (!prod_ref) { // object does not exist, default is top
        return;
      }
      const object_info_t obj_info_ref = odi_map_t::object_info_val(*prod_ref);
      const bool is_stored = obj_info_ref.cache_reg_stored_val();
      const small_range num_refs = obj_info_ref.refcount_val();
      assert(!num_refs.is_zero());
      if (num_refs.is_zero()) { // FIXME: handle the cast where region cast from
                                // an unknown region to a valid region since it
                                // requires to support unknown region, ignore
                                // load / store for now
        return;
      }

      // use odi map
      // keep the equality between rgn == val as a memory store
      // a. generate or reuse a symbolic variable for val if val is variable
      // b. keep that symbolic variable assigned to rgn
      boost::optional<usymb_t> reg_symb = boost::none;
      boost::optional<std::pair<usymb_t, usymb_t>> reg_offset_size_symb =
          boost::none;
      if (val.is_variable()) { // obtain symbols of the val
        ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
        reg_symb = get_symbol_or_fresh(m_eq_regs_dom, val_gvars.get_var());
        if (val_gvars.has_offset_and_size()) {
          reg_offset_size_symb = std::make_pair(
              get_symbol_or_fresh(m_eq_regs_dom,
                                  val_gvars.get_offset_and_size().get_offset()),
              get_symbol_or_fresh(m_eq_regs_dom,
                                  val_gvars.get_offset_and_size().get_size()));
        }
      }
      bool is_hit = ensure_mru_cache((*id_opt), ref, rgn, reg_symb,
                                     reg_offset_size_symb, true);
      // retrieve an abstract object odi again since cache may be flushed
      prod_ref = m_odi_map.find((*id_opt));
      if (!prod_ref) {
        // the pre-access reduction inside the cache update discovered
        // infeasibility: the odi map became bottom, and so is the state
        set_to_bottom();
        return;
      }
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
        // forget the write region rgn_w in the equality field domain
        variable_t rgn_w = get_or_insert_write_region(rgn);
        m_ghost_var_eq_man.forget(rgn_w, *eq_flds_dom);
        out_info.cache_reg_stored_val() = is_hit ? is_stored : false;
        // forget the region rgn in the equality field domain
        m_ghost_var_eq_man.forget(rgn, *eq_flds_dom);
      } else { // val is a variable (i.e. register)
        ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
        boost::optional<ghost_variables_eq_t> val_eq_gvars =
            get_eq_gvars(val.get_variable());
        // Optimize: assign variable with constant value directly
        unsigned fully_assigned = 0;
        auto rgn_eq_gvars = get_or_insert_eq_gvars(rgn);
        auto rgnw_eq_gvars =
            get_or_insert_eq_gvars(get_or_insert_write_region(rgn));
        if (auto val_const = m_base_dom.at(val_gvars.get_var()).singleton()) {
          flds_dom.assign(rgn_gvars.get_var(), *val_const);
          (*eq_flds_dom) -= rgnw_eq_gvars.get_var();
          (*eq_flds_dom) -= rgn_eq_gvars.get_var();
          fully_assigned++;
        } else {
          if (crab_domain_params_man::get().reduction_level() ==
              object_domain_params::reduction_level_t::NO_REDUCTION) {
            // if there is no reduction, we need to make sure soundness
            // without providing rgn_w == val
            (*flds_dom) -= rgn_eq_gvars.get_var();
            (*eq_flds_dom) -= rgnw_eq_gvars.get_var();
          } else {
            // assigning register with symbolic variable
            m_eq_regs_dom.set(val_eq_gvars.value().get_var(), *reg_symb);
          }
        }
        if (val_eq_gvars.value().has_offset_and_size()) {
          auto offset_base_var = val_gvars.get_offset_and_size().get_offset();
          auto size_base_var = val_gvars.get_offset_and_size().get_size();
          if (auto offset_const = m_base_dom.at(offset_base_var).singleton()) {
            flds_dom.assign(rgn_gvars.get_offset_and_size().get_offset(),
                            *offset_const);
            (*eq_flds_dom) -= rgnw_eq_gvars.get_offset_and_size().get_offset();
            (*eq_flds_dom) -= rgn_eq_gvars.get_offset_and_size().get_offset();
            fully_assigned++;
          } else {
            if (crab_domain_params_man::get().reduction_level() ==
                object_domain_params::reduction_level_t::NO_REDUCTION) {
              (*flds_dom) -= rgn_gvars.get_offset_and_size().get_offset();
              (*eq_flds_dom) -=
                  rgnw_eq_gvars.get_offset_and_size().get_offset();
            } else {
              // assigning register with symbolic variable
              m_eq_regs_dom.set(
                  val_eq_gvars.value().get_offset_and_size().get_offset(),
                  std::get<0>(*reg_offset_size_symb));
            }
          }
          if (auto size_const = m_base_dom.at(size_base_var).singleton()) {
            flds_dom.assign(rgn_gvars.get_offset_and_size().get_size(),
                            *size_const);
            (*eq_flds_dom) -= rgnw_eq_gvars.get_offset_and_size().get_size();
            (*eq_flds_dom) -= rgn_eq_gvars.get_offset_and_size().get_size();
            fully_assigned++;
          } else {
            if (crab_domain_params_man::get().reduction_level() ==
                object_domain_params::reduction_level_t::NO_REDUCTION) {
              (*flds_dom) -= rgn_gvars.get_offset_and_size().get_size();
              (*eq_flds_dom) -= rgnw_eq_gvars.get_offset_and_size().get_size();
            } else {
              // assigning register with symbolic variable
              m_eq_regs_dom.set(
                  val_eq_gvars.value().get_offset_and_size().get_size(),
                  std::get<1>(*reg_offset_size_symb));
            }
          }
        }
        // the field's old cache value is kept until the pending store
        // (rgn_w == #symb) is committed by the next reduction; if the
        // symbol's register dies first, the reduction's dead-symbol case
        // drops the field instead (see apply_reduction_from_base_to_object)
        bool is_val_ref = val_eq_gvars.value().has_offset_and_size();
        if ((!is_val_ref && fully_assigned == 1) ||
            (is_val_ref && fully_assigned == 3)) {
          out_info.cache_reg_stored_val() = is_hit ? is_stored : false;
        }
      }

      // if rgn is reference, we add rgn_base == val
      if (val.is_variable() && val.get_type().is_reference()) {
        assign_base_addr(get_or_insert_base_addr(val.get_variable()),
                         get_or_insert_base_addr(rgn));
      }

      // update object info
      object_info_t out_obj_info = object_info_t(
          num_refs, boolean_value::get_true() /*Object is inited*/,
          out_info.sumpresence_val(),
          boolean_value::get_true() /*Cache is used*/,
          boolean_value::get_true() /*Cache is dirty*/,
          out_info.cache_reg_loaded_val(), out_info.cache_reg_stored_val());
      m_odi_map.set(*id_opt, odi_domain_product_t(std::move(out_obj_info),
                                                  std::move(out_prod)));
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

    auto eval = [this](const linear_expression_t &e) {
      interval_t r = e.constant();
      for (const auto &p : e) {
        ERROR_IF_NOT_INT(p.second, __LINE__);
        r += p.first *
             m_base_dom.operator[](get_or_insert_gvars(p.second).get_var());
      }
      return r;
    };

    if (is_bottom()) {
      return;
    }

    apply_reduction_based_on_flags();

    ghost_variables_t ref1_gvars = get_or_insert_gvars(ref1);
    ghost_variables_t ref2_gvars = get_or_insert_gvars(ref2);

    if (ref1_gvars.has_offset_and_size() && ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().assign(
          m_base_dom, ref1_gvars.get_offset_and_size(), offset);
    } else if (ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().forget(m_base_dom);
    }

    if (rgn1 == rgn2 && (eval(offset) != (number_t(0)))) {
      // ref2 steps to another array element inside the same region:
      //
      //   p := make_ref(V, 4);      // array of ints, V holds ALL elements
      //   store_ref(V, p, 1);       // p[0] := 1
      //   q := gep_ref(V, p + 4);   // q = &p[1], same region V
      //   store_ref(V, q, 2);       // p[1] := 2 => V must be [1,2]
      //   x := load_ref(V, p);      // must read [1,2], not exactly 2
      //
      // FIXME: this branch is a workaround. Ideally p and q share the same
      // base address (same allocation), but the cache does STRONG updates
      // for the reference owning it: with p_base == q_base, the store
      // through q would set the cache to V = 2 and the load through p
      // would wrongly read exactly 2. So we promote the object
      // (change_object_status below) and give q NO base equality, which
      // forces p's accesses to the committed summary ([1,2]). The clean
      // fix: keep the base equality and let stores to array-abstracted
      // regions update the summary weakly, bypassing the cache.
      // "NO base equality" must also hold for a REUSED variable: drop any
      // stale binding ref2 carries (e.g. from a previous gep of the same
      // crab variable), or its store would still strong-update the cache
      m_addrs_dom -= get_or_insert_base_addr(ref2);
      if (auto id_opt = get_obj_id(rgn1)) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
        if (!prod_ref) { // object goes to top
          return;
        }
        const object_info_t obj_prod_info =
            odi_map_t::object_info_val(*prod_ref);
        const small_range &num_refs = obj_prod_info.refcount_val();
        if (num_refs.is_one()) {
          change_object_status(*id_opt);
        }
      }
    } else {
      // assign equality: ref2 == ref1
      // In C memory model,
      // pointer arithmetic cannot be performed from different memory objects.
      // Thus, a precondition is ref2 and ref1 are belongs to a same memory obj.
      assign_base_addr(get_or_insert_base_addr(ref1),
                       get_or_insert_base_addr(ref2));
    }

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

      apply_reduction_based_on_flags(true, false);

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
        auto size_lin_csts = m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
            ref_cst, ghost_variable_kind::SIZE);
        m_base_dom += size_lin_csts;
        m_is_bottom = m_base_dom.is_bottom();
      }
      if (!m_is_bottom && ref_cst.is_equality()) {
        if (ref_cst.is_binary()) { // ref_cst is lhs == rhs
          const variable_t &lhs = ref_cst.lhs();
          const variable_t &rhs = ref_cst.rhs();
          assign_base_addr(get_or_insert_base_addr(rhs),
                           get_or_insert_base_addr(lhs));
        }
      }
      apply_reduction_based_on_flags(true);
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
      apply_reduction_based_on_flags();

      // We represent reference as numerical in domain
      m_base_dom.assign(get_or_insert_gvars(int_var).get_var(),
                        get_or_insert_gvars(ref_var).get_var());
      m_addrs_dom -= get_or_insert_base_addr(ref_var);
      // int_var is overwritten: any symbol linking it to an object field is
      // stale now (otherwise a later reduction would re-impose the field's
      // value on an address-valued integer)
      m_ghost_var_eq_man.forget(int_var, m_eq_regs_dom);
    }
  }

  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    OBJECT_DOMAIN_SCOPED_STATS(".int_to_ref");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      apply_reduction_based_on_flags();

      ghost_variables_t ref_gvars = get_or_insert_gvars(ref_var);
      ghost_variables_t int_gvars = get_or_insert_gvars(int_var);
      if (ref_gvars.has_offset_and_size()) {
        ref_gvars.get_offset_and_size().forget(m_base_dom);
      }

      ref_gvars.assign(m_base_dom, int_gvars);
      // ref_var is re-targeted to an arbitrary address: its old base-address
      // equalities (e.g. ownership of an object's MRU cache) no longer hold
      m_addrs_dom -= get_or_insert_base_addr(ref_var);
      m_ghost_var_eq_man.forget(ref_var, m_eq_regs_dom);
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
      // read-only probes go through the const view: the mutable raw
      // accessors detach the COW reference, deep-copying a shared value
      // just to answer a query
      const object_value_t &c_prod_value_ref = prod_value_ref;

      const small_range &num_refs = obj_info_ref.refcount_val();

      // In crab IR, the number of references cannot be zero
      // So zero case should not exist
      assert(!num_refs.is_zero());

      // num_refs > 1, non-singleton, use odi map
      // copy field
      if (!odi_map_t::object_sum_raw_val(c_prod_value_ref).is_top()) {
        base_lhs.forget(odi_map_t::object_sum_raw_val(prod_value_ref));
        base_lhs.assign(odi_map_t::object_sum_raw_val(prod_value_ref),
                        base_rhs);
        base_rhs.forget(odi_map_t::object_sum_raw_val(prod_value_ref));
      }
      if (!odi_map_t::object_cache_raw_val(c_prod_value_ref).is_top()) {
        base_lhs.forget(odi_map_t::object_cache_raw_val(prod_value_ref));
        base_lhs.assign(odi_map_t::object_cache_raw_val(prod_value_ref),
                        base_rhs);
        base_rhs.forget(odi_map_t::object_cache_raw_val(prod_value_ref));
      }
      boost::optional<usymb_t> t_opt =
          odi_map_t::object_eq_raw_val(c_prod_value_ref).get_class_id(rhs_rgn);
      if (t_opt) {
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
      // the copied value still carries rhs's OTHER fields: they are not
      // fields of lhs's object, so nothing can ever read or forget them
      // through lhs -- drop them before installing, or they remain as
      // unreachable dimensions in every later join/inclusion on lhs
      variable_vector_t rhs_flds;
      get_obj_flds(*rhs_id_opt, rhs_flds);
      for (auto &fld : rhs_flds) {
        if (fld == rhs_rgn || is_unknown_region(fld)) {
          continue;
        }
        m_ghost_var_num_man.forget(
            fld, odi_map_t::object_sum_raw_val(prod_value_ref));
        m_ghost_var_num_man.forget(
            fld, odi_map_t::object_cache_raw_val(prod_value_ref));
        m_ghost_var_eq_man.forget(fld,
                                  odi_map_t::object_eq_raw_val(prod_value_ref));
      }

      m_odi_map.set(lhs_id, odi_domain_product_t(std::move(obj_info_ref),
                                                 std::move(prod_value_ref)));
    }

    CRAB_LOG("object-region-copy",
             crab::outs() << "After region_copy(" << lhs_rgn << ":"
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
      // fresh object info: same initial state as region_init
      m_odi_map.set(
          id, odi_domain_product_t(
                  object_info_t(small_range::zero(),
                                /*obj_init=*/boolean_value::get_false(),
                                /*sum_presence=*/boolean_value::get_false(),
                                /*cache_used=*/boolean_value::get_false(),
                                /*cache_dirty=*/boolean_value::get_false(),
                                /*is_loaded=*/false, /*is_stored=*/false),
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
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    OBJECT_DOMAIN_SCOPED_STATS(".select_ref");

    auto compute_lhs_ref = [&lhs_ref, &lhs_rgn](
                               const variable_or_constant_t &ref,
                               const boost::optional<base_dom_variable_t> &rgn,
                               object_domain_t &out) {
      if (ref.is_reference_null()) {
        out -= lhs_ref;
        out.ref_assume(reference_constraint_t::mk_null(lhs_ref));
      } else {
        assert(ref.is_variable());
        assert(rgn);
        linear_expression_t zero_offset(number_t(0));
        out.ref_gep(ref.get_variable(), *rgn, lhs_ref, lhs_rgn, zero_offset);
      }
    };

    if (!is_bottom()) {
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
      // for address domain, we check boolean value of cond.
      // if it is determined, we insert alias information on m_addrs_dom.
      // NOTE: the following operations assume m_base_dom is using
      // flat_bool_domain
      boolean_value cond_val =
          object_domain_impl::base_is_instance_of_Flat_Boolean<
              base_abstract_domain_t>::get_bool_val_by_var(m_base_dom,
                                                           cond_gvars
                                                               .get_var());
      if (cond_val.is_true()) {
        compute_lhs_ref(ref1, rgn1, *this);
      } else if (cond_val.is_false()) {
        compute_lhs_ref(ref2, rgn2, *this);
      } else {
        // unknown condition: apply each branch on its own copy and join the
        // FULL states. Joining only some subdomains would keep branch 1's
        // odi map and register equalities as-is, which is not an upper bound
        // of branch 2 (compute_lhs_ref -> ref_gep can update the odi map).
        object_domain_t inv1(*this);
        compute_lhs_ref(ref1, rgn1, *this);
        compute_lhs_ref(ref2, rgn2, inv1);
        *this |= inv1;
      }
    }
  }

  /**************************** Numerical operations *************************/
  // x := y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    OBJECT_DOMAIN_SCOPED_STATS(".apply");

    if (!is_bottom()) {
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
      apply_reduction_based_on_flags();

      auto b_e = m_ghost_var_num_man.rename_linear_expr(e);
      m_base_dom.assign(get_or_insert_gvars(x).get_var(), b_e);
      m_eq_regs_dom -= x;
    }
  }

  // join(*this, copy_of_this(x := e)); delegated to the base domain, whose
  // weak_assign has exactly that semantics on the ghost variables. x may no
  // longer equal a field afterwards, so its register equality is dropped.
  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    OBJECT_DOMAIN_SCOPED_STATS(".weak_assign");

    if (!is_bottom()) {
      apply_reduction_based_on_flags();

      auto b_e = m_ghost_var_num_man.rename_linear_expr(e);
      m_base_dom.weak_assign(get_or_insert_gvars(x).get_var(), b_e);
      m_eq_regs_dom -= x;
    }
  }

  // Note: object_domain uses a fixed naming scheme for its ghost variables
  // for the base domain (i.e. variable_t and the base domain's variable_t
  // are the same type and no renaming is required), so rhs can be checked
  // against m_base_dom directly. The cache may still hold facts not yet
  // reduced into m_base_dom, so reduce on a copy first.
  bool entails(const linear_constraint_t &rhs) const override {
    if (is_bottom()) {
      return true;
    }
    auto reduced_opt = cow_apply_reduction();
    const object_domain_t &reduced = boost::get_optional_value_or(reduced_opt, *this);
    return reduced.m_base_dom.entails(rhs);
  }

  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {
    OBJECT_DOMAIN_SCOPED_STATS(".weak_assign_bool_cst");

    if (!is_bottom()) {
      apply_reduction_based_on_flags();

      auto b_rhs = m_ghost_var_num_man.rename_linear_cst(rhs);
      m_base_dom.weak_assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                      b_rhs);
      m_eq_regs_dom -= lhs;
    }
  }

  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {
    OBJECT_DOMAIN_SCOPED_STATS(".weak_assign_bool_var");

    if (!is_bottom()) {
      apply_reduction_based_on_flags();

      m_base_dom.weak_assign_bool_var(get_or_insert_gvars(lhs).get_var(),
                                      get_or_insert_gvars(rhs).get_var(),
                                      is_not_rhs);
      m_eq_regs_dom -= lhs;
    }
  }

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    OBJECT_DOMAIN_SCOPED_STATS(".add_constraints");

    if (!is_bottom()) {
      apply_reduction_based_on_flags(true, false);

      auto b_csts = m_ghost_var_num_man.rename_linear_cst_sys(csts);
      m_base_dom += b_csts;
      m_is_bottom = m_base_dom.is_bottom();
      apply_reduction_based_on_flags(true);
    }
  }

  /********************** Boolean operations **********************/

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    OBJECT_DOMAIN_SCOPED_STATS(".assign_bool_cst");

    if (!is_bottom()) {
      if (!rhs.is_tautology() && !rhs.is_contradiction()) {
        // Reduction required for the following code pattern
        // lhs := a < b;
        // assert(lhs);
        // The assertion requires to obtain the value of lhs for assertion check
        apply_reduction_based_on_flags(true);
      } else {
        apply_reduction_based_on_flags();
      }

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
      apply_reduction_based_on_flags(true, false);

      m_base_dom.assume_bool(get_or_insert_gvars(v).get_var(), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
      m_eq_regs_dom -= v;
      apply_reduction_based_on_flags(true);
    }
  }

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {

    OBJECT_DOMAIN_SCOPED_STATS(".select_bool");

    if (!is_bottom()) {
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
      apply_reduction_based_on_flags();

      auto rhs_lin_cst = m_ghost_var_num_man.ghosting_ref_cst_to_linear_cst(
          rhs, ghost_variable_kind::ADDRESS);
      m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(),
                                 rhs_lin_cst);
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

  // Allocation sites and tags are not tracked: decline (return false) so
  // downstream checkers report "unknown" instead of aborting the analysis.
  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
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

    apply_reduction_based_on_flags();

    forget_var(v);
  }

private:
  /// @brief forget one variable; shared core of \c operator-= and \c forget
  /// (which is the fold of this routine over its argument list, following
  /// region_domain's design). Pre: not bottom/top, reduction applied.
  ///
  /// To forget a variable, it is easier to forget if variable is not a rgn
  /// However, if it is a region, we need to perform the followings:
  /// 1. the abstract object that region var belongs to is a singleton,
  ///    forget it in base dom if we keep singletons in the base domain.
  /// 2. object is not singleton, forget it in odi map
  void forget_var(const variable_t &v) {
    if (v.get_type().is_region()) {
      // for now skip analysis for unknown region
      if (is_unknown_region(v)) {
        return;
      }
      if (auto id_opt = get_obj_id(v)) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
        if (!prod_ref) {
          return;
        }
        object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
        object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);

        const small_range &num_refs = prod_info_ref.refcount_val();
        variable_vector_t obj_flds;
        get_obj_flds(*id_opt, obj_flds);

        if (obj_flds.size() ==
            1) { // The region is part of an object that has only one region
          m_odi_map -= *id_opt;
        } else {
          // use odi map
          m_ghost_var_num_man.forget(
              v, odi_map_t::object_sum_raw_val(prod_value_ref));
          m_ghost_var_num_man.forget(
              v, odi_map_t::object_cache_raw_val(prod_value_ref));
          m_ghost_var_eq_man.forget(
              v, odi_map_t::object_eq_raw_val(prod_value_ref));
          if (auto v_w = get_write_region(v)) {
            m_ghost_var_eq_man.forget(
                *v_w, odi_map_t::object_eq_raw_val(prod_value_ref));
          }

          m_odi_map.set(*id_opt,
                        odi_domain_product_t(std::move(prod_info_ref),
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

public:
  /// @brief forget a set of variables
  /// @param variables a set of program variables whose properties are dropped
  /// @note the batched fold of \c operator-=: forget(S) behaves exactly like
  /// applying `*this -= v` for each v in S (the law region_domain follows),
  /// with the guards and the reduction performed once up front.
  void forget(const variable_vector_t &variables) override {
    OBJECT_DOMAIN_SCOPED_STATS(".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    apply_reduction_based_on_flags();

    CRAB_LOG("object-forget", crab::outs() << "Forgetting (";
             object_domain_impl::print_vector(crab::outs(), variables);
             crab::outs() << ")=" << *this << "\n";);

    for (auto &v : variables) {
      forget_var(v);
    }

    CRAB_LOG("object-forget", crab::outs()
                                  << "After Forget:" << *this << "\n";);
  }

  void project(const variable_vector_t &variables) override {
    OBJECT_DOMAIN_SCOPED_STATS(".project");

    if (is_bottom() || is_top()) {
      return;
    }

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
    for (const auto &kv : flds_by_id_map) {
      const obj_id_t &id = kv.first;
      const variable_vector_t &flds = kv.second;
      ref_vars.push_back(get_or_insert_base_addr(id));
      // also keep the MRU base address: the cache's ownership is recorded as
      // ref_base == <id>_mru_base (see test_ref_refer_mru_object), so
      // dropping it would turn every later cache probe into a miss
      ref_vars.push_back(get_or_insert_base_addr(id, true /*is_cache*/));
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (!prod_ref) {
        continue;
      }
      object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
      object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
      const small_range &num_refs = prod_info_ref.refcount_val();
      // non-singleton object
      // project based on corresponding field(s)
      m_ghost_var_num_man.project(
          flds, odi_map_t::object_sum_raw_val(prod_value_ref));
      m_ghost_var_num_man.project(
          flds, odi_map_t::object_cache_raw_val(prod_value_ref));
      m_ghost_var_eq_man.project(flds,
                                 odi_map_t::object_eq_raw_val(prod_value_ref));

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
    // pending write-region ghosts of renamed fields (renamed in the
    // field-equality domain only, where they live)
    std::unordered_map<obj_id_t,
                       std::pair<variable_vector_t, variable_vector_t>>
        renamed_wflds_by_id_map;

    // filtering incoming input variable vectors based on their types
    // since we split them into several subdomains
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &old_v = from[i];
      const variable_t &new_v = to[i];
      if (old_v.get_type() != new_v.get_type()) {
        CRAB_ERROR(domain_name(), "::", __func__, " ", old_v, " and ", new_v,
                   " must preserve the same type");
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
          // register the new name as a field of the same object, so
          // get_obj_id(new_v) resolves after the rename. The old entry is
          // kept: the field->id map is shared across states and other
          // states may still refer to the old name.
          update_fields_id_map(new_v, *id_opt);
          // derived ghosts are renamed along with the field: a pending
          // store's write-region ghost (field-equality domain) and, for
          // reference regions, the base-address ghost (address domain)
          if (get_write_region(old_v) != boost::none) {
            auto &wflds = renamed_wflds_by_id_map[*id_opt];
            wflds.first.push_back(get_or_insert_write_region(old_v));
            wflds.second.push_back(get_or_insert_write_region(new_v));
          }
          if (old_v.get_type().is_reference_region()) {
            from_ref_vars.push_back(get_or_insert_base_addr(old_v));
            to_ref_vars.push_back(get_or_insert_base_addr(new_v));
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
    // rename registers/references in the base domain (their ghost
    // variables), not only their equalities
    m_ghost_var_num_man.rename(from_non_odi_vars, to_non_odi_vars, m_base_dom);
    m_ghost_var_eq_man.rename(from_non_odi_vars, to_non_odi_vars,
                              m_eq_regs_dom);

    // renaming fields
    for (const auto &kv : renamed_flds_by_id_map) {
      const obj_id_t &id = kv.first;
      const variable_vector_t &from_flds = std::get<0>(kv.second);
      const variable_vector_t &to_flds = std::get<1>(kv.second);
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (!prod_ref) {
        continue;
      }
      object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
      object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
      const small_range &num_refs = prod_info_ref.refcount_val();
      // update base_dom or odi_map
      // fields live in the odi map
      m_ghost_var_num_man.rename(from_flds, to_flds,
                                 odi_map_t::object_sum_raw_val(prod_value_ref));
      m_ghost_var_num_man.rename(
          from_flds, to_flds, odi_map_t::object_cache_raw_val(prod_value_ref));
      m_ghost_var_eq_man.rename(from_flds, to_flds,
                                odi_map_t::object_eq_raw_val(prod_value_ref));
      auto wit = renamed_wflds_by_id_map.find(id);
      if (wit != renamed_wflds_by_id_map.end()) {
        m_ghost_var_eq_man.rename(wit->second.first, wit->second.second,
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
      // apply any pending cache/base reduction first (as entails() does), so
      // facts still sitting in a dirty cache are exported too
      auto opt_res = cow_apply_reduction();
      const object_domain_t &reduced =
          boost::get_optional_value_or(opt_res, *this);
      return reduced.m_base_dom.to_linear_constraint_system();
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

    if (name == "is_dereferenceable") {
      apply_reduction_based_on_flags(true);
    } else {
      apply_reduction_based_on_flags();
    }

    if (name == "regions_from_memory_object") {
      // pass region variables into object field map
      // the intrinsics is only added in clam if the object has more than one
      // region.
      assert(inputs.size() >= 1);
      unsigned obj_id_idx = 0;
      for (unsigned sz = inputs.size(); obj_id_idx < sz; ++obj_id_idx) {
        error_if_not_rgn(inputs[obj_id_idx]); // this should not happen
        if (inputs[obj_id_idx].get_type().is_unknown_region()) {
          continue;
        } else {
          break;
        }
      }
      // error_if_not_rgn(inputs[0]);
      if (obj_id_idx == inputs.size()) {
        return;
      }
      auto obj_id_opt = get_obj_id(inputs[obj_id_idx].get_variable());
      obj_id_t obj_id =
          obj_id_opt ? *obj_id_opt : inputs[obj_id_idx].get_variable();
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          continue;
        }
        // Note that, region initialization is before the intrinsic calls
        // Any obj info set up are not for object id will be removed.
        auto old_id_opt = get_obj_id(inputs[i].get_variable());
        if (old_id_opt && obj_id != *old_id_opt) {
          m_odi_map -= *old_id_opt;
          m_addrs_dom -= get_or_insert_base_addr(*old_id_opt);
        }
        update_fields_id_map(inputs[i].get_variable(), obj_id);
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
      std::unordered_set<obj_id_t> id_set;
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // an object's field list may mix known and unknown regions
          // (regions_from_memory_object skips unknowns the same way);
          // the known fields still identify their objects to commit
          continue;
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
        obj_info.sumpresence_val() = obj_info.objinit_val();
        m_odi_map.set(id, odi_domain_product_t(std::move(obj_info),
                                               std::move(obj_value)));
      }
    } else if (name == "copy_memory_object") {
      assert(inputs.size() >= 1);
      assert(outputs.size() >= 1);
      std::vector<variable_t> input_vars;
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        // unknown regions are passed through: object_copy filters
        // src/dst pairs of unknown regions itself, so one unknown field
        // no longer discards the whole copy
        input_vars.push_back(inputs[i].get_variable());
      }
      object_copy(input_vars, outputs, *this, *this);
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

  std::string domain_name() const override { return "Object"; }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      CRAB_LOG("object-print", o << "("
                                 << "Flds-id map=";
               print_flds_id_map(o); o << ",\n"
                                       << "BaseDom=";
               m_ghost_var_num_man.write(o, m_base_dom);
               o << ",\nAddrs=" << m_addrs_dom;
               o << ",\nRegs=" << m_eq_regs_dom; o << ","; object_write(o);
               o << ")\n"; return;);
      o << "{}";
    } else {
      CRAB_LOG("object-print", o << "("
                                 << "Flds-id map=";
               print_flds_id_map(o); o << ",\n"
                                       << "BaseDom=";
               m_ghost_var_num_man.write(o, m_base_dom);
               o << "\nAddrs=" << m_addrs_dom; o << "\nRegs=" << m_eq_regs_dom;
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

  /// @brief check whether \p rgn is a region belonging to some equivalence
  /// class; when \p which is non-null, additionally mark the index of the
  /// matched class
  static bool rgn_in_group(const variable_t &rgn, const classes_t &cls,
                           std::vector<bool> *which = nullptr) {
    if (!rgn.get_type().is_region()) {
      return false;
    }
    for (unsigned i = 0, sz = cls.size(); i < sz; ++i) {
      if (std::find(cls[i].begin(), cls[i].end(), rgn) != cls[i].end()) {
        if (which) {
          (*which)[i] = true;
        }
        return true;
      }
    }
    return false;
  }

  /// @brief copy object based on one dsa node to another
  /// @param[in] src_rgns a set of regions from src dsa node
  /// @param[in] dst_rgns a set of regions from dst dsa node
  /// @param[in] src_dom the source abstract state
  /// @param[in,out] dst_dom the destination abstract state
  static void object_copy(const equiv_class_regions_t &src_rgns,
                          const equiv_class_regions_t &dst_rgns,
                          const object_domain_t &src_dom,
                          object_domain_t &dst_dom) {
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
    obj_id_t dst_id = dst_id_opt
                          ? *dst_id_opt
                          : dst_dom.create_new_obj_id(dst_rgns_no_unknown[0]);
    // register only the known regions: an unknown region must not become a
    // "field" of dst's object in the shared field->id map (it would, e.g.,
    // change forget()'s all-fields-listed drop rule for this object forever)
    for (auto &fld : dst_rgns_no_unknown) {
      dst_dom.update_fields_id_map(fld, dst_id);
    }
    const odi_domain_product_t *prod_ref = src_dom.m_odi_map.find(src_id);
    if (!prod_ref) {
      // default is top
      return;
    }
    object_info_t prod_info_ref = odi_map_t::object_info_val(*prod_ref);
    object_value_t prod_value_ref = odi_map_t::object_odi_val(*prod_ref);
    // const view for read-only probes (see region_copy)
    const object_value_t &c_prod_value_ref = prod_value_ref;
    const small_range &num_refs = prod_info_ref.refcount_val();

    assert(!num_refs.is_zero());

    // num_refs > 1
    for (unsigned i = 0, len = src_rgns_no_unknown.size(); i < len; ++i) {
      auto src_rgn_gvars_opt = src_dom.get_num_gvars(src_rgns_no_unknown[i]);
      auto src_rgn_eq_gvars_opt = src_dom.get_eq_gvars(src_rgns_no_unknown[i]);
      if (src_rgn_gvars_opt == boost::none ||
          src_rgn_eq_gvars_opt == boost::none) {
        continue;
      }
      ghost_variables_t dst_rgn_gvars =
          dst_dom.get_or_insert_gvars(dst_rgns_no_unknown[i]);
      auto dst_rgn_eq_gvars =
          dst_dom.get_or_insert_eq_gvars(dst_rgns_no_unknown[i]);
      if (!odi_map_t::object_sum_raw_val(c_prod_value_ref).is_top()) {
        dst_rgn_gvars.forget(odi_map_t::object_sum_raw_val(prod_value_ref));
        src_rgn_gvars_opt.value().expand(
            odi_map_t::object_sum_raw_val(prod_value_ref), dst_rgn_gvars);
      }
      if (!odi_map_t::object_cache_raw_val(c_prod_value_ref).is_top()) {
        dst_rgn_gvars.forget(odi_map_t::object_cache_raw_val(prod_value_ref));
        src_rgn_gvars_opt.value().expand(
            odi_map_t::object_cache_raw_val(prod_value_ref), dst_rgn_gvars);
      }
      if (!odi_map_t::object_eq_raw_val(c_prod_value_ref).is_top()) {
        dst_rgn_eq_gvars.forget(odi_map_t::object_eq_raw_val(prod_value_ref));
        src_rgn_eq_gvars_opt.value().expand(
            odi_map_t::object_eq_raw_val(prod_value_ref), dst_rgn_eq_gvars);
      }
    }
    dst_dom.m_ghost_var_num_man.project(
        dst_rgns_no_unknown, odi_map_t::object_sum_raw_val(prod_value_ref));
    dst_dom.m_ghost_var_num_man.project(
        dst_rgns_no_unknown, odi_map_t::object_cache_raw_val(prod_value_ref));
    dst_dom.m_ghost_var_eq_man.project(
        dst_rgns_no_unknown, odi_map_t::object_eq_raw_val(prod_value_ref));

    if (auto src_base_addr_opt = src_dom.get_base_addr(src_id)) {
      dst_dom.m_addrs_dom.expand(*src_base_addr_opt,
                                 dst_dom.get_or_insert_base_addr(dst_id));
    }
    // NOTE: self-copy (src_dom == dst_dom, see the copy_memory_object
    // intrinsic) is allowed BECAUSE this set is the last access: no read
    // through prod_ref may follow it (odi_map::find's pointer is
    // invalidated by set -- see S16)
    dst_dom.m_odi_map.set(dst_id,
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

    CRAB_LOG("inter-restrict1", crab::outs()
                                    << "Compute the entry state at callsite "
                                    << callsite.get_function()
                                    << ". Caller: " << caller_dom << "\n";);

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
        // the region in the callee cannot be mapped to two distinct regions in
        // the caller. This mean there is no case that callee's regions belong
        // to a dsa node but caller's do not (picture 8 in SEADSA-SAS17 paper).
        // For the case that a src region belonging to a dsa node but dst does
        // not This differs from copying objects, we implemented it in
        // region_copy. Here, we use object_copy for both src regions and dst
        // regions belong to some dsa nodes.
        if (rgn_in_group(formal, callee_formal_cls) &&
            rgn_in_group(actual, caller_actual_cls)) {
          continue;
        }
        CRAB_LOG("inter-restrict", crab::outs()
                                       << "\t" << formal << ":"
                                       << formal.get_type() << " and " << actual
                                       << ":" << actual.get_type() << "\n";);
        inter_transformers_impl::typed_assign(caller, formal, actual);
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
      object_copy(actual_rgns, formal_rgns, caller_dom, caller);
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

    CRAB_LOG("inter-extend1", crab::outs()
                                  << "Compute the return state at callsite "
                                  << callsite.get_function()
                                  << ". Callee: " << callee_dom << "\n";);

    if (is_bottom()) {
      return;
    }

    if (callee_dom.is_bottom()) {
      set_to_bottom();
      return;
    }

    object_domain_t callee_at_exit(callee_dom);

    if (callee_at_exit.is_top()) {
      callee_at_exit = make_top();
    }

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

    // 1. forget output parameters at the callsite
    // e.g reassignment in case
    forget(caller_out_params);

    // 2. propagate from callee's outputs to caller's lhs of the callsite
    for (unsigned i = 0, e = callee_out_params.size(); i < e; ++i) {
      const variable_t &out_formal = callee_out_params[i];
      const variable_t &out_actual = caller_out_params[i];
      if (!(out_formal == out_actual)) {
        if (rgn_in_group(out_actual, caller_lhs_cls) &&
            rgn_in_group(out_formal, callee_ret_cls)) {
          continue;
        }
        CRAB_LOG("inter-extend", crab::outs()
                                     << "Unifying output " << out_actual
                                     << ":= " << out_formal << "\n";);
        inter_transformers_impl::typed_assign(callee_at_exit, out_actual,
                                              out_formal);
      }
    }

    // 3. propagate from callee's inputs to caller's inputs at callsite
    // This is needed to propagate up new input-output relationships:
    // (1) e.g. V1, V2 = call function_f(V2), where V1 is a new variable,
    // V2 is input and output variable;
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
        tmp.begin(), tmp.end()); // linear because tmp is sorted

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
        if (rgn_in_group(in_actual, caller_actual_cls, &cls_check)) {
          continue;
        }
        inter_transformers_impl::typed_assign(callee_at_exit, in_actual,
                                              in_formal);
      }
    }

    // 4. propagate abstract objects
    // unify out objects
    for (unsigned i = 0, len = caller_lhs_cls.size(); i < len; ++i) {
      const equiv_class_regions_t &lhs_rgns = caller_lhs_cls[i];
      const equiv_class_regions_t &ret_rgns = callee_ret_cls[i];
      object_copy(ret_rgns, lhs_rgns, callee_at_exit, *this);
    }

    // unify in objects in case (3)
    for (unsigned i = 0, len = caller_actual_cls.size(); i < len; ++i) {
      if (cls_check[i]) {
        const equiv_class_regions_t &formal_rgns = callee_formal_cls[i];
        const equiv_class_regions_t &actual_rgns = caller_actual_cls[i];
        object_copy(formal_rgns, actual_rgns, callee_at_exit, *this);
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
  /// @brief the symbol bound to \p variable in \p eq_dom, or a fresh
  /// (still unbound) symbol if it has none
  /// @param eq_dom an equality domain value
  /// @param variable a field
  /// @return the class id of the variable's equivalence class, such as #id1
  /// @note the fresh symbol is NOT inserted into \p eq_dom; the odi-map
  /// operation receiving it performs the binding
  static usymb_t get_symbol_or_fresh(const eq_register_domain_t &eq_dom,
                                     const flds_dom_variable_t &variable) {
    if (boost::optional<usymb_t> symb_opt = eq_dom.get_class_id(variable)) {
      return *symb_opt;
    } else {
      return eq_domain_value_t::fresh_class_id();
    }
  }

  /***************** Cache operations *****************/
  /// @brief commit cache contents into summary if it is dirty, perform
  /// reduction to a specific base domain and reset cache
  /// @param base_dom the base domain the cache is reduced against
  /// @param prod object value stores a product of <SUM dom, cache dom, EQ dom>
  /// @param obj_info_ref the object info stores reference count & cache status
  /// @param id an object id
  void commit_cache_if_dirty(base_abstract_domain_t &base_dom,
                             object_value_t &prod,
                             const object_info_t &obj_info_ref,
                             const obj_id_t &id) const {
    // (1) perform reduction by transferring regs-flds' constraints between
    // cache
    //  domain and base domain
    if (obj_info_ref.cache_reg_loaded_val()) {
      apply_reduction_from_object_to_base(base_dom, prod, id);
    }
    if (obj_info_ref.cache_reg_stored_val()) {
      base_abstract_domain_t regs_only_base =
          project_onto_eq_vars(base_dom, m_eq_regs_dom);
      apply_reduction_from_base_to_object(regs_only_base, prod, id);
    }
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

  /// @brief a wrapper method for adapting the interface of domain reduction
  /// from odi map to base
  /// @param[in] prod the abstract value corresponding to object \p id
  /// @param[in] id the object id
  void apply_reduction_from_object_to_base(const object_value_t &prod,
                                           const obj_id_t &id) {
    apply_reduction_from_object_to_base(m_base_dom, prod, id);
  }

  /// @brief a wrapper method for adapting the interface of domain reduction
  /// from base to odi map
  /// @param[in,out] prod the abstract value corresponding to object \p id
  /// @param[in] id the object id
  void apply_reduction_from_base_to_object(object_value_t &prod,
                                           const obj_id_t &id) const {
    base_abstract_domain_t regs_only_base =
        project_onto_eq_vars(m_base_dom, m_eq_regs_dom);
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