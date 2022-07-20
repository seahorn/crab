#pragma once

#include <boost/optional.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/uf_domain.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/varname_factory.hpp>

#include <crab/domains/object/object_info.hpp>
#include <crab/domains/region/ghost_variable_manager.hpp>
#include <crab/domains/region/tags.hpp>

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

template <class AbsDom>
void assign_interval(AbsDom &dom, const typename AbsDom::variable_t v,
                     typename AbsDom::interval_t ival) {
  dom -= v;
  boost::optional<typename AbsDom::number_t> lb = ival.lb().number();
  boost::optional<typename AbsDom::number_t> ub = ival.ub().number();
  if (lb) {
    dom += (*lb <= v);
  }
  if (ub) {
    dom += (v <= *ub);
  }
}

template <class AbsDom> struct object_equal_to {
  bool operator()(const AbsDom &v1, const AbsDom &v2) const {
    return &v1 == &v2;
  }
};

struct term_op_hash {
  std::size_t operator()(const term::term_operator_t &k) const {
    return std::hash<uint32_t>()(k.value());
  }
};

template <class EUFDom> term::term_operator_t make_fresh_symbol(EUFDom &dom) {
  term::term_operator_t t =
      dom.make_uf(term::term_op_val_generator_t::get_next_val());
  return t;
}

template <class Num, class Ftor>
boost::optional<Ftor>
extract_uf_symbol_from_term(const term::term<Num, Ftor> *t) {
  if (t && t->kind() == term::term_kind::TERM_APP) {
    auto fterm = static_cast<const term::ftor_term<Num, Ftor> *>(t);
    if (fterm->args.empty())
      return term::term_ftor(t);
  }
  return boost::none;
}

template <class EUFDom>
term::term_operator_t
get_or_make_term_symbol(EUFDom &dom, const typename EUFDom::variable_t &v) {
  using term_t = typename EUFDom::term_t;
  const term_t *t = dom.get_term(v);
  if (t != nullptr && t->kind() == term::term_kind::TERM_APP) {
    return term::term_ftor(t);
  }
  return make_fresh_symbol(dom);
}

template <class Var, class Num, class Ftor>
bool is_two_vars_have_same_term(const Var &v1,
                                const term::term<Num, Ftor> *t1_ptr,
                                const Var &v2,
                                const term::term<Num, Ftor> *t2_ptr) {
  using term_op_opt_t = boost::optional<term::term_operator_t>;
  if (!t1_ptr || !t2_ptr) {
    return false;
  }
  term_op_opt_t top1_opt = extract_uf_symbol_from_term(t1_ptr);
  term_op_opt_t top2_opt = extract_uf_symbol_from_term(t2_ptr);
  if (top1_opt == boost::none || top1_opt == boost::none) {
    return false;
  }
  return (*top1_opt).value() == (*top2_opt).value();
}

template <typename TType>
void print_vector(crab::crab_os &o, const std::vector<TType> &vec) {
  typename std::vector<TType>::const_iterator it;
  o << "[";
  for (it = vec.begin(); it != vec.end(); it++) {
    if (it != vec.begin())
      o << ",";
    o << (*it);
  }
  o << "]\n";
}
} // end namespace object_domain_impl

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

  // type name for object domains
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

  // type name for address and register domains
  using uf_domain_t = uf_domain<number_t, varname_t>;
  using usymb_t = typename uf_domain_t::term_t;
  using address_abstract_domain_t = uf_domain_t;
  using uf_register_abstract_domain_t = uf_domain_t;
  using uf_fields_abstract_domain_t = uf_domain_t;

  using object_info_t = typename object_domain_impl::object_info;
  using symb_flds_regs_map_t = std::unordered_map<
      term::term_operator_t,
      std::pair<base_dom_variable_vector_t, base_dom_variable_vector_t>,
      object_domain_impl::term_op_hash>;

  // FIXME: the current solution does not support different dom operations
  static_assert(
      std::is_same<base_abstract_domain_t, field_abstract_domain_t>::value,
      "base_abstract_domain_t and field_abstract_domain_t must be the same");

  static_assert(std::is_same<uf_register_abstract_domain_t,
                             uf_fields_abstract_domain_t>::value,
                "uf_register_abstract_domain_t and uf_fields_abstract_domain_t "
                "must be the same");

  /**------------------ Begin type definitions ------------------**/
  // This domain is based on the region based memory model,
  // an object, more specifically an abstract object, is an object represents
  // a number of concrete objects with the same pattern.
  // An object is spawn from a set of fields (or regions).
  // In addition, we model the memory architecture with a cache.
  // The cache tracks the most recently used (MRU) object by memory load / store
  // where the MRU object is one of the summarzied concrete objects.
  // The operations such as memory load or store are precisely abstracted
  // if the requested properties can be found in the cache.

  // Object identifier / object id
  using obj_id_t = variable_t;

  // ODI: <DOM1, DOM2, uf_DOM>
  // A product domain tracks the followings:
  // Summary domain: An domain keeps properties for all summarized objects
  // Cache domain: An domain for most recetnly used object
  // uf fields domain:
  //    keeps what uninterpreted symbols assign to fields
  //    NOTE that: we could keep a map from each field to
  //               an uninterpreted symbol
  using odi_domain_product_t =
      basic_domain_product2<field_abstract_domain_t,
                            basic_domain_product2<field_abstract_domain_t,
                                                  uf_fields_abstract_domain_t>>;

  // ODI map: object id -> ODI
  // Map from an object id to object's abstract value.
  // An object id indicates which kinds of object is.
  // the two domains must maintain the same id before the lattice
  // operations
  using odi_map_t = ikos::separate_domain<
      obj_id_t, odi_domain_product_t,
      object_domain_impl::object_equal_to<odi_domain_product_t>>;

  // Object infos: object id -> <small range, boolean value, boolean value>
  // Map an object id to finite reduced product domains
  // for inferring object info:
  //  object counting, object cache used, object cache dirty
  using obj_info_env_t =
      ikos::separate_domain<obj_id_t, object_domain_impl::object_info>;

  // Object fields to id map
  // Map region variables to object's id
  // This map is constructed during intrinsic call and share in common
  using obj_flds_id_map_t =
      std::shared_ptr<std::unordered_map<variable_t, obj_id_t>>;

  // References to corresponding base addresses
  // Map reference variables to variables used in EUF domain
  // This map is constructued for address domain, we keep track
  // the variables created for references.
  // The map is constructed dynamically but share accress all states
  using refs_base_addrs_map_t =
      std::shared_ptr<std::unordered_map<variable_t, variable_t>>;

  // Management of ghost variables
  // We model two types of variables in object domain
  // One is for region variable, another is for reference variable.
  // The former one is used to perform the reduction between region and register
  // , the later one is used to infer sea.is_deref.
  // Note that, we are using fixed naming where
  // variable_t and subdomain::variable_t have the same type
  // TODO: if the type of the domain for object is different from
  // the type for base domain, the ghost variable manager should split into two.
  using ghost_var_man_t =
      typename region_domain_impl::ghost_variable_manager_with_fixed_naming<
          object_domain_t, base_abstract_domain_t>;
  using ghost_variables_t = typename ghost_var_man_t::ghost_variables_t;
  using ghost_variable_vector_t = typename std::vector<ghost_variables_t>;
  using ghost_variable_kind = typename ghost_variables_t::ghost_variable_kind;
  /**------------------ End type definitions ------------------**/

  /**------------------ Begin class field definitions ------------------**/

  // a special symbol to represent bottom state of the object domain
  // object domain is bottom if this flag is set or all subdomains are bottom.
  bool m_is_bottom;

  // The abstract state definition:
  // Base domain:
  // Domain for register, singleton object
  base_abstract_domain_t m_base_dom;

  // A map from each region variable to corresponding object domain
  odi_map_t m_odi_map;

  // A domain to keep equalities over base addresses of objects.
  address_abstract_domain_t m_addrs_dom;

  // A domain to keep what uninterpreted symbols assign to registers if those
  // are equal to some objects' fields
  uf_register_abstract_domain_t m_uf_regs_dom;

  // Map object ids to a tuple of (Count,Init,Type)
  // In addition, any regions are not part of an object treat introduced by
  // the intrinsic method, we treat them as objects with single field.
  //
  // Count: count how many allocations are owned by an abstract object.
  // Each call to ref_make models a new allocation.
  // This allows us to decide when only if one address per object
  // (i.e., singleton object).
  //
  // Used: indicate whether the cache domain is used.
  //
  // Dirty: indicate whether the cache domain is dirty.
  obj_info_env_t m_obj_info_env;

  // Map region variables to an object id for an abstract object
  // To determine the object id,
  // we choose the first region variable passed by the intrinsic method.
  // If any regions that not indicating by the intrinsic, we add those regions
  // into this method during the evaluation of make_ref.
  // TODO: for now, we did not cover unknown regions.
  obj_flds_id_map_t m_flds_id_map;

  // Map reference variables to base addresses variables
  // The base addresses
  // In addition, each cache that used for an abstract object
  // has a special base address variable for the MRU object
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
  ghost_var_man_t m_ghost_var_man;
  /**------------------ End class field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  // Constructor Definition
  object_domain(base_abstract_domain_t &&base_dom, odi_map_t &&odi_map,
                address_abstract_domain_t &&addrs_dom,
                uf_register_abstract_domain_t &&uf_regs_dom,
                obj_info_env_t &&obj_info_env, obj_flds_id_map_t &&flds_id_map,
                refs_base_addrs_map_t &&refs_base_addrs_map,
                ghost_var_man_t &&ghost_var_man)
      : m_is_bottom(base_dom.is_bottom()), m_base_dom(base_dom),
        m_odi_map(odi_map), m_addrs_dom(addrs_dom), m_uf_regs_dom(uf_regs_dom),
        m_obj_info_env(obj_info_env), m_flds_id_map(flds_id_map),
        m_refs_base_addrs_map(refs_base_addrs_map),
        m_ghost_var_man(ghost_var_man) {}

  std::function<variable_type(const variable_t &)> get_type() const {
    auto fn = [](const variable_t &v) { return v.get_type(); };
    return fn;
  }

  ghost_variables_t get_or_insert_gvars(const variable_t &v) {
    return m_ghost_var_man.get_or_insert(v);
  }

  boost::optional<ghost_variables_t> get_gvars(const variable_t &v) const {
    return m_ghost_var_man.get(v);
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

  base_dom_linear_expression_t
  rename_linear_expr(const linear_expression_t &e) {
    return m_ghost_var_man.rename_linear_expr(e);
  }

  base_dom_linear_constraint_t
  rename_linear_cst(const linear_constraint_t &cst) {
    return m_ghost_var_man.rename_linear_cst(cst);
  }

  base_dom_linear_constraint_system_t
  rename_linear_cst_sys(const linear_constraint_system_t &csts) {
    return m_ghost_var_man.rename_linear_cst_sys(csts);
  }

  base_dom_linear_constraint_t
  convert_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst,
                                ghost_variable_kind kind) {
    return m_ghost_var_man.ghosting_ref_cst_to_linear_cst(ref_cst, kind);
  }

  // Self join
  // Perform *this = join(*this, right)
  void self_join(const object_domain_t &right) {

    // The join operation is sound if:
    // Invariant 1: the singleton object in base domain, the object info implies
    // the number of reference is exactly one
    // Invariant 2: the non-singleton obj in odi map, the object info implies
    // the number of reference is > 1
    // I.e. for each abstract object,
    //  there is only three states of object reference: 0, 1, and [1, +oo]
    // If those invariants, holds, the following join operation is sound:

    // For each abstract object, generally, we have four cases at the join:
    // The following example only shows properties related to
    // just one abstract object:

    // Case 1: not singleton, odi map contains a field domain for that object
    // Solution: we join the base dom and odi map
    // E.g.
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 = 3, V_2 <= 5 }
    //  s1 U s2 : Base in s1 U Base in s2,
    //            Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }

    // Case 2: singleton object in base domain
    //         odi map does not contain a field domain for that object
    // Solution: we join the base dom and odi map
    // E.g.
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { V_1 = 2, V_2 = 4 }, Object (V_1, V_2) = empty
    //  s1 U s2: Base in s1 U Base in s2
    //           Object (V_1, V_2) = empty

    // Case 3:
    // one state: singleton object in base domain,
    //  odi map does not contain a field domain for that object
    // another: not singleton, odi map contains a field domain for that object
    // Solution: project singleton object, join with non-singleton in the odi
    // map
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
    //    Base = (Base in s1).project(V_1, V_2) U Base in s2
    //    Object (V_1, V_2) = { V_1 <= 3, V_2 <= 6 }

    // Other cases (simple):
    // one state: exists such an abstract object
    // another: does not exist that abstract object
    // Solution: we join the base dom and odi map
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

    // Overall, for both cases, the only third case requires move singleton
    // into odi map. Other cases require join operation on each subdomain as
    // before.

    // 1. join the parts of each state that can use default join operation
    obj_info_env_t out_obj_info_env(m_obj_info_env | right.m_obj_info_env);
    address_abstract_domain_t out_addrs_dom(m_addrs_dom | right.m_addrs_dom);

    // m_flds_id_map and m_refs_base_addrs_map do not require join,
    // they share accross all states. No need to update
    assert(m_flds_id_map == right.m_flds_id_map);
    if (right.m_refs_base_addrs_map && !m_refs_base_addrs_map) {
      m_refs_base_addrs_map = right.m_refs_base_addrs_map;
    }

    // 2. join the odi map
    base_abstract_domain_t &l_base_dom = m_base_dom;
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    odi_map_t out_odi_map;
    for (auto kv : m_obj_info_env) {
      const obj_id_t &id = kv.first;
      const object_info_t &l_obj_info = kv.second;
      const small_range &l_num_refs = l_obj_info.refcount_val();
      auto r_obj_info_ref = right.m_obj_info_env.find(id);
      if (!r_obj_info_ref) {
        continue;
      }
      const object_info_t &r_obj_info = (*r_obj_info_ref);
      const small_range &r_num_refs = r_obj_info.refcount_val();
      const odi_domain_product_t *l_prod_ref = m_odi_map.find(id);
      const odi_domain_product_t *r_prod_ref = right.m_odi_map.find(id);
      // Step 1. Join the odi map if the object is non-singleton in both states
      // For current abstract object, when both states refer odi maps,
      // We need to check if cache(s) is (are) used before:
      // (1) for summary domain, we just perform the join.
      // (2) for cache domain, we can join them if both caches
      //     refer the same memory object
      // if yes, we perform the join;
      // otherwise, we have to commit cache(s) before (1) and (2).
      // To determine whether both caches refer the same object, we can only
      // know if those caches are loaded from the common state.
      // In current implementation, we determine it by checking whether
      // uninterpreted symbols of the mru object are the same.
      if (l_num_refs == small_range::oneOrMore() &&
          r_num_refs == small_range::oneOrMore()) {
        const boolean_value l_used = l_obj_info.cacheused_val();
        const boolean_value r_used = r_obj_info.cacheused_val();
        odi_domain_product_t l_prod =
            l_prod_ref ? (*l_prod_ref) : odi_domain_product_t();
        odi_domain_product_t r_prod =
            r_prod_ref ? (*r_prod_ref) : odi_domain_product_t();
        if (l_used.is_true() && r_used.is_true()) {
          if (mru_refer_same_object(id, m_addrs_dom, right.m_addrs_dom) ==
              false) {
            commit_cache_if_dirty(l_base_dom, m_uf_regs_dom, l_prod,
                                  &l_obj_info, id);
            commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                  &r_obj_info, id);
          }
        } else if (l_used.is_true()) {
          commit_cache_if_dirty(l_base_dom, m_uf_regs_dom, l_prod, &l_obj_info,
                                id);
        } else if (r_used.is_true()) {
          commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                &r_obj_info, id);
        }
        // join the summary, cache and flds part
        odi_domain_product_t o_prod = l_prod | r_prod;
        out_odi_map.set(id, o_prod);
      }
      // Step 2. Handle the case for joining singleton with non-singleton
      // In this case, we need to commit the cache if the cache is used
      // before joining singleton with non-singleton.
      else if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
        // left state: singleton
        // right state: non-singleton
        join_or_widen_singleton_with_non_singleton(
            id, *this, right, r_base_dom, out_odi_map, true /* is_join*/);
      } else if (r_num_refs.is_one() &&
                 l_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(
            id, right, *this, l_base_dom, out_odi_map, true /* is_join*/);
      }
    }

    // 3. join the base domain
    base_abstract_domain_t out_base_dom(l_base_dom | r_base_dom);

    // 4. join the regs domain
    uf_register_abstract_domain_t out_uf_regs_dom(m_uf_regs_dom |
                                                  right.m_uf_regs_dom);

    // update current state
    std::swap(m_obj_info_env, out_obj_info_env);
    std::swap(m_addrs_dom, out_addrs_dom);
    std::swap(m_uf_regs_dom, out_uf_regs_dom);
    std::swap(m_base_dom, out_base_dom);
    std::swap(m_odi_map, out_odi_map);
  }

  // join / widening two abstract states, return a new join / widening state
  object_domain_t join_or_widening(const object_domain_t &left,
                                   const object_domain_t &right,
                                   const bool is_join) const {

    // 1. join the parts of each state that can use default join operation
    obj_info_env_t out_obj_info_env(
        is_join ? left.m_obj_info_env | right.m_obj_info_env
                : left.m_obj_info_env || right.m_obj_info_env);
    address_abstract_domain_t out_addrs_dom(
        is_join ? left.m_addrs_dom | right.m_addrs_dom
                : left.m_addrs_dom || right.m_addrs_dom);

    // m_flds_id_map and m_refs_base_addrs_map do not require join,
    // they share accross all states. No need to update
    assert(left.m_flds_id_map == right.m_flds_id_map);

    // 2. join the odi map
    // We only perform the common join
    base_abstract_domain_t l_base_dom = left.m_base_dom;
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    odi_map_t out_odi_map;
    for (auto kv : left.m_obj_info_env) {
      const obj_id_t &id = kv.first;
      const object_info_t &l_obj_info = kv.second;
      const small_range &l_num_refs = l_obj_info.refcount_val();
      auto r_obj_info_ref = right.m_obj_info_env.find(id);
      if (!r_obj_info_ref) {
        continue;
      }
      const object_info_t &r_obj_info = (*r_obj_info_ref);
      const small_range &r_num_refs = r_obj_info.refcount_val();
      const odi_domain_product_t *l_prod_ref = left.m_odi_map.find(id);
      const odi_domain_product_t *r_prod_ref = right.m_odi_map.find(id);
      // Step 1. Join the odi map if the object is non-singleton in both states
      if (l_num_refs == small_range::oneOrMore() &&
          r_num_refs == small_range::oneOrMore()) {
        const boolean_value l_used = l_obj_info.cacheused_val();
        const boolean_value r_used = r_obj_info.cacheused_val();
        odi_domain_product_t l_prod =
            l_prod_ref ? (*l_prod_ref) : odi_domain_product_t();
        odi_domain_product_t r_prod =
            r_prod_ref ? (*r_prod_ref) : odi_domain_product_t();
        if (l_used.is_true() && r_used.is_true()) {
          if (mru_refer_same_object(id, left.m_addrs_dom, right.m_addrs_dom) ==
              false) {
            commit_cache_if_dirty(l_base_dom, left.m_uf_regs_dom, l_prod,
                                  &l_obj_info, id);
            commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                  &r_obj_info, id);
          }
        } else if (l_used.is_true()) {
          commit_cache_if_dirty(l_base_dom, left.m_uf_regs_dom, l_prod,
                                &l_obj_info, id);
        } else if (r_used.is_true()) {
          commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                &r_obj_info, id);
        }
        // join the summary, cache and flds part
        odi_domain_product_t o_prod =
            is_join ? l_prod | r_prod : l_prod || r_prod;
        out_odi_map.set(id, o_prod);
      }
      // Step 2. Handle the case for joining singleton with non-singleton
      else if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
        // left state: singleton
        // right state: non-singleton
        join_or_widen_singleton_with_non_singleton(id, left, right, r_base_dom,
                                                   out_odi_map, is_join);
      } else if (r_num_refs.is_one() &&
                 l_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(id, right, left, l_base_dom,
                                                   out_odi_map, is_join);
      }
    }

    // 3. join the base domain
    base_abstract_domain_t out_base_dom(is_join ? l_base_dom | r_base_dom
                                                : l_base_dom || r_base_dom);

    // 4. join the regs domain
    uf_register_abstract_domain_t out_uf_regs_dom(
        is_join ? left.m_uf_regs_dom | right.m_uf_regs_dom
                : left.m_uf_regs_dom || right.m_uf_regs_dom);

    obj_flds_id_map_t out_flds_id_map = left.m_flds_id_map;
    refs_base_addrs_map_t out_refs_base_addrs_map =
        left.m_refs_base_addrs_map ? left.m_refs_base_addrs_map
                                   : right.m_refs_base_addrs_map;
    ghost_var_man_t out_ghost_var_man(left.m_ghost_var_man);

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_addrs_dom), std::move(out_uf_regs_dom),
                        std::move(out_obj_info_env), std::move(out_flds_id_map),
                        std::move(out_refs_base_addrs_map),
                        std::move(out_ghost_var_man));
    return res;
  }

  // meet / narrowing two abstract states, return a new meet / narrowing state
  object_domain_t meet_or_narrowing(const object_domain_t &left,
                                    const object_domain_t &right,
                                    const bool is_meet) const {

    // For each abstract object, generally, we have four cases at the meet:
    // The following example only shows properties related to
    // just one abstract object:

    // Case 1: not singleton, odi map contains a field domain for that object
    // Solution: we meet the base dom and odi map
    // E.g.
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = { V_1 = 3, V_2 <= 5 }
    //  s1 meet s2 : Base in s1 meet Base in s2,
    //            Object1 (V_1, V_2) = meet s1 with s2

    // Case 2: singleton object in base domain
    //         odi map does not contain a field domain for that object
    // Solution: we meet the base dom and odi map
    // E.g.
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { V_1 = 2, V_2 = 4 }, Object (V_1, V_2) = empty
    //  s1 meet s2:
    //    Base = { bot }, Object (V_1, V_2) = empty

    // Case 3:
    // one state: singleton object in base domain,
    //  odi map does not contain a field domain for that object
    // another: not singleton, odi map contains a field domain for that object

    // Solution: meet non-singleton in odi map with singleton in base domain
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
    //    b = Base in s1 meet { V_1 <= V_2, V_2 <= 5 }
    //    Base = b meet s2 in Base
    //    Object (V_1, V_2) = empty

    // Other cases:
    // one state: exists such an abstract object
    // another: does not exist that abstract object
    // Example1:
    //  s1: Base = { V_1 = 3, V_2 = 6 }, Object (V_1, V_2) = empty
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 meet s2 :
    //    Base = Base^s1 meet Base^s2, Object (V_1, V_2) = empty

    // Example2:
    //  s1: Base = { ... }, Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }
    //  s2: Base = { ... }, Object (V_1, V_2) = empty
    //  s1 meet s2 :
    //   Base = Base^s1 meet Base^s2,
    //   Object (V_1, V_2) = { V_1 <= V_2, V_2 <= 5 }

    // 1. meet the parts of each state that can use default meet operation
    obj_info_env_t out_obj_info_env(
        is_meet ? left.m_obj_info_env & right.m_obj_info_env
                : left.m_obj_info_env && right.m_obj_info_env);
    address_abstract_domain_t out_addrs_dom(
        is_meet ? left.m_addrs_dom & right.m_addrs_dom
                : left.m_addrs_dom && right.m_addrs_dom);

    // m_flds_id_map and m_refs_base_addrs_map do not require meet,
    // they share accross all states. No need to update
    assert(m_flds_id_map == right.m_flds_id_map);

    // 2. meet the odi map
    base_abstract_domain_t l_base_dom = left.m_base_dom;
    base_abstract_domain_t r_base_dom = right.m_base_dom;
    odi_map_t out_odi_map;
    for (auto kv : left.m_obj_info_env) {
      const obj_id_t &id = kv.first;
      const object_info_t &l_obj_info = kv.second;
      const small_range &l_num_refs = l_obj_info.refcount_val();
      auto r_obj_info_ref = right.m_obj_info_env.find(id);
      if (!r_obj_info_ref) { // only on left
        const odi_domain_product_t *l_prod_ref = left.m_odi_map.find(id);
        if (l_prod_ref) {
          out_odi_map.set(id, (*l_prod_ref));
        }
        continue;
      }
      const object_info_t &r_obj_info = (*r_obj_info_ref);
      const small_range &r_num_refs = r_obj_info.refcount_val();
      const odi_domain_product_t *l_prod_ref = left.m_odi_map.find(id);
      const odi_domain_product_t *r_prod_ref = right.m_odi_map.find(id);
      // Step 1. Join the odi map if the object is non-singleton in both states
      if (l_num_refs == small_range::oneOrMore() &&
          r_num_refs == small_range::oneOrMore()) {
        const boolean_value l_used = l_obj_info.cacheused_val();
        const boolean_value r_used = r_obj_info.cacheused_val();
        odi_domain_product_t l_prod =
            l_prod_ref ? (*l_prod_ref) : odi_domain_product_t();
        odi_domain_product_t r_prod =
            r_prod_ref ? (*r_prod_ref) : odi_domain_product_t();
        if (l_used.is_true() && r_used.is_true()) {
          if (mru_refer_same_object(id, left.m_addrs_dom, right.m_addrs_dom) ==
              false) {
            commit_cache_if_dirty(l_base_dom, left.m_uf_regs_dom, l_prod,
                                  &l_obj_info, id);
            commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                  &r_obj_info, id);
          }
        } else if (l_used.is_true()) {
          commit_cache_if_dirty(l_base_dom, left.m_uf_regs_dom, l_prod,
                                &l_obj_info, id);
        } else if (r_used.is_true()) {
          commit_cache_if_dirty(r_base_dom, right.m_uf_regs_dom, r_prod,
                                &r_obj_info, id);
        }
        // join the summary, cache and flds part
        odi_domain_product_t o_prod =
            is_meet ? l_prod & r_prod : l_prod && r_prod;
        out_odi_map.set(id, o_prod);
      }
      if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
        // left state: singleton
        // right state: non-singleton
        meet_or_narrow_non_singleton_with_singleton(id, l_base_dom, right,
                                                    r_base_dom, is_meet);
        out_odi_map -= id; // remove non singleton in odi map since meet in odi
                           // map will keep that
      } else if (r_num_refs.is_one() &&
                 l_num_refs == small_range::oneOrMore()) {
        meet_or_narrow_non_singleton_with_singleton(id, r_base_dom, left,
                                                    l_base_dom, is_meet);
        out_odi_map -= id; // remove non singleton in odi map since meet in odi
                           // map will keep that
      }
    }

    // add missing odi maps from the right state
    for (auto kv : right.m_obj_info_env) {
      const obj_id_t &id = kv.first;
      auto l_obj_info_ref = left.m_obj_info_env.find(id);
      if (!l_obj_info_ref) { // only on right
        const odi_domain_product_t *r_prod_ref = right.m_odi_map.find(id);
        if (r_prod_ref) {
          out_odi_map.set(id, (*r_prod_ref));
        }
      }
    }

    // 3. meet the base domain
    base_abstract_domain_t out_base_dom(is_meet ? l_base_dom & r_base_dom
                                                : l_base_dom && l_base_dom);

    // 4. meet the regs domain
    uf_register_abstract_domain_t out_uf_regs_dom(
        is_meet ? left.m_uf_regs_dom & right.m_uf_regs_dom
                : left.m_uf_regs_dom && right.m_uf_regs_dom);

    obj_flds_id_map_t out_flds_id_map = left.m_flds_id_map;

    refs_base_addrs_map_t out_refs_base_addrs_map =
        left.m_refs_base_addrs_map ? left.m_refs_base_addrs_map
                                   : right.m_refs_base_addrs_map;

    ghost_var_man_t out_ghost_var_man(left.m_ghost_var_man);

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_addrs_dom), std::move(out_uf_regs_dom),
                        std::move(out_obj_info_env), std::move(out_flds_id_map),
                        std::move(out_refs_base_addrs_map),
                        std::move(out_ghost_var_man));
    return res;
  }

  // compare two abstract states, return boolean
  bool less_than_eq(const object_domain_t &left,
                    const object_domain_t &right) const {

    if (!(left.m_obj_info_env <= right.m_obj_info_env)) {
      CRAB_LOG("object", crab::outs() << "Result3=0\n";);
      return false;
    }

    bool res = left.m_base_dom <= right.m_base_dom;
    CRAB_LOG("object", crab::outs() << "Result4=" << res << "\n";);

    res &= (left.m_odi_map <= right.m_odi_map);
    CRAB_LOG("object", crab::outs() << "Result5=" << res << "\n";);

    res &= (left.m_addrs_dom <= right.m_addrs_dom);
    CRAB_LOG("object", crab::outs() << "Result6=" << res << "\n";);

    res &= (left.m_uf_regs_dom <= right.m_uf_regs_dom);
    CRAB_LOG("object", crab::outs() << "Result7=" << res << "\n";);
    return res;
  }

  bool is_unknown_region(const variable_t &v) {
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

  /***************** Fields and Object id operations *****************/
  boost::optional<obj_id_t> get_obj_id(const variable_t &rgn) {
    if (!m_flds_id_map) {
      return boost::none;
    }
    auto it = (*m_flds_id_map).find(rgn);
    if (it == (*m_flds_id_map).end()) {
      return boost::none;
    }
    return it->second;
  }

  obj_id_t get_obj_id_or_fail(variable_t rgn) {
    assert(m_flds_id_map);
    auto it = (*m_flds_id_map).find(rgn);
    if (it == (*m_flds_id_map).end()) {
      CRAB_ERROR(domain_name(),
                 "::get_odi, the odi map does not include the region ", rgn,
                 " belonging to an abstract object");
    }
    return it->second;
  }

  void get_obj_flds(const obj_id_t &id, variable_vector_t &obj_flds) const {
    assert(m_flds_id_map);
    for (auto kv : (*m_flds_id_map)) {
      if (kv.second == id) {
        obj_flds.push_back(kv.first);
      }
    }
  }

  void get_obj_ghost_flds(const obj_id_t &id,
                          ghost_variable_vector_t &obj_ghost_flds) const {
    variable_vector_t obj_flds;
    get_obj_flds(id, obj_flds);
    for (auto v : obj_flds) {
      auto gvars_opt = get_gvars(v);
      obj_ghost_flds.push_back(*gvars_opt);
    }
  }

  void get_obj_dom_flds(const obj_id_t &id,
                        base_dom_variable_vector_t &obj_dom_flds) const {
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

  bool is_rgn_obj_id(variable_t rgn, obj_id_t id) { return id == rgn; }

  /***************** Base address operations *****************/
  // get or create a variable representing the base address by given a reference
  // this also works for base address of mru object.
  // For mru object, we use a variable named `<obj_id>_mru_base`.
  variable_t get_or_insert_base_addr(const variable_t &v) {
    // if v is a reference variable,
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
    varname_t base_name =
        v_var_name.get_var_factory().get(v_var_name, str_base_name);
    variable_t v_base(base_name, crab::INT_TYPE, 32);
    (*m_refs_base_addrs_map).insert({v, v_base});
    return v_base;
  }

  variable_t get_base_addr_or_fail(const variable_t &v) const {
    assert(v.get_type().is_reference() || v.get_type().is_region());
    if (!m_refs_base_addrs_map) {
      CRAB_ERROR(domain_name(),
                 "::get_base_addr_or_fail: base address map is not created");
    }
    auto it = (*m_refs_base_addrs_map).find(v);
    if (it == (*m_refs_base_addrs_map).end()) {
      CRAB_ERROR(domain_name(), "::get_base_addr_or_fail: ", v, "not found");
    }
    return it->second;
  }

  // check whether mru object in the left state refers the same as the mru
  // object in the right state
  bool
  mru_refer_same_object(const obj_id_t &id,
                        const address_abstract_domain_t &l_addrs_dom,
                        const address_abstract_domain_t &r_addrs_dom) const {
    variable_t mru_base = get_base_addr_or_fail(id);
    const usymb_t *l_term_ptr = l_addrs_dom.get_term(mru_base);
    const usymb_t *r_term_ptr = r_addrs_dom.get_term(mru_base);
    // TODO: so far we use term operator as a constant uninterpreted symbol.
    // This is different from the constant term in crab's uf term.
    // For now, two caches from two states refer the same MRU object
    // if the base addresses of caches refer the same symbol in the uf domain
    // So the way to check two symbols are same is finding operator values
    // Check term/term_operator.hpp for more details.
    return object_domain_impl::is_two_vars_have_same_term(id, l_term_ptr, id,
                                                          r_term_ptr);
  }

  bool test_two_addrs_equality(const variable_t &x, const variable_t &y) {
    const usymb_t *x_term_ptr =
        m_addrs_dom.get_term(get_or_insert_base_addr(x));
    const usymb_t *y_term_ptr =
        m_addrs_dom.get_term(get_or_insert_base_addr(y));
    return object_domain_impl::is_two_vars_have_same_term(x, x_term_ptr, y,
                                                          y_term_ptr);
  }

  bool test_ref_refer_mru_object(const variable_t &ref, const obj_id_t &id) {
    return test_two_addrs_equality(ref, id);
  }

  /***************** Reg-fld domain operations *****************/
  void reg_dom_forget(const variable_vector_t &variables) {
    base_dom_variable_vector_t gvars_vec;
    for (auto const &v : variables) {
      if (auto v_gvars_opt = get_gvars(v)) {
        (*v_gvars_opt).add(gvars_vec);
      }
    }
    m_uf_regs_dom.forget(gvars_vec);
  }

  void reg_dom_project(const variable_vector_t &variables) {
    base_dom_variable_vector_t gvars_vec;
    for (auto const &v : variables) {
      if (auto v_gvars_opt = get_gvars(v)) {
        (*v_gvars_opt).add(gvars_vec);
      }
    }
    m_uf_regs_dom.project(gvars_vec);
  }

  void reg_dom_rename(const variable_vector_t &old_vars,
                      const variable_vector_t &new_vars) {
    base_dom_variable_vector_t old_ghost_vars, new_ghost_vars;
    for (unsigned i = 0, sz = old_vars.size(); i < sz; ++i) {
      variable_t old_var = old_vars[i];
      if (boost::optional<ghost_variables_t> old_gvars = get_gvars(old_var)) {
        (*old_gvars).add(old_ghost_vars);
        variable_t new_var = new_vars[i];
        ghost_variables_t new_gvars = get_or_insert_gvars(new_var);
        new_gvars.add(new_ghost_vars);
      } else {
        continue;
      }
    }
    m_uf_regs_dom.rename(old_ghost_vars, new_ghost_vars);
  }

  void flds_dom_forget(uf_fields_abstract_domain_t &uf_fld_dom,
                       const variable_vector_t &variables) const {
    base_dom_variable_vector_t gvars_vec;
    for (auto const &v : variables) {
      if (auto v_gvars_opt = get_gvars(v)) {
        (*v_gvars_opt).add(gvars_vec);
      }
    }
    uf_fld_dom.forget(gvars_vec);
  }

  void flds_dom_project(uf_fields_abstract_domain_t &uf_fld_dom,
                        const variable_vector_t &variables) const {
    base_dom_variable_vector_t gvars_vec;
    for (auto const &v : variables) {
      if (auto v_gvars_opt = get_gvars(v)) {
        (*v_gvars_opt).add(gvars_vec);
      }
    }
    uf_fld_dom.project(gvars_vec);
  }

  void flds_dom_rename(uf_fields_abstract_domain_t &uf_fld_dom,
                       const variable_vector_t &old_vars,
                       const variable_vector_t &new_vars) {
    base_dom_variable_vector_t old_ghost_vars, new_ghost_vars;
    for (unsigned i = 0, sz = old_vars.size(); i < sz; ++i) {
      variable_t old_var = old_vars[i];
      if (boost::optional<ghost_variables_t> old_gvars = get_gvars(old_var)) {
        (*old_gvars).add(old_ghost_vars);
        variable_t new_var = new_vars[i];
        ghost_variables_t new_gvars = get_or_insert_gvars(new_var);
        new_gvars.add(new_ghost_vars);
      } else {
        continue;
      }
    }
    uf_fld_dom.rename(old_ghost_vars, new_ghost_vars);
  }

  void forget_a_field(const variable_t &v) {
    // for now skip analysis for unknown region
    if (is_unknown_region(v)) {
      return;
    }
    if (auto id_opt = get_obj_id(v)) {
      // retrieve an abstract object info
      auto old_obj_info = m_obj_info_env.at(*id_opt);

      const small_range &num_refs = old_obj_info.refcount_val();

      // Check number of references for an abstract object
      if (num_refs.is_one()) { // singleton object
        m_ghost_var_man.forget(v, m_base_dom);
      } else if (num_refs == small_range::oneOrMore()) { // non-singleton
        const odi_domain_product_t *prod_ref = m_odi_map.find(*id_opt);
        if (prod_ref) {
          odi_domain_product_t out_prod = *prod_ref;
          m_ghost_var_man.forget(v, out_prod.first());
          m_ghost_var_man.forget(v, out_prod.second().first());
          flds_dom_forget(out_prod.second().second(), {v});
          m_odi_map.set(*id_opt, out_prod);
        }
      }
    }
  }

  /***************** ODI map operations *****************/
  void move_singleton_to_odi_map(const obj_id_t &id) {
    crab::CrabStats::count(domain_name() + ".count.move_object");
    crab::ScopedCrabStats __st__(domain_name() + ".move_object");
    base_dom_variable_vector_t flds_vec;
    get_obj_dom_flds(id, flds_vec);

    base_abstract_domain_t singleton_base(m_base_dom);
    singleton_base.project(flds_vec);
    // TODO: performance may matter if base dom contains huge dimensions
    // Might need a projection that only copies flds_vec.

    odi_domain_product_t res_prod; // initial value is top
    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    // In this case, the odi map must not exist that object domain
    assert(!prod_ref);

    // Be careful of the following operation if type of base dom and object dom
    // are different
    res_prod.first() = singleton_base;

    m_odi_map.set(id, res_prod);
    m_base_dom.forget(flds_vec);
  }

  void join_or_widen_singleton_with_non_singleton(
      const obj_id_t &id, const object_domain_t &s_single,
      const object_domain_t &s_non_single,
      base_abstract_domain_t &non_single_base, odi_map_t &odi_map,
      const bool is_join) const {
    // variable_vector_t flds_vec;
    base_dom_variable_vector_t flds_vec;
    s_single.get_obj_dom_flds(id, flds_vec);

    base_abstract_domain_t singleton_base(s_single.m_base_dom);
    singleton_base.project(flds_vec);
    // singleton_base.project(flds_vec);

    odi_domain_product_t res_prod; // initial value is top
    const odi_domain_product_t *prod_ref = s_non_single.m_odi_map.find(id);

    const object_info_t *obj_info_ref = s_non_single.m_obj_info_env.find(id);
    assert(obj_info_ref);
    // Be careful of the following operation if type of base dom and object
    // dom are different
    if (prod_ref) {
      res_prod = (*prod_ref);
      boolean_value c_used = (*obj_info_ref).cacheused_val();
      if (c_used.is_true()) { // commit cache if dirty
        commit_cache_if_dirty(non_single_base, s_non_single.m_uf_regs_dom,
                              res_prod, obj_info_ref, id);
      }
      if (is_join) {
        res_prod.first() |= singleton_base;
      } else {
        res_prod.first() = res_prod.first() || singleton_base;
      }
    }

    odi_map.set(id, res_prod);
  }

  void meet_or_narrow_non_singleton_with_singleton(
      const obj_id_t &id, base_abstract_domain_t &single_base,
      const object_domain_t &s_non_single,
      base_abstract_domain_t &non_single_base, const bool is_meet) const {

    const odi_domain_product_t *prod_ref = s_non_single.m_odi_map.find(id);

    const object_info_t *obj_info_ref = s_non_single.m_obj_info_env.find(id);
    assert(obj_info_ref);

    // Be careful of the following operation if type of base dom and object dom
    // are different
    if (prod_ref) {
      odi_domain_product_t res_prod = (*prod_ref);
      boolean_value c_used = (*obj_info_ref).cacheused_val();
      if (c_used.is_true()) { // commit cache if dirty
        commit_cache_if_dirty(non_single_base, s_non_single.m_uf_regs_dom,
                              res_prod, obj_info_ref, id);
      }
      base_abstract_domain_t &summary = res_prod.first();
      single_base = is_meet ? single_base & summary : single_base && summary;
    }
  }

  /***************** Cache operations *****************/
  void invalidate_cache_if_miss(
      obj_id_t &id, const variable_t &ref, const variable_t &rgn,
      boost::optional<term::term_operator_t> &reg_symb,
      boost::optional<std::pair<term::term_operator_t, term::term_operator_t>>
          &offset_size_symb) {
    variable_t mru_obj_base = get_or_insert_base_addr(id);
    // retrieve an abstract object info
    auto obj_info_ref = m_obj_info_env.find(id);
    assert(obj_info_ref); // The object info must exsits
    boolean_value cache_used = (*obj_info_ref).cacheused_val();
    const odi_domain_product_t *prod_ref = m_odi_map.find(id);
    odi_domain_product_t out_prod =
        prod_ref ? (*prod_ref) : odi_domain_product_t();
    CRAB_LOG("object-entailment",
             crab::outs() << mru_obj_base
                          << " == " << get_or_insert_base_addr(ref) << "?\n"
                          << "Addrs = " << m_addrs_dom << "\n"
                          << "Is mru cached? "
                          << test_two_addrs_equality(id, ref) << "\n";);
    if (cache_used.is_false() || test_two_addrs_equality(id, ref) == false) {
      /* Cache missed condition:
         cache is empty or current ref does not refer to the mru object */
      // Step1: commit cache if the cache is dirty
      commit_cache_if_dirty(m_base_dom, m_uf_regs_dom, out_prod, obj_info_ref,
                            id);
      // Step2: update cache for new MRU object
      update_cache(out_prod.first(), out_prod.second().first());
      // Step3: update address dom and object info
      m_addrs_dom.assign(mru_obj_base, get_or_insert_base_addr(ref));
      m_obj_info_env.set(
          id, object_domain_impl::object_info((*obj_info_ref).refcount_val(),
                                              // Cache is used
                                              boolean_value::get_true(),
                                              // Cache is not dirty
                                              boolean_value::get_false()));
    }
    // Step4: update cache_reg if a field equals to some register
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);
    uf_fields_abstract_domain_t &uf_flds_dom = out_prod.second().second();
    if (reg_symb == boost::none) {
      reg_symb = object_domain_impl::get_or_make_term_symbol(
          uf_flds_dom, rgn_gvars.get_var());
    }
    uf_flds_dom.set(rgn_gvars.get_var(), *reg_symb);
    if (rgn_gvars.has_offset_and_size()) {
      if (offset_size_symb == boost::none) {
        offset_size_symb = std::make_pair(
            object_domain_impl::get_or_make_term_symbol(
                uf_flds_dom, rgn_gvars.get_offset_and_size().get_offset()),
            object_domain_impl::get_or_make_term_symbol(
                uf_flds_dom, rgn_gvars.get_offset_and_size().get_size()));
      }
      uf_flds_dom.set(rgn_gvars.get_offset_and_size().get_offset(),
                      std::get<0>(*offset_size_symb));
      uf_flds_dom.set(rgn_gvars.get_offset_and_size().get_size(),
                      std::get<1>(*offset_size_symb));
    }

    // Step5: update odi map
    m_odi_map.set(id, out_prod);
  }

  // transfer regs-flds' constriants between MRU and base
  void reduce_regs_flds_between_cache_and_base(
      base_abstract_domain_t &base_dom,
      const uf_register_abstract_domain_t &uf_regs_dom,
      odi_domain_product_t &prod,
      const base_dom_variable_vector_t &obj_flds) const {
    // E.g.
    // pre state:
    //      uf_regs_dom = { x == t1 ; z == t2; k = t3 }
    //      base_dom = { 2 <= z; z <= 3; k = 0; }
    //      Object = (
    //                cache_dom = { 0 <= y; y <= 10 },
    //                uf_flds_dom = { y == t1; w == t2 }
    //               )
    // After reduction:
    //      uf_regs_dom = { x == t1 ; z == t2; k = t3 }
    //      base_dom = { k = 0; 2 <= z; z <= 3; 0 <= x; x <= 10; }
    //      Object = (
    //                cache_dom = { 0 <= y; y <= 10; 2 <= w; w <= 3; },
    //                uf_flds_dom = { y == t1; w == t2 }
    //               )
    reduce_regs_from_cache_to_base(base_dom, uf_regs_dom, prod, obj_flds);
    reduce_flds_from_base_to_cache(base_dom, uf_regs_dom, prod, obj_flds);
  }

  void perform_reduction() {
    if (!crab_domain_params_man::get().reduce_everywhere()) {
      // if flag is false, no reduction is performed.
      return;
    }
    CRAB_LOG("object-reduction", crab::outs() << "State Before Reduction:\n"
                                              << *this << "\n";);
    for (auto kv : m_obj_info_env) {
      const obj_id_t &id = kv.first;
      base_dom_variable_vector_t flds_vec;
      get_obj_dom_flds(id, flds_vec);
      const object_info_t &obj_info = kv.second;
      const small_range &num_refs = obj_info.refcount_val();
      if (num_refs == small_range::oneOrMore()) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(id);
        if (prod_ref) {
          reduce_regs_from_cache_to_base(m_base_dom, m_uf_regs_dom, *prod_ref,
                                         flds_vec);
        }
        // m_odi_map.set(id, (out_prod));
      }
    }

    for (auto kv : m_obj_info_env) {
      const obj_id_t &id = kv.first;
      base_dom_variable_vector_t flds_vec;
      get_obj_dom_flds(id, flds_vec);
      const object_info_t &obj_info = kv.second;
      const small_range &num_refs = obj_info.refcount_val();
      if (num_refs == small_range::oneOrMore()) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(id);
        if (prod_ref) {
          odi_domain_product_t out_prod = *prod_ref;
          reduce_flds_from_base_to_cache(m_base_dom, m_uf_regs_dom, out_prod,
                                         flds_vec);
          m_odi_map.set(id, (out_prod));
        }
      }
    }
    CRAB_LOG("object-reduction", crab::outs() << "State After Reduction:\n"
                                              << *this << "\n";);
  }

  void perform_reduction_by_object(const obj_id_t &id,
                                   bool is_from_cache_to_base) {
    base_dom_variable_vector_t flds_vec;
    get_obj_dom_flds(id, flds_vec);
    auto obj_info_ref = m_obj_info_env.find(id);
    if (obj_info_ref) {
      const object_info_t &obj_info = *obj_info_ref;
      const small_range &num_refs = obj_info.refcount_val();
      if (num_refs == small_range::oneOrMore()) {
        const odi_domain_product_t *prod_ref = m_odi_map.find(id);
        if (prod_ref) {
          if (is_from_cache_to_base) {
            // reduce from cache to base
            reduce_regs_from_cache_to_base(m_base_dom, m_uf_regs_dom, *prod_ref,
                                           flds_vec);
          } else {
            // reduce from base to cache
            odi_domain_product_t out_prod = *prod_ref;
            reduce_flds_from_base_to_cache(m_base_dom, m_uf_regs_dom, out_prod,
                                           flds_vec);
            m_odi_map.set(id, (out_prod));
          }
        }
      }
    }
  }

  void build_map(symb_flds_regs_map_t &map,
                 const uf_register_abstract_domain_t &uf_regs_dom,
                 const uf_fields_abstract_domain_t &uf_flds_dom,
                 const base_dom_variable_vector_t &flds) const {
    // it is not allowed to pass a term pointer to different
    // domain value because they use different term managers.
    // Worst case about performance is O(mn) where m and n
    // are number of variables in each domain.
    for (auto &v : flds) {
      const usymb_t *t_ptr = uf_flds_dom.get_term(v);
      auto symb_opt = object_domain_impl::extract_uf_symbol_from_term(t_ptr);
      if (symb_opt == boost::none) {
        continue;
      }
      auto it = map.find(*symb_opt);
      if (it == map.end()) {
        for (auto it_regs = uf_regs_dom.begin(); it_regs != uf_regs_dom.end();
             ++it_regs) {
          if (object_domain_impl::is_two_vars_have_same_term(
                  v, t_ptr, it_regs->first, it_regs->second)) {
            if (it != map.end()) {
              std::get<1>(it->second).push_back(it_regs->first);
            } else {
              map.insert({*symb_opt, {{}, {it_regs->first}}});
            }
          }
        }
      }
      auto it_2 = map.find(*symb_opt);
      if (it_2 != map.end()) {
        std::get<0>(it_2->second).push_back(v);
      }
    }
  }

  // Reduce equalities between registers and fields from base to cache
  void reduce_regs_from_cache_to_base(
      base_abstract_domain_t &base_dom,
      const uf_register_abstract_domain_t &uf_regs_dom,
      const odi_domain_product_t &prod,
      const base_dom_variable_vector_t &obj_flds) const {

    crab::CrabStats::count(domain_name() + ".count.reduction");
    crab::ScopedCrabStats __st__(domain_name() + ".reduction");
    // The following reduction is performed the followings:
    // a. get equalities constraints from uf_regs_dom and uf domain for fields
    // b. rename a field from cache domain with a register if there is
    // a register equals to the same symbol as that field;
    // if there are multiple fields equal to a field, add equality constriants
    // into cache domain.
    // c. forget all fields in cache domain.
    // d. meet reduced and renamed cache domain with base domain
    // E.g.
    // pre state:
    //      uf_regs_dom = { k == t1; w == t2; v == t4; m == t2 }
    //      base_dom = { 3 <= a }
    //      Object = (
    //                summary_dom = {...},
    //                cache_dom = { x <= y; y == z },
    //                uf_flds_dom = { x == t1 ; y == t2; z = t3 }
    //               )
    // After reduction:
    //      uf_regs_dom = { k == t1; w == t2; v == t4; m == t2 }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { x <= y; y == z },
    //                uf_flds_dom = { x == t1 ; y == t2; z = t3 }
    //               )
    CRAB_LOG("object-reduce", crab::outs()
                                  << "Before Reduction from cache to base:\n"
                                  << "base = " << base_dom << "\n"
                                  << "uf_regs = " << uf_regs_dom << "\n"
                                  << "odi = ";
             odi_write(crab::outs(), prod); crab::outs() << "\n";);
    const uf_fields_abstract_domain_t &uf_flds_dom = prod.second().second();
    const base_abstract_domain_t &cache_dom = prod.second().first();
    base_abstract_domain_t reduced_cache = cache_dom;

    symb_flds_regs_map_t map;
    build_map(map, uf_regs_dom, uf_flds_dom, obj_flds);

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t flds = std::get<0>(it->second);
      const base_dom_variable_vector_t regs = std::get<1>(it->second);
      assert(flds.size() > 0);
      assert(regs.size() > 0);
      const base_dom_variable_t &fld = flds[0];
      if (regs.size() == 1) {
        // reduced_cache += ikos::operator==(fld, regs[0]);
        reduced_cache.rename({fld}, {regs[0]});
      } else { // regs.size() > 1
        auto it_regs = regs.begin();
        const base_dom_variable_t &var = *it_regs;
        while (it_regs != regs.end()) {
          if (it_regs == regs.begin()) {
            reduced_cache.rename({fld}, {var});
          } else {
            reduced_cache += ikos::operator==(var, *it_regs);
          }
          ++it_regs;
        }
      }
    }
    map.clear();
    reduced_cache.forget(obj_flds);
    base_dom = base_dom & reduced_cache;
    CRAB_LOG("object-reduce", crab::outs()
                                  << "After Reduction:\n"
                                  << "base = " << base_dom << "\n"
                                  << "uf_regs = " << uf_regs_dom << "\n"
                                  << "odi = ";
             odi_write(crab::outs(), prod); crab::outs() << "\n";);
  }

  // Precondition: if register(s) equal to some field,
  //  the register(s) are reduced from cache to base
  void reduce_flds_from_base_to_cache(
      const base_abstract_domain_t &base_dom,
      const uf_register_abstract_domain_t &uf_regs_dom,
      odi_domain_product_t &prod,
      const base_dom_variable_vector_t &obj_flds) const {

    crab::CrabStats::count(domain_name() + ".count.reduction");
    crab::ScopedCrabStats __st__(domain_name() + ".reduction");
    // The following reduction is performed the followings:
    // a. get equalities constraints from uf_regs_dom and uf domain for fields
    // b. rename register(s) from base domain with a field if they equal to
    // the same symbol as that field;
    // c. project those regs in base domain.
    // d. meet reduced and renamed base domain with cache domain
    // E.g.
    // pre state:
    //      uf_regs_dom = { k == t1; w == t2; v == t4; m == t2; }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { z = 3 },
    //                uf_flds_dom = { x == t1 ; y == t2; }
    //               )
    // After reduction:
    //      uf_regs_dom = { k == t1; w == t2; v == t4; m == t2; }
    //      base_dom = { 3 <= a; k <= w; w == m; }
    //      Object = (
    //                cache_dom = { z = 3, x <= y },
    //                uf_flds_dom = { x == t1 ; y == t2; }
    //               )
    CRAB_LOG("object-reduce", crab::outs()
                                  << "Before Reduction from base to cache:\n"
                                  << "base = " << base_dom << "\n"
                                  << "uf_regs = " << uf_regs_dom << "\n"
                                  << "odi = ";
             odi_write(crab::outs(), prod); crab::outs() << "\n";);
    const uf_fields_abstract_domain_t &uf_flds_dom = prod.second().second();
    base_abstract_domain_t &cache_dom = prod.second().first();
    variable_vector_t remained_regs;
    variable_vector_t renamed_flds;

    symb_flds_regs_map_t map;
    build_map(map, uf_regs_dom, uf_flds_dom, obj_flds);

    // Note that, the regs-flds equalities are one-to-many
    // A field could equal to multiple registers but
    // multiple fields cannot equal to a single register

    base_abstract_domain_t reduced_base = base_dom;

    for (auto it = map.begin(); it != map.end(); ++it) {
      const base_dom_variable_vector_t flds = std::get<0>(it->second);
      const base_dom_variable_vector_t regs = std::get<1>(it->second);
      assert(flds.size() > 0);
      assert(regs.size() > 0);
      const base_dom_variable_t &reg = regs[0];
      if (flds.size() == 1) {
        // reduced_base += ikos::operator==(reg, flds[0]);
        reduced_base.rename({reg}, {flds[0]});
      } else { // flds.size() > 1
        auto it_flds = flds.begin();
        const base_dom_variable_t &var = *it_flds;
        while (it_flds != flds.end()) {
          if (it_flds == flds.begin()) {
            reduced_base.rename({reg}, {var});
            // reduced_base += ikos::operator==(var, reg);
          } else {
            reduced_base += ikos::operator==(var, *it_flds);
          }
          ++it_flds;
        }
      }
    }
    // reduced_base.project(obj_flds);
    reduced_base.project(obj_flds);
    cache_dom = cache_dom & reduced_base;
    CRAB_LOG("object-reduce", crab::outs()
                                  << "After Reduction:\n"
                                  << "base = " << base_dom << "\n"
                                  << "uf_regs = " << uf_regs_dom << "\n"
                                  << "odi = ";
             odi_write(crab::outs(), prod); crab::outs() << "\n";);
  }

  // commit cache contents into summary if cache is dirty
  // reset cache:
  //  (1) perform reduction by transfering regs-flds' constriants between cache
  //  domain and base domain (2) clear UF field domain (3) clear cache domain
  void
  commit_cache_if_dirty(base_abstract_domain_t &base_dom,
                        const uf_register_abstract_domain_t &uf_regs_dom,
                        odi_domain_product_t &prod,
                        const object_domain_impl::object_info *obj_info_ref,
                        const obj_id_t &id) const {
    assert(obj_info_ref);
    base_dom_variable_vector_t obj_flds;
    get_obj_dom_flds(id, obj_flds);
    reduce_regs_flds_between_cache_and_base(base_dom, uf_regs_dom, prod,
                                            obj_flds);
    if ((*obj_info_ref).cachedirty_val().is_true()) {
      // commit cache if the cache is dirty
      commit_cache(prod.first(), prod.second().first());
    }
    prod.second().first().set_to_top();
    prod.second().second().set_to_top();
  }

  // commit cache contents into summary
  void commit_cache(field_abstract_domain_t &summary,
                    field_abstract_domain_t &cache) const {
    // join summary with cache
    summary |= cache;
  }

  // update cache contents from summary to cache
  void update_cache(field_abstract_domain_t &summary,
                    field_abstract_domain_t &cache) const {
    // copy summary to cache
    cache = summary;
  }

  /***************** Print operations *****************/
  void print_flds_id_map(crab_os &o) const {
    o << "Fields -> id map: ";
    if (!m_flds_id_map) {
      o << "not created \n";
    } else if ((*m_flds_id_map).empty()) {
      o << "empty \n";
    } else {
      std::unordered_map<obj_id_t, variable_vector_t> print_map;
      for (auto kv : (*m_flds_id_map)) {
        variable_t field = kv.first;
        obj_id_t id = kv.second;
        auto it = print_map.find(id);
        if (it == print_map.end()) {
          print_map.insert({id, {field}});
        } else {
          it->second.push_back(field);
        }
      }

      // print map
      for (auto kv : print_map) {
        o << "Object (";
        for (auto &v : kv.second) {
          o << v << " ";
        }
        o << ")";
      }
      o << "\n";
    }
  }

  void odi_write(crab_os &o, const odi_domain_product_t &prod) const {
    o << "(";
    o << "summary: ";
    prod.first().write(o);
    o << ", cache: ";
    prod.second().first().write(o);
    o << ", uf_fields: ";
    prod.second().second().write(o);
    o << ")";
  }

  void object_write(crab_os &o) const { // a special output for object domain
    // not using api from seperate domain
    if (m_odi_map.is_bottom()) {
      o << "Object = _|_";
    } else if (m_odi_map.is_top()) {
      o << "Object = {}";
    } else {
      for (auto it = m_odi_map.begin(); it != m_odi_map.end();) {
        o << "Object( ";
        obj_id_t id = it->first;
        variable_vector_t vars;
        get_obj_flds(id, vars);
        std::sort(vars.begin(), vars.end());
        for (auto &v : vars)
          o << v << " ";
        o << ") = ";
        auto prod = it->second;
        odi_write(o, prod);
        ++it;
        if (it != m_odi_map.end()) {
          o << ";\n";
        }
      }
    }
  }

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
      : m_is_bottom(!is_top), m_ghost_var_man(get_type()) {}

  object_domain(const object_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(o.m_base_dom),
        m_odi_map(o.m_odi_map), m_addrs_dom(o.m_addrs_dom),
        m_uf_regs_dom(o.m_uf_regs_dom), m_obj_info_env(o.m_obj_info_env),
        m_flds_id_map(o.m_flds_id_map),
        m_refs_base_addrs_map(o.m_refs_base_addrs_map),
        m_ghost_var_man(o.m_ghost_var_man) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  object_domain(object_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(std::move(o.m_base_dom)),
        m_odi_map(std::move(o.m_odi_map)),
        m_addrs_dom(std::move(o.m_addrs_dom)),
        m_uf_regs_dom(std::move(o.m_uf_regs_dom)),
        m_obj_info_env(std::move(o.m_obj_info_env)),
        m_flds_id_map(std::move(o.m_flds_id_map)),
        m_refs_base_addrs_map(std::move(o.m_refs_base_addrs_map)),
        m_ghost_var_man(std::move(o.m_ghost_var_man)) {
    crab::CrabStats::count(domain_name() + ".count.move");
    crab::ScopedCrabStats __st__(domain_name() + ".move");
  }

  object_domain_t &operator=(const object_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_base_dom = o.m_base_dom;
      m_odi_map = o.m_odi_map;
      m_addrs_dom = o.m_addrs_dom;
      m_uf_regs_dom = o.m_uf_regs_dom;
      m_obj_info_env = o.m_obj_info_env;
      m_flds_id_map = o.m_flds_id_map;
      m_refs_base_addrs_map = o.m_refs_base_addrs_map;
      m_ghost_var_man = o.m_ghost_var_man;
    }
    return *this;
  }

  object_domain_t &operator=(object_domain_t &&o) {
    crab::CrabStats::count(domain_name() + ".count.move");
    crab::ScopedCrabStats __st__(domain_name() + ".move");
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_base_dom = std::move(o.m_base_dom);
      m_odi_map = std::move(o.m_odi_map);
      m_addrs_dom = std::move(o.m_addrs_dom);
      m_uf_regs_dom = std::move(o.m_uf_regs_dom);
      m_obj_info_env = std::move(o.m_obj_info_env);
      m_flds_id_map = std::move(o.m_flds_id_map);
      m_refs_base_addrs_map = std::move(o.m_refs_base_addrs_map);
      m_ghost_var_man = std::move(o.m_ghost_var_man);
    };
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    crab::CrabStats::count(domain_name() + ".count.is_top");
    crab::ScopedCrabStats __st__(domain_name() + ".is_top");

    bool res = (!is_bottom() && m_base_dom.is_top() &&
                m_obj_info_env.is_top() && m_odi_map.is_top());
    return res;
  }

  bool operator<=(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    CRAB_LOG("object", crab::outs() << "Inclusion test:\n\t" << *this << "\n\t"
                                    << o << "\n";);
    if (is_bottom() || o.is_top()) {
      CRAB_LOG("object", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) {
      CRAB_LOG("object", crab::outs() << "Result2=0\n";);
      return false;
    }

    return less_than_eq(*this, o);
  }

  void operator|=(const object_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) { // this is bot, assign this by o
      *this = o;
      return;
    } else if (o.is_bottom()) { // o is bot, nothing change
      return;
    } else if (is_top() | o.is_top()) { // one is top, set to top
      set_to_top();
      return;
    }

    CRAB_LOG("object", crab::outs()
                           << "Join " << *this << " and " << o << " =\n");

    self_join(o);

    CRAB_LOG("object", crab::outs() << "Result=" << *this << "\n");
  }

  object_domain_t operator|(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) { // this is bot, return o
      return o;
    } else if (o.is_bottom()) { // o is bot, return this
      return *this;
    } else if (is_top() | o.is_top()) { // one is top, set to top
      object_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("object", crab::outs()
                           << "Join " << *this << " and " << o << " =\n");

    object_domain_t res(
        std::move(join_or_widening(*this, o, true /*is join*/)));
    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t operator&(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    }

    CRAB_LOG("object", crab::outs()
                           << "Meet " << *this << " and " << o << " =\n");

    object_domain_t res(
        std::move(meet_or_narrowing(*this, o, true /*is meet*/)));

    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  void operator&=(const object_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      *this = o;
      return;
    }

    CRAB_LOG("object", crab::outs()
                           << "Meet " << *this << " and " << o << " =\n");

    // TODO: improve this by avoiding the copy of the left operand.
    *this = std::move(meet_or_narrowing(*this, o, true /*is meet*/));
    CRAB_LOG("object", crab::outs() << *this << "\n");
  }

  object_domain_t operator||(const object_domain_t &o) const override {
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

    CRAB_LOG("object", crab::outs() << "Widening with threshold " << *this
                                    << " and " << o << " =\n");

    object_domain_t res(
        std::move(join_or_widening(*this, o, false /*is widen*/)));

    CRAB_LOG("object", crab::outs() << res << "\n");
    return res;
  }

  object_domain_t operator&&(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

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

  // Initialize a region
  void region_init(const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_init");
    crab::ScopedCrabStats __st__(domain_name() + ".region_init");

    ERROR_IF_NOT_REGION(rgn, __LINE__);

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    if (get_obj_id(rgn) == boost::none) {
      // if a region does not belong to an object, treat it as an object
      // treat region as a field, as well as an object id
      update_fields_id_map(rgn, rgn);
      // initialize information for an abstract object for object
      m_obj_info_env.set(rgn, object_domain_impl::object_info(
                                  // No references owned by the object
                                  small_range::zero(),
                                  // Cache is not used
                                  boolean_value::get_false(),
                                  // Cache is not dirty
                                  boolean_value::get_false()));
    }
    // Assign ghost variables to rgn for modeling its content
    get_or_insert_gvars(rgn);

    CRAB_LOG("object", crab::outs() << "After region_init(" << rgn
                                    << ")=" << *this << "\n";);
  }

  // Create a new reference ref associated with as within region
  void ref_make(const variable_t &ref, const variable_t &rgn,
                /* size of the allocation in bytes */
                const variable_or_constant_t &size,
                /* identifier for the allocation site */
                const allocation_site &as) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    ERROR_IF_NOT_REGION(rgn, __LINE__);
    ERROR_IF_NOT_REF(ref, __LINE__);

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    // if (is_unknown_region(rgn)) {
    //   return;
    // }

    // REDUCTION: perform reduction
    perform_reduction();

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
      auto obj_info_ref = m_obj_info_env.find(*id_opt);
      if (!obj_info_ref) { // object goes to top
        return;
      }
      auto old_obj_info = *obj_info_ref;

      const small_range &num_refs = old_obj_info.refcount_val();

      // Check number of references for an abstract object
      if (num_refs.is_one()) {
        // if the abstract object is a singleton object,
        // now the number of references is increasing,
        // need to move fields' properties into odi map
        move_singleton_to_odi_map(*id_opt);
      }
      // create a fresh symbol to represent a base address
      // at current allocation site
      m_addrs_dom.set(get_or_insert_base_addr(ref),
                      object_domain_impl::make_fresh_symbol(m_addrs_dom));

      m_obj_info_env.set(*id_opt, object_domain_impl::object_info(
                                      old_obj_info.refcount_val().increment(),
                                      old_obj_info.cacheused_val(),
                                      old_obj_info.cachedirty_val()));
    }

    CRAB_LOG("object", crab::outs() << "After ref_make(" << ref << "," << rgn
                                    << ":" << rgn.get_type() << "," << size
                                    << "," << as << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within rgn. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    // Perform a read.
    // ref_load(ref, rgn, res);
    // Definition: assign abstract value stored in rgn to res
    //  rgn: a region variable used in object domain
    //  res: a register variable used in base domain
    // The transformer is sound for ref_load required the followings:
    //  1. The type of region variable is consistent with the type of register
    //  E.g. if region represents integer, the register must be an integer type
    //  2. The reference ref is not evaluated as null pointer.
    //  3. The object domain of the abstract object is either:
    //  3.1 The number of reference = 1: singleton object remained in base dom
    //        Perform assignment on base dom for res := rgn
    //  or the number of reference > 1:
    //  3.2 if object domain is Top, the overapproximate res value as top
    //  3.3 if object domain exists, extract interval for rgn in object domain
    //    assign the interval into base domain
    //
    //  post condition: the odi map not changed, the base domain is updated

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
    // At this point, the requirement for 1 is satisfied

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    const ghost_variables_t &res_gvars = get_or_insert_gvars(res);
    // use ghost variable for field
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_load: reference ", ref,
                                   " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      res_gvars.forget(m_base_dom);
      reg_dom_forget({res});
      return;
    }
    // At this point, the requirement for 2 is satisfied

    if (auto id_opt = get_obj_id(rgn)) {
      // the object id must exist for each region if the
      // region belongs to an abstract object.

      // retrieve an abstract object info
      auto obj_info_ref = m_obj_info_env.find(*id_opt);
      if (!obj_info_ref) { // object goes to top
        return;
      }
      assert(obj_info_ref); // The object info must exsits

      const small_range &num_refs = (*obj_info_ref).refcount_val();

      // In crab IR, the number of references cannot be zero
      //  if ref_load access a not null reference.
      // So zero case should not exist
      assert(!num_refs.is_zero());

      if (num_refs.is_one()) { // singleton object
        // m_base_dom.assign(res, rgn);
        res_gvars.assign(m_base_dom, rgn_gvars);
        if (rgn_gvars.has_offset_and_size() &&
            res_gvars.has_offset_and_size()) {
          rgn_gvars.get_offset_and_size().assign(
              m_base_dom, res_gvars.get_offset_and_size());
        } else if (rgn_gvars.has_offset_and_size()) {
          rgn_gvars.get_offset_and_size().forget(m_base_dom);
        }
        reg_dom_forget({res});
        // the requirment 3.1 is satisfied
      } else { // num_refs > 1, use odi map
        boost::optional<term::term_operator_t> reg_symb = boost::none;
        boost::optional<std::pair<term::term_operator_t, term::term_operator_t>>
            reg_offset_size_symb = boost::none;
        // read from odi map
        const odi_domain_product_t *prod_ref = m_odi_map.find((*id_opt));
        if (prod_ref != nullptr) {
          invalidate_cache_if_miss((*id_opt), ref, rgn, reg_symb,
                                   reg_offset_size_symb);
          if (res.get_type().is_reference()) {
            // res is a reference, assign a fresh symbol into address dom
            // create a fresh symbol to represent its base address
            m_addrs_dom.set(get_or_insert_base_addr(res),
                            object_domain_impl::make_fresh_symbol(m_addrs_dom));
          }
          // assigning register with symbol
          m_uf_regs_dom.set(res_gvars.get_var(), *reg_symb);
          if (res_gvars.has_offset_and_size()) {
            m_uf_regs_dom.set(res_gvars.get_offset_and_size().get_offset(),
                              std::get<0>(*reg_offset_size_symb));
            m_uf_regs_dom.set(res_gvars.get_offset_and_size().get_size(),
                              std::get<1>(*reg_offset_size_symb));
          }
        }
        res_gvars.forget(m_base_dom);
      }
      // post condition meets
    }

    CRAB_LOG("object", crab::outs()
                           << "After " << res << ":="
                           << "ref_load(" << rgn << ":" << rgn.get_type() << ","
                           << ref << ":" << ref.get_type() << ")=" << *this
                           << "\n";);
  }

  // Write the content of val to the address pointed by ref in region.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    // Perform a store.
    // ref_store(ref, rgn, val);
    // Definition: assign abstract value in val to rgn
    //  rgn: a region variable used in object domain
    //  res: a register variable used in base domain
    // The transformer is sound for ref_store required the followings:
    //  1. The type of region variable is consistent with the type of register
    //  E.g. if region represents integer, the register must be an integer type
    //  2. The reference ref is not evaluated as null pointer.
    //  3. The number of references implicitly implies:
    //  a. if is one, it means abstract object is singleton,
    //     perform strong update: update the contents of an object in
    //     abstraction is the same as the concrete
    //  b. if is more than one, it means abstract object is a summarized object,
    //     perform a weak update: update the contents of an object in
    //     abstraction is an over-approximation here: the region is updated for
    //     new val as combined the old properties (through join op)
    //
    // For a. the property of this abstract object is in the base domain:
    //    perform assign in base dom, except if val is not in base domain
    //
    // For b. the object domain of the abstract object is either:
    //     object domain is Top (i.e. the object domain is not found in odi map)
    //  or object domain exists
    //
    //  if object domain is Top, create a new object domain and set into odi map
    //  if object domain exists, create a new one by copying old one
    //
    //  post condition: the odi map is updated, the base domain is not changed

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

    // for now skip analysis for unknown region
    if (is_unknown_region(rgn)) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    // use ghost variable for field
    ghost_variables_t rgn_gvars = get_or_insert_gvars(rgn);
    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_store: reference ",
                                   ref, " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      operator-=(rgn);
      return;
    }

    if (auto id_opt = get_obj_id(rgn)) {
      // the object id must exist for each region if the
      // region belongs to an abstract object.

      // retrieve an abstract object info
      auto obj_info_ref = m_obj_info_env.find(*id_opt);
      if (!obj_info_ref) { // object goes to top
        return;
      }
      assert(obj_info_ref); // The object info must exsits

      const small_range &num_refs = (*obj_info_ref).refcount_val();
      assert(!num_refs.is_zero());

      if (num_refs.is_one()) { // singleton object, perform a strong update
        if (val.is_constant()) {
          m_base_dom.assign(rgn_gvars.get_var(), val.get_constant());
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
          // m_base_dom.assign(rgn_gvars.get_var(), val.get_variable());
        }
      } else { // num_refs > 1, use odi map
        // read from odi map
        const odi_domain_product_t *prod_ref = m_odi_map.find((*id_opt));
        boost::optional<term::term_operator_t> reg_symb = boost::none;
        boost::optional<std::pair<term::term_operator_t, term::term_operator_t>>
            reg_offset_size_symb = boost::none;
        if (val.is_variable()) {
          ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
          reg_symb = object_domain_impl::get_or_make_term_symbol(
              m_uf_regs_dom, val_gvars.get_var());
          if (val_gvars.has_offset_and_size()) {
            reg_offset_size_symb = std::make_pair(
                object_domain_impl::get_or_make_term_symbol(
                    m_uf_regs_dom,
                    val_gvars.get_offset_and_size().get_offset()),
                object_domain_impl::get_or_make_term_symbol(
                    m_uf_regs_dom, val_gvars.get_offset_and_size().get_size()));
          }
        }
        invalidate_cache_if_miss((*id_opt), ref, rgn, reg_symb,
                                 reg_offset_size_symb);
        prod_ref = m_odi_map.find((*id_opt));
        odi_domain_product_t out_prod = (*prod_ref);
        if (val.is_constant()) {
          out_prod.second().first().assign(rgn_gvars.get_var(),
                                           val.get_constant());
          flds_dom_forget(out_prod.second().second(), {rgn});
        } else { // val is a variable (i.e. register)
          ghost_variables_t val_gvars = get_or_insert_gvars(val.get_variable());
          // assigning register with symbol
          m_uf_regs_dom.set(val_gvars.get_var(), *reg_symb);
          if (rgn_gvars.has_offset_and_size() &&
              val_gvars.has_offset_and_size()) {
            m_uf_regs_dom.set(val_gvars.get_offset_and_size().get_offset(),
                              std::get<0>(*reg_offset_size_symb));
            m_uf_regs_dom.set(val_gvars.get_offset_and_size().get_size(),
                              std::get<1>(*reg_offset_size_symb));
          }
          rgn_gvars.forget(out_prod.second().first());
        }
        // update object info
        m_obj_info_env.set(*id_opt, object_domain_impl::object_info(
                                        (*obj_info_ref).refcount_val(),
                                        // Cache is used
                                        boolean_value::get_true(),
                                        // Cache is dirty
                                        boolean_value::get_true()));
        // update odi map
        m_odi_map.set(*id_opt, out_prod);
      }

      // post condition meets
    }

    CRAB_LOG("object", crab::outs()
                           << "After ref_store(" << rgn << ":" << rgn.get_type()
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

    if (is_bottom()) {
      return;
    }

    // for now skip analysis for unknown region or a region not belonging to an
    // abstract object
    // if (is_unknown_region(rgn1) || is_unknown_region(rgn2)) {
    //   return;
    // }

    // REDUCTION: perform reduction
    perform_reduction();

    ghost_variables_t ref1_gvars = get_or_insert_gvars(ref1);
    ghost_variables_t ref2_gvars = get_or_insert_gvars(ref2);

    if (ref1_gvars.has_offset_and_size() && ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().assign(
          m_base_dom, ref1_gvars.get_offset_and_size(), offset);
    } else if (ref2_gvars.has_offset_and_size()) {
      ref2_gvars.get_offset_and_size().forget(m_base_dom);
    }

    if (auto id1_opt = get_obj_id(rgn1)) {
      if (auto id2_opt = get_obj_id(rgn2)) {
        if ((*id1_opt) == (*id2_opt)) { // both regions refer the same object
          // assign equality: ref2 == ref1
          m_addrs_dom.assign(get_or_insert_base_addr(ref2),
                             get_or_insert_base_addr(ref1));
        }
      }
    }

    m_base_dom.assign(ref2_gvars.get_var(), ref1_gvars.get_var() + offset);

    CRAB_LOG("object", crab::outs()
                           << "After (" << rgn2 << "," << ref2
                           << ") := ref_gep(" << rgn1 << "," << ref1 << " + "
                           << offset << ")=" << *this << "\n";);
  }

  // Add constraints between references
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

      // REDUCTION: perform reduction
      perform_reduction();

      auto lin_cst =
          convert_ref_cst_to_linear_cst(ref_cst, ghost_variable_kind::ADDRESS);
      m_base_dom += lin_cst;
      m_is_bottom = m_base_dom.is_bottom();
      if (!m_is_bottom) {
        auto offset_lin_csts =
            convert_ref_cst_to_linear_cst(ref_cst, ghost_variable_kind::OFFSET);
        m_base_dom += offset_lin_csts;
        m_is_bottom = m_base_dom.is_bottom();
      }
    }

    CRAB_LOG("object",
             crab::outs() << "ref_assume(" << ref_cst << ")" << *this << "\n";);
  }

  // Convert a reference to an integer variable
  void ref_to_int(const variable_t &rgn, const variable_t &ref_var,
                  const variable_t &int_var) override {

    crab::CrabStats::count(domain_name() + ".count.ref_to_int");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_to_int");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      // We represent reference as numerical in domain
      m_base_dom.assign(int_var, ref_var);
      m_addrs_dom -= get_or_insert_base_addr(ref_var);
    }
  }

  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    crab::CrabStats::count(domain_name() + ".count.int_to_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".int_to_ref");

    ERROR_IF_NOT_REF(ref_var, __LINE__);
    ERROR_IF_NOT_INT(int_var, __LINE__);

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      ghost_variables_t ref_gvars = get_or_insert_gvars(ref_var);
      ghost_variables_t int_gvars = get_or_insert_gvars(int_var);
      if (ref_gvars.has_offset_and_size()) {
        ref_gvars.get_offset_and_size().forget(m_base_dom);
      }

      ref_gvars.assign(m_base_dom, int_gvars);
      // ref_var is a reference, assign a fresh symbol into address dom
      // create a fresh symbol to represent its base address
      m_addrs_dom.set(get_or_insert_base_addr(ref_var),
                      object_domain_impl::make_fresh_symbol(m_addrs_dom));
    }
  }

  // Make a copy of a region
  // copy from rhs to lhs
  void region_copy(const variable_t &lhs_rgn,
                   const variable_t &rhs_rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_copy");
    crab::ScopedCrabStats __st__(domain_name() + ".region_copy");

    /* These are ensured by well-typed Crab CFGs */
    ERROR_IF_NOT_REGION(lhs_rgn, __LINE__);
    ERROR_IF_NOT_REGION(rhs_rgn, __LINE__);
    if (lhs_rgn.get_type() != rhs_rgn.get_type()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn, ":=", rhs_rgn,
                 " with different types");
    }

    if (is_bottom()) {
      return;
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
    crab::CrabStats::count(domain_name() + ".count.region_cast");
    crab::ScopedCrabStats __st__(domain_name() + ".region_cast");

    CRAB_LOG("object", crab::outs()
                           << "After region_cast(" << src_rgn << ":"
                           << src_rgn.get_type() << "," << dst_rgn << ":"
                           << dst_rgn.get_type() << ")=" << *this << "\n";);
  }

  // Remove a reference ref within region reg
  void ref_free(const variable_t &reg, const variable_t &ref) override {}

  // This default implementation is expensive because it will call the
  // join.
  // DEFAULT_SELECT_REF(object_domain_t)
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    crab::CrabStats::count(domain_name() + ".count.select_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".select_ref");
    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

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
      // ref_var is a reference, assign a fresh symbol into address dom
      // create a fresh symbol to represent its base address
      m_addrs_dom.set(get_or_insert_base_addr(lhs_ref),
                      object_domain_impl::make_fresh_symbol(m_addrs_dom));
    }
  }

  /**************************** Numerical operations *************************/
  // x := y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());
      m_uf_regs_dom -= x;
    }
  }

  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), k);
      m_uf_regs_dom -= x;
    }
  }

  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(),
                       get_or_insert_gvars(z).get_var());
    }
  }

  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply(op, get_or_insert_gvars(x).get_var(),
                       get_or_insert_gvars(y).get_var(), k);
      m_uf_regs_dom -= x;
    }
  }

  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply(op, get_or_insert_gvars(dst).get_var(),
                       get_or_insert_gvars(src).get_var());
      m_uf_regs_dom -= dst;
    }
  }

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    crab::CrabStats::count(domain_name() + ".count.select");
    crab::ScopedCrabStats __st__(domain_name() + ".select");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      auto b_e1 = rename_linear_expr(e1);
      auto b_e2 = rename_linear_expr(e2);
      auto b_cond = rename_linear_cst(cond);
      m_base_dom.select(get_or_insert_gvars(lhs).get_var(), b_cond, b_e1, b_e2);
      m_uf_regs_dom -= lhs;
    }
  }

  // x := e
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      auto b_e = rename_linear_expr(e);
      m_base_dom.assign(get_or_insert_gvars(x).get_var(), b_e);
      m_uf_regs_dom -= x;
    }
  }

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      auto b_csts = rename_linear_cst_sys(csts);
      m_base_dom += b_csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  /********************** Boolean operations **********************/

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      auto b_rhs = rename_linear_cst(rhs);
      m_base_dom.assign_bool_cst(get_or_insert_gvars(lhs).get_var(), b_rhs);
      m_uf_regs_dom -= lhs;
    }
  }

  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {

    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.assign_bool_var(get_or_insert_gvars(lhs).get_var(),
                                 get_or_insert_gvars(rhs).get_var(),
                                 is_not_rhs);
      m_uf_regs_dom -= lhs;
    }
  }

  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.apply_binary_bool(op, get_or_insert_gvars(x).get_var(),
                                   get_or_insert_gvars(y).get_var(),
                                   get_or_insert_gvars(z).get_var());
      m_uf_regs_dom -= x;
    }
  }

  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.assume_bool(get_or_insert_gvars(v).get_var(), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
      m_uf_regs_dom -= v;
    }
  }

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      m_base_dom.select_bool(get_or_insert_gvars(lhs).get_var(),
                             get_or_insert_gvars(cond).get_var(),
                             get_or_insert_gvars(b1).get_var(),
                             get_or_insert_gvars(b2).get_var());
      m_uf_regs_dom -= lhs;
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {

    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      // REDUCTION: perform reduction
      perform_reduction();

      auto rhs_lin_cst =
          convert_ref_cst_to_linear_cst(rhs, ghost_variable_kind::ADDRESS);
      m_base_dom.assign_bool_cst(lhs, rhs_lin_cst);
      m_uf_regs_dom -= lhs;
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
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    CRAB_ERROR(domain_name(), "::expand not implemented");
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const object_domain_t &invariant) override {}

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::to_disjunctive_linear_constraint_system not "
                              "implemented");
  }

  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {

    CRAB_ERROR(domain_name(), "::get_allocation_sites not implemented");
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {

    CRAB_ERROR(domain_name(), "::get_tags not implemented");
    return false;
  }

  // FIXME: The above methods are UNDEFINED METHODS

  // Forget v
  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.-=");
    crab::ScopedCrabStats __st__(domain_name() + ".-=");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    // To forget a variable, it is easier to forget if variable is not a rgn
    // However, if it is a region, we need to perform the followings:
    // 1. the abstract object that region var belongs to is a singleton,
    //    forget it in base dom
    // 2. object is not singleton, forget it in odi map
    if (v.get_type().is_region()) {
      forget_a_field(v);
    } else { // forget a non region variable
      m_ghost_var_man.forget(v, m_base_dom);
      if (v.get_type().is_reference()) {
        m_addrs_dom.operator-=(get_or_insert_base_addr(v));
      }
      reg_dom_forget({v});
    }
  }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        forget_a_field(v);
      } else {
        m_ghost_var_man.forget(v, m_base_dom);
        if (v.get_type().is_reference()) {
          m_addrs_dom.operator-=(get_or_insert_base_addr(v));
        }
        reg_dom_forget({v});
      }
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    // REDUCTION: perform reduction
    perform_reduction();

    variable_vector_t non_odi_vars;
    variable_vector_t ref_vars;
    non_odi_vars.reserve(variables.size());
    ref_vars.reserve(variables.size());
    // The following map keeps fields that need to be remained
    std::unordered_map<variable_t, variable_vector_t> flds_by_id_map;

    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        if (auto id_opt = get_obj_id(v)) {
          // retrieve an abstract object info
          auto old_obj_info = m_obj_info_env.at(*id_opt);

          const small_range &num_refs = old_obj_info.refcount_val();
          if (num_refs.is_one()) { // singleton object
            non_odi_vars.push_back(v);
          } else { // non-singleton object
            auto it = flds_by_id_map.find(*id_opt);
            if (it != flds_by_id_map.end()) {
              it->second.push_back(v);
            } else {
              flds_by_id_map.insert({(*id_opt), {v}});
            }
          }
        }
      } else {
        non_odi_vars.push_back(v);
        if (v.get_type().is_reference()) {
          ref_vars.push_back(get_or_insert_base_addr(v));
        }
      }
    }

    m_ghost_var_man.project(variables, m_base_dom);
    m_addrs_dom.project(ref_vars);
    reg_dom_project(non_odi_vars);

    // projecting fields need to reconstruct odi map
    // keep odi that only remained
    for (auto kv : m_obj_info_env) {
      const obj_id_t &id = kv.first;
      auto it = flds_by_id_map.find(id);
      if (it == flds_by_id_map.end()) {
        // if not found, remove that object in odi map
        m_odi_map.operator-=(id);
        continue;
      }
      const variable_vector_t &flds = it->second;
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (prod_ref) {
        odi_domain_product_t o_prod = *prod_ref;
        // project based on corresponding field(s)
        m_ghost_var_man.project(flds, o_prod.first());
        m_ghost_var_man.project(flds, o_prod.second().first());
        flds_dom_project(o_prod.second().second(), flds);
        m_odi_map.set(id, o_prod);
      }
    }
  }

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

    variable_vector_t from_non_odi_vars;
    variable_vector_t to_non_odi_vars;
    variable_vector_t from_ref_vars;
    variable_vector_t to_ref_vars;
    from_non_odi_vars.reserve(from.size());
    to_non_odi_vars.reserve(from.size());
    from_ref_vars.reserve(from.size());
    to_ref_vars.reserve(from.size());
    // The following map keeps fields that need to be renamed
    std::unordered_map<variable_t,
                       std::pair<variable_vector_t, variable_vector_t>>
        renamed_flds_by_id_map;
    for (int i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &old_v = from[i];
      const variable_t &new_v = to[i];
      if (old_v.get_type() != new_v.get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve same type");
      }
      if (old_v.get_type().is_region()) {
        if (auto id_opt = get_obj_id(old_v)) {
          // retrieve an abstract object info
          auto old_obj_info = m_obj_info_env.at(*id_opt);

          const small_range &num_refs = old_obj_info.refcount_val();
          if (num_refs.is_one()) { // singleton object
            from_non_odi_vars.push_back(old_v);
            to_non_odi_vars.push_back(new_v);
          } else { // non-singleton object
            auto it = renamed_flds_by_id_map.find(*id_opt);
            if (it != renamed_flds_by_id_map.end()) {
              std::get<0>(it->second).push_back(old_v);
              std::get<1>(it->second).push_back(new_v);
            } else {
              renamed_flds_by_id_map.insert({(*id_opt), {{old_v}, {new_v}}});
            }
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
    m_ghost_var_man.rename(from_non_odi_vars, to_non_odi_vars, m_base_dom);
    m_addrs_dom.rename(from_ref_vars, to_ref_vars);
    reg_dom_rename(from_non_odi_vars, to_non_odi_vars);
    // renaming fields need to reconstruct odi map
    for (auto kv : m_obj_info_env) {
      const obj_id_t &id = kv.first;
      auto it = renamed_flds_by_id_map.find(id);
      if (it == renamed_flds_by_id_map.end()) {
        continue;
      }
      const variable_vector_t &from_flds = std::get<0>(it->second);
      const variable_vector_t &to_flds = std::get<1>(it->second);
      const odi_domain_product_t *prod_ref = m_odi_map.find(id);
      if (prod_ref) {
        odi_domain_product_t o_prod = *prod_ref;
        // project based on corresponding field(s)
        m_ghost_var_man.rename(from_flds, to_flds, o_prod.first());
        m_ghost_var_man.rename(from_flds, to_flds, o_prod.second().first());
        flds_dom_rename(o_prod.second().second(), from_flds, to_flds);
        m_odi_map.set(id, o_prod);
      }
    }
  }

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.to_interval");
    crab::ScopedCrabStats __st__(domain_name() + ".to_interval");
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom[v];
    }
  }

  interval_t at(const variable_t &v) const override {
    crab::CrabStats::count(domain_name() + ".count.to_interval");
    crab::ScopedCrabStats __st__(domain_name() + ".to_interval");
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
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
    crab::CrabStats::count(domain_name() + ".count.intrinsic");
    crab::ScopedCrabStats __st__(domain_name() + ".intrinsic");

    //=================================================================//
    //       Special intrinsics supported by the object domain
    //=================================================================//
    // ---DSA region analysis---
    //      This analysis indicates which regions might belong to the
    //      same memory object. The intrinstic is added only if object
    //      has more than one field.
    // TODO: need to support other analysis
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
    auto error_if_not_rgn = [&name](const variable_or_constant_t &x) {
      if (!x.is_variable() || !x.get_type().is_region()) {
        // the input vector should only contains region variables
        CRAB_ERROR("Intrinsics ", name, " parameter ", x,
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

    // REDUCTION: perform reduction
    perform_reduction();

    if (name == "regions_from_memory_object") {
      // pass region varaibles into object field map
      // the intrinsics is only added in clam if the object has more than one
      // region.
      assert(inputs.size() >= 1);
      error_if_not_rgn(inputs[0]);
      obj_id_t obj_id = inputs[0].get_variable();
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        update_fields_id_map(inputs[i].get_variable(), obj_id);
        // Note that, region initialization is before the intrinsic calls
        // Any obj info set up are not for object id will be removed.
        m_obj_info_env -= inputs[i].get_variable();
      }

      // initialize information for an abstract object for object
      m_obj_info_env.set(obj_id, object_domain_impl::object_info(
                                     // No references owned by the object
                                     small_range::zero(),
                                     // Cache is not used
                                     boolean_value::get_false(),
                                     // Cache is not dirty
                                     boolean_value::get_false()));
      // No need to update odi map
      // because top for an abstract object in seperate domain means not exist
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
          // current reference refers mru object
          perform_reduction_by_object(*id_opt, reduce_direction);
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
        if (auto ref_gvars_opt = get_gvars(ref)) {

          if ((*ref_gvars_opt).has_offset_and_size()) {
            // CRAB_LOG("region-domain-is-deref",
            // 	     crab::outs() << "\t" << "Found ghost variables for " << ref
            // << "\n"
            // 	     << "\toffset=" <<
            // (*ref_gvars_opt).get_offset_and_size().get_offset() << "\n"
            // 	     << "\tsize=" <<
            // (*ref_gvars_opt).get_offset_and_size().get_size() << "\n";);
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
    }
  }

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

  std::string domain_name() const override {
    field_abstract_domain_t fld_dom;
    return "ObjectDomain(Base:" + m_base_dom.domain_name() +
           ", Object:" + fld_dom.domain_name() +
           ", Addrs:" + m_addrs_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("object-print", o << "("
                                 << "ObjectInfo=" << m_obj_info_env;
               o << ","
                 << "BaseDom=";
               m_ghost_var_man.write(o, m_base_dom);
               o << ","
                 << "Addrs=" << m_addrs_dom;
               o << ","
                 << "Regs=" << m_uf_regs_dom;
               o << ","; object_write(o); o << ")\n"; return;);
      o << "Base = ";
      m_ghost_var_man.write(o, m_base_dom);
      o << ",\nuf_addrs = ";
      m_addrs_dom.write(o);
      o << ",\nuf_regs = ";
      m_uf_regs_dom.write(o);
      o << ",\n";
      object_write(o);
    }
  }
  /**------------- End domain API definitions -------------------**/

}; // class object_domain

template <typename Params>
struct abstract_domain_traits<object_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab