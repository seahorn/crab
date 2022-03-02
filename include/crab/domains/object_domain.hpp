#pragma once

#include <boost/optional.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/varname_factory.hpp>

#include "region/region_info.hpp"
#include "region/tags.hpp"

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

  // FIXME: the current solution does not support different dom operations
  static_assert(
      std::is_same<base_abstract_domain_t, field_abstract_domain_t>::value,
      "base_abstract_domain_t and field_abstract_domain_t must be the same");

  /**------------------ Begin type definitions ------------------**/

  // Object identifier / object id
  using obj_id_t = variable_t;

  // ODI map
  // Map from an object id to object's abstract value.
  // An object id indicates which kinds of object is.
  // the two domains must maintain the same id before the lattice
  // operations
  using odi_map_t = ikos::separate_domain<
      obj_id_t, field_abstract_domain_t,
      object_domain_impl::object_equal_to<field_abstract_domain_t>>;

  // Object infos
  // Map an object id to finite reduced product domains
  // for inferring object info:
  //  object counting, object initial flag, <unused>
  using obj_info_env_t =
      ikos::separate_domain<obj_id_t, region_domain_impl::region_info>;

  // Object fields to id map
  // Map region variables to object's id
  // This map is constructed during intrinsic call and share in common
  using obj_flds_id_map_t =
      std::shared_ptr<std::unordered_map<variable_t, obj_id_t>>;
  /**------------------ End type definitions ------------------**/

  /**------------------ Begin field definitions ------------------**/
  // a special symbol to represent bottom state of the object domain
  // object domain is bottom if this flag is set or all subdomains are bottom.
  bool m_is_bottom;

  // The abstract state definition:
  // Base domain:
  // Domain for register, singleton object
  base_abstract_domain_t m_base_dom;

  // A map from each region variable to corresponding object domain
  odi_map_t m_odi_map;

  // Map region variables to a tuple of (Count,Init,Type)
  //
  // Count: count how many allocations are owned by an abstract object.
  // Each call to ref_make models a new allocation.
  // This allows us to decide when only if one address per object
  // (i.e., singleton object).
  //
  // Init: unused.
  //
  // Type: unused.
  obj_info_env_t m_obj_info_env;

  // Map region variables to an object id for an abstract object
  // To determine the object id,
  // now is using the first region variable passed by the intrinsic method
  obj_flds_id_map_t m_flds_id_map;
  /**------------------ End field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  // Constructor Definition
  object_domain(base_abstract_domain_t &&base_dom, odi_map_t &&odi_map,
                obj_info_env_t &&obj_info_env, obj_flds_id_map_t &&flds_id_map)
      : m_is_bottom(base_dom.is_bottom()), m_base_dom(base_dom),
        m_odi_map(odi_map), m_obj_info_env(obj_info_env),
        m_flds_id_map(flds_id_map) {}

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

    // 1. join the parts of each state that do not require common renaming
    obj_info_env_t out_obj_info_env(m_obj_info_env | right.m_obj_info_env);

    // 2. join the subdomain for base and objects
    // 2.1 join the base domain
    base_abstract_domain_t out_base_dom =
        base_abstract_domain_t(m_base_dom | right.m_base_dom);

    // 2.2 join the objects' domain over odi map
    odi_map_t out_odi_map(m_odi_map | right.m_odi_map);

    // Handle case 3
    for (auto it = m_obj_info_env.begin(); it != m_obj_info_env.end(); it++) {
      obj_id_t rep = it->first;
      const small_range &left_num_refs = it->second.refcount_val();
      auto right_obj_info_ref = right.m_obj_info_env.find(rep);
      if (!right_obj_info_ref) {
        continue;
      }
      const small_range &right_num_refs = (*right_obj_info_ref).refcount_val();
      if (left_num_refs.is_one() &&
          right_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(rep, m_base_dom, out_odi_map,
                                                   true /* is_join*/);
      } else if (right_num_refs.is_one() &&
                 left_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(
            rep, right.m_base_dom, out_odi_map, true /* is_join*/);
      }
    }

    // m_flds_id_map is fixed after intrinsic calls. No need to update
    assert(m_flds_id_map == right.m_flds_id_map);

    // update current state
    std::swap(m_obj_info_env, out_obj_info_env);
    std::swap(m_base_dom, out_base_dom);
    std::swap(m_odi_map, out_odi_map);
  }

  // join / widening two abstract states, return a new join / widening state
  object_domain_t join_or_widening(const object_domain_t &left,
                                   const object_domain_t &right,
                                   const bool is_join) const {

    // 1. join / widening the parts of each state that do not require common
    // renaming
    obj_info_env_t out_obj_info_env(
        is_join ? left.m_obj_info_env | right.m_obj_info_env
                : left.m_obj_info_env || right.m_obj_info_env);

    // 2. join / widening the subdomain for base and objects
    // 2.1 join / widening the base domain
    base_abstract_domain_t out_base_dom =
        base_abstract_domain_t(is_join ? left.m_base_dom | right.m_base_dom
                                       : left.m_base_dom || right.m_base_dom);

    // 2.2 join / widening the objects' domain over odi map
    // In the real operation,
    // the following operates on two odi maps containg the same keys
    // even if it can performed without such assumption.
    odi_map_t out_odi_map(is_join ? left.m_odi_map | right.m_odi_map
                                  : left.m_odi_map || right.m_odi_map);

    // Handle case 3, see comments on self_join method.
    for (auto it = left.m_obj_info_env.begin(); it != left.m_obj_info_env.end();
         it++) {
      obj_id_t rep = it->first;
      const small_range &left_num_refs = it->second.refcount_val();
      auto right_obj_info_ref = right.m_obj_info_env.find(rep);
      if (!right_obj_info_ref) {
        continue;
      }
      const small_range &right_num_refs = (*right_obj_info_ref).refcount_val();
      if (left_num_refs.is_one() &&
          right_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(rep, left.m_base_dom,
                                                   out_odi_map, is_join);
      } else if (right_num_refs.is_one() &&
                 left_num_refs == small_range::oneOrMore()) {
        join_or_widen_singleton_with_non_singleton(rep, right.m_base_dom,
                                                   out_odi_map, is_join);
      }
    }

    // m_flds_id_map is fixed after intrinsic calls. No need to update
    assert(left.m_flds_id_map == right.m_flds_id_map);

    obj_flds_id_map_t out_flds_id_map = left.m_flds_id_map;

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_obj_info_env),
                        std::move(out_flds_id_map));
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

    // 1. meet / narrowing the parts of each state that do not require common
    // renaming
    obj_info_env_t out_obj_info_env(
        is_meet ? left.m_obj_info_env & right.m_obj_info_env
                : left.m_obj_info_env && right.m_obj_info_env);
    // 2. meet / narrowing the subdomain for base and objects
    // 2.1 meet / narrowing the base subdomain
    base_abstract_domain_t out_base_dom =
        base_abstract_domain_t(is_meet ? left.m_base_dom & right.m_base_dom
                                       : left.m_base_dom && right.m_base_dom);

    // 2.2 meet / narrowing the objects' domain over odi map
    // In the real operation,
    // the following operates on two odi maps containg the same keys
    // even if it can performed without such assumption.
    odi_map_t out_odi_map(is_meet ? left.m_odi_map & right.m_odi_map
                                  : left.m_odi_map && right.m_odi_map);

    // Handle case 3
    for (auto it = left.m_obj_info_env.begin(); it != left.m_obj_info_env.end();
         it++) {
      obj_id_t rep = it->first;
      const small_range &left_num_refs = it->second.refcount_val();
      auto right_obj_info_ref = right.m_obj_info_env.find(rep);
      if (!right_obj_info_ref) {
        continue;
      }
      const small_range &right_num_refs = (*right_obj_info_ref).refcount_val();
      if (left_num_refs.is_one() &&
          right_num_refs == small_range::oneOrMore()) {
        meet_or_narrow_non_singleton_with_singleton(rep, out_base_dom,
                                                    right.m_odi_map, is_meet);
        out_odi_map -= rep; // remove non singleton in odi map since meet in odi
                            // map will keep that
      } else if (right_num_refs.is_one() &&
                 left_num_refs == small_range::oneOrMore()) {
        meet_or_narrow_non_singleton_with_singleton(rep, out_base_dom,
                                                    left.m_odi_map, is_meet);
        out_odi_map -= rep; // remove non singleton in odi map since meet in odi
                            // map will keep that
      }
    }

    // m_flds_id_map is fixed after intrinsic calls. No need to update
    assert(left.m_flds_id_map == right.m_flds_id_map);

    obj_flds_id_map_t out_flds_id_map = left.m_flds_id_map;

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_obj_info_env),
                        std::move(out_flds_id_map));
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

    res = left.m_odi_map <= right.m_odi_map;

    CRAB_LOG("object", crab::outs() << "Result5=" << res << "\n";);
    return res;
  }

  linear_constraint_t
  convert_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst) {
    if (ref_cst.is_tautology()) {
      return linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        variable_t x = ref_cst.lhs();
        if (ref_cst.is_equality()) {
          return linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
        variable_t x = ref_cst.lhs();
        variable_t y = ref_cst.rhs();
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return linear_constraint_t(x == y + offset);
        } else if (ref_cst.is_disequality()) {
          return linear_constraint_t(x != y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return linear_constraint_t(x <= y + offset);
        } else if (ref_cst.is_less_than()) {
          return linear_constraint_t(x < y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return linear_constraint_t(x >= y + offset);
        } else if (ref_cst.is_greater_than()) {
          return linear_constraint_t(x > y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
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

  boost::optional<obj_id_t> get_obj_id(variable_t rgn) {
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

  void get_obj_flds(obj_id_t id, variable_vector_t &obj_flds) const {
    assert(m_flds_id_map);
    for (auto kv : (*m_flds_id_map)) {
      if (kv.second == id) {
        obj_flds.push_back(kv.first);
      }
    }
  }

  bool is_rgn_obj_id(variable_t rgn, obj_id_t id) { return id == rgn; }

  void object_write(crab_os &o) const { // a special output for object domain
    // not using api from seperate domain
    if (m_odi_map.is_bottom()) {
      o << "Object = _|_";
    } else if (m_odi_map.is_top()) {
      o << "Object = {}";
    } else {
      assert(m_flds_id_map != nullptr);
      for (auto it = m_odi_map.begin(); it != m_odi_map.end();) {
        o << "Object( ";
        obj_id_t id = it->first;
        variable_vector_t vars;
        get_obj_flds(id, vars);
        std::sort(vars.begin(), vars.end());
        for (auto &v : vars)
          o << v << " ";
        o << ") = ";
        auto dom = it->second;
        dom.write(o);
        ++it;
        if (it != m_odi_map.end()) {
          o << "; ";
        }
        o << "\n";
      }
    }
  }

  void move_singleton_to_odi_map(variable_t rep) {
    variable_vector_t flds_vec;
    get_obj_flds(rep, flds_vec);

    base_abstract_domain_t tmp_base(m_base_dom);
    tmp_base.project(flds_vec);

    field_abstract_domain_t res_obj_dom;
    const field_abstract_domain_t *obj_dom_ref = m_odi_map.find(rep);
    // In this case, the odi map must not exist that object domain
    assert(!obj_dom_ref);

    // Be careful of the following operation if type of base dom and object dom
    // are different
    res_obj_dom = tmp_base;

    m_odi_map.set(rep, res_obj_dom);
    m_base_dom.forget(flds_vec);
  }

  void join_or_widen_singleton_with_non_singleton(
                        variable_t rep, const base_abstract_domain_t &base_dom,
                        odi_map_t &odi_map, const bool is_join) const {
    variable_vector_t flds_vec;
    get_obj_flds(rep, flds_vec);

    base_abstract_domain_t tmp_base(base_dom);
    tmp_base.project(flds_vec);

    field_abstract_domain_t res_obj_dom;
    const field_abstract_domain_t *obj_dom_ref = odi_map.find(rep);

    // Be careful of the following operation if type of base dom and object dom
    // are different
    if (obj_dom_ref) {
      res_obj_dom =
          is_join ? tmp_base | *obj_dom_ref : tmp_base || *obj_dom_ref;
    }

    odi_map.set(rep, res_obj_dom);
  }

  void meet_or_narrow_non_singleton_with_singleton(
                        variable_t rep, base_abstract_domain_t &base_dom,
                        const odi_map_t &odi_map, const bool is_meet) const {

    const field_abstract_domain_t *obj_dom_ref = odi_map.find(rep);

    // Be careful of the following operation if type of base dom and object dom
    // are different
    if (obj_dom_ref) {
      base_dom = is_meet ? base_dom & *obj_dom_ref : base_dom && *obj_dom_ref;
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

  object_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  object_domain(const object_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(o.m_base_dom),
        m_odi_map(o.m_odi_map), m_obj_info_env(o.m_obj_info_env),
        m_flds_id_map(o.m_flds_id_map) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  object_domain(object_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(std::move(o.m_base_dom)),
        m_odi_map(std::move(o.m_odi_map)),
        m_obj_info_env(std::move(o.m_obj_info_env)),
        m_flds_id_map(std::move(o.m_flds_id_map)) {}

  object_domain_t &operator=(const object_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_base_dom = o.m_base_dom;
      m_odi_map = o.m_odi_map;
      m_obj_info_env = o.m_obj_info_env;
      m_flds_id_map = o.m_flds_id_map;
    }
    return *this;
  }

  object_domain_t &operator=(object_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_base_dom = std::move(o.m_base_dom);
      m_odi_map = std::move(o.m_odi_map);
      m_obj_info_env = std::move(o.m_obj_info_env);
      m_flds_id_map = std::move(o.m_flds_id_map);
    };
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
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

    if (is_bottom() || o.is_top()) { // bot | o = top
      *this = o;
      return;
    } else if (o.is_bottom() || is_top()) { // this = top | o = bot
      return;
    }

    CRAB_LOG("object", crab::outs()
                           << "Join " << *this << " and " << o << " =\n");

    self_join(o);

    CRAB_LOG("object", crab::outs() << *this << "\n");
  }

  object_domain_t operator|(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top()) { // bot | o = top
      return o;
    } else if (o.is_bottom() || is_top()) { // this = top | o = bot
      return *this;
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

  object_domain_t widening_thresholds(
      const object_domain_t &o,
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

  // x := e operates on register dom
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      m_base_dom.assign(x, e);
    }
  }

  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      m_base_dom += csts;
    }
  }

  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.assign_bool_cst(lhs, rhs);
    }
  }

  // if region variable is not in odi_map_t, create a new odi for it. That
  // means the region variable belongs to an object that only contains a single
  // region.
  /***************** Regions and reference operations *****************/
  // Initialize a region
  void region_init(const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_init");
    crab::ScopedCrabStats __st__(domain_name() + ".region_init");

    ERROR_IF_NOT_REGION(rgn, __LINE__);

    if (is_bottom()) {
      return;
    }

    CRAB_LOG("object", crab::outs()
                           << "After region_init(" << rgn << ":"
                           << rgn.get_type() << ")=" << *this << "\n";);
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
    if (is_unknown_region(rgn)) {
      return;
    }

    if (auto id_opt = get_obj_id(rgn)) {
      // After instrisic calls, the object id must exist for each region if the
      // region belongs to an abstract object

      // retrieve an abstract object info
      auto old_obj_info = m_obj_info_env.at(*id_opt);

      const small_range &num_refs = old_obj_info.refcount_val();

      // Check number of references for an abstract object
      if (num_refs.is_one()) {
        // if the abstract object is a singleton object,
        // now the number of references is increasing,
        // need to move fields' properties into odi map
        const field_abstract_domain_t *obj_dom_ref = m_odi_map.find(*id_opt);
        move_singleton_to_odi_map(*id_opt);
      }

      m_obj_info_env.set(*id_opt,
                         region_domain_impl::region_info(
                             old_obj_info.refcount_val().increment(),
                             old_obj_info.init_val(), old_obj_info.type_val()));
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

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object", CRAB_WARN(domain_name(), "::ref_load: reference ", ref,
                                   " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      m_base_dom -= res;
      return;
    }
    // At this point, the requirement for 2 is satisfied

    if (auto id_opt = get_obj_id(rgn)) {
      // After instrisic calls, the object id must exist for each region if the
      // region belongs to an abstract object

      // retrieve an abstract object info
      auto obj_info_ref = m_obj_info_env.find(*id_opt);
      assert(obj_info_ref); // The object info must exsits

      const small_range &num_refs = (*obj_info_ref).refcount_val();

      // In crab IR, the number of references cannot be zero
      //  if ref_load access a not null reference.
      // So zero case should not exist
      assert(!num_refs.is_zero());

      if (num_refs.is_one()) { // singleton object
        m_base_dom.assign(res, rgn);
        // the requirment 3.1 is satisfied
      } else { // num_refs > 1, use odi map
        // read from odi map
        const field_abstract_domain_t *obj_dom_ref = m_odi_map.find(*id_opt);
        if (!obj_dom_ref) {  // not found in odi map, means object domain is top
          m_base_dom -= res; // forget res means res == top
          // the requirment 3.2 is satisfied
        } else {
          interval_t ival = (*obj_dom_ref).at(rgn);
          object_domain_impl::assign_interval(m_base_dom, res,
                                              ival); // res := ival
          // the requirment 3.3 is satisfied
        }
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

    if (auto id_opt = get_obj_id(rgn)) {
      // After instrisic calls, the object id must exist for each region if the
      // region belongs to an abstract object

      // retrieve an abstract object info
      auto(*obj_info_ref) = m_obj_info_env.find(*id_opt);
      assert((*obj_info_ref)); // The object info must exsits

      const small_range &num_refs = (*obj_info_ref).refcount_val();
      assert(!num_refs.is_zero());

      if (num_refs.is_one()) { // singleton object, perform a strong update
        if (val.is_constant()) {
          m_base_dom.assign(rgn, val.get_constant());
        } else {
          m_base_dom.assign(rgn, val.get_variable());
        }
      } else { // num_refs > 1, use odi map
        // read from odi map
        auto obj_dom_ptr = m_odi_map.find(*id_opt);
        // if not found, it means that domain is top, no need to update

        if (obj_dom_ptr) {
          field_abstract_domain_t res_obj_dom =
              field_abstract_domain_t((*obj_dom_ptr));
          if (val.is_constant()) {
            res_obj_dom.assign(rgn, val.get_constant());
          } else { // val is a variable
            interval_t ival = m_base_dom[val.get_variable()];
            object_domain_impl::assign_interval(res_obj_dom, rgn, ival);
          }

          res_obj_dom |= *obj_dom_ptr; // perform a weak update
          m_odi_map.set(*id_opt, res_obj_dom);
        }
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
    if (is_unknown_region(rgn1) || is_unknown_region(rgn2)) {
      return;
    }

    m_base_dom.assign(ref2, ref1 + offset);

    CRAB_LOG("object", crab::outs()
                           << "After (" << rgn2 << "," << ref2
                           << ") := ref_gep(" << rgn1 << "," << ref1 << " + "
                           << offset << ")=" << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, z);
    }
  }

  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, k);
    }
  }

  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, z);
    }
  }

  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, k);
    }
  }

  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, dst, src);
    }
  }

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {

    if (!is_bottom()) {
      m_base_dom.select(lhs, cond, e1, e2);
    }
  }

  /********************** Boolean operations **********************/

  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {

    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.assign_bool_var(lhs, rhs, is_not_rhs);
    }
  }

  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.apply_binary_bool(op, x, y, z);
    }
  }

  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_bottom()) {
      m_base_dom.assume_bool(v, is_negated);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {

    if (!is_bottom()) {
      m_base_dom.select_bool(lhs, cond, b1, b2);
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {

    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      auto rhs_lin_cst = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.assign_bool_cst(lhs, rhs_lin_cst);
    }
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

      auto lin_cst = convert_ref_cst_to_linear_cst(ref_cst);
      m_base_dom += lin_cst;
      m_is_bottom = m_base_dom.is_bottom();
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
      // We represent reference as numerical in domain
      m_base_dom.assign(int_var, ref_var);
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
      m_base_dom.assign(ref_var, int_var);
    }
  }

  // This default implementation is expensive because it will call the
  // join.
  DEFAULT_SELECT_REF(object_domain_t)

  // FIXME: The followings are UNDEFINED METHODS
  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
  }

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
  // Make a copy of a region
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {}
  // Cast between regions of different types
  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {}
  // Remove a reference ref within region reg
  void ref_free(const variable_t &reg, const variable_t &ref) override {}
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

  /********************** Backward array operations **********************/
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const object_domain_t &invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const object_domain_t &invariant) override {}
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const object_domain_t &invariant) override {}
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const object_domain_t &invariant) override {}
  void backward_array_assign(const variable_t &a, const variable_t &b,
                             const object_domain_t &invariant) override {}

  /********************** Miscellaneous operations **********************/

  // Normalize the abstract domain if such notion exists.
  void normalize() override {}

  // Reduce the size of the abstract domain representation.
  void minimize() override {}

  // Make a new copy of var without relating var with new_var
  void expand(const variable_t &var, const variable_t &new_var) override {}

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const object_domain_t &invariant) override {}

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

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(domain_name(), "::to_disjunctive_linear_constraint_system not "
                              "implemented");
  }

  // FIXME: The above methods are UNDEFINED METHODS

  // Forget v
  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
  }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_bottom() || is_top()) {
      return;
    }
  }

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom[v];
    }
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom.at(v);
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
    auto error_if_not_rgn = [&name](const variable_or_constant_t &x) {
      if (!x.is_variable() || !x.get_type().is_region()) {
        // the input vector should only contains region variables
        CRAB_ERROR("Intrinsics ", name, " parameter ", x,
                   " should be a region");
      }
    };

    if (is_bottom()) {
      return;
    }

    if (name == "regions_from_memory_object") {
      // pass region varaibles into object field map
      // the intrinsics is only added in clam if the object has more than one
      // region.
      assert(inputs.size() >= 1);
      error_if_not_rgn(inputs[0]);
      obj_id_t obj_id = inputs[0].get_variable();
      if (!m_flds_id_map) { // create a object fields - id map if it is null
        m_flds_id_map =
            std::make_shared<std::unordered_map<variable_t, variable_t>>();
      }
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        (*m_flds_id_map).insert({inputs[i].get_variable(), obj_id});
      }

      // initialize information for an abstract object
      m_obj_info_env.set(obj_id, region_domain_impl::region_info(
                                     // No references owned by the object
                                     small_range::zero(),
                                     // unused
                                     boolean_value::get_false(),
                                     // unused
                                     obj_id.get_type()));

      // No need to update odi map
      // because top for an abstract object in seperate domain is not exists
    }
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    if (is_bottom()) {
      return boolean_value::bottom();
    }

    if (!ref.get_type().is_reference()) {
      return boolean_value::get_false();
    }

    interval_t ival = m_base_dom[ref];
    number_t zero(0);

    if (!(interval_t(zero) <= ival)) {
      return boolean_value::get_false();
    }

    boost::optional<number_t> x = ival.lb().number();
    boost::optional<number_t> y = ival.ub().number();
    if (x && y && *x == zero && *y == zero) { // ref == 0
      return boolean_value::get_true();
    }

    return boolean_value::top();
  }

  std::string domain_name() const override {
    field_abstract_domain_t fld_dom;
    return "ObjectDomain(Regs:" + m_base_dom.domain_name() +
           ", Object:" + fld_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      o << "Base = ";
      m_base_dom.write(o);
      o << ", ";
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