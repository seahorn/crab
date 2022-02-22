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

#include "object/region_object.hpp"
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

  /**------------------ Begin type definitions ------------------**/

  // ODI map
  // Map from a representative to object's abstract value.
  // A representative region indicates which kinds of object is.
  // the two domains must maintain the same representative before the lattice
  // operations
  using odi_map_t = ikos::separate_domain<
      variable_t, object_domain_impl::region_object<field_abstract_domain_t>>;

  // Object infos
  // Map representative region (i.e. object) to finite reduced product domains
  // for inferring object info:
  //  object counting, object initial flag, <unused>
  using obj_info_env_t =
      ikos::separate_domain<variable_t, region_domain_impl::region_info>;

  // Region representative map
  // Map regions to object's region representative
  using flds_rep_map_t = std::unordered_map<variable_t, variable_t>;
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

  // Map region variables to a representative region for an abstract object
  // To determine the representative region,
  // now is using the first region variable passed by the intrinsic method
  flds_rep_map_t m_flds_rep_map;
  /**------------------ End field definitions ------------------**/

  /**------------------ Begin helper method definitions ------------------**/
  // Constructor Definition
  object_domain(base_abstract_domain_t &&base_dom, odi_map_t &&odi_map,
                obj_info_env_t &&obj_info_env, flds_rep_map_t &&flds_rep_map)
      : m_is_bottom(base_dom.is_bottom()), m_base_dom(base_dom),
        m_odi_map(odi_map), m_obj_info_env(obj_info_env),
        m_flds_rep_map(flds_rep_map) {}

  // Self join
  // Perform *this = join(*this, right)
  void self_join(const object_domain_t &right) {

    // 1. join the parts of each state that do not require common renaming
    obj_info_env_t out_obj_info_env(m_obj_info_env | right.m_obj_info_env);

    // 2. join the subdomain for base and objects
    // 2.1 join the base domain
    base_abstract_domain_t out_base_dom =
        base_abstract_domain_t(m_base_dom | right.m_base_dom);

    // 2.2 join the objects' domain over odi map
    // In the real join operation,
    // the following join operates on two odi maps containg the same keys
    // even if it can performed without such assumption.
    odi_map_t out_odi_map(m_odi_map | right.m_odi_map);

    // The following map update is not necessary because m_flds_rep_map
    // is fixed after intrinsic calls.
    // Here, in case to keep operation sound.
    flds_rep_map_t out_flds_rep_map;
    for (auto &kv : m_flds_rep_map) {
      const variable_t &v = kv.first;
      auto it = right.m_flds_rep_map.find(v);
      if (it != right.m_flds_rep_map.end()) {
        assert(it->second == kv.second);
        out_flds_rep_map.insert({v, kv.second});
      }
    }

    // update current state
    std::swap(m_obj_info_env, out_obj_info_env);
    std::swap(m_base_dom, out_base_dom);
    std::swap(m_odi_map, out_odi_map);
    std::swap(m_flds_rep_map, out_flds_rep_map);

    CRAB_LOG("object-domain", crab::outs() << *this << "\n");
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

    // The following map update is not necessary because m_flds_rep_map
    // is fixed after intrinsic calls.
    // Here, in case to keep operation sound.
    flds_rep_map_t out_flds_rep_map;
    for (auto &kv : left.m_flds_rep_map) {
      const variable_t &v = kv.first;
      auto it = right.m_flds_rep_map.find(v);
      if (it != right.m_flds_rep_map.end()) {
        assert(it->second == kv.second);
        out_flds_rep_map.insert({v, kv.second});
      }
    }

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_obj_info_env),
                        std::move(out_flds_rep_map));
    return res;
  }

  // meet / narrowing two abstract states, return a new meet / narrowing state
  object_domain_t meet_or_narrowing(const object_domain_t &left,
                                    const object_domain_t &right,
                                    const bool is_meet) const {

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

    // The following map update is not necessary because m_flds_rep_map
    // is fixed after intrinsic calls.
    // Here, in case to keep operation sound.
    flds_rep_map_t out_flds_rep_map(left.m_flds_rep_map);
    for (auto &kv : right.m_flds_rep_map) {
      const variable_t &v = kv.first;
      auto it = out_flds_rep_map.find(v);
      if (it == out_flds_rep_map.end()) {
        out_flds_rep_map.insert({v, kv.second});
      }
    }

    object_domain_t res(std::move(out_base_dom), std::move(out_odi_map),
                        std::move(out_obj_info_env),
                        std::move(out_flds_rep_map));
    return res;
  }

  // compare two abstract states, return boolean
  bool less_than_eq(const object_domain_t &left,
                    const object_domain_t &right) const {

    if (!(left.m_obj_info_env <= right.m_obj_info_env)) {
      CRAB_LOG("object-leq", crab::outs() << "Result3=0\n";);
      return false;
    }

    bool res = left.m_base_dom <= right.m_base_dom;
    CRAB_LOG("object-leq", crab::outs() << "Result4=" << res << "\n";);

    res = left.m_odi_map <= right.m_odi_map;

    CRAB_LOG("object-leq", crab::outs() << "Result5=" << res << "\n";);
    return res;
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

  variable_t get_object_rep_by_rgn(variable_t rgn) {
    auto it = m_flds_rep_map.find(rgn);
    if (it == m_flds_rep_map.end()) {
      CRAB_ERROR(
          domain_name(),
          "::get_object_rep_by_rgn, the odi map does not include the region ",
          rgn, " belonging to an abstract object");
    }
    return it->second;
  }

  void get_rgn_vec(variable_t rep, variable_vector_t &rgn_vec) {
    for (auto kv : m_flds_rep_map) {
      if (kv.second == rep) {
        rgn_vec.push_back(kv.second);
      }
    }
  }

  template <class DOM>
  void assign_interval(DOM &dom, const variable_t v, interval_t ival) {
    dom -= v;
    boost::optional<number_t> lb = ival.lb().number();
    boost::optional<number_t> ub = ival.ub().number();
    if (lb) {
      dom += (*lb <= v);
    }
    if (ub) {
      dom += (v <= *ub);
    }
  }

  void move_singleton_to_odi_map(variable_t rep) {
    variable_vector_t rgn_vec;
    base_abstract_domain_t tmp_base(m_base_dom);
    tmp_base.project(rgn_vec);
    // Be careful of the following operation if base dom and object dom are
    // different numerical domain
    auto obj_dom = m_odi_map[rep];
    obj_dom.set_dom(tmp_base);
    m_odi_map.set(rep, obj_dom);
  }

  void object_write(crab_os &o) const { // a special output for object domain
    // not using api from seperate domain
    if (m_odi_map.is_bottom()) {
      o << "_|_";
    } else if (m_odi_map.is_top()) {
      o << "{}";
    } else {
      for (auto it = m_odi_map.begin(); it != m_odi_map.end();) {
        o << "Object( ";
        variable_t rep = it->first;
        variable_vector_t vars;
        for (auto &kv : m_flds_rep_map) {
          if (rep == kv.second) {
            vars.push_back(kv.first);
          }
        }
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
        m_flds_rep_map(o.m_flds_rep_map) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  object_domain(object_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_base_dom(std::move(o.m_base_dom)),
        m_odi_map(std::move(o.m_odi_map)),
        m_obj_info_env(std::move(o.m_obj_info_env)),
        m_flds_rep_map(std::move(o.m_flds_rep_map)) {}

  object_domain_t &operator=(const object_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_base_dom = o.m_base_dom;
      m_odi_map = o.m_odi_map;
      m_obj_info_env = o.m_obj_info_env;
      m_flds_rep_map = o.m_flds_rep_map;
    }
    return *this;
  }

  object_domain_t &operator=(object_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_base_dom = std::move(o.m_base_dom);
      m_odi_map = std::move(o.m_odi_map);
      m_obj_info_env = std::move(o.m_obj_info_env);
      m_flds_rep_map = std::move(o.m_flds_rep_map);
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
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top()) { // bot | o = top
      *this = o;
      return;
    } else if (o.is_bottom() || is_top()) { // this = top | o = bot
      return;
    }

    CRAB_LOG("object", crab::outs() << "Join " << *this << " and " << o << "=");

    self_join(o);
  }

  object_domain_t operator|(const object_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top()) { // bot | o = top
      return o;
    } else if (o.is_bottom() || is_top()) { // this = top | o = bot
      return *this;
    }

    CRAB_LOG("object", crab::outs() << "Join " << *this << " and " << o << "=");

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

    CRAB_LOG("object", crab::outs() << "Meet " << *this << " and " << o << "=");

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
                           << "Widening " << *this << " and " << o << "=");

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
                                    << " and " << o << "=");

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
                           << "Narrowing " << *this << " and " << o << "=");

    object_domain_t res(
        std::move(meet_or_narrowing(*this, o, false /*is narrow*/)));

    crab::outs() << res << "\n";
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

  void region_assign(const variable_t &x, const linear_expression_t &e) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    std::string name = domain_name();
    auto error_if_not_rgn = [&name](const variable_t &x) {
      if (!x.get_type().is_region()) {
        // the lhs should be a region variable
        CRAB_ERROR(name, "::region_assign parameter ", x,
                   " should be a region");
      }
    };
    error_if_not_rgn(x);

    if (!is_bottom()) {
      auto it = m_flds_rep_map.find(x);
      if (it == m_flds_rep_map.end()) {
        CRAB_ERROR(domain_name(),
                   "::region_assign, the odi map does not include the region ",
                   x, " belonging to an abstract object");
      }
      variable_t rep = it->second;
      auto obj_dom = m_odi_map[rep];
      obj_dom.assign(x, e);
      m_odi_map.set(rep, obj_dom);
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

  // add constraint into object domain
  void add_cons_for_objects(const linear_constraint_t &cst) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    // assume cst always assign within an object
    if (!is_bottom()) {
      for (auto v : cst.variables()) {
        if (v.get_type().is_region()) {
          auto it = m_flds_rep_map.find(v);
          if (it == m_flds_rep_map.end()) {
            CRAB_ERROR(domain_name(),
                       "::add_cons_for_objects, the odi map does not include "
                       "the region ",
                       v, " belonging to an object");
          }
          variable_t rep = it->second;
          auto obj_dom = m_odi_map[rep];
          obj_dom += cst;
          m_odi_map.set(rep, obj_dom);
          break;
        }
      }
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

    variable_t rep = get_object_rep_by_rgn(rgn);

    if (rep == rgn) {
      // retrieve an abstract object info
      auto obj_info = m_obj_info_env[rep];
      // get region counting from m_obj_info_env
      const small_range &count_num = obj_info.count_dom();
      // if count >= 1 report error, this case means it already initialized
      if (count_num <= small_range::oneOrMore()) {
        CRAB_ERROR("region_domain::init_region: ", rgn,
                   " cannot be initialized twice");
      }
      // update m_obj_info_env
      m_obj_info_env.set(rep, region_domain_impl::region_info(
                                  // No references owned by the object
                                  small_range::zero(),
                                  // unused
                                  boolean_value::get_false(),
                                  // unused
                                  rgn.get_type()));
    }

    CRAB_LOG("object-region-init", crab::outs() << "After region_init(" << rgn
                                                << ":" << rgn.get_type()
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

    variable_t rep = get_object_rep_by_rgn(rgn);
    auto old_obj_info = m_obj_info_env[rep];
    const small_range &num_refs = old_obj_info.count_dom();

    if (num_refs.is_one()) {
      // if the counting is already one, use odi map in the future
      move_singleton_to_odi_map(rep);
    }

    m_obj_info_env.set(rgn,
                       region_domain_impl::region_info(
                           old_obj_info.count_dom().increment(),
                           old_obj_info.init_dom(), old_obj_info.type_dom()));

    CRAB_LOG("object-make", crab::outs()
                                << "After ref_make(" << ref << "," << rgn << ":"
                                << rgn.get_type() << "," << size << "," << as
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

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("object-load", CRAB_WARN(domain_name(), "::ref_load: reference ",
                                        ref, " is null."););
      // set_to_bottom();
      crab::CrabStats::count(domain_name() +
                             ".count.ref_load.skipped.null_reference");
      m_base_dom -= res;
      return;
    }

    variable_t rep = get_object_rep_by_rgn(rgn);
    // retrieve an abstract object info
    auto obj_info = m_obj_info_env[rep];

    const small_range &num_refs = obj_info.count_dom();
    if (num_refs.is_zero()) {
      m_base_dom -= res;
    } else if (num_refs.is_one()) { // singleton object
      // strong read
      m_base_dom.assign(res, rgn);
    } else {
      // read from odi map
      // weak read
      auto obj_dom = m_odi_map[rep];
      interval_t ival = obj_dom[rgn];
      assign_interval(m_base_dom, res, ival);
    }

    CRAB_LOG("object-load", crab::outs()
                                 << "After " << res << ":="
                                 << "ref_load(" << rgn << ":" << rgn.get_type()
                                 << "," << ref << ":" << ref.get_type()
                                 << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in region.
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

    variable_t rep = get_object_rep_by_rgn(rgn);
    // retrieve an abstract object info
    auto obj_info = m_obj_info_env[rep];

    const small_range &num_refs = obj_info.count_dom();
    if (num_refs.is_zero()) {
      m_base_dom -= rgn;
    } else if (num_refs.is_one()) { // singleton object, strong update
      if (val.is_constant())
        m_base_dom.assign(rgn, val.get_constant());
      else
        m_base_dom.assign(rgn, val.get_variable());
    } else { // use odi map, weak update
      auto obj_dom = m_odi_map[rep];
      if (val.is_constant()) {
        obj_dom.assign(rgn, val.get_constant());
      } else { // val is variable, this should not be a region variable
        interval_t ival = m_base_dom[val.get_variable()];
        assign_interval(obj_dom, rgn, ival);
      }
      m_odi_map.set(rep, obj_dom);
    }

    CRAB_LOG("object-store",
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

    if (is_bottom()) {
      return;
    }

    m_base_dom.assign(ref2, ref1 + offset);

    CRAB_LOG("object-gep", crab::outs()
                               << "After (" << rgn2 << "," << ref2
                               << ") := ref_gep(" << rgn1 << "," << ref1
                               << " + " << offset << ")=" << *this << "\n";);
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

  // FIXME: The followings are UNDEFINED METHODS
  bool get_allocation_sites(const variable_t &ref,
                            std::vector<allocation_site> &out) override {
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {

    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    // if (!is_bottom()) {
    //   m_base_dom.assign_bool_cst(lhs, rhs);
    // }
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
    return linear_constraint_t::get_false();
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
  void operator-=(const variable_t &v) override {}
  void forget(const variable_vector_t &variables) override {}

  void project(const variable_vector_t &variables) override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {}

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain.
  interval_t operator[](const variable_t &v) override {

    return interval_t::bottom();
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
      for (int i = 0, sz = inputs.size(); i < sz; ++i) {
        error_if_not_rgn(inputs[i]); // this should not happen
        if (inputs[i].get_type().is_unknown_region()) {
          // get_or_insert_gvars does not support unknown regions so we bail
          // out. The base domain shouldn't care about regions anyway.
          return;
        }
        m_flds_rep_map.insert(
            {inputs[i].get_variable(), inputs[0].get_variable()});
      }

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