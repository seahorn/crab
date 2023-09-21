#pragma once

#include <boost/optional.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/object/cow_domain_ref.hpp>
#include <crab/domains/object/object_info.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/term/term_operators.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {
namespace object_domain_impl {

/// @brief A simple comparison function to determine if two values are equal
/// @tparam AbsDom a template for domain value
/// @note this structure is used for patricia tree. Although comparing
/// two values should use operator <= (inclusion test).
/// Here we make a faster comparison by checking address only.
template <class AbsDom> struct object_equal_to {
  bool operator()(const AbsDom &v1, const AbsDom &v2) const {
    bool res = v1 == v2;
    return res;
  }
};

#define ODI_DOMAIN_SCOPED_STATS(NAME) CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define ODI_COUNT_STATS(NAME) CRAB_DOMAIN_COUNT_STATS(NAME, 0)

/* Environment from Key to Value with all lattice operations */

/// @brief Environment from Key to Object value with all lattice operations
/// @tparam Key the key to determine which abstract object (DSA node) is
///         now we represent it as object id
/// @tparam ObjectDom the type name for object domain class
/// @tparam BaseAbsDom the type name for base domain such as Octagons
template <typename Key, typename ObjectDom, typename BaseAbsDom>
class odi_map_domain {
public:
  using object_domain_t = ObjectDom;
  using key_t = Key;
  using base_domain_t = BaseAbsDom;
  using odi_map_domain_t =
      odi_map_domain<key_t, object_domain_t, base_domain_t>;

private:
  // Domain types
  using number_t = typename ObjectDom::number_t;
  using varname_t = typename ObjectDom::varname_t;
  using eq_domain_value_t = typename ObjectDom::eq_domain_value_t;
  using ghost_variables_eq_t = typename ObjectDom::ghost_variables_eq_t;
  using usymb_t = typename eq_domain_value_t::domain_t;
  using base_domain_ref_t = abstract_domain_ref<typename BaseAbsDom::variable_t,
                                              BaseAbsDom>;
  using eq_domain_ref_t = abstract_domain_ref<typename BaseAbsDom::variable_t,
                                              eq_domain_value_t>;

public:
  // Map types
  using summary_domain_t = base_domain_ref_t;
  using cache_domain_t = base_domain_ref_t;
  using eq_domain_t = eq_domain_ref_t;
  using odi_info_t = typename object_domain_impl::object_info;
  using odi_value_t =
      // basic_domain_product2 is a product domain without reduction
      basic_domain_product2<summary_domain_t,
                            basic_domain_product2<cache_domain_t, eq_domain_t>>;
  using map_raw_value_t = basic_domain_product2<odi_info_t, odi_value_t>;
  using map_value_t = std::shared_ptr<map_raw_value_t>;

private:
  using patricia_tree_t =
      ikos::patricia_tree<Key, map_value_t, object_equal_to<map_value_t>>;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using unary_op_t = typename patricia_tree_t::unary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

  using base_dom_variable_t = typename BaseAbsDom::variable_t;
  using base_dom_variable_vector_t = std::vector<base_dom_variable_t>;

  enum cache_status_t { LEFT_ONLY, RIGHT_ONLY, LEFT_AND_RIGHT, NONE };

public:
  using iterator = typename patricia_tree_t::iterator;
  using mapped_type = map_value_t;

private:
  bool m_is_bottom;

  // @def The map contains all objects invariants including values and infos
  // The type of this map is:
  //   object id --> <info, value>
  //   where info is a tuple of <ref counting, cache is used?, cache is dirty?>
  //         value is a tuple of <SUM_DOM, CACHE_DOM, EQ_DOM>
  // If the map does not contain a binding, this means we never see this
  // abstract object.
  // If the map contains info to show the object is singleton,
  // the value is top if we keep object in base domain;
  // Otherwise, we keep the singleton object in CACHE_DOM
  // If the map contains an object is not singleton, ...
  // NOTE: the operation on this map is depended on other domains
  // such as address domain (checking the mru object in CACHE_DOM) eq_domain for
  // registers (used for reduction) basedomain (stored properties for registers,
  // references, etc.)
  patricia_tree_t m_odi_map;

  // Implementation details:
/* clang-format off */
  //       ┌──────────────────────────────────────────────────┐ ╔══════════════╗
  //       │                                                  │ ║shared object:║
  //       │   ┌─────────────────┐                            │ ╚══════════════╝
  //       │   │                 │                            │ ┌─────────────┐ 
  //       │   │equality:        │                            │ │             │ 
  //       │   │                 ├────────────────────────────┼▶│             │ 
  //       │   │                 │                            │ └─────────────┘ 
  //       │   │─────────────────│                            │ ┌─────────────┐ 
  //       │   │                 │                            │ │             │ 
  //       │   │cache:           │────────────────────────────┼▶│             │ 
  //       │   │                 │                            │ └─────────────┘ 
  //       │   │                 │                            │
  //       │   │─────────────────│                            │ ┌─────────────┐ 
  //       │   │                 │                            │ │             │ 
  //       │   │summary:         ├────────────────────────────┼▶│             │ 
  // key:  │   │                 │                            │ └─────────────┘ 
  //  ID ─▶│   └─────────────────┘                            │
  //       │                                                  │
  //       │                                                  │
  //       │┌────────────────────────────────────────────────┐│
  //       ││               │                 │              ││
  //       ││reference count│   cache used    │    cache     ││
  //       ││               │      flag       │  dirty flag  ││
  //       ││────────────────────────────────────────────────││
  //       ││  cache is loaded by a  │  cache is stored by a ││
  //       ││          reg           │          reg          ││
  //       ││   flag for reduction   │   flag for reduction  ││
  //       │└────────────────────────────────────────────────┘│
  //       │╔══════════════╗                                  │
  //       │║shared object:║                                  │
  //       │╚══════════════╝                                  │
  //       └──────────────────────────────────────────────────┘
/* clang-format on */

  /// @brief A special class to compute join / widening when merging two trees
  class join_or_widening_op : public binary_op_t {
  private:
    const object_domain_t &m_left;
    const object_domain_t &m_right;
    base_domain_t &m_l_base_dom;
    base_domain_t &m_r_base_dom;
    bool m_is_join;

    /// @brief A special operation for the case where
    ///       state one contains a singleton object in the base domain
    ///       state two contains a non-singleton object in the odi map
    /// @param key the object id
    /// @param s_single the state containing singleton object
    /// @param non_single_odi the state containing non-singleton
    /// @return a new non-singleton by join two states
    odi_value_t join_or_widen_singleton_with_non_singleton(
        const key_t &key, const object_domain_t &s_single,
        const odi_value_t &non_single_odi) const {
      odi_value_t res_prod = non_single_odi;
      // Project on fields only
      base_domain_t singleton_base = s_single.project_singleton(key);

      if (m_is_join) {
        *res_prod.first() |= singleton_base;
      } else {
        *res_prod.first() = *res_prod.first() || singleton_base;
      }
      return res_prod;
    }

    /// @brief the join or widen odi
    /// @param key the object id
    /// @param l the odi stored on the left state
    /// @param r the odi stored on the right state
    /// @return a new odi after join or widening
    map_value_t join_or_widen_odi(const key_t &key, const map_value_t &l,
                                  const map_value_t &r) {
      ODI_DOMAIN_SCOPED_STATS(".join");

      const odi_info_t &l_obj_info = l->first();
      const odi_value_t &l_odi_val = l->second();
      const odi_info_t &r_obj_info = r->first();
      const odi_value_t &r_odi_val = r->second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      map_raw_value_t out_val;
      cache_status_t status = cache_status_t::NONE;
      // commit cache if dirty
      if (l_num_refs == small_range::oneOrMore() &&
          l_obj_info.cachedirty_val().is_true()) {
        status = cache_status_t::LEFT_ONLY;
      }
      if (r_num_refs == small_range::oneOrMore() &&
          r_obj_info.cachedirty_val().is_true()) {
        status = status == cache_status_t::LEFT_ONLY
                     ? cache_status_t::LEFT_AND_RIGHT
                     : cache_status_t::RIGHT_ONLY;
      }

      if (!crab_domain_params_man::get().singletons_in_base()) {
        switch (status) {
        case LEFT_ONLY: {
          odi_value_t new_l_odi_val = l_odi_val;
          m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val, l_obj_info,
                                       key);
          // pairwise join
          out_val.second() = m_is_join ? new_l_odi_val | r_odi_val
                                       : new_l_odi_val || r_odi_val;
        } break;
        case RIGHT_ONLY: {
          odi_value_t new_r_odi_val = r_odi_val;
          m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val, r_obj_info,
                                        key);
          // pairwise join
          out_val.second() = m_is_join ? l_odi_val | new_r_odi_val
                                       : l_odi_val || new_r_odi_val;
        } break;
        case LEFT_AND_RIGHT: {
          odi_value_t new_l_odi_val = l_odi_val;
          m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val, l_obj_info,
                                       key);
          odi_value_t new_r_odi_val = r_odi_val;
          m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val, r_obj_info,
                                        key);
          // pairwise join
          out_val.second() = m_is_join ? new_l_odi_val | new_r_odi_val
                                       : new_l_odi_val || new_r_odi_val;
        } break;
        case NONE: {
          // pairwise join
          out_val.second() =
              m_is_join ? l_odi_val | r_odi_val : l_odi_val || r_odi_val;
        } break;
        default:
          // this is an unreachable default clause
          CRAB_ERROR(typeid(join_or_widening_op).name(), "::", __func__,
                     "Unhandled special enum constant!");
          break;
        }
      } else {
        // NOTE: if we keep singleton object in the base, the join operation
        // should cover the case for singleton with non singleton
        if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
          if (r_obj_info.cachedirty_val().is_true()) {
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            out_val.second() = join_or_widen_singleton_with_non_singleton(
                key, m_left, new_r_odi_val);
          } else {
            out_val.second() = join_or_widen_singleton_with_non_singleton(
                key, m_left, r_odi_val);
          }
        } else if (r_num_refs.is_one() &&
                   l_num_refs == small_range::oneOrMore()) {
          if (l_obj_info.cachedirty_val().is_true()) {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            out_val.second() = join_or_widen_singleton_with_non_singleton(
                key, m_right, new_l_odi_val);
          } else {
            out_val.second() = join_or_widen_singleton_with_non_singleton(
                key, m_right, l_odi_val);
          }
        } else {
          switch (status) {
          case LEFT_ONLY: {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            // pairwise join
            out_val.second() = m_is_join ? new_l_odi_val | r_odi_val
                                         : new_l_odi_val || r_odi_val;
          } break;
          case RIGHT_ONLY: {
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            // pairwise join
            out_val.second() = m_is_join ? l_odi_val | new_r_odi_val
                                         : l_odi_val || new_r_odi_val;
          } break;
          case LEFT_AND_RIGHT: {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            // pairwise join
            out_val.second() = m_is_join ? new_l_odi_val | new_r_odi_val
                                         : new_l_odi_val || new_r_odi_val;
          } break;
          case NONE: {
            // pairwise join
            out_val.second() =
                m_is_join ? l_odi_val | r_odi_val : l_odi_val || r_odi_val;
          } break;
          default:
            // this is an unreachable default clause
            CRAB_ERROR(typeid(join_or_widening_op).name(), "::", __func__,
                       "Unhandled special enum constant!");
            break;
          }
        }
      }
      // After the join / widening, update the object info
      out_val.first() =
          odi_info_t(l_num_refs | r_num_refs,
                     // Cache is not used
                     boolean_value::get_false(),
                     // Cache is not dirty
                     boolean_value::get_false(), false, false);
      return std::make_shared<map_raw_value_t>(out_val);
    }

    std::string domain_name() const { return "ODI"; }

  protected:
    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key &key, const map_value_t &x, const map_value_t &y) override {
      if (x == y) {
        return {false, boost::optional<map_value_t>(x)};
      }
      map_value_t z = join_or_widen_odi(key, x, y);
      if (z->is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    }

    virtual bool default_is_absorbing() override { return true; }

  public:
    join_or_widening_op(const object_domain_t &left,
                        const object_domain_t &right, base_domain_t &l_base_dom,
                        base_domain_t &r_base_dom, bool is_join)
        : m_left(left), m_right(right), m_l_base_dom(l_base_dom),
          m_r_base_dom(r_base_dom), m_is_join(is_join) {}
  }; // class join_op

  /// @brief A special class to compute meet / narrowing when merging two trees
  class meet_or_narrow_op : public binary_op_t {
  private:
    const object_domain_t &m_left;
    const object_domain_t &m_right;
    base_domain_t &m_l_base_dom;
    base_domain_t &m_r_base_dom;
    bool m_is_meet;

    /// @brief A special operation for the case where
    ///       state one contains a singleton object in the base domain
    ///       state two contains a non-singleton object in the odi map
    /// @param key the object id
    /// @param single the state containing singleton object
    /// @param non_single the state containing non-singleton
    void meet_or_narrow_non_singleton_with_singleton(
        const key_t &key, base_domain_t &single,
        const odi_value_t &non_single) const {
      const base_domain_t &summary = *non_single.first();

      // Be careful of the following operation if type of base dom and object
      // dom are different
      if (m_is_meet) {
        single &= summary;
      } else {
        single = single && summary;
      }
    }

    /// @brief the meet or narrowing odi
    /// @param key the object id
    /// @param l the odi stored on the left state
    /// @param r the odi stored on the right state
    /// @return a new odi after meet or narrowing
    map_value_t meet_or_narrow_odi(const key_t &key, const map_value_t &l,
                                   const map_value_t &r) {
      const odi_info_t &l_obj_info = l->first();
      const odi_value_t &l_odi_val = l->second();
      const odi_info_t &r_obj_info = r->first();
      const odi_value_t &r_odi_val = r->second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      map_raw_value_t out_val;
      cache_status_t status = cache_status_t::NONE;
      // commit cache if dirty
      if (l_num_refs == small_range::oneOrMore() &&
          l_obj_info.cachedirty_val().is_true()) {
        status = cache_status_t::LEFT_ONLY;
      }
      if (r_num_refs == small_range::oneOrMore() &&
          r_obj_info.cachedirty_val().is_true()) {
        status = status == cache_status_t::LEFT_ONLY
                     ? cache_status_t::LEFT_AND_RIGHT
                     : cache_status_t::RIGHT_ONLY;
      }
      if (!crab_domain_params_man::get().singletons_in_base()) {
        switch (status) {
        case LEFT_ONLY: {
          odi_value_t new_l_odi_val = l_odi_val;
          m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val, l_obj_info,
                                       key);
          // pairwise meet
          out_val.second() = m_is_meet ? new_l_odi_val & r_odi_val
                                       : new_l_odi_val && r_odi_val;
        } break;
        case RIGHT_ONLY: {
          odi_value_t new_r_odi_val = r_odi_val;
          m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val, r_obj_info,
                                        key);
          // pairwise meet
          out_val.second() = m_is_meet ? l_odi_val & new_r_odi_val
                                       : l_odi_val && new_r_odi_val;
        } break;
        case LEFT_AND_RIGHT: {
          odi_value_t new_l_odi_val = l_odi_val;
          m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val, l_obj_info,
                                       key);
          odi_value_t new_r_odi_val = r_odi_val;
          m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val, r_obj_info,
                                        key);
          // pairwise meet
          out_val.second() = m_is_meet ? new_l_odi_val & new_r_odi_val
                                       : new_l_odi_val && new_r_odi_val;
        } break;
        case NONE: {
          // pairwise meet
          out_val.second() =
              m_is_meet ? l_odi_val & r_odi_val : l_odi_val && r_odi_val;
        } break;
        default:
          // this is an unreachable default clause
          CRAB_ERROR(typeid(meet_or_narrow_op).name(), "::", __func__,
                     "Unhandled special enum constant!");
          break;
        }
      } else {
        // NOTE: if we keep singleton object in the base, the meet operation
        // should cover the case for singleton with non singleton
        if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
          if (r_obj_info.cachedirty_val().is_true()) {
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            meet_or_narrow_non_singleton_with_singleton(key, m_l_base_dom,
                                                        new_r_odi_val);
          } else {
            meet_or_narrow_non_singleton_with_singleton(key, m_l_base_dom,
                                                        r_odi_val);
          }
        } else if (r_num_refs.is_one() &&
                   l_num_refs == small_range::oneOrMore()) {
          if (l_obj_info.cachedirty_val().is_true()) {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            meet_or_narrow_non_singleton_with_singleton(key, m_r_base_dom,
                                                        new_l_odi_val);
          } else {
            meet_or_narrow_non_singleton_with_singleton(key, m_r_base_dom,
                                                        l_odi_val);
          }
        } else {
          switch (status) {
          case LEFT_ONLY: {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            // pairwise meet
            out_val.second() = m_is_meet ? new_l_odi_val & r_odi_val
                                         : new_l_odi_val && r_odi_val;
          } break;
          case RIGHT_ONLY: {
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            // pairwise meet
            out_val.second() = m_is_meet ? l_odi_val & new_r_odi_val
                                         : l_odi_val && new_r_odi_val;
          } break;
          case LEFT_AND_RIGHT: {
            odi_value_t new_l_odi_val = l_odi_val;
            m_left.commit_cache_if_dirty(m_l_base_dom, new_l_odi_val,
                                         l_obj_info, key);
            odi_value_t new_r_odi_val = r_odi_val;
            m_right.commit_cache_if_dirty(m_r_base_dom, new_r_odi_val,
                                          r_obj_info, key);
            // pairwise meet
            out_val.second() = m_is_meet ? new_l_odi_val & new_r_odi_val
                                         : new_l_odi_val && new_r_odi_val;
          } break;
          case NONE: {
            // pairwise meet
            out_val.second() =
                m_is_meet ? l_odi_val & r_odi_val : l_odi_val && r_odi_val;
          } break;
          default:
            // this is an unreachable default clause
            CRAB_ERROR(typeid(meet_or_narrow_op).name(), "::", __func__,
                       "Unhandled special enum constant!");
            break;
          }
        }
      }
      out_val.first() = odi_info_t(l_num_refs & r_num_refs,
                        // Cache is not used
                        boolean_value::get_false(),
                        // Cache is not dirty
                        boolean_value::get_false(), false, false);
      return std::make_shared<map_raw_value_t>(out_val);
    }

  protected:
    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key &key, const map_value_t &x, const map_value_t &y) override {
      if (x == y) {
        return {false, boost::optional<map_value_t>(x)};
      }
      map_value_t z = meet_or_narrow_odi(key, x, y);
      if (z->is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    }

    virtual bool default_is_absorbing() override { return false; }

  public:
    meet_or_narrow_op(const object_domain_t &left, const object_domain_t &right,
                      base_domain_t &l_base_dom, base_domain_t &r_base_dom,
                      bool is_meet)
        : m_left(left), m_right(right), m_l_base_dom(l_base_dom),
          m_r_base_dom(r_base_dom), m_is_meet(is_meet) {}
  }; // class meet_or_narrow_op

  /// @brief A special class to compute widening with thresholds when merging
  /// two trees
  /// TODO: this class is not used because
  ///       object domain does not apply widening with thresholds
  template <typename Thresholds>
  class widening_thresholds_op : public binary_op_t {
    const Thresholds &m_ts;

  public:
    widening_thresholds_op(const Thresholds &ts) : m_ts(ts) {}

    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key & /*key*/, const map_value_t &x,
          const map_value_t &y) override {
      map_value_t z = x->widening_thresholds(*y, m_ts);
      if (z->is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    };

    virtual bool default_is_absorbing() override { return true; }
  }; // class widening_thresholds_op

  /// @brief A special class to operator \subseteq when merging two trees
  class inclusion_test_op : public partial_order_t {
    virtual bool leq(const map_value_t &x, const map_value_t &y) override {
      if (x == y) {
        return true;
      }
      const odi_info_t &l_obj_info = x->first();
      const odi_value_t &l_odi_val = x->second();
      const odi_info_t &r_obj_info = y->first();
      const odi_value_t &r_odi_val = y->second();
      bool res = l_obj_info <= r_obj_info;
      res &= l_odi_val.first() <= r_odi_val.first();
      res &= l_odi_val.second().first() <= r_odi_val.second().first();
      return res;
    }

    virtual bool default_is_top() override { return true; }
  }; // class inclusion_test_op

  /// @brief apply operation when merge two trees
  /// @param o the binary operation for merging two values with the same key
  /// @param t1 a tree
  /// @param t2 a tree
  /// @param is_bottom a bool value indicates the result is bottom or not
  /// @return a tree after merged
  static patricia_tree_t apply_operation(binary_op_t &o, patricia_tree_t t1,
                                         const patricia_tree_t &t2,
                                         bool &is_bottom) {
    is_bottom = t1.merge_with(t2, o);
    return t1;
  }

  /// @brief apply operation when merge two trees (self-apply)
  /// @param o the binary operation for merging two values with the same key
  /// @param t1 a tree, the merging result will update this tree
  /// @param t2 a tree
  /// @param is_bottom a bool value indicates the result is bottom or not
  static void apply_self_operation(binary_op_t &o, patricia_tree_t &t1,
                                   const patricia_tree_t &t2, bool &is_bottom) {
    is_bottom = t1.merge_with(t2, o);
  }

  /// @brief a simple log method
  /// @param o crab ostream
  void write(crab::crab_os &o) {
    if (is_bottom()) {
      o << "_|_";
    } else {
      o << "{";
      auto it = m_odi_map.begin();
      auto end = m_odi_map.end();
      while (it != end) {
        key_t k = it->first;
        k.write(o);
        o << " -> ";
        map_value_t v = it->second;
        v->write(o);
        ++it;
        if (it != end) {
          o << "; ";
        }
      }
      o << "}";
    }
  }

  odi_map_domain(patricia_tree_t &&t)
      : m_is_bottom(false), m_odi_map(std::move(t)) {}

  odi_map_domain(bool b) : m_is_bottom(!b) {}

public:
  /**------------------ Begin domain APIs ------------------**/
  odi_map_domain() : m_is_bottom(false) {}
  // NOTE: The copy constructor is a shallow copy to keep subtrees to be
  // sharable
  odi_map_domain(const odi_map_domain_t &o) = default;
  odi_map_domain(odi_map_domain_t &&o) = default;
  odi_map_domain_t &operator=(const odi_map_domain_t &o) = default;
  odi_map_domain_t &operator=(odi_map_domain_t &&o) = default;

  static odi_map_domain_t top() { return odi_map_domain_t(); }
  static odi_map_domain_t bottom() { return odi_map_domain_t(false); }

  bool is_bottom() const { return m_is_bottom; }

  bool is_top() const { return (!is_bottom() && m_odi_map.size() == 0); }

  void set_to_bottom() {
    m_is_bottom = true;
    m_odi_map.clear();
  }

  void set_to_top() {
    m_is_bottom = false;
    m_odi_map.clear();
  }

  // Inclusion test
  bool operator<=(const odi_map_domain_t &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    } else {
      inclusion_test_op leq_op;
      return m_odi_map.leq(o.m_odi_map, leq_op);
    }
  }

  // Join
  odi_map_domain_t join(const odi_map_domain_t &o, const object_domain_t &left,
                        const object_domain_t &right, base_domain_t &l_base_dom,
                        base_domain_t &r_base_dom) const {
    ODI_DOMAIN_SCOPED_STATS(".join");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() || o.is_top()) {
      odi_map_domain_t abs;
      abs.set_to_top();
      return abs;
    } else {
      join_or_widening_op jop(left, right, l_base_dom, r_base_dom, true);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Self join
  void compound_join(const odi_map_domain_t &o, const object_domain_t &left,
                     const object_domain_t &right, base_domain_t &l_base_dom,
                     base_domain_t &r_base_dom) {
    ODI_DOMAIN_SCOPED_STATS(".join");

    if (is_bottom()) { // this is bot, assign this by o
      *this = o;
    } else if (o.is_bottom()) { // o is bot, nothing change
    } else if (is_top() || o.is_top()) {
      set_to_top();
    } else {
      join_or_widening_op jop(left, right, l_base_dom, r_base_dom, true);
      bool is_bottom = false /*unused*/;
      apply_self_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
    }
  }

  // Meet
  odi_map_domain_t meet(const odi_map_domain_t &o, const object_domain_t &left,
                        const object_domain_t &right, base_domain_t &l_base_dom,
                        base_domain_t &r_base_dom) const {
    ODI_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    } else {
      meet_or_narrow_op mop(left, right, l_base_dom, r_base_dom, true);
      bool is_bottom;
      patricia_tree_t res =
          apply_operation(mop, m_odi_map, o.m_odi_map, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return odi_map_domain_t(std::move(res));
      }
    }
  }

  // Self meet
  void compound_meet(const odi_map_domain_t &o, const object_domain_t &left,
                     const object_domain_t &right, base_domain_t &l_base_dom,
                     base_domain_t &r_base_dom) {
    ODI_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      *this = o;
    } else {
      meet_or_narrow_op mop(left, right, l_base_dom, r_base_dom, true);
      bool is_bottom = false /*unused*/;
      apply_self_operation(mop, m_odi_map, o.m_odi_map, is_bottom);
      if (is_bottom) {
        set_to_bottom();
      }
    }
  }

  // Widening
  odi_map_domain_t widening(const odi_map_domain_t &o,
                            const object_domain_t &left,
                            const object_domain_t &right,
                            base_domain_t &l_base_dom,
                            base_domain_t &r_base_dom) const {
    ODI_DOMAIN_SCOPED_STATS(".widening");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      join_or_widening_op wop(left, right, l_base_dom, r_base_dom, false);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(wop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Widening with thresholds
  template <typename Thresholds>
  odi_map_domain_t widening_thresholds(const odi_map_domain_t &o,
                                       const Thresholds &ts) const {
    ODI_DOMAIN_SCOPED_STATS(".widening");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      widening_thresholds_op<Thresholds> wtop;
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(wtop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Narrowing
  odi_map_domain_t narrowing(const odi_map_domain_t &o,
                             const object_domain_t &left,
                             const object_domain_t &right,
                             base_domain_t &l_base_dom,
                             base_domain_t &r_base_dom) const {
    ODI_DOMAIN_SCOPED_STATS(".narrowing");
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    } else {
      meet_or_narrow_op nop(left, right, l_base_dom, r_base_dom, false);
      bool is_bottom;
      patricia_tree_t res =
          apply_operation(nop, m_odi_map, o.m_odi_map, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return odi_map_domain_t(std::move(res));
      }
    }
  }

  // Project onto
  // NOTE: this operation is on key not on values
  void project(const std::vector<key_t> &keys) {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }
    const double factor = 0.6;
    const size_t totals = m_odi_map.size();
    const size_t remained = keys.size();
    if (remained >= (size_t)(totals * factor)) {
      // most keys remained, remove the rest
      // binary search on a sorting key array --- O(nlgn)
      std::vector<key_t> remove_keys;
      remove_keys.reserve(totals - remained);
      std::vector<key_t> sorted_keys(keys);
      std::sort(sorted_keys.begin(), sorted_keys.end());
      for (auto it = m_odi_map.begin(); it != m_odi_map.end(); ++it) {
        key_t &key = it->first;
        if (!std::binary_search(sorted_keys.begin(), sorted_keys.end(), key)) {
          remove_keys.push_back(key);
        }
      }
      forget(remove_keys);
    } else {
      odi_map_domain_t new_val;
      // less keys remained, construct a new one
      for (auto &key : keys) {
        new_val.set(key, at(keys));
      }
      std::swap(*this, new_val);
    }
  }

  // Forget
  // NOTE: this operation is on key not on values
  void operator-=(const key_t &k) {
    ODI_DOMAIN_SCOPED_STATS(".forget");

    if (!is_bottom()) {
      m_odi_map.remove(k);
    }
  }

  // NOTE: this operation is on key not on values
  void forget(const std::vector<key_t> &keys) {
    ODI_DOMAIN_SCOPED_STATS(".forget");

    if (!is_bottom()) {
      for (auto &key : keys) {
        m_odi_map.remove(key);
      }
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const odi_map_domain_t &d) {
    d.write(o);
    return o;
  }

  std::string domain_name() const { return "ODIMapDomain"; }
  /**------------------ End domain APIs ------------------**/
  /**------------------ Begin Map APIs ------------------**/
  iterator begin() const {
    if (is_bottom()) {
      CRAB_ERROR(domain_name(), "::", __func__,
                 " trying to invoke iterator on bottom");
    } else {
      return m_odi_map.begin();
    }
  }
  iterator end() const {
    if (is_bottom()) {
      CRAB_ERROR(domain_name(), "::", __func__,
                 " trying to invoke iterator on bottom");
    } else {
      return m_odi_map.end();
    }
  }

  /// @brief check an object exists or not
  /// @param key the object id
  /// @return true is it exists
  bool contains(const key_t &key) const {
    return m_odi_map.lookup(key) != boost::none;
  }

  /// @brief update odi map by a new <obj_id, <info, value>>
  /// @param key object id
  /// @param v the new value by a product of <info, value>
  void set(const key_t &key, map_raw_value_t &&v) {
    ODI_DOMAIN_SCOPED_STATS(".set");
    if (!is_bottom()) {
      if (v.is_bottom()) {
        set_to_bottom();
      } else if (v.is_top()) {
        m_odi_map.remove(key);
      } else {
        map_value_t v_ptr = std::make_shared<map_raw_value_t>(std::move(v));
        m_odi_map.insert(key, v_ptr);
      }
    }
  }

  /// @brief remove an object from the map
  /// @param key object id
  void remove(const key_t &key) { m_odi_map.remove(key); }

  /// @brief find a value by giving a key
  /// @param k object id
  /// @return a pointer to const if the value exists
  /// The value to the map is a constant object which should not be modified
  const map_raw_value_t *find(const key_t &k) const {
    assert(!is_bottom());
    auto val_ptr = m_odi_map.find(k);
    if (val_ptr) {
      return val_ptr->get();
    } else {
      return nullptr;
    }
  }

  /// @brief find a value by giving a key, otherwise, return a top value
  /// @param k object id
  /// @return a value if the key exists or a top if it does not
  map_raw_value_t at(const key_t &k) const {
    if (is_bottom()) {
      map_raw_value_t res;
      return res.make_bottom();
    } else {
      boost::optional<map_value_t> v = m_odi_map.lookup(k);
      if (v) {
        map_value_t val_ptr = *v;
        return *val_ptr;
      } else {
        map_raw_value_t res;
        return res.make_top();
      }
    }
  }

  /// @brief get the size of the tree
  /// @return the size
  std::size_t size() const {
    if (is_bottom()) {
      return 0;
    } else if (is_top()) {
      CRAB_ERROR(domain_name(), "::", __func__, " is undefined if top");
    } else {
      return m_odi_map.size();
    }
  }

  /// @brief get all object ids from the map
  /// @return a vector stored all keys
  std::vector<key_t> keys() const {
    if (is_bottom()) {
      return {};
    } else if (is_top()) {
      CRAB_ERROR(domain_name(), "::", __func__, " is undefined if top");
    } else {
      std::vector<key_t> keys;
      for (auto kv : m_odi_map) {
        const key_t &key = kv.first;
        keys.push_back(key);
      }
      return keys;
    }
  }

  /// @brief a special log method for <SUM_DOM, CACHE_DOM, EQ_DOM>
  /// @param o crab ostream
  /// @param prod a odi value
  void odi_val_write(crab_os &o, const odi_value_t &prod) const {
    o << "summary: ";
    prod.first().write(o);
    o << ", cache: ";
    prod.second().first().write(o);
    o << ", eq_fields: ";
    prod.second().second().write(o);
    o << ")";
  }

  /// @brief a special log method for <info, value>
  /// @param o crab ostream
  /// @param prod the value stored on the map
  void odi_write(crab_os &o, const map_raw_value_t &prod) const {
    const odi_info_t &obj_info = prod.first();
    const odi_value_t &obj_val = prod.second();
    o << "(";
    o << "info: ";
    obj_info.write(o);
    o << ", ";
    auto &num_refs = obj_info.refcount_val();
    if (num_refs.is_zero()) {
      o << "not init";
    } else if (num_refs.is_one() &&
               crab_domain_params_man::get().singletons_in_base()) {
      o << "values in base_dom";
    } else {
      odi_val_write(o, obj_val);
    }
  }

  void odi_write(crab_os &o, const map_value_t &prod) const {
    odi_write(o, *prod);
  }

  void dump() const { write(crab::outs()); }

  /**------------------ End Map APIs ------------------**/

  /**------------------ Begin Cache APIs ------------------**/
  // commit cache contents into summary
  void commit_cache(summary_domain_t &summary, cache_domain_t &cache) const {
    // join summary with cache
    summary |= cache;
  }

  // update cache contents from summary to cache
  void update_cache(summary_domain_t &summary, cache_domain_t &cache) const {
    // copy summary to cache
    cache = summary;
  }

  void commit_cache_if_dirty(const odi_info_t &obj_info,
                             odi_value_t &obj_val) const {
    ODI_DOMAIN_SCOPED_STATS(".commit_cache");
    summary_domain_t &sum = odi_map_domain_t::object_sum_val(obj_val);
    cache_domain_t &cache = odi_map_domain_t::object_cache_val(obj_val);
    eq_domain_t &eq_fld = odi_map_domain_t::object_eq_val(obj_val);
    if (obj_info.cachedirty_val().is_true()) {
      // commit cache if the cache is dirty
      commit_cache(sum, cache);
    }
    // (2) clear cache domain
    cache.set_to_top();
    // (3) clear eq field domain
    eq_fld.set_to_top();
  }

  bool invalidate_cache_if_miss(
      key_t &key, object_domain_t &abs_state,
      ghost_variables_eq_t &&rgn_eq_gvars, boost::optional<usymb_t> &reg_symb,
      boost::optional<std::pair<usymb_t, usymb_t>> &offset_size_symb,
      bool is_ref_mru) {
    const map_raw_value_t *obj_prod_ref = find(key);

    bool update_new_mru = false;
    // retrieve an abstract object info
    if (!obj_prod_ref) {
      CRAB_ERROR(domain_name(), "::", __func__, ": accessing ", key,
                 " is not found on the odi map");
    }

    odi_info_t out_obj_info = obj_prod_ref->first();
    auto cache_used = out_obj_info.cacheused_val();

    // NOTE: copy the object value is required
    odi_value_t out_prod = obj_prod_ref->second();
    eq_domain_t &eq_flds_dom = out_prod.second().second();

    /* Cache missed condition:
        cache is empty or current reference does not refer to the mru object
    */
    if (cache_used.is_false() || is_ref_mru == false) {
      if (out_obj_info.refcount_val() == small_range::oneOrMore()) {
        // Step1: commit cache if the cache is dirty
        abs_state.commit_cache_if_dirty(out_prod, out_obj_info, key);
      }
      else if (out_obj_info.cache_reg_loaded_val()) {
        out_obj_info.cache_reg_loaded_val() = false;
        abs_state.apply_reduction_from_object_to_base(out_prod, key);
      }
      // Step2: update cache for new MRU object)
      update_cache(out_prod.first(), out_prod.second().first());
      // Step3: update address dom and object info
      update_new_mru = true;

      out_obj_info = odi_info_t(
          out_obj_info.refcount_val(),
          // Cache is used
          boolean_value::get_true(),
          // Cache is not dirty
          boolean_value::get_false(),
          // Cache is loaded to a reg?
          reg_symb == boost::none,
          // Cache is stored from a reg?
          reg_symb != boost::none);
    } else {
      if (reg_symb) {
        if (out_obj_info.cache_reg_loaded_val()) {
          out_obj_info.cache_reg_loaded_val() = false;
          abs_state.apply_reduction_from_object_to_base(out_prod, key);
        }
      } else {
        if (out_obj_info.cache_reg_stored_val()) {
          out_obj_info.cache_reg_stored_val() = false;
          abs_state.apply_reduction_from_base_to_object(out_prod, key);
        }
      }
      out_obj_info.cache_reg_loaded_val() = reg_symb == boost::none;
      out_obj_info.cache_reg_stored_val() = reg_symb != boost::none;
    }
    // Step4: update field == reg, the equality represents either:
    //        reg := load_ref(ref, field), or
    //        store_ref(ref, field, reg)
    if (reg_symb == boost::none) { // it means load_ref
      // assigning the same symbol as rgn
      reg_symb =
          abs_state.get_or_insert_symbol(*eq_flds_dom, rgn_eq_gvars.get_var());
    }
    (*eq_flds_dom).set(rgn_eq_gvars.get_var(), *reg_symb);
    if (rgn_eq_gvars.has_offset_and_size()) {
      if (offset_size_symb == boost::none) {
        offset_size_symb = std::make_pair(
            abs_state.get_or_insert_symbol(
                *eq_flds_dom, rgn_eq_gvars.get_offset_and_size().get_offset()),
            abs_state.get_or_insert_symbol(
                *eq_flds_dom, rgn_eq_gvars.get_offset_and_size().get_size()));
      }
      (*eq_flds_dom)
          .set(rgn_eq_gvars.get_offset_and_size().get_offset(),
               std::get<0>(*offset_size_symb));
      (*eq_flds_dom)
          .set(rgn_eq_gvars.get_offset_and_size().get_size(),
               std::get<1>(*offset_size_symb));
    }

    // Step5: update odi map
    set(key, map_raw_value_t(std::move(out_obj_info), std::move(out_prod)));

    return update_new_mru;
  }
  /**------------------ End Cache APIs ------------------**/

  /**------------------ Begin Product Domain APIs ------------------**/
  static const odi_info_t &object_info_val(const map_raw_value_t &prod) {
    return prod.first();
  }

  static odi_info_t &object_info_val(map_raw_value_t &prod) {
    return prod.first();
  }

  static const odi_value_t &object_odi_val(const map_raw_value_t &prod) {
    return prod.second();
  }

  static odi_value_t &object_odi_val(map_raw_value_t &prod) {
    return prod.second();
  }

  static const summary_domain_t &object_sum_val(const map_raw_value_t &prod) {
    return prod.second().first();
  }

  static const base_domain_t &object_sum_raw_val(const map_raw_value_t &prod) {
    return *prod.second().first();
  }

  static summary_domain_t &object_sum_val(map_raw_value_t &prod) {
    return prod.second().first();
  }

  static base_domain_t &object_sum_raw_val(map_raw_value_t &prod) {
    return *prod.second().first();
  }

  static const summary_domain_t &object_sum_val(const odi_value_t &odi_val) {
    return odi_val.first();
  }

  static const base_domain_t &object_sum_raw_val(const odi_value_t &odi_val) {
    return *odi_val.first();
  }

  static summary_domain_t &object_sum_val(odi_value_t &odi_val) {
    return odi_val.first();
  }

  static base_domain_t &object_sum_raw_val(odi_value_t &odi_val) {
    return *odi_val.first();
  }

  static const cache_domain_t &object_cache_val(const map_raw_value_t &prod) {
    return prod.second().second().first();
  }

  static const base_domain_t &
  object_cache_raw_val(const map_raw_value_t &prod) {
    return *prod.second().second().first();
  }

  static cache_domain_t &object_cache_val(map_raw_value_t &prod) {
    return prod.second().second().first();
  }

  static base_domain_t &object_cache_raw_val(map_raw_value_t &prod) {
    return *prod.second().second().first();
  }

  static const cache_domain_t &object_cache_val(const odi_value_t &odi_val) {
    return odi_val.second().first();
  }

  static const base_domain_t &object_cache_raw_val(const odi_value_t &odi_val) {
    return *odi_val.second().first();
  }

  static cache_domain_t &object_cache_val(odi_value_t &odi_val) {
    return odi_val.second().first();
  }

  static base_domain_t &object_cache_raw_val(odi_value_t &odi_val) {
    return *odi_val.second().first();
  }

  static const eq_domain_t &object_eq_val(const map_raw_value_t &prod) {
    return prod.second().second().second();
  }

  static const eq_domain_value_t &
  object_eq_raw_val(const map_raw_value_t &prod) {
    return *prod.second().second().second();
  }

  static eq_domain_t &object_eq_val(map_raw_value_t &prod) {
    return prod.second().second().second();
  }

  static eq_domain_value_t &object_eq_raw_val(map_raw_value_t &prod) {
    return *prod.second().second().second();
  }

  static const eq_domain_t &object_eq_val(const odi_value_t &odi_val) {
    return odi_val.second().second();
  }

  static const eq_domain_value_t &
  object_eq_raw_val(const odi_value_t &odi_val) {
    return *odi_val.second().second();
  }

  static eq_domain_t &object_eq_val(odi_value_t &odi_val) {
    return odi_val.second().second();
  }

  static eq_domain_value_t &object_eq_raw_val(odi_value_t &odi_val) {
    return *odi_val.second().second();
  }

  /**------------------ End Product Domain APIs ------------------**/
}; // class odi_map_domain
} // end namespace object_domain_impl
} // end namespace domains
} // end namespace crab