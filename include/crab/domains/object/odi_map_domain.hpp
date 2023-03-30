#pragma once

#include <boost/optional.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/object/object_info.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/domains/term/term_operators.hpp>
#include <crab/domains/uf_domain.hpp>
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
    // return v1 <= v2 && v2 <= v1;
    return &v1 == &v2;
  }
};

/* Environment from Key to Value with all lattice operations */

/// @brief Environment from Key to Object value with all lattice operations
/// @tparam Key the key to determine which abstract object (DSA node) is
///         now we represent it as object id
/// @tparam ObjectDom the type name for object domain class
/// @tparam BaseAbsDom the type name for base domain such as Octagons
template <typename Key, typename ObjectDom, typename BaseAbsDom>
class odi_map_domain {
private:
  // Domain types
  using number_t = typename ObjectDom::number_t;
  using varname_t = typename ObjectDom::varname_t;
  using eq_domain_t = typename ObjectDom::eq_fields_abstract_domain_t;
  // Map types
  using odi_info_t = typename object_domain_impl::object_info;
  using odi_value_t =
      // basic_domain_product2 is a product domain without reduction
      basic_domain_product2<BaseAbsDom,
                            basic_domain_product2<BaseAbsDom, eq_domain_t>>;
  using map_value_t = basic_domain_product2<odi_info_t, odi_value_t>;
  using patricia_tree_t =
      ikos::patricia_tree<Key, map_value_t, object_equal_to<map_value_t>>;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using unary_op_t = typename patricia_tree_t::unary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

  using base_dom_variable_t = typename BaseAbsDom::variable_t;
  using base_dom_variable_vector_t = std::vector<base_dom_variable_t>;

public:
  using odi_map_domain_t = odi_map_domain<Key, ObjectDom, BaseAbsDom>;
  using object_domain_t = ObjectDom;
  using iterator = typename patricia_tree_t::iterator;
  using key_t = Key;
  using mapped_type = map_value_t;

private:
  bool m_is_bottom;

  // @def The map contains all objects invariants including values and infos
  // The type of this map is:
  //   object id --> <info, value>
  //   where info is a tuple of <ref counting, cache is used?, cache is dirty?>
  //         value is a tuple of <SUM_DOM, CACHE_DOM, UF_DOM>
  // If the map does not contain a binding, this means we never see this
  // abstract object.
  // If the map contains info to show the object is singleton,
  // the value is top if we keep object in base domain;
  // Otherwise, we keep the singleton object in CACHE_DOM
  // If the map contains an object is not singleton, ...
  // NOTE: the operation on this map is depended on other domains
  // such as address domain (checking the mru object in CACHE_DOM) ufdomain for
  // registers (used for reduction) basedomain (stored properties for registers,
  // references, etc.) Sol1: impement the odi map inside object domain, do not
  // create a new domain Sol2: during domain operation such as join, meet or
  // inclusion test, passing the abstract state of the object domain.
  patricia_tree_t m_odi_map;

  /// @brief A special class to compute join / widening when merging two trees
  class join_or_widening_op : public binary_op_t {
  private:
    const object_domain_t &m_left;
    const object_domain_t &m_right;
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
      BaseAbsDom singleton_base = s_single.project_singleton(key);

      if (m_is_join) {
        res_prod.first() |= singleton_base;
      } else {
        res_prod.first() = res_prod.first() || singleton_base;
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
      const odi_info_t &l_obj_info = l.first();
      const odi_value_t &l_odi_val = l.second();
      const odi_info_t &r_obj_info = r.first();
      const odi_value_t &r_odi_val = r.second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      map_value_t out_val;
      if (!crab_domain_params_man::get().singletons_in_base()) {
        // pairwise join
        out_val.second() =
            m_is_join ? l_odi_val | r_odi_val : l_odi_val || r_odi_val;
      } else {
        // NOTE: if we keep singleton object in the base, the join operation
        // should cover the case for singleton with non singleton
        if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
          out_val.second() = join_or_widen_singleton_with_non_singleton(
              key, m_left, r_odi_val);
        } else if (r_num_refs.is_one() &&
                   l_num_refs == small_range::oneOrMore()) {
          out_val.second() = join_or_widen_singleton_with_non_singleton(
              key, m_right, l_odi_val);
        } else {
          out_val.second() =
              m_is_join ? l_odi_val | r_odi_val : l_odi_val || r_odi_val;
        }
      }
      // After the join / widening, update the object info
      out_val.first() = odi_info_t(l_num_refs | r_num_refs,
                                   // Cache is not used
                                   boolean_value::get_false(),
                                   // Cache is not dirty
                                   boolean_value::get_false());
      return out_val;
    }

  protected:
    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key &key, const map_value_t &x, const map_value_t &y) override {
      map_value_t z = join_or_widen_odi(key, x, y);
      if (z.is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    }

    virtual bool default_is_absorbing() override { return true; }

  public:
    join_or_widening_op(const object_domain_t &left,
                        const object_domain_t &right, bool is_join)
        : m_left(left), m_right(right), m_is_join(is_join) {}
  }; // class join_op

  /// @brief A special class to compute meet / narrowing when merging two trees
  class meet_or_narrow_op : public binary_op_t {
  private:
    object_domain_t &m_left;
    object_domain_t &m_right;
    bool m_is_meet;

    /// @brief A special operation for the case where
    ///       state one contains a singleton object in the base domain
    ///       state two contains a non-singleton object in the odi map
    /// @param key the object id
    /// @param single the state containing singleton object
    /// @param non_single the state containing non-singleton
    void meet_or_narrow_non_singleton_with_singleton(
        const key_t &key, object_domain_t &single,
        const odi_value_t &non_single) const {
      const BaseAbsDom &summary = non_single.first();

      // Be careful of the following operation if type of base dom and object
      // dom are different
      single.meet_or_narrow_base(summary, m_is_meet);
    }

    /// @brief the meet or narrowing odi
    /// @param key the object id
    /// @param l the odi stored on the left state
    /// @param r the odi stored on the right state
    /// @return a new odi after meet or narrowing
    map_value_t meet_or_narrow_odi(const key_t &key, const map_value_t &l,
                                   const map_value_t &r) {
      const odi_info_t &l_obj_info = l.first();
      const odi_value_t &l_odi_val = l.second();
      const odi_info_t &r_obj_info = r.first();
      const odi_value_t &r_odi_val = r.second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      map_value_t out_val;
      if (!crab_domain_params_man::get().singletons_in_base()) {
        // pairwise join
        out_val.second() =
            m_is_meet ? l_odi_val & r_odi_val : l_odi_val && r_odi_val;
      } else {
        // NOTE: if we keep singleton object in the base, the meet operation
        // should cover the case for singleton with non singleton
        if (l_num_refs.is_one() && r_num_refs == small_range::oneOrMore()) {
          meet_or_narrow_non_singleton_with_singleton(key, m_left, r_odi_val);
        } else if (r_num_refs.is_one() &&
                   l_num_refs == small_range::oneOrMore()) {
          meet_or_narrow_non_singleton_with_singleton(key, m_right, l_odi_val);
        } else {
          out_val.second() =
              m_is_meet ? l_odi_val & r_odi_val : l_odi_val && r_odi_val;
        }
      }
      out_val.first() = odi_info_t(l_num_refs & r_num_refs,
                                   // Cache is not used
                                   boolean_value::get_false(),
                                   // Cache is not dirty
                                   boolean_value::get_false());
      return out_val;
    }

  protected:
    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key &key, const map_value_t &x, const map_value_t &y) override {
      map_value_t z = meet_or_narrow_odi(key, x, y);
      if (z.is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    }

    virtual bool default_is_absorbing() override { return false; }

  public:
    meet_or_narrow_op(object_domain_t &left, object_domain_t &right,
                      bool is_meet)
        : m_left(left), m_right(right), m_is_meet(is_meet) {}
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
      map_value_t z = x.widening_thresholds(y, m_ts);
      if (z.is_top()) {
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
      const odi_info_t &l_obj_info = x.first();
      const odi_value_t &l_odi_val = x.second();
      const odi_info_t &r_obj_info = y.first();
      const odi_value_t &r_odi_val = y.second();
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
        v.write(o);
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
                        const object_domain_t &right) const {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      join_or_widening_op jop(left, right, true);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Self join
  void compound_join(const odi_map_domain_t &o, const object_domain_t &left,
                     const object_domain_t &right) {
    if (is_bottom()) { // this is bot, assign this by o
      *this = o;
    } else if (o.is_bottom()) { // o is bot, nothing change
      return;
    } else {
      join_or_widening_op jop(left, right, true);
      bool is_bottom = false /*unused*/;
      apply_self_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
    }
  }

  // Meet
  odi_map_domain_t meet(const odi_map_domain_t &o, object_domain_t &left,
                        object_domain_t &right) const {
    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    } else {
      meet_or_narrow_op mop(left, right, true);
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
  void compound_meet(const odi_map_domain_t &o, object_domain_t &left,
                     object_domain_t &right) {
    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      *this = o;
    } else {
      meet_or_narrow_op mop(left, right, true);
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
                            const object_domain_t &right) const {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      join_or_widening_op wop(left, right, false);
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
  odi_map_domain_t narrowing(const odi_map_domain_t &o, object_domain_t &left,
                             object_domain_t &right) const {
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    } else {
      meet_or_narrow_op nop(left, right, false);
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
    if (!is_bottom()) {
      m_odi_map.remove(k);
    }
  }

  // NOTE: this operation is on key not on values
  void forget(const std::vector<key_t> &keys) {
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
  void set(const key_t &key, const map_value_t &v) {
    if (!is_bottom()) {
      if (v.is_bottom()) {
        set_to_bottom();
      } else if (v.is_top()) {
        m_odi_map.remove(key);
      } else {
        m_odi_map.insert(key, v);
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
  const map_value_t *find(const key_t &k) const {
    assert(!is_bottom());
    return m_odi_map.find(k);
  }

  /// @brief find a value by giving a key, otherwise, return a top value
  /// @param k object id
  /// @return a value if the key exists or a top if it does not
  map_value_t at(const key_t &k) const {
    if (is_bottom()) {
      map_value_t res;
      return res.make_bottom();
    } else {
      boost::optional<map_value_t> v = m_odi_map.lookup(k);
      if (v) {
        return *v;
      } else {
        map_value_t res;
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
    std::vector<key_t> keys;
    for (auto kv : m_odi_map) {
      const key_t &key = kv.first;
      keys.push_back(key);
    }
    return keys;
  }

  /// @brief a special log method for <SUM_DOM, CACHE_DOM, UF_DOM>
  /// @param o crab ostream
  /// @param prod a odi value
  void odi_val_write(crab_os &o, const odi_value_t &prod) const {
    o << "summary: ";
    prod.first().write(o);
    o << ", cache: ";
    prod.second().first().write(o);
    o << ", uf_fields: ";
    prod.second().second().write(o);
    o << ")";
  }

  /// @brief a special log method for <info, value>
  /// @param o crab ostream
  /// @param prod the value stored on the map
  void odi_write(crab_os &o, const map_value_t &prod) const {
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

  /**------------------ End Map APIs ------------------**/
}; // class odi_map_domain
} // end namespace object_domain_impl
} // end namespace domains
} // end namespace crab