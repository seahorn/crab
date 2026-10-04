#pragma once

#include <boost/optional.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/object/cow_domain_ref.hpp>
#include <crab/domains/object/object_info.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/domains/small_range.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {
namespace object_domain_impl {

/// @brief Equality functor consulted by the patricia tree to detect no-op
/// updates: when an insert/merge produces a value equal to the existing leaf,
/// the tree keeps the OLD node, preserving structural sharing across states.
/// @tparam MapValuePtr shared_ptr to the odi product <info, <sum, <cache, eq>>>
/// @note Semantic equality (mutual inclusion) would be exact but as costly as
/// the domain operation itself. Instead we decide equality in O(1) by identity:
/// two values are equal if they are the same allocation, or if their infos
/// coincide field-by-field and their three subdomain COW references share the
/// same underlying allocations (set()/merge always rewrap results in a fresh
/// outer shared_ptr, so the outer pointers alone almost never match).
template <class MapValuePtr> struct object_equal_to {
  bool operator()(const MapValuePtr &v1, const MapValuePtr &v2) const {
    if (v1 == v2) { // same allocation (both null included)
      return true;
    }
    if (!v1 || !v2) {
      return false;
    }
    return v1->first().equals(v2->first()) &&
           v1->second().first().same_absval(v2->second().first()) &&
           v1->second().second().first().same_absval(
               v2->second().second().first()) &&
           v1->second().second().second().same_absval(
               v2->second().second().second());
  }
};

/// @brief selector for the binary lattice operations, shared by the odi
/// map's merge operators and the object domain's general combine; widening
/// with thresholds is WIDENING plus a non-null thresholds pointer
enum class combine_kind { JOIN, WIDENING, MEET, NARROWING };

#define ODI_DOMAIN_SCOPED_STATS(NAME) CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)

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
  using usymb_t = typename eq_domain_value_t::class_id_t;
  using base_domain_ref_t =
      abstract_domain_ref<typename BaseAbsDom::variable_t, BaseAbsDom>;
  using eq_domain_ref_t =
      abstract_domain_ref<typename BaseAbsDom::variable_t, eq_domain_value_t>;

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

  // The map contains all objects' invariants including values and infos:
  //   object id --> <info, value>
  //   where info is a tuple of
  //     <ref count, obj init?, summary present?, cache used?, cache dirty?>
  //     plus two reduction flags (cache loaded by reg / stored from reg),
  //   and value is a tuple of <SUM_DOM, CACHE_DOM, EQ_DOM>.
  // No binding for an id means the object has not been seen (top).
  // A singleton object (refcount == 1) keeps its properties in CACHE_DOM;
  // once more allocations arrive, the cache is committed into SUM_DOM.
  // NOTE: the operations on this map depend on other domains, such as the
  // address domain (checking the MRU object in CACHE_DOM), the eq_domain
  // for registers (used for reduction), and the base domain (stored
  // properties for registers, references, etc.)
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
  //       ││ reference │  object  │ summary  │ cache │cache ││
  //       ││   count   │init flag │ presence │ used  │dirty ││
  //       ││           │          │   flag   │ flag  │flag  ││
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

  /// @brief A special class to compute join / widening / meet / narrowing
  /// when merging two trees
  class combine_op : public binary_op_t {
  private:
    const odi_map_domain_t &m_l_odi_map;
    const odi_map_domain_t &m_r_odi_map;
    combine_kind m_kind;
    // non-null only for widening with thresholds (m_kind == WIDENING)
    const thresholds<number_t> *m_ts;

    // join and widening grow the value; meet and narrowing shrink it
    bool is_grow() const {
      return m_kind == combine_kind::JOIN || m_kind == combine_kind::WIDENING;
    }

    /// @brief combine two odis pointwise, according to m_kind
    map_value_t combine_odi(const map_value_t &l, const map_value_t &r) {
      const odi_info_t &l_obj_info = l->first();
      const odi_value_t &c_l_odi_val = l->second();
      const odi_info_t &r_obj_info = r->first();
      const odi_value_t &c_r_odi_val = r->second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      boolean_value l_sum_presence = l_obj_info.sumpresence_val();
      boolean_value r_sum_presence = r_obj_info.sumpresence_val();
      const boolean_value l_is_init = l_obj_info.objinit_val();
      const boolean_value r_is_init = r_obj_info.objinit_val();
      map_raw_value_t out_val;

      // The copies below are cheap: each subdomain is copy-on-write.
      odi_value_t l_odi_val = odi_value_t(c_l_odi_val);
      odi_value_t r_odi_val = odi_value_t(c_r_odi_val);

      // commit cache if dirty
      m_l_odi_map.commit_cache_if_dirty(l_obj_info, l_odi_val);
      if (l_sum_presence.is_false()) {
        l_sum_presence = l_obj_info.cacheused_val();
      }

      m_r_odi_map.commit_cache_if_dirty(r_obj_info, r_odi_val);
      if (r_sum_presence.is_false()) {
        r_sum_presence = r_obj_info.cacheused_val();
      }

      switch (m_kind) {
      case combine_kind::JOIN: // pairwise join
        out_val.second() = l_odi_val | r_odi_val;
        break;
      case combine_kind::WIDENING:
        if (m_ts == nullptr) { // pairwise widening
          out_val.second() = l_odi_val || r_odi_val;
        } else {
          // pairwise widening with thresholds, component-wise:
          // basic_domain_product2 has no widening_thresholds, the COW
          // wrappers forward it, and the equality component's
          // widening_thresholds is its join (thresholds are meaningless for
          // equalities)
          odi_value_t out_v;
          odi_map_domain_t::object_sum_val(out_v) =
              odi_map_domain_t::object_sum_val(l_odi_val).widening_thresholds(
                  odi_map_domain_t::object_sum_val(r_odi_val), *m_ts);
          odi_map_domain_t::object_cache_val(out_v) =
              odi_map_domain_t::object_cache_val(l_odi_val).widening_thresholds(
                  odi_map_domain_t::object_cache_val(r_odi_val), *m_ts);
          odi_map_domain_t::object_eq_val(out_v) =
              odi_map_domain_t::object_eq_val(l_odi_val).widening_thresholds(
                  odi_map_domain_t::object_eq_val(r_odi_val), *m_ts);
          out_val.second() = std::move(out_v);
        }
        break;
      case combine_kind::MEET: // pairwise meet
        out_val.second() = l_odi_val & r_odi_val;
        break;
      case combine_kind::NARROWING: // pairwise narrowing
        out_val.second() = l_odi_val && r_odi_val;
        break;
      }

      // After the combine, update the object info
      const bool grow = is_grow();
      out_val.first() =
          odi_info_t(grow ? l_num_refs | r_num_refs : l_num_refs & r_num_refs,
                     grow ? l_is_init | r_is_init : l_is_init & r_is_init,
                     grow ? l_sum_presence | r_sum_presence
                          : l_sum_presence & r_sum_presence,
                     /*cache_used=*/boolean_value::get_false(),
                     /*cache_dirty=*/boolean_value::get_false(),
                     /*is_loaded=*/false, /*is_stored=*/false);
      return std::make_shared<map_raw_value_t>(out_val);
    }

  protected:
    virtual std::pair<bool, boost::optional<map_value_t>>
    apply(const Key &key, const map_value_t &x, const map_value_t &y) override {
      if (x == y) {
        return {false, boost::optional<map_value_t>(x)};
      }
      map_value_t z = combine_odi(x, y);
      if (!is_grow() && z->is_bottom()) {
        // a contradictory per-object meet makes the whole map bottom
        // (patricia_trees.hpp: first == true means the result is bottom;
        // a join/widening of two non-bottom odis cannot be bottom)
        return {true, boost::optional<map_value_t>()};
      } else if (z->is_top()) {
        return {false, boost::optional<map_value_t>()};
      } else {
        return {false, boost::optional<map_value_t>(z)};
      }
    }

    // join/widening: an entry missing on one side means top, and top
    // absorbs; meet/narrowing: the present entry is kept
    virtual bool default_is_absorbing() override { return is_grow(); }

  public:
    combine_op(const odi_map_domain_t &left, const odi_map_domain_t &right,
               combine_kind kind, const thresholds<number_t> *ts = nullptr)
        : m_l_odi_map(left), m_r_odi_map(right), m_kind(kind), m_ts(ts) {}
  }; // class combine_op

  /// @brief A special class to perform the inclusion test (<=) when merging
  /// two trees
  class inclusion_test_op : public partial_order_t {
  private:
    const odi_map_domain_t &m_l_odi_map;
    const odi_map_domain_t &m_r_odi_map;

  protected:
    // Only refcounts and summaries decide the inclusion. The caches are
    // deliberately NOT compared against each other: each state's cache
    // describes that state's own MRU object, and two states may have
    // cached DIFFERENT objects (the identity lives in the address domain,
    // invisible here). Dirty caches are folded into the summary copies
    // first; the summary is the only identity-agnostic component.
    virtual bool leq(const map_value_t &x, const map_value_t &y) override {
      if (x == y) {
        return true;
      }
      const odi_info_t &l_obj_info = x->first();
      const odi_value_t &c_l_odi_val = x->second();
      const small_range &l_num_refs = l_obj_info.refcount_val();
      const odi_info_t &r_obj_info = y->first();
      const odi_value_t &c_r_odi_val = y->second();
      const small_range &r_num_refs = r_obj_info.refcount_val();
      //
      odi_value_t l_odi_val = odi_value_t(c_l_odi_val);
      odi_value_t r_odi_val = odi_value_t(c_r_odi_val);
      m_l_odi_map.commit_cache_if_dirty(l_obj_info, l_odi_val);
      m_r_odi_map.commit_cache_if_dirty(r_obj_info, r_odi_val);
      bool res = l_num_refs <= r_num_refs;
      res &= l_odi_val.first() <= r_odi_val.first();
      res &= l_odi_val.second().first() <= r_odi_val.second().first();
      return res;
    }

    virtual bool default_is_top() override { return true; }

  public:
    inclusion_test_op(const odi_map_domain_t &left,
                      const odi_map_domain_t &right)
        : m_l_odi_map(left), m_r_odi_map(right) {}
  }; // class inclusion_test_op

  /// @brief apply an operation when merging two trees
  /// @param o the binary operation for merging two values with the same key
  /// @param t1 a tree
  /// @param t2 a tree
  /// @param is_bottom a bool value indicates whether the result is bottom
  /// @return the merged tree
  static patricia_tree_t apply_operation(binary_op_t &o, patricia_tree_t t1,
                                         const patricia_tree_t &t2,
                                         bool &is_bottom) {
    is_bottom = t1.merge_with(t2, o);
    return t1;
  }

  /// @brief apply an operation when merging two trees (self-apply)
  /// @param o the binary operation for merging two values with the same key
  /// @param t1 a tree, the merging result will update this tree
  /// @param t2 a tree
  /// @param is_bottom a bool value indicates whether the result is bottom
  static void apply_self_operation(binary_op_t &o, patricia_tree_t &t1,
                                   const patricia_tree_t &t2, bool &is_bottom) {
    is_bottom = t1.merge_with(t2, o);
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
      inclusion_test_op leq_op(*this, o);
      return m_odi_map.leq(o.m_odi_map, leq_op);
    }
  }

  // Join
  odi_map_domain_t join(const odi_map_domain_t &o) const {
    ODI_DOMAIN_SCOPED_STATS(".join");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() || o.is_top()) {
      return odi_map_domain_t();
    } else {
      combine_op jop(*this, o, combine_kind::JOIN);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Self join
  void compound_join(const odi_map_domain_t &o) {
    ODI_DOMAIN_SCOPED_STATS(".join");

    if (is_bottom()) { // this is bot, assign this by o
      *this = o;
    } else if (o.is_bottom()) { // o is bot, nothing change
    } else if (is_top() || o.is_top()) {
      set_to_top();
    } else {
      combine_op jop(*this, o, combine_kind::JOIN);
      bool is_bottom = false /*unused*/;
      apply_self_operation(jop, m_odi_map, o.m_odi_map, is_bottom);
    }
  }

  // Meet
  odi_map_domain_t meet(const odi_map_domain_t &o) const {
    ODI_DOMAIN_SCOPED_STATS(".meet");

    if (is_bottom() || o.is_top()) { // bot & o or this & top, return this
      return *this;
    } else if (o.is_bottom() || is_top()) { // this & bot or top & o, return o
      return o;
    } else {
      combine_op mop(*this, o, combine_kind::MEET);
      bool is_bottom = false;
      patricia_tree_t res =
          apply_operation(mop, m_odi_map, o.m_odi_map, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return odi_map_domain_t(std::move(res));
      }
    }
  }

  // Widening
  odi_map_domain_t widening(const odi_map_domain_t &o) const {
    ODI_DOMAIN_SCOPED_STATS(".widening");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      combine_op wop(*this, o, combine_kind::WIDENING);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(wop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Widening with thresholds
  odi_map_domain_t widening_thresholds(const odi_map_domain_t &o,
                                       const thresholds<number_t> &ts) const {
    ODI_DOMAIN_SCOPED_STATS(".widening");

    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      combine_op wop(*this, o, combine_kind::WIDENING, &ts);
      bool is_bottom = false /*unused*/;
      patricia_tree_t res =
          apply_operation(wop, m_odi_map, o.m_odi_map, is_bottom);
      return odi_map_domain_t(std::move(res));
    }
  }

  // Narrowing
  odi_map_domain_t narrowing(const odi_map_domain_t &o) const {
    ODI_DOMAIN_SCOPED_STATS(".narrowing");
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    } else {
      combine_op nop(*this, o, combine_kind::NARROWING);
      bool is_bottom = false;
      patricia_tree_t res =
          apply_operation(nop, m_odi_map, o.m_odi_map, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return odi_map_domain_t(std::move(res));
      }
    }
  }

  /// @brief general binary lattice operation, selected by \p kind
  /// @param o the other odi map
  /// @param kind which operation to perform
  /// @param ts non-null only for widening with thresholds
  odi_map_domain_t combine(const odi_map_domain_t &o, combine_kind kind,
                           const thresholds<number_t> *ts = nullptr) const {
    if (kind == combine_kind::JOIN) {
      return join(o);
    } else if (kind == combine_kind::WIDENING) {
      return ts ? widening_thresholds(o, *ts) : widening(o);
    } else if (kind == combine_kind::MEET) {
      return meet(o);
    } else {
      return narrowing(o);
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

  /// @brief update odi map by a new <obj_id, <info, value>>
  /// @param key object id
  /// @param v the new value by a product of <info, value>
  void set(const key_t &key, const map_raw_value_t &v) {
    ODI_DOMAIN_SCOPED_STATS(".set");
    if (!is_bottom()) {
      if (v.is_bottom()) {
        set_to_bottom();
      } else if (v.is_top()) {
        m_odi_map.remove(key);
      } else {
        map_value_t v_ptr = std::make_shared<map_raw_value_t>(v);
        m_odi_map.insert(key, v_ptr);
      }
    }
  }

  /// @brief find a value by giving a key
  /// @param k object id
  /// @return a pointer to const if the value exists
  /// The value in the map is a constant object which should not be modified.
  /// @warning the returned pointer (and any reference derived from it) is
  /// invalidated by any set()/operator-= on that key: set() replaces the
  /// tree leaf and, when this state is the entry's sole owner, frees the
  /// old value. Copy the scalars you need BEFORE mutating the map.
  const map_raw_value_t *find(const key_t &k) const {
    assert(!is_bottom());
    auto val_ptr = m_odi_map.find(k);
    if (val_ptr) {
      return val_ptr->get();
    } else {
      return nullptr;
    }
  }

  /// @brief a special log method for <SUM_DOM, CACHE_DOM, EQ_DOM>
  /// @param o crab ostream
  /// @param prod an odi value
  void odi_val_write(crab_os &o, const odi_value_t &prod) const {
    o << "summary: ";
    prod.first().write(o);
    o << ", cache: ";
    prod.second().first().write(o);
    o << ", eq_fields: ";
    prod.second().second().write(o);
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
    } else {
      odi_val_write(o, obj_val);
    }
    o << ")";
  }

  void odi_write(crab_os &o, const map_value_t &prod) const {
    odi_write(o, *prod);
  }

  /**------------------ End Map APIs ------------------**/

  /**------------------ Begin Cache APIs ------------------**/
  // commit cache contents into summary
  void commit_cache(summary_domain_t &summary,
                    const cache_domain_t &cache) const {
    summary |= cache;
  }

  // update cache contents from summary to cache
  void update_cache(const summary_domain_t &summary,
                    cache_domain_t &cache) const {
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
      bool is_summary_absence = obj_info.sumpresence_val().is_false();
      const small_range &num_refs = obj_info.refcount_val();
      if (is_summary_absence || num_refs.is_one()) {
        // Even if an object is not singleton, the summary is still absent
        // if no cache flush has occurred
        sum = cache;
      } else {
        commit_cache(sum, cache);
      }
    }
    // skip the reset when already top: set_to_top on a shared value
    // allocates a replacement, and a fresh pointer defeats the
    // leaf-sharing equality (object_equal_to) on later merges
    if (!cache.is_top()) {
      cache.set_to_top();
    }
    if (!eq_fld.is_top()) {
      eq_fld.set_to_top();
    }
  }

  bool invalidate_cache_if_miss(
      key_t &key, object_domain_t &abs_state,
      ghost_variables_eq_t &&rgn_eq_gvars, boost::optional<usymb_t> &reg_symb,
      boost::optional<std::pair<usymb_t, usymb_t>> &offset_size_symb,
      bool is_ref_mru, bool is_store) {
    const map_raw_value_t *obj_prod_ref = find(key);

    bool update_new_mru = false;
    if (!obj_prod_ref) {
      CRAB_ERROR(domain_name(), "::", __func__, ": accessing ", key,
                 " is not found on the odi map");
    }

    odi_info_t out_obj_info = obj_prod_ref->first();
    auto cache_used = out_obj_info.cacheused_val();

    // NOTE: copying the object value is required
    odi_value_t out_prod = obj_prod_ref->second();
    eq_domain_t &eq_flds_dom = out_prod.second().second();

    /* Cache missed condition:
        cache is empty or current reference does not refer to the mru object
    */
    if (cache_used.is_false() || is_ref_mru == false) { // cache is missed
      // if (out_obj_info.refcount_val() == small_range::oneOrMore()) {
      // Step1: commit cache if the cache is dirty
      abs_state.commit_cache_if_dirty(out_prod, out_obj_info, key);
      if (out_obj_info.sumpresence_val().is_false()) {
        out_obj_info.sumpresence_val() = cache_used;
      }
      // }
      // Step2: update cache for new MRU object
      update_cache(out_prod.first(), out_prod.second().first());
      // Step3: update address dom and object info
      update_new_mru = true;

      out_obj_info =
          odi_info_t(out_obj_info.refcount_val(), out_obj_info.objinit_val(),
                     out_obj_info.sumpresence_val(),
                     // Cache is used
                     boolean_value::get_true(),
                     // Cache is not dirty
                     boolean_value::get_false(),
                     // Cache will be loaded to a reg?
                     is_store == false,
                     // Cache will be stored from a reg?
                     is_store == true && reg_symb != boost::none);
    } else {          // cache is still hit
      if (is_store) { // call from ref_store
        // if cache needs to be stored from a reg1, but it used
        // to be loaded by a reg2, perform reduction from object to base
        // and then update the cache, without losing precision and while
        // remaining sound.
        if (reg_symb == boost::none && out_obj_info.cache_reg_loaded_val()) {
          out_obj_info.cache_reg_loaded_val() = false;
          abs_state.apply_reduction_from_object_to_base(out_prod, key);
        }
        // delay the reduction: keep the equality rgn_w == reg for the
        // pending store
        out_obj_info.cache_reg_stored_val() =
            is_store == true && reg_symb != boost::none;
      } else { // call from ref_load
        // if cache needs to be loaded by a reg1, but it used to be stored
        // from a reg2, perform reduction from base to object
        if (out_obj_info.cache_reg_stored_val()) {
          out_obj_info.cache_reg_stored_val() = false;
          if (out_obj_info.cache_reg_loaded_val()) {
            abs_state.apply_reduction_from_object_to_base(out_prod, key);
          }
          abs_state.apply_reduction_from_base_to_object(out_prod, key);
        }
        out_obj_info.cache_reg_loaded_val() = true;
      }
    }
    if (!is_store ||
        reg_symb != boost::none) { // only for load_ref and store_ref by reg
      // Step4: update field == reg, the equality represents either:
      //        reg := load_ref(ref, field), or
      //        store_ref(ref, field, reg)
      if (!is_store) { // it means caller is from load_ref
        // assigning the same symbol as rgn
        reg_symb =
            abs_state.get_symbol_or_fresh(*eq_flds_dom, rgn_eq_gvars.get_var());
      }
      (*eq_flds_dom).set(rgn_eq_gvars.get_var(), *reg_symb);
      if (rgn_eq_gvars.has_offset_and_size()) {
        if (offset_size_symb == boost::none) {
          offset_size_symb = std::make_pair(
              abs_state.get_symbol_or_fresh(
                  *eq_flds_dom,
                  rgn_eq_gvars.get_offset_and_size().get_offset()),
              abs_state.get_symbol_or_fresh(
                  *eq_flds_dom, rgn_eq_gvars.get_offset_and_size().get_size()));
        }
        (*eq_flds_dom)
            .set(rgn_eq_gvars.get_offset_and_size().get_offset(),
                 std::get<0>(*offset_size_symb));
        (*eq_flds_dom)
            .set(rgn_eq_gvars.get_offset_and_size().get_size(),
                 std::get<1>(*offset_size_symb));
      }
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

  static summary_domain_t &object_sum_val(odi_value_t &odi_val) {
    return odi_val.first();
  }

  static base_domain_t &object_sum_raw_val(odi_value_t &odi_val) {
    return *odi_val.first();
  }

  // const overload: reads through the COW reference without detaching
  static const base_domain_t &object_sum_raw_val(const odi_value_t &odi_val) {
    return *odi_val.first();
  }

  static const cache_domain_t &object_cache_val(const odi_value_t &odi_val) {
    return odi_val.second().first();
  }

  static cache_domain_t &object_cache_val(odi_value_t &odi_val) {
    return odi_val.second().first();
  }

  static base_domain_t &object_cache_raw_val(odi_value_t &odi_val) {
    return *odi_val.second().first();
  }

  // const overload: reads through the COW reference without detaching
  static const base_domain_t &object_cache_raw_val(const odi_value_t &odi_val) {
    return *odi_val.second().first();
  }

  static const eq_domain_t &object_eq_val(const odi_value_t &odi_val) {
    return odi_val.second().second();
  }

  static eq_domain_t &object_eq_val(odi_value_t &odi_val) {
    return odi_val.second().second();
  }

  static eq_domain_value_t &object_eq_raw_val(odi_value_t &odi_val) {
    return *odi_val.second().second();
  }

  // const overload: reads through the COW reference without detaching
  static const eq_domain_value_t &
  object_eq_raw_val(const odi_value_t &odi_val) {
    return *odi_val.second().second();
  }

  /**------------------ End Product Domain APIs ------------------**/
}; // class odi_map_domain
} // end namespace object_domain_impl
} // end namespace domains
} // end namespace crab