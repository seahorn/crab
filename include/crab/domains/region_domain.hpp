#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/varname_factory.hpp>

#include <algorithm>
#include <functional>
#include <unordered_map>

namespace crab {
namespace domains {

/*
 * Class that represents the range of intervals
 * {[0,0], [1,1], [0,1], [0,+oo], [1,+oo]}
 */
class small_range {

  using kind_t = enum {
    Bottom,
    ExactlyZero,
    ExactlyOne,
    ZeroOrOne,
    ZeroOrMore,
    OneOrMore
  };

  kind_t m_value;

  small_range(kind_t v) : m_value(v){};

  /*
      [0,0] | [0,0] = [0,0]
      [0,0] | [1,1] = [0,1]
      [0,0] | [0,1] = [0,1]
      [0,0] | [0,+oo] = [0,+oo]
      [0,0] | [1,+oo] = [0,+oo]
  */
  small_range join_zero_with(const small_range &other) const {
    if (other.m_value == ExactlyZero) {
      return other;
    } else if (other.m_value == ExactlyOne) {
      return small_range(ZeroOrOne);
    } else if (other.m_value == ZeroOrOne) {
      return other;
    } else {
      return small_range(ZeroOrMore);
    }
  }

  /*
      [1,1] | [0,0] = [0,1]
      [1,1] | [1,1] = [1,+oo]
      [1,1] | [0,1] = [0,+oo]
      [1,1] | [0,+oo] = [0,+oo]
      [1,1] | [1,+oo] = [1,+oo]
  */
  small_range join_one_with(const small_range &other) const {
    if (other.m_value == ExactlyZero) {
      return small_range(ZeroOrOne);
    } else if (other.m_value == ExactlyOne) {
      return small_range(OneOrMore);
    } else if (other.m_value == OneOrMore) {
      return other;
    } else {
      return small_range(ZeroOrMore);
    }
  }

  /*
      [0,1] | [0,0] = [0,1]
      [0,1] | [1,1] = [0,+oo]
      [0,1] | [0,1] = [0,+oo]
      [0,1] | [0,+oo] = [0,+oo]
      [0,1] | [1,+oo] = [0,+oo]
  */
  small_range join_zero_or_one_with(const small_range &other) const {
    if (other.m_value == ExactlyZero) {
      return small_range(ZeroOrOne);
    } else {
      return small_range(ZeroOrMore);
    }
  }

  /*
      [1,+oo] | [0,0] = [0,+oo]
      [1,+oo] | [1,1] = [1,+oo]
      [1,+oo] | [0,1] = [0,+oo]
      [1,+oo] | [0,+oo] = [0,+oo]
      [1,+oo] | [1,+oo] = [1,+oo]
  */
  small_range join_one_or_more_with(const small_range &other) const {
    if (other.m_value == ExactlyOne || other.m_value == OneOrMore) {
      return small_range(OneOrMore);
    } else {
      return small_range(ZeroOrMore);
    }
  }

  /*
      [0,0] & [0,0] = [0,0]
      [0,0] & [1,1] = _|_
      [0,0] & [0,1] = [0,0]
      [0,0] & [0,+oo] = [0,0]
      [0,0] & [1,+oo] = _|_
  */

  small_range meet_zero_with(const small_range &other) const {
    assert(is_zero());
    if (other.m_value == ExactlyOne || other.m_value == OneOrMore) {
      return small_range::bottom();
    } else {
      return *this;
    }
  }

  /*
      [1,1] & [0,0] = _|_
      [1,1] & [1,1] = [1,1]
      [1,1] & [0,1] = [1,1]
      [1,1] & [0,+oo] = [1,1]
      [1,1] & [1,+oo] = [1,1]
  */
  small_range meet_one_with(const small_range &other) const {
    assert(is_one());
    if (other.m_value == ExactlyZero) {
      return small_range::bottom();
    } else {
      return *this;
    }
  }

  /*
      [0,1] & [0,0] = [0,0]
      [0,1] & [1,1] = [1,1]
      [0,1] & [0,1] = [0,1]
      [0,1] & [0,+oo] = [0,1]
      [0,1] & [1,+oo] = [1,1]
  */
  small_range meet_zero_or_one_with(const small_range &other) const {
    assert(m_value == ZeroOrOne);
    if (other.is_zero() || other.is_one()) {
      return other;
    } else if (other.m_value == OneOrMore) {
      return one();
    } else {
      return *this;
    }
  }

  /*
      [1,+oo] & [0,0] = _|_
      [1,+oo] & [1,1] = [1,1]
      [1,+oo] & [0,1] = [1,1]
      [1,+oo] & [0,+oo] = [1,+oo]
      [1,+oo] & [1,+oo] = [1,+oo]
  */

  small_range meet_one_or_more_with(const small_range &other) const {
    if (other.is_zero()) {
      return small_range::bottom();
    } else if (other.is_one()) {
      return other;
    } else if (other.m_value == ZeroOrOne) {
      return one();
    } else {
      assert(other.is_top() || other.m_value == OneOrMore);
      return *this;
    }
  }

public:
  small_range() : m_value(ZeroOrMore) {}

  static small_range bottom() { return small_range(Bottom); }
  static small_range top() { return small_range(ZeroOrMore); }
  static small_range zero() { return small_range(ExactlyZero); }
  static small_range one() { return small_range(ExactlyOne); }
  static small_range oneOrMore() { return small_range(OneOrMore); }

  small_range(const small_range &other) = default;
  small_range(small_range &&other) = default;
  small_range &operator=(const small_range &other) = default;
  small_range &operator=(small_range &&other) = default;

  bool is_bottom() const { return (m_value == Bottom); }

  bool is_top() const { return (m_value == ZeroOrMore); }

  bool is_zero() const { return (m_value == ExactlyZero); }

  bool is_one() const { return (m_value == ExactlyOne); }

  /*
     [0,+oo]
       |   \
       |    \
      [0,1] [1,+oo]
       / \ /
      0   1
   */
  bool operator<=(small_range other) const {
    if (m_value == other.m_value) {
      return true;
    } else if (m_value == Bottom || other.m_value == ZeroOrMore) {
      return true;
    } else if (m_value == ExactlyZero) {
      return other.m_value != ExactlyOne && other.m_value != OneOrMore;
    } else if (m_value == ExactlyOne) {
      return other.m_value != ExactlyZero;
    } else if (m_value == ZeroOrOne || m_value == OneOrMore) {
      return other.m_value == ZeroOrMore;
    } else if (m_value == ZeroOrMore) {
      assert(other.m_value != ZeroOrMore);
      return false;
    }
    // should be unreachable
    return false;
  }

  bool operator==(small_range other) const { return m_value == other.m_value; }

  small_range operator|(small_range other) const {
    if (is_bottom()) {
      return other;
    } else if (other.is_bottom()) {
      return *this;
    } else if (is_zero()) {
      return join_zero_with(other);
    } else if (other.is_zero()) {
      return join_zero_with(*this);
    } else if (is_one()) {
      return join_one_with(other);
    } else if (other.is_one()) {
      return join_one_with(*this);
    } else if (m_value == ZeroOrOne) {
      return join_zero_or_one_with(other);
    } else if (other.m_value == ZeroOrOne) {
      return join_zero_or_one_with(*this);
    } else if (m_value == OneOrMore) {
      return join_one_or_more_with(other);
    } else if (other.m_value == OneOrMore) {
      return join_one_or_more_with(*this);
    } else {
      return small_range(ZeroOrMore);
    }
  }

  small_range operator||(small_range other) const { return *this | other; }

  small_range operator&(small_range other) const {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    }

    if (is_zero()) {
      return meet_zero_with(other);
    } else if (other.is_zero()) {
      return other.meet_zero_with(*this);
    } else if (is_one()) {
      return meet_one_with(other);
    } else if (other.is_one()) {
      return other.meet_one_with(*this);
    } else if (m_value == ZeroOrOne) {
      return meet_zero_or_one_with(other);
    } else if (other.m_value == ZeroOrOne) {
      return other.meet_zero_or_one_with(*this);
    } else if (m_value == OneOrMore) {
      return meet_one_or_more_with(other);
    } else if (other.m_value == OneOrMore) {
      return other.meet_one_or_more_with(*this);
    } else { // unreachable because top cases handled above
      CRAB_ERROR("unexpected small_range::meet operands");
    }
  }

  small_range operator&&(small_range other) const { return *this & other; }

  small_range increment(void) {
    if (!is_bottom()) {
      if (m_value == ExactlyZero) {
        m_value = ExactlyOne;
      } else if (m_value == ExactlyOne || m_value == ZeroOrMore ||
                 m_value == ZeroOrOne || m_value == OneOrMore) {
        m_value = OneOrMore;
      } else {
        CRAB_ERROR("small_range::increment unreachable");
      }
    }
    return *this;
  }

  void write(crab_os &o) const {
    switch (m_value) {
    case Bottom:
      o << "_|_";
      break;
    case ExactlyZero:
      o << "[0,0]";
      break;
    case ExactlyOne:
      o << "[1,1]";
      break;
    case ZeroOrOne:
      o << "[0,1]";
      break;
    case ZeroOrMore:
      o << "[0,+oo]";
      break;
    case OneOrMore:
      o << "[1,+oo]";
      break;
    default:
      CRAB_ERROR("unexpected small_range value");
    }
  }

  friend crab_os &operator<<(crab_os &o, const small_range &v) {
    v.write(o);
    return o;
  }
}; // end class small_range

namespace region_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;
  using base_abstract_domain_t = BaseAbsDom;
  using base_varname_t = typename BaseAbsDom::varname_t;

  // Reason about deallocation
  enum { deallocation = 1 };
  // Improve precision but it might not be sound
  enum { refine_uninitialized_regions = 0 };

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
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_variable_vector_t =
      typename base_abstract_domain_t::variable_vector_t;
  using base_variable_t = typename base_abstract_domain_t::variable_t;
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

  bool m_is_bottom; // special symbol for bottom
  // To create synthetic variables for the base domain.
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
  // Abstract domain to count how many addresses are in a region.
  // This allows us to decide when strong update is sound: only if
  // one address per region (i.e., singleton).
  rgn_counting_env_t m_rgn_counting_dom;
  // Whether the summarized region has some object (or slice) that has
  // been initialized
  rgn_bool_env_t m_rgn_init_dom;
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
  //   ref_init(rgn): create a singleton equivalence class with rgn.
  //   ref_gep(reg1, rgn1, ref2, rgn2, o): join together the
  //                                       equivalence classes of rgn1
  //                                       and rgn2.
  rgn_dealloc_t m_rgn_dealloc_dom;
  // The base abstract domain: all the heavy lifting is done here.
  // m_base_dom does not have any variable of reference type.
  base_abstract_domain_t m_base_dom;

  region_domain(base_varname_allocator_t &&alloc, var_map_t &&var_map,
                rev_var_map_t &&rev_var_map,
                rgn_counting_env_t &&rgn_counting_dom,
                rgn_bool_env_t &&rgn_init_dom, rgn_dealloc_t &&rgn_dealloc_dom,
                base_abstract_domain_t &&base_dom)
      : m_is_bottom(base_dom.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
        m_base_dom(std::move(base_dom)) {}

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

    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);

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
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);

    CRAB_LOG("region", crab::outs() << *this << "\n");
  }

  region_domain_t do_join_or_widening(const region_domain_t &left,
                                      const region_domain_t &right,
                                      base_dom_binop_t base_dom_op) const {

    // rgn_counting_dom does not require common renaming
    // the domain is finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    // rgn_init_dom does not require common renaming
    // the domain is finite
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom | right.m_rgn_init_dom);
    // rgn_dealloc_dom does not require common renaming
    // the domain is finite
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);

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
        std::move(out_rgn_init_dom), std::move(out_rgn_dealloc_dom),
        std::move(out_base_dom));
    return res;
  }

  region_domain_t do_meet_or_narrowing(const region_domain_t &left,
                                       const region_domain_t &right,
                                       base_dom_binop_t base_dom_op) const {

    // these two domains are finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom &
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom & right.m_rgn_init_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom &
                                      right.m_rgn_dealloc_dom);

    // This shouldn't happen but just in case ...
    if (out_rgn_counting_dom.is_bottom() || out_rgn_init_dom.is_bottom() ||
        out_rgn_dealloc_dom.is_bottom()) {
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
        std::move(out_rgn_init_dom), std::move(out_rgn_dealloc_dom),
        std::move(out_base_dom));
    return res;
  }

  // Create a fresh variable in the base domain to shadow v
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
        CRAB_ERROR(domain_name() + "::make_base_variable: unreachable");
      }
      base_variable_t bv(var_allocator.next(), ty, bitwidth);
      return bv;
    } else {
      base_variable_t bv(var_allocator.next(), v.get_type());
      return bv;
    }
  }

  // Get the synthetic variable in the base domain that represents v
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

  // Return true if ref is definitely null. Ask the base numerical
  // domain for that.
  bool is_null_ref(const variable_t &ref) {
    if (!ref.get_type().is_reference()) {
      return false;
    }
    if (auto var_opt = get_var(ref)) {
      interval_t ival = m_base_dom[*var_opt];
      boost::optional<number_t> x = ival.lb().number();
      boost::optional<number_t> y = ival.ub().number();
      number_t zero(0);
      return (x && y && *x == zero && *y == zero);
    }
    return false;
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
        m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom), m_base_dom(o.m_base_dom) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  region_domain(region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
        m_rgn_dealloc_dom(std::move(o.m_rgn_dealloc_dom)),
        m_base_dom(std::move(o.m_base_dom)) {}

  region_domain_t &operator=(const region_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_rgn_init_dom = o.m_rgn_init_dom;
      m_rgn_dealloc_dom = o.m_rgn_dealloc_dom;
      m_base_dom = o.m_base_dom;
    }
    return *this;
  }

  region_domain_t &operator=(region_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_alloc = std::move(o.m_alloc);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_rgn_init_dom = std::move(o.m_rgn_init_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_base_dom = std::move(o.m_base_dom);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    return (!is_bottom() && m_base_dom.is_top() && m_rgn_counting_dom.is_top());
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

    if (Params::deallocation) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result5=0\n";);
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
    CRAB_LOG("region-leq", crab::outs() << "Result6=" << res << "\n";);
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
    if (Params::refine_uninitialized_regions) {
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
    if (Params::refine_uninitialized_regions) {
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
        if (Params::deallocation) {
          m_rgn_dealloc_dom.forget(v);
        }
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

    if (Params::deallocation) {
      // No deallocated objects in the region
      m_rgn_dealloc_dom.set(rgn, boolean_value::get_false());
    }

    if (!rgn.get_type().is_unknown_region()) {
      // Assign a synthetic variable to rgn for modeling its content.
      base_variable_t v = make_base_variable(m_alloc, rgn);
    }
    CRAB_LOG("region", crab::outs() << "After region_init(" << rgn
                                    << ")=" << *this << "\n";);
  }

  void region_copy(const variable_t &lhs_rgn,
                   const variable_t &rhs_rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_copy");
    crab::ScopedCrabStats __st__(domain_name() + ".region_copy");

    if (is_bottom()) {
      return;
    }

    if (!lhs_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn,
                 " must have a region type");
    }
    if (!rhs_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_copy ", rhs_rgn,
                 " must have a region type");
    }

    m_rgn_counting_dom.set(lhs_rgn, m_rgn_counting_dom[rhs_rgn]);
    m_rgn_init_dom.set(lhs_rgn, m_rgn_init_dom[rhs_rgn]);
    if (Params::deallocation) {
      m_rgn_dealloc_dom.forget(lhs_rgn);
      m_rgn_dealloc_dom.expand(rhs_rgn, lhs_rgn);
    }

    if (lhs_rgn.get_type().is_unknown_region() ||
        rhs_rgn.get_type().is_unknown_region()) {
      return;
    }

    if (lhs_rgn.get_type() != rhs_rgn.get_type()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn, ":=", rhs_rgn,
                 " with different types");
    }

    const base_variable_t &base_lhs = rename_var(lhs_rgn);
    const base_variable_t &base_rhs = rename_var(rhs_rgn);
    m_base_dom.forget({base_lhs});
    m_base_dom.expand(base_rhs, base_lhs);

    // auto ty = lhs_rgn.get_type();
    // if (ty.is_bool_region()) {
    //   m_base_dom.assign_bool_var(base_lhs, base_rhs, false);
    // } else if (ty.is_integer_region() || ty.is_real_region() ||
    // ty.is_reference_region()) {
    //   m_base_dom.assign(base_lhs, base_rhs);
    // } else if (ty.is_bool_array_region() ||
    // 	       ty.is_int_array_region() ||
    // 	       ty.is_real_array_region()) {
    //   m_base_dom.array_assign(base_lhs, base_rhs);
    // } else {
    //   CRAB_ERROR(domain_name() + "::region_copy with unexpected region
    //   type");
    // }
  }

  // Create a new reference ref to region rgn.
  void ref_make(const variable_t &ref, const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    if (is_bottom()) {
      return;
    }

    if (!ref.get_type().is_reference()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_make: ", ref,
                                   " must be a reference"););
      return;
    }

    // Update region counting
    auto num_refs = m_rgn_counting_dom[rgn];
    m_rgn_counting_dom.set(rgn, num_refs.increment());

    // Assign a base domain variable to ref
    rename_var(ref);

    CRAB_LOG("region", crab::outs() << "After ref_make(" << ref << "," << rgn
                                    << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within rgn. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    if (is_bottom()) {
      return;
    }

    const base_variable_t &base_res = rename_var(res);

    // TODO: handle unknown regions
    if (rgn.get_type().is_unknown_region()) {
      CRAB_LOG(
          "region-load",
          CRAB_WARN("region_domain::ref_load: skip unknown region ", rgn););
      m_base_dom -= base_res;
      return;
    }

    // Non-scalar regions should be handled by ref_load_from_array
    if (!rgn.get_type().is_scalar_region()) {
      CRAB_LOG("region-load", CRAB_WARN("region_domain::ref_load: ", rgn,
                                        " must be a scalar region"););
      m_base_dom -= base_res;
      return;
    } else {
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
    }

    if (!ref.get_type().is_reference()) {
      // This shouldn't happen
      CRAB_LOG("region-load", CRAB_WARN("region_domain::ref_load: ", ref,
                                        " must be a reference"););
      m_base_dom -= base_res;
      return;
    }

    if (is_null_ref(ref)) {
      CRAB_LOG("region-load",
               CRAB_WARN("region_domain::ref_load: reference ", ref, " is null."
                         //, " Set to bottom ..."
               ););
      // set_to_bottom();
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
      CRAB_LOG("region-load", crab::outs() << "Reading from singleton\n";);
      if (auto region_var_opt = get_var(rgn)) {
        ref_load(*region_var_opt);
      } else {
        m_base_dom -= base_res;
      }
    } else {
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

  // Write the content of val to the address pointed by ref in region
  // rgn.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    if (is_bottom()) {
      return;
    }

    // TODO: handle unknown regions
    if (rgn.get_type().is_unknown_region()) {
      CRAB_LOG(
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      return;
    }

    bool is_uninitialized_rgn = m_rgn_init_dom[rgn].is_false();

    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    base_variable_t rgn_var = rename_var(rgn);

    // Non-scalar regions should be handled by ref_store_to_array
    if (!rgn.get_type().is_scalar_region()) {
      CRAB_LOG("region-store", CRAB_WARN("region_domain::ref_store: ", rgn,
                                         " must be a scalar region"););
      m_base_dom -= rgn_var;
      return;
    } else {
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
    }

    if (!ref.get_type().is_reference()) {
      // This shouldn't happen
      CRAB_LOG("region-store", CRAB_WARN("region_domain::ref_store: ", ref,
                                         " must be a reference"););
      m_base_dom -= rgn_var;
      return;
    }

    if (is_null_ref(ref)) {
      CRAB_LOG("region-store",
               CRAB_WARN("region_domain::ref_store: reference ", ref,
                         " is null. " //, " Set to bottom ..."
               ););
      // set_to_bottom();
      m_base_dom -= rgn_var;
      return;
    }

    auto ref_store = [&rgn_var, &val, this](base_abstract_domain_t &base_dom) {
      // At this point region of references has been translated to region of
      // integers. val can be any scalar including a reference.
      if (val.get_type().is_bool()) {
        if (val.is_constant()) {
          base_dom.assign_bool_cst(
              rgn_var,
              (val.is_bool_true() ? base_linear_constraint_t::get_true()
                                  : base_linear_constraint_t::get_false()));
        } else {
          base_dom.assign_bool_var(rgn_var, rename_var(val.get_variable()),
                                   false);
        }
      } else if (val.get_type().is_integer() || val.get_type().is_real() ||
                 val.get_type().is_reference()) {
        if (val.is_constant()) {
          base_dom.assign(rgn_var, val.get_constant());
        } else {
          base_dom.assign(rgn_var, rename_var(val.get_variable()));
        }
      } else {
        CRAB_ERROR(domain_name(), "::ref_store: unsupported type ",
                   val.get_type());
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
      ref_store(m_base_dom);
    } else {
      /* weak update */
      CRAB_LOG("region-store", crab::outs() << "Performing weak update\n";);
      base_abstract_domain_t tmp(m_base_dom);
      ref_store(tmp);
      m_base_dom |= tmp;
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
      if (Params::deallocation) {
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
        } else {
          if (found1) {
            // forget about the whole equivalence class to play safe
            m_rgn_dealloc_dom.remove_equiv_class(rgn1);
          }
          if (found2) {
            // forget about the whole equivalence class to play safe.
            m_rgn_dealloc_dom.remove_equiv_class(rgn2);
          }
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

    assert(ref_var.get_type().is_reference());
    assert(int_var.get_type().is_integer());

    if (!is_bottom()) {
      base_variable_t src_var = rename_var(ref_var);
      base_variable_t dst_var = rename_var(int_var);
      m_base_dom.assign(dst_var, src_var);
    }
  }

  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    crab::CrabStats::count(domain_name() + ".count.int_to_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".int_to_ref");

    assert(ref_var.get_type().is_reference());
    assert(int_var.get_type().is_integer());

    if (!is_bottom()) {
      base_variable_t src_var = rename_var(int_var);
      base_variable_t dst_var = rename_var(ref_var);

      m_base_dom.assign(dst_var, src_var);

      // Update region counting
      auto num_refs = m_rgn_counting_dom[rgn];
      m_rgn_counting_dom.set(rgn, num_refs.increment());
    }
  }

  // TODO: this default implementation is expensive because it will
  // call the join
  DEFAULT_SELECT_REF(region_domain_t)
  
  void ref_remove(const variable_t &ref, const variable_t &rgn) {
    if (!Params::deallocation) {
      return;
    }

    if (is_bottom()) {
      return;
    }
    assert(ref.get_type().is_reference());
    assert(rgn.get_type().is_region());

    // We conservatively mark the region's equivalence class as
    // possibly deallocated
    if (Params::deallocation) {
      if (!m_rgn_dealloc_dom.contains(rgn)) {
        CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn,
                                     " in ", m_rgn_dealloc_dom));
      } else {
        m_rgn_dealloc_dom.set(rgn, boolean_value::top());
      }
    }

    CRAB_LOG("region", crab::outs() << "After ref_remove(" << ref << "," << rgn
                                    << ")=" << *this << "\n";);
  }

  // arithmetic operations
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
    }
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), k);
    }
  }
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_base_dom.assign(rename_var(x), b_e);
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    if (!is_bottom()) {
      auto b_e1 = rename_linear_expr(e1);
      auto b_e2 = rename_linear_expr(e2);
      auto b_cond = rename_linear_cst(cond);
      m_base_dom.select(rename_var(lhs), b_cond, b_e1, b_e2);
    }
  }
  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    if (!is_bottom()) {
      auto b_e = rename_linear_expr(e);
      m_base_dom.backward_assign(rename_var(x), b_e, invariant.m_base_dom);
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
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), z);
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
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      auto b_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.assign_bool_cst(rename_var(lhs), b_rhs);
    }
  }

  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.assign_bool_var(rename_var(lhs), rename_var(rhs), is_not_rhs);
    }
  }
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.apply_binary_bool(op, rename_var(x), rename_var(y),
                                   rename_var(z));
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
    }
  }
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".array_assign");

    if (!is_bottom()) {
      m_base_dom.array_assign(rename_var(lhs), rename_var(rhs));
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
    }
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_assign");

    if (!is_bottom()) {
      m_base_dom.backward_array_assign(rename_var(lhs), rename_var(rhs),
                                       invariant.m_base_dom);
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    for (auto &v : variables) {
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;
        if (Params::deallocation) {
          m_rgn_dealloc_dom.forget(v);
        }
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
    if (Params::deallocation) {
      m_rgn_dealloc_dom.project(sorted_variables);
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

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    CRAB_ERROR("region_domain::rename not implemented");
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    //=================================================================//
    //       Special intrinsics supported by the region domain
    //=================================================================//
    // free(rgn, ref)
    //       deallocate the memory object pointed by the reference
    // b := is_unfreed_or_null(rgn, ref)
    //       b is true if the reference does not point to a freed
    //       memory object or it is null.
    // unfreed_or_null(rgn, ref)
    //       ensures that the reference does not point to a freed memory object
    // nonnull(ref)
    //       ensures that the reference is not null
    //=================================================================//
    auto error_if_not_arity = [&name, &inputs, &outputs](unsigned num_inputs,
                                                         unsigned num_outputs) {
      if (inputs.size() != num_inputs || outputs.size() != num_outputs) {
        CRAB_ERROR("Intrinsics ", name, " unexpected number of parameters");
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

    if (is_bottom()) {
      return;
    }

    if (Params::deallocation && name == "free") {
      error_if_not_arity(2, 0);
      variable_t rgn = inputs[0];
      variable_t ref = inputs[1];
      error_if_not_rgn(rgn);
      error_if_not_ref(ref);

      ref_remove(ref, rgn);
    } else if (Params::deallocation && name == "is_unfreed_or_null") {
      error_if_not_arity(2, 1);
      variable_t rgn = inputs[0];
      variable_t ref = inputs[1];
      variable_t bv = outputs[0];
      error_if_not_rgn(rgn);
      error_if_not_ref(ref);
      error_if_not_bool(bv);

      auto set_bool_var_to_true = [this](const variable_t &bool_var) {
        /// Require that the base domain can reason about booleans
        operator-=(bool_var);
        assume_bool(bool_var, false /*not negated*/);
      };

      // the reference is definitely null
      if (is_null_ref(ref)) {
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
    } else if (Params::deallocation && name == "unfreed_or_null") {
      error_if_not_arity(2, 0);
      variable_t rgn = inputs[0];
      variable_t ref = inputs[1];
      error_if_not_rgn(rgn);
      error_if_not_ref(ref);

      if (!is_null_ref(ref)) {
        if (Params::deallocation) {
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
      variable_t ref = inputs[0];
      error_if_not_ref(ref);

      auto nonnull_cst = reference_constraint_t::mk_gt_null(ref);
      ref_assume(nonnull_cst);
    } else {
      // pass the intrinsics to the base domain
      //=== base domain ===/
      std::vector<base_variable_t> base_inputs, base_outputs;
      base_inputs.reserve(inputs.size());
      base_outputs.reserve(outputs.size());
      for (unsigned i = 0, sz = inputs.size(); i < sz; ++i) {
        base_inputs.push_back(rename_var(inputs[i]));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
        base_outputs.push_back(rename_var(outputs[i]));
      }
      m_base_dom.intrinsic(name, base_inputs, base_outputs);
    }
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const region_domain_t &invariant) override {
    if (!is_bottom()) {
      //=== base domain ===/
      std::vector<base_variable_t> base_inputs, base_outputs;
      base_inputs.reserve(inputs.size());
      base_outputs.reserve(outputs.size());
      for (unsigned i = 0, sz = inputs.size(); i < sz; ++i) {
        base_inputs.push_back(rename_var(inputs[i]));
      }
      for (unsigned i = 0, sz = outputs.size(); i < sz; ++i) {
        base_outputs.push_back(rename_var(outputs[i]));
      }
      m_base_dom.backward_intrinsic(name, base_inputs, base_outputs,
                                    invariant.m_base_dom);
    }
  }
  /* end intrinsics operations */

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
      CRAB_LOG("region-print", o << "(RgnCounter=" << m_rgn_counting_dom << ","
                                 << "RgnInit=" << m_rgn_init_dom << ","
                                 << "RgnDealloc=" << m_rgn_dealloc_dom
                                 << ")\n";);
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
