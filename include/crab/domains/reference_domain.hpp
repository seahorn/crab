#pragma once

#include <crab/cfg/var_factory.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/killgen_domain.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/memory_regions.hpp>
#include <crab/types/reference_constraints.hpp>

#include <functional>
#include <unordered_map>

namespace crab {
namespace domains {

/*
 * Class that represents the range of intervals
 * {[0,0], [1,1], [0,1], [0,+oo], [1,+oo]}
 */
class small_range {

  typedef enum {
    Bottom, /*unused*/
    ExactlyZero,
    ExactlyOne,
    ZeroOrOne,
    ZeroOrMore,
    OneOrMore
  } kind_t;

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
      return small_range(OneOrMore);
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
    if (other.m_value == ExactlyZero) {
      return small_range(ZeroOrOne);
    } else if (other.m_value == ExactlyOne || other.m_value == OneOrMore) {
      return small_range(OneOrMore);
    } else {
      return small_range(ZeroOrMore);
    }
  }

public:
  small_range() : m_value(ZeroOrMore) {}

  static small_range bottom() { return small_range(Bottom); }
  static small_range top() { return small_range(ZeroOrMore); }
  static small_range zero() { return small_range(ExactlyZero); }
  static small_range one() { return small_range(ExactlyOne); }

  small_range(const small_range &other) = default;
  small_range(small_range &&other) = default;
  small_range &operator=(const small_range &other) = default;
  small_range &operator=(small_range &&other) = default;

  bool is_bottom() { return (m_value == Bottom); }

  bool is_top() { return (m_value == ZeroOrMore); }

  bool is_zero() const { return (m_value == ExactlyZero); }

  bool is_one() const { return (m_value == ExactlyOne); }

  bool operator<=(small_range other) {
    if (m_value == Bottom || other.m_value == ZeroOrMore) {
      return true;
    } else if (m_value == ExactlyZero) {
      return other.m_value != ExactlyOne && other.m_value != OneOrMore;
    } else if (m_value == ExactlyOne) {
      return other.m_value != ExactlyZero;
    } else if (m_value == ZeroOrOne || m_value == OneOrMore ||
               m_value == ZeroOrMore) {
      assert(other.m_value != ZeroOrMore);
      return false;
    }

    // should be unreachable
    return false;
  }

  bool operator==(small_range other) { return m_value == other.m_value; }

  small_range operator|(small_range other) {
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

  small_range operator||(small_range other) { return *this | other; }

  small_range operator&(small_range other) {
    // TODO
    return *this;
  }

  small_range operator&&(small_range other) { return *this & other; }

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

namespace reference_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::cfg::var_factory_impl::str_var_alloc_col;
  using base_abstract_domain_t = BaseAbsDom;
  using base_varname_t = typename BaseAbsDom::varname_t;

  static_assert(std::is_same<Number, typename BaseAbsDom::number_t>::value,
                "Number type and BaseAbsDom::number_t must be the same");
  // This is a strong requirement
  static_assert(
      std::is_same<base_varname_t,
                   typename varname_allocator_t::varname_t>::value,
      "BaseAbsDom::varname_t and allocator_varname_t must be the same");
};
} // end namespace reference_domain_impl

// Abstract domain for references
template <typename Params>
class reference_domain final
    : public abstract_domain<reference_domain<Params>> {
  using reference_domain_t = reference_domain<Params>;
  using abstract_domain_t = abstract_domain<reference_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::interval_t;
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_variable_t = typename base_abstract_domain_t::variable_t;
  using base_linear_expression_t =
      typename base_abstract_domain_t::linear_expression_t;
  using base_linear_constraint_t =
      typename base_abstract_domain_t::linear_constraint_t;
  using base_linear_constraint_system_t =
      typename base_abstract_domain_t::linear_constraint_system_t;
  using base_varname_allocator_t = typename Params::varname_allocator_t;

  // Abstract domains
  using ref_counting_domain_t = ikos::separate_domain<memory_region, small_range>;
  using regions_powerset_t = flat_killgen_domain<memory_region>;
  using regions_domain_t =
      separate_killgen_domain<variable_t, regions_powerset_t>;
  // Variable map: region to base domain's variable.
  using region_to_content_map_t =
      std::unordered_map<memory_region, base_variable_t>;
  using rev_region_to_content_map_t =
      std::unordered_map<base_variable_t, memory_region>;
  // Variable map: variable to base domain's variable.
  using var_map_t = std::unordered_map<variable_t, base_variable_t>;
  using rev_var_map_t = std::unordered_map<base_variable_t, variable_t>;

  bool m_is_bottom;

  // To create synthetic variables for the base domain.
  base_varname_allocator_t m_alloc;

  // Map region to (scalar or array) variable to represent its whole
  // content. This allows us to model region contents
  region_to_content_map_t m_region_to_content;
  rev_region_to_content_map_t m_rev_region_to_content;
  // Map variable_to to base domain's variable:
  //  - int/real/bool variable is mapped while preserving its type
  //  - ref variables are mapped to **integer** variables to ensure
  //    that the base domain is completely unaware of references.
  var_map_t m_var_map;
  rev_var_map_t m_rev_var_map;
  // Abstract domain to count unique references per region.  This
  // allows us to decide when strong update is sound: only if exactly
  // one reference to the region.
  ref_counting_domain_t m_ref_counting_dom;
  // Abstract domain to keep track all possible regions for a
  // reference variable.
  //
  // I think a reference can only point to one region but we check
  // that ref_load or ref_store use references and regions which are
  // matched together via ref_make.
  regions_domain_t m_regions_dom;
  // The base abstract domain: all the heavy lifting is done here.
  // m_base_dom does not have any variable of reference type.
  base_abstract_domain_t m_base_dom;

  reference_domain(base_varname_allocator_t &&alloc,
                   region_to_content_map_t &&region_to_content,
                   rev_region_to_content_map_t &&rev_region_to_content,
                   var_map_t &&var_map, rev_var_map_t &&rev_var_map,
                   ref_counting_domain_t &&ref_counting_dom,
                   regions_domain_t &&regions_dom,
                   base_abstract_domain_t &&base_dom)
      : m_is_bottom(base_dom.is_bottom()), m_alloc(std::move(alloc)),
        m_region_to_content(std::move(region_to_content)),
        m_rev_region_to_content(std::move(rev_region_to_content)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_ref_counting_dom(std::move(ref_counting_dom)),
        m_regions_dom(std::move(regions_dom)), m_base_dom(std::move(base_dom)) {
  }

  using ref_counting_binop_t = std::function<ref_counting_domain_t(
      ref_counting_domain_t, ref_counting_domain_t)>;
  using regions_dom_binop_t =
      std::function<regions_domain_t(regions_domain_t, regions_domain_t)>;
  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;

  reference_domain_t do_join_or_widening(reference_domain_t &left,
                                         reference_domain_t &right,
                                         ref_counting_binop_t ref_counting_op,
                                         regions_dom_binop_t regions_dom_op,
                                         base_dom_binop_t base_dom_op) {
    ref_counting_domain_t out_ref_counting_dom(
        ref_counting_op(left.m_ref_counting_dom, right.m_ref_counting_dom));

    regions_domain_t out_regions_dom(
        regions_dom_op(left.m_regions_dom, right.m_regions_dom));

    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);
    base_abstract_domain_t left_dom(left.m_base_dom);
    base_abstract_domain_t right_dom(right.m_base_dom);
    region_to_content_map_t out_region_to_content;
    rev_region_to_content_map_t out_rev_region_to_content;
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;

    // -- Compute common renamings
    variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings =
        left.m_region_to_content.size() + left.m_var_map.size();
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

    for (auto &kv : left.m_region_to_content) {
      const memory_region &mem = kv.first;
      auto it = right.m_region_to_content.find(mem);
      if (it != right.m_region_to_content.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, mem)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_region_to_content.insert({mem, out_v});
        out_rev_region_to_content.insert({out_v, mem});
      }
    }

    left_dom.rename(left_vars, out_vars);
    right_dom.rename(right_vars, out_vars);
    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    reference_domain_t res(
        std::move(out_alloc), std::move(out_region_to_content),
        std::move(out_rev_region_to_content), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_ref_counting_dom),
        std::move(out_regions_dom), std::move(out_base_dom));
    return res;
  }

  reference_domain_t do_meet_or_narrowing(reference_domain_t &left,
                                          reference_domain_t &right,
                                          ref_counting_binop_t ref_counting_op,
                                          regions_dom_binop_t regions_dom_op,
                                          base_dom_binop_t base_dom_op) {

    ref_counting_domain_t out_ref_counting_dom(
        ref_counting_op(left.m_ref_counting_dom, right.m_ref_counting_dom));

    regions_domain_t out_regions_dom(
        regions_dom_op(left.m_regions_dom, right.m_regions_dom));

    base_varname_allocator_t out_alloc(left.m_alloc, right.m_alloc);
    base_abstract_domain_t left_dom(left.m_base_dom);
    base_abstract_domain_t right_dom(right.m_base_dom);
    region_to_content_map_t out_region_to_content;
    rev_region_to_content_map_t out_rev_region_to_content;
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;

    variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings =
        left.m_region_to_content.size() + left.m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);

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
        // only on left: no renaming needed
        out_var_map.insert(kv);
        out_rev_var_map.insert({kv.second, kv.first});
      }
    }

    // add missing maps from right operand
    for (auto &kv : right.m_var_map) {
      const variable_t &v = kv.first;
      auto it = left.m_var_map.find(v);
      if (it == left.m_var_map.end()) {
        // only right left: no renaming needed
        out_var_map.insert(kv);
        out_rev_var_map.insert({kv.second, kv.first});
      }
    }

    for (auto &kv : left.m_region_to_content) {
      const memory_region &mem = kv.first;
      auto it = right.m_region_to_content.find(mem);
      if (it != right.m_region_to_content.end()) {
        // common renaming
        base_variable_t out_v(std::move(make_base_variable(out_alloc, mem)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_region_to_content.insert({mem, out_v});
        out_rev_region_to_content.insert({out_v, mem});
      } else {
        // only on left: no renaming needed
        out_region_to_content.insert(kv);
        out_rev_region_to_content.insert({kv.second, kv.first});
      }
    }

    // add missing maps from right operand
    for (auto &kv : right.m_region_to_content) {
      const memory_region &mem = kv.first;
      auto it = left.m_region_to_content.find(mem);
      if (it == left.m_region_to_content.end()) {
        // only right left: no renaming needed
        out_region_to_content.insert(kv);
        out_rev_region_to_content.insert({kv.second, kv.first});
      }
    }

    left_dom.rename(left_vars, out_vars);
    right_dom.rename(right_vars, out_vars);

    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    reference_domain_t res(
        std::move(out_alloc), std::move(out_region_to_content),
        std::move(out_rev_region_to_content), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_ref_counting_dom),
        std::move(out_regions_dom), std::move(out_base_dom));
    return res;
  }

  // Create a fresh variable in the base domain to shadow reg
  base_variable_t make_base_variable(base_varname_allocator_t &var_allocator,
                                     const memory_region &reg) {
    typename base_variable_t::type_t ty = UNK_TYPE;
    unsigned bitwidth = 0;
    switch (reg.get_type()) {
    case memory_region::type_t::BOOL:
      ty = BOOL_TYPE;
      break;
    case memory_region::type_t::INT:
      ty = INT_TYPE;
      bitwidth = reg.get_bitwidth();
      break;
    case memory_region::type_t::REAL:
      ty = REAL_TYPE;
      break;
    case memory_region::type_t::REF:
      ty = INT_TYPE;
      bitwidth = 32;
      break;
    default:;
      CRAB_ERROR("reference_domain::make_base_variable: unreachable");
    }
    base_variable_t v(var_allocator.next(), ty, bitwidth);
    return v;
  }

  // Create a fresh variable in the base domain to shadow v
  base_variable_t make_base_variable(base_varname_allocator_t &var_allocator,
                                     const variable_t &v) {
    if (v.is_ref_type()) {
      base_variable_t v(var_allocator.next(), crab::INT_TYPE, 32);
      return v;
    } else {
      base_variable_t nv(var_allocator.next(), v.get_type(), v.get_bitwidth());
      return nv;
    }
  }

  // Get the synthetic variable in the base domain that represents
  // reg's content.
  boost::optional<base_variable_t>
  get_region_var(const memory_region &reg) const {
    auto it = m_region_to_content.find(reg);
    if (it != m_region_to_content.end()) {
      return it->second;
    } else {
      return boost::optional<base_variable_t>();
    }
  }

  bool has_region_var(const memory_region &reg) const {
    return (m_region_to_content.find(reg) != m_region_to_content.end());
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
    auto it = m_var_map.find(v);
    if (it != m_var_map.end()) {
      return it->second;
    } else {
      base_variable_t bv = make_base_variable(m_alloc, v);
      auto res = m_var_map.insert({v, bv});
      if (res.second) {
        m_rev_var_map.insert({bv, v});
      }
      return res.first->second;
    }
  }

  // used only for pretty-printing and to_linear_constraint_system()
  boost::optional<variable_t> rev_rename_var(const base_variable_t &v,
                                             bool ignore_references) {
    auto it = m_rev_var_map.find(v);
    if (it != m_rev_var_map.end()) {
      if (!ignore_references || !it->second.is_ref_type()) {
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
                         bool ignore_references) {
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
    return base_linear_constraint_t(
        rename_linear_expr(cst.expression()),
        (typename base_linear_constraint_t::kind_t)cst.kind());
  }

  // used only for pretty-printing and to_linear_constraint_system()
  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const base_linear_constraint_t &cst,
                        bool ignore_references) {
    if (boost::optional<linear_expression_t> e =
            rev_rename_linear_expr(cst.expression(), ignore_references)) {
      return linear_constraint_t(
          *e, (typename linear_constraint_t::kind_t)cst.kind());
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

  // Return true if ref is definitely null. Ask the base numerical
  // domain for that.
  bool is_null_ref(const variable_t &ref) {
    if (!ref.is_ref_type()) {
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
  void set_to_top() override {
    reference_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    reference_domain_t abs(false);
    std::swap(*this, abs);
  }

  reference_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  reference_domain(const reference_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_alloc(o.m_alloc),
        m_region_to_content(o.m_region_to_content),
        m_rev_region_to_content(o.m_rev_region_to_content),
        m_var_map(o.m_var_map), m_rev_var_map(o.m_rev_var_map),
        m_ref_counting_dom(o.m_ref_counting_dom),
        m_regions_dom(o.m_regions_dom), m_base_dom(o.m_base_dom) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }
  reference_domain(reference_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_region_to_content(std::move(o.m_region_to_content)),
        m_rev_region_to_content(std::move(o.m_rev_region_to_content)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_ref_counting_dom(std::move(o.m_ref_counting_dom)),
        m_regions_dom(std::move(o.m_regions_dom)),
        m_base_dom(std::move(o.m_base_dom)) {}

  reference_domain_t &operator=(const reference_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_region_to_content = o.m_region_to_content;
      m_rev_region_to_content = o.m_rev_region_to_content;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_ref_counting_dom = o.m_ref_counting_dom;
      m_regions_dom = o.m_regions_dom;
      m_base_dom = o.m_base_dom;
    }
    return *this;
  }

  reference_domain_t &operator=(reference_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_alloc = std::move(o.m_alloc);
      m_region_to_content = std::move(o.m_region_to_content);
      m_rev_region_to_content = std::move(o.m_rev_region_to_content);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_ref_counting_dom = std::move(o.m_ref_counting_dom);
      m_regions_dom = std::move(o.m_regions_dom);
      m_base_dom = std::move(o.m_base_dom);
    }
    return *this;
  }

  bool is_bottom() override { return m_is_bottom; }

  bool is_top() override { return m_base_dom.is_top(); }

  bool operator<=(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    }

    if (!(m_ref_counting_dom <= o.m_ref_counting_dom)) {
      return false;
    }

    if (!(m_regions_dom <= o.m_regions_dom)) {
      return false;
    }

    // perform the common renaming for "region" variables
    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);
    base_abstract_domain_t left_dom(m_base_dom);
    base_abstract_domain_t right_dom(o.m_base_dom);
    variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings = m_region_to_content.size() + m_var_map.size();
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

    for (auto &kv : m_region_to_content) {
      const memory_region &mem = kv.first;
      auto it = o.m_region_to_content.find(mem);
      if (it != o.m_region_to_content.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, mem)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
      }
    }

    left_dom.rename(left_vars, out_vars);
    right_dom.rename(right_vars, out_vars);
    return left_dom <= right_dom;
  }

  void operator|=(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

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

    CRAB_LOG("reference", crab::outs()
                              << "Join " << *this << " and " << o << "=");

    ref_counting_domain_t out_ref_counting_dom(m_ref_counting_dom |
                                               o.m_ref_counting_dom);
    regions_domain_t out_regions_dom(m_regions_dom | o.m_regions_dom);
    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);
    base_abstract_domain_t left_dom(m_base_dom);
    base_abstract_domain_t right_dom(o.m_base_dom);
    region_to_content_map_t out_region_to_content;
    rev_region_to_content_map_t out_rev_region_to_content;
    var_map_t out_var_map;
    rev_var_map_t out_rev_var_map;
    variable_vector_t left_vars, right_vars, out_vars;
    // upper bound to avoid reallocations
    size_t num_renamings = m_region_to_content.size() + m_var_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    // perform the common renamings

    for (auto &kv : m_var_map) {
      const variable_t &v = kv.first;
      auto it = o.m_var_map.find(v);
      if (it != o.m_var_map.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, v)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_var_map.insert({v, out_v});
        out_rev_var_map.insert({out_v, v});
      }
    }

    for (auto &kv : m_region_to_content) {
      const memory_region &mem = kv.first;
      auto it = o.m_region_to_content.find(mem);
      if (it != o.m_region_to_content.end()) {
        base_variable_t out_v(std::move(make_base_variable(out_alloc, mem)));
        left_vars.push_back(kv.second);
        right_vars.push_back(it->second);
        out_vars.push_back(out_v);
        out_region_to_content.insert({mem, out_v});
        out_rev_region_to_content.insert({out_v, mem});
      }
    }

    left_dom.rename(left_vars, out_vars);
    right_dom.rename(right_vars, out_vars);
    base_abstract_domain_t out_base_dom(left_dom | right_dom);

    m_is_bottom = out_base_dom.is_bottom();
    std::swap(m_alloc, out_alloc);
    std::swap(m_ref_counting_dom, out_ref_counting_dom);
    std::swap(m_region_to_content, out_region_to_content);
    std::swap(m_rev_region_to_content, out_rev_region_to_content);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_regions_dom, out_regions_dom);
    std::swap(m_base_dom, out_base_dom);

    CRAB_LOG("reference", crab::outs() << *this << "\n");
  }

  reference_domain_t operator|(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    // Trivial cases first
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() | o.is_top()) {
      reference_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Join " << *this << " and " << o << "=");

    auto ref_counting_op = [](ref_counting_domain_t v1,
                              ref_counting_domain_t v2) { return v1 | v2; };
    auto region_op = [](regions_domain_t v1, regions_domain_t v2) {
      return v1 | v2;
    };
    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return v1 | v2; };
    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, ref_counting_op, region_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator&(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Meet " << *this << " and " << o << "=");

    auto ref_counting_op = [](ref_counting_domain_t v1,
                              ref_counting_domain_t v2) { return (v1 & v2); };
    auto region_op = [](regions_domain_t v1, regions_domain_t v2) {
      return (v1 & v2);
    };
    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return (v1 & v2); };
    reference_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, ref_counting_op, region_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator||(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Widening " << *this << " and " << o << "=");

    auto ref_counting_op = [](ref_counting_domain_t v1,
                              ref_counting_domain_t v2) { return v1 || v2; };
    auto region_op = [](regions_domain_t v1, regions_domain_t v2) {
      return v1 | v2;
    };
    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return v1 || v2; };
    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, ref_counting_op, region_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t widening_thresholds(
      reference_domain_t o,
      const iterators::thresholds<number_t> &thresholds) override {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Widening " << *this << " and " << o << "=");

    auto ref_counting_op = [](ref_counting_domain_t v1,
                              ref_counting_domain_t v2) { return v1 || v2; };
    auto region_op = [](regions_domain_t v1, regions_domain_t v2) {
      return v1 | v2;
    };
    auto base_dom_op = [&thresholds](base_abstract_domain_t v1,
                                     base_abstract_domain_t v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, ref_counting_op, region_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator&&(reference_domain_t o) override {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Meet " << *this << " and " << o << "=");

    auto ref_counting_op = [](ref_counting_domain_t v1,
                              ref_counting_domain_t v2) { return (v1 && v2); };
    auto region_op = [](regions_domain_t v1, regions_domain_t v2) {
      return (v1 & v2);
    };
    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return (v1 && v2); };
    reference_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, ref_counting_op, region_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  void operator-=(variable_t v) override {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (!is_bottom()) {
      if (v.is_ref_type()) {
        m_regions_dom -= v;
      }
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        m_base_dom -= it->second;
        m_rev_var_map.erase(it->second);
        m_var_map.erase(it);
      }
    }
  }

  // Initialize region reg.
  void region_init(memory_region reg) override {
    crab::CrabStats::count(getDomainName() + ".count.region_init");
    crab::ScopedCrabStats __st__(getDomainName() + ".region_init");

    if (is_bottom()) {
      return;
    }

    if (has_region_var(reg)) {
      CRAB_ERROR("reference_domain::init_region: ", reg,
                 " cannot be initialized twice");
    }

    // Set to zero the number of references
    m_ref_counting_dom.set(reg, small_range::zero());

    // Assign a synthetic variable to reg for modeling its content.
    base_variable_t v = make_base_variable(m_alloc, reg);
    m_region_to_content.insert({reg, v});
    m_rev_region_to_content.insert({v, reg});

    CRAB_LOG("reference", crab::outs() << "After region_init(" << reg
                                       << ")=" << *this << "\n";);
  }

  // Create a new reference ref to region reg.
  void ref_make(variable_t ref, memory_region reg) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_make");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_make");

    if (is_bottom()) {
      return;
    }

    if (!has_region_var(reg)) {
      CRAB_ERROR("reference_domain::ref_make: ", reg, " must be initialized");
    }
    if (!ref.is_ref_type()) {
      CRAB_ERROR("reference_domain::ref_make: ", ref, " must be a reference");
    }

    // Update reference counting
    auto num_refs = m_ref_counting_dom[reg];
    m_ref_counting_dom.set(reg, num_refs.increment());

    // Assign a base domain variable to ref
    rename_var(ref);

    // Remember that ref points to reg
    auto region_set = m_regions_dom[ref];
    m_regions_dom.set(ref, region_set | reg);
    CRAB_LOG("reference", crab::outs() << "After ref_make(" << ref << "," << reg
                                       << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within reg. The content is
  // stored in res.
  void ref_load(variable_t ref, memory_region reg, variable_t res) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_load");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_load");

    if (is_bottom()) {
      return;
    }

    if (!has_region_var(reg)) {
      CRAB_ERROR("reference_domain::ref_load: ", reg, " must be initialized");
    }
    if (!ref.is_ref_type()) {
      CRAB_ERROR("reference_domain::ref_load: ", ref, " must be a reference");
    }

    if (is_null_ref(ref)) {
      CRAB_WARN("reference_domain::ref_load: rerefence ", ref, " is null. ",
                " Set to bottom ...");
      set_to_bottom();
      return;
    }

    const base_variable_t &base_res = rename_var(res);
    auto num_refs = m_ref_counting_dom[reg];
    if (num_refs.is_zero()) {
      CRAB_WARN("reference_domain::ref_load: TODO region has no references");
      m_base_dom -= base_res;
    } else {
      auto ref_load = [&ref, &reg, &base_res,
                       this](const base_variable_t &reg_var) {
        auto regions_set = this->m_regions_dom[ref];
        if (!(regions_powerset_t(reg) == regions_set)) {
          CRAB_WARN("reference_domain::ref_load: reference ", ref,
                    " is not a reference of region ", reg);
          this->m_base_dom -= base_res;
          return;
        }

        if (reg_var.get_type() != base_res.get_type()) {
          CRAB_ERROR("reference_domain::ref_load: ", "Type of region ", reg,
                     " does not match with ", base_res);
        }

        switch (base_res.get_type()) {
        case BOOL_TYPE:
          this->m_base_dom.assign_bool_var(base_res, reg_var, false);
          break;
        case INT_TYPE:
        case REAL_TYPE:
          this->m_base_dom.assign(base_res, reg_var);
          break;
        default:
          // variables of type base_variable_t cannot be REF_TYPE
          CRAB_ERROR("reference_domain::ref_load: unsupported type in ",
                     base_res);
        }
      };

      if (num_refs.is_one()) {
        CRAB_LOG("reference", crab::outs() << "Reading from singleton\n";);
        base_variable_t region_var = *(get_region_var(reg));
        ref_load(region_var);
      } else {
        CRAB_LOG("reference", crab::outs() << "Reading from non-singleton\n";);
        base_variable_t region_var = *(get_region_var(reg));
        base_variable_t fresh_region_var = make_base_variable(m_alloc, reg);
        m_base_dom.expand(region_var, fresh_region_var);
        ref_load(fresh_region_var);
      }
    }
    CRAB_LOG("reference", crab::outs() << "After " << res << ":="
                                       << "ref_load(" << ref << "," << reg
                                       << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in region
  // reg.
  void ref_store(variable_t ref, memory_region reg,
                 linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_store");

    if (is_bottom()) {
      return;
    }

    if (!has_region_var(reg)) {
      CRAB_ERROR("reference_domain::ref_store: ", reg, " must be initialized");
    }
    if (!ref.is_ref_type()) {
      CRAB_ERROR("reference_domain::ref_store: ", ref, " must be a reference");
    }

    if (is_null_ref(ref)) {
      CRAB_LOG("reference",
               CRAB_WARN("reference_domain::ref_store: rerefence ", ref,
                         " is null. ", " Set to bottom ..."););
      set_to_bottom();
      return;
    }

    // const base_variable_t &base_val = rename_var(val);
    auto num_refs = m_ref_counting_dom[reg];
    if (num_refs.is_zero()) {
      CRAB_WARN("reference_domain::ref_store: TODO region has no references");
    } else {
      auto ref_store = [&ref, &reg, &val,
                        this](base_abstract_domain_t &base_dom) {
        auto regions_set = this->m_regions_dom[ref];
        if (!(regions_powerset_t(reg) == regions_set)) {
          CRAB_WARN("reference_domain::ref_store: reference ", ref,
                    " is not a reference of region ", reg);
          return;
        }

        base_variable_t reg_var = *(get_region_var(reg));
        switch (reg_var.get_type()) {
        case BOOL_TYPE:
          if (val.is_constant()) {
            if (val.constant() >= number_t(1)) {
              base_dom.assign_bool_cst(reg_var,
                                       base_linear_constraint_t::get_true());
            } else {
              base_dom.assign_bool_cst(reg_var,
                                       base_linear_constraint_t::get_false());
            }
          } else if (auto var = val.get_variable()) {
            base_dom.assign_bool_var(reg_var, rename_var(*var), false);
          } else {
            CRAB_ERROR("reference_domain::ref_store: unsupported boolean ",
                       val);
          }
          break;
        case INT_TYPE:
        case REAL_TYPE:
          base_dom.assign(reg_var, rename_linear_expr(val));
          break;
        default:
          // variables of type base_variable_t cannot be REF_TYPE
          CRAB_ERROR("reference_domain::ref_store: unsupported type in ",
                     reg_var);
        }
      };

      if (num_refs.is_one()) {
        /* strong update */
        CRAB_LOG("reference", crab::outs() << "Performing strong update\n";);
        ref_store(m_base_dom);
      } else {
        /* weak update */
        CRAB_LOG("reference", crab::outs() << "Performing weak update\n";);
        base_abstract_domain_t tmp(m_base_dom);
        ref_store(tmp);
        m_base_dom |= tmp;
      }
    }
    CRAB_LOG("reference", crab::outs()
                              << "After ref_store(" << ref << "," << reg << ","
                              << val << ")=" << *this << "\n";);
  }

  // Create a new reference ref2 to region reg2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(variable_t ref1, memory_region reg1, variable_t ref2,
               memory_region reg2, linear_expression_t offset) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_gep");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_gep");

    if (!is_bottom()) {
      ref_make(ref2, reg2);
    }
    CRAB_LOG("reference", crab::outs()
                              << "After ref_gep(" << ref1 << "," << reg1 << ","
                              << ref2 << "," << reg2 << "," << offset
                              << ")=" << *this << "\n";);
  }

  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(variable_t lhs, variable_t ref, memory_region region,
                           linear_expression_t index,
                           linear_expression_t elem_size) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_load_from_array");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_load_from_array");

    CRAB_WARN("reference_domain::ref_load_from_array not implemented");
  }

  // Treat region as an array and perform an array store.
  void ref_store_to_array(variable_t ref, memory_region region,
                          linear_expression_t index,
                          linear_expression_t elem_size,
                          linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_store_to_array");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_store_to_array");

    CRAB_WARN("reference_domain::ref_store_array not implemented");
  }

  void ref_assume(reference_constraint_t cst) override {
    crab::CrabStats::count(getDomainName() + ".count.ref_assume");
    crab::ScopedCrabStats __st__(getDomainName() + ".ref_assume");

    if (!is_bottom()) {
      if (cst.is_tautology()) {
        return;
      }
      if (cst.is_contradiction()) {
        set_to_bottom();
        return;
      }
      if (cst.is_unary()) {
        assert(cst.lhs().is_ref_type());
        if (auto ref_opt = get_var(cst.lhs())) {
          if (cst.is_equality()) {
            m_base_dom += linear_constraint_t(*ref_opt == number_t(0));
          } else if (cst.is_disequality()) {
            // m_base_dom += linear_constraint_t(*ref_opt != number_t(0));
            m_base_dom += linear_constraint_t(*ref_opt > number_t(0));
          } else if (cst.is_less_or_equal_than()) {
            m_base_dom += linear_constraint_t(*ref_opt <= number_t(0));
          } else if (cst.is_less_than()) {
            m_base_dom += linear_constraint_t(*ref_opt < number_t(0));
          }
        } else {
          CRAB_LOG("reference",
                   CRAB_WARN("reference_domain::ref_assume: no found ",
                             cst.lhs(), " as reference"));
        }
      } else {
        assert(cst.lhs().is_ref_type());
        assert(cst.rhs().is_ref_type());
        auto ref1_opt = get_var(cst.lhs());
        auto ref2_opt = get_var(cst.rhs());
        number_t offset = cst.offset();
        if (ref1_opt && ref2_opt) {
          if (cst.is_equality()) {
            m_base_dom += linear_constraint_t(*ref1_opt == *ref2_opt + offset);
          } else if (cst.is_disequality()) {
            m_base_dom += linear_constraint_t(*ref1_opt != *ref2_opt + offset);
          } else if (cst.is_less_or_equal_than()) {
            m_base_dom += linear_constraint_t(*ref1_opt <= *ref2_opt + offset);
          } else if (cst.is_less_than()) {
            m_base_dom += linear_constraint_t(*ref1_opt < *ref2_opt + offset);
          }
        } else {
          if (!ref1_opt) {
            CRAB_LOG("reference",
                     CRAB_WARN("reference_domain::ref_assume: no found ",
                               cst.lhs(), " as reference"));
          }
          if (!ref2_opt) {
            CRAB_LOG("reference",
                     CRAB_WARN("reference_domain::ref_assume: no found ",
                               cst.rhs(), " as reference"));
          }
        }
      }
    }
    m_is_bottom = m_base_dom.is_bottom();
    CRAB_LOG("reference", crab::outs()
                              << "ref_assume(" << cst << ")" << *this << "\n";);
  }

  // arithmetic operations
  void apply(arith_operation_t op, variable_t x, variable_t y,
             variable_t z) override {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
    }
  }
  void apply(arith_operation_t op, variable_t x, variable_t y, number_t k) override {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), k);
    }
  }
  void assign(variable_t x, linear_expression_t e) override {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (!is_bottom()) {
      m_base_dom.assign(rename_var(x), rename_linear_expr(e));
    }
  }
  void backward_assign(variable_t x, linear_expression_t e,
                       reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign");

    if (!is_bottom()) {
      m_base_dom.backward_assign(rename_var(x), rename_linear_expr(e),
                                 invariant.m_base_dom);
    }
  }
  void backward_apply(arith_operation_t op, variable_t x, variable_t y, number_t z,
                      reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, rename_var(x), rename_var(y), z,
                                invariant.m_base_dom);
    }
  }
  void backward_apply(arith_operation_t op, variable_t x, variable_t y, variable_t z,
                      reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, rename_var(x), rename_var(y), rename_var(z),
                                invariant.m_base_dom);
    }
  }
  void operator+=(linear_constraint_system_t csts) override {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

    if (!is_bottom()) {
      m_base_dom += rename_linear_cst_sys(csts);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  // not part of the numerical_domains api but it should be
  void set(variable_t x, interval_t intv) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (!is_bottom()) {
      if (intv.is_bottom()) {
        set_to_bottom();
      } else {
        m_base_dom.set(rename_var(x), intv);
      }
    }
  }
  interval_t operator[](variable_t x) {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom[rename_var(x)];
    }
  }

  // int cast operations and bitwise operations
  void apply(int_conv_operation_t op, variable_t dst, variable_t src) override {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(dst), rename_var(src));
    }
  }
  void apply(bitwise_operation_t op, variable_t x, variable_t y,
             variable_t z) override {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), rename_var(z));
    }
  }
  void apply(bitwise_operation_t op, variable_t x, variable_t y,
             number_t z) override {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, rename_var(x), rename_var(y), z);
    }
  }

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) override {
    crab::CrabStats::count(getDomainName() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.assign_bool_cst(rename_var(lhs), rename_linear_cst(rhs));
    }
  }
  void assign_bool_var(variable_t lhs, variable_t rhs,
                       bool is_not_rhs) override {
    crab::CrabStats::count(getDomainName() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.assign_bool_var(rename_var(lhs), rename_var(rhs), is_not_rhs);
    }
  }
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) override {
    crab::CrabStats::count(getDomainName() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.apply_binary_bool(op, rename_var(x), rename_var(y),
                                   rename_var(z));
    }
  }
  void assume_bool(variable_t v, bool is_negated) override {
    crab::CrabStats::count(getDomainName() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(getDomainName() + ".assume_bool");

    if (!is_bottom()) {
      m_base_dom.assume_bool(rename_var(v), is_negated);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_cst(
          rename_var(lhs), rename_linear_cst(rhs), invariant.m_base_dom);
    }
  }
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_var(rename_var(lhs), rename_var(rhs),
                                          is_not_rhs, invariant.m_base_dom);
    }
  }
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() +
                           ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".backward_apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.backward_apply_binary_bool(op, rename_var(x), rename_var(y),
                                            rename_var(z),
                                            invariant.m_base_dom);
    }
  }

  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.array_init");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_init");

    if (!is_bottom()) {
      m_base_dom.array_init(rename_var(a), rename_linear_expr(elem_size),
                            rename_linear_expr(lb_idx),
                            rename_linear_expr(ub_idx),
                            rename_linear_expr(val));
    }
  }
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) override {
    crab::CrabStats::count(getDomainName() + ".count.array_load");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_load");

    if (!is_bottom()) {
      m_base_dom.array_load(rename_var(lhs), rename_var(a),
                            rename_linear_expr(elem_size),
                            rename_linear_expr(i));
    }
  }
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store");

    if (!is_bottom()) {
      m_base_dom.array_store(rename_var(a), rename_linear_expr(elem_size),
                             rename_linear_expr(i), rename_linear_expr(v),
                             is_strong_update);
    }
  }
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store");

    if (!is_bottom()) {
      m_base_dom.array_store(
          rename_var(a_new), rename_var(a_old), rename_linear_expr(elem_size),
          rename_linear_expr(i), rename_linear_expr(v), is_strong_update);
    }
  }
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store_range");

    if (!is_bottom()) {
      m_base_dom.array_store_range(rename_var(a), rename_linear_expr(elem_size),
                                   rename_linear_expr(i), rename_linear_expr(j),
                                   rename_linear_expr(v));
    }
  }
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j,
                         linear_expression_t v) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store_range");

    if (!is_bottom()) {
      m_base_dom.array_store_range(rename_var(a_new), rename_var(a_old),
                                   rename_linear_expr(elem_size), i, j, v);
    }
  }
  void array_assign(variable_t lhs, variable_t rhs) override {
    crab::CrabStats::count(getDomainName() + ".count.array_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_assign");

    if (!is_bottom()) {
      m_base_dom.array_assign(rename_var(lhs), rename_var(rhs));
    }
  }
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_init");

    if (!is_bottom()) {
      m_base_dom.backward_array_init(
          rename_var(a), rename_linear_expr(elem_size),
          rename_linear_expr(lb_idx), rename_linear_expr(ub_idx),
          rename_linear_expr(val), invariant.m_base_dom);
    }
  }
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_load");

    if (!is_bottom()) {
      m_base_dom.backward_array_load(
          rename_var(lhs), rename_var(a), rename_linear_expr(elem_size),
          rename_linear_expr(i), invariant.m_base_dom);
    }
  }
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_store");

    if (!is_bottom()) {
      m_base_dom.backward_array_store(
          rename_var(a), rename_linear_expr(elem_size), rename_linear_expr(i),
          rename_linear_expr(v), is_strong_update, invariant.m_base_dom);
    }
  }
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_store");

    if (!is_bottom()) {
      m_base_dom.backward_array_store(
          rename_var(a_new), rename_var(a_old), rename_linear_expr(elem_size),
          rename_linear_expr(i), rename_linear_expr(v), is_strong_update,
          invariant.m_base_dom);
    }
  }
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() +
                           ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".backward_array_store_range");

    if (!is_bottom()) {
      m_base_dom.backward_array_store_range(
          rename_var(a), rename_linear_expr(elem_size), rename_linear_expr(i),
          rename_linear_expr(j), rename_linear_expr(v), invariant.m_base_dom);
    }
  }
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() +
                           ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".backward_array_store_range");

    if (!is_bottom()) {
      m_base_dom.backward_array_store_range(
          rename_var(a_new), rename_var(a_old), rename_linear_expr(elem_size),
          rename_linear_expr(i), rename_linear_expr(j), rename_linear_expr(v),
          invariant.m_base_dom);
    }
  }
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             reference_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_assign");

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
    std::transform(variables.begin(), variables.end(),
                   std::back_inserter(base_vars),
                   [this](const variable_t &v) { return this->rename_var(v); });
    m_base_dom.forget(base_vars);
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    std::vector<base_variable_t> base_vars;
    base_vars.reserve(variables.size());
    regions_domain_t out_regions_dom;
    for (auto const &v : variables) {
      if (v.is_ref_type()) {
        out_regions_dom.set(v, m_regions_dom[v]);
      }
      base_vars.push_back(rename_var(v));
    }
    m_base_dom.project(base_vars);
    std::swap(m_regions_dom, out_regions_dom);
  }

  void normalize() override {}
  void minimize() override {}

  void expand(variable_t x, variable_t new_x) override {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_ERROR("reference_domain::expand not implemented");
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    CRAB_ERROR("reference_domain::rename not implemented");
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          reference_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }
  /* end intrinsics operations */

  // Only integer or boolean variables are exposed.
  //
  // Variables that shadow memory regions and reference variables are
  // ignored.
  linear_constraint_system_t to_linear_constraint_system() override {
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
  to_disjunctive_linear_constraint_system() override {
    CRAB_ERROR("reference_domain::to_disjunctive_linear_constraint_system not "
               "implemented");
  }

  static std::string getDomainName() {
    return "ReferenceDomain(" + base_abstract_domain_t::getDomainName() + ")";
  }

  void write(crab_os &o) override {
    if (is_bottom()) {
      o << "_|_";
    } else {
      o << "(RefCounter=" << m_ref_counting_dom << ",";
      o << "RefToReg=" << m_regions_dom << ",";
#if 1
      o << "MapVar={";
      for (auto &kv : m_var_map) {
        o << kv.first << "->" << kv.second << ";";
      }
      o << "},";
      o << "MapRegToVar={";
      for (auto &kv : m_region_to_content) {
        o << kv.first << "->" << kv.second << ";";
      }
      o << "},";
#endif
      o << "BaseDom=" << m_base_dom << ")";
    }
  }

}; // class reference_domain

template <typename Params>
struct abstract_domain_traits<reference_domain<Params>> {
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;
};

} // end namespace domains
} // end namespace crab
