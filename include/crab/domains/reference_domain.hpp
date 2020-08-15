#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
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
    Bottom, /*unused*/
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
       |  \ 
       |  [1,+oo]
      [0,1] /   
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

  bool operator==(small_range other) { return m_value == other.m_value; }

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
    // TODO
    CRAB_WARN("small_range::meet not implemented");
    return *this;
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

namespace reference_domain_impl {
template <class Number, class VariableName, class BaseAbsDom> class Params {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using varname_allocator_t = crab::var_factory_impl::str_var_alloc_col;
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
    : public abstract_domain_api<reference_domain<Params>> {
  using reference_domain_t = reference_domain<Params>;
  using abstract_domain_t = abstract_domain_api<reference_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using number_t = typename Params::number_t;
  using varname_t = typename Params::varname_t;

private:
  using base_abstract_domain_t = typename Params::base_abstract_domain_t;
  using base_variable_vector_t = typename base_abstract_domain_t::variable_vector_t;
  using base_variable_t = typename base_abstract_domain_t::variable_t;
  using base_linear_expression_t =
      typename base_abstract_domain_t::linear_expression_t;
  using base_linear_constraint_t =
      typename base_abstract_domain_t::linear_constraint_t;
  using base_linear_constraint_system_t =
      typename base_abstract_domain_t::linear_constraint_system_t;
  using base_varname_allocator_t = typename Params::varname_allocator_t;

  // Abstract domains
  using rgn_counting_domain_t =
      ikos::separate_domain<variable_t, small_range>;
  // Variable map: variable to base domain's variable.
  using var_map_t = std::unordered_map<variable_t, base_variable_t>;
  using rev_var_map_t = std::unordered_map<base_variable_t, variable_t>;

  bool m_is_bottom;

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
  // exactly one address per region (i.e., singleton).
  rgn_counting_domain_t m_rgn_counting_dom;
  // Abstract domain to keep track all possible regions for a
  // reference variable.
  //
  // The base abstract domain: all the heavy lifting is done here.
  // m_base_dom does not have any variable of reference type.
  base_abstract_domain_t m_base_dom;

  reference_domain(base_varname_allocator_t &&alloc,
                   var_map_t &&var_map, rev_var_map_t &&rev_var_map,
                   rgn_counting_domain_t &&rgn_counting_dom,
                   base_abstract_domain_t &&base_dom)
      : m_is_bottom(base_dom.is_bottom()), m_alloc(std::move(alloc)),
        m_var_map(std::move(var_map)), m_rev_var_map(std::move(rev_var_map)),
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_base_dom(std::move(base_dom)) {
  }

  using rgn_counting_binop_t = std::function<rgn_counting_domain_t(
      rgn_counting_domain_t, rgn_counting_domain_t)>;
  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;

  reference_domain_t do_join_or_widening(const reference_domain_t &left,
                                         const reference_domain_t &right,
                                         rgn_counting_binop_t rgn_counting_op,
                                         base_dom_binop_t base_dom_op) const {
    rgn_counting_domain_t out_rgn_counting_dom(
        rgn_counting_op(left.m_rgn_counting_dom, right.m_rgn_counting_dom));

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

    // JN: project might necessary to avoid keep variables that exist
    // in the base domain but they doen't exist on m_var_map.
    // 
    // If such a variable exists only on either left_dom or right_dom
    // the join removes it.  However, if we have the same variable in
    // both left_dom and righ_dom then the join will preserve it.
    // 
    left_dom.project(left_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.project(right_vars);
    right_dom.rename(right_vars, out_vars);
    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    reference_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_rgn_counting_dom),
	std::move(out_base_dom));
    return res;
  }

  reference_domain_t do_meet_or_narrowing(const reference_domain_t &left,
                                          const reference_domain_t &right,
                                          rgn_counting_binop_t rgn_counting_op,
                                          base_dom_binop_t base_dom_op) const {
    rgn_counting_domain_t out_rgn_counting_dom(
        rgn_counting_op(left.m_rgn_counting_dom, right.m_rgn_counting_dom));

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


    left_dom.rename(left_vars, out_vars);
    left_dom.rename(only_left_vars, only_left_out_vars);    
    right_dom.rename(right_vars, out_vars);
    right_dom.rename(only_right_vars, only_right_out_vars);
    
    base_abstract_domain_t out_base_dom(base_dom_op(left_dom, right_dom));

    reference_domain_t res(
        std::move(out_alloc), std::move(out_var_map),
        std::move(out_rev_var_map), std::move(out_rgn_counting_dom),
        std::move(out_base_dom));
    return res;
  }

  // Create a fresh variable in the base domain to shadow v
  base_variable_t make_base_variable(base_varname_allocator_t &var_allocator,
                                     const variable_t &v) const {
    if (v.is_ref_type()) {
      base_variable_t bv(var_allocator.next(), crab::INT_TYPE, 32);
      return bv;
    } else if (v.is_region_type()) {
      typename base_variable_t::type_t ty;
      unsigned bitwidth = 0;
      switch (v.get_type()) {
      case REG_BOOL_TYPE:
	ty = BOOL_TYPE;
	bitwidth = 1;
	break;
      case REG_ARR_BOOL_TYPE:
	ty = ARR_BOOL_TYPE;
	bitwidth = 1;
	break;
      case REG_INT_TYPE:
	ty = INT_TYPE;
	bitwidth = v.get_bitwidth();
	break;
      case REG_ARR_INT_TYPE:
	ty = ARR_INT_TYPE;
	bitwidth = v.get_bitwidth();
	break;
      case REG_REF_TYPE:
	ty = INT_TYPE;
	bitwidth = 32;
	break;	
      case REG_REAL_TYPE:
	ty = REAL_TYPE;
	break;
      case REG_ARR_REAL_TYPE:
	ty = ARR_REAL_TYPE;
	break;
      default:;
	CRAB_ERROR(domain_name() + "::make_base_variable: unreachable");
      }
      base_variable_t bv(var_allocator.next(), ty, bitwidth);
      return bv;
    } else {
      base_variable_t bv(var_allocator.next(), v.get_type(), v.get_bitwidth());
      return bv;
    }
  }

  bool is_initialized_region(const variable_t &rgn) const {
    small_range count_num = m_rgn_counting_dom[rgn];
    return (count_num <= small_range::oneOrMore());
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
                                             bool ignore_references) const {
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
    return base_linear_constraint_t(
        rename_linear_expr(cst.expression()),
        (typename base_linear_constraint_t::kind_t)cst.kind());
  }

  // used only for pretty-printing and to_linear_constraint_system()
  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const base_linear_constraint_t &cst,
                        bool ignore_references) const {
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
  reference_domain_t make_top() const override {
    return reference_domain_t(true);
  }

  reference_domain_t make_bottom() const override {
    return reference_domain_t(false);
  }

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
        m_var_map(o.m_var_map), m_rev_var_map(o.m_rev_var_map),
        m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_base_dom(o.m_base_dom) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  reference_domain(reference_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), m_alloc(std::move(o.m_alloc)),
        m_var_map(std::move(o.m_var_map)),
        m_rev_var_map(std::move(o.m_rev_var_map)),
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_base_dom(std::move(o.m_base_dom)) {}

  reference_domain_t &operator=(const reference_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_alloc = o.m_alloc;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_base_dom = o.m_base_dom;
    }
    return *this;
  }

  reference_domain_t &operator=(reference_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_alloc = std::move(o.m_alloc);
      m_var_map = std::move(o.m_var_map);
      m_rev_var_map = std::move(o.m_rev_var_map);
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_base_dom = std::move(o.m_base_dom);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    return m_base_dom.is_top() && m_rgn_counting_dom.is_top();
  }

  bool operator<=(const reference_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    }

    CRAB_LOG("reference",
	     crab::outs() << "Inclusion test:\n\t" << *this << "\n\t" << o << "\n";);
    
    if (!(m_rgn_counting_dom <= o.m_rgn_counting_dom)) {
      CRAB_LOG("reference", crab::outs() << "Result=0\n";);
      return false;
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

    left_dom.rename(left_vars, out_vars);
    right_dom.rename(right_vars, out_vars);
    bool res = left_dom <= right_dom;
    CRAB_LOG("reference", crab::outs() << "Result=" << res << "\n";);
    return res;
  }

  void operator|=(const reference_domain_t &o) override {
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

    CRAB_LOG("reference", crab::outs()
                              << "Join " << *this << " and " << o << "=");

    rgn_counting_domain_t out_rgn_counting_dom(m_rgn_counting_dom |
                                               o.m_rgn_counting_dom);
    base_varname_allocator_t out_alloc(m_alloc, o.m_alloc);
    base_abstract_domain_t left_dom(m_base_dom);
    base_abstract_domain_t right_dom(o.m_base_dom);
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

    // See comments in do_join_or_widening
    // 
    left_dom.project(left_vars);
    left_dom.rename(left_vars, out_vars);
    right_dom.project(right_vars);
    right_dom.rename(right_vars, out_vars);
    
    base_abstract_domain_t out_base_dom(left_dom | right_dom);

    m_is_bottom = out_base_dom.is_bottom();
    std::swap(m_alloc, out_alloc);
    std::swap(m_rgn_counting_dom, out_rgn_counting_dom);
    std::swap(m_var_map, out_var_map);
    std::swap(m_rev_var_map, out_rev_var_map);
    std::swap(m_base_dom, out_base_dom);

    CRAB_LOG("reference", crab::outs() << *this << "\n");
  }

  reference_domain_t operator|(const reference_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

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

    auto rgn_counting_op = [](const rgn_counting_domain_t &v1,
                              const rgn_counting_domain_t &v2) {
      return v1 | v2;
    };
    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 | v2; };
    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, rgn_counting_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator&(const reference_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Meet " << *this << " and " << o << "=");

    auto rgn_counting_op = [](const rgn_counting_domain_t &v1,
                              const rgn_counting_domain_t &v2) {
      return (v1 & v2);
    };
    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return (v1 & v2);
    };
    reference_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, rgn_counting_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator||(const reference_domain_t &o) const override {
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

    CRAB_LOG("reference", crab::outs()
                              << "Widening " << *this << " and " << o << "=");

    auto rgn_counting_op = [](const rgn_counting_domain_t &v1,
                              const rgn_counting_domain_t &v2) {
      return v1 || v2;
    };
    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 || v2;
    };
    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, rgn_counting_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t widening_thresholds(
      const reference_domain_t &o,
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

    CRAB_LOG("reference", crab::outs()
                              << "Widening " << *this << " and " << o << "=");

    auto rgn_counting_op = [](const rgn_counting_domain_t &v1,
                              const rgn_counting_domain_t &v2) {
      return v1 || v2;
    };
    auto base_dom_op = [&thresholds](const base_abstract_domain_t &v1,
                                     const base_abstract_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    reference_domain_t res(std::move(do_join_or_widening(
        *this, o, rgn_counting_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  reference_domain_t operator&&(const reference_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("reference", crab::outs()
                              << "Meet " << *this << " and " << o << "=");

    auto rgn_counting_op = [](rgn_counting_domain_t v1,
                              rgn_counting_domain_t v2) { return (v1 && v2); };
    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return (v1 && v2); };
    reference_domain_t res(std::move(do_meet_or_narrowing(
        *this, o, rgn_counting_op, base_dom_op)));

    CRAB_LOG("reference", crab::outs() << res << "\n");
    return res;
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (!is_bottom()) {
      if (v.is_region_type()) {
	m_rgn_counting_dom -= v;
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

    if (is_initialized_region(rgn)) {
      CRAB_ERROR("reference_domain::init_region: ", rgn,
                 " cannot be initialized twice");
    }

    // Set to zero the number of references
    m_rgn_counting_dom.set(rgn, small_range::zero());

    // Assign a synthetic variable to rgn for modeling its content.
    base_variable_t v = make_base_variable(m_alloc, rgn);
    CRAB_LOG("reference", crab::outs() << "After region_init(" << rgn
                                       << ")=" << *this << "\n";);
  }

  void region_assign(const variable_t &lhs_rgn, const variable_t &rhs_rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".region_assign");
    
    if (is_bottom()) {
      return;
    }

    if (!lhs_rgn.is_region_type()) {
      CRAB_ERROR(domain_name() + "::region_assign ", lhs_rgn, " must have a region type");      
    }
    if (!rhs_rgn.is_region_type()) {
      CRAB_ERROR(domain_name() + "::region_assign ", rhs_rgn, " must have a region type");            
    } 
    if (!lhs_rgn.same_type_and_bitwidth(rhs_rgn)) {
      CRAB_ERROR(domain_name() + "::region_assign ", lhs_rgn, ":=", rhs_rgn, " with different types");
    }
    
    m_rgn_counting_dom.set(lhs_rgn, m_rgn_counting_dom[rhs_rgn]);
    
    const base_variable_t &base_lhs = rename_var(lhs_rgn);
    const base_variable_t &base_rhs = rename_var(rhs_rgn);    
    switch (lhs_rgn.get_type()) {
    case REG_BOOL_TYPE:
      m_base_dom.assign_bool_var(base_lhs, base_rhs, false);
      break;
    case REG_INT_TYPE:
    case REG_REAL_TYPE:
    case REG_REF_TYPE:  
      m_base_dom.assign(base_lhs, base_rhs);
      break;
    case REG_ARR_BOOL_TYPE:
    case REG_ARR_INT_TYPE:
    case REG_ARR_REAL_TYPE:
      m_base_dom.array_assign(base_lhs, base_rhs);
      break;
    default:;
      CRAB_ERROR(domain_name() + 
          "::region_assign with unexpected region  type");
    }
  }
  
  // Create a new reference ref to region rgn.
  void ref_make(const variable_t &ref, const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    if (is_bottom()) {
      return;
    }

    if (!ref.is_ref_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_make: ", ref, " must be a reference"););
      return;
    }

    // Update reference counting
    auto num_refs = m_rgn_counting_dom[rgn];
    m_rgn_counting_dom.set(rgn, num_refs.increment());

    // Assign a base domain variable to ref
    rename_var(ref);

    CRAB_LOG("reference", crab::outs() << "After ref_make(" << ref << "," << rgn
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

    if (!rgn.is_scalar_region_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_load: ", rgn, " must be a scalar region"););
      m_base_dom -= base_res;
      return;
    }
    
    if (!ref.is_ref_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_load: ", ref, " must be a reference"););
      m_base_dom -= base_res;
      return;
    }

    if (is_null_ref(ref)) {
      CRAB_WARN("reference_domain::ref_load: reference ", ref, " is null. ",
                " Set to bottom ...");
      set_to_bottom();
      return;
    }

    auto ref_load = [&rgn, &base_res,this](const base_variable_t &rgn_var) {
        if (rgn_var.get_type() != base_res.get_type()) {
          CRAB_ERROR("reference_domain::ref_load: ", "Type of region ", rgn,
                     " does not match with ", base_res);
        }
        switch (base_res.get_type()) {
        case BOOL_TYPE:
          this->m_base_dom.assign_bool_var(base_res, rgn_var, false);
          break;
        case INT_TYPE:
        case REAL_TYPE:
          this->m_base_dom.assign(base_res, rgn_var);
          break;
        default:
          // variables of type base_variable_t cannot be REF_TYPE or
          // ARR_TYPE
          CRAB_ERROR("reference_domain::ref_load: unsupported type in ",
                     base_res);
        }
      };

    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_one()) {
      CRAB_LOG("reference", crab::outs() << "Reading from singleton\n";);
      if (auto region_var_opt = get_var(rgn)) {
	ref_load(*region_var_opt);
      } else {
	m_base_dom -= base_res;
      }
    } else {
      CRAB_LOG("reference", crab::outs() << "Reading from non-singleton\n";);
      if (auto region_var_opt = get_var(rgn)) {
	base_variable_t fresh_region_var = make_base_variable(m_alloc, rgn);
	m_base_dom.expand(*region_var_opt, fresh_region_var);
	ref_load(fresh_region_var);
      } else {
	m_base_dom -= base_res;
      }
    }
    CRAB_LOG("reference", crab::outs() << "After " << res << ":="
                                       << "ref_load(" << ref << "," << rgn
                                       << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in region
  // rgn.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    if (is_bottom()) {
      return;
    }

    if (!rgn.is_scalar_region_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_store: ", rgn, " must be a scalar region"););
      return;
    }

    if (!ref.is_ref_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_store: ", ref, " must be a reference"););
      return;
    }

    if (is_null_ref(ref)) {
      CRAB_LOG("reference",
               CRAB_WARN("reference_domain::ref_store: rerefence ", ref,
                         " is null. ", " Set to bottom ..."););
      set_to_bottom();
      return;
    }
    
    auto ref_store = [&rgn, &val,this](base_abstract_domain_t &base_dom) {
	base_variable_t rgn_var = rename_var(rgn);
        switch (rgn_var.get_type()) {
        case BOOL_TYPE:
          if (val.is_constant()) {
            if (val.constant() >= number_t(1)) {
              base_dom.assign_bool_cst(rgn_var,
                                       base_linear_constraint_t::get_true());
            } else {
              base_dom.assign_bool_cst(rgn_var,
                                       base_linear_constraint_t::get_false());
            }
          } else if (auto var = val.get_variable()) {
            base_dom.assign_bool_var(rgn_var, rename_var(*var), false);
          } else {
            CRAB_ERROR("reference_domain::ref_store: unsupported boolean ",
                       val);
          }
          break;
        case INT_TYPE:
        case REAL_TYPE:
          base_dom.assign(rgn_var, rename_linear_expr(val));
          break;
        default:
          // variables of type base_variable_t cannot be REF_TYPE or
          // ARR_TYPE
          CRAB_ERROR("reference_domain::ref_store: unsupported type in ",
                     rgn_var);
        }
      };

    auto num_refs = m_rgn_counting_dom[rgn];    
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
    CRAB_LOG("reference", crab::outs()
                              << "After ref_store(" << ref << "," << rgn << ","
                              << val << ")=" << *this << "\n";);
  }

  // Create a new reference ref2 to region rgn2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &rgn1,
               const variable_t &ref2, const variable_t &rgn2,
               const linear_expression_t &offset) override {
    crab::CrabStats::count(domain_name() + ".count.ref_gep");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_gep");
    
    if (!is_bottom()) {
      m_base_dom.assign(rename_var(ref2),
			rename_var(ref1) + rename_linear_expr(offset));

      if (!(rgn1 == rgn2 && offset.equal(number_t(0)))) {
	// Update reference counting
	auto num_refs = m_rgn_counting_dom[rgn2];
	m_rgn_counting_dom.set(rgn2, num_refs.increment());
      }
    }

    CRAB_LOG("reference", crab::outs()
	     << "After (" << rgn2 << "," << ref2 <<") := ref_gep(" << rgn1 << "," << ref1
	     << " + " << offset << ")=" << *this << "\n";);
  }

  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref /*unused*/,
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
    if (!lhs.is_int_type() &&
	!lhs.is_bool_type() &&
	!lhs.is_real_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_load_from_array: ", lhs, " must be bool/int/real type"););
      m_base_dom -= base_lhs;

    }
    
    if (!rgn.is_array_region_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_load_from_array: ", rgn, " must be an array region"););
      m_base_dom -= base_lhs;
      return;
    }

    // TODO: check that the array element matches the type of lhs

    // if (!ref.is_ref_type()) {
    //   CRAB_LOG("reference",
    // 	       CRAB_WARN("reference_domain::ref_load_from_array: ", ref, " must be a reference"););
    //   m_base_dom -= base_lhs;
    //   return;
    // }

    // if (is_null_ref(ref)) {
    //   CRAB_WARN("reference_domain::ref_load_from_array: reference ", ref, " is null. ",
    //             " Set to bottom ...");
    //   set_to_bottom();
    //   return;
    // }

    if (boost::optional<base_variable_t> arr_var_opt = get_var(rgn)) {
      auto num_refs = m_rgn_counting_dom[rgn];
      if (num_refs.is_one()) {
	CRAB_LOG("reference", crab::outs() << "Reading from singleton\n";);
	m_base_dom.array_load(base_lhs, *arr_var_opt, rename_linear_expr(elem_size),
			      rename_linear_expr(index));
      } else {
	// TODO: get the bitwidth from index
	base_variable_t unknown_index(m_alloc.next(), crab::INT_TYPE, 32);
	m_base_dom.array_load(base_lhs, *arr_var_opt, rename_linear_expr(elem_size), unknown_index);
      }
    } else {
      m_base_dom -= base_lhs;
    }
  }

  // Treat region as an array and perform an array store.
  void ref_store_to_array(const variable_t &ref /*unused*/, const variable_t &rgn,
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
    if (!rgn.is_array_region_type()) {
      CRAB_LOG("reference",
	       CRAB_WARN("reference_domain::ref_store_to_array: ", rgn, " must be an array region"););
      return;
    }

    // TODO: check that the array element matches the type of val

    // if (!ref.is_ref_type()) {
    //   CRAB_LOG("reference",
    // 	       CRAB_WARN("reference_domain::ref_store_to_array: ", ref, " must be a reference"););
    //   return;
    // }

    // if (is_null_ref(ref)) {
    //   CRAB_WARN("reference_domain::ref_store_to_array: reference ", ref, " is null. ",
    //             " Set to bottom ...");
    //   set_to_bottom();
    //   return;
    // }

    auto num_refs = m_rgn_counting_dom[rgn];
    base_variable_t arr_var = rename_var(rgn);
    if (num_refs.is_one()) {
      CRAB_LOG("reference", crab::outs() << "Reading from singleton\n";);
      m_base_dom.array_store(arr_var, rename_linear_expr(elem_size),
			     rename_linear_expr(index), rename_linear_expr(val), false);
    } else {
      // TODO: get the bitwidth from index
      base_variable_t unknown_index(m_alloc.next(), crab::INT_TYPE, 32);
      m_base_dom.array_store(arr_var, rename_linear_expr(elem_size),
			     unknown_index, rename_linear_expr(val), false);
    }
  }

  void ref_assume(const reference_constraint_t &cst) override {
    crab::CrabStats::count(domain_name() + ".count.ref_assume");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_assume");

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
	base_variable_t x = rename_var(cst.lhs());
	if (cst.is_equality()) {
	  m_base_dom += base_linear_constraint_t(x == number_t(0));
	} else if (cst.is_disequality()) {
	  m_base_dom += base_linear_constraint_t(x > number_t(0));
	} else if (cst.is_less_or_equal_than()) {
            m_base_dom += base_linear_constraint_t(x <= number_t(0));
	} else if (cst.is_less_than()) {
	  m_base_dom += base_linear_constraint_t(x < number_t(0));
	} else if (cst.is_greater_or_equal_than()) {
	  m_base_dom += base_linear_constraint_t(x >= number_t(0));
	} else if (cst.is_greater_than()) {
	  m_base_dom += base_linear_constraint_t(x > number_t(0));
	}
      } else {
        assert(cst.lhs().is_ref_type());
        assert(cst.rhs().is_ref_type());
	base_variable_t x = rename_var(cst.lhs());
	base_variable_t y = rename_var(cst.rhs());
        number_t offset = cst.offset();
	if (cst.is_equality()) {
	  m_base_dom += base_linear_constraint_t(x == y + offset);
	} else if (cst.is_disequality()) {
	  m_base_dom += base_linear_constraint_t(x != y + offset);
	} else if (cst.is_less_or_equal_than()) {
	  m_base_dom += base_linear_constraint_t(x <= y + offset);
	} else if (cst.is_less_than()) {
	  m_base_dom += base_linear_constraint_t(x < y + offset);
	} else if (cst.is_greater_or_equal_than()) {
	  m_base_dom += base_linear_constraint_t(x >= y + offset);
	} else if (cst.is_greater_than()) {
	  m_base_dom += base_linear_constraint_t(x > y + offset);
	}
      }
    }
    m_is_bottom = m_base_dom.is_bottom();
    CRAB_LOG("reference", crab::outs()
	     << "ref_assume(" << cst << ")" << *this << "\n";);
  }

  void ref_to_int(const variable_t &rgn, const variable_t &ref_var,
		  const variable_t &int_var) override {
    crab::CrabStats::count(domain_name() + ".count.ref_to_int");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_to_int");
    if (!is_bottom()) {
      CRAB_WARN(domain_name() + "::ref_to_int not implemented");
    }

  }
  void int_to_ref(const variable_t &int_var,
		  const variable_t &rgn, const variable_t &ref_var) override {
    crab::CrabStats::count(domain_name() + ".count.int_to_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".int_to_ref");
    if (!is_bottom()) {
      CRAB_WARN(domain_name() + "::int_to_ref not implemented");
    }
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
      m_base_dom.assign(rename_var(x), rename_linear_expr(e));
    }
  }
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    if (!is_bottom()) {
      m_base_dom.backward_assign(rename_var(x), rename_linear_expr(e),
                                 invariant.m_base_dom);
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, rename_var(x), rename_var(y), z,
                                invariant.m_base_dom);
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const reference_domain_t &invariant) override {
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
      m_base_dom += rename_linear_cst_sys(csts);
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
      m_base_dom.assign_bool_cst(rename_var(lhs), rename_linear_cst(rhs));
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
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_cst(
          rename_var(lhs), rename_linear_cst(rhs), invariant.m_base_dom);
    }
  }
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_var(rename_var(lhs), rename_var(rhs),
                                          is_not_rhs, invariant.m_base_dom);
    }
  }
  void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const reference_domain_t &invariant) override {
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
      m_base_dom.array_init(rename_var(a), rename_linear_expr(elem_size),
                            rename_linear_expr(lb_idx),
                            rename_linear_expr(ub_idx),
                            rename_linear_expr(val));
    }
  }
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".array_load");

    if (!is_bottom()) {
      m_base_dom.array_load(rename_var(lhs), rename_var(a),
                            rename_linear_expr(elem_size),
                            rename_linear_expr(i));
    }
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store");

    if (!is_bottom()) {
      m_base_dom.array_store(rename_var(a), rename_linear_expr(elem_size),
                             rename_linear_expr(i), rename_linear_expr(v),
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
      m_base_dom.array_store_range(rename_var(a), rename_linear_expr(elem_size),
                                   rename_linear_expr(i), rename_linear_expr(j),
                                   rename_linear_expr(v));
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
                           const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_init");

    if (!is_bottom()) {
      m_base_dom.backward_array_init(
          rename_var(a), rename_linear_expr(elem_size),
          rename_linear_expr(lb_idx), rename_linear_expr(ub_idx),
          rename_linear_expr(val), invariant.m_base_dom);
    }
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_load");

    if (!is_bottom()) {
      m_base_dom.backward_array_load(
          rename_var(lhs), rename_var(a), rename_linear_expr(elem_size),
          rename_linear_expr(i), invariant.m_base_dom);
    }
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store");

    if (!is_bottom()) {
      m_base_dom.backward_array_store(
          rename_var(a), rename_linear_expr(elem_size), rename_linear_expr(i),
          rename_linear_expr(v), is_strong_update, invariant.m_base_dom);
    }
  }
  void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const reference_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store_range");

    if (!is_bottom()) {
      m_base_dom.backward_array_store_range(
          rename_var(a), rename_linear_expr(elem_size), rename_linear_expr(i),
          rename_linear_expr(j), rename_linear_expr(v), invariant.m_base_dom);
    }
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const reference_domain_t &invariant) override {
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
    for (auto &v: variables) {
      if (v.is_region_type()) {
	m_rgn_counting_dom -= v;
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
      base_vars.push_back(rename_var(v));
    }
    m_base_dom.project(base_vars);

    // -- update m_var_map and m_rev_var_map
    std::vector<variable_t> var_map_to_remove;
    var_map_to_remove.reserve(m_var_map.size());
    for (auto &kv: m_var_map) {
      if (!std::binary_search(sorted_variables.begin(), sorted_variables.end(), kv.first)) {
	var_map_to_remove.push_back(kv.first);
	m_rev_var_map.erase(kv.second);
      }
    }
    for (auto &v: var_map_to_remove) {
      m_var_map.erase(v);
    }

    // -- update m_rgn_counting_dom
    if (!m_rgn_counting_dom.is_top() && !m_rgn_counting_dom.is_bottom()) {
      std::vector<variable_t> rgn_counting_vars;
      rgn_counting_vars.reserve(m_rgn_counting_dom.size());
      for (auto kv: m_rgn_counting_dom) {
	if (!std::binary_search(sorted_variables.begin(), sorted_variables.end(), kv.first)) {            
	  rgn_counting_vars.push_back(kv.first);
	}
      }
      for (auto &v: rgn_counting_vars) {
	m_rgn_counting_dom -= v;
      }
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
    CRAB_ERROR("reference_domain::expand not implemented");
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    CRAB_ERROR("reference_domain::rename not implemented");
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const reference_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
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
    CRAB_ERROR("reference_domain::to_disjunctive_linear_constraint_system not "
               "implemented");
  }

  std::string domain_name() const override {
    return "ReferenceDomain(" + m_base_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("reference-print",      
	       o << "(RgnCounter=" << m_rgn_counting_dom << ","
	       << "MapVar={";
	       for (auto &kv : m_var_map) {
		 o << kv.first << "->" << kv.second << ";";
	       }
	       o << "}," << "BaseDom=" << m_base_dom << ")\n";
	       );
      std::unordered_map<std::string, std::string> renaming_map;
      for (auto &kv: m_rev_var_map) {
	renaming_map[kv.first.name().str()] = kv.second.name().str();
      }
      m_alloc.add_renaming_map(renaming_map);
      o << m_base_dom;
      m_alloc.clear_renaming_map();

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
