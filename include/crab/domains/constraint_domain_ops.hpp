/**
 * Operations at the boundary between an abstract domain and its linear
 * constraint representation. They go in opposite directions:
 *
 *   constraint_simplification: a constraint, sharpened by asking the
 *     domain, becomes a tighter constraint.
 *   constraint_extraction: a domain, queried about one variable, yields
 *     the constraints mentioning it.
 *
 * Neither is an abstract domain operation. Both are customization
 * points: a generic default here, specialized by the few domains that
 * can do better, and called by *other* domains in the middle of their
 * own operations.
 **/

#pragma once

#include <boost/optional.hpp>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>

namespace crab {
namespace domains {

// Only named by the static_asserts below, which reject the generic
// (type-erased) domains. Declaring them is enough for std::is_same.
template <typename Variable> class abstract_domain;
template <typename Variable> class abstract_domain_ref;

template <typename Domain> class constraint_simplification {
public:
  static_assert(
      !std::is_same<Domain, abstract_domain<typename Domain::variable_t>>::value,
      "constraint_simplification not supported for generic domain");
  static_assert(
      !std::is_same<Domain,
                    abstract_domain_ref<typename Domain::variable_t>>::value,
      "constraint_simplification not supported for generic domain");

  using number_t = typename Domain::number_t;
  using variable_t = typename Domain::variable_t;
  using linear_expression_t = typename Domain::linear_expression_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  using linear_constraint_system_t =
      typename Domain::linear_constraint_system_t;

  // Convert a disequality into a strict inequality:
  // - if cst is x!=y and abs_val |= x <= y then add x < y
  // - if cst is x!=y and abs_val |= x >= y then add x > y
  static void lower_disequality(const Domain &abs_val,
                                const linear_constraint_t &cst,
                                linear_constraint_system_t &out_csts) {

    // TODO: we could use abs_val to infer more disequalities from cst.
    auto get_binary_operands =
        [](const linear_constraint_t &c)
        -> boost::optional<std::pair<variable_t, variable_t>> {
      if (c.is_disequation()) {
        if (c.size() == 2 && c.constant() == 0) {
          auto it = c.begin();
          auto nx = it->first;
          auto vx = it->second;
          ++it;
          assert(it != c.end());
          auto ny = it->first;
          auto vy = it->second;
          if (nx == (ny * -1)) {
            return std::make_pair(vx, vy);
          }
        }
      }
      return boost::none;
    };

    if (auto pair = get_binary_operands(cst)) {
      variable_t x = (*pair).first;
      variable_t y = (*pair).second;
      linear_constraint_t x_le_y(x <= y);
      linear_constraint_t x_lt_y(linear_expression_t(x) < linear_expression_t(y));
      if (abs_val.entails(x_le_y)) {
        out_csts += x_lt_y;
      } else {
        linear_constraint_t x_ge_y(x >= y);
        linear_constraint_t x_gt_y(linear_expression_t(x) > linear_expression_t(y));
        if (abs_val.entails(x_ge_y)) {
          out_csts += x_gt_y;
        }
      }
    }
  }
};

// Extract from dom the linear constraints that mention x.
template <typename Domain> class constraint_extraction {
public:
  static_assert(
      !std::is_same<Domain, abstract_domain<typename Domain::variable_t>>::value,
      "constraint_extraction not supported for generic domain");
  static_assert(
      !std::is_same<Domain,
                    abstract_domain_ref<typename Domain::variable_t>>::value,
      "constraint_extraction not supported for generic domain");

  using variable_t = typename Domain::variable_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  using linear_constraint_system_t =
      typename Domain::linear_constraint_system_t;

  static void extract(Domain &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    auto all_csts = dom.to_linear_constraint_system();
    for (auto const &cst : all_csts) {
      if (only_equalities && (!cst.is_equality())) {
        continue;
      }
      if (std::find(std::begin(cst.variables()), std::end(cst.variables()), x) !=
          std::end(cst.variables())) {
        csts += cst;
      }
    }
  }
};

} // end namespace domains
} // end namespace crab
