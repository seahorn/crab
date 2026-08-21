#pragma once

/**
 * JSON serialization for reference constraints.
 *
 * The reference counterpart of linear_constraints_to_json.hpp. A reference
 * constraint relates two reference variables, or a reference variable and the
 * null reference:
 *
 *     {"op": "==", "lhs": {var}, "rhs": {var}, "offset": "4"}   p == q + 4
 *     {"op": "!=", "lhs": {var}}                                p != null
 *     {"op": "true"}                                            null == null
 *     {"op": "false"}                                           null != null
 *
 * `op` is one of "==", "!=", "<", "<=", ">", ">=", or the nullary "true" /
 * "false" for tautologies and contradictions.
 *
 * A missing "rhs" means the null reference, which is why it is omitted rather
 * than emitted as null: `p == null` and `p == q` are different constraints, and
 * the absence of the key is what distinguishes them. "offset" is only present
 * on binary constraints, where it defaults to "0".
 */

#include <crab/support/debug.hpp>
#include <crab/support/json.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/variable_to_json.hpp>

namespace crab {
namespace json {

/** A reference constraint. */
template <typename Number, typename VariableName>
void write(writer &w,
           const crab::reference_constraint<Number, VariableName> &c) {
  w.begin_object();
  if (c.is_tautology()) {
    w.kv_string("op", "true");
    w.end_object();
    return;
  }
  if (c.is_contradiction()) {
    w.kv_string("op", "false");
    w.end_object();
    return;
  }

  const char *op = nullptr;
  if (c.is_equality()) {
    op = "==";
  } else if (c.is_disequality()) {
    op = "!=";
  } else if (c.is_less_than()) {
    op = "<";
  } else if (c.is_less_or_equal_than()) {
    op = "<=";
  } else if (c.is_greater_than()) {
    op = ">";
  } else if (c.is_greater_or_equal_than()) {
    op = ">=";
  } else {
    CRAB_ERROR("reference_constraints_to_json: unexpected constraint kind");
  }
  w.kv_string("op", op);

  // Neither operand is emitted when it is null: is_unary() means the rhs is
  // the null reference, and a constraint with no lhs is a tautology or a
  // contradiction, both handled above.
  w.key("lhs");
  write(w, c.lhs());
  if (c.is_binary()) {
    w.key("rhs");
    write(w, c.rhs());
    w.kv_number("offset", c.offset().get_str());
  }
  w.end_object();
}

} // end namespace json
} // end namespace crab
