#pragma once

/**
 * JSON serialization for linear expressions, constraints, and (disjunctive)
 * constraint systems.
 *
 * These are the form in which every Crab abstract domain can export its
 * invariants (`to_linear_constraint_system` /
 * `to_disjunctive_linear_constraint_system`), so this is the domain-independent
 * way to get an abstract state out of Crab in machine-readable form.
 *
 * Encoding
 * --------
 * A linear expression is emitted as its variable terms plus a constant:
 *
 *     {"terms": [["2", "x"], ["-1", "y"]], "const": "3"}     <=>   2x - y + 3
 *
 * A constraint is emitted with the constant moved to the right-hand side, so
 * that it reads directly as `<terms> <op> <const>`, together with the type of
 * its variables:
 *
 *     {"op": "<=", "type": {"kind": "int", "bitwidth": 32},
 *      "terms": [["1", "y"]], "const": "9"}                  <=>   y <= 9
 *
 * `op` is one of "<=", "<", "=", "!=", or the nullary "true" / "false" for
 * tautologies and contradictions (in which case the other keys are absent).
 *
 * Crab type-checks constraints, so every variable occurring in one has the same
 * type and a single "type" per constraint suffices. It is `null` in the corner
 * case of a constraint with no variables that is neither a tautology nor a
 * contradiction.
 *
 * All numbers are emitted as decimal *strings*: Crab's numbers are
 * arbitrary-precision and JSON numbers are not.
 */

#include <crab/support/json.hpp>
#include <crab/types/linear_constraints.hpp>
#include <crab/types/variable_to_json.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace json {

/**
 * The common type of the variables in a linear expression, or none if the
 * expression is a constant. Crab type-checks expressions, so the first
 * variable's type is the type of all of them.
 */
template <typename Number, typename VariableName>
boost::optional<crab::variable_type>
type_of(const ikos::linear_expression<Number, VariableName> &e) {
  auto it = e.begin();
  if (it == e.end()) {
    return boost::optional<crab::variable_type>();
  }
  return boost::optional<crab::variable_type>(it->second.get_type());
}

/** A linear expression: {"terms": [[coef, var], ...], "const": c}. */
template <typename Number, typename VariableName>
void write(writer &w, const ikos::linear_expression<Number, VariableName> &e) {
  w.begin_object();
  w.key("type");
  if (auto ty = type_of(e)) {
    write(w, *ty);
  } else {
    w.value_null();
  }
  w.key("terms");
  w.begin_array();
  for (auto it = e.begin(), et = e.end(); it != et; ++it) {
    w.begin_array(true /*compact*/);
    w.value_number(it->first.get_str());
    w.value_string(to_string(it->second));
    w.end_array();
  }
  w.end_array();
  w.kv_number("const", e.constant().get_str());
  w.end_object();
}

/**
 * A linear constraint, normalized to `<terms> <op> <const>`.
 *
 * Crab stores a constraint as `expr <op> 0` with the constant folded into
 * `expr`; we move it to the right-hand side, which is both more readable and
 * closer to what consumers want.
 */
template <typename Number, typename VariableName>
void write(writer &w, const ikos::linear_constraint<Number, VariableName> &c) {
  using linear_constraint_t = ikos::linear_constraint<Number, VariableName>;

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
  switch (c.kind()) {
  case linear_constraint_t::INEQUALITY:
    op = "<=";
    break;
  case linear_constraint_t::STRICT_INEQUALITY:
    op = "<";
    break;
  case linear_constraint_t::EQUALITY:
    op = "=";
    break;
  case linear_constraint_t::DISEQUATION:
    op = "!=";
    break;
  }
  w.kv_string("op", op);
  w.key("type");
  if (auto ty = type_of(c.expression())) {
    write(w, *ty);
  } else {
    w.value_null();
  }
  w.key("terms");
  w.begin_array();
  for (auto it = c.begin(), et = c.end(); it != et; ++it) {
    w.begin_array(true /*compact*/);
    w.value_number(it->first.get_str());
    w.value_string(to_string(it->second));
    w.end_array();
  }
  w.end_array();
  // Crab's representation is `expr <op> 0` with the constant inside expr, so
  // the right-hand side is -constant.
  w.kv_number("const", (-c.expression().constant()).get_str());
  w.end_object();
}

/** A conjunction of constraints: [cst, ...]. */
template <typename Number, typename VariableName>
void write(writer &w,
           const ikos::linear_constraint_system<Number, VariableName> &csts) {
  w.begin_array();
  for (auto const &c : csts) {
    write(w, c);
  }
  w.end_array();
}

/**
 * A disjunctive constraint system, i.e. an abstract state:
 *
 *   {"kind": "false"}                          bottom
 *   {"kind": "true"}                           top
 *   {"kind": "disj", "disjuncts": [[cst,...], ...]}
 *
 * Bottom and top are tagged explicitly rather than being represented by an
 * empty list, because "no disjuncts" and "no constraints" would otherwise be
 * indistinguishable.
 */
template <typename Number, typename VariableName>
void write(writer &w, const ikos::disjunctive_linear_constraint_system<
                          Number, VariableName> &dcsts) {
  w.begin_object();
  if (dcsts.is_false()) {
    w.kv_string("kind", "false");
  } else if (dcsts.is_true()) {
    w.kv_string("kind", "true");
  } else {
    w.kv_string("kind", "disj");
    w.key("disjuncts");
    w.begin_array();
    for (auto const &csts : dcsts) {
      write(w, csts);
    }
    w.end_array();
  }
  w.end_object();
}

} // end namespace json
} // end namespace crab
