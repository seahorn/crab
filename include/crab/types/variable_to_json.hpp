#pragma once

/**
 * JSON serialization for Crab variable types and variables.
 *
 * A type is emitted as `{"kind": <string>}`, with an extra "bitwidth" for the
 * kinds that carry one:
 *
 *     {"kind": "int", "bitwidth": 32}
 *     {"kind": "bool"}
 *     {"kind": "int_array"}
 *
 * A variable is emitted as `{"name": <string>, "type": <type>}`.
 */

#include <crab/support/json.hpp>
#include <crab/types/variable.hpp>

namespace crab {
namespace json {

/** A variable type: {"kind": ..., ["bitwidth": ...]}. */
inline void write(writer &w, const crab::variable_type &ty) {
  w.begin_object(true /*compact*/);
  if (ty.is_bool()) {
    w.kv_string("kind", "bool");
  } else if (ty.is_integer()) {
    w.kv_string("kind", "int");
    w.kv_unsigned("bitwidth", ty.get_integer_bitwidth());
  } else if (ty.is_real()) {
    w.kv_string("kind", "real");
  } else if (ty.is_reference()) {
    w.kv_string("kind", "ref");
  } else if (ty.is_bool_array()) {
    w.kv_string("kind", "bool_array");
  } else if (ty.is_integer_array()) {
    w.kv_string("kind", "int_array");
  } else if (ty.is_real_array()) {
    w.kv_string("kind", "real_array");
  } else if (ty.is_unknown_region()) {
    w.kv_string("kind", "unknown_region");
  } else if (ty.is_bool_region()) {
    w.kv_string("kind", "bool_region");
  } else if (ty.is_integer_region()) {
    w.kv_string("kind", "int_region");
    w.kv_unsigned("bitwidth", ty.get_integer_region_bitwidth());
  } else if (ty.is_real_region()) {
    w.kv_string("kind", "real_region");
  } else if (ty.is_reference_region()) {
    w.kv_string("kind", "ref_region");
  } else if (ty.is_bool_array_region()) {
    w.kv_string("kind", "bool_array_region");
  } else if (ty.is_int_array_region()) {
    w.kv_string("kind", "int_array_region");
  } else if (ty.is_real_array_region()) {
    w.kv_string("kind", "real_array_region");
  } else {
    w.kv_string("kind", "unknown");
  }
  w.end_object();
}

/** A variable: {"name": ..., "type": ...}. */
template <typename Number, typename VariableName>
void write(writer &w, const crab::variable<Number, VariableName> &v) {
  w.begin_object(true /*compact*/);
  w.kv_string("name", to_string(v));
  w.key("type");
  write(w, v.get_type());
  w.end_object();
}

/**
 * A variable or a constant, discriminated by "kind":
 *
 *     {"kind": "var",   "name": "x", "type": {...}}
 *     {"kind": "const", "value": "4", "type": {...}}
 *
 * Both carry a type, so a consumer can treat the two uniformly.
 */
template <typename Number, typename VariableName>
void write(writer &w,
           const crab::variable_or_constant<Number, VariableName> &vc) {
  w.begin_object(true /*compact*/);
  if (vc.is_variable()) {
    w.kv_string("kind", "var");
    w.kv_string("name", to_string(vc.get_variable()));
  } else {
    w.kv_string("kind", "const");
    w.kv_number("value", vc.get_constant().get_str());
  }
  w.key("type");
  write(w, vc.get_type());
  w.end_object();
}

} // end namespace json
} // end namespace crab
