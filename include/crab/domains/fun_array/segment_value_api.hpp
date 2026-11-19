#pragma once

#include <crab/support/os.hpp>

#include <utility>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Paper's A (POPL'11, §9): what holds for every element of one segment.
 *
 * Use of Curiously Recurring Template Pattern (CRTP): Value is the
 * implementation, ScalarDom the scalar domain used by the conversions
 * (§10) and the reduction.
 *
 * Required static member of Value (static members cannot be virtual):
 *
 *   static Value make_top_for(const variable_t &a);
 *
 * It returns the top value for the segments of array a, used when a is
 * first accessed.
 **/
template <class Value, class ScalarDom> class segment_value_api {
public:
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;

  virtual ~segment_value_api() {}

  // Lattice operations
  virtual Value make_bottom() const = 0;
  virtual Value make_top() const = 0;
  virtual void set_to_bottom() = 0;
  virtual void set_to_top() = 0;
  virtual bool is_bottom() const = 0;
  virtual bool is_top() const = 0;
  virtual bool operator<=(const Value &o) const = 0;
  virtual Value operator|(const Value &o) const = 0;
  virtual Value operator&(const Value &o) const = 0;
  // Widening
  virtual Value operator||(const Value &o) const = 0;
  // Narrowing
  virtual Value operator&&(const Value &o) const = 0;

  // Changes to the segmentation
  // [l,h) is cut into [l,e) and [e,h)
  virtual std::pair<Value, Value> split() const = 0;
  // [l,m) and [m,h) become [l,h)
  virtual Value merge(const Value &right) const = 0;
  // The one-element segment [i, i+1) cut out of this segment: the index
  // placeholder is fixed to i
  virtual Value singleton(const linear_expression_t &i) const = 0;

  // Conversions with the scalar domain (§10). Both are strong: they must be
  // called on a one-element segment [i, i+1) (see singleton). For a weak
  // update, join the segment with a copy that received the store; for a weak
  // read, load from a copy and discard it.
  // A[i] := val: the element becomes val, with what s knows about val
  virtual void store(const variable_t &a, const linear_expression_t &val,
                     const ScalarDom &s) = 0;
  // lhs := A[i]: s gets the element's constraints, and the segment records
  // that the element equals lhs. The caller forgets lhs first, in s and in
  // every segment value.
  virtual void load(const variable_t &a, const variable_t &lhs,
                    ScalarDom &s) = 0;

  // Does the value mention scalar x?
  virtual bool mentions(const variable_t &x) const = 0;
  // Drop the value placeholder of array a
  virtual void forget_array(const variable_t &a) = 0;

  // Scalar changes reaching the value
  virtual void forget(const variable_vector_t &scalars) = 0;
  // Forget every scalar not in scalars; arrays and the index are kept
  virtual void project(const variable_vector_t &scalars) = 0;
  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) = 0;
  // Meet with what s knows about the scalars this value mentions
  virtual void reduce(const ScalarDom &s) = 0;

  virtual void write(crab_os &o) const = 0;
};

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
