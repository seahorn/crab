#pragma once

/*
   Implementation of the abstract transfer functions by reducing them
   to abstract domain operations.

   These are the main Crab statements for which we define their abstract
   transfer functions:

   ARITHMETIC and BOOLEAN
     x := y bin_op z;
     x := y;
     assume(cst)
     assert(cst);
     x := select(cond, y, z);

   ARRAYS
     a[l...u] := v (a,b are arrays and v can be bool or integer)
     a[i] := v;
     v := a[i];
     a := b

   REGIONS and REFERENCES
     TODO

   FUNCTIONS
     x := foo(arg1,...,argn);
     return r;

   havoc(x);

 */

#include <crab/cfg/cfg.hpp>
#include <crab/domains/abstract_domain_operators.hpp>
#include <boost/optional.hpp>

namespace crab {
namespace analyzer {

/**
 * API abstract transformer
 **/
template <typename BasicBlockLabel, typename Number, typename VariableName>
class abs_transformer_api
    : public crab::cfg::const_statement_visitor<BasicBlockLabel, Number,
                                                VariableName> {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using bb_label_t = BasicBlockLabel;

  using var_t = variable<number_t, varname_t>;
  using var_or_cst_t = variable_or_constant<number_t, varname_t>;
  using lin_exp_t = ikos::linear_expression<number_t, varname_t>;
  using lin_cst_t = ikos::linear_constraint<number_t, varname_t>;
  using lin_cst_sys_t = ikos::linear_constraint_system<number_t, varname_t>;
  using ref_cst_t = reference_constraint<number_t, varname_t>;

  using havoc_t = crab::cfg::havoc_stmt<bb_label_t, number_t, varname_t>;
  using unreach_t =
      crab::cfg::unreachable_stmt<bb_label_t, number_t, varname_t>;

  using bin_op_t = crab::cfg::binary_op<bb_label_t, number_t, varname_t>;
  using assign_t = crab::cfg::assignment<bb_label_t, number_t, varname_t>;
  using assume_t = crab::cfg::assume_stmt<bb_label_t, number_t, varname_t>;
  using select_t = crab::cfg::select_stmt<bb_label_t, number_t, varname_t>;
  using assert_t = crab::cfg::assert_stmt<bb_label_t, number_t, varname_t>;
  using int_cast_t = crab::cfg::int_cast_stmt<bb_label_t, number_t, varname_t>;

  using callsite_t = crab::cfg::callsite_stmt<bb_label_t, number_t, varname_t>;
  using intrinsic_t =
      crab::cfg::intrinsic_stmt<bb_label_t, number_t, varname_t>;

  using arr_init_t =
      crab::cfg::array_init_stmt<bb_label_t, number_t, varname_t>;
  using arr_store_t =
      crab::cfg::array_store_stmt<bb_label_t, number_t, varname_t>;
  using arr_load_t =
      crab::cfg::array_load_stmt<bb_label_t, number_t, varname_t>;
  using arr_assign_t =
      crab::cfg::array_assign_stmt<bb_label_t, number_t, varname_t>;

  using region_init_t =
      crab::cfg::region_init_stmt<bb_label_t, number_t, varname_t>;
  using region_copy_t =
      crab::cfg::region_copy_stmt<bb_label_t, number_t, varname_t>;
  using region_cast_t =
      crab::cfg::region_cast_stmt<bb_label_t, number_t, varname_t>;  
  using make_ref_t = crab::cfg::make_ref_stmt<bb_label_t, number_t, varname_t>;
  using remove_ref_t = crab::cfg::remove_ref_stmt<bb_label_t, number_t, varname_t>;  
  using load_from_ref_t =
      crab::cfg::load_from_ref_stmt<bb_label_t, number_t, varname_t>;
  using store_to_ref_t =
      crab::cfg::store_to_ref_stmt<bb_label_t, number_t, varname_t>;
  using gep_ref_t = crab::cfg::gep_ref_stmt<bb_label_t, number_t, varname_t>;
  using assume_ref_t =
      crab::cfg::assume_ref_stmt<bb_label_t, number_t, varname_t>;
  using assert_ref_t =
      crab::cfg::assert_ref_stmt<bb_label_t, number_t, varname_t>;
  using select_ref_t =
      crab::cfg::ref_select_stmt<bb_label_t, number_t, varname_t>;
  using ref_to_int_t =
      crab::cfg::ref_to_int_stmt<bb_label_t, number_t, varname_t>;
  using int_to_ref_t =
      crab::cfg::int_to_ref_stmt<bb_label_t, number_t, varname_t>;
  using bool_bin_op_t =
      crab::cfg::bool_binary_op<bb_label_t, number_t, varname_t>;
  using bool_assign_cst_t =
      crab::cfg::bool_assign_cst<bb_label_t, number_t, varname_t>;
  using bool_assign_var_t =
      crab::cfg::bool_assign_var<bb_label_t, number_t, varname_t>;
  using bool_assume_t =
      crab::cfg::bool_assume_stmt<bb_label_t, number_t, varname_t>;
  using bool_select_t =
      crab::cfg::bool_select_stmt<bb_label_t, number_t, varname_t>;
  using bool_assert_t =
      crab::cfg::bool_assert_stmt<bb_label_t, number_t, varname_t>;

protected:
  virtual void exec(const havoc_t &) {}
  virtual void exec(const unreach_t &) {}
  virtual void exec(const bin_op_t &) {}
  virtual void exec(const assign_t &) {}
  virtual void exec(const assume_t &) {}
  virtual void exec(const select_t &) {}
  virtual void exec(const assert_t &) {}
  virtual void exec(const int_cast_t &) {}
  virtual void exec(const callsite_t &) {}
  virtual void exec(const intrinsic_t &) {}
  virtual void exec(const arr_init_t &) {}
  virtual void exec(const arr_store_t &) {}
  virtual void exec(const arr_load_t &) {}
  virtual void exec(const arr_assign_t &) {}
  virtual void exec(const region_init_t &) {}
  virtual void exec(const region_copy_t &) {}
  virtual void exec(const region_cast_t &) {}  
  virtual void exec(const make_ref_t &) {}
  virtual void exec(const remove_ref_t &) {}  
  virtual void exec(const load_from_ref_t &) {}
  virtual void exec(const store_to_ref_t &) {}
  virtual void exec(const gep_ref_t &) {}
  virtual void exec(const assume_ref_t &) {}
  virtual void exec(const assert_ref_t &) {}
  virtual void exec(const select_ref_t &) {}
  virtual void exec(const int_to_ref_t &) {}
  virtual void exec(const ref_to_int_t &) {}
  virtual void exec(const bool_bin_op_t &) {}
  virtual void exec(const bool_assign_cst_t &) {}
  virtual void exec(const bool_assign_var_t &) {}
  virtual void exec(const bool_assume_t &) {}
  virtual void exec(const bool_select_t &) {}
  virtual void exec(const bool_assert_t &) {}

public: /* visitor api */
  virtual void visit(const havoc_t &s) override { exec(s); }
  virtual void visit(const unreach_t &s) override { exec(s); }
  virtual void visit(const bin_op_t &s) override { exec(s); }
  virtual void visit(const assign_t &s) override { exec(s); }
  virtual void visit(const assume_t &s) override { exec(s); }
  virtual void visit(const select_t &s) override { exec(s); }
  virtual void visit(const assert_t &s) override { exec(s); }
  virtual void visit(const int_cast_t &s) override { exec(s); }
  virtual void visit(const callsite_t &s) override { exec(s); }
  virtual void visit(const intrinsic_t &s) override { exec(s); }
  virtual void visit(const arr_init_t &s) override { exec(s); }
  virtual void visit(const arr_store_t &s) override { exec(s); }
  virtual void visit(const arr_load_t &s) override { exec(s); }
  virtual void visit(const arr_assign_t &s) override { exec(s); }
  virtual void visit(const region_init_t &s) override { exec(s); }
  virtual void visit(const region_copy_t &s) override { exec(s); }
  virtual void visit(const region_cast_t &s) override { exec(s); }  
  virtual void visit(const make_ref_t &s) override { exec(s); }
  virtual void visit(const remove_ref_t &s) override { exec(s); }  
  virtual void visit(const load_from_ref_t &s) override { exec(s); }
  virtual void visit(const store_to_ref_t &s) override { exec(s); }
  virtual void visit(const gep_ref_t &s) override { exec(s); }
  virtual void visit(const assume_ref_t &s) override { exec(s); }
  virtual void visit(const assert_ref_t &s) override { exec(s); }
  virtual void visit(const select_ref_t &s) override { exec(s); }
  virtual void visit(const ref_to_int_t &s) override { exec(s); }
  virtual void visit(const int_to_ref_t &s) override { exec(s); }
  virtual void visit(const bool_bin_op_t &s) override { exec(s); }
  virtual void visit(const bool_assign_cst_t &s) override { exec(s); }
  virtual void visit(const bool_assign_var_t &s) override { exec(s); }
  virtual void visit(const bool_assume_t &s) override { exec(s); }
  virtual void visit(const bool_select_t &s) override { exec(s); }
  virtual void visit(const bool_assert_t &s) override { exec(s); }
};

/**
 * Convert CFG operations into abstract domain operations
 **/

template <typename T>
inline boost::optional<T> conv_op(crab::cfg::binary_operation_t op);
template <typename T>
inline boost::optional<T> conv_op(crab::cfg::bool_binary_operation_t op);
template <typename T>
inline boost::optional<T> conv_op(crab::cfg::cast_operation_t op);

template <>
inline boost::optional<domains::arith_operation_t>
conv_op(crab::cfg::binary_operation_t op) {
  switch (op) {
  case crab::cfg::BINOP_ADD:
    return domains::OP_ADDITION;
  case crab::cfg::BINOP_SUB:
    return domains::OP_SUBTRACTION;
  case crab::cfg::BINOP_MUL:
    return domains::OP_MULTIPLICATION;
  case crab::cfg::BINOP_SDIV:
    return domains::OP_SDIV;
  case crab::cfg::BINOP_UDIV:
    return domains::OP_UDIV;
  case crab::cfg::BINOP_SREM:
    return domains::OP_SREM;
  case crab::cfg::BINOP_UREM:
    return domains::OP_UREM;
  default:
    return boost::optional<domains::arith_operation_t>();
  }
}

template <>
inline boost::optional<domains::bitwise_operation_t>
conv_op(crab::cfg::binary_operation_t op) {
  switch (op) {
  case crab::cfg::BINOP_AND:
    return domains::OP_AND;
  case crab::cfg::BINOP_OR:
    return domains::OP_OR;
  case crab::cfg::BINOP_XOR:
    return domains::OP_XOR;
  case crab::cfg::BINOP_SHL:
    return domains::OP_SHL;
  case crab::cfg::BINOP_LSHR:
    return domains::OP_LSHR;
  case crab::cfg::BINOP_ASHR:
    return domains::OP_ASHR;
  default:
    return boost::optional<domains::bitwise_operation_t>();
  }
}

template <>
inline boost::optional<domains::int_conv_operation_t>
conv_op(crab::cfg::cast_operation_t op) {
  switch (op) {
  case crab::cfg::CAST_TRUNC:
    return domains::OP_TRUNC;
  case crab::cfg::CAST_SEXT:
    return domains::OP_SEXT;
  default:
    //case crab::cfg::CAST_ZEXT:
    return domains::OP_ZEXT;
  }
}

template <>
inline boost::optional<domains::bool_operation_t>
conv_op(crab::cfg::bool_binary_operation_t op) {
  switch (op) {
  case crab::cfg::BINOP_BAND:
    return domains::OP_BAND;
  case crab::cfg::BINOP_BOR:
    return domains::OP_BOR;
  default:
    //case crab::cfg::BINOP_BXOR:
    return domains::OP_BXOR;
  }
}
} // namespace analyzer
} // namespace crab

