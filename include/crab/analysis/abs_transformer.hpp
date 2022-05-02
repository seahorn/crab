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
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace analyzer {

/**
 * API abstract transformer
 **/
template <typename BasicBlockLabel, typename Number, typename VariableName>
class abs_transformer_api
    : public crab::cfg::statement_visitor<BasicBlockLabel, Number,
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
  virtual void exec(havoc_t &) {}
  virtual void exec(unreach_t &) {}
  virtual void exec(bin_op_t &) {}
  virtual void exec(assign_t &) {}
  virtual void exec(assume_t &) {}
  virtual void exec(select_t &) {}
  virtual void exec(assert_t &) {}
  virtual void exec(int_cast_t &) {}
  virtual void exec(callsite_t &) {}
  virtual void exec(intrinsic_t &) {}
  virtual void exec(arr_init_t &) {}
  virtual void exec(arr_store_t &) {}
  virtual void exec(arr_load_t &) {}
  virtual void exec(arr_assign_t &) {}
  virtual void exec(region_init_t &) {}
  virtual void exec(region_copy_t &) {}
  virtual void exec(region_cast_t &) {}  
  virtual void exec(make_ref_t &) {}
  virtual void exec(remove_ref_t &) {}  
  virtual void exec(load_from_ref_t &) {}
  virtual void exec(store_to_ref_t &) {}
  virtual void exec(gep_ref_t &) {}
  virtual void exec(assume_ref_t &) {}
  virtual void exec(assert_ref_t &) {}
  virtual void exec(select_ref_t &) {}
  virtual void exec(int_to_ref_t &) {}
  virtual void exec(ref_to_int_t &) {}
  virtual void exec(bool_bin_op_t &) {}
  virtual void exec(bool_assign_cst_t &) {}
  virtual void exec(bool_assign_var_t &) {}
  virtual void exec(bool_assume_t &) {}
  virtual void exec(bool_select_t &) {}
  virtual void exec(bool_assert_t &) {}

public: /* visitor api */
  virtual void visit(havoc_t &s) override { exec(s); }
  virtual void visit(unreach_t &s) override { exec(s); }
  virtual void visit(bin_op_t &s) override { exec(s); }
  virtual void visit(assign_t &s) override { exec(s); }
  virtual void visit(assume_t &s) override { exec(s); }
  virtual void visit(select_t &s) override { exec(s); }
  virtual void visit(assert_t &s) override { exec(s); }
  virtual void visit(int_cast_t &s) override { exec(s); }
  virtual void visit(callsite_t &s) override { exec(s); }
  virtual void visit(intrinsic_t &s) override { exec(s); }
  virtual void visit(arr_init_t &s) override { exec(s); }
  virtual void visit(arr_store_t &s) override { exec(s); }
  virtual void visit(arr_load_t &s) override { exec(s); }
  virtual void visit(arr_assign_t &s) override { exec(s); }
  virtual void visit(region_init_t &s) override { exec(s); }
  virtual void visit(region_copy_t &s) override { exec(s); }
  virtual void visit(region_cast_t &s) override { exec(s); }  
  virtual void visit(make_ref_t &s) override { exec(s); }
  virtual void visit(remove_ref_t &s) override { exec(s); }  
  virtual void visit(load_from_ref_t &s) override { exec(s); }
  virtual void visit(store_to_ref_t &s) override { exec(s); }
  virtual void visit(gep_ref_t &s) override { exec(s); }
  virtual void visit(assume_ref_t &s) override { exec(s); }
  virtual void visit(assert_ref_t &s) override { exec(s); }
  virtual void visit(select_ref_t &s) override { exec(s); }
  virtual void visit(ref_to_int_t &s) override { exec(s); }
  virtual void visit(int_to_ref_t &s) override { exec(s); }
  virtual void visit(bool_bin_op_t &s) override { exec(s); }
  virtual void visit(bool_assign_cst_t &s) override { exec(s); }
  virtual void visit(bool_assign_var_t &s) override { exec(s); }
  virtual void visit(bool_assume_t &s) override { exec(s); }
  virtual void visit(bool_select_t &s) override { exec(s); }
  virtual void visit(bool_assert_t &s) override { exec(s); }
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

/**
 * Abstract forward transformer for all statements. Function calls
 * can be redefined by derived classes. By default, all function
 * calls are ignored in a sound manner (by havoc'ing all outputs).
 **/
template <class BasicBlock, class AbsD>
class intra_abs_transformer
    : public abs_transformer_api<typename BasicBlock::basic_block_label_t,
                                 typename AbsD::number_t,
                                 typename AbsD::varname_t> {

public:
  using abs_dom_t = AbsD;
  using number_t = typename abs_dom_t::number_t;
  using varname_t = typename abs_dom_t::varname_t;
  using variable_t = typename abs_dom_t::variable_t;
  using variable_or_constant_t = variable_or_constant<number_t, varname_t>;
  using basic_block_label_t = typename BasicBlock::basic_block_label_t;

  static_assert(std::is_same<typename AbsD::number_t,
                             typename BasicBlock::number_t>::value,
                "Basic block and abstract domain must have same number type");
  static_assert(
      std::is_same<typename AbsD::varname_t,
                   typename BasicBlock::varname_t>::value,
      "Basic block and abstract domain must have same variable name type");

public:
  using abs_transform_api_t =
      abs_transformer_api<basic_block_label_t, number_t, varname_t>;

  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::ref_cst_t;
  using typename abs_transform_api_t::var_t;

  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_ref_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_ref_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::gep_ref_t;
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::int_to_ref_t;
  using typename abs_transform_api_t::intrinsic_t;
  using typename abs_transform_api_t::load_from_ref_t;
  using typename abs_transform_api_t::make_ref_t;
  using typename abs_transform_api_t::remove_ref_t;  
  using typename abs_transform_api_t::ref_to_int_t;
  using typename abs_transform_api_t::region_copy_t;
  using typename abs_transform_api_t::region_cast_t;  
  using typename abs_transform_api_t::region_init_t;
  using typename abs_transform_api_t::select_ref_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::store_to_ref_t;
  using typename abs_transform_api_t::unreach_t;

protected:
  abs_dom_t m_inv;
  bool m_ignore_assert;

private:
  template <typename NumOrVar>
  void apply(abs_dom_t &inv, crab::cfg::binary_operation_t op,
             const variable_t &x, const variable_t &y, NumOrVar z) {
    if (auto top = conv_op<domains::arith_operation_t>(op)) {
      inv.apply(*top, x, y, z);
    } else if (auto top = conv_op<domains::bitwise_operation_t>(op)) {
      inv.apply(*top, x, y, z);
    } else {
      CRAB_ERROR("unsupported binary operator", op);
    }
  }

public:
  intra_abs_transformer(abs_dom_t inv, bool ignore_assert = false)
      : m_inv(inv), m_ignore_assert(ignore_assert) {}

  virtual ~intra_abs_transformer() {}

  void set_abs_value(abs_dom_t &&inv) { m_inv = std::move(inv); }

  abs_dom_t get_abs_value() const { return m_inv; }
  abs_dom_t &get_abs_value() { return m_inv; }

  virtual void exec(bin_op_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= crab::cfg::BINOP_SDIV &&
           stmt.op() <= crab::cfg::BINOP_UREM))) {
      pre_bot = m_inv.is_bottom();
    }

    const lin_exp_t &op1 = stmt.left();
    const lin_exp_t &op2 = stmt.right();
    if (op1.get_variable() && op2.get_variable()) {
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            (*op2.get_variable()));
    } else {
      assert(op1.get_variable() && op2.is_constant());
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            op2.constant());
    }

    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= crab::cfg::BINOP_SDIV &&
           stmt.op() <= crab::cfg::BINOP_UREM))) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);
  }

  virtual void exec(select_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.select(stmt.lhs(), stmt.cond(), stmt.left(), stmt.right());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(assign_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.assign(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(assume_t &stmt) override {
    m_inv.operator+=(stmt.constraint());
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(assert_t &stmt) override {
    if (m_ignore_assert) {
      return;
    }

    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.operator+=(stmt.constraint());

    if (::crab::CrabSanityCheckFlag) {
      if (!stmt.constraint().is_contradiction()) {
        bool post_bot = m_inv.is_bottom();
        if (!(pre_bot || !post_bot)) {
          CRAB_WARN("Invariant became bottom after ", stmt, ".",
                    " This might indicate that the assertion is violated");
        }
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(int_cast_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    if (auto op = conv_op<crab::domains::int_conv_operation_t>(stmt.op())) {
      m_inv.apply(*op, stmt.dst(), stmt.src());
    } else {
      CRAB_ERROR("unsupported cast operator ", stmt.op());
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_assign_cst_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (stmt.is_rhs_linear_constraint()) {
      m_inv.assign_bool_cst(stmt.lhs(), stmt.rhs_as_linear_constraint());
    } else {
      m_inv.assign_bool_ref_cst(stmt.lhs(), stmt.rhs_as_reference_constraint());
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_assign_var_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }
    m_inv.assign_bool_var(stmt.lhs(), stmt.rhs(), stmt.is_rhs_negated());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_bin_op_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (auto op = conv_op<domains::bool_operation_t>(stmt.op())) {
      m_inv.apply_binary_bool(*op, stmt.lhs(), stmt.left(), stmt.right());
    } else {
      CRAB_WARN("Unsupported statement ", stmt);
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_assume_t &stmt) override {
    m_inv.assume_bool(stmt.cond(), stmt.is_negated());
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_select_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.select_bool(stmt.lhs(), stmt.cond(), stmt.left(), stmt.right());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(bool_assert_t &stmt) override {
    if (m_ignore_assert) {
      return;
    }

    m_inv.assume_bool(stmt.cond(), false);
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(havoc_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.operator-=(stmt.get_variable());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(unreach_t &stmt) override {
    m_inv.set_to_bottom();
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(arr_init_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_init(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                     stmt.ub_index(), stmt.val());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }
  
  virtual void exec(arr_store_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    if (stmt.lb_index().equal(stmt.ub_index())) {
      m_inv.array_store(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                        stmt.value(), stmt.is_strong_update());
    } else {
      m_inv.array_store_range(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                              stmt.ub_index(), stmt.value());
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(arr_load_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_load(stmt.lhs(), stmt.array(), stmt.elem_size(), stmt.index());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(arr_assign_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.array_assign(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(make_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_make(stmt.lhs(), stmt.region(), stmt.size(), stmt.alloc_site());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(remove_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_free(stmt.region(), stmt.ref());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }
  
  virtual void exec(region_init_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.region_init(stmt.region());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(region_copy_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.region_copy(stmt.lhs_region(), stmt.rhs_region());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(region_cast_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.region_cast(stmt.src(), stmt.dst());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }
  
  virtual void exec(load_from_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_load(stmt.ref(), stmt.region(), stmt.lhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(store_to_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_store(stmt.ref(), stmt.region(), stmt.val());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(gep_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_gep(stmt.rhs(), stmt.rhs_region(), stmt.lhs(), stmt.lhs_region(),
                  stmt.offset());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(assume_ref_t &stmt) override {
    m_inv.ref_assume(stmt.constraint());
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(assert_ref_t &stmt) override {
    if (m_ignore_assert) {
      return;
    }
    m_inv.ref_assume(stmt.constraint());
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(select_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.select_ref(stmt.lhs_ref(), stmt.lhs_rgn(),
		     stmt.cond(),
		     stmt.left_ref(), stmt.left_rgn(),
		     stmt.right_ref(), stmt.right_rgn());
    
    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(int_to_ref_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.int_to_ref(stmt.int_var(), stmt.region(), stmt.ref_var());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(ref_to_int_t &stmt) override {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.ref_to_int(stmt.region(), stmt.ref_var(), stmt.int_var());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << stmt << " :" << m_inv <<"\n";);    
  }

  virtual void exec(intrinsic_t &cs) override {
    if (cs.get_intrinsic_name() == "print_invariants") {
      // Note that we don't call the abstract transformer "intrinsic".
      // Instead, we directly print the projected invariants here.
      typename abs_dom_t::variable_vector_t vars;
      auto const&inputs = cs.get_args(); 
      for (auto in: inputs) {
	if (in.is_variable()) {
	  vars.push_back(in.get_variable());
	}
      }
      abs_dom_t copy(m_inv);
      copy.project(vars);
      crab::outs() << cs << "\n" << "\t" << copy << "\n";
    } else {
      m_inv.intrinsic(cs.get_intrinsic_name(), cs.get_args(), cs.get_lhs());
      CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << cs << " :" << m_inv <<"\n";);    
    }
  }

  virtual void exec(callsite_t &cs) override {
    for (const variable_t &vt : cs.get_lhs()) {
      m_inv.operator-=(vt); // havoc
    }
    CRAB_VERBOSE_IF(5, crab::outs() << "EXECUTED " << cs << " :" << m_inv <<"\n";);    
  }

};

///////////////////////////////////////
/// For inter-procedural analysis
///////////////////////////////////////

template <typename AbsDom> class inter_transformer_helpers {
public:
  using linear_expression_t = typename AbsDom::linear_expression_t;
  using reference_constraint_t = typename AbsDom::reference_constraint_t;
  using variable_t = typename AbsDom::variable_t;
  using number_t = typename AbsDom::number_t;

  static void unify(AbsDom &inv, const variable_t &lhs, const variable_t &rhs) {
    assert(lhs.get_type() == rhs.get_type());
    auto ty = lhs.get_type();
    if (ty.is_bool()) {
      inv.assign_bool_var(lhs, rhs, false);
    } else if (ty.is_integer() || ty.is_real()) {
      inv.assign(lhs, rhs);
    } else if (ty.is_reference()) {
      inv -= lhs;
      inv.ref_assume(reference_constraint_t::mk_eq(lhs, rhs, number_t(0)));
    } else if (ty.is_region()) {
      inv.region_copy(lhs, rhs);
    } else if (ty.is_array()) {
      inv.array_assign(lhs, rhs);
    } else {
      CRAB_ERROR("abs_transformer::unify unsupported type");
    }
  }
};

/////////////////////////////////
/// For backward analysis
/////////////////////////////////

/**
 * Abstract transformer to compute necessary preconditions.
 **/
template <class BasicBlock, class AbsD, class InvT>
class intra_necessary_preconditions_abs_transformer final
    : public abs_transformer_api<typename BasicBlock::basic_block_label_t,
                                 typename AbsD::number_t,
                                 typename AbsD::varname_t> {
public:
  using abs_dom_t = AbsD;
  using number_t = typename abs_dom_t::number_t;
  using varname_t = typename abs_dom_t::varname_t;
  using variable_t = typename abs_dom_t::variable_t;
  using basic_block_label_t = typename BasicBlock::basic_block_label_t;
  using statement_t =
      crab::cfg::statement<basic_block_label_t, number_t, varname_t>;
  using abs_transform_api_t =
      abs_transformer_api<basic_block_label_t, number_t, varname_t>;
  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_ref_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_ref_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::gep_ref_t;
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::int_to_ref_t;
  using typename abs_transform_api_t::intrinsic_t;
  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::load_from_ref_t;
  using typename abs_transform_api_t::make_ref_t;
  using typename abs_transform_api_t::remove_ref_t;  
  using typename abs_transform_api_t::ref_to_int_t;
  using typename abs_transform_api_t::region_copy_t;
  using typename abs_transform_api_t::region_cast_t;  
  using typename abs_transform_api_t::region_init_t;
  using typename abs_transform_api_t::select_ref_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::store_to_ref_t;
  using typename abs_transform_api_t::unreach_t;
  using typename abs_transform_api_t::var_t;

private:
  static_assert(std::is_same<number_t, typename BasicBlock::number_t>::value,
                "Basic block and abstract domain must have same number type");
  static_assert(
      std::is_same<varname_t, typename BasicBlock::varname_t>::value,
      "Basic block and abstract domain must have same variable name type");

  // used to compute the (necessary) preconditions
  abs_dom_t m_pre;
  // used to refine the preconditions: map from statement_t to abs_dom_t.
  InvT *m_invariants;
  // ignore assertions
  bool m_ignore_assert;
  // if m_ignore_assert is false then enable compute preconditions
  // from good states, otherwise from bad states (by negating the
  // conditions of the assert statements).
  bool m_good_states;

  abs_dom_t make_top() const { return m_pre.make_top(); }

  abs_dom_t get_forward_invariant(const statement_t *stmt) const {
    assert(m_invariants);
    assert(stmt);

    auto it = m_invariants->find(stmt);
    if (it != m_invariants->end()) {
      return it->second;
    } else {
      return make_top();
    }
  }

public:
  intra_necessary_preconditions_abs_transformer(abs_dom_t post, InvT *invars,
                                                bool good_states,
                                                bool ignore_assert = false)
      : m_pre(post), m_invariants(invars), m_ignore_assert(ignore_assert),
        m_good_states(good_states) {

    if (!m_invariants) {
      CRAB_ERROR("Invariant table cannot be null");
    }
  }

  ~intra_necessary_preconditions_abs_transformer() = default;

  abs_dom_t preconditions() { return m_pre; }

  virtual void exec(bin_op_t &stmt) override {
    auto op = conv_op<domains::arith_operation_t>(stmt.op());
    if (!op || op >= domains::OP_UDIV) {
      // ignore UDIV, SREM, UREM
      // CRAB_WARN("backward operation ", stmt.op(), " not implemented");
      m_pre -= stmt.lhs();
      return;
    }

    const lin_exp_t &op1 = stmt.left();
    const lin_exp_t &op2 = stmt.right();
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt.lhs() << " := " << op1 << " "
                                << *op << " " << op2 << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");

    if (op1.get_variable() && op2.get_variable()) {
      m_pre.backward_apply(*op, stmt.lhs(), (*op1.get_variable()),
                           (*op2.get_variable()), std::move(invariant));
    } else {
      assert(op1.get_variable() && op2.is_constant());
      m_pre.backward_apply(*op, stmt.lhs(), (*op1.get_variable()),
                           op2.constant(), std::move(invariant));
    }
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // select(x := cond ? e1: e2, post) can be reduced to
  //   pre: goto b_then;
  //   pre: goto b_else;
  //   b_then:
  //     assume(cond);
  //     x := e1;
  //     goto post;
  //   b_else:
  //     assume(not(cond));
  //     x := e2;
  //     goto post;
  //   post: ....
  virtual void exec(select_t &stmt) override {
    abs_dom_t old_pre = get_forward_invariant(&stmt);

    // -- one of the two branches is false
    abs_dom_t then_inv(old_pre);
    then_inv += stmt.cond();
    if (then_inv.is_bottom()) {
      m_pre.backward_assign(stmt.lhs(), stmt.right(), std::move(old_pre));
      m_pre += stmt.cond().negate();
      return;
    }

    abs_dom_t else_inv(old_pre);
    else_inv += stmt.cond().negate();
    if (else_inv.is_bottom()) {
      m_pre.backward_assign(stmt.lhs(), stmt.left(), std::move(old_pre));
      m_pre += stmt.cond();
      return;
    }

    // -- both branches can be possible so we join them
    abs_dom_t pre_then(m_pre);
    pre_then.backward_assign(stmt.lhs(), stmt.left(), old_pre);
    pre_then += stmt.cond();

    abs_dom_t pre_else(m_pre);
    pre_else.backward_assign(stmt.lhs(), stmt.right(), old_pre);
    pre_else += stmt.cond().negate();

    m_pre = pre_then | pre_else;
  }

  // x := e
  virtual void exec(assign_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", auto rhs = stmt.rhs();
             crab::outs() << "** " << stmt.lhs() << " := " << rhs << "\n"
                          << "\tFORWARD INV=" << invariant << "\n"
                          << "\tPOST=" << m_pre << "\n");

    m_pre.backward_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assume(c)
  // the precondition must contain c so forward and backward are the same.
  virtual void exec(assume_t &stmt) override {
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre += stmt.constraint();
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assert(c)
  virtual void exec(assert_t &stmt) override {
    if (!m_ignore_assert) {
      CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                           << "\tPOST=" << m_pre << "\n");
      if (m_good_states) {
        // similar to assume(c)
        m_pre += stmt.constraint();
      } else {
        // here we are interested in computing preconditions of the
        // error states. Thus, we propagate backwards "not c" which
        // represents the error states.
        abs_dom_t error = m_pre.make_top();
        error += stmt.constraint().negate();
        // we need to join to consider all possible preconditions to
        // error. Otherwise, if we would have two assertions
        // "assert(x >= -2); assert(x <= 2);" we would have
        // incorrectly contradictory constraints.
        m_pre |= error;
      }

      CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
    }
  }

  // similar to assume(false)
  virtual void exec(unreach_t &stmt) override {
    m_pre.set_to_bottom();
  }

  // x := *
  // x can be anything before the assignment
  virtual void exec(havoc_t &stmt) override {
    m_pre -= stmt.get_variable();
  }

  virtual void exec(int_cast_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.backward_assign(stmt.dst(), stmt.src(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  virtual void exec(bool_assign_cst_t &stmt) override {
    m_pre -= stmt.lhs();
  }
  
  virtual void exec(bool_assign_var_t &stmt) override {
    m_pre -= stmt.lhs();
  }
  
  virtual void exec(bool_bin_op_t &stmt) override {
    m_pre -= stmt.lhs();
  }
  
  virtual void exec(bool_select_t &stmt) override {
    m_pre -= stmt.lhs();
  }

  virtual void exec(bool_assume_t &stmt) override {
    // same as forward
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.assume_bool(stmt.cond(), stmt.is_negated());
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  virtual void exec(bool_assert_t &stmt) override {
    if (!m_ignore_assert) {
      CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                           << "\tPOST=" << m_pre << "\n");
      if (m_good_states) {
        // similar to bool_assume(c)
        m_pre.assume_bool(stmt.cond(), false /*non-negated*/);
      } else {
        abs_dom_t error = m_pre.make_top();
        error.assume_bool(stmt.cond(), true /*negated*/);
        m_pre |= error;
      }
      CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
    }
  }

  virtual void exec(arr_init_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_init(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                              stmt.ub_index(), stmt.val(),
                              std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  virtual void exec(arr_load_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_load(stmt.lhs(), stmt.array(), stmt.elem_size(),
                              stmt.index(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  virtual void exec(arr_store_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");

    if (stmt.lb_index().equal(stmt.ub_index())) {
      m_pre.backward_array_store(stmt.array(), stmt.elem_size(),
                                 stmt.lb_index(), stmt.value(),
                                 stmt.is_strong_update(), std::move(invariant));
    } else {
      m_pre.backward_array_store_range(stmt.array(), stmt.elem_size(),
                                       stmt.lb_index(), stmt.ub_index(),
                                       stmt.value(), std::move(invariant));
    }
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  virtual void exec(arr_assign_t &stmt) override {
    abs_dom_t invariant = get_forward_invariant(&stmt);
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // NOT IMPLEMENTED
  virtual void exec(region_init_t &stmt) override {}
  virtual void exec(region_copy_t &stmt) override {}
  virtual void exec(region_cast_t &stmt) override {}  
  virtual void exec(make_ref_t &stmt) override {}
  virtual void exec(remove_ref_t &stmt) override {}  
  virtual void exec(load_from_ref_t &stmt) override {}
  virtual void exec(store_to_ref_t &stmt) override {}
  virtual void exec(gep_ref_t &stmt) override {}
  virtual void exec(assume_ref_t &stmt) override {}
  virtual void exec(assert_ref_t &stmt) override {}
  virtual void exec(select_ref_t &stmt) override {}
  virtual void exec(int_to_ref_t &stmt) override {}
  virtual void exec(ref_to_int_t &stmt) override {}

  /// -- Call and return can be redefined by derived classes

  virtual void exec(callsite_t &cs) override {
    for (const variable_t &vt : cs.get_lhs()) {
      m_pre -= vt;
    }
  }

  virtual void exec(intrinsic_t &cs) override {
    abs_dom_t invariant = get_forward_invariant(&cs);
    m_pre.backward_intrinsic(cs.get_intrinsic_name(), cs.get_args(),
                             cs.get_lhs(), std::move(invariant));
  }
};

} // namespace analyzer
} // namespace crab
