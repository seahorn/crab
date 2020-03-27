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
     a[l...u] := v (a,b are arrays and v can be bool/integer/pointer)
     a[i] := v;
     v := a[i];
     a := b

   POINTERS
     *p = q;
     p = *q;
     p := q+n
     p := &obj;
     p := &fun
     p := null;

   FUNCTIONS
     x := foo(arg1,...,argn);
     return r;

   havoc(x);

 */

#include <crab/cfg/cfg.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/abstract_domain_operators.hpp>
#include <crab/domains/linear_constraints.hpp>

namespace crab {
namespace analyzer {

/**
 * API abstract transformer
 **/
template <typename Number, typename VariableName>
class abs_transformer_api
    : public crab::cfg::statement_visitor<Number, VariableName> {
public:
  typedef Number number_t;
  typedef VariableName varname_t;

  typedef ikos::variable<number_t, VariableName> var_t;
  typedef ikos::linear_expression<number_t, VariableName> lin_exp_t;
  typedef ikos::linear_constraint<number_t, VariableName> lin_cst_t;
  typedef ikos::linear_constraint_system<number_t, VariableName> lin_cst_sys_t;

  typedef crab::cfg::havoc_stmt<number_t, VariableName> havoc_t;
  typedef crab::cfg::unreachable_stmt<number_t, VariableName> unreach_t;

  typedef crab::cfg::binary_op<number_t, VariableName> bin_op_t;
  typedef crab::cfg::assignment<number_t, VariableName> assign_t;
  typedef crab::cfg::assume_stmt<number_t, VariableName> assume_t;
  typedef crab::cfg::select_stmt<number_t, VariableName> select_t;
  typedef crab::cfg::assert_stmt<number_t, VariableName> assert_t;
  typedef crab::cfg::int_cast_stmt<number_t, VariableName> int_cast_t;
  
  typedef crab::cfg::callsite_stmt<number_t, VariableName> callsite_t;
  typedef crab::cfg::return_stmt<number_t, VariableName> return_t;
  typedef crab::cfg::intrinsic_stmt<number_t, VariableName> intrinsic_t;  

  typedef crab::cfg::array_init_stmt<number_t, VariableName> arr_init_t;
  typedef crab::cfg::array_store_stmt<number_t, VariableName> arr_store_t;
  typedef crab::cfg::array_load_stmt<number_t, VariableName> arr_load_t;
  typedef crab::cfg::array_assign_stmt<number_t, VariableName> arr_assign_t;
  typedef crab::cfg::ptr_store_stmt<number_t, VariableName> ptr_store_t;
  typedef crab::cfg::ptr_load_stmt<number_t, VariableName> ptr_load_t;

  typedef crab::cfg::ptr_assign_stmt<number_t, VariableName> ptr_assign_t;
  typedef crab::cfg::ptr_object_stmt<number_t, VariableName> ptr_object_t;
  typedef crab::cfg::ptr_function_stmt<number_t, VariableName> ptr_function_t;
  typedef crab::cfg::ptr_null_stmt<number_t, VariableName> ptr_null_t;
  typedef crab::cfg::ptr_assume_stmt<number_t, VariableName> ptr_assume_t;
  typedef crab::cfg::ptr_assert_stmt<number_t, VariableName> ptr_assert_t;

  typedef crab::cfg::bool_binary_op<number_t, VariableName> bool_bin_op_t;
  typedef crab::cfg::bool_assign_cst<number_t, VariableName> bool_assign_cst_t;
  typedef crab::cfg::bool_assign_var<number_t, VariableName> bool_assign_var_t;
  typedef crab::cfg::bool_assume_stmt<number_t, VariableName> bool_assume_t;
  typedef crab::cfg::bool_select_stmt<number_t, VariableName> bool_select_t;
  typedef crab::cfg::bool_assert_stmt<number_t, VariableName> bool_assert_t;

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
  virtual void exec(return_t &) {}
  virtual void exec(intrinsic_t &) {}  
  virtual void exec(arr_init_t &) {}
  virtual void exec(arr_store_t &) {}
  virtual void exec(arr_load_t &) {}
  virtual void exec(arr_assign_t &) {}
  virtual void exec(ptr_store_t &) {}
  virtual void exec(ptr_load_t &) {}
  virtual void exec(ptr_assign_t &) {}
  virtual void exec(ptr_object_t &) {}
  virtual void exec(ptr_function_t &) {}
  virtual void exec(ptr_null_t &) {}
  virtual void exec(ptr_assume_t &) {}
  virtual void exec(ptr_assert_t &) {}
  virtual void exec(bool_bin_op_t &) {}
  virtual void exec(bool_assign_cst_t &) {}
  virtual void exec(bool_assign_var_t &) {}
  virtual void exec(bool_assume_t &) {}
  virtual void exec(bool_select_t &) {}
  virtual void exec(bool_assert_t &) {}

public: /* visitor api */
  void visit(havoc_t &s) { exec(s); }
  void visit(unreach_t &s) { exec(s); }
  void visit(bin_op_t &s) { exec(s); }
  void visit(assign_t &s) { exec(s); }
  void visit(assume_t &s) { exec(s); }
  void visit(select_t &s) { exec(s); }
  void visit(assert_t &s) { exec(s); }
  void visit(int_cast_t &s) { exec(s); }
  void visit(callsite_t &s) { exec(s); }
  void visit(return_t &s) { exec(s); }
  void visit(intrinsic_t &s) { exec(s); }  
  void visit(arr_init_t &s) { exec(s); }
  void visit(arr_store_t &s) { exec(s); }
  void visit(arr_load_t &s) { exec(s); }
  void visit(arr_assign_t &s) { exec(s); }
  void visit(ptr_store_t &s) { exec(s); }
  void visit(ptr_load_t &s) { exec(s); }
  void visit(ptr_assign_t &s) { exec(s); }
  void visit(ptr_object_t &s) { exec(s); }
  void visit(ptr_function_t &s) { exec(s); }
  void visit(ptr_null_t &s) { exec(s); }
  void visit(ptr_assume_t &s) { exec(s); }
  void visit(ptr_assert_t &s) { exec(s); }
  void visit(bool_bin_op_t &s) { exec(s); }
  void visit(bool_assign_cst_t &s) { exec(s); }
  void visit(bool_assign_var_t &s) { exec(s); }
  void visit(bool_assume_t &s) { exec(s); }
  void visit(bool_select_t &s) { exec(s); }
  void visit(bool_assert_t &s) { exec(s); }
};

/**
 * Abstract forward transformer for all statements. Function calls
 * can be redefined by derived classes. By default, all function
 * calls are ignored in a sound manner (by havoc'ing all outputs).
 **/
template <class AbsD>
class intra_abs_transformer
    : public abs_transformer_api<typename AbsD::number_t,
                                 typename AbsD::varname_t> {

public:
  typedef AbsD abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::variable_t variable_t;

public:
  typedef abs_transformer_api<number_t, varname_t> abs_transform_api_t;
  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::intrinsic_t;  
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::ptr_assert_t;
  using typename abs_transform_api_t::ptr_assign_t;
  using typename abs_transform_api_t::ptr_assume_t;
  using typename abs_transform_api_t::ptr_function_t;
  using typename abs_transform_api_t::ptr_load_t;
  using typename abs_transform_api_t::ptr_null_t;
  using typename abs_transform_api_t::ptr_object_t;
  using typename abs_transform_api_t::ptr_store_t;
  using typename abs_transform_api_t::return_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::unreach_t;
  using typename abs_transform_api_t::var_t;

protected:
  abs_dom_t m_inv;
  bool m_ignore_assert;

private:
  template <typename NumOrVar>
  void apply(abs_dom_t &inv, binary_operation_t op, variable_t x, variable_t y,
             NumOrVar z) {
    if (auto top = conv_op<ikos::operation_t>(op)) {
      inv.apply(*top, x, y, z);
    } else if (auto top = conv_op<ikos::bitwise_operation_t>(op)) {
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

  abs_dom_t &get_abs_value() { return m_inv; }

  void exec(bin_op_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= BINOP_SDIV && stmt.op() <= BINOP_UREM))) {
      pre_bot = m_inv.is_bottom();
    }

    auto op1 = stmt.left();
    auto op2 = stmt.right();
    if (op1.get_variable() && op2.get_variable()) {
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            (*op2.get_variable()));
    } else {
      assert(op1.get_variable() && op2.is_constant());
      apply(m_inv, stmt.op(), stmt.lhs(), (*op1.get_variable()),
            op2.constant());
    }

    if (::crab::CrabSanityCheckFlag &&
        (!(stmt.op() >= BINOP_SDIV && stmt.op() <= BINOP_UREM))) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(select_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    abs_dom_t inv1(m_inv);
    abs_dom_t inv2(m_inv);

    inv1 += stmt.cond();
    inv2 += stmt.cond().negate();

    if (::crab::CrabSanityCheckFlag) {
      if (!pre_bot && (inv1.is_bottom() && inv2.is_bottom())) {
        CRAB_ERROR(
            "select condition and its negation cannot be false simultaneously ",
            stmt);
      }
    }

    if (inv2.is_bottom()) {
      inv1.assign(stmt.lhs(), stmt.left());
      m_inv = inv1;
    } else if (inv1.is_bottom()) {
      inv2.assign(stmt.lhs(), stmt.right());
      m_inv = inv2;
    } else {
      inv1.assign(stmt.lhs(), stmt.left());
      inv2.assign(stmt.lhs(), stmt.right());
      m_inv = inv1 | inv2;
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(assign_t &stmt) {
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
  }

  void exec(assume_t &stmt) { m_inv.operator+=(stmt.constraint()); }

  void exec(assert_t &stmt) {
    if (m_ignore_assert)
      return;

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
  }

  void exec(int_cast_t &stmt) {
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
  }

  void exec(bool_assign_cst_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.assign_bool_cst(stmt.lhs(), stmt.rhs());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assign_var_t &stmt) {
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
  }

  void exec(bool_bin_op_t &stmt) {
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
  }

  void exec(bool_assume_t &stmt) {
    m_inv.assume_bool(stmt.cond(), stmt.is_negated());
  }

  void exec(bool_select_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    abs_dom_t inv1(m_inv);
    abs_dom_t inv2(m_inv);
    const bool negate = true;
    inv1.assume_bool(stmt.cond(), !negate);
    inv2.assume_bool(stmt.cond(), negate);
    if (inv2.is_bottom()) {
      inv1.assign_bool_var(stmt.lhs(), stmt.left(), !negate);
      m_inv = inv1;
    } else if (inv1.is_bottom()) {
      inv2.assign_bool_var(stmt.lhs(), stmt.right(), !negate);
      m_inv = inv2;
    } else {
      inv1.assign_bool_var(stmt.lhs(), stmt.left(), !negate);
      inv2.assign_bool_var(stmt.lhs(), stmt.right(), !negate);
      m_inv = inv1 | inv2;
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(bool_assert_t &stmt) {
    if (m_ignore_assert)
      return;

    m_inv.assume_bool(stmt.cond(), false);
  }

  void exec(havoc_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    m_inv.operator-=(stmt.variable());

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(unreach_t &stmt) { m_inv.set_to_bottom(); }

  void exec(arr_init_t &stmt) {
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
  }

  void exec(arr_store_t &stmt) {
    bool pre_bot = false;
    if (::crab::CrabSanityCheckFlag) {
      pre_bot = m_inv.is_bottom();
    }

    auto new_arr_v = stmt.new_array();
    if (stmt.lb_index().equal(stmt.ub_index())) {
      if (new_arr_v) {
        m_inv.array_store(*new_arr_v, stmt.array(), stmt.elem_size(),
                          stmt.lb_index(), stmt.value(),
                          stmt.is_strong_update());
      } else {
        m_inv.array_store(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                          stmt.value(), stmt.is_strong_update());
      }
    } else {
      if (new_arr_v) {
        m_inv.array_store_range(*new_arr_v, stmt.array(), stmt.elem_size(),
                                stmt.lb_index(), stmt.ub_index(), stmt.value());
      } else {
        m_inv.array_store_range(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                                stmt.ub_index(), stmt.value());
      }
    }

    if (::crab::CrabSanityCheckFlag) {
      bool post_bot = m_inv.is_bottom();
      if (!(pre_bot || !post_bot)) {
        CRAB_ERROR("Invariant became bottom after ", stmt);
      }
    }
  }

  void exec(arr_load_t &stmt) {
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
  }

  void exec(arr_assign_t &stmt) {
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
  }

  void exec(ptr_null_t &stmt) { m_inv.pointer_mk_null(stmt.lhs()); }

  void exec(ptr_object_t &stmt) {
    m_inv.pointer_mk_obj(stmt.lhs(), stmt.rhs());
  }

  void exec(ptr_assign_t &stmt) {
    m_inv.pointer_assign(stmt.lhs(), stmt.rhs(), stmt.offset());
  }

  void exec(ptr_function_t &stmt) {
    m_inv.pointer_function(stmt.lhs(), stmt.rhs());
  }

  void exec(ptr_load_t &stmt) { m_inv.pointer_load(stmt.lhs(), stmt.rhs()); }

  void exec(ptr_store_t &stmt) { m_inv.pointer_store(stmt.lhs(), stmt.rhs()); }

  void exec(ptr_assume_t &stmt) { m_inv.pointer_assume(stmt.constraint()); }

  void exec(ptr_assert_t &stmt) {
    if (m_ignore_assert)
      return;
    m_inv.pointer_assert(stmt.constraint());
  }

  void exec(intrinsic_t &cs) {
    m_inv.intrinsic(cs.get_intrinsic_name(), cs.get_args(), cs.get_lhs());
  }
  
  virtual void exec(callsite_t &cs) {
    for (auto vt : cs.get_lhs()) {
      m_inv.operator-=(vt); // havoc
    }
  }
  
  virtual void exec(return_t &ret) {}
  
};

///////////////////////////////////////
/// For inter-procedural analysis
///////////////////////////////////////

template <typename AbsDom> class inter_transformer_helpers {
public:
  typedef typename AbsDom::linear_expression_t linear_expression_t;
  typedef typename AbsDom::variable_t variable_t;
  typedef typename AbsDom::number_t number_t;

  static void unify(AbsDom &inv, variable_t lhs, variable_t rhs) {
    assert(lhs.get_type() == rhs.get_type());
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      inv.assign_bool_var(lhs, rhs, false);
      break;
    case INT_TYPE:
    case REAL_TYPE:
      inv.assign(lhs, rhs);
      break;
    case PTR_TYPE:
      inv.pointer_assign(lhs, rhs, number_t(0));
      break;
    case ARR_BOOL_TYPE:
    case ARR_INT_TYPE:
    case ARR_REAL_TYPE:
    case ARR_PTR_TYPE:
      inv.array_assign(lhs, rhs);
      break;
    default:
      CRAB_ERROR("unsuported type");
    }
  }
};

/////////////////////////////////
/// For backward analysis
/////////////////////////////////

/**
 * Abstract transformer to compute necessary preconditions.
 **/
template <class AbsD, class InvT>
class intra_necessary_preconditions_abs_transformer final
    : public abs_transformer_api<typename AbsD::number_t,
                                 typename AbsD::varname_t> {
public:
  typedef AbsD abs_dom_t;
  typedef typename abs_dom_t::number_t number_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::variable_t variable_t;
  typedef crab::cfg::statement<number_t, varname_t> statement_t;
  typedef abs_transformer_api<number_t, varname_t> abs_transform_api_t;
  using typename abs_transform_api_t::arr_assign_t;
  using typename abs_transform_api_t::arr_init_t;
  using typename abs_transform_api_t::arr_load_t;
  using typename abs_transform_api_t::arr_store_t;
  using typename abs_transform_api_t::assert_t;
  using typename abs_transform_api_t::assign_t;
  using typename abs_transform_api_t::assume_t;
  using typename abs_transform_api_t::bin_op_t;
  using typename abs_transform_api_t::bool_assert_t;
  using typename abs_transform_api_t::bool_assign_cst_t;
  using typename abs_transform_api_t::bool_assign_var_t;
  using typename abs_transform_api_t::bool_assume_t;
  using typename abs_transform_api_t::bool_bin_op_t;
  using typename abs_transform_api_t::bool_select_t;
  using typename abs_transform_api_t::callsite_t;
  using typename abs_transform_api_t::intrinsic_t;  
  using typename abs_transform_api_t::havoc_t;
  using typename abs_transform_api_t::int_cast_t;
  using typename abs_transform_api_t::lin_cst_sys_t;
  using typename abs_transform_api_t::lin_cst_t;
  using typename abs_transform_api_t::lin_exp_t;
  using typename abs_transform_api_t::ptr_assert_t;
  using typename abs_transform_api_t::ptr_assign_t;
  using typename abs_transform_api_t::ptr_assume_t;
  using typename abs_transform_api_t::ptr_function_t;
  using typename abs_transform_api_t::ptr_load_t;
  using typename abs_transform_api_t::ptr_null_t;
  using typename abs_transform_api_t::ptr_object_t;
  using typename abs_transform_api_t::ptr_store_t;
  using typename abs_transform_api_t::return_t;
  using typename abs_transform_api_t::select_t;
  using typename abs_transform_api_t::unreach_t;
  using typename abs_transform_api_t::var_t;

private:
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

  void exec(bin_op_t &stmt) {
    auto op = conv_op<ikos::operation_t>(stmt.op());
    if (!op || op >= ikos::OP_UDIV) {
      // ignore UDIV, SREM, UREM
      // CRAB_WARN("backward operation ", stmt.op(), " not implemented");
      m_pre -= stmt.lhs();
      return;
    }

    auto op1 = stmt.left();
    auto op2 = stmt.right();
    abs_dom_t invariant = (*m_invariants)[&stmt];

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
  void exec(select_t &stmt) {
    abs_dom_t old_pre = (*m_invariants)[&stmt];

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
  void exec(assign_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];

    CRAB_LOG("backward-tr", auto rhs = stmt.rhs();
             crab::outs() << "** " << stmt.lhs() << " := " << rhs << "\n"
                          << "\tFORWARD INV=" << invariant << "\n"
                          << "\tPOST=" << m_pre << "\n");

    m_pre.backward_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assume(c)
  // the precondition must contain c so forward and backward are the same.
  void exec(assume_t &stmt) {
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre += stmt.constraint();
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // assert(c)
  void exec(assert_t &stmt) {
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
        abs_dom_t error;
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
  void exec(unreach_t &stmt) { m_pre.set_to_bottom(); }

  // x := *
  // x can be anything before the assignment
  void exec(havoc_t &stmt) { m_pre -= stmt.variable(); }

  void exec(int_cast_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.backward_assign(stmt.dst(), stmt.src(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(bool_assign_cst_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_assign_var_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_bin_op_t &stmt) { m_pre -= stmt.lhs(); }
  void exec(bool_select_t &stmt) { m_pre -= stmt.lhs(); }

  void exec(bool_assume_t &stmt) {
    // same as forward
    CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                         << "\tPOST=" << m_pre << "\n");
    m_pre.assume_bool(stmt.cond(), stmt.is_negated());
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(bool_assert_t &stmt) {
    if (!m_ignore_assert) {
      CRAB_LOG("backward-tr", crab::outs() << "** " << stmt << "\n"
                                           << "\tPOST=" << m_pre << "\n");
      if (m_good_states) {
        // similar to bool_assume(c)
        m_pre.assume_bool(stmt.cond(), false /*non-negated*/);
      } else {
        abs_dom_t error;
        error.assume_bool(stmt.cond(), true /*negated*/);
        m_pre |= error;
      }
      CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
    }
  }

  void exec(arr_init_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_init(stmt.array(), stmt.elem_size(), stmt.lb_index(),
                              stmt.ub_index(), stmt.val(),
                              std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_load_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];

    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_load(stmt.lhs(), stmt.array(), stmt.elem_size(),
                              stmt.index(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_store_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");

    auto new_arr_v = stmt.new_array();
    if (stmt.lb_index().equal(stmt.ub_index())) {
      if (new_arr_v) {
        m_pre.backward_array_store(
            *new_arr_v, stmt.array(), stmt.elem_size(), stmt.lb_index(),
            stmt.value(), stmt.is_strong_update(), std::move(invariant));
      } else {
        m_pre.backward_array_store(
            stmt.array(), stmt.elem_size(), stmt.lb_index(), stmt.value(),
            stmt.is_strong_update(), std::move(invariant));
      }
    } else {
      if (new_arr_v) {
        m_pre.backward_array_store_range(
            *new_arr_v, stmt.array(), stmt.elem_size(), stmt.lb_index(),
            stmt.ub_index(), stmt.value(), std::move(invariant));
      } else {
        m_pre.backward_array_store_range(stmt.array(), stmt.elem_size(),
                                         stmt.lb_index(), stmt.ub_index(),
                                         stmt.value(), std::move(invariant));
      }
    }
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  void exec(arr_assign_t &stmt) {
    abs_dom_t invariant = (*m_invariants)[&stmt];
    CRAB_LOG("backward-tr", crab::outs()
                                << "** " << stmt << "\n"
                                << "\tFORWARD INV=" << invariant << "\n"
                                << "\tPOST=" << m_pre << "\n");
    m_pre.backward_array_assign(stmt.lhs(), stmt.rhs(), std::move(invariant));
    CRAB_LOG("backward-tr", crab::outs() << "\tPRE=" << m_pre << "\n");
  }

  // NOT IMPLEMENTED
  void exec(ptr_null_t &stmt) {}
  void exec(ptr_object_t &stmt) {}
  void exec(ptr_assign_t &stmt) {}
  void exec(ptr_function_t &stmt) {}
  void exec(ptr_load_t &stmt) {}
  void exec(ptr_store_t &stmt) {}
  void exec(ptr_assume_t &stmt) {}
  void exec(ptr_assert_t &stmt) {}

  /// -- Call and return can be redefined by derived classes

  virtual void exec(callsite_t &cs) {
    for (auto vt : cs.get_lhs()) {
      m_pre -= vt;
    }
  }
  virtual void exec(return_t &stmt) {}

  void exec(intrinsic_t &cs) {
    abs_dom_t invariant = (*m_invariants)[&cs];    
    m_pre.backward_intrinsic(cs.get_intrinsic_name(), cs.get_args(), cs.get_lhs(),
			     std::move(invariant));    
  }
  
};

} // namespace analyzer
} // namespace crab
