#pragma once

#include <crab/analysis/abs_transformer_api.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace analyzer {

/**
 * Abstract forward transformer for all statements. Function calls
 * can be redefined by derived classes. By default, all function
 * calls are ignored in a sound manner (by havocking all outputs).
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

  abs_dom_t &get_abs_value() {
    return m_inv;
  }

  const abs_dom_t &get_abs_value() const {
    return m_inv;
  }

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
} // namespace analyzer
} // namespace crab
