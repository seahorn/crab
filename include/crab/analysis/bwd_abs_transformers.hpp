#pragma once

#include <crab/analysis/abs_transformer_api.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace analyzer {

/**
 * Abstract transformers to compute necessary preconditions.
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

  // current necessary preconditions
  abs_dom_t m_pre;
  // top abstract value
  abs_dom_t m_top_absval;
  // used to refine the preconditions: map from statement_t to abs_dom_t.
  InvT *m_invariants;
  // ignore assertions
  bool m_ignore_assert;
  // if m_ignore_assert is false then enable compute preconditions
  // from good states, otherwise from bad states (by negating the
  // conditions of the assert statements).
  bool m_good_states;

  abs_dom_t get_forward_invariant(const statement_t *stmt) const {
    assert(m_invariants);
    assert(stmt);

    auto it = m_invariants->find(stmt);
    if (it != m_invariants->end()) {
      return it->second;
    } else {
      return m_top_absval;
    }
  }

public:
  intra_necessary_preconditions_abs_transformer(abs_dom_t post, InvT *invars,
                                                bool good_states,
                                                bool ignore_assert = false)
    : m_pre(post), m_top_absval(post.make_top()),
      m_invariants(invars), m_ignore_assert(ignore_assert),
      m_good_states(good_states) {
    assert(m_top_absval.is_top());
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
        abs_dom_t error(m_top_absval);
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
        abs_dom_t error(m_top_absval);
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
