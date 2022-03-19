#pragma once

#include <crab/cfg/cfg.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>

namespace crab {
namespace cfg {

// forward declaration
template <typename CFG> class type_checker_visitor;

template <class CFG> class type_checker {
  CFG m_cfg;

public:
  type_checker(CFG cfg) : m_cfg(cfg) {}

  void run() {
    CRAB_LOG("type-check", crab::outs() << "Type checking CFG ...\n";);

    // some sanity checks about the CFG
    if (m_cfg.size() == 0) {
      CRAB_ERROR("CFG must have at least one basic block");
    }

    // -- LLVM does not enforce having a return instruction so a CFG
    //    might not have an exit block.
    // if (!m_cfg.has_exit())
    //   CRAB_ERROR("CFG must have exit block");
    // if (m_cfg.size() == 1) {
    //   if (!(m_cfg.exit() == m_cfg.entry()))
    //     CRAB_ERROR("CFG entry and exit must be the same");
    // }

    // check all statement are well typed
    type_checker_visitor<CFG> vis;
    for (auto &b : boost::make_iterator_range(m_cfg.begin(), m_cfg.end())) {
      b.accept(&vis);
    }

    CRAB_LOG("type-check", crab::outs() << "CFG is well-typed!\n";);
  }
}; // end class type_checker

template <class CFG>
class type_checker_visitor
    : public statement_visitor<typename CFG::basic_block_label_t,
                               typename CFG::number_t,
                               typename CFG::varname_t> {

  using B = typename CFG::basic_block_label_t;
  using V = typename CFG::varname_t;
  using N = typename CFG::number_t;
  using statement_t = typename CFG::statement_t;

  using bin_op_t = typename statement_visitor<B, N, V>::bin_op_t;
  using assign_t = typename statement_visitor<B, N, V>::assign_t;
  using assume_t = typename statement_visitor<B, N, V>::assume_t;
  using assert_t = typename statement_visitor<B, N, V>::assert_t;
  using int_cast_t = typename statement_visitor<B, N, V>::int_cast_t;
  using select_t = typename statement_visitor<B, N, V>::select_t;
  using havoc_t = typename statement_visitor<B, N, V>::havoc_t;
  using unreach_t = typename statement_visitor<B, N, V>::unreach_t;
  using callsite_t = typename statement_visitor<B, N, V>::callsite_t;
  using intrinsic_t = typename statement_visitor<B, N, V>::intrinsic_t;
  using arr_init_t = typename statement_visitor<B, N, V>::arr_init_t;
  using arr_store_t = typename statement_visitor<B, N, V>::arr_store_t;
  using arr_load_t = typename statement_visitor<B, N, V>::arr_load_t;
  using arr_assign_t = typename statement_visitor<B, N, V>::arr_assign_t;
  using make_ref_t = typename statement_visitor<B, N, V>::make_ref_t;
  using remove_ref_t = typename statement_visitor<B, N, V>::remove_ref_t;
  using region_init_t = typename statement_visitor<B, N, V>::region_init_t;
  using region_copy_t = typename statement_visitor<B, N, V>::region_copy_t;
  using region_cast_t = typename statement_visitor<B, N, V>::region_cast_t;
  using load_from_ref_t = typename statement_visitor<B, N, V>::load_from_ref_t;
  using store_to_ref_t = typename statement_visitor<B, N, V>::store_to_ref_t;
  using gep_ref_t = typename statement_visitor<B, N, V>::gep_ref_t;
  using assume_ref_t = typename statement_visitor<B, N, V>::assume_ref_t;
  using assert_ref_t = typename statement_visitor<B, N, V>::assert_ref_t;
  using select_ref_t = typename statement_visitor<B, N, V>::select_ref_t;
  using int_to_ref_t = typename statement_visitor<B, N, V>::int_to_ref_t;
  using ref_to_int_t = typename statement_visitor<B, N, V>::ref_to_int_t;
  using bool_bin_op_t = typename statement_visitor<B, N, V>::bool_bin_op_t;
  using bool_assign_cst_t =
      typename statement_visitor<B, N, V>::bool_assign_cst_t;
  using bool_assign_var_t =
      typename statement_visitor<B, N, V>::bool_assign_var_t;
  using bool_assume_t = typename statement_visitor<B, N, V>::bool_assume_t;
  using bool_assert_t = typename statement_visitor<B, N, V>::bool_assert_t;
  using bool_select_t = typename statement_visitor<B, N, V>::bool_select_t;

  using lin_exp_t = ikos::linear_expression<N, V>;
  using lin_cst_t = ikos::linear_constraint<N, V>;
  using variable_t = variable<N, V>;
  using variable_or_constant_t = variable_or_constant<N, V>;
  using variable_type_t = typename variable_t::type_t;
  using variable_bitwidth_t = typename variable_t::bitwidth_t;

  std::unordered_map<std::string, variable_t> m_varname_map;

  // check that there is no two variables with same name but
  // different types.
  void check_varname(const variable_t &v) {
    // XXX: v.name() which is of type V is not necessarily a
    // string.  We want to check the actual string name but V
    // does not have a method that returns the variable name as
    // a string. The solution is to do this hack using
    // crab_string_os
    crab::crab_string_os v1_os;
    v1_os << v;
    std::string v1_strname = v1_os.str();
    auto it = m_varname_map.find(v1_strname);
    if (it == m_varname_map.end()) {
      m_varname_map.insert({v1_strname, v});
    } else {
      crab::crab_string_os v2_os;
      v2_os << it->second.name();
      std::string v2_strname = v2_os.str();
      if (v1_strname == v2_strname && (it->second.get_type() != v.get_type())) {
        crab::crab_string_os os;
        os << "(type checking) variable name " << v.name()
           << " is used with different types ";
        it->second.dump(os);
        os << " and ";
        v.dump(os);
        CRAB_ERROR(os.str());
      }
    }
  }

  // check variable is a number
  void check_num(const variable_t &v, std::string msg, statement_t &s) {
    if (!v.get_type().is_integer() && !v.get_type().is_real()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is an integer or boolean
  void check_int_or_bool(const variable_t &v, std::string msg, statement_t &s) {
    if (!v.get_type().is_integer() && !v.get_type().is_bool()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is an integer
  void check_int(const variable_t &v, std::string msg, statement_t &s) {
    if (!v.get_type().is_integer()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is a boolean
  void check_bool(const variable_t &v, std::string msg, statement_t &s) {
    if (!v.get_type().is_bool()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is region
  void check_region(const variable_t &v, statement_t &s) {
    if (!v.get_type().is_region()) {
      crab::crab_string_os os;
      os << "(type checking) " << v << " is not a region variable in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is a reference
  void check_ref(const variable_t &v, std::string msg, statement_t &s) {
    if (!v.get_type().is_reference()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  void check_ref(const variable_or_constant_t &v, std::string msg,
                 statement_t &s) {
    if (!v.get_type().is_reference()) {
      crab::crab_string_os os;
      os << "(type checking) " << v << " is not a reference";
      if (msg != "") {
        os << ". " << msg;
      }
      os << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  void check_region_consistent_with_data(const variable_t &rgn,
                                         const variable_type &data_type,
                                         statement_t &s) {
    if (!rgn.get_type().is_region()) {
      return;
    }

    if (rgn.get_type().is_unknown_region()) {
      // anything is consistent with an unknown region
      return;
    } else if (rgn.get_type().is_bool_region() && data_type.is_bool()) {
      return;
    } else if (rgn.get_type().is_reference_region() &&
               data_type.is_reference()) {
      return;
    } else if (rgn.get_type().is_integer_region() &&
               data_type.is_integer(
                   rgn.get_type().get_integer_region_bitwidth())) {
      return;
    } else if (rgn.get_type().is_real_region() && data_type.is_real()) {
      return;
    } else {
      // TODO: other cases
      return;
    }
    crab::crab_string_os os;
    os << "(type checking) the type of " << rgn
       << " must be consistent with the type of data " << data_type << " in "
       << s;
    CRAB_ERROR(os.str());
  }

  // check two variables have same types
  void check_same_type(const variable_t &v1, const variable_t &v2,
                       std::string msg, statement_t &s) {
    if (v1.get_type() != v2.get_type()) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check two variables have different names
  void check_different_name(const variable_t &v1, const variable_t &v2,
                            std::string msg, statement_t &s) {
    if (v1 == v2) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check linear expressinon is just a number or variable
  void check_num_or_var(const lin_exp_t &e, std::string msg, statement_t &s) {
    if (!(e.is_constant() || e.get_variable())) {
      crab::crab_string_os os;
      os << "(type checking) " << msg << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // check variable is an array
  void check_is_array(const variable_t &v, statement_t &s) {
    if (!v.get_type().is_array()) {
      crab::crab_string_os os;
      os << "(type checking) " << v << " must be an array variable in " << s;
      CRAB_ERROR(os.str());
    }
  }

  // v1 is array type and v2 is a scalar type consistent with v1
  void check_array_and_scalar_type(const variable_t &v1, const variable_t &v2,
                                   statement_t &s) {
    if (v1.get_type().is_bool_array()) {
      if (v2.get_type().is_bool())
        return;
    } else if (v1.get_type().is_integer_array()) {
      if (v2.get_type().is_integer())
        return;
    } else if (v1.get_type().is_real_array()) {
      if (v2.get_type().is_real())
        return;
    } else {
      crab::crab_string_os os;
      os << "(type checking) " << v1 << " must be an array variable in " << s;
      CRAB_ERROR(os.str());
    }

    crab::crab_string_os os;
    os << "(type checking) " << v1 << " and " << v2
       << " do not have consistent types in " << s;
    CRAB_ERROR(os.str());
  }

public:
  type_checker_visitor() {}

  virtual void visit(bin_op_t &s) override {
    const variable_t &lhs = s.lhs();
    const lin_exp_t &op1 = s.left();
    const lin_exp_t &op2 = s.right();

    check_varname(lhs);
    check_num(lhs, "lhs must be integer or real", s);

    if (boost::optional<variable_t> v1 = op1.get_variable()) {
      check_varname(*v1);
      check_same_type(lhs, *v1,
                      "first operand cannot have different type from lhs", s);
    } else {
      CRAB_ERROR("(type checking) first binary operand must be a variable in ",
                 s);
    }
    if (boost::optional<variable_t> v2 = op2.get_variable()) {
      check_varname(*v2);
      check_same_type(lhs, *v2,
                      "second operand cannot have different type from lhs", s);
    } else {
      // TODO: we can still check that we use z_number
      // (q_number) of INT_TYPE (REAL_TYPE)
    }
  }

  virtual void visit(assign_t &s) override {
    const variable_t &lhs = s.lhs();
    const lin_exp_t &rhs = s.rhs();

    check_varname(lhs);
    check_num(lhs, "lhs must be integer or real", s);

    for (auto const &v : rhs.variables()) {
      check_varname(v);
      check_same_type(lhs, v, "variable cannot have different type from lhs",
                      s);
    }
  }

  virtual void visit(assume_t &s) override {
    bool first = true;
    const variable_t *first_var = nullptr;
    for (auto const &v : s.constraint().variables()) {
      check_varname(v);
      check_num(v, "assume variables must be integer or real", s);
      if (first) {
        first_var = &v;
        first = false;
      }
      assert(first_var);
      check_same_type(*first_var, v, "inconsistent types in assume variables",
                      s);
    }
  }

  virtual void visit(assert_t &s) override {
    bool first = true;
    const variable_t *first_var;
    for (auto const &v : s.constraint().variables()) {
      check_varname(v);
      check_num(v, "assert variables must be integer or real", s);
      if (first) {
        first_var = &v;
        first = false;
      }
      assert(first_var);
      check_same_type(*first_var, v, "inconsistent types in assert variables",
                      s);
    }
  }

  virtual void visit(select_t &s) override {
    check_num(s.lhs(), "lhs must be integer or real", s);
    check_varname(s.lhs());
    for (auto const &v : s.left().variables()) {
      check_varname(v);
      check_same_type(s.lhs(), v, "inconsistent types in select variables", s);
    }
    for (auto const &v : s.right().variables()) {
      check_varname(v);
      check_same_type(s.lhs(), v, "inconsistent types in select variables", s);
    }

    // -- The condition can have different type from lhs/left/right
    //    operands.
    bool first = true;
    const variable_t *first_var;
    for (auto const &v : s.cond().variables()) {
      check_varname(v);
      check_num(v, "assume variables must be integer or real", s);
      if (first) {
        first_var = &v;
        first = false;
      }
      assert(first_var);
      check_same_type(*first_var, v,
                      "inconsistent types in select condition variables", s);
    }
  }

  virtual void visit(int_cast_t &s) override {
    const variable_t &src = s.src();
    const variable_t &dst = s.dst();

    auto get_bitwidth = [](const variable_t &v) {
      return (v.get_type().is_integer() ? v.get_type().get_integer_bitwidth()
                                        : 1);
    };

    check_varname(src);
    check_varname(dst);
    switch (s.op()) {
    case CAST_TRUNC:
      check_int(src, "source operand must be integer", s);
      check_int_or_bool(dst, "destination must be integer or bool", s);
      if (get_bitwidth(src) <= get_bitwidth(dst)) {
        CRAB_ERROR("(type checking) bitwidth of source operand must be "
                   "greater than destination in ",
                   s);
      }
      break;
    default:
      // case CAST_SEXT:
      // case CAST_ZEXT:
      check_int(dst, "destination operand must be integer", s);
      check_int_or_bool(src, "source must be integer or bool", s);
      if (get_bitwidth(dst) <= get_bitwidth(src)) {
        CRAB_ERROR("(type checking) bitwidth of destination must be greater "
                   "than source in ",
                   s);
      }
    }
  }

  virtual void visit(havoc_t &) override {}

  virtual void visit(unreach_t &) override {}

  virtual void visit(bool_bin_op_t &s) override {
    check_varname(s.lhs());
    check_varname(s.left());
    check_varname(s.right());
    check_bool(s.lhs(), "lhs must be boolean", s);
    check_bool(s.left(), "first operand must be boolean", s);
    check_bool(s.right(), "second operand must be boolean", s);
  };

  virtual void visit(bool_assign_cst_t &s) override {
    check_bool(s.lhs(), "lhs must be boolean", s);
    check_varname(s.lhs());
    bool first = true;
    const variable_t *first_var;

    if (s.is_rhs_linear_constraint()) {
      for (auto const &v : s.rhs_as_linear_constraint().variables()) {
        check_varname(v);
        check_num(v, "rhs variables must be integer or real", s);
        if (first) {
          first_var = &v;
          first = false;
        }
        assert(first_var);
        check_same_type(*first_var, v, "inconsistent types in rhs variables",
                        s);
      }
    } else {
      for (auto const &v : s.rhs_as_reference_constraint().variables()) {
        check_varname(v);
        check_ref(v, "rhs variables must be reference", s);
        if (first) {
          first_var = &v;
          first = false;
        }
        assert(first_var);
        check_same_type(*first_var, v, "inconsistent types in rhs variables",
                        s);
      }
    }
  };

  virtual void visit(bool_assign_var_t &s) override {
    check_bool(s.lhs(), "lhs must be boolean", s);
    check_bool(s.rhs(), "rhs must be boolean", s);
    check_varname(s.lhs());
    check_varname(s.rhs());
  };

  virtual void visit(bool_assume_t &s) override {
    check_bool(s.cond(), "condition must be boolean", s);
    check_varname(s.cond());
  };

  virtual void visit(bool_assert_t &s) override {
    check_bool(s.cond(), "condition must be boolean", s);
    check_varname(s.cond());
  };

  virtual void visit(bool_select_t &s) override {
    check_bool(s.lhs(), "lhs must be boolean", s);
    check_bool(s.lhs(), "condition must be boolean", s);
    check_bool(s.left(), "first operand must be boolean", s);
    check_bool(s.right(), "second operand must be boolean", s);
    check_varname(s.lhs());
    check_varname(s.left());
    check_varname(s.right());
  };

  virtual void visit(arr_init_t &s) override {
    // TODO: check that e_sz is the same number that v's bitwidth
    const variable_t &a = s.array();
    const lin_exp_t &e_sz = s.elem_size();
    const lin_exp_t &lb = s.lb_index();
    const lin_exp_t &ub = s.ub_index();
    const lin_exp_t &v = s.val();

    check_is_array(a, s);
    check_varname(a);
    check_num_or_var(lb, "array lower bound index must be number or variable",
                     s);
    check_num_or_var(ub, "array upper bound index must be number or variable",
                     s);
    check_num_or_var(e_sz, "array element size must be number or variable", s);
    check_num_or_var(v, "array value must be number or variable", s);
    // TODO: check_varname lb and ub if variables
    if (boost::optional<variable_t> vv = v.get_variable()) {
      check_varname(*vv);
      check_array_and_scalar_type(a, *vv, s);
    }
  }

  virtual void visit(arr_store_t &s) override {
    // TODO: check that e_sz is the same number that v's bitwidth
    const variable_t &a = s.array();
    check_is_array(a, s);
    check_varname(a);

    const lin_exp_t &e_sz = s.elem_size();
    const lin_exp_t &v = s.value();
    if (s.is_strong_update()) {
      if (!(s.lb_index().equal(s.ub_index()))) {
        crab::crab_string_os os;
        os << "(type checking) "
           << "array lower and upper indexes must be equal because "
              "strong_update array in "
           << s;
        CRAB_ERROR(os.str());
      }
    }
    for (auto const &iv : s.lb_index().variables()) {
      check_varname(iv);
      check_int_or_bool(
          iv, "array index must contain only integer or boolean variables", s);
    }
    for (auto const &iv : s.lb_index().variables()) {
      check_varname(iv);
      check_int_or_bool(
          iv, "array index must contain only integer or boolean variables", s);
    }
    check_num_or_var(e_sz, "array element size must be number or variable", s);
    check_num_or_var(v, "array value must be number or variable", s);
    if (boost::optional<variable_t> vv = v.get_variable()) {
      check_varname(*vv);
      check_array_and_scalar_type(a, *vv, s);
    }
  }

  virtual void visit(arr_load_t &s) override {
    // TODO: check that e_sz is the same number that lhs's bitwidth
    const variable_t &a = s.array();
    const lin_exp_t &e_sz = s.elem_size();
    const variable_t &lhs = s.lhs();
    check_varname(lhs);
    check_is_array(a, s);
    check_varname(a);
    for (auto const &iv : s.index().variables()) {
      check_varname(iv);
      check_int_or_bool(
          iv, "array index must contain only integer or boolean variables", s);
    }
    check_num_or_var(e_sz, "array element size must be number or variable", s);
    check_array_and_scalar_type(a, lhs, s);
  }

  virtual void visit(arr_assign_t &s) override {
    const variable_t &lhs = s.lhs();
    const variable_t &rhs = s.rhs();
    check_is_array(lhs, s);
    check_is_array(rhs, s);
    check_varname(lhs);
    check_varname(rhs);
    check_same_type(lhs, rhs, "array assign variables must have same type", s);
  }

  virtual void visit(callsite_t &s) override {
    // The type consistency with the callee parameters is done
    // elsewhere.
    for (const variable_t &v : s.get_lhs()) {
      check_varname(v);
    }
    for (const variable_t &v : s.get_args()) {
      check_varname(v);
    }
  }

  virtual void visit(intrinsic_t &s) override {
    for (const variable_t &v : s.get_lhs()) {
      check_varname(v);
    }
    for (const variable_or_constant_t &vc : s.get_args()) {
      if (vc.is_variable()) {
        check_varname(vc.get_variable());
      }
    }
  }

  virtual void visit(region_init_t &s) override { check_region(s.region(), s); }

  virtual void visit(region_copy_t &s) override {
    check_region(s.lhs_region(), s);
    check_region(s.rhs_region(), s);
    check_same_type(s.lhs_region(), s.rhs_region(),
                    "region_copy must have same types", s);
  }

  virtual void visit(region_cast_t &s) override {
    check_region(s.src(), s);
    check_region(s.dst(), s);
    if (!(s.src().get_type().is_unknown_region() ^
          s.dst().get_type().is_unknown_region())) {
      crab::crab_string_os os;
      os << "(type checking) "
         << "one operand must be unknown region and the other not "
         << " in " << s;
      CRAB_ERROR(os.str());
    }
  }

  virtual void visit(make_ref_t &s) override {
    check_region(s.region(), s);
    check_ref(s.lhs(), "", s);
  }

  virtual void visit(remove_ref_t &s) override {
    check_region(s.region(), s);
    check_ref(s.ref(), "", s);
  }

  virtual void visit(load_from_ref_t &s) override {
    check_region(s.region(), s);
    check_ref(s.ref(), "", s);
    check_region_consistent_with_data(s.region(), s.lhs().get_type(), s);
  }

  virtual void visit(store_to_ref_t &s) override {
    check_region(s.region(), s);
    check_ref(s.ref(), "", s);
    check_region_consistent_with_data(s.region(), s.val().get_type(), s);
  }

  virtual void visit(gep_ref_t &s) override {
    check_region(s.lhs_region(), s);
    check_region(s.rhs_region(), s);
    check_ref(s.lhs(), "", s);
    check_ref(s.rhs(), "", s);
  }

  virtual void visit(assume_ref_t &s) override {}
  virtual void visit(assert_ref_t &s) override {}
  virtual void visit(select_ref_t &s) override {
    // TODO: check region operands
    check_ref(s.lhs_ref(), "lhs must be reference", s);
    check_bool(s.cond(), "condition must be boolean", s);
    check_ref(s.left_ref(), "first operand must be reference", s);
    check_ref(s.right_ref(), "second operand must be reference", s);
    check_varname(s.lhs_ref());
    if (s.left_ref().is_variable())
      check_varname(s.left_ref().get_variable());
    if (s.right_ref().is_variable())
      check_varname(s.right_ref().get_variable());
  }
  virtual void visit(int_to_ref_t &s) override {
    check_region(s.region(), s);
    check_ref(s.ref_var(), "", s);
    check_num(s.int_var(), "first input must be a number", s);
  }
  virtual void visit(ref_to_int_t &s) override {
    check_region(s.region(), s);
    check_ref(s.ref_var(), "", s);
    check_num(s.int_var(), "first input must be a number", s);
  }

}; // end class type_checker_visitor

} // end namespace cfg
} // end namespace crab
