#pragma once

/*
   Property checker for division by zero
 */

#include <crab/checkers/base_property.hpp>
#include <crab/common/types.hpp>

namespace crab {

namespace checker {

template <typename Analyzer>
class div_zero_property_checker : public property_checker<Analyzer> {

  typedef typename Analyzer::varname_t varname_t;
  typedef typename Analyzer::number_t number_t;

  typedef ikos::interval<number_t> interval_t;
  typedef typename Analyzer::abs_dom_t abs_dom_t;
  typedef property_checker<Analyzer> base_checker_t;
  using typename base_checker_t::bin_op_t;
  using typename base_checker_t::lin_cst_sys_t;
  using typename base_checker_t::lin_cst_t;
  using typename base_checker_t::lin_exp_t;
  using typename base_checker_t::var_t;

public:
  using analyzer_t = Analyzer;

  div_zero_property_checker(int verbose = 0) : base_checker_t(verbose) {}

  std::string get_property_name() const override {
    return "integer division by zero checker";
  }

  void check(bin_op_t &s) override {
    if (!this->m_abs_tr)
      return;

    if (s.op() == BINOP_SDIV || s.op() == BINOP_UDIV || s.op() == BINOP_SREM ||
        s.op() == BINOP_UREM) {

      auto &inv = this->m_abs_tr->get_abs_value();
      if (inv.is_bottom()) {
        this->m_db.add(_UNREACH);
        return;
      }

      auto divisor_expr = s.right();
      if (divisor_expr.is_constant()) {
        number_t divisor = divisor_expr.constant();
        if (divisor == number_t(0)) {
          this->m_db.add(_ERR);
        } else {
          this->m_db.add(_SAFE);
        }
      } else if (auto var = divisor_expr.get_variable()) {
        interval_t divisor_intv = inv[(*var)];
        if (auto divisor = divisor_intv.singleton()) {
          if (*divisor == number_t(0)) {
            lin_cst_t e_cst(*var != number_t(0));
            crab::crab_string_os os;
            if (this->m_verbose >= 1) {
              os << "Property : " << e_cst << "\n";
              os << "Invariant: " << inv;
            }
            this->add_error(os.str(), &s);
          } else {
            this->m_db.add(_SAFE);
          }
        } else if (interval_t(number_t(0)) <= divisor_intv) {
          lin_cst_t w_cst(*var != number_t(0));
          crab::crab_string_os os;
          if (this->m_verbose >= 2) {
            os << "Property : " << w_cst << "\n";
            os << "Invariant: " << inv;
          }
          this->add_warning(os.str(), &s);
        } else {
          this->m_db.add(_SAFE);
        }
      } else
        CRAB_ERROR("DivZero only supports constant or single var as divisor.");
    }

    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }
};
} // namespace checker
} // namespace crab
