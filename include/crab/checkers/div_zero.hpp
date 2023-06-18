#pragma once

/**
 * Property checker for division by zero
 **/

#include <crab/checkers/base_property.hpp>

namespace crab {

namespace checker {

template <typename Analyzer>
class div_zero_property_checker : public property_checker<Analyzer> {

  using varname_t = typename Analyzer::varname_t;
  using number_t = typename Analyzer::number_t;

  using interval_t = ikos::interval<number_t>;
  using abs_dom_t = typename Analyzer::abs_dom_t;
  using base_checker_t = property_checker<Analyzer>;
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

    if (s.op() == cfg::BINOP_SDIV || s.op() == cfg::BINOP_UDIV ||
        s.op() == cfg::BINOP_SREM || s.op() == cfg::BINOP_UREM) {

      auto &inv = this->m_abs_tr->get_abs_value();
      if (inv.is_bottom()) {
        this->m_db.add(check_kind::CRAB_UNREACH, s.get_debug_info());
        return;
      }

      auto divisor_expr = s.right();
      if (divisor_expr.is_constant()) {
        number_t divisor = divisor_expr.constant();
        if (divisor == number_t(0)) {
          this->m_db.add(check_kind::CRAB_ERR, s.get_debug_info());
        } else {
          this->m_db.add(check_kind::CRAB_SAFE, s.get_debug_info());
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
            this->m_db.add(check_kind::CRAB_SAFE, s.get_debug_info());
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
          this->m_db.add(check_kind::CRAB_SAFE, s.get_debug_info());
        }
      } else
        CRAB_ERROR("DivZero only supports constant or single var as divisor.");
    }

    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }
};
} // namespace checker
} // namespace crab
