#ifndef DIVISION_ZERO_PROPERTY_CHECKER_HPP
#define DIVISION_ZERO_PROPERTY_CHECKER_HPP

/* 
   Property checker for division by zero
 */

#include <crab/common/types.hpp>
#include <crab/checkers/base_property.hpp>

namespace crab {

  namespace checker {

    template<typename Analyzer>
    class div_zero_property_checker: public property_checker <Analyzer> {
      
      typedef typename Analyzer::varname_t varname_t;
      typedef typename Analyzer::number_t number_t;
      
      typedef ikos::interval<number_t> interval_t;
      typedef typename Analyzer::abs_dom_t abs_dom_t;
      typedef property_checker<Analyzer> base_checker_t;
      using typename base_checker_t::var_t;
      using typename base_checker_t::lin_exp_t;
      using typename base_checker_t::lin_cst_t;
      using typename base_checker_t::lin_cst_sys_t;
      using typename base_checker_t::bin_op_t;

     public:
      
      div_zero_property_checker (int verbose = 0)
          : base_checker_t (verbose) { }
      
      std::string get_property_name () const override {
        return "integer division by zero checker";
      }

      void check (bin_op_t &s) override { 
        if (!this->m_abs_tr) return;        

        if (s.op () == BINOP_SDIV || s.op () == BINOP_UDIV ||
            s.op () == BINOP_SREM || s.op () == BINOP_UREM) {
         
          auto &inv = this->m_abs_tr->inv ();
          if (inv.is_bottom ()) {
            this->m_db.add (_UNREACH);
            return;
          }

          auto divisor_expr = s.right ();
          if (divisor_expr.is_constant ()) {
            number_t divisor = divisor_expr.constant ();
            if (divisor == number_t (0)) {
              this->m_db.add (_ERR);             
            } else {
              this->m_db.add (_SAFE);
            }
          }
          else if (auto var = divisor_expr.get_variable ()) {
            num_dom_detail::checker_ops<abs_dom_t> num_inv (inv);
            interval_t divisor_intv = num_inv [(*var).name()];
            if (auto divisor = divisor_intv.singleton ()) {
              if (*divisor == number_t (0)) {
		lin_cst_t e_cst(*var != number_t (0));
                LOG_ERR(this->m_verbose, inv, e_cst, s.get_debug_info());
              } else {
                this->m_db.add (_SAFE);
              }
            } else if (interval_t (number_t (0)) <= divisor_intv) {
	      lin_cst_t w_cst(*var != number_t (0));
	      LOG_WARN(this->m_verbose, inv, w_cst, s.get_debug_info());
                       
            } else {
              this->m_db.add (_SAFE);
            }
          }
          else 
            CRAB_ERROR ("DivZero only supports constant or single var as divisor.");
        }

        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }      

  }; 
  } // end namespace
} // end namespace
#endif 
