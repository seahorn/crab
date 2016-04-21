#ifndef ASSERT_PROPERTY_CHECKER_HPP
#define ASSERT_PROPERTY_CHECKER_HPP

/* 
   User-definable assertion checker
 */

#include <crab/common/types.hpp>
#include <crab/checkers/BaseProperty.hpp>

namespace crab {

  namespace checker {

    template<typename Analyzer>
    class AssertPropertyChecker: public PropertyChecker <Analyzer> {
      
      typedef typename Analyzer::varname_t varname_t;
      typedef crab::domains::interval<z_number> interval_t;
      typedef typename Analyzer::abs_dom_t abs_dom_t;
      typedef PropertyChecker<Analyzer> base_checker_t;

      using typename base_checker_t::z_var_t;
      using typename base_checker_t::z_lin_exp_t;
      using typename base_checker_t::z_lin_cst_t;
      using typename base_checker_t::z_lin_cst_sys_t;
      using typename base_checker_t::z_bin_op_t;
      using typename base_checker_t::z_assign_t;
      using typename base_checker_t::z_assume_t;
      using typename base_checker_t::z_assert_t;

     public:
      
      AssertPropertyChecker (int verbose = 0): base_checker_t (verbose) { }
      
      virtual std::string get_property_name () const override {
        return "user-defined assertion checker";
      }

      virtual void check (z_assert_t& s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        if (inv.is_bottom ()) {
          this->m_db.add (_UNREACH);
          return;
        }
        
        num_dom_detail::GetAs<abs_dom_t,varname_t> wrapper_inv (inv);
        z_lin_cst_t cst = s.constraint ();
        if (wrapper_inv.entails(cst)) {
          LOG_SAFE(this->m_verbose, inv, cst, s.getDebugInfo());
        } else if (wrapper_inv.intersect (cst)) {
          LOG_WARN(this->m_verbose, inv, cst, s.getDebugInfo());
        } else {
          LOG_ERR(this->m_verbose, inv, cst, s.getDebugInfo());
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }

    }; 
  } // end namespace
} // end namespace
#endif 
