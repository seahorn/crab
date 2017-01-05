#ifndef ASSERT_PROPERTY_CHECKER_HPP
#define ASSERT_PROPERTY_CHECKER_HPP

/* 
   User-definable assertion checker
 */

#include <crab/common/types.hpp>
#include <crab/checkers/base_property.hpp>

namespace crab {

  namespace checker {

    template<typename Analyzer>
    class assert_property_checker: public property_checker <Analyzer> {
      
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
      using typename base_checker_t::assign_t;
      using typename base_checker_t::assume_t;
      using typename base_checker_t::assert_t;
      using typename base_checker_t::bool_assert_t;      

     public:
      
      assert_property_checker (int verbose = 0): base_checker_t (verbose) { }
      
      virtual std::string get_property_name () const override {
        return "user-defined assertion checker using " + abs_dom_t::getDomainName ();
      }

      virtual void check (assert_t& s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        lin_cst_t cst = s.constraint ();

        // Answering a reachability question
        if (cst.is_contradiction()) {
          if (inv.is_bottom()) {
            LOG_SAFE(this->m_verbose, inv, cst, s.get_debug_info());
          } else {
            LOG_WARN(this->m_verbose, inv, cst, s.get_debug_info());
          }
          return;
        }
          
        if (inv.is_bottom ()) {
          this->m_db.add (_UNREACH);
          return;
        }
        
        num_dom_detail::checker_ops<abs_dom_t> cinv (inv);
        if (cinv.entails(cst)) {
          LOG_SAFE(this->m_verbose, inv, cst, s.get_debug_info());
        } else if (cinv.intersect (cst)) {
          LOG_WARN(this->m_verbose, inv, cst, s.get_debug_info());
        } else {
          LOG_ERR(this->m_verbose, inv, cst, s.get_debug_info());
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }


      virtual void check (bool_assert_t& s) override { 
        if (!this->m_abs_tr) return;        
        
	CRAB_WARN ("TODO: checker with bool assertions");
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }
      
    }; 
  } // end namespace
} // end namespace
#endif 
