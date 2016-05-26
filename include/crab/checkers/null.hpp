#ifndef NULL_PROPERTY_CHECKER_HPP
#define NULL_PROPERTY_CHECKER_HPP

/* 
   Property checker for null-dereference
 */

#include <crab/checkers/base_property.hpp>
#include <crab/domains/nullity.hpp>

namespace crab {

  namespace checker {

     namespace null_detail {

       // To avoid compilation errors if the Analyzer's abstract domain
       // is not compatible with nullity information.
       template<typename Domain, typename VariableName>
       struct get_as {
         typedef crab::domains::nullity_value nullity_value_t;
         get_as (Domain&) { }
         nullity_value_t operator[](VariableName v){
           return nullity_value_t::top ();
         }
       };
     
       template<typename VariableName>
       struct get_as <crab::domains::nullity_domain<VariableName>, VariableName> {
         typedef crab::domains::nullity_value nullity_value_t;
         typedef crab::domains::nullity_domain<VariableName> nullity_domain_t;
         nullity_domain_t& m_inv;
         get_as (nullity_domain_t& inv) : m_inv (inv) { }
         nullity_value_t operator[](VariableName v){
           return m_inv [v];
         }
       };
     } // end null_detail namespace

    template<typename Analyzer>
    class null_property_checker: public property_checker <Analyzer> {
      
      typedef typename Analyzer::varname_t varname_t;
      typedef crab::domains::nullity_domain<varname_t> nullity_domain_t;
      typedef typename Analyzer::abs_dom_t abs_dom_t;
      typedef property_checker<Analyzer> base_checker_t;
      using typename base_checker_t::z_var_t;
      using typename base_checker_t::z_lin_exp_t;
      using typename base_checker_t::z_lin_cst_t;
      using typename base_checker_t::z_lin_cst_sys_t;
      using typename base_checker_t::z_ptr_load_t;
      using typename base_checker_t::z_ptr_store_t;

     public:
      
      null_property_checker (int verbose = 0): base_checker_t (verbose) { }
      
      std::string get_property_name () const override {
        return "null-dereference checker";
      }
            
      void check (z_ptr_store_t &s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        auto ptr = s.lhs ();
        null_detail::get_as<abs_dom_t,varname_t> null_inv (inv);
        crab::domains::nullity_value val = null_inv [ptr];

        if (val.is_bottom ()) {
          this->m_db.add (_UNREACH);
        } else if (val.is_non_null ()) {
          this->m_db.add (_SAFE);
        } else if (val.is_null ()) {
          this->m_db.add (_ERR);
        } else {
          this->m_db.add (_WARN);
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }
      
      void check (z_ptr_load_t &s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        auto ptr = s.lhs ();
        null_detail::get_as<abs_dom_t, varname_t> null_inv (inv);
        crab::domains::nullity_value val = null_inv [ptr];

        if (val.is_bottom ()) {
          this->m_db.add (_UNREACH);
        } else if (val.is_non_null ()) {
          this->m_db.add (_SAFE);
        } else if (val.is_null ()) {
          this->m_db.add (_ERR);
        } else {
          this->m_db.add (_WARN);
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }

      
  }; 
  } // end namespace
} // end namespace
#endif 
