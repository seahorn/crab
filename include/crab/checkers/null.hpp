#ifndef NULL_PROPERTY_CHECKER_HPP
#define NULL_PROPERTY_CHECKER_HPP

/* 
   Property checker for null-dereference
 */

#include <crab/checkers/base_property.hpp>
#include <crab/domains/domain_products.hpp>

namespace crab {

  namespace checker {

     namespace null_detail {

       // arbitrary abstract domain without nullity information
       template<typename Domain>
       struct get_as {
         typedef typename Domain::varname_t varname_t;
         typedef crab::domains::nullity_value nullity_value_t;
         get_as (Domain&) { }
         nullity_value_t operator[](varname_t v)
         { return nullity_value_t::top (); }
       };
  
       // nullity domain
       template<typename Number, typename VariableName>
       class get_as <crab::domains::nullity_domain<Number, VariableName> > {
         typedef crab::domains::nullity_value nullity_value_t;
         typedef crab::domains::nullity_domain<Number, VariableName> nullity_domain_t;

         nullity_domain_t &m_inv;

        public:

         get_as (nullity_domain_t &inv) : m_inv (inv) { }

         nullity_value_t operator[](VariableName v)
         { return m_inv [v]; }
       };

       // Reduced product of a numerical abstract domain with nullity
       template<typename Dom, typename Number, typename VariableName>
       class get_as <ikos::numerical_domain_product2<Number, 
                      VariableName, Dom, crab::domains::nullity_domain<Number, VariableName> > > {

         typedef crab::domains::nullity_value nullity_value_t;
         typedef crab::domains::nullity_domain<Number, VariableName> nullity_domain_t;
         typedef ikos::numerical_domain_product2<Number, VariableName, Dom, nullity_domain_t>
         numerical_domain_product2_t;
         
         numerical_domain_product2_t &m_inv;

        public:

         get_as (numerical_domain_product2_t& inv) : m_inv (inv) { }

         nullity_value_t operator[](VariableName v)
         { return m_inv.second() [v]; }
       };

     } // end null_detail namespace

    template<typename Analyzer>
    class null_property_checker: public property_checker <Analyzer> {

      typedef typename Analyzer::abs_dom_t abs_dom_t;      
      typedef typename abs_dom_t::varname_t varname_t;
      typedef typename abs_dom_t::number_t number_t;
      typedef crab::domains::nullity_domain<number_t, varname_t> nullity_domain_t;
      typedef property_checker<Analyzer> base_checker_t;
      using typename base_checker_t::z_var_t;
      using typename base_checker_t::z_lin_exp_t;
      using typename base_checker_t::z_lin_cst_t;
      using typename base_checker_t::z_lin_cst_sys_t;
      using typename base_checker_t::ptr_load_t;
      using typename base_checker_t::ptr_store_t;

      std::string checked_prop_str (varname_t p) {
        std::string res = p.str () + " != 0";
        return res;
      }

     public:
      
      null_property_checker (int verbose = 0): base_checker_t (verbose) { }
      
      std::string get_property_name () const override {
        return "null-dereference checker";
      }
            
      void check (ptr_store_t &s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        auto ptr = s.lhs ();
        null_detail::get_as<abs_dom_t> null_inv (inv);
        crab::domains::nullity_value val = null_inv [ptr];

        if (val.is_bottom ()) {
          this->m_db.add (_UNREACH);
        } else if (val.is_non_null ()) {
          LOG_SAFE(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
          //this->m_db.add (_SAFE);
        } else if (val.is_null ()) {
          LOG_ERR(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
          //this->m_db.add (_ERR);
        } else {
          LOG_WARN(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
          //this->m_db.add (_WARN);
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }
      
      void check (ptr_load_t &s) override { 
        if (!this->m_abs_tr) return;        
        
        auto &inv = this->m_abs_tr->inv ();
        auto ptr = s.rhs ();
        null_detail::get_as<abs_dom_t> null_inv (inv);
        crab::domains::nullity_value val = null_inv [ptr];
        
        if (val.is_bottom ()) {
          this->m_db.add (_UNREACH);
        } else if (val.is_non_null ()) {
          LOG_SAFE(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
        } else if (val.is_null ()) {
          LOG_ERR(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
        } else {
          LOG_WARN(this->m_verbose, inv, checked_prop_str(ptr), s.get_debug_info());
        }
        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
      }

      
  }; 
  } // end namespace
} // end namespace
#endif 
