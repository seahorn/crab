#ifndef BOXES_DOMAIN_HPP
#define BOXES_DOMAIN_HPP

/* 
   Disjunctive intervals based on "Boxes: A Symbolic Abstract Domain
   of Boxes" by A. Gurfinkel and S. Chaki in SAS'10
*/


// Uncomment for enabling debug information
// #include <crab/common/dbg.hpp>

#include <crab/config.h>

#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/domain_traits_impl.hpp>

using namespace boost;
using namespace ikos;

#define LDD_NOT_FOUND "No LDD. Run cmake with -DUSE_LDD=ON"

#ifndef HAVE_LDD
namespace crab {
   namespace domains {
      template<typename Number, typename VariableName>
      class boxes_domain: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName > {
              
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;
        typedef boxes_domain <Number, VariableName> boxes_domain_t;

        boxes_domain(): ikos::writeable() { }    

        static boxes_domain_t top() { CRAB_ERROR (LDD_NOT_FOUND); }

        static boxes_domain_t bottom() { CRAB_ERROR (LDD_NOT_FOUND); }

        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable() { }
        
        bool is_bottom() { CRAB_ERROR (LDD_NOT_FOUND); }

        bool is_top() { CRAB_ERROR (LDD_NOT_FOUND); }

        bool operator<=(boxes_domain_t other) { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator|(boxes_domain_t other)
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator&(boxes_domain_t other) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator||(boxes_domain_t other)
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator&& (boxes_domain_t other) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void operator-=(VariableName var) 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        void operator += (linear_constraint_system_t csts) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void assign (VariableName x, linear_expression_t e) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
          
        void apply (operation_t op, VariableName x, VariableName y, Number z) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(operation_t op, VariableName x, Number k) 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        linear_constraint_system_t to_linear_constraint_system ()
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void write(ostream& o) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
          
        const char* getDomainName () const {return "BOXES";}  
      }; 
   
   } // namespace domains
}// namespace crab

#else

// Real implementation starts here
namespace crab {

   namespace domains {

      template<typename Number, typename VariableName>
      class boxes_domain: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName > {
              
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;

        typedef boxes_domain <Number, VariableName> boxes_domain_t;
        typedef interval <Number> interval_t;
        
       private:
        
        typedef bound <Number> bound_t; 
        
        
       public:
        
        boxes_domain(): ikos::writeable() { }    
        
        static boxes_domain_t top() { 
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        static boxes_domain_t bottom() {
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable() { }
        
        // boxes_domain_t& operator=(boxes_domain_t other) {
        //   _inv = other._inv;
        //   return *this;
        // }
        
        bool is_bottom() { 
          return false;
        }
        
        bool is_top() { 
          return true;
        }
        
        bool operator<=(boxes_domain_t other) {
          return true;
        }
        
        boxes_domain_t operator|(boxes_domain_t other) {
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        boxes_domain_t operator&(boxes_domain_t other) {
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        boxes_domain_t operator||(boxes_domain_t other) {
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        boxes_domain_t operator&& (boxes_domain_t other) {
          CRAB_ERROR("TODO: BOXES domain operation");
        }
        
        void operator-=(VariableName var) {
          CRAB_WARN("TODO: BOXES domain operation");
        }

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          CRAB_WARN("TODO: BOXES domain operation");
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          CRAB_WARN("TODO: BOXES domain operation");
        }

        void operator += (linear_constraint_system_t csts) {
          CRAB_WARN("TODO: BOXES domain operation");
        }
        
        void assign (VariableName x, linear_expression_t e) {
          CRAB_WARN("TODO: BOXES domain operation");

          CRAB_DEBUG("apply ", x, " := ", e, *this);
        }
        
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          CRAB_WARN("TODO: BOXES domain operation");
          
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          CRAB_WARN("TODO: BOXES domain operation");
          
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          CRAB_WARN("TODO: BOXES domain operation");
          
          CRAB_DEBUG("apply ", x, " := ", x, " ", op, " ", k, *this);
        }

        // bitwise_operators_api
        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
          //_inv.apply (op, x, y, width);
        }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          //_inv.apply (op, x, k, width);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          CRAB_WARN("Bitwise operation ", op, " not implemented");
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          CRAB_WARN("Bitwise operation ", op, " not implemented");
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, *this);
        }
        
        // division_operators_api
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          //_inv.apply (op, x, y, z);
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          //_inv.apply (op, x, y, k);
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, *this);
        }
        
        linear_constraint_system_t to_linear_constraint_system (){
          CRAB_ERROR ("TODO: BOXES domain operation");
        }
        
        void write(ostream& o) {
          CRAB_ERROR ("TODO: BOXES domain operation");
        }
        
        const char* getDomainName () const {return "BOXES ";}  
        
      }; 
   
   } // namespace domains

   namespace domain_traits {
    
     template <typename VariableName, typename Number, typename Iterator >
     void forget (boxes_domain<Number, VariableName>& inv, 
                  Iterator it, Iterator end) {
       //inv.forget (it, end);
     }
   
     template <typename VariableName, typename Number, typename Iterator >
     void project (boxes_domain<Number, VariableName>& inv, 
                   Iterator it, Iterator end) {
       //inv.project (it, end);
     }

   } // namespace domain_traits
}// namespace crab
#endif /* HAVE_LDD */
#endif 
