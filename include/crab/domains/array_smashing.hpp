/*******************************************************************************
 * Array smashing domain
 * 
 * FIXME: assume all array accesses are aligned wrt to the size of the
 * array element (e.g., if the size of the array element is 4 bytes
 * then all array accesses must be multiple of 4). Note that this
 * assumption does not hold in real programs.
 ******************************************************************************/

#ifndef ARRAY_SMASHING_HPP
#define ARRAY_SMASHING_HPP

// Uncomment for enabling debug information
// #include <crab/common/dbg.hpp>

#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/domain_traits_impl.hpp>

using namespace boost;
using namespace ikos;

namespace crab {

   namespace domains {

      //! Abstract domain to reason about summarized variables. All array
      //  elements are `smashed` into a single cell.
      template<typename NumDomain>
      class array_smashing: 
         public ikos::writeable, 
         public numerical_domain<typename NumDomain::number_t,
                                 typename NumDomain::varname_t>,
         public bitwise_operators<typename NumDomain::number_t, 
                                  typename NumDomain::varname_t>, 
         public division_operators<typename NumDomain::number_t,
                                   typename NumDomain::varname_t> {
              
       public:

        typedef typename NumDomain::number_t Number;
        typedef typename NumDomain::varname_t VariableName;

        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;

        typedef array_smashing <NumDomain> array_smashing_t;
        typedef NumDomain content_domain_t;

        typedef interval <Number> interval_t;
        
       private:
        
        typedef bound <Number> bound_t; 
        
        //! scalar and summarized array variables        
        NumDomain _inv; 
        
        array_smashing (NumDomain inv): 
            ikos::writeable (), 
            _inv (inv) { }
        
        void strong_update (VariableName lhs, linear_expression_t rhs ) {
          _inv.assign (lhs, rhs);
        }
        
        void weak_update (VariableName lhs, linear_expression_t rhs) {
          NumDomain other (_inv);
          other.assign (lhs, rhs);
          _inv = _inv | other;
        }
        
       public:
        
        array_smashing(): ikos::writeable(), _inv (NumDomain::top()) { }    
        
        static array_smashing_t top() { 
          return array_smashing (NumDomain::top ()); 
        }
        
        static array_smashing_t bottom() {
          return array_smashing (NumDomain::bottom ());
        }
        
        array_smashing (const array_smashing_t& other): 
            ikos::writeable(), 
            _inv (other._inv) { }
        
        array_smashing_t& operator=(array_smashing_t other) {
          _inv = other._inv;
          return *this;
        }
        
        bool is_bottom() { 
          return (_inv.is_bottom ());
        }
        
        bool is_top() { 
          return (_inv.is_top());
        }
        
        bool operator<=(array_smashing_t other) {
          return (_inv <= other._inv);
        }

        void operator|=(array_smashing_t other) {
          *this = *this | other;
        }
        
        array_smashing_t operator|(array_smashing_t other) {
          return array_smashing_t (_inv | other._inv);
        }
        
        array_smashing_t operator&(array_smashing_t other) {
          return array_smashing_t (_inv & other._inv);
        }
        
        array_smashing_t operator||(array_smashing_t other) {
          return array_smashing_t (_inv || other._inv);
        }

        template<typename Thresholds>
        array_smashing_t widening_thresholds (array_smashing_t other, 
                                              const Thresholds &ts) {
          return array_smashing_t (_inv.widening_thresholds (other._inv, ts));
        }
        
        array_smashing_t operator&& (array_smashing_t other) {
          return array_smashing_t (_inv && other._inv);
        }
        
        void operator-=(VariableName var) {
          _inv -= var;
        }

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          crab::domain_traits::forget (_inv, begin, end);
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          crab::domain_traits::project (_inv, begin, end);
        }

        void operator += (linear_constraint_system_t csts) {
          _inv += csts;
        }
        
        void assign (VariableName x, linear_expression_t e) {
          _inv.assign (x, e);
          
          CRAB_DEBUG("apply ", x, " := ", e, *this);
        }
        
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          _inv.apply (op, x, y, z);
          
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);
          
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          _inv.apply (op, x, k);
          
          CRAB_DEBUG("apply ", x, " := ", x, " ", op, " ", k, *this);
        }

        // bitwise_operators_api
        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
          _inv.apply (op, x, y, width);
        }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          _inv.apply (op, x, k, width);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);

          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          _inv.apply (op, x, y, k);

          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, *this);
        }
        
        // division_operators_api
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);

          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          _inv.apply (op, x, y, k);

          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, *this);
        }
        
        void array_init (VariableName a, 
                         const vector<ikos::z_number>& values) {
          if (values.empty ()) return;
          
          interval_t init = interval_t::bottom ();
          for (auto const &v: values) {
            // assume automatic conversion from z_number to bound_t
            init = init | interval_t (bound_t (v)); 
          }
          _inv.set (a, init);
          CRAB_DEBUG("Array init: ",*this);
        }
        
        // All the array elements are initialized to val
        void assume_array (VariableName a, interval_t val) {
          _inv.set (a, val);
          
          CRAB_DEBUG("Assume array: ",*this);
        }
        
        void load (VariableName lhs, VariableName a, 
                   VariableName /*i*/, z_number /*n_bytes*/) {
          
          // We need to be careful when assigning a summarized variable a
          // into a non-summarized variable lhs. Simply _inv.assign (lhs,
          // a) is not sound.
          /* ask for a temp var */
          VariableName a_prime = a.getVarFactory().get(); 
          crab::domain_traits::expand (_inv, a, a_prime);
          _inv.assign (lhs, linear_expression_t (a_prime));
          _inv -= a_prime; 
          
          CRAB_DEBUG("Load: ",*this);
        }
        
        
        void store (VariableName a, VariableName /*i*/,
                    linear_expression_t val, z_number /*n_bytes*/,
                    bool is_singleton) {
          
          if (is_singleton)
            strong_update (a, val);
          else 
            weak_update (a, val);
          
          CRAB_DEBUG("Store: ",*this);
        }
        
        linear_constraint_system_t to_linear_constraint_system (){
          return _inv.to_linear_constraint_system ();
        }
        
        NumDomain  get_content_domain () const {      
          return _inv;
        }

        void write(ostream& o) {
          o << _inv;
        }
        
        static const char* getDomainName () {
          std::string name ("ArraySmashing(" + 
                            std::string(NumDomain::getDomainName ()) + 
                            ")");
          return name.c_str ();
        }  
        
      }; // end array_smashing
   
   } // namespace domains

   namespace domain_traits {
    
     using namespace domains;            

     template <typename BaseDomain>
     void array_init (array_smashing<BaseDomain>& inv, 
                      typename BaseDomain::varname_t a, 
                      const vector<ikos::z_number> &values) {
       inv.array_init (a, values);
     }
   
     template <typename BaseDomain>
     void assume_array (array_smashing<BaseDomain>& inv, 
                        typename BaseDomain::varname_t a, 
                        typename BaseDomain::number_t val) {
       inv.assume_array (a, interval<typename BaseDomain::number_t> 
                         (bound <typename BaseDomain::number_t> (val)));
     }
   
     template <typename BaseDomain>
     void assume_array (array_smashing<BaseDomain>& inv, 
                        typename BaseDomain::varname_t a, 
                        interval<typename BaseDomain::number_t> val) {
       inv.assume_array (a, val);
     }
   
     template <typename BaseDomain>
     void array_load (array_smashing<BaseDomain>& inv, 
                      typename BaseDomain::varname_t lhs, 
                      typename BaseDomain::varname_t a, 
                      typename BaseDomain::varname_t i, 
                      z_number n_bytes) {
       inv.load (lhs, a, i, n_bytes);
     }
   
     template <typename BaseDomain>
     void array_store (array_smashing<BaseDomain>& inv, 
                       typename BaseDomain::varname_t a, 
                       typename BaseDomain::varname_t i,
                       typename BaseDomain::linear_expression_t val,
                       z_number n_bytes, bool is_singleton) {
       inv.store (a, i, val, n_bytes, is_singleton);
     }

     template <typename BaseDomain, typename Iterator >
     void forget (array_smashing<BaseDomain>& inv, 
                  Iterator it, Iterator end) {
       inv.forget (it, end);
     }
   
     template <typename BaseDomain, typename Iterator >
     void project (array_smashing<BaseDomain>& inv, 
                   Iterator it, Iterator end) {
       inv.project (it, end);
     }

   } // namespace domain_traits
}// namespace crab
#endif 
