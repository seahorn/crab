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


#ifndef HAVE_LDD

#define LDD_NOT_FOUND "No LDD. Run cmake with -DUSE_LDD=ON"
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
#include <crab/domains/ldd/ldd.hpp>

using namespace crab::domains::ldd;

namespace crab {

   namespace domains {
       
      class LddManagerWrapper {
        public:
         static LddManager* get () {
           if (!m_ldd_man) {
             const unsigned Sz = 100; // FIXME: make me template parameter
             DdManager* cudd = Cudd_Init (0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
             theory_t* theory = tvpi_create_boxz_theory (Sz);
             m_ldd_man = Ldd_Init (cudd, theory);

             Cudd_AutodynEnable (cudd, CUDD_REORDER_GROUP_SIFT);
           }
           return m_ldd_man;
         }
        
        ~LddManagerWrapper () {
          DdManager *cudd = NULL;
          theory_t *theory = NULL;
          if (m_ldd_man) {
            cudd = Ldd_GetCudd (m_ldd_man);
            theory = Ldd_GetTheory (m_ldd_man);
            Ldd_Quit(m_ldd_man);
          }
          
          if (theory) tvpi_destroy_theory(theory);
          if (cudd) Cudd_Quit(cudd);
        }

        private:
         LddManagerWrapper(){};
         LddManagerWrapper(LddManagerWrapper const&);  
         LddManagerWrapper& operator=(LddManagerWrapper const&); 
         static LddManager* m_ldd_man;
      };

      #define LDD_MAN LddManagerWrapper::get()
   
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

        LddNodePtr m_ldd;
        boxes_domain (LddNodePtr ldd): m_ldd (ldd) { }

        LddNodePtr approx (LddNodePtr v) {
          return lddPtr (LDD_MAN, Ldd_TermMinmaxApprox(LDD_MAN, &*v));
        }

        void dump () {
          LddManager *ldd_man = getLddManager (m_ldd);
          DdManager *cudd = Ldd_GetCudd (ldd_man);
          FILE *fp = Cudd_ReadStdout(cudd);
          Cudd_SetStdout(cudd,stderr);
          if (m_ldd.get () == Ldd_GetTrue (ldd_man)) cout << "true\n";
          else if (m_ldd.get () == Ldd_GetFalse (ldd_man)) cout << "false\n";	
          else Ldd_PrintMinterm(ldd_man, m_ldd.get ());
          Cudd_SetStdout(cudd,fp);      
        }
        
       public:
        
        boxes_domain(): ikos::writeable() { 
          m_ldd = lddPtr (LDD_MAN, Ldd_GetTrue (LDD_MAN));
        }
        
        ~boxes_domain () { }
                
        static boxes_domain_t top() { 
          return boxes_domain_t (lddPtr (LDD_MAN, Ldd_GetTrue (LDD_MAN)));
                                         
        }
        
        static boxes_domain_t bottom() {
          return boxes_domain_t (lddPtr (LDD_MAN, Ldd_GetFalse (LDD_MAN)));
                                         
        }
        
        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable(), m_ldd (other.m_ldd) { }
        
        boxes_domain_t& operator=(boxes_domain_t other) {
           m_ldd = other.m_ldd;
           return *this;
        }
        
        bool is_bottom() { 
          return &*m_ldd == Ldd_GetFalse (LDD_MAN);
        }
        
        bool is_top() { 
          return &*m_ldd == Ldd_GetTrue (LDD_MAN);
        }
        
        bool operator<=(boxes_domain_t other) {
          return Ldd_TermLeq (LDD_MAN, &(*m_ldd), &(*other.m_ldd));
        }
        
        boxes_domain_t operator|(boxes_domain_t other) {
          return boxes_domain_t (approx (lddPtr (LDD_MAN, 
                                                 Ldd_Or (LDD_MAN, &*m_ldd, &*other.m_ldd))));
        }
        
        boxes_domain_t operator&(boxes_domain_t other) {
          return boxes_domain_t (lddPtr (LDD_MAN, 
                                         Ldd_And (LDD_MAN, &*m_ldd, &*other.m_ldd)));
        }
        
        boxes_domain_t operator||(boxes_domain_t other) {
          return boxes_domain_t (lddPtr (LDD_MAN, 
                                         Ldd_IntervalWiden (LDD_MAN, &*m_ldd, &*other.m_ldd)));
        }
        
        boxes_domain_t operator&& (boxes_domain_t other) {
          CRAB_WARN("TODO: BOXES narrowing operation");
          return *this;
        }
        
        void operator-=(VariableName var) {
          CRAB_WARN("TODO: BOXES forget operation");
        }

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          CRAB_WARN("TODO: BOXES forget operation");
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          CRAB_WARN("TODO: BOXES project operation");
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
          CRAB_WARN("TODO: BOXES domain operation");
        }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          CRAB_WARN("TODO: BOXES domain operation");
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
          CRAB_WARN("TODO: BOXES domain operation");
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", z, *this);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          CRAB_WARN("TODO: BOXES domain operation");
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, *this);
        }
        
        linear_constraint_system_t to_linear_constraint_system (){
          CRAB_ERROR ("TODO: BOXES domain operation");
        }
        
        void write(ostream& o) {
          auto csts = to_linear_constraint_system ();
          o << csts;
        }
        
        const char* getDomainName () const {return "BOXES ";}  
        
      }; 
   
   } // namespace domains

   namespace domain_traits {
    
     template <typename VariableName, typename Number, typename Iterator >
     void forget (boxes_domain<Number, VariableName>& inv, 
                  Iterator it, Iterator end) {
       inv.forget (it, end);
     }
   
     template <typename VariableName, typename Number, typename Iterator >
     void project (boxes_domain<Number, VariableName>& inv, 
                   Iterator it, Iterator end) {
       inv.project (it, end);
     }

   } // namespace domain_traits
}// namespace crab
#endif /* HAVE_LDD */
#endif 
