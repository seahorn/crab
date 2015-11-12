#ifndef APRON_DOMAINS_HPP
#define APRON_DOMAINS_HPP

/// Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/domain_traits_impl.hpp>

using namespace boost;
using namespace ikos;

namespace crab {
   namespace domains {
      typedef enum { INT, OCT, PK } apron_domain_id_t;
   }
}

#ifndef HAVE_APRON
/*
 * Dummy implementation if Apron not found 
 */
#define APRON_NOT_FOUND "No Apron. Run cmake with -DUSE_APRON=ON"

namespace crab {
   namespace domains {
      template<typename Number, typename VariableName, apron_domain_id_t ApronDom>
      class apron_domain: 
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
        typedef apron_domain <Number, VariableName, ApronDom> apron_domain_t;
        typedef interval <Number> interval_t;

        apron_domain(): ikos::writeable() { }    

        static void addTrackVar (VariableName v) 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        static void resetTrackVars () 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        static apron_domain_t top() { CRAB_ERROR (APRON_NOT_FOUND); }

        static apron_domain_t bottom() { CRAB_ERROR (APRON_NOT_FOUND); }

        apron_domain (const apron_domain_t& other): 
            ikos::writeable() { }
        
        bool is_bottom() { CRAB_ERROR (APRON_NOT_FOUND); }

        bool is_top() { CRAB_ERROR (APRON_NOT_FOUND); }

        bool operator<=(apron_domain_t other) { CRAB_ERROR (APRON_NOT_FOUND); }
        
        apron_domain_t operator|(apron_domain_t other)
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        apron_domain_t operator&(apron_domain_t other) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        apron_domain_t operator||(apron_domain_t other)
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        apron_domain_t operator&& (apron_domain_t other) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void operator-=(VariableName var) 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        interval_t operator[](VariableName v) 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        void set(VariableName v, interval_t ival) 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        void operator += (linear_constraint_system_t csts) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void assign (VariableName x, linear_expression_t e) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
          
        void apply (operation_t op, VariableName x, VariableName y, Number z) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(operation_t op, VariableName x, Number k) 
        { CRAB_ERROR (APRON_NOT_FOUND); }

        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        linear_constraint_system_t to_linear_constraint_system ()
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
        void write(ostream& o) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
          
        const char* getDomainName () const {return "Dummy Apron";}  
      }; 
   } // namespace domains
}// namespace crab
#else

/* 
 *  Real implementation starts here 
 */

#include <crab/domains/apron/apron.hpp>
#include <boost/bimap.hpp>

namespace crab {
   namespace domains {

      template<typename Number, typename VariableName, apron_domain_id_t ApronDom>
      class apron_domain: 
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
        typedef apron_domain <Number, VariableName, ApronDom> apron_domain_t;
        typedef interval <Number> interval_t;

       private:
        typedef interval_domain <Number, VariableName> interval_domain_t;
        typedef boost::bimap< VariableName , ap_dim_t > var_bimap_t;
        typedef boost::shared_ptr<var_bimap_t> var_map_ptr;
        typedef typename var_bimap_t::value_type binding_t;
        
        ap_state_ptr m_apstate;
        static ap_manager_t* m_apman;
        static var_map_ptr m_var_map;
        
        static ap_manager_t* get_man () {
          if (!m_apman) {
            if (ApronDom == INT)
              m_apman = box_manager_alloc ();
            else if (ApronDom == OCT)
              m_apman = oct_manager_alloc ();
            else if (ApronDom == PK)
              m_apman = pk_manager_alloc (false);
            else
              CRAB_ERROR("Unknown apron domain");
          }
          return m_apman;
        }

        static var_map_ptr get_var_map () {
          if (!m_var_map) {
            m_var_map = var_map_ptr (new var_bimap_t ());            
          }
          return m_var_map;
        }

        static size_t get_dims () {
          if (!m_var_map || m_var_map->empty ())  {
            CRAB_WARN ("Apron state with 0 dimensions. ",
                       "Call addTrackVar before apron starts");
            return 0;
          }
          else
            return m_var_map->size ();                               
        }
        
        ap_dim_t get_var_dim (VariableName v) const {
          auto it = get_var_map()->left.find (v);
          if (it != get_var_map()->left.end ())
            return it->second;
          CRAB_ERROR ("Apron could not find ", v, ". ",
                      "Call addTrackVar before apron starts");
        }

        VariableName get_var_name (ap_dim_t i) const {
          auto it = get_var_map()->right.find (i);
          if (it != get_var_map()->right.end ())
            return it->second;
          CRAB_ERROR ("Apron dimension ", i, " out of range!");
        }


        // --- build apron binary operations
        inline ap_texpr0_t* ADD(ap_texpr0_t* a, ap_texpr0_t*b) const {
          return ap_texpr0_binop(AP_TEXPR_ADD,a,b,AP_RTYPE_INT,AP_RDIR_NEAREST);
        }
        
        inline ap_texpr0_t* SUB(ap_texpr0_t* a, ap_texpr0_t*b) const {
          return ap_texpr0_binop(AP_TEXPR_SUB,a,b,AP_RTYPE_INT,AP_RDIR_NEAREST);
        }
        
        inline ap_texpr0_t* MUL(ap_texpr0_t* a, ap_texpr0_t*b) const {
          return ap_texpr0_binop(AP_TEXPR_MUL,a,b,AP_RTYPE_INT,AP_RDIR_NEAREST);
        }
        
        inline ap_texpr0_t* DIV(ap_texpr0_t* a, ap_texpr0_t*b) const {
          return ap_texpr0_binop(AP_TEXPR_DIV,a,b,AP_RTYPE_INT,AP_RDIR_NEAREST);
        }

        // --- from crab to apron

        inline ap_texpr0_t* var2texpr (VariableName v) const { 
          return ap_texpr0_dim( get_var_dim (v));
        }

        inline ap_texpr0_t* num2texpr (Number i) const {  
          return ap_texpr0_cst_scalar_int ((int) i); 
        }

        inline ap_texpr0_t* intv2texpr (Number a, Number b) const { 
          return ap_texpr0_cst_interval_int ((int) a,(int) b); 
        }
        
        inline ap_texpr0_t* expr2texpr (linear_expression_t e) const {
          Number cst = e.constant ();
          ap_texpr0_t* res = num2texpr (cst);
          for (auto p: e) {
            ap_texpr0_t* term = MUL (num2texpr (p.first), var2texpr (p.second.name ()));
            res = ADD (res, term); 
          }
          return res;
        }

        inline ap_tcons0_t const2tconst (linear_constraint_t cst) const {

          assert (!cst.is_tautology ());
          assert (!cst.is_contradiction ());

          linear_expression_t exp = cst.expression();
          if (cst.is_equality ()) {
            return ap_tcons0_make (AP_CONS_EQ, expr2texpr (exp), NULL);            
          }
          else if (cst.is_inequality ()) {
            return ap_tcons0_make (AP_CONS_SUPEQ, expr2texpr (-exp), NULL);
          }
          else  { 
            assert (cst.is_disequation ());
            return ap_tcons0_make (AP_CONS_DISEQ, expr2texpr (exp), NULL);            
          }
        }
          
        // --- from apron to crab 
         
        Number coeff2Num (ap_coeff_t* coeff) {
          assert (coeff->discr == AP_COEFF_SCALAR);
          ap_scalar_t* scalar = coeff->val.scalar;
          assert (scalar->discr == AP_SCALAR_MPQ);
          mpq_ptr c = scalar->val.mpq;
          return Number ((mpz_class) mpq_class(c));
        }

        linear_expression_t term2expr (ap_coeff_t* coeff, ap_dim_t i) {
          return variable_t (get_var_name (i)) * coeff2Num(coeff) ;
        }

        linear_constraint_t tconst2const (ap_lincons0_t cons) {

          assert (cons.scalar == NULL); // Not modulo form
          ap_linexpr0_t* linexp = cons.linexpr0;
          assert (ap_linexpr0_is_linear (linexp));

          linear_expression_t e (0);
          for (unsigned i=0; i < get_dims (); ++i) {
            ap_coeff_t* coeff = ap_linexpr0_coeffref (linexp, i);
            if (ap_coeff_zero (coeff)) continue;
            e = e + term2expr ( coeff, i);
          }
          ap_coeff_t* cst = ap_linexpr0_cstref (linexp);

          linear_constraint_t res;
          switch (cons.constyp) {
            case AP_CONS_EQ:
              // e == k <--> e -k == 0
              if (!ap_coeff_zero (cst)) e = e - coeff2Num(cst);
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::EQUALITY);
              break;
            case AP_CONS_SUPEQ:
              // e >= k <--> -e <= -k <--> -e + k <= 0
              e = Number (0) - e;
              if (!ap_coeff_zero (cst)) e = e + coeff2Num(cst);
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY);
              break;
            case AP_CONS_SUP:
              // e > k <--> -e < -k <--> -e + k < 0 <---> -e + 1 + k <= 0 (only if integers)
              e = Number (0) - e + 1;
              if (!ap_coeff_zero (cst)) e = e + coeff2Num(cst);
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY);
              break;
            case AP_CONS_EQMOD:
              res = T ();
            case AP_CONS_DISEQ:
              // e != k <--> e - k != 0
              if (!ap_coeff_zero (cst)) e = e - coeff2Num(cst);
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::DISEQUATION);
              break;
          }
          return res;
        }

        inline linear_constraint_t T () const {
          return linear_constraint_t (linear_expression_t (Number(1)) == 
                                      linear_expression_t (Number(0)));          
        }

        inline linear_constraint_t F () const {
          return linear_constraint_t (linear_expression_t (Number(1)) == 
                                      linear_expression_t (Number(1)));          
        }

        void dump () {  
          vector<char*> names;
          for (unsigned i=0; i < get_dims (); i++) {
            string varname = get_var_name (i).str ();
            char* name = new char [varname.length () + 1];
            strcpy (name, varname.c_str ());
            names.push_back (name);
          }

          ap_abstract0_fprint (stdout, 
                               get_man (), &*m_apstate, &names[0]);
        }

        apron_domain (ap_state_ptr apState): 
            ikos::writeable (), m_apstate (apState)  { }

       public:

        static void addTrackVar (VariableName v) {
          auto it = get_var_map()->left.find (v);
          if (it == get_var_map()->left.end ()) {
            ap_dim_t id = get_var_map ()->size ();
            get_var_map ()->insert (binding_t (v, id));
          }
        }
            
        static void resetTrackVars () { get_var_map ()->clear (); }

       public:

        apron_domain (): 
            ikos::writeable (),
            m_apstate (apPtr (get_man(), ap_abstract0_top (get_man(), get_dims (), 0))){
        }

        static apron_domain_t top() { 
          return apron_domain_t (apPtr (get_man(), 
                                        ap_abstract0_top (get_man(), get_dims (), 0)));
        }

        static apron_domain_t bottom() { 
          return apron_domain_t (apPtr (get_man(), 
                                        ap_abstract0_bottom (get_man(), get_dims (), 0)));
        }

        apron_domain (const apron_domain_t& other): 
            ikos::writeable(), m_apstate (other.m_apstate) {  }
        
        
        apron_domain_t operator=(const apron_domain_t& o) {
          if (this != &o) {
            m_apstate = o.m_apstate;
          }
          return *this;
        }
        
        bool is_bottom() { 
          return ap_abstract0_is_bottom (get_man(), &*m_apstate);
        }

        bool is_top() { 
          return ap_abstract0_is_top (get_man(), &*m_apstate);
        }

        bool operator<=(apron_domain_t other) { 
          return ap_abstract0_is_leq (get_man(), &*m_apstate, &*other.m_apstate);
        }
        
        apron_domain_t operator|(apron_domain_t other) {
          return apron_domain_t (apPtr (get_man(), 
                                        ap_abstract0_join (get_man(), 
                                                           false,
                                                           &*m_apstate, &*other.m_apstate)));
        }        
        
        apron_domain_t operator&(apron_domain_t other) {
          return apron_domain_t (apPtr (get_man(), 
                                        ap_abstract0_meet (get_man(), 
                                                           false,
                                                           &*m_apstate, &*other.m_apstate)));
        }        
        
        apron_domain_t operator||(apron_domain_t other) {
          return apron_domain_t (apPtr (get_man(), 
                                        ap_abstract0_widening (get_man(), 
                                                               &*m_apstate, &*other.m_apstate)));
        }        
        
        apron_domain_t operator&&(apron_domain_t other) {
          // apron does not define narrowing operators.
          // It might not terminate so make sure that the fixpoint
          // runs a finite number of descending iterations
          return (*this & other);
        }        

        template<typename Range>
        void forget (const Range &vars) {
          vector<ap_dim_t> dims;
          for (auto v: vars) 
            dims.push_back (get_var_dim (v));
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_forget_array (get_man (), 
                                                        false, 
                                                        &*m_apstate, 
                                                        &dims[0], dims.size(), 
                                                        false));
        }

        void operator-=(VariableName var) {
          vector<ap_dim_t> dims;
          dims.push_back (get_var_dim (var));
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_forget_array (get_man (), 
                                                        false, 
                                                        &*m_apstate, 
                                                        &dims[0], dims.size(), 
                                                        false));
        }

        interval_t operator[](VariableName v) {
          ap_interval_t* intv = ap_abstract0_bound_dimension (get_man (),
                                                              &*m_apstate, get_var_dim (v));
          ap_scalar_t* lb = intv->inf;
          ap_scalar_t* ub = intv->sup;

          assert (lb->discr == AP_SCALAR_MPQ);
          mpq_ptr lb_c = lb->val.mpq;
          assert (ub->discr == AP_SCALAR_MPQ);
          mpq_ptr ub_c = ub->val.mpq;
          return interval_t (Number ((mpz_class) mpq_class(lb_c)),
                             Number ((mpz_class) mpq_class(ub_c)));
        }

        void set(VariableName v, interval_t ival) {

          // -- forget v
          *this -= v;

          // -- add constraints v >= lb and v <= ub
          linear_constraint_system_t csts;
          auto lb = ival.lb ();
          if (lb.is_finite ())  {
            // v >= lb <--> -v + lb <= 0
            assert (lb.number ());
            linear_expression_t e = (Number(-1) * linear_expression_t (v)) + *(lb.number ());
            csts += (linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY));
          }
          auto ub = ival.ub ();
          if (ub.is_finite ()) {
            // v <= ub <--> v - ub <= 0
            assert (ub.number ());
            linear_expression_t e = (linear_expression_t (v) - *(ub.number ()));
            csts += (linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY));
          }

          if (csts.size () > 0)
            *this += csts;
        }

        void operator += (linear_constraint_system_t csts) {
          if(is_bottom()) return;

          ap_tcons0_array_t array = ap_tcons0_array_make (csts.size ());
          unsigned i=0;
          for (auto cst : csts) { 
            ap_tcons0_t tcons = const2tconst (cst);
            array.p[i] = tcons;
            ++i;
          }

          m_apstate = apPtr (get_man (), 
                             ap_abstract0_meet_tcons_array (get_man (), false, 
                                                            &*m_apstate, &array));

          //cout << "Added " << csts << " --> "; dump (); cout << "\n";
        }
       
        void assign (VariableName x, linear_expression_t e) {
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_assign_texpr(get_man (), false, 
                                                       &*m_apstate, 
                                                       get_var_dim (x), expr2texpr (e), 
                                                       NULL));

          //cout << x << ":=" << e << " --> "; dump (); cout << "\n";
        }
          
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          ap_texpr0_t* a = var2texpr (y);
          ap_texpr0_t* b = num2texpr (z);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("Apron: unreachable");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim (x), 
                                                                   res, 
                                                                   NULL));

          // cout << x << ":=" << y << op << z << " --> "; dump (); cout << "\n";
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          ap_texpr0_t* a = var2texpr (y);
          ap_texpr0_t* b = var2texpr (z);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("Apron: unreachable");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim (x), 
                                                                   res, 
                                                                   NULL));

          // cout << x << ":=" << y << op << z << " --> "; dump (); cout << "\n";
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          ap_texpr0_t* a = var2texpr (x);
          ap_texpr0_t* b = num2texpr (k);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("Apron: unreachable");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim (x), 
                                                                   res, 
                                                                   NULL));

          // cout << x << ":=" << x << op << k << " --> "; dump (); cout << "\n";
        }

        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
          // since reasoning about infinite precision we simply assign and
          // ignore the width.
          assign(x, linear_expression_t(y));
        }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          // since reasoning about infinite precision we simply assign and
          // ignore the width.
          assign(x, k);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi = operator[](z);
          interval_t xi = interval_t::top();
          switch (op) {
            case OP_AND: xi = yi.And(zi); break;
            case OP_OR: xi = yi.Or(zi); break;
            case OP_XOR: xi = yi.Xor(zi); break; 
            case OP_SHL: xi = yi.Shl(zi); break; 
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
            default: CRAB_ERROR("Apron: unreachable");
          }
          set(x, xi);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi(k);
          interval_t xi = interval_t::top();
          switch (op) {
            case OP_AND: xi = yi.And(zi); break;
            case OP_OR: xi = yi.Or(zi); break;
            case OP_XOR: xi = yi.Xor(zi); break; 
            case OP_SHL: xi = yi.Shl(zi); break; 
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
            default: CRAB_ERROR("Apron: unreachable");
          }
          set(x, xi);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          if (op == OP_SDIV){
            apply(OP_DIVISION, x, y, z);
          }
          else {
            // Convert to intervals and perform the operation
            interval_t yi = operator[](y);
            interval_t zi = operator[](z);
            interval_t xi = interval_t::top ();
            
            switch (op) {
              case OP_UDIV: xi = yi.UDiv(zi); break;
              case OP_SREM: xi = yi.SRem(zi); break;
              case OP_UREM: xi = yi.URem(zi); break;
              default: CRAB_ERROR("Apron: unreachable");
            }
            set(x, xi);
          }
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          if (op == OP_SDIV){
            apply(OP_DIVISION, x, y, k);
          }
          else {
            // Convert to intervals and perform the operation
            interval_t yi = operator[](y);
            interval_t zi(k);
            interval_t xi = interval_t::top ();
            switch (op) {
              case OP_UDIV: xi = yi.UDiv(zi); break;
              case OP_SREM: xi = yi.SRem(zi); break;
              case OP_UREM: xi = yi.URem(zi); break;
              default: CRAB_ERROR("Apron: unreachable");
            }
            set(x, xi);
          }
        }
        
        linear_constraint_system_t to_linear_constraint_system (){
          linear_constraint_system_t csts;
          if(is_bottom ())  {
            csts += F ();
          }
          else if(is_top ()) {
            csts += T ();
          }
          else {
            ap_lincons0_array_t lcons_arr = ap_abstract0_to_lincons_array (get_man (), &*m_apstate);
            for (unsigned i=0 ; i < lcons_arr.size; i++)
              csts += tconst2const (lcons_arr.p[i]);
          }
          return csts;
        }
        
        void write(ostream& o) {

          ap_abstract0_canonicalize (get_man (), &*m_apstate);

          if(is_bottom()){
            o << "_|_";
            return;
          }
          else if (is_top()){
            o << "{}";
            return;
          }
          else {
            linear_constraint_system_t inv = to_linear_constraint_system ();
            o << inv;
          }
        }          

        const char* getDomainName () const {
          if (ApronDom == INT) 
            return "Apron Intervals"; 
          else if (ApronDom == OCT) 
            return "Apron Octagon"; 
          else if (ApronDom == PK) 
            return "Apron NewPolka";
          else 
            CRAB_ERROR("Unknown apron domain");
        }
      }; 
   
     template<typename N, typename V, apron_domain_id_t D>
     ap_manager_t* apron_domain<N,V,D>::m_apman = nullptr;

     template<typename N, typename V, apron_domain_id_t D>
     typename apron_domain<N,V,D>::var_map_ptr apron_domain<N,V,D>::m_var_map = nullptr;

   } // namespace domains

   namespace domain_traits {
       using namespace domains;

       // template <typename Number, typename VariableName>
       // void expand (DBM<Number,VariableName>& inv, 
       //              VariableName x, VariableName new_x) {
       //   inv.expand (x, new_x);
       // }
    
       // template <typename Number, typename VariableName>
       // void normalize (DBM<Number,VariableName>& inv) {
       //   inv.normalize ();
       // }
    
       // template <typename Number, typename VariableName, typename Iterator >
       // void forget (DBM<Number,VariableName>& inv, Iterator it, Iterator end) {
       //   inv.forget (it, end);
       // }

       // template <typename Number, typename VariableName, typename Iterator >
       // void project (DBM<Number,VariableName>& inv, Iterator it, Iterator end) {
       //   inv.project (it, end);
       // }
   } // namespace domain_traits
}// namespace crab
#endif 
#endif 
