#ifndef APRON_DOMAINS_HPP
#define APRON_DOMAINS_HPP

/// Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>

using namespace boost;
using namespace ikos;

namespace crab {
   namespace domains {
      typedef enum { APRON_INT, APRON_OCT, APRON_PK } apron_domain_id_t;
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

     using namespace apron;

      // TODO: resize (increase/decrease number of dimensions)
      template<typename Number, typename VariableName, apron_domain_id_t ApronDom, unsigned Dims = 50>
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
        typedef apron_domain <Number, VariableName, ApronDom, Dims> apron_domain_t;
        typedef interval <Number> interval_t;

       private:
        typedef interval_domain <Number, VariableName> interval_domain_t;
        typedef boost::bimap< VariableName , ap_dim_t > var_bimap_t;
        typedef boost::shared_ptr<var_bimap_t> var_map_ptr;
        typedef typename var_bimap_t::value_type binding_t;

        static ap_manager_t* m_apman;
        
        ap_state_ptr m_apstate;
        var_map_ptr m_var_map;

        static ap_manager_t* get_man () {
          if (!m_apman) {
            if (ApronDom == APRON_INT)
              m_apman = box_manager_alloc ();
            else if (ApronDom == APRON_OCT)
              m_apman = oct_manager_alloc ();
            else if (ApronDom == APRON_PK)
              m_apman = pk_manager_alloc (false);
            else
              CRAB_ERROR("Unknown apron domain");
          }
          return m_apman;
        }

        size_t get_dims (ap_state_ptr s) const {
          ap_dimension_t dims = _ap_abstract0_dimension (&*s);
          return dims.intdim;
        }

        size_t get_dims () const { return get_dims (m_apstate); }

        ap_dim_t get_var_dim (VariableName v) const {
          auto it = m_var_map->left.find (v);
          if (it != m_var_map->left.end ())
            return it->second;
          else {
            unsigned dim = m_var_map->size ();
            if (dim > get_dims ())
              CRAB_ERROR ("TODO: apron needs to resize\n");
            m_var_map->insert (binding_t (v, dim));
            return dim;
          }
        }

        bool has_var_name (const var_bimap_t& m, ap_dim_t i) const {
          return m.right.find (i) != m.right.end ();
        }

        VariableName get_var_name (const var_bimap_t& m, ap_dim_t i) const {
          auto it = m.right.find (i);
          if (it != m.right.end ())
            return it->second;            
          CRAB_ERROR ("Apron dimension ", i, " is not used!");
        }

        bool has_var_name (ap_dim_t i) const {
          return has_var_name (*m_var_map, i);
        }

        VariableName get_var_name (ap_dim_t i) const {
          return get_var_name (*m_var_map, i);
        }

        void add_dimensions (ap_state_ptr& s, size_t dims) const {
          if (dims <= 0) return;

          ap_dimchange_t* dimchange =  ap_dimchange_alloc (dims, 0);
          for (unsigned i=0 ; i<dims ; i++)
            dimchange->dim[i] = get_dims (s); // add dimension at the end

          s = apPtr (get_man (), 
                     ap_abstract0_add_dimensions(get_man (), false, 
                                                 &*s, dimchange, false));
          ap_dimchange_free (dimchange);
          
        }

        var_map_ptr merge_var_map (const var_bimap_t& m_x, ap_state_ptr& s_x,
                                   const var_bimap_t& m_y, ap_state_ptr& s_y) {

          size_t max_dim = std::max (get_dims (s_x), get_dims (s_y));
          add_dimensions (s_x, max_dim - get_dims(s_x));
          add_dimensions (s_y, max_dim - get_dims(s_y));

          assert (get_dims (s_x) == get_dims (s_y));
          
          // -- collect all vars from the two maps
          std::set<VariableName> vars;
          for (auto px: m_x.left) vars.insert (px.first);
          for (auto py: m_y.left) vars.insert (py.first);
            

          // -- create a fresh map
          var_map_ptr res = var_map_ptr (new var_bimap_t ());
          unsigned dim = 0;
          for (auto v: vars) {
            res->left.insert (make_pair (v, dim));
            ++dim;
          }

          if (res->size () > get_dims ())
            CRAB_ERROR ("TODO: apron needs to resize!");
          
          // build the permutations maps
          ap_dimperm_t* perm_x = ap_dimperm_alloc (get_dims ());
          ap_dimperm_t* perm_y = ap_dimperm_alloc (get_dims ());
          ap_dimperm_set_id (perm_x); // identity permutation
          ap_dimperm_set_id (perm_y); // identity permutation

          for (auto const &p: m_x.left)
            perm_x->dim [p.second] = res->left.at (p.first); 
          for (auto const &p: m_y.left)
            perm_y->dim [p.second] = res->left.at (p.first);

          // cout << "Before permutations: \n";
          // dump (m_x, s_x);
          // dump (m_y, s_y);

          // cout << "Permutations \n";
          // ap_dimperm_fprint(stdout, perm_x);          
          // cout << "Permutations \n";
          // ap_dimperm_fprint(stdout, perm_y);          

          // apply the permutations
          s_x = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_x, perm_x));
          s_y = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_y, perm_y));

          // cout << "After permutations: \n";
          // dump (*res, s_x);
          // dump (*res, s_y);

          ap_dimperm_free (perm_x);
          ap_dimperm_free (perm_y);

          return res;
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
          // LEAKING!
          return ap_texpr0_dim( get_var_dim (v));
        }

        inline ap_texpr0_t* num2texpr (Number i) const {  
          // LEAKING!
          return ap_texpr0_cst_scalar_int ((int) i); 
        }

        inline ap_texpr0_t* intv2texpr (Number a, Number b) const { 
          // LEAKING!
          return ap_texpr0_cst_interval_int ((int) a,(int) b); 
        }
        
        inline ap_texpr0_t* expr2texpr (linear_expression_t e) const {
          // LEAKING!
          Number cst = e.constant ();
          ap_texpr0_t* res = num2texpr (cst);
          for (auto p: e) {
            ap_texpr0_t* term = MUL (num2texpr (p.first), var2texpr (p.second.name ()));
            res = ADD (res, term); 
          }
          return res;
        }

        inline ap_tcons0_t const2tconst (linear_constraint_t cst) const {
          // LEAKING!
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

          // add possible constant
          ap_coeff_t* cst = ap_linexpr0_cstref (linexp);
          if (!ap_coeff_zero (cst)) 
            e = e + coeff2Num(cst);
          linear_constraint_t res;
          switch (cons.constyp) {
            case AP_CONS_EQ:
              // e == k 
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::EQUALITY);
              break;
            case AP_CONS_SUPEQ:
              // e >= k 
              e = -e;
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY);
              break;
            case AP_CONS_SUP:
              // e > k 
              e = -e + 1;
              res =  linear_constraint_t (e, linear_constraint_t::kind_t::INEQUALITY);
              break;
            case AP_CONS_EQMOD:
              res = T ();
              break;
            case AP_CONS_DISEQ:
              // e != k 
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

        void dump (const var_bimap_t& m, ap_state_ptr apstate ) {  
          vector<char*> names;
          for (unsigned i=0; i < get_dims () ; i++){
            string varname;
            if (has_var_name (m, i))
              varname = get_var_name (m, i).str ();
            else // unused dimension
              varname = string ("x") + std::to_string (i);
            char* name = new char [varname.length () + 1];
            strcpy (name, varname.c_str ());
            names.push_back (name);
          }
          ap_abstract0_fprint (stdout, get_man (), &*apstate, &names[0]);
          for (auto n : names) { delete n; }
        }

        void dump () { dump (*m_var_map, m_apstate); }

        apron_domain (ap_state_ptr apState, var_map_ptr varMap): 
            ikos::writeable (), 
            m_apstate (apState), 
            m_var_map (varMap) { }

        bool is_var_unconstrained (VariableName x) const {
          return ap_abstract0_is_dimension_unconstrained (get_man (), 
                                                          &*m_apstate,
                                                          get_var_dim (x));
        }

       public:

        apron_domain (): 
            ikos::writeable (),
            m_apstate (apPtr (get_man(), ap_abstract0_top (get_man(), Dims, 0))),
            m_var_map (var_map_ptr (new var_bimap_t ())) { }

        ~apron_domain () { }

        static apron_domain_t top() { 
          return apron_domain_t (apPtr (get_man(), ap_abstract0_top (get_man(), Dims, 0)),
                                 var_map_ptr (new var_bimap_t ()));
        }

        static apron_domain_t bottom() { 
          return apron_domain_t (apPtr (get_man(), ap_abstract0_bottom (get_man(), Dims, 0)),
                                 var_map_ptr (new var_bimap_t ()));
        }

        apron_domain (const apron_domain_t& other): 
            ikos::writeable(), 
            m_apstate (other.m_apstate), 
            m_var_map (other.m_var_map) {  }
        
        
        apron_domain_t operator=(const apron_domain_t& o) {
          if (this != &o) {
            m_apstate = o.m_apstate;
            m_var_map = o.m_var_map;
          }
          return *this;
        }
        
        bool is_bottom() { 
          return ap_abstract0_is_bottom (get_man(), &*m_apstate);
        }

        bool is_top() { 
          return ap_abstract0_is_top (get_man(), &*m_apstate);
        }

        bool operator<=(apron_domain_t o) { 
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom()) 
            return true;
          else if(o.is_bottom())
            return false;
          else if (o.is_top ())
            return true;
          else if (is_top () && !o.is_top ())
            return false;
          else if (is_top () && o.is_top ())
            return true;
          else { 
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            ap_state_ptr y = apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)));
            merge_var_map (*m_var_map, x, *(o.m_var_map), y);
            return ap_abstract0_is_leq (get_man(), &*x, &*y);
          }
        }
        
        apron_domain_t operator|(apron_domain_t o) {
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom() || o.is_top ())
            return o;
          else if (is_top () || o.is_bottom())
            return *this;
          else {
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            ap_state_ptr y = apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), y);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_join (get_man(), false, &*x, &*y)), m);
          }
        }        
        
        apron_domain_t operator&(apron_domain_t o) {
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom() || o.is_bottom())
            return bottom();
          else if (is_top())
            return o;
          else if (o.is_top())
            return *this;
          else{
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            ap_state_ptr y = apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), y);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_meet (get_man(), false, &*x, &*y)), m);
          }
        }        
        
        apron_domain_t operator||(apron_domain_t o) {
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom())
            return o;
          else if (o.is_bottom())
            return *this;
          else {
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            ap_state_ptr y = apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), y);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_widening (get_man(), &*x, &*y)), m);
          }
        }        
        
        apron_domain_t operator&&(apron_domain_t o) {
          // apron does not define narrowing operators.
          // It might not terminate so make sure that the fixpoint
          // runs a finite number of descending iterations
          return (*this & o);
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

        // dual of forget: remove all variables except vars
        template<typename Range>
        void project (const Range& vars) {
          if (is_bottom ()) return;
          std::set<VariableName> s1,s2,s3;
          for (auto p: m_var_map->left) s1.insert (p.first);
          s2.insert (vars.begin (), vars.end ());
          boost::set_difference (s1,s2,std::inserter (s3, s3.end ()));
          forget (s3);
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

          ap_tcons0_array_clear(&array);
          //cout << "Added " << csts << " --> "; dump (); cout << "\n";
        }
       
        void assign (VariableName x, linear_expression_t e) {
          ap_texpr0_t* t = expr2texpr (e);
          assert (t);
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_assign_texpr(get_man (), false, 
                                                       &*m_apstate, 
                                                       get_var_dim (x),t, 
                                                       NULL));

          ap_texpr0_free (t);
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

          
          ap_texpr0_free (res);
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

          ap_texpr0_free (res);
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

          ap_texpr0_free (res);
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

            ap_lincons0_array_clear (&lcons_arr);
          }
          return csts;
        }
        
        void expand (VariableName x, VariableName dup) {
          if (m_var_map->left.find (dup) != m_var_map->left.end ())
            CRAB_ERROR ("Apron tried to exxpand to an existing variable");
          
          // --- increases number of dimensions by one
          m_apstate = apPtr (get_man(),
                             ap_abstract0_expand(get_man (), false, &* m_apstate, get_var_dim (x), 1));

          // --- the additional dimension is put at the end of integer
          //     dimensions.
          m_var_map->insert (binding_t (dup, get_dims () - 1));            
        }

        void normalize () {
          ap_abstract0_canonicalize (get_man (), &*m_apstate);
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
            //dump ();
            linear_constraint_system_t inv = to_linear_constraint_system ();
            o << inv;

          }
        }          

        const char* getDomainName () const {
          if (ApronDom == APRON_INT) 
            return "Apron Intervals"; 
          else if (ApronDom == APRON_OCT) 
            return "Apron Octagon"; 
          else if (ApronDom == APRON_PK) 
            return "Apron NewPolka";
          else 
            CRAB_ERROR("Unknown apron domain");
        }
      }; 
   
      template<typename N, typename V, apron_domain_id_t D, unsigned S>
      ap_manager_t* apron_domain<N,V,D,S>::m_apman = nullptr;

   } // namespace domains

   namespace domain_traits {
      using namespace domains;

      template <typename Number, typename VariableName, apron_domain_id_t ApronDom>
      void expand (apron_domain<Number,VariableName, ApronDom>& inv,
                   VariableName x, VariableName new_x) {
        inv.expand (x, new_x);
      }
    
      template <typename Number, typename VariableName, apron_domain_id_t ApronDom>
      void normalize (apron_domain<Number,VariableName, ApronDom>& inv) {
         inv.normalize ();
      }
    
      template <typename Number, typename VariableName, apron_domain_id_t ApronDom, typename Iterator>
      void forget (apron_domain<Number,VariableName, ApronDom>& inv, Iterator it, Iterator end) {
         inv.forget (boost::make_iterator_range (it, end));
      }

      template <typename Number, typename VariableName, apron_domain_id_t ApronDom, typename Iterator>
      void project (apron_domain<Number,VariableName, ApronDom>& inv, Iterator it, Iterator end) {
        inv.project (boost::make_iterator_range (it, end));
      }

   } // namespace domain_traits
}// namespace crab
#endif 
#endif 
