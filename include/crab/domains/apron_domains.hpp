#ifndef APRON_DOMAINS_HPP
#define APRON_DOMAINS_HPP

/// Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/domains/numerical_domains_api.hpp>

#include "boost/range/algorithm/set_algorithm.hpp"

using namespace boost;
using namespace ikos;

namespace crab {
   namespace domains {
      typedef enum { APRON_INT, APRON_OCT, APRON_OPT_OCT, APRON_PK } apron_domain_id_t;
   }
}

#ifndef HAVE_APRON
/*
 * Dummy implementation if Apron not found 
 */
#define APRON_NOT_FOUND "No Apron. Run cmake with -DUSE_APRON=ON"

namespace crab {
   namespace domains {
      template<typename Number, typename VariableName, apron_domain_id_t ApronDom, unsigned Dims = 50u>
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

        apron_domain(): ikos::writeable() { }    

        static apron_domain_t top() { CRAB_ERROR (APRON_NOT_FOUND); }

        static apron_domain_t bottom() { CRAB_ERROR (APRON_NOT_FOUND); }

        apron_domain (const apron_domain_t& other): 
            ikos::writeable() { }
        
        bool is_bottom() { CRAB_ERROR (APRON_NOT_FOUND); }

        bool is_top() { CRAB_ERROR (APRON_NOT_FOUND); }

        bool operator<=(apron_domain_t other) { CRAB_ERROR (APRON_NOT_FOUND); }

        void operator|=(apron_domain_t other)
        { CRAB_ERROR (APRON_NOT_FOUND); }
        
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

      template<typename Number, typename VariableName, apron_domain_id_t ApronDom, unsigned Dims = 50u>
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
        typedef bound <Number> bound_t;
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
            else if (ApronDom == APRON_OPT_OCT)
              m_apman = opt_oct_manager_alloc ();
            else if (ApronDom == APRON_PK)
              m_apman = pk_manager_alloc (false);
            else
              CRAB_ERROR("ERROR: unknown apron domain");

            //cout << "Create apron manager: " << ap_manager_get_library (m_apman) << "\n";
          }
          return m_apman;
        }

        size_t get_dims (ap_state_ptr s) const {
          ap_dimension_t dims = _ap_abstract0_dimension (&*s);
          return dims.intdim;
        }

        size_t get_dims () const { return get_dims (m_apstate); }

        // If v is in the map then it maps v to a dimension, otherwise null
        boost::optional<ap_dim_t> get_var_dim (const var_bimap_t& m, VariableName v) const {
          auto it = m.left.find (v);
          if (it != m.left.end ())
            return it->second;
          else
            return boost::optional<ap_dim_t> ();
        }

        boost::optional<ap_dim_t> get_var_dim (VariableName v) const {
          return get_var_dim (*m_var_map, v);
        }


        ap_dim_t map_insert (var_bimap_t& m, VariableName v){
          unsigned dim = m.size ();
          // if (dim >= get_dims ())
          //   CRAB_ERROR ("Apron needs more dimensions!\n");
          m.insert (binding_t (v, dim));
          return dim;
        }

        // If v is in m then it maps v to a dimension, otherwise it
        // inserts v in m and returns a fresh dimension.
        ap_dim_t get_var_dim_insert (var_bimap_t& m, VariableName v)  {
          if (auto dim = get_var_dim (m, v))
            return *dim;
          else
            return map_insert (m, v);
        }

        ap_dim_t get_var_dim_insert (VariableName v) {
          return get_var_dim_insert (*m_var_map, v);
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

        void remove_dimensions (ap_state_ptr& s, vector<ap_dim_t> dims) const {
#if 1
          if ( (get_dims (s) - dims.size ()) <= Dims) return;

          // --- sort dimensions in descending order to avoid
          //     unnecessary shuffling in case remove_dimensions does
          //     not sort.
          std::sort(dims.begin(), dims.end(), std::greater<int>());

          ap_dimchange_t* dimchange =  ap_dimchange_alloc (dims.size (), 0);
          for (unsigned i=0; i<dims.size () ; i++) {
            // remove dimension dims[i] and shift to the left all the
            // dimensions greater than dims[i]
            dimchange->dim[i] = dims[i]; 
          }

          s = apPtr (get_man (), 
                     ap_abstract0_remove_dimensions(get_man (), false, 
                                                    &*s, dimchange));
          ap_dimchange_free (dimchange);
#endif 
        }

        bool enlarge (const size_t size, ap_state_ptr& state) {
#if 1
          const double load_factor = 0.75;
          unsigned threshold = (unsigned) ((double) get_dims (state) * load_factor);
          const double grow_factor = 1.5;
          if (size >= threshold) {
            // cout << "# Used dims= " << size  << "\n";
            // cout << "Threshold=" << threshold << "\n";
            // cout << "Enlarging from " << get_dims (state) << " to " 
            //      << get_dims (state) + (Dims * grow_factor) << "\n";
            unsigned new_dims = (size - threshold) + (unsigned ) ((double) Dims * grow_factor);
            add_dimensions (state, new_dims);
            return true;
          }
#endif 
          return false;
        }

        // unused 
        bool is_unconstrained_var (const ap_state_ptr s, VariableName v) const {
          if (auto dim = get_var_dim (v))
            return ap_abstract0_is_dimension_unconstrained (get_man (), &*s, *dim);
          else
            return false;
        }

        bool check_perm (ap_dimperm_t* perm, size_t size){
#if 0
          // it does not check injectivity
          if (perm->size != size) return false;
          for (unsigned i=0; i<perm->size; i++){
            if (perm->dim[i]>=size){
              return false;
            }
          }
#endif 
          return true;

        }

        // var_map_ptr merge_var_map (var_bimap_t m_x, ap_state_ptr& s_x,
        //                            var_bimap_t m_y, ap_state_ptr& s_y) {

        //   // -- ensure both states have the same number of dimensions,
        //   // -- otherwise it will crash
        //   size_t max_dim = std::max (get_dims (s_x), get_dims (s_y));
        //   add_dimensions (s_x, max_dim - get_dims(s_x));
        //   add_dimensions (s_y, max_dim - get_dims(s_y));

        //   // Add variables until occur in both maps
        //   for (auto const& px: m_x.left)
        //     get_var_dim_insert (m_y, px.first);
        //   for (auto const& py: m_y.left) 
        //     get_var_dim_insert (m_x, py.first);
          
        //   if (get_dims (s_x) != get_dims (s_y))
        //     CRAB_ERROR ("Merging two apron states with different dimensions",
        //                 get_dims (s_x), " != ", get_dims (s_y));

        //   size_t dims = get_dims (s_x);
          
        //   // -- initialize the permutations maps
        //   ap_dimperm_t* perm_x = ap_dimperm_alloc (dims);
        //   ap_dimperm_t* perm_y = ap_dimperm_alloc (dims);
          

        //   for (unsigned i=0; i < dims; i++){
        //     // Use AP_DIM_MAX as a dont-care symbol
        //     perm_x->dim [i] = AP_DIM_MAX; 
        //     perm_y->dim [i] = AP_DIM_MAX; 
        //   }
          
        //   // -- create a common map
        //   var_map_ptr res = var_map_ptr (new var_bimap_t ());
        //   for (auto const& px: m_x.left) {
        //     ap_dim_t i = get_var_dim_insert (*res, px.first);
        //     perm_x->dim [px.second] = i;
        //   }
        //   for (auto const& py: m_y.left) {
        //     ap_dim_t i = get_var_dim_insert (*res, py.first);
        //     perm_y->dim [py.second] = i;
        //   }

        //   unsigned unwanted_dim_x = res->size ();
        //   unsigned unwanted_dim_y = res->size ();
        //   for (unsigned i=0; i < dims; i++) {
        //     if (perm_x->dim [i] == AP_DIM_MAX)
        //       perm_x->dim [i] = unwanted_dim_x++;
        //     if (perm_y->dim [i] == AP_DIM_MAX)
        //       perm_y->dim [i] = unwanted_dim_y++;
        //   }

        //   if (!check_perm(perm_x, dims))
        //     CRAB_ERROR ("Permutation x is not well formed!");

        //   if (!check_perm(perm_y, dims))
        //     CRAB_ERROR ("Permutation y is not well formed!");

        //   // // s_x and s_y can have disjoint set of variables
        //   // if (res->left.size () >= dims) {
        //   //   enlarge (res->left.size (), s_x);
        //   //   enlarge (res->left.size (), s_y);
        //   // }

        //   #if 0        
        //   cout << "Permutations X: \n";
        //   ap_dimperm_fprint(stdout, perm_x);          
        //   cout << "Permutations Y: \n";
        //   ap_dimperm_fprint(stdout, perm_y);          
        //   #endif 

        //   // apply the permutations
        //   s_x = apPtr (get_man (), 
        //                ap_abstract0_permute_dimensions(get_man (), false, &*s_x, perm_x));
        //   s_y = apPtr (get_man (), 
        //                ap_abstract0_permute_dimensions(get_man (), false, &*s_y, perm_y));

        //   ap_dimperm_free (perm_x);
        //   ap_dimperm_free (perm_y);

        //   // remove unused dimensions

        //   vector<ap_dim_t> rem_dims;
        //   for (unsigned i=0; i< dims; i++) {
        //     if (!has_var_name (*res, i))
        //       rem_dims.push_back (i);
        //   }            

        //   remove_dimensions (s_x, rem_dims);
        //   remove_dimensions (s_y, rem_dims);

        //   return res;
        // }

        var_map_ptr merge_var_map (const var_bimap_t& m_x, ap_state_ptr& s_x,
                                   const var_bimap_t& m_y, ap_state_ptr& s_y) {
          
          // -- ensure both states have the same number of dimensions,
          // -- otherwise it will crash
          size_t max_dim = std::max (get_dims (s_x), get_dims (s_y));
          add_dimensions (s_x, max_dim - get_dims(s_x));
          add_dimensions (s_y, max_dim - get_dims(s_y));

          assert (get_dims (s_x) == get_dims (s_y));

          // -- collect all vars from the two maps
          std::set<VariableName> vars;
          for (auto const& px: m_x.left) {  
            vars.insert (px.first);
          }
          for (auto const& py: m_y.left) {
            vars.insert (py.first);
          }
            
          // -- create a fresh map
          var_map_ptr res = var_map_ptr (new var_bimap_t ());
          for (auto v: vars) {
            ap_dim_t dim = res->size ();
            res->insert (binding_t (v, dim));
          }

          // s_x and s_y can have disjoint set of variables
          if (res->left.size () >= get_dims (s_x)) {
            enlarge (res->left.size (), s_x);
            enlarge (res->left.size (), s_y);
          }
          
          // build the permutations maps
          ap_dimperm_t* perm_x = ap_dimperm_alloc (get_dims (s_x));
          ap_dimperm_t* perm_y = ap_dimperm_alloc (get_dims (s_x));
          char * xmap1 = (char *)calloc(get_dims (s_x),sizeof(char));
          char * xmap2 = (char *)calloc(get_dims (s_x),sizeof(char));
          for (auto const &px: m_x.left) {
            ap_dim_t ind = res->left.at (px.first);
            perm_x->dim [px.second] = ind;
            // This sets 1 if the index that has been assigned
            xmap1[px.second] = 1;
            // This sets 1 if the value has been assigned
            xmap2[ind] = 1;
            //perm_x->dim [ind] = tmp; 
          }
          ap_dim_t i, counter = 0;
          for(i=0; i < get_dims (s_x); i++){
            // If the index has beena assigned, skip
            if(xmap1[i]){
              continue;
            }
            // Find the next available element that has not been assigned
            while(xmap2[counter]){
			counter++;
            }
            perm_x->dim[i] = counter;
            counter++;
          }
          free(xmap1);
          free(xmap2);
          char * ymap1 = (char *)calloc(get_dims (s_x),sizeof(char));
          char * ymap2 = (char *)calloc(get_dims (s_x),sizeof(char));
          for (auto const &py: m_y.left) {
            //if (!is_unconstrained_var (s_y, py.first))
            ap_dim_t ind = res->left.at (py.first);
            perm_y->dim [py.second] = ind;
            //perm_y->dim [ind] = tmp;
            ymap1[py.second] = 1;
            ymap2[ind] = 1; 
          }
          counter = 0;
          for(i=0; i < get_dims (s_x); i++){
            if(ymap1[i]){
              continue;
            }
            while(ymap2[counter]){
              counter++;
            }
            perm_y->dim[i] = counter;
            counter++;
          }
          free(ymap1);
          free(ymap2);

          #if 0          
          cout << "Permutations \n";
          ap_dimperm_fprint(stdout, perm_x);          
          cout << "Permutations \n";
          ap_dimperm_fprint(stdout, perm_y);          
          #endif 

          // apply the permutations
          s_x = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_x, perm_x));
          s_y = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_y, perm_y));

          ap_dimperm_free (perm_x);
          ap_dimperm_free (perm_y);

          // remove unused dimensions
          vector<ap_dim_t> rem_dims;
          for (unsigned i=0; i< get_dims (s_x); i++) {
            if (!has_var_name (*res, i))
              rem_dims.push_back (i);
          }            

          remove_dimensions (s_x, rem_dims);
          remove_dimensions (s_y, rem_dims);

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

        inline ap_texpr0_t* var2texpr (VariableName v) { 
          return ap_texpr0_dim (get_var_dim_insert (v));
        }

        inline ap_texpr0_t* num2texpr (Number i) const {  
          mpq_class n = ((mpz_class) i); 
          return ap_texpr0_cst_scalar_mpq (n.get_mpq_t ());
        }

        inline ap_texpr0_t* intv2texpr (Number a, Number b) const { 
          mpq_class n1 = ((mpz_class) a); 
          mpq_class n2 = ((mpz_class) b); 
          return ap_texpr0_cst_interval_mpq (n1.get_mpq_t (), n2.get_mpq_t ());
        }
        
        inline ap_texpr0_t* expr2texpr (linear_expression_t e)  {
          Number cst = e.constant ();
          ap_texpr0_t* res = num2texpr (cst);
          for (auto p: e) {
            ap_texpr0_t* term = MUL (num2texpr (p.first), var2texpr (p.second.name ()));
            res = ADD (res, term); 
          }
          return res;
        }

        inline ap_tcons0_t const2tconst (linear_constraint_t cst)  {
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
          if (scalar->discr == AP_SCALAR_DOUBLE) { // elina uses double
            return Number ((long) scalar->val.dbl);
          }
          else if (scalar->discr == AP_SCALAR_MPQ) {
            mpq_class n (scalar->val.mpq);
            return Number ((mpz_class) n);
          }
          else
            CRAB_ERROR ("ERROR: apron translation only covers double or mpq scalars");
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
                                      linear_expression_t (Number(1)));          
        }

        inline linear_constraint_t F () const {
          return linear_constraint_t (linear_expression_t (Number(1)) == 
                                      linear_expression_t (Number(0)));          
        }

        void dump (const var_bimap_t& m, ap_state_ptr apstate ) {  
          vector<char*> names;
          for (unsigned i=0; i < get_dims (apstate) ; i++){
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

       public:
        void print_stats () { ap_abstract0_fprint (stdout, get_man (), &*m_apstate, NULL); }

       private:
        void dump_vars () {
          #if 0
          cout << "vars={";
          for (auto const& p: m_var_map->left)
            cout << p.first <<";";
          cout << "}\n";
          #endif 
        }

        apron_domain (ap_state_ptr apState, var_map_ptr varMap): 
            ikos::writeable (), 
            m_apstate (apState), 
            m_var_map (varMap) { 
          
          if (ap_abstract0_is_top (get_man(), &*m_apstate) ||
              ap_abstract0_is_bottom (get_man(), &*m_apstate)) {
            m_var_map = var_map_ptr (new var_bimap_t ());
          }

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

        // XXX
        apron_domain (const apron_domain_t& o): 
            ikos::writeable(), 
            m_apstate (apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)))),
            m_var_map (var_map_ptr (new var_bimap_t (*o.m_var_map))) 
            //m_apstate (o.m_apstate),
            //m_var_map (o.m_var_map)
        {  }
        
        // XXX        
        apron_domain_t operator=(const apron_domain_t& o) {
          //m_apstate  = apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)));
          //std::swap (*m_var_map, *o.m_var_map);
          m_apstate = o.m_apstate;
          m_var_map = o.m_var_map;
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
            // XXX JN: copy this->m_apstate because merge_var_map will modify it
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            merge_var_map (*m_var_map, x, *(o.m_var_map), o.m_apstate);
            return ap_abstract0_is_leq (get_man(), &*x, &*o.m_apstate);
          }
        }

        void operator|=(apron_domain_t o) {
          if (is_bottom() || o.is_top ())
            *this = o;
          else if (is_top () || o.is_bottom())
            return ;
          else {
              m_var_map = merge_var_map (*m_var_map, m_apstate, *(o.m_var_map), o.m_apstate);
              m_apstate = apPtr (get_man(), ap_abstract0_join (get_man(), false, 
                                                               &*m_apstate, &*o.m_apstate));
          }
        }
        
        apron_domain_t operator|(apron_domain_t o) {
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom() || o.is_top ())
            return o;
          else if (is_top () || o.is_bottom())
            return *this;
          else {

            // XXX JN: Copy this->m_apstate because merge_var_map will
            // modify it and the join is a function so this should not
            // be modified.
            
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_join (get_man(), false, 
                                                             &*x, &*o.m_apstate)), m);
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

            // XXX JN: Copy this->m_apstate because merge_var_map will
            // modify it and the meet is a function so this should not
            // be modified.
            
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_meet (get_man(), false, 
                                                             &*x, &*o.m_apstate)), m);
          }
        }        
        
        apron_domain_t operator||(apron_domain_t o) {
          // cover all trivial cases to avoid permutating dimensions
          if (is_bottom())
            return o;
          else if (o.is_bottom())
            return *this;
          else {

            // XXX JN: Copy this->m_apstate because merge_var_map will
            // modify it and the widening is a function so this should not
            // be modified.
            
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_ptr  m = merge_var_map (*m_var_map, x, *(o.m_var_map), o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_widening (get_man(), 
                                                                 &*x, &*o.m_apstate)), m);
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
          for (auto v: vars)  {
            if (auto dim = get_var_dim (v)) {
              dims.push_back (*dim);
              m_var_map->left.erase (v);
            }
          }

          if (dims.empty ()) return;

          m_apstate = apPtr (get_man (), 
                             ap_abstract0_forget_array (get_man (), 
                                                        false, 
                                                        &*m_apstate, 
                                                        &dims[0], dims.size(), 
                                                        false));
          remove_dimensions (m_apstate, dims);
        }

        // remove all variables except vars
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
          if (auto dim = get_var_dim (var)) {
            dims.push_back (*dim);
            m_var_map->left.erase (var);
            m_apstate = apPtr (get_man (), 
                               ap_abstract0_forget_array (get_man (), 
                                                          false, 
                                                          &*m_apstate, 
                                                          &dims[0], dims.size(), 
                                                          false));
            remove_dimensions (m_apstate, dims);
          }
        }

        interval_t operator[](VariableName v) {

          if (is_bottom ()) 
            return interval_t::bottom ();

          if (auto dim = get_var_dim (v)) {

            ap_interval_t* intv = ap_abstract0_bound_dimension (get_man (),
                                                                &*m_apstate, 
                                                                *dim);
            if (ap_interval_is_top (intv))
              return interval_t::top ();

            ap_scalar_t* lb = intv->inf;
            ap_scalar_t* ub = intv->sup;
            
            if (lb->discr == AP_SCALAR_DOUBLE && ub->discr == AP_SCALAR_DOUBLE) { 

              if (ap_scalar_infty(lb) == -1) {     // [-oo, k]
                return interval_t (bound_t::minus_infinity (),
                                   Number ((long) ub->val.dbl));

              }
              else if (ap_scalar_infty(ub) == 1) { // [k, +oo]
                return interval_t (Number ((long) lb->val.dbl),
                                   bound_t::plus_infinity ());
              }
              else { 
                assert (ap_scalar_infty(lb) == 0); // lb is finite
                assert (ap_scalar_infty(ub) == 0); // ub is finite
                return interval_t (Number ((long) lb->val.dbl),
                                   Number ((long) ub->val.dbl));
              }

            }
            else if (lb->discr == AP_SCALAR_MPQ && ub->discr == AP_SCALAR_MPQ ) {

              if (ap_scalar_infty(lb) == -1) {     // [-oo, k]
                mpq_class sup (ub->val.mpq);
                return interval_t (bound_t::minus_infinity (),
                                   Number ((mpz_class) sup));

              }
              else if (ap_scalar_infty(ub) == 1) { // [k, +oo]
                mpq_class inf (lb->val.mpq);
                return interval_t (Number ((mpz_class) inf),
                                   bound_t::plus_infinity ());
              }
              else {
                assert (ap_scalar_infty(lb) == 0); // lb is finite
                assert (ap_scalar_infty(ub) == 0); // ub is finite

                mpq_class inf (lb->val.mpq);
                mpq_class sup (ub->val.mpq);
                return interval_t (Number ((mpz_class) inf), 
                                   Number ((mpz_class) sup));
              }

            }
            else 
              CRAB_ERROR ("ERROR: apron translation only covers double or mpq scalars");
          }
          else 
            return interval_t::top ();
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

          if (enlarge (m_var_map->left.size(), m_apstate))
            dump_vars ();

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
          CRAB_DEBUG("Added ",csts, " --> ", *this);
        }
       
        void assign (VariableName x, linear_expression_t e) {
          if(is_bottom()) return;

          if (enlarge (m_var_map->left.size(), m_apstate))
            dump_vars();
                      
          ap_texpr0_t* t = expr2texpr (e);
          assert (t);
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_assign_texpr(get_man (), false, 
                                                       &*m_apstate, 
                                                       get_var_dim_insert (x),t, 
                                                       NULL));

          ap_texpr0_free (t);
          CRAB_DEBUG(x, ":=", e , " --> ", *this);
        }
          
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          if(is_bottom()) return;

          if (enlarge (m_var_map->left.size(), m_apstate))
            dump_vars ();

          ap_texpr0_t* a = var2texpr (y);
          ap_texpr0_t* b = num2texpr (z);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("ERROR apron: operation not supported");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim_insert (x), 
                                                                   res, 
                                                                   NULL));

          
          ap_texpr0_free (res);
          CRAB_DEBUG(x, ":=", y, op, z, " --> ", *this);
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          if(is_bottom()) return;

          if (enlarge (m_var_map->left.size(), m_apstate))
            dump_vars ();

          ap_texpr0_t* a = var2texpr (y);
          ap_texpr0_t* b = var2texpr (z);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("ERROR apron: operation not supported");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim_insert (x), 
                                                                   res, 
                                                                   NULL));

          ap_texpr0_free (res);
          CRAB_DEBUG(x, ":=", y, op, z, " --> ", *this);
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          if(is_bottom()) return;

          if (enlarge (m_var_map->left.size(), m_apstate))
            dump_vars ();

          ap_texpr0_t* a = var2texpr (x);
          ap_texpr0_t* b = num2texpr (k);
          ap_texpr0_t* res = nullptr;
          
          switch (op){
            case OP_ADDITION: res = ADD (a, b); break;
            case OP_SUBTRACTION: res = SUB (a, b); break;
            case OP_MULTIPLICATION: res = MUL (a, b); break;
            case OP_DIVISION: res = DIV (a, b); break;
            default: CRAB_ERROR("ERROR apron: operation not supported");
          }
          assert (res);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   get_var_dim_insert (x), 
                                                                   res, 
                                                                   NULL));

          ap_texpr0_free (res);
          CRAB_DEBUG(x, ":=", x, op, k, " --> ", *this);
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
            default: CRAB_ERROR("ERROR apron: operation not supported");
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
            default: CRAB_ERROR("ERROR apron: operation not supported");
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
              default: CRAB_ERROR("ERROR apron: operation not supported");
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
              default: CRAB_ERROR("ERROR apron: operation not supported");
            }
            set(x, xi);
          }
        }
        
        linear_constraint_system_t to_linear_constraint_system () {
          linear_constraint_system_t csts;
          if(is_bottom ())  {
            csts += F ();
          }
          else if(is_top ()) {
            csts += T ();
          }
          else {
            normalize ();

            ap_lincons0_array_t lcons_arr = ap_abstract0_to_lincons_array (get_man (), &*m_apstate);
            for (unsigned i=0 ; i < lcons_arr.size; i++)
              csts += tconst2const (lcons_arr.p[i]);

            ap_lincons0_array_clear (&lcons_arr);
          }
          return csts;
        }
        
        void expand (VariableName x, VariableName dup) {
          *this -= dup;
          // --- increases number of dimensions by one
          m_apstate = apPtr (get_man(),
                             ap_abstract0_expand(get_man (), false, &* m_apstate, 
                                                 get_var_dim_insert (x), 1));

          // --- the additional dimension is put at the end of integer
          //     dimensions.
          m_var_map->insert (binding_t (dup, get_dims () - 1));            
        }

        void normalize () {
          ap_abstract0_canonicalize (get_man (), &*m_apstate);
        }

        void write(ostream& o) {
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
          else if (ApronDom == APRON_OPT_OCT) 
            return "Apron Optimized Octagon"; 
          else if (ApronDom == APRON_PK) 
            return "Apron NewPolka";
          else 
            CRAB_ERROR("Unknown apron domain");
        }
      }; 

      // --- global datastructures

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

      // temporary
      template <typename Number, typename VariableName, apron_domain_id_t ApronDom>
      void print_stats (apron_domain<Number,VariableName, ApronDom>& inv) {
        inv.print_stats ();
      }


   } // namespace domain_traits
}// namespace crab
#endif 
#endif 
