#ifndef APRON_DOMAINS_HPP
#define APRON_DOMAINS_HPP

#include "boost/range/algorithm/set_algorithm.hpp"

#include <crab/config.h>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/operators_api.hpp>
#include <crab/domains/domain_traits.hpp>

using namespace boost;
using namespace ikos;

namespace crab {
   namespace domains {
      typedef enum { APRON_INT, 
                     APRON_OCT, 
                     APRON_OPT_OCT, 
                     APRON_PK } apron_domain_id_t;
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
         public division_operators< Number, VariableName >,
         public array_operators< Number, VariableName >,
         public pointer_operators< Number, VariableName > {
              
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
            ikos::writeable() { 
        }
        
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
        
        template<typename Thresholds>
        apron_domain_t widening_thresholds (apron_domain_t e, const Thresholds &ts) 
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
        
        void write(crab_os& o) 
        { CRAB_ERROR (APRON_NOT_FOUND); }
          
        static std::string getDomainName () {
          return "Dummy Apron";
        }  
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

      template<typename Number, typename VariableName, apron_domain_id_t ApronDom>
      class apron_domain_: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName >,
         public array_operators< Number, VariableName >,
         public pointer_operators< Number, VariableName > {
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;
        typedef interval <Number> interval_t;

       private:
        typedef apron_domain_ <Number, VariableName, ApronDom> apron_domain_t;
        typedef interval_domain <Number, VariableName> interval_domain_t;
        typedef bound <Number> bound_t;
        typedef boost::bimap< VariableName , ap_dim_t > var_map_t;
        typedef typename var_map_t::value_type binding_t;

        static ap_manager_t* m_apman;
        
        ap_state_ptr m_apstate; 
        var_map_t m_var_map;

        static ap_manager_t* get_man () {
          if (!m_apman) {
            switch (ApronDom) {
              case APRON_INT: m_apman = box_manager_alloc (); break;
              case APRON_OCT: m_apman = oct_manager_alloc (); break;
              case (APRON_OPT_OCT): m_apman = opt_oct_manager_alloc (); break;
              case APRON_PK: m_apman = pk_manager_alloc (false); break;
              default: CRAB_ERROR("ERROR: unknown apron domain");
            }
          }
          return m_apman;
        }

        size_t get_dims (ap_state_ptr s) const {
          ap_dimension_t dims = _ap_abstract0_dimension (&*s);
          return dims.intdim;
        }

        size_t get_dims () const { return get_dims (m_apstate); }

        // If v is in the map then it maps v to a dimension, otherwise null
        boost::optional<ap_dim_t> get_var_dim (const var_map_t& m, VariableName v) const {
          auto it = m.left.find (v);
          if (it != m.left.end ())
            return it->second;
          else
            return boost::optional<ap_dim_t> ();
        }

        boost::optional<ap_dim_t> get_var_dim (VariableName v) const {
          return get_var_dim (m_var_map, v);
        }

        ap_dim_t get_var_dim_insert (VariableName v) {
          assert (m_var_map.size () == get_dims ());
          if (auto dim = get_var_dim (v))
            return *dim;
          else {
            ap_dim_t i = m_var_map.size ();
            m_var_map.insert (binding_t (v, i));
            add_dimensions (m_apstate, 1);
            assert (m_var_map.size () == get_dims ());
            return i;
          }
        }
        
        bool has_var_name (const var_map_t& m, ap_dim_t i) const {
          return m.right.find (i) != m.right.end ();
        }

        bool has_var_name (ap_dim_t i) const {
          return has_var_name (m_var_map, i);
        }

        VariableName get_var_name (const var_map_t& m, ap_dim_t i) const {
          auto it = m.right.find (i);
          if (it != m.right.end ())
            return it->second;            
          CRAB_ERROR ("Apron dimension ", i, " is not used!");
        }

        VariableName get_var_name (ap_dim_t i) const {
          return get_var_name (m_var_map, i);
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
          if (dims.empty ()) return;

          // Apron assumption: make sure that the removing dimensions
          //                   are in ascending order.
          std::sort(dims.begin(), dims.end());

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

          #if 0
          crab::outs() << "Removed " << dims.size () << " dimensions\n";
          crab::outs() << "Size = " << get_dims (s) << "\n";
          #endif           
        }

        bool check_perm (ap_dimperm_t* perm, size_t size){
          // it does not check injectivity
          if (perm->size != size) return false;
          for (unsigned i=0; i<perm->size; i++){
            if (perm->dim[i]>=size){
              return false;
            }
          }
          return true;
        }

        var_map_t merge_var_map (const var_map_t& m_x, ap_state_ptr& s_x,
                                   const var_map_t& m_y, ap_state_ptr& s_y) {

          assert (m_x.size () == get_dims (s_x));
          assert (m_y.size () == get_dims (s_y));

          // -- collect all vars from the two maps
          std::set<VariableName> vars;
          for (auto const& px: m_x.left)
            vars.insert (px.first);
          for (auto const& py: m_y.left)
            vars.insert (py.first);

          assert (vars.size () >= get_dims (s_x));
          assert (vars.size () >= get_dims (s_y));

          add_dimensions (s_x, vars.size () - get_dims (s_x));
          add_dimensions (s_y, vars.size () - get_dims (s_y));

          assert (get_dims (s_x) == get_dims (s_y));

          // -- create a fresh map 
          var_map_t res;
          for (auto v: vars) {
            ap_dim_t i = res.size ();
            assert (i < get_dims (s_x));
            res.insert (binding_t (v, i));
          }
          
          // build the permutations maps
          ap_dimperm_t* perm_x = ap_dimperm_alloc (get_dims (s_x));
          ap_dimperm_t* perm_y = ap_dimperm_alloc (get_dims (s_x));
          char * xmap1 = (char *)calloc (get_dims (s_x), sizeof(char));
          if (!xmap1) CRAB_ERROR ("calloc does not have more available memory");
          char * xmap2 = (char *)calloc (get_dims (s_x), sizeof(char));
          if (!xmap2) CRAB_ERROR ("calloc does not have more available memory");
          for (auto const &px: m_x.left) {
            ap_dim_t ind = res.left.at (px.first);
            perm_x->dim [px.second] = ind;
            // This sets 1 if the index that has been assigned
            assert (px.second < get_dims(s_x));
            xmap1[px.second] = 1;
            // This sets 1 if the value has been assigned
            assert (ind < get_dims(s_x));
            xmap2[ind] = 1;
          }
          ap_dim_t i, counter = 0;
          for(i=0; i < get_dims (s_x); i++){
            // If the index has been assigned, skip
            if(xmap1[i]) continue;
            // Find the next available element that has not been assigned
            while(xmap2[counter])
              counter++;
            perm_x->dim[i] = counter;
            counter++;
          }
          free (xmap1);
          free (xmap2);

          char * ymap1 = (char *)calloc (get_dims (s_x), sizeof(char));
          if (!ymap1) CRAB_ERROR ("calloc does not have more available memory");
          char * ymap2 = (char *)calloc (get_dims (s_x), sizeof(char));
          if (!ymap2) CRAB_ERROR ("calloc does not have more available memory");
          for (auto const &py: m_y.left) {
            ap_dim_t ind = res.left.at (py.first);
            perm_y->dim [py.second] = ind;
            assert (py.second < get_dims(s_x));
            ymap1[py.second] = 1;
            assert (ind < get_dims(s_x));
            ymap2[ind] = 1; 
          }

          counter = 0;
          for(i=0; i < get_dims (s_x); i++){
            if(ymap1[i]) continue;
            while(ymap2[counter])
              counter++;
            perm_y->dim[i] = counter;
            counter++;
          }

          free (ymap1);
          free (ymap2);

          #if 0          
          crab::outs() << "Permutations \n";
          ap_dimperm_fprint(stdout, perm_x);          
          crab::outs() << "Permutations \n";
          ap_dimperm_fprint(stdout, perm_y);          
          #endif 

          assert (check_perm (perm_x, get_dims (s_x)));
          assert (check_perm (perm_y, get_dims (s_x)));

          // apply the permutations
          s_x = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_x, perm_x));
          s_y = apPtr (get_man (), 
                       ap_abstract0_permute_dimensions(get_man (), false, &*s_y, perm_y));

          ap_dimperm_free (perm_x);
          ap_dimperm_free (perm_y);

          assert (res.size () == get_dims (s_x));
          assert (res.size () == get_dims (s_y));

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
            
            if (!has_var_name (i)) continue; // unused dimension

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

        void dump (const var_map_t& m, ap_state_ptr apstate ) {  
          crab::outs() << "\nNumber of dimensions=" << get_dims (apstate) << "\n";
          crab::outs() << "variable map ["; 
          vector<char*> names;
          for (unsigned i=0; i < get_dims (apstate) ; i++){
            string varname;
            if (has_var_name (m, i))
              varname = get_var_name (m, i).str ();
            else // unused dimension
              varname = string ("_x") + std::to_string (i);
            crab::outs() << i << " -> " << varname << ";";
            char* name = new char [varname.length () + 1];
            strcpy (name, varname.c_str ());
            names.push_back (name);
          }
          crab::outs() << "]\n";
          ap_abstract0_fprint (stdout, get_man (), &*apstate, &names[0]);
          for (auto n : names) { delete n; }
        }

        void dump () { dump (m_var_map, m_apstate); }

       public:
        void print_stats () { ap_abstract0_fprint (stdout, get_man (), &*m_apstate, NULL); }

       private:

        apron_domain_ (ap_state_ptr apState, var_map_t varMap): 
            ikos::writeable (), 
            m_apstate (apState), 
            m_var_map (varMap) { 

          vector<ap_dim_t> dims;
          var_map_t res;
          /// XXX: we must iterate on the dimension id's to preserve
          /// order between them
          for (auto const& p: m_var_map.right) {  
            if (ap_abstract0_is_dimension_unconstrained (get_man (),
                                                         &*m_apstate, 
                                                         p.first)) {
              dims.push_back (p.first);
            }
            else {
              ap_dim_t i = res.size ();
              res.insert (binding_t (p.second, i));
            }
          }
          remove_dimensions (m_apstate, dims);
          std::swap (m_var_map, res);

          assert (m_var_map.size () == get_dims ());
        }


        apron_domain_ (ap_state_ptr&& apState, var_map_t&& varMap): 
            ikos::writeable (), 
            m_apstate (std::move (apState)), 
            m_var_map (std::move (varMap)) { 

          vector<ap_dim_t> dims;
          var_map_t res;
          /// XXX: we must iterate on the dimension id's to preserve
          /// order between them
          for (auto const& p: m_var_map.right) {  
            if (ap_abstract0_is_dimension_unconstrained (get_man (),
                                                         &*m_apstate, 
                                                         p.first)) {
              dims.push_back (p.first);
            }
            else {
              ap_dim_t i = res.size ();
              res.insert (binding_t (p.second, i));
            }
          }
          remove_dimensions (m_apstate, dims);
          std::swap (m_var_map, res);

          assert (m_var_map.size () == get_dims ());
        }


       public:

        apron_domain_ (bool isBot = false): 
            ikos::writeable (),
            m_apstate (apPtr (get_man(), 
                              (isBot ? 
                               ap_abstract0_bottom (get_man(), 0, 0) : 
                               ap_abstract0_top (get_man(), 0, 0))))
        { }

        ~apron_domain_ () { }

        apron_domain_ (const apron_domain_t& o): 
            ikos::writeable(), 
            m_apstate (apPtr (get_man (), ap_abstract0_copy (get_man (), &*(o.m_apstate)))),
            m_var_map (o.m_var_map)
        {  
          crab::CrabStats::count (getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
        }

        apron_domain_ (apron_domain_t&& o): 
            ikos::writeable(), 
            m_apstate (std::move (o.m_apstate)), 
            m_var_map (std::move (o.m_var_map)) { }
        
        apron_domain_t& operator=(const apron_domain_t& o) {
          crab::CrabStats::count (getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
          if (this != &o) {
            m_apstate = o.m_apstate;
            m_var_map = o.m_var_map;
          }
          return *this;
        }

        apron_domain_t& operator=(apron_domain_t&& o) {
          if (this != &o) {
            m_apstate = std::move (o.m_apstate);
            m_var_map = std::move (o.m_var_map);
          }
          return *this;
        }
        
        static apron_domain_t top() { 
          return apron_domain_t (false);
        }

        static apron_domain_t bottom() { 
          return apron_domain_t (true);
        }

        bool is_bottom() { 
          return ap_abstract0_is_bottom (get_man(), &*m_apstate);
        }

        bool is_top() { 
          return ap_abstract0_is_top (get_man(), &*m_apstate);
        }

        bool operator<=(apron_domain_t o) { 
          crab::CrabStats::count (getDomainName() + ".count.leq");
          crab::ScopedCrabStats __st__(getDomainName() + ".leq");

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
            merge_var_map (m_var_map, x, o.m_var_map, o.m_apstate);
            return ap_abstract0_is_leq (get_man(), &*x, &*o.m_apstate);
          }
        }

        void operator|=(apron_domain_t o) {
          crab::CrabStats::count (getDomainName() + ".count.join");
          crab::ScopedCrabStats __st__(getDomainName() + ".join");

          if (is_bottom() || o.is_top ())
            *this = o;
          else if (is_top () || o.is_bottom())
            return ;
          else {
            m_var_map = merge_var_map (m_var_map, m_apstate, o.m_var_map, o.m_apstate);
            m_apstate = apPtr (get_man(), 
                               ap_abstract0_join (get_man(), false, 
                                                  &*m_apstate, &*o.m_apstate));
          }
        }
        
        apron_domain_t operator|(apron_domain_t o) {
          crab::CrabStats::count (getDomainName() + ".count.join");
          crab::ScopedCrabStats __st__(getDomainName() + ".join");

          if (is_bottom() || o.is_top ())
            return o;
          else if (is_top () || o.is_bottom())
            return *this;
          else {
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_t  m = merge_var_map (m_var_map, x, o.m_var_map, o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_join (get_man(), false, 
                                                             &*x, &*o.m_apstate)), m);
          }
        }        
        
        apron_domain_t operator&(apron_domain_t o) {
          crab::CrabStats::count (getDomainName() + ".count.meet");
          crab::ScopedCrabStats __st__(getDomainName() + ".meet");

          if (is_bottom() || o.is_bottom())
            return bottom();
          else if (is_top())
            return o;
          else if (o.is_top())
            return *this;
          else{
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_t  m = merge_var_map (m_var_map, x, o.m_var_map, o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_meet (get_man(), false, 
                                                             &*x, &*o.m_apstate)), m);
          }
        }        
        
        apron_domain_t operator||(apron_domain_t o) {
          crab::CrabStats::count (getDomainName() + ".count.widening");
          crab::ScopedCrabStats __st__(getDomainName() + ".widening");

          if (is_bottom())
            return o;
          else if (o.is_bottom())
            return *this;
          else {
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_t  m = merge_var_map (m_var_map, x, o.m_var_map, o.m_apstate);
            return apron_domain_t (apPtr (get_man(), 
                                          ap_abstract0_widening (get_man(), 
                                                                 &*x, &*o.m_apstate)), m);
          }
        }        
        
        template<typename Thresholds>
        apron_domain_t widening_thresholds (apron_domain_t o, const Thresholds &ts) {
          // TODO: use thresholds
          return (*this || o);
        }

        apron_domain_t operator&&(apron_domain_t o) {
          crab::CrabStats::count (getDomainName() + ".count.narrowing");
          crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

          if (is_bottom() || o.is_bottom())
            return bottom();
          else if (is_top())
            return o;
          else if (o.is_top())
            return *this;
          else{
            ap_state_ptr x = apPtr (get_man (), ap_abstract0_copy (get_man (), &*m_apstate));
            var_map_t  m = merge_var_map (m_var_map, x, o.m_var_map, o.m_apstate);
            switch (ApronDom) {
              case APRON_OCT:
                return apron_domain_t (apPtr (get_man(), 
                                              ap_abstract0_oct_narrowing (get_man(),
                                                                          &*x, &*o.m_apstate)), m);
              case APRON_OPT_OCT:
                return apron_domain_t (apPtr (get_man(), 
                                              ap_abstract0_opt_oct_narrowing (get_man(),
                                                                              &*x, &*o.m_apstate)), m);
              case APRON_INT:
              case APRON_PK:
              default:
                CRAB_WARN (" used meet instead of narrowing: \n",
                           "make sure only a finite number of descending iterations are run.");
                return apron_domain_t (apPtr (get_man(), 
                                              ap_abstract0_meet (get_man(), false, 
                                                               &*x, &*o.m_apstate)), m);
            }
          }
        }        

        template<typename Range>
        void forget (const Range &vars) {
          crab::CrabStats::count (getDomainName() + ".count.forget");
          crab::ScopedCrabStats __st__(getDomainName() + ".forget");

          vector<ap_dim_t> vector_dims;
          std::set<ap_dim_t> set_dims;

          for (auto v: vars)  {
            if (auto dim = get_var_dim (v)) {
              vector_dims.push_back (*dim);
              set_dims.insert (*dim);
            }
          }

          if (vector_dims.empty ()) return;

          m_apstate = apPtr (get_man (), 
                             ap_abstract0_forget_array (get_man (), 
                                                        false, 
                                                        &*m_apstate, 
                                                        &vector_dims[0], vector_dims.size(), 
                                                        false));

          // -- Remove forgotten dimensions while compacting
          var_map_t res;
          /// XXX: we must iterate on the dimension id's to preserve
          /// order between them
          for (auto const& p: m_var_map.right) {  
             if (set_dims.count (p.first) <= 0) {
               ap_dim_t i = res.size ();
               res.insert (binding_t (p.second, i));
             }
          }

          remove_dimensions (m_apstate, vector_dims);
          std::swap (m_var_map, res);
        }

        void operator-=(VariableName var) {
          vector<ap_dim_t> vector_dims;
          if (auto dim = get_var_dim (var)) {
            vector_dims.push_back (*dim);
            m_apstate = apPtr (get_man (), 
                               ap_abstract0_forget_array (get_man (), 
                                                          false, 
                                                          &*m_apstate, 
                                                          &vector_dims[0], vector_dims.size(), 
                                                          false));
            // -- Remove forgotten dimensions while compacting
            var_map_t res;
            /// XXX: we must iterate on the dimension id's to preserve
            /// order between them
            for (auto const& p: m_var_map.right) {  
              if (p.first != *dim) {
                ap_dim_t i = res.size ();
                res.insert (binding_t (p.second, i));
              }
            }
            remove_dimensions (m_apstate, vector_dims);
            std::swap (m_var_map, res);
          }
        }

        // remove all variables except vars
        template<typename Range>
        void project (const Range& vars) {
          crab::CrabStats::count (getDomainName() + ".count.project");
          crab::ScopedCrabStats __st__(getDomainName() + ".project");

          if (is_bottom ()) return;
          std::set<VariableName> s1,s2,s3;
          for (auto p: m_var_map.left) s1.insert (p.first);
          s2.insert (vars.begin (), vars.end ());
          boost::set_difference (s1,s2,std::inserter (s3, s3.end ()));
          forget (s3);
        }

        interval_t operator[](VariableName v) {
          crab::CrabStats::count (getDomainName() + ".count.to_intervals");
          crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");

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
          crab::CrabStats::count (getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");

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
          crab::CrabStats::count (getDomainName() + ".count.add_constraints");
          crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

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
          CRAB_LOG("apron", 
                   crab::outs() << "--- "<< "Assume "<<csts<< " --> "<< *this<<"\n";);
        }
       
        void assign (VariableName x, linear_expression_t e) {
          crab::CrabStats::count (getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");

          if(is_bottom()) return;

          ap_texpr0_t* t = expr2texpr (e);
          assert (t);
          auto dim_x = get_var_dim_insert (x);
          m_apstate = apPtr (get_man (), 
                             ap_abstract0_assign_texpr(get_man (), false, 
                                                       &*m_apstate, 
                                                       dim_x, t, NULL));
                                                       
          ap_texpr0_free (t);
          CRAB_LOG("apron",
                   crab::outs() << "--- "<< x<< ":="<< e << " --> "<< *this<<"\n";);
        }
          
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if(is_bottom()) return;

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

          auto dim_x = get_var_dim_insert (x);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   dim_x, res, NULL));
                                                                   
          ap_texpr0_free (res);
          CRAB_LOG("apron",
                   crab::outs() << "--- "<< x<< ":="<< y<< op<< z<< " --> "<< *this<<"\n";);
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if(is_bottom()) return;

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

          auto dim_x = get_var_dim_insert (x);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   dim_x, res, NULL));

          ap_texpr0_free (res);
          CRAB_LOG("apron",
                   crab::outs() << "--- "<< x<< ":="<< y<< op<< z<< " --> "<< *this<<"\n";);
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if(is_bottom()) return;

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

          auto dim_x = get_var_dim_insert (x);
          m_apstate = apPtr (get_man (), ap_abstract0_assign_texpr(get_man (), false, 
                                                                   &*m_apstate, 
                                                                   dim_x, res, NULL));

          ap_texpr0_free (res);
          CRAB_LOG("apron",
                   crab::outs() << "--- "<< x<< ":="<< x<< op<< k<< " --> "<< *this<<"\n";);
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
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

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
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

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
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

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
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

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
          auto dim_x = get_var_dim_insert (x);
          m_apstate = apPtr (get_man(),
                             ap_abstract0_expand(get_man (), false, &* m_apstate, 
                                                 dim_x, 1));
          // --- the additional dimension is put at the end of integer
          //     dimensions.
          m_var_map.insert (binding_t (dup, get_dims () - 1));            
        }

        void normalize () {
          ap_abstract0_canonicalize (get_man (), &*m_apstate);
        }

        void write(crab_os& o) {
          if(is_bottom()){
            o << "_|_";
            return;
          }
          else if (is_top()){
            o << "{}";
            return;
          }
          else {
            // dump ();
            linear_constraint_system_t inv = to_linear_constraint_system ();
            o << inv;
          }
        }          

        static std::string getDomainName () {
          switch (ApronDom) {
            case APRON_INT:     return "ApronIntervals"; 
            case APRON_OCT:     return "ApronOctagon"; 
            case APRON_OPT_OCT: return "ApronOptimizedOctagon"; 
            case APRON_PK:      return "ApronNewPolka";
            default: CRAB_ERROR("Unknown apron domain");
          }
        }
      }; 

      // Quick wrapper which uses shared references with copy-on-write.
      template<class Number, class VariableName, apron_domain_id_t ApronDom>
      class apron_domain : public ikos::writeable,
                           public numerical_domain<Number, VariableName >,
                           public bitwise_operators<Number,VariableName >,
                           public division_operators<Number, VariableName >,
                           public array_operators< Number, VariableName >,
                           public pointer_operators< Number, VariableName > {
       public:
        using typename numerical_domain< Number, VariableName >::linear_expression_t;
        using typename numerical_domain< Number, VariableName >::linear_constraint_t;
        using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName >::variable_t;
        using typename numerical_domain< Number, VariableName >::number_t;
        using typename numerical_domain< Number, VariableName >::varname_t;
        typedef typename linear_constraint_t::kind_t constraint_kind_t;
        typedef interval<Number>  interval_t;
        
        typedef apron_domain_<Number, VariableName, ApronDom> apron_domain_impl_t;
        typedef std::shared_ptr<apron_domain_impl_t> apron_domain_ref_t;
        typedef apron_domain<Number, VariableName, ApronDom> apron_domain_t;

        apron_domain(apron_domain_ref_t _ref) : norm_ref(_ref) { }
        
        apron_domain(apron_domain_ref_t _base, apron_domain_ref_t _norm) 
            : base_ref(_base), norm_ref(_norm)
        { }
        
        apron_domain_t create(apron_domain_impl_t&& t)
        {
          return std::make_shared<apron_domain_impl_t>(std::move(t));
        }
        
        apron_domain_t create_base(apron_domain_impl_t&& t)
        {
          apron_domain_ref_t base = std::make_shared<apron_domain_impl_t>(t);
          apron_domain_ref_t norm = std::make_shared<apron_domain_impl_t>(std::move(t));  
          return apron_domain_t(base, norm);
        }
        
        void lock(void)
        {
          // Allocate a fresh copy.
          if(!norm_ref.unique())
            norm_ref = std::make_shared<apron_domain_impl_t>(*norm_ref);
          base_ref.reset();
        }
        
       public:
        
        static apron_domain_t top() { return apron_domain(false); }
        
        static apron_domain_t bottom() { return apron_domain(true); }
        
        apron_domain(bool is_bottom = false)
            : norm_ref(std::make_shared<apron_domain_impl_t>(is_bottom)) { }
        
        apron_domain(const apron_domain_t& o)
            : base_ref(o.base_ref), norm_ref(o.norm_ref)
        { }
        
        apron_domain& operator=(const apron_domain_t& o) {
          base_ref = o.base_ref;
          norm_ref = o.norm_ref;
          return *this;
        }
        
        apron_domain_impl_t& base(void) {
          if(base_ref)
            return *base_ref;
          else
            return *norm_ref;
        }
        apron_domain_impl_t& norm(void) { return *norm_ref; }
        
        bool is_bottom() { return norm().is_bottom(); }
        bool is_top() { return norm().is_top(); }
        bool operator<=(apron_domain_t& o) { return norm() <= o.norm(); }
        void operator|=(apron_domain_t o) { lock(); norm() |= o.norm(); }
        apron_domain_t operator|(apron_domain_t o) { return create(norm() | o.norm()); }
        apron_domain_t operator||(apron_domain_t o) { return create_base(base() || o.norm()); }
        apron_domain_t operator&(apron_domain_t o) { return create(norm() & o.norm()); }
        apron_domain_t operator&&(apron_domain_t o) { return create(norm() && o.norm()); }
        
        template<typename Thresholds>
        apron_domain_t widening_thresholds (apron_domain_t o, const Thresholds &ts) {
          return create_base(base().template widening_thresholds<Thresholds>(o.norm(), ts));
        }
        
        void normalize() { norm(); }
        void operator+=(linear_constraint_system_t csts) { lock(); norm() += csts; } 
        void operator-=(VariableName v) { lock(); norm() -= v; }
        interval_t operator[](VariableName x) { return norm()[x]; }
        void set(VariableName x, interval_t intv) { lock(); norm().set(x, intv); }
        
        template<typename Range>
        void forget (Range vs) { lock(); norm().forget(vs); }
        void assign(VariableName x, linear_expression_t e) { lock(); norm().assign(x, e); }
        void apply(operation_t op, VariableName x, VariableName y, Number k) {
          lock(); norm().apply(op, x, y, k);
        }
        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
          lock(); norm().apply(op, x, y, width);
        }
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          lock(); norm().apply(op, x, k, width);
        }
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          lock(); norm().apply(op, x, y, k);
        }
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          lock(); norm().apply(op, x, y, z);
        }
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          lock(); norm().apply(op, x, y, z);
        }
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          lock(); norm().apply(op, x, y, z);
        }
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          lock(); norm().apply(op, x, y, k);
        }
        void expand (VariableName x, VariableName y) { lock(); norm().expand(x, y); }
        
        template<typename Range>
        void project (Range vs) { lock(); norm().project(vs); }
        
        template <typename NumDomain>
        void push (const VariableName& x, NumDomain&inv){ lock(); norm().push(x, inv); }
        
        void write(crab_os& o) { norm().write(o); }
        
        linear_constraint_system_t to_linear_constraint_system () {
          return norm().to_linear_constraint_system();
        }
        static std::string getDomainName () { return apron_domain_impl_t::getDomainName(); }

       protected:  
        apron_domain_ref_t base_ref;  
        apron_domain_ref_t norm_ref;
      };

      // --- global datastructures

      template<typename N, typename V, apron_domain_id_t D>
      ap_manager_t* apron_domain_<N,V,D>::m_apman = nullptr;

      // -- domain traits
      template<typename Number, typename VariableName, apron_domain_id_t ApronDom>
      class domain_traits <apron_domain<Number,VariableName, ApronDom> > {
       public:
        
        typedef apron_domain<Number,VariableName, ApronDom> apron_domain_t;
        
        static void normalize (apron_domain_t& inv) { 
          inv.normalize ();
        }
        
        template <typename Iter>
        static void forget (apron_domain_t& inv, Iter it, Iter end) {
          inv.forget (boost::make_iterator_range (it, end));
        }
        
        template <typename Iter >
        static void project (apron_domain_t& inv, Iter it, Iter end) {
          inv.project (boost::make_iterator_range (it, end));
        }
        
        static void expand (apron_domain_t& inv, VariableName x, VariableName new_x) {
          inv.expand (x, new_x);
        }
        
      };
   
   } // namespace domains
}// namespace crab
#endif 
#endif 
