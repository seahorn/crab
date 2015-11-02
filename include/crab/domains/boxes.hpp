#ifndef BOXES_DOMAIN_HPP
#define BOXES_DOMAIN_HPP

/*** 
   A disjunctive domain of intervals based on "Boxes: A Symbolic
   Abstract Domain of Boxes" by A. Gurfinkel and S. Chaki published
   SAS'10.
***/


// Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

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

#include <boost/bimap.hpp>
#include <boost/shared_ptr.hpp>

using namespace crab::domains::ldd;

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
         // --- map from crab variable index to ldd term index
         typedef boost::bimap< VariableName , int > var_bimap_t;
         typedef boost::shared_ptr<var_bimap_t> VarMapPtr;
         typedef typename var_bimap_t::value_type binding_t;

        public:

        LddNodePtr m_ldd;
        static LddManager* m_ldd_man;
        static VarMapPtr m_var_map;

        static LddManager* get_ldd_man () {
          if (!m_ldd_man) {
            const unsigned Sz = 100; // FIXME: make me template parameter
            DdManager* cudd = Cudd_Init (0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
            theory_t* theory = tvpi_create_boxz_theory (Sz);
             m_ldd_man = Ldd_Init (cudd, theory);
             
             Cudd_AutodynEnable (cudd, CUDD_REORDER_GROUP_SIFT);
          }
          return m_ldd_man;
        }

        static VarMapPtr get_var_map () {
          if (!m_var_map) {
            m_var_map = VarMapPtr (new var_bimap_t ());            
          }
          return m_var_map;
        }

        boxes_domain (LddNodePtr ldd): m_ldd (ldd) { }

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

        int getVarId (VariableName v) {
          auto it = get_var_map()->left.find (v);
          if (it != get_var_map()->left.end ())
             return it->second;
          else {
            int id = get_var_map ()->size ();
            get_var_map ()->insert (binding_t (v, ++id));
            return id;
          }
        }

        VariableName getVarName (int v) {
          auto it = get_var_map ()->right.find (v);
          if (it != get_var_map ()->right.end ())
             return it->second;
          else {
             CRAB_ERROR ("Index", v, "cannot be mapped back to a variable name");
          }
        }

        // convex approximation
        LddNodePtr approx (LddNodePtr v) {
          return lddPtr (get_ldd_man(), Ldd_TermMinmaxApprox(get_ldd_man(), &*v));
        }

        LddNodePtr join (LddNodePtr v1, LddNodePtr v2) {
          return lddPtr (get_ldd_man(), Ldd_Or (get_ldd_man(), &*v1, &*v2));
        }

        /** return term for variable v, neg for negation of variable */
        linterm_t termForVal(VariableName v, bool neg = false) 
        {
          int varId = getVarId (v);
          int sgn = neg ? -1 : 1;
          linterm_t term = 
              Ldd_GetTheory (get_ldd_man())->create_linterm_sparse_si (&varId, &sgn, 1);
          return term; 
        }
        
        /** v := a * x + k, where a, k are constants and x variable */
        void apply (VariableName v, VariableName x, Number a, Number k) {
          if (is_top () || is_bottom ()) return ;

          linterm_t t = termForVal (v);
          linterm_t r = termForVal (x);
          
          mpq_class aa = a.get (); // FIXME: implicit cast from mpz_class to mpq_class
          mpq_class kk = k.get (); // FIXME: implicit cast from mpz_class to mpq_class
          
          constant_t aaa = (constant_t) tvpi_create_cst (aa.get_mpq_t ());
          constant_t kkk = (constant_t) tvpi_create_cst (kk.get_mpq_t ());        
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, r, aaa, kkk, kkk));
              
          Ldd_GetTheory (get_ldd_man())->destroy_cst(aaa);
          Ldd_GetTheory (get_ldd_man())->destroy_cst(kkk);
        }

        Number numFromLddCst (constant_t cst, theory_t *theory) {
          mpq_class v;
          // XXX We know that the theory is tvpi, use its method direclty.
          tvpi_cst_set_mpq (v.get_mpq_t (), (tvpi_cst_t)cst);

          // FIXME: it's assuming that Number has a constructor that
          // takes mpz_class
          mpz_class n = static_cast<mpz_class> (v);
          return Number (n);
        }

        linear_expression_t exprFromLddTerm (linterm_t term, theory_t *theory) {
          linear_expression_t e (0);
          for(size_t i = 0;i < theory->term_size (term); i++) {
            Number k = numFromLddCst (theory->term_get_coeff (term,i), theory);
            VariableName v =  getVarName (theory->term_get_var(term,i));
            e = e + (k * linear_expression_t (v));
          }
          return e;
        }

        linear_constraint_t cstFromLddCons(lincons_t lincons, theory_t *theory
                                           /*, bool isNegated*/) {

          linear_expression_t lhs = exprFromLddTerm(theory->get_term(lincons), theory);
          Number rhs = numFromLddCst(theory->get_constant(lincons), theory);

          // if (!isNegated) {
          if (theory->is_strict (lincons)) {
            // lhs < rhs <-> lhs+1 <= rhs but only for integers!
            linear_expression_t e = lhs +1;
            e = e - rhs;
            return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
          }
          else  {
            // lhs <= rhs
            linear_expression_t e = lhs - rhs;
            return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
          }
          // }
          // else {
          //   if (theory->is_strict (lincons)) {
          //     // lhs >= rhs <-> rhs <= lhs
          //     linear_expression_t e = rhs - lhs;
          //     return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
          //   }
          //   else  {
          //     // lhs > rhs <-> rhs +1 <= lhs
          //     linear_expression_t e = rhs +1;
          //     e = e - lhs;
          //     return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
          //   }
          // }
        }

        // Pre: n is convex and cannot be bottom
        void toLinCstSysRecur(LddNode* n, 
                              linear_constraint_system_t& csts) {

          LddNode *N = Ldd_Regular (n);
          if (N == Ldd_GetTrue (get_ldd_man ())) return;


          auto cst = cstFromLddCons (Ldd_GetCons(get_ldd_man (), N), 
                                     Ldd_GetTheory (get_ldd_man ()));

          if (Ldd_Regular (Ldd_T(N)) == Ldd_GetTrue (get_ldd_man ()) &&
              Ldd_Regular (Ldd_E(N)) == Ldd_GetTrue (get_ldd_man ()))
            csts += cst;
          else if (Ldd_Regular (Ldd_T(N)) == Ldd_GetTrue (get_ldd_man ()))
            csts += cst.negate ();
          else
            csts += cst;
            
          toLinCstSysRecur (Ldd_T (N), csts);
          toLinCstSysRecur (Ldd_E (N), csts);

        }

       public:
        
        boxes_domain(): ikos::writeable() { 
          m_ldd = lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man()));
        }
        
        ~boxes_domain () { }
                
        static boxes_domain_t top() { 
          return boxes_domain_t (lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man())));
        }
        
        static boxes_domain_t bottom() {
          return boxes_domain_t (lddPtr (get_ldd_man(), Ldd_GetFalse (get_ldd_man())));
        }
        
        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable(), m_ldd (other.m_ldd) { }
        
        boxes_domain_t& operator=(boxes_domain_t other) {
           m_ldd = other.m_ldd;
           return *this;
        }
        
        bool is_bottom() { 
          return &*m_ldd == Ldd_GetFalse (get_ldd_man());
        }
        
        bool is_top() { 
          return &*m_ldd == Ldd_GetTrue (get_ldd_man());
        }
        
        bool operator<=(boxes_domain_t other) {
          bool res = Ldd_TermLeq (get_ldd_man(), &(*m_ldd), &(*other.m_ldd));

          CRAB_DEBUG ("Check if ", *this, " <= ", other, " ---> ", res);
          return res;
        }
        
        boxes_domain_t operator|(boxes_domain_t other) {
          return boxes_domain_t (join (m_ldd, other.m_ldd));
        }
        
        boxes_domain_t operator&(boxes_domain_t other) {
          return boxes_domain_t (lddPtr (get_ldd_man(), 
                                         Ldd_And (get_ldd_man(), &*m_ldd, &*other.m_ldd)));
        }

        boxes_domain_t operator||(boxes_domain_t other) {
          // It is not necessarily true that the new value is bigger
          // than the old value so we apply 
          // widen(old, new) = widen (old, (join (old,new)))
          LddNodePtr v = join (m_ldd, other.m_ldd); 
          LddNodePtr w = lddPtr (get_ldd_man (), 
                                 Ldd_BoxWiden2 (get_ldd_man (), &*m_ldd, &*v));

#if 1
          /** ensure that 'w' is only the fronteer of the compution. 
              Not sure whether this is still a widening though 
          */
          w = lddPtr (get_ldd_man (), 
                      Ldd_And (get_ldd_man (), &*w, Ldd_Not (&*m_ldd)));
          /** ensure the output is at least as big as newV */
          w = lddPtr (get_ldd_man (), Ldd_Or (get_ldd_man (), &*w, &*other.m_ldd));
#endif 
          boxes_domain_t res (w); 

          CRAB_DEBUG ("Widening ", *this, " and ", other, "=", res);
          return res;
        }
        
        boxes_domain_t operator&& (boxes_domain_t other) {
          // TODO: narrowing
          boxes_domain_t res(*this);

          CRAB_DEBUG ("Narrowing ", *this, " and ", other, "=", res);
          return res;
        }
        
        void operator-=(VariableName var) {
          int id = getVarId (var);
          m_ldd =  lddPtr (get_ldd_man(), Ldd_ExistsAbstract (get_ldd_man(), &*m_ldd, id));
        }

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          for (auto v: boost::make_iterator_range (begin, end))
            operator-= (v);
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          CRAB_WARN("TODO: BOXES project operation");
        }

        /** 
            Restricts n by a given constraint e.
            Currently, e can be one of 
               TRUE | FALSE | x < k | x <= k | -x < k | -x <= k | x == k | x != k 
            where x is a variable and k is a constant
        */	
        void operator+= (linear_constraint_t cst)
        {
          if (is_bottom () || cst.is_tautology ())  
            return;

          if (cst.is_contradiction ())  {
            m_ldd = lddPtr (get_ldd_man(), Ldd_GetFalse (get_ldd_man()));
            return;
          }
    
          linear_expression_t exp = cst.expression();
          if (exp.size() > 1) {
            CRAB_WARN ("BOXES only supports constraints with at most one variable.");
            return;
          }
          assert (exp.size () > 0);

          Number k = -exp.constant(); 
          typename linear_expression_t::iterator it=exp.begin();
          Number coef = it->first;
          variable_t x  = it->second;      

          if (!(coef == 1 || coef == -1)) {
            CRAB_WARN ("BOXES only support constraint with unit coefficients");
            return;
          }

          if (coef == -1) k = -k;

          mpq_class kk = k.get (); // FIXME: implicit cast from mpz_class to mpq_class
          constant_t kkk = (constant_t) tvpi_create_cst (kk.get_mpq_t ());
          linterm_t term = termForVal (x.name ());

          auto theory = Ldd_GetTheory (get_ldd_man());
          if (cst.is_equality ()) {  // x == k <-> x<=k and !(x<k)
            // x<=k
            lincons_t cons1 = theory->create_cons (term, 0 /*non-strict*/, kkk);
            LddNodePtr res = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons1));
            theory->destroy_lincons (cons1);
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*res));
              // !(x<k)
            lincons_t cons2 = theory->create_cons (term, 1 /*strict*/, kkk);
            res = lddPtr (get_ldd_man(), Ldd_Not (theory->to_ldd (get_ldd_man(), cons2)));
            theory->destroy_lincons (cons2);
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*res));
          }
          else if (cst.is_inequality()) {
            if (coef == 1) {
              // x <= k;
              lincons_t cons = theory->create_cons (term, 0 /*non-strict*/, kkk);
              LddNodePtr res = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons));
              theory->destroy_lincons (cons);
              m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*res));
            }
            else {
              // -x <= k <-> x >= -k <-> ! (x < -k)
              lincons_t cons = theory->create_cons (term, 1 /*strict*/, kkk);
              LddNodePtr res = lddPtr (get_ldd_man(), Ldd_Not (theory->to_ldd (get_ldd_man(), cons)));
              theory->destroy_lincons (cons);
              m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*res));
            }
          }
          else { // x != k or x!=-k
            CRAB_WARN("TODO BOXES disequalities not implemented");
          }

          CRAB_DEBUG("Assume(", cst, ") --> ", *this);
        }    
                    

        void operator += (linear_constraint_system_t csts) {
          for (auto cst : csts)
            operator += (cst);
        }

        void set (VariableName v, interval_t ival) {
          
          if (is_bottom ()) return ;
          
          constant_t kmin = NULL, kmax = NULL;
          
          if (boost::optional <Number> l = ival.lb ().number ()) {
            // FIXME: implicit cast from mpz_class to mpq_class            
            mpq_class val = (*l).get ();
            kmin = (constant_t) tvpi_create_cst (val.get_mpq_t ());
          }
          
          if (boost::optional <Number> u = ival.ub ().number ()) {
            // FIXME: implicit cast from mpz_class to mpq_class
            mpq_class val = (*u).get ();
            kmax = (constant_t)tvpi_create_cst (val.get_mpq_t ());
          }
          
          linterm_t t = termForVal(v);
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, kmin, kmax));
              
          if (kmin) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmin);
          if (kmax) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmax);
        }

        interval_t operator[](VariableName v) { 
          CRAB_WARN("Projection to intervals not implemented");
          return interval_t::top ();
        }
                
        void assign (VariableName x, linear_expression_t e) {
          if (is_bottom ()) return;

          if (e.is_constant ()) {
            Number k = e.constant ();
            // FIXME: implicit cast from mpz_class to mpq_class
            mpq_class kk = k.get (); 
            constant_t kkk = (constant_t) tvpi_create_cst (kk.get_mpq_t ());
            linterm_t t = termForVal (x);
            m_ldd = lddPtr(get_ldd_man(), 
                           Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, kkk, kkk));
            Ldd_GetTheory (get_ldd_man())->destroy_cst(kkk);
            
          }
          else if  (optional<variable_t> v = e.get_variable()){
            VariableName y = (*v).name();
            if (!(x==y))
              apply (x, y, 1, 0);      
          }
          else {
            CRAB_WARN("BOXES only supports a cst or var on the rhs of assignment");
            this->operator-=(x);
          }
          CRAB_DEBUG("---", x, ":=", e,"\n",*this);
        }
        
        void apply (operation_t op, VariableName x, VariableName y, Number k) {
          switch(op){
            case OP_ADDITION:
              apply (x, y, 1, k);
              break;
            case OP_SUBTRACTION:
              apply (x, y, 1, -k);
              break;
            case OP_MULTIPLICATION:
              apply (x, y, k, 0);
              break;
            case OP_DIVISION:
              // x = y / k <-> x = 1/k *y + 0
            default: 
              operator-= (x);
          }
          
          CRAB_DEBUG("apply ", x, " := ", y, " ", op, " ", k, " --- ", *this);
        }

        void apply(operation_t op, VariableName x, Number k) {
          switch(op){
            case OP_ADDITION:
              apply (x, x, 1, k);
              break;
            case OP_SUBTRACTION:
              apply (x, x, 1, -k);
              break;
            case OP_MULTIPLICATION:
              apply (x, x, k, 0);
              break;
            case OP_DIVISION:
              // x = x / k <-> x = 1/k *x + 0
            default: 
              operator-= (x);
          }
          CRAB_DEBUG("apply ", x, " := ", x, " ", op, " ", k, "---", *this);
        }

        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          CRAB_WARN("BOXES arithmetic operation with two variables on the rhs not implemented");
          // TODO: fall back to intervals and perform operation there          
          operator-=(x);
        }
        
        // bitwise_operators_api
        void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
          // since reasoning about infinite precision we simply assign and
          // ignore the width.
          assign(x, linear_expression_t(y));
        }
        
        void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
          // since reasoning about infinite precision we simply assign
          // and ignore the width.
          assign(x, k);
        }      
  
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          CRAB_WARN("BOXES bitwise operation ", op, " not implemented");
          // TODO: fall back to intervals and perform operation there
          operator-=(x);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          CRAB_WARN("BOXES bitwise operation ", op, " not implemented");
          // TODO: fall back to intervals and perform operation there
          operator-=(x);
        }
        
        // division_operators_api
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          CRAB_WARN("BOXES division operation", op, " not implemented");
          // TODO: fall back to intervals and perform operation there
          operator-=(x);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          CRAB_WARN("BOXES division operation", op, " not implemented");
          // TODO: fall back to intervals and perform operation there
          operator-=(x);
        }
        
        linear_constraint_system_t to_linear_constraint_system () {
          linear_constraint_system_t csts;
    
          if(is_bottom ()) {
            csts += linear_constraint_t (linear_expression_t (Number(1)) == 
                                         linear_expression_t (Number(0)));
            return csts;
          }

          if(is_top ()) {
            csts += linear_constraint_t (linear_expression_t (Number(1)) == 
                                         linear_expression_t (Number(1)));
            return csts;
          }

          // --- produce convex approximation
          LddNodePtr v = approx (m_ldd);

          // --- extract linear inequalities from the convex ldd
          toLinCstSysRecur(&*v, csts);
          
          return csts;
        }
        
        void write(ostream& o) {
          auto csts = to_linear_constraint_system ();
          o << csts;
#if 1
          // -- debugging
          o.flush ();
          dump ();
#endif 
        }
        
        const char* getDomainName () const {return "BOXES";}  
        
      }; 

     template<typename N, typename V>
     LddManager* boxes_domain<N, V>::m_ldd_man = nullptr;

     template<typename N, typename V>
     typename boxes_domain<N,V>::VarMapPtr boxes_domain<N, V>::m_var_map = nullptr;
      
   } // namespace domains

   namespace domain_traits {

     using namespace domains;

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
