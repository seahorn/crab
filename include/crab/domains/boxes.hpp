#ifndef BOXES_DOMAIN_HPP
#define BOXES_DOMAIN_HPP

/************************************************************************** 
 * A disjunctive domain of intervals based on "Boxes: A Symbolic
 * Abstract Domain of Boxes" by A. Gurfinkel and S. Chaki published
 * SAS'10.
 **************************************************************************/

#include <crab/config.h>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/operators_api.hpp>

using namespace boost;
using namespace ikos;


#ifndef HAVE_LDD
/*
 * Dummy implementation if ldd not found 
 */
#define LDD_NOT_FOUND "No LDD. Run cmake with -DUSE_LDD=ON"
namespace crab {
   namespace domains {
      template<typename Number, typename VariableName, size_t LddSize = 100>
      class boxes_domain: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName >,
         public array_operators< Number, VariableName>,
         public pointer_operators< Number, VariableName> {
              
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;
        typedef boxes_domain <Number, VariableName, LddSize> boxes_domain_t;
        typedef interval <Number> interval_t;
        typedef int LddNodePtr;

        boxes_domain(): ikos::writeable() { }    

        LddNodePtr getLdd () 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        VariableName getVarName (int v) const 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        static boxes_domain_t top() { CRAB_ERROR (LDD_NOT_FOUND); }

        static boxes_domain_t bottom() { CRAB_ERROR (LDD_NOT_FOUND); }

        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable() { }
        
        bool is_bottom() { CRAB_ERROR (LDD_NOT_FOUND); }

        bool is_top() { CRAB_ERROR (LDD_NOT_FOUND); }

        bool operator<=(boxes_domain_t other) { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void operator|=(boxes_domain_t other)
        { CRAB_ERROR (LDD_NOT_FOUND); }

        boxes_domain_t operator|(boxes_domain_t other)
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator&(boxes_domain_t other) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator||(boxes_domain_t other)
        { CRAB_ERROR (LDD_NOT_FOUND); }

        template<typename Thresholds>
        boxes_domain_t widening_thresholds (boxes_domain_t other, 
                                            const Thresholds &ts)
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        boxes_domain_t operator&& (boxes_domain_t other) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
        
        void operator-=(VariableName var) 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        interval_t operator[](VariableName v) 
        { CRAB_ERROR (LDD_NOT_FOUND); }

        void set(VariableName v, interval_t ival) 
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
        
        void write(crab_os& o) 
        { CRAB_ERROR (LDD_NOT_FOUND); }
          
        static std::string getDomainName () {
          return "Dummy Boxes";
        }  
      }; 

   } // namespace domains
}// namespace crab

#else
/* 
 *  Real implementation starts here 
 */

#include <crab/domains/ldd/ldd.hpp>
#include <crab/domains/ldd/ldd_print.hpp>
#include <crab/domains/domain_traits.hpp>
#include <boost/bimap.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>

namespace crab {

   namespace domains {

      using namespace crab::domains::ldd;       

      /*
       * The wrapper has two global datastructures:
       * 1) a ldd manager and 
       * 2) a map from VariableName to ldd dimension.
       *
       * FIXME: since the ldd manager is shared we need to fix a
       * single size for all ldds. Since ldds are sparse we can fix a
       * size big enough for our programs.
       */
      template<typename Number, typename VariableName, size_t LddSize = 3000>
      class boxes_domain: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName >,
         public array_operators< Number, VariableName>,
         public pointer_operators< Number, VariableName> {

        typedef interval_domain <Number, VariableName> interval_domain_t;
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;

        typedef boxes_domain <Number, VariableName, LddSize> boxes_domain_t;
        typedef interval <Number> interval_t;

       private:
        // --- map from crab variable index to ldd term index
        typedef boost::bimap< VariableName , int > var_bimap_t;
        typedef boost::shared_ptr<var_bimap_t> var_map_ptr;
        typedef typename var_bimap_t::value_type binding_t;

        LddNodePtr m_ldd;
        static LddManager* m_ldd_man;
        static var_map_ptr m_var_map;

        static LddManager* get_ldd_man () {
          if (!m_ldd_man) {
            DdManager* cudd = Cudd_Init (0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
            theory_t* theory = tvpi_create_boxz_theory (LddSize);
            CRAB_LOG ("boxes",crab::outs() << "Created a ldd of size " << LddSize <<"\n";);
            m_ldd_man = Ldd_Init (cudd, theory);
            Cudd_AutodynEnable (cudd, CUDD_REORDER_GROUP_SIFT);
          }
          return m_ldd_man;
        }

        static int num_of_vars () {
          if (!m_var_map) return 0;
          
          return (*m_var_map).left.size ();
        }

        static var_map_ptr get_var_map () {
          if (!m_var_map) {
            m_var_map = var_map_ptr (new var_bimap_t ());            
          }
          return m_var_map;
        }
       
        int get_var_dim (VariableName v) const {
          auto it = get_var_map()->left.find (v);
          if (it != get_var_map()->left.end ()) {
            return it->second;
          } else {
            int id = get_var_map ()->size ();
            if (id >= LddSize) {
              CRAB_ERROR ("The Ldd size of ", LddSize, " needs to be larger");
            }
            get_var_map ()->insert (binding_t (v, id));
            return id;
          }
        }

        inline constant_t mkCst (Number k) {
          mpq_class kk ((mpz_class) k); 
          return (constant_t) tvpi_create_cst (kk.get_mpq_t ());
        }
        
        // convex approximation
        LddNodePtr convex_approx () const {
          return lddPtr (get_ldd_man(), Ldd_TermMinmaxApprox(get_ldd_man(), &*m_ldd));
        }

        // pre: ldd is not either true or false
        void project (LddNodePtr& ldd, VariableName v) const {
          for (auto p: (*m_var_map).left) {
            if (!(p.first == v)) {
              int id = get_var_dim (p.first);
              ldd = lddPtr (get_ldd_man(), 
                            Ldd_ExistsAbstract (get_ldd_man(), &*ldd, id));
              
            }
          }
        }

        // non-relational approximation (but still non-convex)
        // expensive operation
        LddNodePtr non_relational_approx () const {
          if (is_top () || is_bottom ())
            return m_ldd;
                                  
          LddNodePtr res = lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man()));
          for (auto p: (*m_var_map).left) {
            LddNodePtr tmp (m_ldd);
            project (tmp, p.first);
            res = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*res, &*tmp));
          }
          return res;
        }

        LddNodePtr join (LddNodePtr v1, LddNodePtr v2) {
          return lddPtr (get_ldd_man(), Ldd_Or (get_ldd_man(), &*v1, &*v2));
        }

        /** return term for variable v, neg for negation of variable */
        linterm_t termForVal(VariableName v, bool neg = false) 
        {
          int varId = get_var_dim (v);
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

          constant_t c1 = mkCst (a);
          constant_t c2 = mkCst (k);
          // FIXME: Ldd_TermReplace seems to leak memory
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, r, c1, c2, c2));

          Ldd_GetTheory (get_ldd_man())->destroy_term(t);
          Ldd_GetTheory (get_ldd_man())->destroy_term(r);              
          Ldd_GetTheory (get_ldd_man())->destroy_cst(c1);
          Ldd_GetTheory (get_ldd_man())->destroy_cst(c2);
        }

        Number numFromLddCst (constant_t cst, theory_t *theory) {
          mpq_class v;
          // XXX We know that the theory is tvpi, use its method direclty.
          tvpi_cst_set_mpq (v.get_mpq_t (), (tvpi_cst_t) cst);
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

        linear_constraint_t cstFromLddCons(lincons_t lincons, theory_t *theory) {

          linear_expression_t lhs = exprFromLddTerm(theory->get_term(lincons), theory);
          Number rhs = numFromLddCst(theory->get_constant(lincons), theory);

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

        // Restricts the ldd by the constraint e
        // where e can be one of 
        //    c*x <= k | c*x == k | c*x != k 
        //    and c is 1 or -1 and k is an integer constant
        typedef typename linear_constraint_t::kind_t kind_t;

        void add_unit_constraint (Number coef, variable_t x, kind_t kind, Number k) {

          assert (coef == 1 || coef == -1);
          
          auto theory = Ldd_GetTheory (get_ldd_man());
          if (kind == kind_t::EQUALITY) {  // x == k <-> x<=k and !(x<k)
            // x<=k
            linterm_t term1 = termForVal (x.name (), (coef == 1 ? false : true));
            constant_t c1 = mkCst (k);
            lincons_t cons1 = theory->create_cons (term1, 0 /*non-strict*/, c1);
            LddNodePtr n1 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons1));
            theory->destroy_lincons (cons1);
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*n1));
            // !(x<k)
            linterm_t term2 = termForVal (x.name (), (coef == 1 ? false : true));
            constant_t c2 = mkCst (k);
            lincons_t cons2 = theory->create_cons (term2, 1 /*strict*/, c2);
            LddNodePtr n2 = lddPtr (get_ldd_man(), Ldd_Not (theory->to_ldd (get_ldd_man(), cons2)));
            theory->destroy_lincons (cons2);
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*n2));
          }
          else if (kind == kind_t::INEQUALITY) {
            // case 1:  x <= k  
            // case 2: -x <= k
            linterm_t term = termForVal (x.name (), (coef == 1 ? false : true));
            constant_t c = mkCst (k);
            lincons_t cons = theory->create_cons (term, 0 /*non-strict*/, c);
            LddNodePtr n = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons));
            theory->destroy_lincons (cons);
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*n));
          }
          else { // assert (kind == kind_t::DISEQUALITY)
            // case 1:  x != k  <->  x < k OR  x > k <->  x < k OR -x < -k
            // case 2: -x != k  <-> -x < k OR -x > k <-> -x < k OR  x < -k
            linterm_t term1 = termForVal (x.name (), (coef == 1 ? false : true));
            constant_t c1 = mkCst (k);
            lincons_t cons1 = theory->create_cons (term1, 1 /*strict*/, c1);
            LddNodePtr n1 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons1));
            theory->destroy_lincons (cons1);
            
            linterm_t term2 = termForVal (x.name (), (coef == 1 ? true : false));
            constant_t c2 = mkCst (-k);
            lincons_t cons2 = theory->create_cons (term2, 1 /*strict*/, c2);
            LddNodePtr n2 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons2));
            theory->destroy_lincons (cons2);
            
            LddNodePtr n3 = lddPtr (get_ldd_man (), Ldd_Or (get_ldd_man(), &*n1, &*n2));
            m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*n3));
          }
        } 

        boxes_domain (LddNodePtr ldd): m_ldd (ldd) {
#if 1
          const int CST_FACTOR = 10; /* JN: some magic constant factor */
          int threshold = num_of_vars () * CST_FACTOR;
          if (threshold > 0) {
            int before_num_paths = Ldd_PathSize (NULL, &*m_ldd);
            // -- we throw away relationships between variables if the
            //    number of ldd paths is not linear with some constant
            //    factor in the size of the diagram.
            if (before_num_paths >= threshold) {
              m_ldd = non_relational_approx ();
              int after_num_paths = Ldd_PathSize (NULL, &*m_ldd);
              if (after_num_paths == before_num_paths) {
                // -- if the previous approximation did not work then
                // -- we further simplify the ldd by making it convex
                m_ldd = convex_approx ();                
                CRAB_WARN ("ldd is growing too fast (number of paths= " , before_num_paths, "). ", 
                           "After making convex the ldd ",
                           "the number of paths is ", Ldd_PathSize (NULL, &*m_ldd));
                           
              }
              else {
                CRAB_WARN ("ldd is growing too fast (number of paths= " , before_num_paths, "). ", 
                           "After deleting relationships between variables ",
                           "the number of paths is ", after_num_paths);
              }
            }
          }
#endif 
        }

       public:

        LddNodePtr getLdd () { 
          return m_ldd; 
        }

        VariableName getVarName (int v) const {
          auto it = get_var_map ()->right.find (v);
          if (it != get_var_map ()->right.end ())
             return it->second;
          else {
             CRAB_ERROR ("Index ", v, " cannot be mapped back to a variable name");
          }
        }

        boxes_domain(): ikos::writeable() { 
          m_ldd = lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man()));
        }
        
        ~boxes_domain () { 
          // DdManager *cudd = nullptr;
          // theory_t *theory = nullptr;
          // if (m_ldd_man)  {
	//   cudd = Ldd_GetCudd (m_ldd_man);
	//   theory = Ldd_GetTheory (m_ldd_man);
	//   Ldd_Quit (m_ldd_man);
	// }
          // if (theory) tvpi_destroy_theory(theory);
          // if (cudd) Cudd_Quit(cudd);
        }
                
        static boxes_domain_t top() { 
          return boxes_domain_t (lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man())));
        }
        
        static boxes_domain_t bottom() {
          return boxes_domain_t (lddPtr (get_ldd_man(), Ldd_GetFalse (get_ldd_man())));
        }
        
        boxes_domain (const boxes_domain_t& other): 
            ikos::writeable(), m_ldd (other.m_ldd) { 
          crab::CrabStats::count (getDomainName () + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
        }
        
        boxes_domain_t& operator=(const boxes_domain_t& other) {
          crab::CrabStats::count (getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
          if (this != &other) 
            m_ldd = other.m_ldd;
          return *this;
        }
        
        bool is_bottom() const { 
          return &*m_ldd == Ldd_GetFalse (get_ldd_man());
        }
        
        bool is_top() const { 
          return &*m_ldd == Ldd_GetTrue (get_ldd_man());
        }
        
        bool operator<=(boxes_domain_t other) {
          crab::CrabStats::count (getDomainName() + ".count.leq");
          crab::ScopedCrabStats __st__(getDomainName() + ".leq");

          bool res = Ldd_TermLeq (get_ldd_man(), &(*m_ldd), &(*other.m_ldd));

          CRAB_LOG ("boxes", 
                    crab::outs() << "Check if " <<  *this << " <= " <<  other 
                              <<  " ---> " <<  res <<"\n";);
          return res;
        }

        void operator|=(boxes_domain_t other) {
          *this = *this | other;
        }
        
        boxes_domain_t operator|(boxes_domain_t other) {
          crab::CrabStats::count (getDomainName() + ".count.join");
          crab::ScopedCrabStats __st__(getDomainName() + ".join");

          return boxes_domain_t (join (m_ldd, other.m_ldd));
        }
        
        boxes_domain_t operator&(boxes_domain_t other) {
          crab::CrabStats::count (getDomainName() + ".count.meet");
          crab::ScopedCrabStats __st__(getDomainName() + ".meet");

          return boxes_domain_t (lddPtr (get_ldd_man(), 
                                         Ldd_And (get_ldd_man(), &*m_ldd, &*other.m_ldd)));
        }

        boxes_domain_t operator||(boxes_domain_t other) {
          crab::CrabStats::count (getDomainName() + ".count.widening");
          crab::ScopedCrabStats __st__(getDomainName() + ".widening");

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

          CRAB_LOG ("boxes",
                    crab::outs() << "Widening " <<  *this << " and " <<  other 
                              <<  "=" << res <<"\n";);
          return res;
        }

        template<typename Thresholds>
        boxes_domain_t widening_thresholds (boxes_domain_t other, 
                                            const Thresholds & /*ts*/) {
          CRAB_WARN (" boxes widening operator with thresholds not implemented");
          return (*this || other);
        }
        
        boxes_domain_t operator&& (boxes_domain_t other) {
          crab::CrabStats::count (getDomainName() + ".count.narrowing");
          crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

          boxes_domain_t res (*this & other);
          CRAB_WARN (" boxes narrowing operator replaced with meet");
          
          CRAB_LOG ("boxes", 
                    crab::outs() << "Narrowing " << *this << " and " << other << "=" << res <<"\n";);
          return res;
        }
        
        void operator-=(VariableName var) {
          crab::CrabStats::count (getDomainName() + ".count.forget");
          crab::ScopedCrabStats __st__(getDomainName() + ".forget");

          if (is_bottom ()) return;

          int id = get_var_dim (var);
          m_ldd =  lddPtr (get_ldd_man(), 
                           Ldd_ExistsAbstract (get_ldd_man(), &*m_ldd, id));
        }

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          if (is_bottom ()) return;

          for (auto v: boost::make_iterator_range (begin, end))
            operator-= (v);
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          crab::CrabStats::count (getDomainName() + ".count.project");
          crab::ScopedCrabStats __st__(getDomainName() + ".project");

          if (is_bottom ()) return;

          std::set<VariableName> s1,s2,s3;
          for (auto p: (*m_var_map).left) s1.insert (p.first);
          s2.insert (begin, end);
          boost::set_difference (s1,s2,std::inserter (s3, s3.end ()));
          forget (s3.begin (), s3.end ());
        }

        void operator+= (linear_constraint_t cst)
        {
          crab::CrabStats::count (getDomainName() + ".count.add_constraints");
          crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

          if (is_bottom () || cst.is_tautology ())  
            return;

          if (cst.is_contradiction ())  {
            m_ldd = lddPtr (get_ldd_man(), Ldd_GetFalse (get_ldd_man()));
            return;
          }

          linear_expression_t exp = cst.expression();    
          unsigned int size = exp.size ();
          if (size == 0) 
            return; // this should not happen
          else if (size == 1) {
            auto it = exp.begin();
            Number cx = it->first;
            if (cx == 1 || cx == -1) {
              Number k = -exp.constant(); 
              variable_t x = it->second;      
              add_unit_constraint (cx, x, cst.kind (), k);
            }
            else 
              CRAB_WARN (" boxes only supports constraints with unit coefficients");
          }
          else if (size >= 2) {
            // TODO: we can always do the operation in the interval
            //       domain and meet the result with the ldd. But it
            //       might be expensive if a basic block has two many
            //       assume's.
            CRAB_WARN (" boxes only supports constraints with at most one variable.");
          }
          
          CRAB_LOG("boxes", crab::outs() << "Assume(" << cst << ") --> " <<  *this <<"\n";);
        }    

        void operator += (linear_constraint_system_t csts) {
          if (is_bottom ()) return;

          for (auto cst : csts) operator += (cst);
        }

        void set (VariableName v, interval_t ival) {
          crab::CrabStats::count (getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");
          
          if (is_bottom ()) return ;

          constant_t kmin = NULL, kmax = NULL;       
          if (boost::optional <Number> l = ival.lb ().number ())
            kmin = mkCst (*l);
          if (boost::optional <Number> u = ival.ub ().number ())
            kmax = mkCst (*u);
          
          linterm_t t = termForVal(v);
          // FIXME: Ldd_TermReplace seems to leak memory
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, kmin, kmax));
              
          if (kmin) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmin);
          if (kmax) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmax);
          Ldd_GetTheory (get_ldd_man())->destroy_term(t);
        }

        // FIXME: expensive operation
        interval_t operator[](VariableName v) { 

          crab::CrabStats::count (getDomainName() + ".count.to_intervals");
          crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");

          if (is_bottom ()) 
            return interval_t::bottom ();

          if (is_top ()) 
            return interval_t::top ();

          // make convex the ldd
          LddNodePtr tmp = convex_approx ();
          // forget any variable that is not v
          for (auto p: (*m_var_map).left) {
            if (! (p.first == v)) {
              int id = get_var_dim (p.first);
              tmp =  lddPtr (get_ldd_man(), 
                             Ldd_ExistsAbstract (get_ldd_man(), &*tmp, id));
            }
          }
          // convert to interval domain
          linear_constraint_system_t csts;
          toLinCstSysRecur (&*tmp, csts);
          interval_domain_t intv = interval_domain_t::top ();
          for (auto cst: csts) { intv += cst; }
          // do projection with intervals
          return intv[v];
        }
                
        void assign (VariableName x, linear_expression_t e) {
          crab::CrabStats::count (getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");

          if (is_bottom ()) 
            return;

          if (e.is_constant ()) {
            constant_t c = mkCst (e.constant ());
            linterm_t t = termForVal (x);
            m_ldd = lddPtr(get_ldd_man(), 
                           Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, c, c));
            Ldd_GetTheory (get_ldd_man())->destroy_cst (c);
            Ldd_GetTheory (get_ldd_man())->destroy_term(t);
          }
          else if (optional<variable_t> v = e.get_variable()){
            VariableName y = (*v).name();
            if (!(x==y))
              apply (x, y, 1, 0);      
          }
          else {
            CRAB_WARN(" boxes only supports cst or var on the rhs of assignment");
            *this -= x;
          }
          CRAB_LOG("boxes", crab::outs() << "---" << x << ":=" << e << "\n" << *this <<"\n";);
        }
        
        void apply (operation_t op, VariableName x, VariableName y, Number k) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          switch(op){
            case OP_ADDITION: apply (x, y, 1, k); break;
            case OP_SUBTRACTION: apply (x, y, 1, -k); break;
            case OP_MULTIPLICATION: apply (x, y, k, 0); break;
            case OP_DIVISION: 
              {
                // Convert to intervals and perform the operation
                interval_t yi = operator[](y);
                interval_t zi (k);
                interval_t xi = yi / zi;
                set (x, xi);              
                break;
              }
            default: CRAB_ERROR ("Boxes unreachable");
          }
          
          CRAB_LOG("boxes", 
                   crab::outs() << "apply " << x << " := " << y << " " <<  op 
                             << " " << k << " --- " <<  *this <<"\n";);
        }

        void apply(operation_t op, VariableName x, Number k) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ())
            return;

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
              {
                // Convert to intervals and perform the operation
                interval_t yi = operator[](x);
                interval_t zi (k);
                interval_t xi = yi / zi;
                set(x, xi);              
                break;
              }
          }
          CRAB_LOG("boxes",
                   crab::outs() << "apply " << x << " := " << x << " " << op << " " 
                                << k <<  "---" <<  *this <<"\n";);
        }

        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          interval_t zi = operator[](z);
          if (auto k = zi.singleton ())
            apply (op, x, y, *k);            
          else {
            // Convert to intervals and perform the operation
            interval_t yi = operator[](y);
            interval_t xi = interval_t::bottom();
            switch (op) {
              case OP_ADDITION: xi = yi + zi; break;
              case OP_SUBTRACTION: xi = yi - zi; break;
              case OP_MULTIPLICATION: xi = yi * zi; break;
              case OP_DIVISION: xi = yi / zi; break;
            }            
            set(x, xi);
          }
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
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi = operator[](z);
          interval_t xi = interval_t::bottom();
          switch (op) {
            case OP_AND:  xi = yi.And(zi); break;
            case OP_OR:   xi = yi.Or(zi) ; break;
            case OP_XOR:  xi = yi.Xor(zi); break;
            case OP_SHL:  xi = yi.Shl(zi); break;
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
          }
          set(x, xi);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi(k);
          interval_t xi = interval_t::bottom();
          switch (op) {
            case OP_AND:  xi = yi.And(zi); break;
            case OP_OR:   xi = yi.Or(zi) ; break;
            case OP_XOR:  xi = yi.Xor(zi); break;
            case OP_SHL:  xi = yi.Shl(zi); break;
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
          }
          set(x, xi);
        }
        
        // division_operators_api
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          if (op == OP_SDIV){
            apply(OP_DIVISION, x, y, z);
          }
          else{
            // Convert to intervals and perform the operation
            interval_t yi = operator[](y);
            interval_t zi = operator[](z);
            interval_t xi = interval_t::bottom();
            switch (op) {
              case OP_UDIV: xi = yi.UDiv(zi); break;
              case OP_SREM: xi = yi.SRem(zi); break;
              case OP_UREM: xi = yi.URem(zi); break;
              default: CRAB_ERROR("Boxes: unreachable");
            }
            set(x, xi);
          }
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          crab::CrabStats::count (getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom ()) 
            return;

          if (op == OP_SDIV){
            apply(OP_DIVISION, x, y, k);
          }
          else{
            // Convert to intervals and perform the operation
            interval_t yi = operator[](y);
            interval_t zi(k);
            interval_t xi = interval_t::bottom();
            switch (op) {
              case OP_UDIV: xi = yi.UDiv(zi); break;
              case OP_SREM: xi = yi.SRem(zi); break;
              case OP_UREM: xi = yi.URem(zi); break;
              default: CRAB_ERROR("Boxes: unreachable");
            }
            set(x, xi);
          }
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
          LddNodePtr v = convex_approx ();

          // --- extract linear inequalities from the convex ldd
          toLinCstSysRecur(&*v, csts);
          
          return csts;
        }

        
        void write (crab_os& o) {
          // TODO: write to o rather than stdout
          LddManager *ldd_man = getLddManager (m_ldd);
          DdManager *cudd = Ldd_GetCudd (ldd_man);
          FILE *fp = Cudd_ReadStdout(cudd);
          Cudd_SetStdout(cudd, stdout);
          if (m_ldd.get () == Ldd_GetTrue (ldd_man)) 
            o << "{}";
          else if (m_ldd.get () == Ldd_GetFalse (ldd_man)) 
            o << "_|_";	
          else 
          {
            // -- build dictionary
            vector<char*> vnames;
            vnames.reserve (num_of_vars ());
            for (unsigned int i=0; i < num_of_vars (); i++) {
              const char * name = getVarName (i).str ().c_str();
              char * cname = (char*) malloc (sizeof(char) * (std::strlen (name)));
              std::strcpy (cname, name);
              vnames.push_back (cname);
            }

            Ldd_PrintMintermSmtLibv1 (ldd_man, m_ldd.get (), &vnames[0]);

            // -- destroy dict
            for (unsigned int i=0; i < num_of_vars (); i++)
              free (vnames[i]);
          }
          Cudd_SetStdout(cudd,fp);      
        }
        
        static std::string getDomainName () {
          return "Boxes";
        }  
        
      }; 

     template<typename N, typename V, size_t S>
     LddManager* boxes_domain<N,V,S>::m_ldd_man = nullptr;

     template<typename N, typename V, size_t S>
     typename boxes_domain<N,V,S>::var_map_ptr boxes_domain<N,V,S>::m_var_map = nullptr;

     template<typename Number, typename VariableName, size_t LddSize>
     class domain_traits <boxes_domain<Number,VariableName, LddSize> > {
      public:
       typedef boxes_domain<Number, VariableName, LddSize> boxes_domain_t;

       static void normalize (boxes_domain_t& inv) { }

       template <typename Iter>
       static void forget (boxes_domain_t& inv, Iter it, Iter end) {
         inv.forget (it, end);
       }
       
       template <typename Iter>
       static void project (boxes_domain_t& inv, Iter it, Iter end) {
         inv.project (it, end);
       }

       static void expand (boxes_domain_t& inv, VariableName x, VariableName new_x) {
         CRAB_WARN ("boxes expand operation not implemented");
       }

     };
   } // namespace domains
}// namespace crab
#endif /* HAVE_LDD */
#endif 
