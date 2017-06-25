#ifndef BOXES_DOMAIN_HPP
#define BOXES_DOMAIN_HPP

/************************************************************************** 
 * A wrapper for a disjunctive domain of intervals based on "Boxes: A
 * Symbolic Abstract Domain of Boxes" by A. Gurfinkel and S. Chaki
 * published in SAS'10.
 **************************************************************************/

#include <crab/config.h>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/operators_api.hpp>

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
         public pointer_operators< Number, VariableName>,
	 public boolean_operators< Number, VariableName> {
              
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
	typedef disjunctive_linear_constraint_system<Number, VariableName>
	disjunctive_linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;
        typedef boxes_domain <Number, VariableName, LddSize> boxes_domain_t;
        typedef interval <Number> interval_t;

        boxes_domain(): ikos::writeable() { }    

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

        disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system ()
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
#include <boost/optional.hpp>

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
      template<typename Number, typename VariableName, size_t LddSize>
      class boxes_domain_: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName >,
         public array_operators< Number, VariableName>,
         public pointer_operators< Number, VariableName>,
	 public boolean_operators< Number, VariableName>{

        typedef interval_domain <Number, VariableName> interval_domain_t;
       public:
        using typename numerical_domain< Number, VariableName>::linear_expression_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_t;
        using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
	typedef disjunctive_linear_constraint_system<Number, VariableName>
	disjunctive_linear_constraint_system_t;
        using typename numerical_domain< Number, VariableName>::variable_t;
        using typename numerical_domain< Number, VariableName>::number_t;
        using typename numerical_domain< Number, VariableName>::varname_t;

        typedef boxes_domain_ <Number, VariableName, LddSize> boxes_domain_t;
        typedef interval <Number> interval_t;

       private:
        // --- map from crab variable index to ldd term index
        typedef boost::bimap< VariableName , int > var_map_t;
        typedef typename var_map_t::value_type binding_t;

        LddNodePtr m_ldd;
        static LddManager* m_ldd_man;
        static var_map_t m_var_map;

        static LddManager* get_ldd_man () {
          if (!m_ldd_man) {
            DdManager* cudd = Cudd_Init (0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
            theory_t* theory = ldd::create_box_theory<Number> (LddSize);
            CRAB_LOG ("boxes",crab::outs() << "Created a ldd of size " << LddSize <<"\n";);
            m_ldd_man = Ldd_Init (cudd, theory);
            //Cudd_AutodynEnable (cudd, CUDD_REORDER_GROUP_SIFT);
	    Ldd_SanityCheck (get_ldd_man ());
          }
          return m_ldd_man;
        }

        static int num_of_vars () {
          return m_var_map.left.size ();
        }

        int get_var_dim (VariableName v) const {
          auto it = m_var_map.left.find (v);
          if (it != m_var_map.left.end ()) {
            return it->second;
          } else {
	    // XXX: reserved dim 0 for SPECIAL variable
            unsigned int id = m_var_map.size () + 1;
            if (id >= LddSize) {
              CRAB_ERROR ("The Ldd size of ", LddSize, " needs to be larger");
            }
            m_var_map.insert (binding_t (v, id));
            return id;
          }
        }

        inline constant_t mk_cst (ikos::z_number k) {
          mpq_class kk ((mpz_class) k); 
          return (constant_t) tvpi_create_cst (kk.get_mpq_t ());
        }

        inline constant_t mk_cst (ikos::q_number k) {
          mpq_class kk ((mpq_class) k); 
          return (constant_t) tvpi_create_cst (kk.get_mpq_t ());
        }
	
        // convex approximation
        LddNodePtr convex_approx () const {
          return lddPtr (get_ldd_man(), Ldd_TermMinmaxApprox(get_ldd_man(), &*m_ldd));
        }

        // pre: ldd is not either true or false
        void project (LddNodePtr& ldd, VariableName v) const {
          for (auto p: m_var_map.left) {
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
          for (auto p: m_var_map.left) {
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
        linterm_t term_from_var(VariableName v, bool neg = false) 
        {
          int dim = get_var_dim (v);
          int sgn = neg ? -1 : 1;
          linterm_t term = 
              Ldd_GetTheory (get_ldd_man())->create_linterm_sparse_si (&dim, &sgn, 1);
          return term; 
        }

        /** return term from SPECIAL variable $0 */
        linterm_t term_from_special_var(bool neg = false) 
        {
          int dim = 0; // reserved for SPECIAL variable
	  int sgn = neg ? -1 : 1;
          return  Ldd_GetTheory (get_ldd_man())->create_linterm_sparse_si (&dim, &sgn, 1);
        }

	void copy_term (VariableName v, boost::optional<VariableName> x)
	{
          if (is_top () || is_bottom ()) return ;

	  this->operator-=(v); // remove v before assigning new term
	  
          linterm_t lhs = term_from_var (v);
          linterm_t rhs;
	  
	  if (x)
	    rhs = term_from_var (*x);
	  else
	    rhs = term_from_special_var ();
	  
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermCopy(get_ldd_man(), &(*m_ldd), lhs, rhs));
          Ldd_GetTheory (get_ldd_man())->destroy_term(lhs);
          Ldd_GetTheory (get_ldd_man())->destroy_term(rhs);

	  if (!x) {
	    // XXX: if we copy the SPECIAL var to v we forget the SPECIAL
	    //      var after we copy.
	    crab::CrabStats::count (getDomainName() + ".count.forget");
	    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
	    int dim = 0; // SPECIAL variable is dim 0
	    m_ldd =  lddPtr (get_ldd_man(), 
			     Ldd_ExistsAbstract (get_ldd_man(), &*m_ldd, dim));
	  }
	  
	}

        /** v := a * x + k, where a, k are constants and x variable */
        void apply (VariableName v, VariableName x, Number a, Number k) {
          if (is_top () || is_bottom ()) return ;

          linterm_t t = term_from_var (v);
          linterm_t r = term_from_var (x);

          constant_t c1 = mk_cst (a);
          constant_t c2 = mk_cst (k);
          // FIXME: Ldd_TermReplace seems to leak memory
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, r, c1, c2, c2));

          Ldd_GetTheory (get_ldd_man())->destroy_term(t);
          Ldd_GetTheory (get_ldd_man())->destroy_term(r);              
          Ldd_GetTheory (get_ldd_man())->destroy_cst(c1);
          Ldd_GetTheory (get_ldd_man())->destroy_cst(c2);
        }

        void num_from_ldd_cst (constant_t cst, theory_t *theory, ikos::z_number& res) {
          mpq_class v;
          // XXX We know that the theory is tvpi, use its method direclty.
          tvpi_cst_set_mpq (v.get_mpq_t (), (tvpi_cst_t) cst);
          res = ikos::z_number (static_cast<mpz_class> (v));
        }

        void num_from_ldd_cst (constant_t cst, theory_t *theory, ikos::q_number& res) {
          mpq_class v;
          // XXX We know that the theory is tvpi, use its method direclty.
          tvpi_cst_set_mpq (v.get_mpq_t (), (tvpi_cst_t) cst);
	  res = ikos::q_number(v);
        }
	
        linear_expression_t expr_from_ldd_term (linterm_t term, theory_t *theory) {
          linear_expression_t e (0);
          for(size_t i = 0;i < (size_t) theory->term_size (term); i++) {
	    Number k (0); // any value
            num_from_ldd_cst (theory->term_get_coeff (term,i), theory, k);
            VariableName v =  getVarName (theory->term_get_var(term,i));
            e = e + (k * linear_expression_t (v));
          }
          return e;
        }

        linear_constraint_t cst_from_ldd_cons(lincons_t lincons, theory_t *theory) {

          linear_expression_t lhs = expr_from_ldd_term(theory->get_term(lincons), theory);
	  Number rhs (0); // any value
          num_from_ldd_cst(theory->get_constant(lincons), theory, rhs);

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
        void to_lin_cst_sys_recur(LddNode* n, 
				  linear_constraint_system_t& csts) {

          LddNode *N = Ldd_Regular (n);
          if (N == Ldd_GetTrue (get_ldd_man ())) return;


          auto cst = cst_from_ldd_cons (Ldd_GetCons(get_ldd_man (), N), 
                                     Ldd_GetTheory (get_ldd_man ()));

          if (Ldd_Regular (Ldd_T(N)) == Ldd_GetTrue (get_ldd_man ()) &&
              Ldd_Regular (Ldd_E(N)) == Ldd_GetTrue (get_ldd_man ()))
            csts += cst;
          else if (Ldd_Regular (Ldd_T(N)) == Ldd_GetTrue (get_ldd_man ()))
            csts += cst.negate ();
          else
            csts += cst;
            
          to_lin_cst_sys_recur (Ldd_T (N), csts);
          to_lin_cst_sys_recur (Ldd_E (N), csts);

        }

        // Restricts the ldd by the constraint e
        // where e can be one of 
        //    c*x <= k | c*x == k | c*x != k 
        //    and c is 1 or -1 and k is an integer constant
        typedef typename linear_constraint_t::kind_t kind_t;

        LddNodePtr gen_unit_constraint (Number coef, linterm_t term, kind_t kind, Number k) {

          assert (coef == 1 || coef == -1);
          
          auto theory = Ldd_GetTheory (get_ldd_man());
          if (kind == kind_t::EQUALITY) {  // x == k <-> x<=k and !(x<k)
            // x<=k  
            constant_t c1 = mk_cst (k);
	    // XXX: make copy of term so that we can free memory using
	    // destroy_lincons without having double free
	    linterm_t copy_term = Ldd_GetTheory (get_ldd_man())->dup_term (term);
            lincons_t cons1 = theory->create_cons (copy_term, 0 /*non-strict*/, c1);
            LddNodePtr n1 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons1));
            // !(x<k)
            constant_t c2 = mk_cst (k);
            lincons_t cons2 = theory->create_cons (term, 1 /*strict*/, c2);
            LddNodePtr n2 = lddPtr(get_ldd_man(),Ldd_Not (theory->to_ldd (get_ldd_man(), cons2)));
            theory->destroy_lincons (cons1);
            theory->destroy_lincons (cons2);
	    return lddPtr (get_ldd_man (), Ldd_And (get_ldd_man(), &*n1, &*n2));
          } else if (kind == kind_t::INEQUALITY) {
            // case 1:  x <= k  
            // case 2: -x <= k
            constant_t c = mk_cst (k);
            lincons_t cons = theory->create_cons (term, 0 /*non-strict*/, c);
            LddNodePtr n = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons));
            theory->destroy_lincons (cons);
	    return n; 
          } else { // assert (kind == kind_t::DISEQUALITY)
            // case 1:  x != k  <->  x < k OR  x > k <->  x < k OR -x < -k
            // case 2: -x != k  <-> -x < k OR -x > k <-> -x < k OR  x < -k
            constant_t c1 = mk_cst (k);
	    // XXX: make copy of term so that we can free memory using
	    // destroy_lincons without having double free
	    linterm_t copy_term = Ldd_GetTheory (get_ldd_man())->dup_term (term);
            lincons_t cons1 = theory->create_cons (copy_term, 1 /*strict*/, c1);
            LddNodePtr n1 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons1));
            constant_t c2 = mk_cst (-k);
            lincons_t cons2 = theory->create_cons (term, 1 /*strict*/, c2);
            LddNodePtr n2 = lddPtr (get_ldd_man(), theory->to_ldd (get_ldd_man(), cons2));
            LddNodePtr n3 = lddPtr (get_ldd_man (), Ldd_Or (get_ldd_man(), &*n1, &*n2));
            theory->destroy_lincons (cons1);
            theory->destroy_lincons (cons2);            
	    return n3;
          }
        } 

	LddNodePtr gen_unit_constraint (Number coef, variable_t var, kind_t kind, Number k) {
	  linterm_t term = term_from_var (var.name (), (coef == 1 ? false : true));
	  return gen_unit_constraint(coef, term, kind, k);
	}

	LddNodePtr gen_unit_constraint (Number coef, kind_t kind, Number k) {
	  linterm_t term = term_from_special_var ((coef == 1 ? false : true));
	  return gen_unit_constraint(coef, term, kind, k);
	}
	
	
        void add_unit_constraint (Number coef, variable_t x, kind_t kind, Number k) {
	  LddNodePtr n = gen_unit_constraint (coef, x, kind, k);
	  m_ldd = lddPtr (get_ldd_man(), Ldd_And (get_ldd_man(), &*m_ldd, &*n));
        } 

        boxes_domain_ (LddNodePtr ldd): m_ldd (ldd) {
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
                CRAB_WARN ("ldd is growing too fast (number of paths= " ,
			   before_num_paths, "). ", 
                           "After making convex the ldd ",
                           "the number of paths is ",
			   Ldd_PathSize (NULL, &*m_ldd));
                           
              }
              else {
                CRAB_WARN ("ldd is growing too fast (number of paths= " ,
			   before_num_paths, "). ", 
                           "After deleting relationships between variables ",
                           "the number of paths is ",
			   after_num_paths);
              }
            }
          }
          #endif 
        }

	void to_disjunctive_linear_constraint_system_aux
	(LddManager *ldd, LddNode *n,
	 disjunctive_linear_constraint_system_t &e,
	 std::vector<int> &list) {
	    /**
	     * the latest negative constraint to be printed. 
	     * is NULL if there isn't one 
	     */
	    lincons_t negc = nullptr;
	    
	    DdNode *N = Cudd_Regular (n);
	    
	    if (cuddIsConstant (N)) {
	      /* n == N here implies that n is one */
	      if (n == N) {

		linear_constraint_system_t csts;
		
		/* for each level */
		for (int i = 0; i < ldd->cudd->size; i++) {
		  /* let p be the index of level i */
		  int p = ldd->cudd->invperm [i];
		  /* let v be the value of p */
		  int v = list [p];
		  
		  /* skip don't care */
		  if (v == 2) continue;
		  
		  if (v == 0 && ldd->ddVars [p] != nullptr)
		    {
		      lincons_t c;
		      c = THEORY->negate_cons (ldd->ddVars [p]);
		
		      if (negc != nullptr) {
			/* consider negative constraint if it is not implied
			   by c
			*/
			if (!THEORY->is_stronger_cons (c, negc)) {
			  csts += cst_from_ldd_cons (negc, Ldd_GetTheory (ldd));
			}
			THEORY->destroy_lincons (negc);
		      }
		      
		      /* store the current constraint to be considered later */
		      negc = c;
		      continue;
		    }
		  
		  /* if there is a negative constraint waiting to be
		     considered, conjoin it now
		  */
		  if (negc != nullptr) {
		    csts += cst_from_ldd_cons (negc, Ldd_GetTheory (ldd)); 
		    THEORY->destroy_lincons (negc);
		    negc = nullptr;
		  }

		  /* if v is not a don't care but p does not correspond to a
		   * constraint, consider it as a Boolean variable */
		  if (v != 2 && ldd->ddVars [p] == nullptr)  {
		    // TODO ????
		    //fprintf (stderr, "%sb%d", (v == 0 ? "!" : " "), p); 
		  }
		  /* v is true */
		  else if (v == 1) {
		    csts += cst_from_ldd_cons (ldd->ddVars [p],Ldd_GetTheory (ldd));
		  }
		} // end for
		
		/* if there is a constraint waiting to be considered, do it
		   now */
		if (negc != nullptr) {
		  csts += cst_from_ldd_cons (negc, Ldd_GetTheory (ldd)); 
		  THEORY->destroy_lincons (negc);	    
		  negc = nullptr;
		}
		e += csts;
	      }
	    } 
	    else {
	      DdNode *Nv = Cudd_NotCond (cuddT(N), N != n);
	      DdNode *Nnv = Cudd_NotCond (cuddE(N), N != n);
	      int index = N->index;
	      list[index] = 0;
	      to_disjunctive_linear_constraint_system_aux(ldd, Nnv, e, list);
	      list[index] = 1;
	      to_disjunctive_linear_constraint_system_aux(ldd, Nv, e, list);
	      list[index] = 2;
	    }
	    return;
       }

        LddNodePtr getLdd () { 
          return m_ldd; 
        }

        VariableName getVarName (int v) const {
          auto it = m_var_map.right.find (v);
          if (it != m_var_map.right.end ())
             return it->second;
          else {
             CRAB_ERROR ("Index ", v, " cannot be mapped back to a variable name");
          }
        }
	
       public:

        boxes_domain_()
	  : ikos::writeable(), m_ldd (lddPtr (get_ldd_man(), Ldd_GetTrue (get_ldd_man()))) {}
        
        ~boxes_domain_ () { 
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
        
        boxes_domain_ (const boxes_domain_t& other): 
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

          // CRAB_LOG ("boxes", 
          //           crab::outs() << "Check if " <<  *this << " <= " <<  other 
          //                     <<  " ---> " <<  res <<"\n";);
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
                    crab::outs() << "Narrowing " << *this << " and " << other << "="
		                 << res <<"\n";);
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
          for (auto p: m_var_map.left) s1.insert (p.first);
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
            kmin = mk_cst (*l);
          if (boost::optional <Number> u = ival.ub ().number ())
            kmax = mk_cst (*u);
          
          linterm_t t = term_from_var(v);
          // FIXME: Ldd_TermReplace seems to leak memory
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, kmin, kmax));
              
          if (kmin) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmin);
          if (kmax) Ldd_GetTheory (get_ldd_man())->destroy_cst(kmax);
          Ldd_GetTheory (get_ldd_man())->destroy_term(t);
        }

        // expensive operation
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
          for (auto p: m_var_map.left) {
            if (! (p.first == v)) {
              int id = get_var_dim (p.first);
              tmp =  lddPtr (get_ldd_man(), 
                             Ldd_ExistsAbstract (get_ldd_man(), &*tmp, id));
            }
          }
          // convert to interval domain
          linear_constraint_system_t csts;
          to_lin_cst_sys_recur (&*tmp, csts);
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
            constant_t c = mk_cst(e.constant ());
            linterm_t t = term_from_var (x);
            m_ldd = lddPtr(get_ldd_man(), 
                           Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL, NULL, c, c));
            Ldd_GetTheory (get_ldd_man())->destroy_cst (c);
            Ldd_GetTheory (get_ldd_man())->destroy_term(t);
          } else if (boost::optional<variable_t> v = e.get_variable()){
            VariableName y = (*v).name();
            if (!(x==y)) copy_term(x,y); //apply (x, y, 1, 0);      
          } else {
	    // Convert e to intervals
	    interval_t r = e.constant();
	    for (auto p : e)
	      r += p.first * this->operator[](p.second.name());
	    set(x, r);
          }

          CRAB_LOG("boxes", crab::outs() << "--- " << x << ":=" << e << "\n" << *this <<"\n";);
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
                   crab::outs() << "--- apply " << x << " := " << y << " " <<  op 
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
                   crab::outs() << "--- apply " << x << " := " << x << " " << op << " " 
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
            csts += linear_constraint_t::get_false();
          }
          else if(is_top ()) {
            csts += linear_constraint_t::get_true();
          } else {
	    // --- produce convex approximation
	    LddNodePtr v = convex_approx ();
	    // --- extract linear inequalities from the convex ldd
	    to_lin_cst_sys_recur(&*v, csts);
	  }
          return csts;
        }

	disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system () {
	  LddManager *ldd_man = getLddManager (m_ldd);
	  
	  std::vector<int> list;
	  list.reserve (ldd_man->cudd->size);
	  for (int i=0; i < ldd_man->cudd->size; i++) list.push_back(2);
	  
	  disjunctive_linear_constraint_system_t r;
	  to_disjunctive_linear_constraint_system_aux (ldd_man, m_ldd.get (),
						       r, list);
	  return r;
	}
	
	void write (crab_os& o) {
          if (is_top ()) {
            o << "{}";
	  } else if (is_bottom ())  {
            o << "_|_";	
	  } else  {
	    auto r = to_disjunctive_linear_constraint_system ();
	    o << r;
          }
        }

        static std::string getDomainName ()
	{ return "Boxes"; }

	////////
	//// boolean_operators_api
	////////
	
      private:
	
	// return x >= 1 or $0 >=1 if !x
	inline LddNode* gen_true_var (boost::optional<VariableName> x)
	{ // x >=1 <--> -x <= -1
	  if (x) {
	    auto r = gen_unit_constraint (Number(-1), variable_t(*x),
					  linear_constraint_t::INEQUALITY,
					  Number(-1));
	    return &*r;
	  }
	  else {
	    auto r = gen_unit_constraint (Number(-1), 
					  linear_constraint_t::INEQUALITY,
					  Number(-1));
	    return &*r;
	  }
	}

	// return x <= 0 or $0 <= 0 if !x
	inline LddNode* gen_false_var (boost::optional<VariableName> x)
	{
	  if (x) {
	    auto r = gen_unit_constraint (Number(1), variable_t(*x),
					  linear_constraint_t::INEQUALITY,
					  Number(0));
	    return &*r;
	  } else {
	    auto r = gen_unit_constraint (Number(1), 
					  linear_constraint_t::INEQUALITY,
					  Number(0));
	    return &*r;
	  }
	}

	// return Ite(y>=1 Op z>=1, x>=1, x <= 0)
	LddNode* gen_binary_bool(bool_operation_t op,
				 boost::optional<VariableName> x,
				 VariableName y, VariableName z)
	{
	  switch (op) {
	  case OP_BAND:
	    return Ldd_Ite(get_ldd_man (),
			   Ldd_And(get_ldd_man(), gen_true_var(y), gen_true_var(z)), 
			   gen_true_var(x), gen_false_var(x));
	    break;
	  case OP_BOR:
	    return Ldd_Ite(get_ldd_man (),
			   Ldd_Or(get_ldd_man(), gen_true_var(y), gen_true_var(z)),  
			   gen_true_var(x), gen_false_var(x));
	    break;
	  case OP_BXOR:
	    return Ldd_Ite(get_ldd_man (),
			   Ldd_Xor(get_ldd_man(), gen_true_var(y), gen_true_var(z)),
			   gen_true_var(x), gen_false_var(x));
	    break;
	  default: CRAB_ERROR ("Unknown boolean operator");
	  }
	}
	
      public:
	
	void assign_bool_cst (VariableName lhs, linear_constraint_t cst) override
	{
          crab::CrabStats::count (getDomainName() + ".count.assign_bool_cst");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_cst");
	  
	  if (is_bottom ()) return;
	  
	  // XXX: lhs should not appear in cst so we can remove lhs
	  // without losing precision
	  this->operator-=(lhs);
	  
	  if (cst.is_tautology ()) {
	    // m_ldd &= lhs >= 1	    
	    m_ldd = lddPtr (get_ldd_man(),
			    Ldd_And (get_ldd_man(), &*m_ldd, gen_true_var(lhs)));
	  } else if (cst.is_contradiction ()) {
	    // m_ldd &= lhs <= 0
	    m_ldd = lddPtr (get_ldd_man(),
			    Ldd_And (get_ldd_man(), &*m_ldd, gen_false_var(lhs)));
	  } else {
	    linear_expression_t exp = cst.expression();    
	    unsigned int size = exp.size ();
	    if (size == 0) return; // this should not happen 
	    else if (size == 1) {
	      auto it = exp.begin();
	      Number cx = it->first;
	      if (cx == 1 || cx == -1) {
		Number k = -exp.constant(); 
		variable_t vx = it->second;
		// m_ldd &= ite (cst, lhs >=1, lhs<=0);
		m_ldd = lddPtr (get_ldd_man (),
				Ldd_And (get_ldd_man(), &*m_ldd,
				   Ldd_Ite(get_ldd_man (),
					   &*(gen_unit_constraint (cx, vx, cst.kind (),k)),
					   gen_true_var(lhs), gen_false_var(lhs))));
						 
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
	  }
	  
	  CRAB_LOG("boxes",
		   crab::outs () << lhs << ":=" << cst << "=\n" << *this << "\n");
	}    
	
	void assign_bool_var (VariableName x, VariableName y) override
	{
          crab::CrabStats::count (getDomainName() + ".count.assign_bool_var");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_var");
	  
	  copy_term (x, y);
	  
	  CRAB_LOG("boxes",
		   crab::outs()  << x << ":=" << y << "=\n" << *this << "\n");
	}
	
	void apply_binary_bool(bool_operation_t op, VariableName x,
			       VariableName y, VariableName z) override
	{
          crab::CrabStats::count (getDomainName() + ".count.apply_bin_bool");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply_bin_bool");
	  
	  // XXX: if *lhs is null then it represents the SPECIAL
	  // variable $0.
	  boost::optional<VariableName> lhs; 
	  
	  if (!(x == y) && !(x == z)) {
	    lhs = boost::optional<VariableName> (x);
	    // XXX: x does not appear on the rhs so we can remove it
	    // without losing precision.
	    this->operator-=(x);
	  }
	  
	  m_ldd = lddPtr (get_ldd_man (),
			  Ldd_And(get_ldd_man(), &*m_ldd,
				  gen_binary_bool(op, lhs, y, z)));
	  
	  if ((x == y) || (x == z)) {
	    // XXX: if we are here we added ite(y op z, $0 >=1, $0 <= 0);
	    // so we still need to assign $0 to x:
	    copy_term(x, boost::optional<VariableName>());
	  }

	  CRAB_LOG("boxes",
		   crab::outs () << x << ":=" << y << " " << op << " " << z << "=\n"
     		                 << *this << "\n");
	  
	}
	
	void assume_bool (VariableName x, bool is_negated) override
	{
          crab::CrabStats::count (getDomainName() + ".count.assume_bool");
          crab::ScopedCrabStats __st__(getDomainName() + ".assume_bool");
	  
	  m_ldd = lddPtr (get_ldd_man(),
			  Ldd_And (get_ldd_man(), &*m_ldd,
				   ((is_negated) ?
				    gen_false_var (x): gen_true_var (x))));

	  CRAB_LOG("boxes",
		   if (!is_negated) 
		     crab::outs () << "assume(" << x << ")=" << *this << "\n";
		   else
		     crab::outs () << "assume(not(" << x << "))=" << *this << "\n";); 
	}
        
      }; 

     template<typename N, typename V, size_t S>
     LddManager* boxes_domain_<N,V,S>::m_ldd_man = nullptr;

     template<typename N, typename V, size_t S>
     typename boxes_domain_<N,V,S>::var_map_t boxes_domain_<N,V,S>::m_var_map;

     template<typename Number, typename VariableName, size_t LddSize>
     class domain_traits <boxes_domain_<Number,VariableName, LddSize> > {
      public:
       typedef boxes_domain_<Number, VariableName, LddSize> boxes_domain_t;

       template<class CFG>
       static void do_initialization (CFG cfg) { }

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

    #if 0
     // Without copy-on-write
    template<typename Number, typename VariableName, size_t LddSize=3000>
    using boxes_domain = boxes_domain_<Number,VariableName,LddSize>;     
    #else 
    // Quick wrapper which uses shared references with copy-on-write.
    template<class Number, class VariableName, size_t LddSize=3000>
    class boxes_domain: public writeable,
               public numerical_domain<Number, VariableName >,
               public bitwise_operators<Number,VariableName >,
               public division_operators<Number, VariableName >,
               public array_operators<Number, VariableName >,
	       public pointer_operators<Number, VariableName >,
               public boolean_operators<Number, VariableName > {
      public:
      using typename numerical_domain< Number, VariableName >::linear_expression_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
      typedef disjunctive_linear_constraint_system<Number, VariableName>
      disjunctive_linear_constraint_system_t;
      using typename numerical_domain< Number, VariableName >::variable_t;
      using typename numerical_domain< Number, VariableName >::number_t;
      using typename numerical_domain< Number, VariableName >::varname_t;
      typedef typename linear_constraint_t::kind_t constraint_kind_t;
      typedef interval<Number>  interval_t;
      typedef boxes_domain <Number, VariableName, LddSize> boxes_domain_t;      

    private:
      
      typedef boxes_domain_ <Number, VariableName, LddSize> boxes_impl_t;
      typedef std::shared_ptr<boxes_impl_t> boxes_ref_t;

      boxes_ref_t _ref;
      
      boxes_domain(boxes_ref_t ref) : _ref(ref) { }

      boxes_domain_t create(boxes_impl_t&& t)
      { return std::make_shared<boxes_impl_t>(std::move(t)); }
        
      void detach(void)
      { 
        if(!_ref || !_ref.unique())
          _ref = std::make_shared<boxes_impl_t>(*_ref);
      }
      
      boxes_impl_t& ref(void) { return *_ref; }

      const boxes_impl_t& ref(void) const { return *_ref; }
      
    public:

      boxes_domain(bool is_bottom = false): _ref(nullptr)
      {
	if (is_bottom) 
	  _ref = std::make_shared<boxes_impl_t>(boxes_impl_t::bottom());
	else 
	  _ref = std::make_shared<boxes_impl_t>(boxes_impl_t::top());
      }
      
      static boxes_domain_t top() { return boxes_domain(false); }
    
      static boxes_domain_t bottom() { return boxes_domain(true); }

      boxes_domain(const boxes_domain_t& o)
        : _ref(o._ref) { }
      
      boxes_domain_t& operator=(const boxes_domain_t& o) 
      {
        _ref = o._ref;
        return *this;
      }

      bool is_bottom() { return ref().is_bottom(); }
      bool is_top() { return ref().is_top(); }
      bool operator<=(boxes_domain_t& o) { return ref() <= o.ref(); }
      void operator|=(boxes_domain_t o) { detach(); ref() |= o.ref(); }
      boxes_domain_t operator|(boxes_domain_t o) { return create(ref() | o.ref()); }
      boxes_domain_t operator||(boxes_domain_t o) { return create(ref() || o.ref()); }
      boxes_domain_t operator&(boxes_domain_t o) { return create(ref() & o.ref()); }
      boxes_domain_t operator&&(boxes_domain_t o) { return create(ref() && o.ref()); }

      template<typename Thresholds>
      boxes_domain_t widening_thresholds (boxes_domain_t o, const Thresholds &ts) {
        return create(ref().template widening_thresholds<Thresholds>(o.ref(), ts));
      }
     
      void operator+=(linear_constraint_system_t csts) { detach(); ref() += csts; } 
      void operator-=(VariableName v) { detach(); ref() -= v; }
      interval_t operator[](VariableName x) { return ref()[x]; }
      void set(VariableName x, interval_t intv) { detach(); ref().set(x, intv); }

      template<typename Iterator>
      void forget (Iterator vIt, Iterator vEt)
      { detach(); ref().forget(vIt, vEt); }
      void assign(VariableName x, linear_expression_t e)
      { detach(); ref().assign(x, e); }
      void apply(operation_t op, VariableName x, VariableName y, Number k)
      { detach(); ref().apply(op, x, y, k); }
      void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width)
      { detach(); ref().apply(op, x, y, width); }
      void apply(conv_operation_t op, VariableName x, Number k, unsigned width)
      { detach(); ref().apply(op, x, k, width); }
      void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k)
      { detach(); ref().apply(op, x, y, k); }
      void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z)
      { detach(); ref().apply(op, x, y, z); }
      void apply(operation_t op, VariableName x, VariableName y, VariableName z)
      { detach(); ref().apply(op, x, y, z); }
      void apply(div_operation_t op, VariableName x, VariableName y, VariableName z)
      { detach(); ref().apply(op, x, y, z); }
      void apply(div_operation_t op, VariableName x, VariableName y, Number k)
      { detach(); ref().apply(op, x, y, k); }

      void assign_bool_cst (VariableName x, linear_constraint_t cst)
      { detach(); ref().assign_bool_cst(x,cst);}
      void assign_bool_var (VariableName x, VariableName y)
      { detach(); ref().assign_bool_var(x,y);}
      void apply_binary_bool(bool_operation_t op,VariableName x, VariableName y, VariableName z)
      { detach(); ref().apply_binary_bool(op,x,y,z);}
      void assume_bool (VariableName x, bool is_negated)
      { detach(); ref().assume_bool(x,is_negated);}
      
      template<typename Iterator>
      void project (Iterator vIt, Iterator vEt) { detach(); ref().project(vIt, vEt); }

      void write(crab_os& o) { ref().write(o); }

      linear_constraint_system_t to_linear_constraint_system ()
      { return ref().to_linear_constraint_system(); }

      disjunctive_linear_constraint_system_t
      to_disjunctive_linear_constraint_system ()
      { return ref().to_disjunctive_linear_constraint_system(); }
      
      static std::string getDomainName () { return boxes_impl_t::getDomainName(); }
      
    };
     
    template<typename Number, typename VariableName, size_t LddSize>
    class domain_traits <boxes_domain<Number,VariableName, LddSize> > {
    public:
      typedef boxes_domain<Number, VariableName, LddSize> boxes_domain_t;
      
      template<class CFG>
      static void do_initialization (CFG cfg) {}

      static void normalize (boxes_domain_t& inv) {}

      template <typename Iter>
      static void forget (boxes_domain_t& inv, Iter it, Iter end)
      { inv.forget (it, end); }
      
      template <typename Iter>
      static void project (boxes_domain_t& inv, Iter it, Iter end)
      { inv.project (it, end); }
      
      static void expand (boxes_domain_t& inv, VariableName x, VariableName new_x)
      { CRAB_WARN ("boxes expand operation not implemented"); }
    };
    #endif

     
   } // namespace domains
}// namespace crab
#endif /* HAVE_LDD */
#endif 
