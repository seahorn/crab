#pragma once

/************************************************************************** 
 * A wrapper for a disjunctive domain of intervals from "Boxes: A
 * Symbolic Abstract Domain of Boxes" by A. Gurfinkel and S. Chaki
 * published in SAS'10.
 **************************************************************************/

#include <crab/config.h>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>

#ifndef HAVE_LDD
/*
 * Dummy implementation if ldd not found 
 */
#define LDD_NOT_FOUND "No LDD. Run cmake with -DCRAB_USE_LDD=ON"
namespace crab {
   namespace domains {
      template<typename Number, typename VariableName,int ConvexReduce=-1, size_t LddSize=100>
      class boxes_domain final:
       public abstract_domain<boxes_domain<Number,VariableName,ConvexReduce,LddSize>> {
        typedef boxes_domain<Number, VariableName, ConvexReduce, LddSize> boxes_domain_t;
	typedef abstract_domain<boxes_domain_t> abstract_domain_t;
	
       public:
	
        using typename abstract_domain_t::linear_expression_t;
        using typename abstract_domain_t::linear_constraint_t;
        using typename abstract_domain_t::linear_constraint_system_t;
	using typename abstract_domain_t::disjunctive_linear_constraint_system_t;	
        using typename abstract_domain_t::variable_t;
	using typename abstract_domain_t::variable_vector_t;
	using typename abstract_domain_t::pointer_constraint_t;
	typedef Number number_t;
	typedef VariableName varname_t;	
        typedef interval<number_t> interval_t;

        boxes_domain() {}    
	
        void set_to_top() { CRAB_ERROR(LDD_NOT_FOUND); }

        void set_to_bottom() { CRAB_ERROR(LDD_NOT_FOUND); }

        boxes_domain(const boxes_domain_t& other) {}
        
	bool is_bottom() { CRAB_ERROR(LDD_NOT_FOUND); }

        bool is_top() { CRAB_ERROR(LDD_NOT_FOUND); }

        bool operator<=(boxes_domain_t other) { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void operator|=(boxes_domain_t other)
        { CRAB_ERROR(LDD_NOT_FOUND); }

        boxes_domain_t operator|(boxes_domain_t other)
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        boxes_domain_t operator&(boxes_domain_t other) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        boxes_domain_t operator||(boxes_domain_t other)
        { CRAB_ERROR(LDD_NOT_FOUND); }

        boxes_domain_t widening_thresholds(boxes_domain_t other, 
					   const iterators::thresholds<number_t> &ts)
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        boxes_domain_t operator&&(boxes_domain_t other) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void operator-=(variable_t var) 
        { CRAB_ERROR(LDD_NOT_FOUND); }

        interval_t operator[](variable_t v) 
        { CRAB_ERROR(LDD_NOT_FOUND); }

        void set(variable_t v, interval_t ival) 
        { CRAB_ERROR(LDD_NOT_FOUND); }

        void operator +=(linear_constraint_system_t csts) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void assign(variable_t x, linear_expression_t e) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
          
        void apply(operation_t op, variable_t x, variable_t y, number_t z) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void apply(operation_t op, variable_t x, variable_t y, variable_t z) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void backward_assign(variable_t x, linear_expression_t e,
			      boxes_domain_t invariant) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
          
        void backward_apply(operation_t op,
			     variable_t x, variable_t y, number_t z,
			     boxes_domain_t invariant) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void backward_apply(operation_t op,
			    variable_t x, variable_t y, variable_t z,
			    boxes_domain_t invariant)
        { CRAB_ERROR(LDD_NOT_FOUND); }
	
        void apply(int_conv_operation_t op, variable_t dst, variable_t src) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
                
        void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
        
        void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) 
        { CRAB_ERROR(LDD_NOT_FOUND); }

	void assign_bool_cst(variable_t lhs, linear_constraint_t rhs)
        { CRAB_ERROR(LDD_NOT_FOUND); }
	
	void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs)
        { CRAB_ERROR(LDD_NOT_FOUND); }

	void apply_binary_bool(bool_operation_t op, variable_t x,variable_t y,variable_t z)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void assume_bool(variable_t v, bool is_negated)
        { CRAB_ERROR(LDD_NOT_FOUND); }
	
	void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
				      boxes_domain_t invariant)
        { CRAB_ERROR(LDD_NOT_FOUND); }	  	  
	
	void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
				      boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	  	  	  
	void backward_apply_binary_bool(crab::domains::bool_operation_t op,
					variable_t x,variable_t y,variable_t z,
					boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	// array ops
	void array_init(variable_t a, linear_expression_t elem_size,
			linear_expression_t lb_idx, linear_expression_t ub_idx, 
			linear_expression_t val)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void array_load(variable_t lhs,
			variable_t a, linear_expression_t elem_size,
			linear_expression_t i)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void array_store(variable_t a, linear_expression_t elem_size,
			 linear_expression_t i, linear_expression_t v, 
			 bool is_strong_update)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void array_store(variable_t a_new, variable_t a_old,
			 linear_expression_t elem_size,
			 linear_expression_t i, linear_expression_t v, 
			 bool is_strong_update)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	
	void array_store_range(variable_t a, linear_expression_t elem_size,
			       linear_expression_t i, linear_expression_t j,
			       linear_expression_t v)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void array_store_range(variable_t a_new, variable_t a_old,
			       linear_expression_t elem_size,
			       linear_expression_t i, linear_expression_t j,
			       linear_expression_t v)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void array_assign(variable_t lhs, variable_t rhs)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	// backward array ops
	void backward_array_init(variable_t a, linear_expression_t elem_size,
				 linear_expression_t lb_idx, linear_expression_t ub_idx, 
				 linear_expression_t val, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void backward_array_load(variable_t lhs,
				 variable_t a, linear_expression_t elem_size,
				 linear_expression_t i, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void backward_array_store(variable_t a, linear_expression_t elem_size,
				  linear_expression_t i, linear_expression_t v, 
				  bool is_strong_update, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void backward_array_store(variable_t a_new, variable_t a_old,
				  linear_expression_t elem_size,
				  linear_expression_t i, linear_expression_t v, 
				  bool is_strong_update, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	
	void backward_array_store_range(variable_t a, linear_expression_t elem_size,
					linear_expression_t i, linear_expression_t j,
					linear_expression_t v, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	void backward_array_store_range(variable_t a_new, variable_t a_old,
					linear_expression_t elem_size,
					linear_expression_t i, linear_expression_t j,
					linear_expression_t v, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	
	void backward_array_assign(variable_t lhs, variable_t rhs, boxes_domain_t invariant)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	// pointer ops
	void pointer_load(variable_t lhs, variable_t rhs)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_store(variable_t lhs, variable_t rhs)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	  
	  
	void pointer_mk_obj(variable_t lhs, ikos::index_t address)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_function(variable_t lhs, varname_t func)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_mk_null(variable_t lhs)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_assume(pointer_constraint_t cst)
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
	void pointer_assert(pointer_constraint_t cst)
	{ CRAB_ERROR(LDD_NOT_FOUND); }	
	
	void forget(const variable_vector_t& variables)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	void project(const variable_vector_t& variables)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	void expand(variable_t var, variable_t new_var)
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	void normalize()
	{ CRAB_ERROR(LDD_NOT_FOUND); }

	void minimize()
	{ CRAB_ERROR(LDD_NOT_FOUND); }
	
        linear_constraint_system_t to_linear_constraint_system()
        { CRAB_ERROR(LDD_NOT_FOUND); }

        disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system()
        { CRAB_ERROR(LDD_NOT_FOUND); }
	
        void write(crab_os& o) 
        { CRAB_ERROR(LDD_NOT_FOUND); }
          
        static std::string getDomainName() {
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
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>

#include <algorithm>
#include <boost/bimap.hpp>
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
       *
       * FIXME: Ldd_TermReplace seems to leak memory sometimes.
       */
      template<typename Number, typename VariableName, int ConvexReduce, size_t LddSize>
      class boxes_domain_ final:
       public abstract_domain<boxes_domain_<Number,VariableName,ConvexReduce,LddSize>> {
	
        typedef interval_domain<Number, VariableName> interval_domain_t;
        typedef boxes_domain_<Number, VariableName, ConvexReduce, LddSize> boxes_domain_t;
	typedef abstract_domain<boxes_domain_t> abstract_domain_t;
	
       public:
	
	typedef Number number_t;
	typedef VariableName varname_t;
        using typename abstract_domain_t::linear_expression_t;
        using typename abstract_domain_t::linear_constraint_t;
        using typename abstract_domain_t::linear_constraint_system_t;
	using typename abstract_domain_t::disjunctive_linear_constraint_system_t;		
        using typename abstract_domain_t::variable_t;
	using typename abstract_domain_t::variable_vector_t;	
	using typename abstract_domain_t::pointer_constraint_t;
        typedef interval<number_t> interval_t;

       private:
	
        typedef typename linear_constraint_t::kind_t kind_t;
	
        // --- map from crab variable index to ldd term index
        typedef boost::bimap<variable_t, int > var_map_t;
        typedef typename var_map_t::value_type binding_t;

        LddNodePtr m_ldd;
        static LddManager* m_ldd_man;
        static var_map_t m_var_map;

	// -- bool reasoning is mostly based on disjunctions so for
	//    efficiency we might want to disable it if precision
	//    gains do not pay off.
	const bool m_bool_reasoning = true;
	
        static LddManager* get_ldd_man() {
          if (!m_ldd_man) {
            DdManager* cudd = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
            theory_t* theory = ldd::create_box_theory<number_t>(LddSize);
            CRAB_LOG("boxes",
		      crab::outs() << "Created a ldd of size " << LddSize <<"\n";);
            m_ldd_man = Ldd_Init(cudd, theory);
            //Cudd_AutodynEnable(cudd, CUDD_REORDER_GROUP_SIFT);
	    Ldd_SanityCheck(get_ldd_man());
          }
          return m_ldd_man;
        }

	inline theory_t* get_theory() {
	  return Ldd_GetTheory(get_ldd_man());
	}
	
        static int num_of_vars() {
          return m_var_map.left.size();
        }

        int get_var_dim(variable_t v) const {
          auto it = m_var_map.left.find(v);
          if (it != m_var_map.left.end()) {
            return it->second;
          } else {
	    // XXX: reserved dim 0 for SPECIAL variable
            unsigned int id = m_var_map.size() + 1;
            if (id >= LddSize) {
              CRAB_ERROR("The Ldd size of ", LddSize, " needs to be larger");
            }
            m_var_map.insert(binding_t(v, id));
            return id;
          }
        }

        inline constant_t mk_cst(ikos::z_number k) {
	  ikos::q_number qk(k);
          return (constant_t) tvpi_create_cst(qk.get_mpq_t());
        }

        inline constant_t mk_cst(ikos::q_number k) {
          return (constant_t) tvpi_create_cst(k.get_mpq_t());
        }
	
        // convex approximation
        void convex_approx() {
          m_ldd = lddPtr(get_ldd_man(),
			 Ldd_TermMinmaxApprox(get_ldd_man(), &*m_ldd));
        }

        LddNodePtr convex_approx(LddNodePtr ldd) const {
          return lddPtr(get_ldd_man(),
			Ldd_TermMinmaxApprox(get_ldd_man(), &*ldd));
        }
	
        // pre: ldd is not either true or false
        void project(LddNodePtr& ldd, variable_t v) const {
          std::vector<int> qvars;
	  // num_of_vars is shared by all ldd's
	  qvars.reserve(num_of_vars()-1); 
          for (auto p: m_var_map.left) 
            if (!(p.first == v))
	      qvars.push_back(get_var_dim(p.first));
	    
  	  ldd =  lddPtr(get_ldd_man(), 
			 Ldd_MvExistAbstract(get_ldd_man(), &*ldd,
					     &qvars[0], qvars.size()));
	}
	
        LddNodePtr join(LddNodePtr v1, LddNodePtr v2) {
          return lddPtr(get_ldd_man(), Ldd_Or(get_ldd_man(), &*v1, &*v2));
        }

        /** return term for variable v, neg for negation of variable */
        linterm_t term_from_var(variable_t v, bool neg = false)  {
          int dim = get_var_dim(v);
          int sgn = neg ? -1 : 1;
          linterm_t term = 
              Ldd_GetTheory(get_ldd_man())->
	    create_linterm_sparse_si(&dim, &sgn, 1);
          return term; 
        }

        /** return term from SPECIAL variable $0 */
        linterm_t term_from_special_var(bool neg = false)  {
          int dim = 0; // reserved for SPECIAL variable
	  int sgn = neg ? -1 : 1;
          return  Ldd_GetTheory(get_ldd_man())->
	    create_linterm_sparse_si(&dim, &sgn, 1);
        }

	void copy_term(variable_t v, boost::optional<variable_t> x) {
          if (is_top() || is_bottom()) return ;

	  this->operator-=(v); // remove v before assigning new term
	  
          linterm_t lhs = term_from_var(v);
          linterm_t rhs;
	  
	  if (x)
	    rhs = term_from_var(*x);
	  else
	    rhs = term_from_special_var();
	  
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermCopy(get_ldd_man(), &(*m_ldd), lhs, rhs));
          Ldd_GetTheory(get_ldd_man())->destroy_term(lhs);
          Ldd_GetTheory(get_ldd_man())->destroy_term(rhs);

	  if (!x) {
	    // XXX: if we copy the SPECIAL var to v we forget the SPECIAL
	    //      var after we copy.
	    crab::CrabStats::count(getDomainName() + ".count.forget");
	    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
	    int dim = 0; // SPECIAL variable is dim 0
	    m_ldd =  lddPtr(get_ldd_man(), 
			     Ldd_ExistsAbstract(get_ldd_man(), &*m_ldd, dim));
	  }
	  
	}

        /** 
	 **  All the expressiveness about numerical operations with boxes 
	 **  is limited to:
	 **     v := a * x + [k.lb(),k.ub()], 
	 **     where a is a constant, k an interval and x variable 
	 **  Each numerical operation is reduced to this form.
	 */
        void apply_ldd(variable_t v, variable_t x, number_t a, interval_t k) {
          if (is_top() || is_bottom()) return ;

          linterm_t t = term_from_var(v);
          linterm_t r = term_from_var(x);

          constant_t c = mk_cst(a);

	  constant_t kmin = NULL;
	  constant_t kmax = NULL;
	  if (k.lb().is_finite())
	    kmin = mk_cst(*(k.lb().number()));
	  if (k.ub().is_finite())
	    kmax = mk_cst(*(k.ub().number()));
	  
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermReplace(get_ldd_man(), &(*m_ldd),
					 t, r, c, kmin, kmax));

          Ldd_GetTheory(get_ldd_man())->destroy_term(t);
          Ldd_GetTheory(get_ldd_man())->destroy_term(r);              
          Ldd_GetTheory(get_ldd_man())->destroy_cst(c);
	  if (kmin)
	    Ldd_GetTheory(get_ldd_man())->destroy_cst(kmin);
	  if (kmax)
	    Ldd_GetTheory(get_ldd_man())->destroy_cst(kmax);
        }

	/* Special case of apply_ldd : v := [lb,ub] */
        LddNodePtr apply_interval(variable_t v, interval_t ival) {
	  assert(!is_bottom());

          constant_t kmin = NULL, kmax = NULL;       
          if (boost::optional <number_t> l = ival.lb().number())
            kmin = mk_cst(*l);
          if (boost::optional <number_t> u = ival.ub().number())
            kmax = mk_cst(*u);
          
          linterm_t t = term_from_var(v);
          LddNodePtr res = lddPtr(get_ldd_man(), 
				  Ldd_TermReplace(get_ldd_man(), &(*m_ldd),
						  t, NULL, NULL, kmin, kmax));
          if (kmin) get_theory()->destroy_cst(kmin);
          if (kmax) get_theory()->destroy_cst(kmax);
          get_theory()->destroy_term(t);
	  return res;
        }

        void num_from_ldd_cst(constant_t cst, ikos::z_number& res) {
          // XXX We know that the theory is tvpi, use its method direclty.
	  q_number q;
          tvpi_cst_set_mpq(q.get_mpq_t(),(tvpi_cst_t) cst);
	  res = q.round_to_lower(); 
        }

        void num_from_ldd_cst(constant_t cst, ikos::q_number& res) {
          // XXX We know that the theory is tvpi, use its method direclty.
          tvpi_cst_set_mpq(res.get_mpq_t(),(tvpi_cst_t) cst);
        }
	
        linear_expression_t expr_from_ldd_term(linterm_t term) {
          linear_expression_t e(0);
          for(size_t i = 0;i < (size_t) get_theory()->term_size(term); i++) {
	    number_t k(0); // any value
            num_from_ldd_cst(get_theory()->term_get_coeff(term,i), k);
            variable_t v(getVarName(get_theory()->term_get_var(term,i)));
            e = e + (k * linear_expression_t(v));
          }
          return e;
        }


	linear_constraint<ikos::z_number,varname_t>
	cst_from_ldd_strict_cons(const linear_expression<ikos::z_number,varname_t> &l,
				 const linear_expression<ikos::z_number,varname_t> &r) {
	  // l < r <-> l+1 <= r
	  linear_expression<ikos::z_number,varname_t> e = l + 1;
	  e = e - r;
	  return linear_constraint<ikos::z_number,varname_t>
	    (e, linear_constraint_t::INEQUALITY);
	}

	linear_constraint<ikos::q_number,varname_t>
	cst_from_ldd_strict_cons(const linear_expression<ikos::q_number,varname_t> &l,
				 const linear_expression<ikos::q_number,varname_t> &r) {
	  return linear_constraint<ikos::q_number,varname_t>
	    (l-r, linear_constraint_t::STRICT_INEQUALITY);
	}
	
        linear_constraint_t cst_from_ldd_cons(lincons_t lincons) {

          linear_expression_t lhs =
	    expr_from_ldd_term(get_theory()->get_term(lincons));
	  number_t rhs(0); // any value
          num_from_ldd_cst(get_theory()->get_constant(lincons), rhs);

          if (get_theory()->is_strict(lincons)) {
	    // lhs < rhs
	    return cst_from_ldd_strict_cons(lhs, rhs);
          } else  {
            // lhs <= rhs
            linear_expression_t e = lhs - rhs;
            return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
          }
        }

        // // Pre: n is convex and cannot be bottom
        // void to_lin_cst_sys_recur(LddNode* n, 
	// 			  linear_constraint_system_t& csts) {
        //   LddNode *N = Ldd_Regular(n);
        //   if (N == Ldd_GetTrue(get_ldd_man())) return;
	//   lincons_t lincons = Ldd_GetCons(get_ldd_man(), N);
        //   if (Ldd_Regular(Ldd_T(N)) == Ldd_GetTrue(get_ldd_man()) &&
        //       Ldd_Regular(Ldd_E(N)) == Ldd_GetTrue(get_ldd_man())) {
	//     csts += cst_from_ldd_cons(lincons);
	//   } else if (Ldd_Regular(Ldd_T(N)) == Ldd_GetTrue(get_ldd_man())) {
	//     csts += cst_from_ldd_cons(get_theory()->negate_cons(lincons));
	//   } else {
	//     csts += cst_from_ldd_cons(lincons);
	//   }
        //   to_lin_cst_sys_recur(Ldd_T(N), csts);
        //   to_lin_cst_sys_recur(Ldd_E(N), csts);
        // }

	
        // Create a new ldd representing the constraint e
        // where e can be one of 
        //    c*x <= k | c*x == k | c*x != k 
        //    and c is 1 or -1 and k is an integer constant
        LddNodePtr make_unit_constraint(number_t coef, linterm_t term, kind_t kind,
					number_t k) {

          assert(coef == 1 || coef == -1);
          
          if (kind == kind_t::EQUALITY) {  // x == k <-> x<=k and !(x<k)
            // x<=k  
            constant_t c1 = mk_cst(k);
	    // XXX: make copy of term so that we can free memory using
	    // destroy_lincons without having double free
	    linterm_t copy_term = get_theory()->dup_term(term);
            lincons_t cons1 = get_theory()->create_cons(copy_term, 0 /*non-strict*/, c1);
            LddNodePtr n1 = lddPtr(get_ldd_man(),
				    get_theory()->to_ldd(get_ldd_man(), cons1));
            // !(x<k)
            constant_t c2 = mk_cst(k);
            lincons_t cons2 = get_theory()->create_cons(term, 1 /*strict*/, c2);
            LddNodePtr n2 = lddPtr(get_ldd_man(),
				   Ldd_Not(get_theory()->to_ldd(get_ldd_man(), cons2)));
            get_theory()->destroy_lincons(cons1);
            get_theory()->destroy_lincons(cons2);
	    return lddPtr(get_ldd_man(),
			   Ldd_And(get_ldd_man(), &*n1, &*n2));
          } else if (kind == kind_t::INEQUALITY) {
            // x <= k  
            constant_t c = mk_cst(k);
            lincons_t cons = get_theory()->create_cons(term, 0 /*non-strict*/, c);
            LddNodePtr n = lddPtr(get_ldd_man(),
				   get_theory()->to_ldd(get_ldd_man(), cons));
            get_theory()->destroy_lincons(cons);
	    return n;	    
          } else if (kind == kind_t::STRICT_INEQUALITY) {
            // x < k  
            constant_t c = mk_cst(k);
            lincons_t cons = get_theory()->create_cons(term, 1 /*strict*/, c);
            LddNodePtr n = lddPtr(get_ldd_man(),
	    			   get_theory()->to_ldd(get_ldd_man(), cons));
	    get_theory()->destroy_lincons(cons);
	    return n; 
           } else {
	    assert(kind == kind_t::DISEQUATION);
            return lddPtr(get_ldd_man(),
			  Ldd_Not(&*make_unit_constraint(coef, term, kind_t::EQUALITY, k)));
          }
        } 

	LddNodePtr make_unit_constraint(number_t coef, variable_t var, kind_t kind,
					 number_t k) {
	  linterm_t term = term_from_var(var,(coef == 1 ? false : true));
	  return make_unit_constraint(coef, term, kind, k);
	}

	LddNodePtr make_unit_constraint(number_t coef, kind_t kind, number_t k) {
	  linterm_t term = term_from_special_var((coef == 1 ? false : true));
	  return make_unit_constraint(coef, term, kind, k);
	}
	
	
        void apply_unit_constraint(number_t coef, variable_t x, kind_t kind, number_t k) {
	  LddNodePtr n = make_unit_constraint(coef, x, kind, k);
	  m_ldd = lddPtr(get_ldd_man(), Ldd_And(get_ldd_man(), &*m_ldd, &*n));
        } 

	
	// Extract interval constraints from linear expression exp
	void unitcsts_of_exp(const linear_expression_t& exp, 
			     std::vector<std::pair<variable_t, number_t>>& lbs,
			     std::vector<std::pair<variable_t, number_t>>& ubs) {
	  
	  number_t unbounded_lbcoeff;
	  number_t unbounded_ubcoeff;
	  boost::optional<variable_t> unbounded_lbvar;
	  boost::optional<variable_t> unbounded_ubvar;
	  number_t exp_ub = - (exp.constant());
	  std::vector<std::pair<std::pair<number_t,variable_t>, number_t>> pos_terms;
	  std::vector<std::pair<std::pair<number_t,variable_t>, number_t>> neg_terms;
	  for(auto p : exp) {
	    number_t coeff(p.first);
	    if(coeff > number_t(0)) {
	      variable_t y(p.second);
	      // evaluate the variable in the domain
	      auto y_lb = this->operator[](y).lb();
	      if(y_lb.is_infinite()) {
		if(unbounded_lbvar)
		  return;
		unbounded_lbvar = y;
		unbounded_lbcoeff = coeff;
	      } else {
		number_t ymin(*(y_lb.number()));
		exp_ub -= ymin*coeff;
		pos_terms.push_back({{coeff, y}, ymin}); 
	      }
	    } else {
	      variable_t y(p.second);
	      // evaluate the variable in the domain	    
	      auto y_ub = this->operator[](y).ub(); 
	      if(y_ub.is_infinite()) {
		if(unbounded_ubvar)
		  return;
		unbounded_ubvar = y;
		unbounded_ubcoeff = -(coeff);
	      } else {
		number_t ymax(*(y_ub.number()));
		exp_ub -= ymax*coeff;
		neg_terms.push_back({{-coeff, y}, ymax});
	      }
	    }
	  }
	  
	  if(unbounded_lbvar) {
	    if(!unbounded_ubvar) {
	      // Add bounds for x
	      variable_t x(*unbounded_lbvar);
	      ubs.push_back({x, exp_ub/unbounded_lbcoeff});
	    }
	  } else {
	    if(unbounded_ubvar) {
	      // Bounds for y
	      variable_t y(*unbounded_ubvar);
	      lbs.push_back({y, -exp_ub/unbounded_ubcoeff});
	    } else {
	      for(auto pl : neg_terms)
		lbs.push_back({pl.first.second, -exp_ub/pl.first.first + pl.second});
	      for(auto pu : pos_terms)
		ubs.push_back({pu.first.second, exp_ub/pu.first.first + pu.second});
	    }
	  }
	}
	
	interval_t eval_interval(const linear_expression_t &e) {
	  interval_t r = e.constant();
	  for (auto p : e) r += p.first * operator[](p.second);
	  return r;
	}

	// Given a constraint a1*x1 + ... + an*xn <= k and pivot xi,
        // it computes the interval:
	//  k - intv (a1*x1 + ... + ai-1*xi-1 + ai+1*xi+1 + ... an*xn)
	interval_t compute_residual(const linear_expression_t &e, variable_t pivot) {
	  interval_t residual(-e.constant());
	  for (typename linear_expression_t::const_iterator it = e.begin();
	       it != e.end(); ++it) {
	    variable_t v = it->second;
	    if (!(v == pivot)) {
	      residual = residual - (interval_t(it->first) * this->operator[](v));
	    }
	  }
	  return residual;
	}

	bool add_linear_leq(const linear_expression_t& e) {
	  std::vector<std::pair<variable_t, number_t>> lbs;
	  std::vector<std::pair<variable_t, number_t>> ubs;

	  unitcsts_of_exp(e, lbs, ubs);
	  CRAB_LOG("boxes-extract-leq",
		   for (auto p: lbs) {
		     crab::outs() << "\t" << p.first<< ">="<< p.second <<"\n";
		   }
		   for (auto p: ubs) {
		       crab::outs() << "\t" << p.first<< "<="<< p.second <<"\n";
		   });
		       
	  for (auto p: lbs) {
	    apply_unit_constraint(-number_t(1) /*coef*/, p.first, kind_t::INEQUALITY,
				  -p.second);
	    if (is_bottom()) {
	      return false;
	    }
	  }
	  for (auto p: ubs) {
	    apply_unit_constraint(number_t(1) /*coef*/, p.first, kind_t::INEQUALITY,
				  p.second);
	    if (is_bottom()) {
	      return false;
	    }
	  } 
	  return true;
	}


	bool add_linear_diseq(const linear_expression_t&e) {
	  assert(!is_bottom());

	  for (typename linear_expression_t::const_iterator it = e.begin(); it != e.end();
	       ++it) {
	    variable_t v = it->second;
	    interval_t res_i = compute_residual(e, v) / interval_t(it->first);
	    boost::optional<number_t> k = res_i.singleton();
	    if (k) {
	      linear_constraint_t cst(v != *k);
	      operator+=(cst);
	      if (is_bottom()) return false;
	    } else {
	      // XXX: we lose precision here because res_i is an
	      // interval which already over-approximates the possible
	      // values of the residual.
	    }
	  }
	  return true;
	}      
	
	#ifdef OLD_NON_UNIT_CONSTRAINTS
	interval_t compute_residual(interval_domain_t intervals,
				    const linear_constraint_t &cst, variable_t pivot) {
	  interval_t residual(cst.constant());
	  for (typename linear_constraint_t::const_iterator it = cst.begin(); 
	       it != cst.end(); ++it) {
	    variable_t v = it->second;
	    if (v.index() != pivot.index()) {
	      residual = residual -(interval_t(it->first) * intervals[v]);
	    }
	  }
	  return residual;
	}
	
	void intvcst_from_lin_const(const linear_constraint_t &cst,
				    linear_constraint_system_t &intvcsts) {
	  interval_domain_t intervals = to_intervals(m_ldd);
	  for (typename linear_constraint_t::const_iterator it = cst.begin(); 
	       it != cst.end(); ++it) {
	    number_t c = it->first;
	    variable_t pivot = it->second;
	    interval_t rhs = compute_residual(intervals, cst, pivot);
	    
	    if (!(rhs.lb().is_finite() && rhs.ub().is_finite())) continue;

	    if (auto k = rhs.singleton()) {
	      linear_expression_t term(c * pivot);
	      linear_expression_t e = term -(*k);
	      intvcsts += linear_constraint_t(e, cst.kind());
	    } else {
	      rhs = rhs / interval_t(c);
	      number_t min = *(rhs.lb().number());
	      number_t max = *(rhs.ub().number());
	      
	      switch (cst.kind()) {
	      case kind_t::EQUALITY: {
		intvcsts += linear_constraint_t(pivot >= min);
		intvcsts += linear_constraint_t(pivot <= max);
		break;
	      }
	      case kind_t::INEQUALITY: {
		if (c < 0) 
		  intvcsts += linear_constraint_t(pivot >= min);
		else
		  intvcsts += linear_constraint_t(pivot <= max);
		break;
	      }
	      case kind_t::STRICT_INEQUALITY: {
		if (c < 0)
		  intvcsts += linear_constraint_t(pivot > min);
		else
		  intvcsts += linear_constraint_t(pivot < max);		  
		break;
	      }
	      default:
		assert(cst.kind() == kind_t::DISEQUATION);
		CRAB_WARN("ignored ", cst);
	      }
	    }
	  }
	}
	#endif 
	
        boxes_domain_(LddNodePtr ldd): m_ldd(ldd) {
	  if (ConvexReduce > 0) {
	    // XXX: the value of ConvexReduce is quite arbitrary. A
	    // good value seems around 1000000
	    unsigned threshold = num_of_vars() * ConvexReduce;
	    //unsigned num_paths = Ldd_PathSize(NULL, &*m_ldd);
	    unsigned num_paths = (unsigned) Cudd_CountPath(&*m_ldd);
	    if (threshold > 0 && num_paths > threshold) {
	      convex_approx();                
	      CRAB_WARN("ldd size was too large: ", num_paths, ". Made ldd convex.");
	    }
	  }
        }

	interval_domain_t to_intervals(LddNodePtr &ldd) {
          crab::CrabStats::count(getDomainName() + ".count.to_intervals");
          crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");

          if (&*ldd == Ldd_GetFalse(get_ldd_man())) 
            return interval_domain_t::bottom();
          if (&*ldd == Ldd_GetTrue(get_ldd_man())) 
            return interval_domain_t::top();

	  LddNodePtr ldd_copy(ldd);
	  ldd_copy = convex_approx(ldd_copy);
	  
	  auto disjs = to_disjunctive_linear_constraint_system(ldd_copy);
	  if (disjs.is_true()) {
	    return interval_domain_t::top();
	  } else if(disjs.is_false()) {
	    return interval_domain_t::bottom();	    
	  } else {
	    if (disjs.size() != 1) {
	      CRAB_ERROR("Boxes::to_intervals: it should not be disjunctive ",
			 disjs);
	    }
	    interval_domain_t intv;	  
	    for (auto c : *(disjs.begin())) {
	      intv += c;
	    }
	    return intv;
	  }
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
	    
	    DdNode *N = Cudd_Regular(n);
	    
	    if (cuddIsConstant(N)) {
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
		      c = ldd->theory->negate_cons(ldd->ddVars [p]);
		
		      if (negc != nullptr) {
			/* consider negative constraint if it is not implied
			   by c
			*/
			if (!ldd->theory->is_stronger_cons(c, negc)) {
			  csts += cst_from_ldd_cons(negc);
			}
			ldd->theory->destroy_lincons(negc);
		      }
		      
		      /* store the current constraint to be considered later */
		      negc = c;
		      continue;
		    }
		  
		  /* if there is a negative constraint waiting to be
		     considered, conjoin it now
		  */
		  if (negc != nullptr) {
		    csts += cst_from_ldd_cons(negc); 
		    ldd->theory->destroy_lincons(negc);
		    negc = nullptr;
		  }

		  /* if v is not a don't care but p does not correspond to a
		   * constraint, consider it as a Boolean variable */
		  if (v != 2 && ldd->ddVars [p] == nullptr)  {
		    // XXX: I don't know what to do here
		    //fprintf(stderr, "%sb%d",(v == 0 ? "!" : " "), p); 
		  }
		  /* v is true */
		  else if (v == 1) {
		    csts += cst_from_ldd_cons(ldd->ddVars [p]);
		  }
		} // end for
		
		/* if there is a constraint waiting to be considered, do it
		   now */
		if (negc != nullptr) {
		  csts += cst_from_ldd_cons(negc); 
		  ldd->theory->destroy_lincons(negc);	    
		  negc = nullptr;
		}

		if (csts.is_true()) {
		  // FIXME: if csts is true then e should become true
		} else {
		  e += csts;
		}
	      }
	    } 
	    else {
	      DdNode *Nv = Cudd_NotCond(cuddT(N), N != n);
	      DdNode *Nnv = Cudd_NotCond(cuddE(N), N != n);
	      int index = N->index;
	      list[index] = 0;
	      to_disjunctive_linear_constraint_system_aux(ldd, Nnv, e, list);
	      list[index] = 1;
	      to_disjunctive_linear_constraint_system_aux(ldd, Nv, e, list);
	      list[index] = 2;
	    }
	    return;
       }

	inline LddNode* mk_true(boost::optional<variable_t> x) {
	  if (x) {
	    // x >=1 <--> -x <= -1
	    LddNodePtr  r = make_unit_constraint(number_t(-1), variable_t(*x),
						 linear_constraint_t::INEQUALITY,
						 number_t(-1));
	    return &*r;
	  }
	  else {
	    LddNodePtr  r = make_unit_constraint(number_t(-1), 
						 linear_constraint_t::INEQUALITY,
						 number_t(-1));
	    return &*r;
	  }
	}

	
	inline LddNode* mk_false(boost::optional<variable_t> x) {
	  
	  if (x) {
	    // return x <= 0 
	    LddNodePtr r = make_unit_constraint(number_t(1), variable_t(*x),
						linear_constraint_t::INEQUALITY,
						number_t(0));
	    return &*r;
	  } else {
	    LddNodePtr r = make_unit_constraint(number_t(1), 
						linear_constraint_t::INEQUALITY,
						number_t(0));
	    return &*r;
	  }
	}

	// return Ite(y>=1 Op z>=1, x>= 1, x <= 0)
	LddNodePtr gen_binary_bool(bool_operation_t op,
				   boost::optional<variable_t> x,
				   variable_t y, variable_t z)
	{
	  switch (op) {
	  case OP_BAND: {
	    LddNodePtr c = lddPtr(get_ldd_man(),
				  Ldd_And(get_ldd_man(), mk_true(y), mk_true(z)));
	    return lddPtr(get_ldd_man(),
			  Ldd_Ite(get_ldd_man(),
				  &*c, mk_true(x), mk_false(x)));
	    break;
	  }
	  case OP_BOR: {
	    LddNodePtr c = lddPtr(get_ldd_man(),
				  Ldd_Or(get_ldd_man(), mk_true(y), mk_true(z)));
	    return lddPtr(get_ldd_man(),
			  Ldd_Ite(get_ldd_man(),
				  &*c, mk_true(x), mk_false(x)));
	    break;
	  }
	  case OP_BXOR: {
	    LddNodePtr c = lddPtr(get_ldd_man(),
				  Ldd_Xor(get_ldd_man(), mk_true(y), mk_true(z)));
	    return lddPtr(get_ldd_man(),
			  Ldd_Ite(get_ldd_man(),
				  &*c, mk_true(x), mk_false(x)));
				  
	    break;
	  }
	  default: CRAB_ERROR("Unknown boolean operator");
	  }
	}
	
        variable_t getVarName(int v) const {
          auto it = m_var_map.right.find (v);
          if (it != m_var_map.right.end())
             return it->second;
          else {
             CRAB_ERROR("Index ", v, " cannot be mapped back to a variable name");
          }
        }

	
       public:

        boxes_domain_():
	  m_ldd(lddPtr(get_ldd_man(), Ldd_GetTrue(get_ldd_man()))) {}
        
        ~boxes_domain_() { 
          // DdManager *cudd = nullptr;
          // theory_t *theory = nullptr;
          // if (m_ldd_man)  {
	  //   cudd = Ldd_GetCudd(m_ldd_man);
	  //   theory = Ldd_GetTheory(m_ldd_man);
	  //   Ldd_Quit(m_ldd_man);
	  // }
          // if (theory) tvpi_destroy_theory(theory);
          // if (cudd) Cudd_Quit(cudd);
        }

	
        void set_to_top() { 
	  boxes_domain_t abs(lddPtr(get_ldd_man(),
				    Ldd_GetTrue(get_ldd_man())));
	  std::swap(*this, abs);
        }
        
        void set_to_bottom() {
	  boxes_domain_t abs(lddPtr(get_ldd_man(),
				    Ldd_GetFalse(get_ldd_man())));
	  std::swap(*this, abs);
        }
        
        boxes_domain_(const boxes_domain_t& other): 
	  // m_ldd(other.m_ldd)
	  m_ldd(lddPtr(get_ldd_man(), &(*other.m_ldd))) { 
          crab::CrabStats::count(getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
        }

        boxes_domain_(boxes_domain_t&& other):
	  m_ldd(std::move(other.m_ldd)) { }  
	
        boxes_domain_t& operator=(const boxes_domain_t& other) {
          crab::CrabStats::count(getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
          if (this != &other) {
	    //m_ldd = other.m_ldd;
            m_ldd = lddPtr(get_ldd_man(), &(*other.m_ldd));
	  }
          return *this;
        }
      	
        bool is_bottom() { 
          return &*m_ldd == Ldd_GetFalse(get_ldd_man());
        }
        
        bool is_top() { 
          return &*m_ldd == Ldd_GetTrue(get_ldd_man());
        }
        
        bool operator<=(boxes_domain_t other) {
          crab::CrabStats::count(getDomainName() + ".count.leq");
          crab::ScopedCrabStats __st__(getDomainName() + ".leq");

          bool res = Ldd_TermLeq(get_ldd_man(), &(*m_ldd), &(*other.m_ldd));

          CRAB_LOG("boxes", 
		   crab::outs() << "Check if " <<  *this << " <= " <<  other 
		                <<  " ---> " <<  res <<"\n";);
          return res;
        }

        void operator|=(boxes_domain_t other) {
          *this = *this | other;
        }
        
        boxes_domain_t operator|(boxes_domain_t other) {
          crab::CrabStats::count(getDomainName() + ".count.join");
          crab::ScopedCrabStats __st__(getDomainName() + ".join");

          return boxes_domain_t(join(m_ldd, other.m_ldd));
        }
        
        boxes_domain_t operator&(boxes_domain_t other) {
          crab::CrabStats::count(getDomainName() + ".count.meet");
          crab::ScopedCrabStats __st__(getDomainName() + ".meet");

          return boxes_domain_t(lddPtr(get_ldd_man(), 
                                         Ldd_And(get_ldd_man(),
						  &*m_ldd, &*other.m_ldd)));
        }

        boxes_domain_t operator||(boxes_domain_t other) {
          crab::CrabStats::count(getDomainName() + ".count.widening");
          crab::ScopedCrabStats __st__(getDomainName() + ".widening");

          // It is not necessarily true that the new value is bigger
          // than the old value so we apply 
          // widen(old, new) = widen(old,(join(old,new)))
          LddNodePtr v = join(m_ldd, other.m_ldd);
          LddNodePtr w = lddPtr(get_ldd_man(), 
                                 Ldd_BoxWiden2(get_ldd_man(), &*m_ldd, &*v));

          #if 0
          /** Trick from ufo (needed for SV-COMP ssh programs):
	      ensure that 'w' is only the fronteer of the computation.
              Not sure whether this is still a widening though 
          */
          w = lddPtr(get_ldd_man(), 
                      Ldd_And(get_ldd_man(), &*w, Ldd_Not(&*m_ldd)));
          /** ensure the output is at least as big as newV */
          w = lddPtr(get_ldd_man(), Ldd_Or(get_ldd_man(), &*w, &*other.m_ldd));
	  #endif

          boxes_domain_t res(w); 
	  
          CRAB_LOG("boxes",
                    crab::outs() << "Performed widening \n"
		                 << "**" << *this  << "\n" 
		                 << "** " << other  << "\n" 
                                 << "= " << res <<"\n";);
          return res;
        }
	
        boxes_domain_t widening_thresholds(boxes_domain_t other, 
					   const iterators::thresholds<number_t>& /*ts*/) {
          //CRAB_WARN(" boxes widening operator with thresholds not implemented");
          return(*this || other);
        }

        boxes_domain_t operator&&(boxes_domain_t other) {
          crab::CrabStats::count(getDomainName() + ".count.narrowing");
          crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

          boxes_domain_t res(*this & other);
          //CRAB_WARN(" boxes narrowing operator replaced with meet");
          CRAB_LOG("boxes",
                    crab::outs() << "Performed narrowing \n"
		                 << "**" << *this  << "\n" 
		                 << "** " << other  << "\n" 
                                 << "= " << res <<"\n";);
          return res;
        }

	boxes_domain_t complement() const {
          crab::CrabStats::count(getDomainName() + ".count.complement");
          crab::ScopedCrabStats __st__(getDomainName() + ".complement");
	  LddNodePtr w = lddPtr(get_ldd_man(), Ldd_Not(&*m_ldd));
	  boxes_domain_t res(w);
          CRAB_LOG("boxes",
		    boxes_domain_t tmp(*this);
                    crab::outs() << "Performed complement \n"
		                 << "**" << tmp  << "\n" 
		                 << "= " << res <<"\n";);
	  
	  return res;
	}
	
        void operator-=(variable_t var) {
          crab::CrabStats::count(getDomainName() + ".count.forget");
          crab::ScopedCrabStats __st__(getDomainName() + ".forget");

          if (is_bottom() || is_top()) return;

          int id = get_var_dim(var);
          m_ldd =  lddPtr(get_ldd_man(), 
                           Ldd_ExistsAbstract(get_ldd_man(), &*m_ldd, id));
        }

        void forget(const variable_vector_t& variables) {
          crab::CrabStats::count(getDomainName() + ".count.forget");
          crab::ScopedCrabStats __st__(getDomainName() + ".forget");
	  
          if (is_bottom() || is_top()) return;

	  std::vector<int> qvars;
	  qvars.reserve(variables.size());
          for (variable_t v: variables) {
	    qvars.push_back(get_var_dim(v));
	  }

  	  m_ldd =  lddPtr(get_ldd_man(), 
			   Ldd_MvExistAbstract(get_ldd_man(), &*m_ldd,
						&qvars[0], qvars.size()));
        }

        void project(const variable_vector_t& variables) {
          crab::CrabStats::count(getDomainName() + ".count.project");
          crab::ScopedCrabStats __st__(getDomainName() + ".project");

          if (is_bottom() || is_top()) return;
	  
	  if (variables.empty()) {
	    set_to_top();
	    return;
	  }
	  
          std::set<variable_t> s1,s2;
	  variable_vector_t s3;
          for (auto p: m_var_map.left) s1.insert(p.first);
          s2.insert(variables.begin(), variables.end());
          std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
			      std::back_inserter(s3));
          forget(s3);
        }

	void expand(variable_t v, variable_t new_v) {
          crab::CrabStats::count(getDomainName() + ".count.expand");
          crab::ScopedCrabStats __st__(getDomainName() + ".expand");

          if (is_top() || is_bottom()) return ;
	  
	  // new_v should be completely unconstrained
	  this->operator-=(new_v); 
	  
          linterm_t lnew_v = term_from_var(new_v);
          linterm_t lv = term_from_var(v);
	  
          m_ldd = lddPtr(get_ldd_man(), 
                         Ldd_TermCopy(get_ldd_man(), &(*m_ldd), lnew_v, lv));
          Ldd_GetTheory(get_ldd_man())->destroy_term(lnew_v);
          Ldd_GetTheory(get_ldd_man())->destroy_term(lv);
	}
	  
        void operator+=(linear_constraint_t cst) {
          crab::CrabStats::count(getDomainName() + ".count.add_constraints");
          crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

          CRAB_LOG("boxes",
		   crab::outs() << "--- assume(" << cst << ") --> ";);
	  
          if (is_bottom() || cst.is_tautology()) {  
            return;
	  }

          if (cst.is_contradiction())  {
	    set_to_bottom();
            return;
          }

	  // XXX: we do nothing with unsigned linear inequalities
	  // TODO: we can express these constraints with more disjunctions.
	  if (cst.is_inequality() && cst.is_unsigned()) {
	    CRAB_WARN("unsigned inequality skipped in boxes domain");	  
	    return;
	  }
	  
          linear_expression_t exp = cst.expression();    
          unsigned int size = exp.size();
          if (size == 0) {
	    return; // this should not happen
	  } else if (size == 1) {
            auto it = exp.begin();
            number_t cx = it->first;
            if (cx == number_t(1) || cx == number_t(-1)) {
              number_t k = -exp.constant(); 
	      variable_t x = it->second;      
              apply_unit_constraint(cx, x, cst.kind(), k);
            } else {
	      // XXX: this is not possible when program is originated by clang/llvm
              CRAB_WARN("non-unit coefficients not implemented in boxes");
	    }
          } else if (size >= 2) {
	    #ifdef OLD_NON_UNIT_CONSTRAINTS
	    linear_constraint_system_t intvcsts;
	    intvcst_from_lin_const(cst, intvcsts);
	    this->operator+=(intvcsts);
	    CRAB_LOG("boxes-lincons-to-int",
		     crab::outs() << cst << " converted to interval constraints : "
		                  << intvcsts << "\n";);
            #else
	    if (cst.is_disequation()) {
	      if (!add_linear_diseq(cst.expression())) {
		CRAB_LOG("boxes", crab::outs() << "_|_\n";);
		return;
	      }
	    } else if(cst.is_inequality()) {
	      if (!add_linear_leq(cst.expression())) {
		CRAB_LOG("boxes", crab::outs() << "_|_\n";);
		return;
	      }
	    } else if(cst.is_equality()) {
	      linear_expression_t exp = cst.expression();
	      // split equality into two inequalities
	      if(!add_linear_leq(exp) || !add_linear_leq(-exp)) {
		CRAB_LOG("boxes", crab::outs() << "_|_\n";);		
		return;
	      }
	    }
	    #endif 
          }
          
          CRAB_LOG("boxes", crab::outs() << *this <<"\n";);
        }    

	void normalize() {}

	void minimize() {}
	
        void operator+=(linear_constraint_system_t csts) {
          if (is_bottom()) return;
          for(auto cst : csts) operator += (cst);
        }
	
        void set(variable_t v, interval_t ival) {
          crab::CrabStats::count(getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");
          
          if (is_bottom()) return ;

	  m_ldd = apply_interval(v, ival);
        }

        interval_t operator[](variable_t v) {
          crab::CrabStats::count(getDomainName() + ".count.to_intervals");
          crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");

          if (is_bottom()) 
            return interval_t::bottom();
          if (is_top()) 
            return interval_t::top();
	  
	  LddNodePtr ldd(m_ldd);
	  project(ldd, v);
	  interval_domain_t intv = to_intervals(ldd);
	  interval_t i = intv[v];
	  CRAB_LOG("boxes-project",
		   crab::outs() << "Before projecting on " << v << ": " << *this << "\n"
		                 << "Projection " << i << "\n";);
	  return i;
        }
                
	// x := e
        void assign(variable_t x, linear_expression_t e) {
          crab::CrabStats::count(getDomainName() + ".count.assign");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign");
	  
          if (is_bottom()) return; 

          if (e.is_constant()) {
            constant_t c = mk_cst(e.constant());
            linterm_t t = term_from_var(x);
            m_ldd = lddPtr(get_ldd_man(), 
                           Ldd_TermReplace(get_ldd_man(), &(*m_ldd),
	    				   t, NULL, NULL, c, c));
            get_theory()->destroy_cst(c);
            get_theory()->destroy_term(t);
          } else if (boost::optional<variable_t> v = e.get_variable()){
            variable_t y =(*v);
            if (!(x==y)) {
	      //copy_term(x,y);
	      apply_ldd(x, y, 1, number_t(0));
	    }
          } else {
	    m_ldd = apply_interval(x, eval_interval(e));
          }

          CRAB_LOG("boxes",
		   crab::outs() << "--- " << x << ":=" << e << "\n"
		                << *this <<"\n";);
        }

	// x := y op k
        void apply(operation_t op, variable_t x, variable_t y, number_t k) {
          crab::CrabStats::count(getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom()) {
            return;
	  }
	  
	  if (op >= OP_ADDITION && op <= OP_MULTIPLICATION) {	  
	    switch(op){
            case OP_ADDITION:
	      apply_ldd(x, y, 1, k);
	      break;
            case OP_SUBTRACTION:
	      apply_ldd(x, y, 1, -k);
	      break;
	    case OP_MULTIPLICATION:
	      apply_ldd(x, y, k, number_t(0));
	      break;
	    default:
	      CRAB_ERROR("Unexpected operator ", op);
	    }

	    CRAB_LOG("boxes", 
		     crab::outs() << "--- " << x << " := " << y << " " <<  op 
                                  << " " << k << "\n" <<  *this <<"\n";);
	  } else {
            // Convert to intervals 
            interval_t yi(this->operator[](y));
            interval_t zi(k);
            interval_t xi(interval_t::bottom());
            switch(op) {
	      case OP_SDIV:
		xi = yi / zi;
		break;
              case OP_UDIV:
		xi = yi.UDiv(zi);
		break;
              case OP_SREM:
		xi = yi.SRem(zi);
		break;
              case OP_UREM:
		xi = yi.URem(zi);
		break;
              default:
		CRAB_ERROR("Unexpected operator ", op);
            }
            set(x, xi);
          }
	}

	// x := y op z
        void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
          crab::CrabStats::count(getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");
	  
          if (is_bottom()) {
	    return;
	  }

	  // --- if z is a singleton we do not lose precision
	  interval_t zi = this->operator[](z);
	  if (auto k = zi.singleton()) {
	    apply(op, x, y, *k);
	    return;
	  }

	  // --- if y or z is top then we give up
	  interval_t yi = this->operator[](y);	    
	  if (yi.is_top() || zi.is_top()) {
	    this->operator-=(x);
	    goto apply_end;
	  }

	  // --- if y is a singleton we do not lose precision
	  if (auto k = yi.singleton()) {
	    if (op == OP_ADDITION || op == OP_MULTIPLICATION) {
	      apply(op, x, z, *k);
	      return;
	    } else if(op == OP_SUBTRACTION) {
	      // x = yi - z <--> x = -z + [yi.lb,yi.ub]
	      apply_ldd(x, z, -1, *k);
	      goto apply_end;
	    }
	  } 	  

	  // XXX: ideally we would like to extract all the interval
	  // constraints from x:= y op z, forget x and, then add all
	  // the interval constraints. However, we would need to
	  // project on y and z while keeping all the disjunctive
	  // information about these variables. Note that the method
	  // operator[] loses all the disjunctive information.
	  
	  
	  // -- We need to abstract either y or z and lose some
	  // precision. We abstract the one with the smaller interval.
	  if (op == OP_ADDITION) {
	    if (yi <= zi) { // abstract y
	      // x = yi + z <--> x = z + yi
	      apply_ldd(x, z, 1, yi);
	    } else {        // abstract z
	      // x = y + zi
	      apply_ldd(x, y, 1, zi);
	    }
	    goto apply_end;
	  }

	  // -- We need to abstract either y or z and lose some
	  // precision. We abstract the one with the smaller interval.	  
	  if (op == OP_SUBTRACTION) {
	    if (yi <= zi) { // abstract y
	      // x = yi - z <-->  x = -z + yi
	      apply_ldd(x, z, -1, yi);
	    } else {        // abstract z
	      // x = y - zi
	      apply_ldd(x, y, 1, zi * number_t(-1));
	    }
	    goto apply_end;
	  }

	  // --- we go to intervals and lose precision
	  switch (op) {
	    case OP_MULTIPLICATION: 
	      set(x, yi * zi);
	      break;
	    case OP_SDIV: 
	      set(x, yi / zi);
	      break;
  	    case OP_UDIV:
	      set(x, yi.UDiv(zi));
	      break;
	    case OP_SREM:
	      set(x, yi.SRem(zi));
	      break;
	    case OP_UREM:
	      set(x, yi.URem(zi));
	      break;
	    default:
	      CRAB_ERROR("Unexpected operator ", op);
	  }

	apply_end:
          CRAB_LOG("boxes", 
                   crab::outs() << "--- " << x << " := " << y << " " <<  op 
		                << " " << z << "\n" <<  *this <<"\n";);
	}
	
        void backward_assign(variable_t x, linear_expression_t e,
			      boxes_domain_t invariant) {
	  crab::CrabStats::count(getDomainName() + ".count.backward_assign");
	  crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign");
	  
	  CRAB_LOG("boxes",
		   crab::outs() << "Backward " << x << ":=" << e
 		                 << "\n\tPOST=" << *this << "\n");
	  BackwardAssignOps<boxes_domain_t>::assign(*this, x, e, invariant);
	  CRAB_LOG("boxes",
		   crab::outs() << "\tPRE=" << *this << "\n");
	}

        void backward_apply(operation_t op,
			     variable_t x, variable_t y, number_t z,
			     boxes_domain_t invariant) {
	  crab::CrabStats::count(getDomainName() + ".count.backward_apply");
	  crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");
	  
	  CRAB_LOG("boxes",
		   crab::outs() << "Backward " << x << ":=" << y << op << z 
 		                 << "\n\tPOST=" << *this << "\n");
	  BackwardAssignOps<boxes_domain_t>::apply(*this, op, x, y, z,
						   invariant);
	  CRAB_LOG("boxes",
		   crab::outs() << "\tPRE=" << *this << "\n");
	  
	}

        void backward_apply(operation_t op,
			    variable_t x, variable_t y, variable_t z,
			    boxes_domain_t invariant)  {
	  crab::CrabStats::count(getDomainName() + ".count.backward_apply");
	  crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");
	  
	  CRAB_LOG("boxes",
		   crab::outs() << "Backward " << x << ":=" << y << op << z 
 		                 << "\n\tPOST=" << *this << "\n");
	  BackwardAssignOps<boxes_domain_t>::apply(*this, op, x, y, z,
						   invariant);
	  CRAB_LOG("boxes",
		   crab::outs() << "\tPRE=" << *this << "\n");
	  
	}	
	
        void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
          // since reasoning about infinite precision we simply assign and
          // ignore the widths.
          assign(dst, src);
        }
	
        void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
          crab::CrabStats::count(getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom()) 
            return;

          // Convert to intervals and perform the operation
          interval_t yi = this->operator[](y);
          interval_t zi = this->operator[](z);
          interval_t xi = interval_t::bottom();
          switch(op) {
            case OP_AND:  xi = yi.And(zi); break;
            case OP_OR:   xi = yi.Or(zi) ; break;
            case OP_XOR:  xi = yi.Xor(zi); break;
            case OP_SHL:  xi = yi.Shl(zi); break;
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
          }
          set(x, xi);
        }
        
        void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
          crab::CrabStats::count(getDomainName() + ".count.apply");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply");

          if (is_bottom()) 
            return;

          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi(k);
          interval_t xi = interval_t::bottom();
          switch(op) {
            case OP_AND:  xi = yi.And(zi); break;
            case OP_OR:   xi = yi.Or(zi) ; break;
            case OP_XOR:  xi = yi.Xor(zi); break;
            case OP_SHL:  xi = yi.Shl(zi); break;
            case OP_LSHR: xi = yi.LShr(zi); break;
            case OP_ASHR: xi = yi.AShr(zi); break;
          }
          set(x, xi);
        }
        
	////////
	//// boolean operations
	////////
	void assign_bool_cst(variable_t lhs, linear_constraint_t cst) override
	{
	  if (!m_bool_reasoning) return;
	  
          crab::CrabStats::count(getDomainName() + ".count.assign_bool_cst");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_cst");
	  
	  if (is_bottom()) return;
	  
	  if (cst.is_tautology()) {
	    // this->operator-=(lhs);
	    // m_ldd = lddPtr(get_ldd_man(),
	    // 		    Ldd_And(get_ldd_man(), &*m_ldd, mk_true(lhs)));
	    
	    assign(lhs, number_t(1));
	  } else if (cst.is_contradiction()) {
	    // this->operator-=(lhs);	    
	    // m_ldd = lddPtr(get_ldd_man(),
	    // 		    Ldd_And(get_ldd_man(), &*m_ldd, mk_false(lhs)));
	    
	    assign(lhs, number_t(0));
	  } else {
	    linear_expression_t exp = cst.expression();    
	    unsigned int size = exp.size();
	    if (size == 0) return; // this should not happen 
	    else if (size == 1) {
	      auto it = exp.begin();
	      number_t cx = it->first;
	      if(cx == 1 || cx == -1) {
		// XXX: lhs should not appear in cst so we can remove lhs
		// without losing precision
		this->operator-=(lhs);
		
		number_t k = -exp.constant(); 
		variable_t vx = it->second;
		// m_ldd &= ite (cst, lhs >= 1, lhs <= 0);
		LddNodePtr ldd_cst = make_unit_constraint(cx, vx, cst.kind(),k);
		LddNodePtr ldd_rhs = lddPtr(get_ldd_man(),
					    Ldd_Ite(get_ldd_man(),
						    &*(ldd_cst),
						    mk_true(lhs), mk_false(lhs)));
		m_ldd = lddPtr(get_ldd_man(),
				Ldd_And(get_ldd_man(), &*m_ldd, &*ldd_rhs));
	      } else {
		CRAB_LOG("boxes",
		crab::outs() << "non-unit coefficients not implemented\n";);
	      }
	    } else if (size >= 2) {
	      CRAB_LOG("boxes",
	      crab::outs() << "non-unary constraints for boolean ops not implemented\n";);
	    }
	  }
	  CRAB_LOG("boxes",
		   crab::outs() << lhs << ":= " << "(" << cst << ")\n" << *this << "\n");
	}    
	
	void assign_bool_var(variable_t x, variable_t y, bool is_not_y) override {
	  if (!m_bool_reasoning) return;
	  
          crab::CrabStats::count(getDomainName() + ".count.assign_bool_var");
          crab::ScopedCrabStats __st__(getDomainName() + ".assign_bool_var");

	  if (is_not_y)
	    apply_ldd(x, y, number_t(-1), number_t(1));
	  else
	    copy_term(x, y);
	  
	  CRAB_LOG("boxes",
		   crab::outs()  << x << ":=";
		   if (is_not_y)
		     crab::outs() << "not(" << y << ")";
		   else
		     crab::outs() << y;		     
		   crab::outs() << "\n" << *this << "\n");
	}
	
	void apply_binary_bool(bool_operation_t op, variable_t x,
			       variable_t y, variable_t z) override {
	  if (!m_bool_reasoning) return;
	  
          crab::CrabStats::count(getDomainName() + ".count.apply_bin_bool");
          crab::ScopedCrabStats __st__(getDomainName() + ".apply_bin_bool");
	  
	  // XXX: if *lhs is null then it represents the SPECIAL
	  // variable $0.
	  boost::optional<variable_t> lhs; 
	  
	  if (!(x == y) && !(x == z)) {
	    lhs = boost::optional<variable_t>(x);
	    // XXX: x does not appear on the rhs so we can remove it
	    // without losing precision.
	    this->operator-=(x);
	  }
	  
	  m_ldd = lddPtr(get_ldd_man(),
			  Ldd_And(get_ldd_man(),
				  &*m_ldd,
				  &*gen_binary_bool(op, lhs, y, z)));
	  
	  if ((x == y) || (x == z)) {
	    // XXX: if we are here we added ite(y op z, $0 >= 1, $0 <= 0);
	    // so we still need to assign $0 to x:
	    copy_term(x, boost::optional<variable_t>());
	  }

	  CRAB_LOG("boxes",
		   crab::outs() << x << ":= " << y << " " << op << " " << z << "\n"
     		                 << *this << "\n");
	  
	}

	void assume_bool(variable_t x, bool is_negated) override {
	  if (!m_bool_reasoning) return;
	  
          crab::CrabStats::count(getDomainName() + ".count.assume_bool");
          crab::ScopedCrabStats __st__(getDomainName() + ".assume_bool");
	  
	  m_ldd = lddPtr(get_ldd_man(),
			  Ldd_And(get_ldd_man(), &*m_ldd,
				  ((is_negated) ?
				    mk_false(x): mk_true(x))));

	  CRAB_LOG("boxes",
	     if (!is_negated)  {
	       crab::outs() << "--- bool_assume(" << x << ")" << "\n" << *this << "\n";
	     } else { 
	       crab::outs() << "--- bool_assume(not(" << x << "))" << "\n" << *this << "\n";
	     });
	}
	
	// Backward boolean operations
	void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
				      boxes_domain_t inv){
	  if (is_bottom()) return;

	  this->operator-=(lhs);
	  CRAB_WARN("boxes backward boolean assignment not implemented");
	}
	
	void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
				      boxes_domain_t inv) {
	  if (is_bottom()) return;
	  
	  this->operator-=(lhs);
	  CRAB_WARN("boxes backward boolean assignment not implemented");
	}
	
	void backward_apply_binary_bool(bool_operation_t op,
					variable_t x, variable_t y, variable_t z,
					boxes_domain_t inv) {
	  if (is_bottom()) return;

	  this->operator-=(x);
	  CRAB_WARN("boxes backward boolean apply not implemented");
	}

	/*
	  Begin unimplemented operations.
	  
	  boxes_domain implements standard abstract operations of a
	  numerical domain plus boolean operations.  The
	  implementation of array and pointer operations is empty
	  because they should never be called.
	*/
	// array operations
	void array_init(variable_t a, linear_expression_t elem_size,
			linear_expression_t lb_idx, linear_expression_t ub_idx, 
			linear_expression_t val) {}      
	void array_load(variable_t lhs,
			variable_t a, linear_expression_t elem_size,
			linear_expression_t i) {}
	void array_store(variable_t a, linear_expression_t elem_size,
			 linear_expression_t i, linear_expression_t v, 
			 bool is_strong_update) {}
	void array_store(variable_t a_new, variable_t a_old,
			 linear_expression_t elem_size,
			 linear_expression_t i, linear_expression_t v, 
			 bool is_strong_update) {}      	
	void array_store_range(variable_t a, linear_expression_t elem_size,
			       linear_expression_t i, linear_expression_t j,
			       linear_expression_t v) {}
	void array_store_range(variable_t a_new, variable_t a_old,
			       linear_expression_t elem_size,
			       linear_expression_t i, linear_expression_t j,
			       linear_expression_t v) {}	
	void array_assign(variable_t lhs, variable_t rhs) {}
	// backward array operations
	void backward_array_init(variable_t a, linear_expression_t elem_size,
				 linear_expression_t lb_idx, linear_expression_t ub_idx, 
				 linear_expression_t val, boxes_domain_t invariant) {}
	void backward_array_load(variable_t lhs,
				 variable_t a, linear_expression_t elem_size,
				 linear_expression_t i, boxes_domain_t invariant) {}
	void backward_array_store(variable_t a, linear_expression_t elem_size,
				  linear_expression_t i, linear_expression_t v, 
				  bool is_strong_update, boxes_domain_t invariant) {}
	void backward_array_store(variable_t a_new, variable_t a_old,
				  linear_expression_t elem_size,
				  linear_expression_t i, linear_expression_t v, 
				  bool is_strong_update, boxes_domain_t invariant) {}	
	void backward_array_store_range(variable_t a, linear_expression_t elem_size,
					linear_expression_t i, linear_expression_t j,
					linear_expression_t v, boxes_domain_t invariant) {}
	void backward_array_store_range(variable_t a_new, variable_t a_old,
					linear_expression_t elem_size,
					linear_expression_t i, linear_expression_t j,
					linear_expression_t v, boxes_domain_t invariant) {}	
	void backward_array_assign(variable_t lhs, variable_t rhs, boxes_domain_t invariant) {}
	// pointer operations
	void pointer_load(variable_t lhs, variable_t rhs)  {}
	void pointer_store(variable_t lhs, variable_t rhs) {} 
	void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
	void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
	void pointer_function(variable_t lhs, varname_t func) {}
	void pointer_mk_null(variable_t lhs) {}
	void pointer_assume(pointer_constraint_t cst) {}
	void pointer_assert(pointer_constraint_t cst) {}
	/* End unimplemented operations */
	
        linear_constraint_system_t to_linear_constraint_system() {
	  crab::CrabStats::count(getDomainName() + ".count.to_linear_constraint_system");
	  crab::ScopedCrabStats __st__(getDomainName() + ".to_linear_constraint_system");
	  
          linear_constraint_system_t csts;
    
          if(is_bottom()) {
            csts += linear_constraint_t::get_false();
          }
          else if(is_top()) {
            csts += linear_constraint_t::get_true();
          } else {
	    // --- produce convex approximation
	    LddNodePtr ldd(m_ldd);	    
	    ldd = convex_approx(ldd);
	    // --- extract linear inequalities from the convex ldd
	    auto disjs = to_disjunctive_linear_constraint_system(ldd);
	    if (disjs.is_false()) {
	      csts += linear_constraint_t::get_false();
	    } else if (disjs.is_true()) {
	      csts += linear_constraint_t::get_true();
	    } else {
	      if (disjs.size() != 1) {
		CRAB_ERROR("Boxes::to_linear_constraint_system: it should not be disjunctive ",
			   disjs);
	      }
	      for (auto c : *(disjs.begin())) {
		csts += c;
	      }
	    }
	  }
	  
          return csts;
        }

	disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system(LddNodePtr &ldd) {
	  LddManager *ldd_man = get_ldd_man();
	  
	  std::vector<int> list;
	  list.reserve (ldd_man->cudd->size);
	  for (int i=0; i < ldd_man->cudd->size; i++) list.push_back(2);
	  
	  disjunctive_linear_constraint_system_t r;
	  to_disjunctive_linear_constraint_system_aux(ldd_man, ldd.get(),
						       r, list);
	  return r;
	}

	disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system() {
	  return to_disjunctive_linear_constraint_system(m_ldd);
	}
	
	void write(crab_os& o) {
	  crab::CrabStats::count(getDomainName() + ".count.write");
	  crab::ScopedCrabStats __st__(getDomainName() + ".write");
	  
          if (is_top()) {
            o << "{}";
	  } else if(is_bottom())  {
            o << "_|_";	
	  } else  {
	    auto r = to_disjunctive_linear_constraint_system();
	    o << r;
          }
        }

        static std::string getDomainName()
	{ return "Boxes"; }        
      }; 

     template<typename N, typename V, int R, size_t S>
     LddManager* boxes_domain_<N,V,R,S>::m_ldd_man = nullptr;

     template<typename N, typename V, int R, size_t S>
     typename boxes_domain_<N,V,R,S>::var_map_t boxes_domain_<N,V,R,S>::m_var_map;

     #if 1
     // Without copy-on-write
     template<typename Number, typename VariableName, int ConvexReduce=-1, size_t LddSize=3000>
     using boxes_domain = boxes_domain_<Number,VariableName,ConvexReduce,LddSize>;     
     #else

     template <typename Number, typename VariableName, int ConvexReduce, size_t LddSize> 
     struct abstract_domain_traits<boxes_domain_<Number, VariableName, ConvexReduce, LddSize>> {
       typedef Number number_t;
       typedef VariableName varname_t;       
     };
     
     // Quick wrapper which uses shared references with copy-on-write.
     template<class Number, class VariableName, int ConvexReduce=-1, size_t LddSize=3000>
     class boxes_domain final: 
       public abstract_domain<boxes_domain<Number,VariableName,ConvexReduce,LddSize>> {
			      
       typedef boxes_domain<Number, VariableName, ConvexReduce, LddSize> boxes_domain_t;
       typedef abstract_domain<boxes_domain_t> abstract_domain_t;
      
      public:
      
       using typename abstract_domain_t::linear_expression_t;
       using typename abstract_domain_t::linear_constraint_t;
       using typename abstract_domain_t::linear_constraint_system_t;
       using typename abstract_domain_t::disjunctive_linear_constraint_system_t;        
       using typename abstract_domain_t::variable_t;
       using typename abstract_domain_t::variable_vector_t;
       using typename abstract_domain_t::pointer_constraint_t;
       typedef Number number_t;
       typedef VariableName varname_t;
       typedef typename linear_constraint_t::kind_t constraint_kind_t;
       typedef interval<number_t>  interval_t;


     private:
       
       typedef boxes_domain_<number_t, varname_t, ConvexReduce, LddSize> boxes_impl_t;
       typedef std::shared_ptr<boxes_impl_t> boxes_ref_t;
       
       boxes_ref_t _ref;
       
       boxes_domain(boxes_ref_t ref) : _ref(ref) { }
       
       boxes_domain_t create(boxes_impl_t&& t) {
	 return std::make_shared<boxes_impl_t>(std::move(t));
       }
       
       void detach(void) {
	 if(!_ref || !_ref.unique()) {
	   _ref = std::make_shared<boxes_impl_t>(*_ref);
	 }
       }
       
       boxes_impl_t& ref(void) { return *_ref; }
       
       const boxes_impl_t& ref(void) const { return *_ref; }
       
     public:
       
       boxes_domain(bool is_bottom = false): _ref(nullptr) {
	 if (is_bottom) {
	   _ref = std::make_shared<boxes_impl_t>(boxes_impl_t::bottom());
	 } else {
	  _ref = std::make_shared<boxes_impl_t>(boxes_impl_t::top());
	 }
       }
       
       void set_to_top() {
	 boxes_domain abs(false);
	 std::swap(*this, abs);
       }
       
       void set_to_bottom() {
	 boxes_domain abs(true);
	 std::swap(*this, abs);
       }
       
       boxes_domain(const boxes_domain_t& o)
	 : _ref(o._ref) { }
      
       boxes_domain_t& operator=(const boxes_domain_t& o) {
	 _ref = o._ref;
	 return *this;
       }
       
       bool is_bottom() { return ref().is_bottom(); }
       bool is_top() { return ref().is_top(); }
       bool operator<=(boxes_domain_t o) { return ref() <= o.ref(); }
       void operator|=(boxes_domain_t o) { detach(); ref() |= o.ref(); }
       boxes_domain_t operator|(boxes_domain_t o) { return create(ref() | o.ref()); }
       boxes_domain_t operator||(boxes_domain_t o) { return create(ref() || o.ref()); }
       boxes_domain_t operator&(boxes_domain_t o) { return create(ref() & o.ref()); }
       boxes_domain_t operator&&(boxes_domain_t o) { return create(ref() && o.ref()); }
       
       boxes_domain_t widening_thresholds(boxes_domain_t o,
					  const iterators::thresholds<number_t> &ts) {
	 return create(ref().widening_thresholds(o.ref(), ts));
       }
       
       void operator+=(linear_constraint_system_t csts) { detach(); ref() += csts; } 
       void operator-=(variable_t v) { detach(); ref() -= v; }
       interval_t operator[](variable_t x) { return ref()[x]; }
       void set(variable_t x, interval_t intv) { detach(); ref().set(x, intv); }
       
       void assign(variable_t x, linear_expression_t e)
       { detach(); ref().assign(x, e); }
       void apply(operation_t op, variable_t x, variable_t y, number_t k)
       { detach(); ref().apply(op, x, y, k); }
       void apply(operation_t op, variable_t x, variable_t y, variable_t z)
       { detach(); ref().apply(op, x, y, z); }
       
       void backward_assign(variable_t x, linear_expression_t e,
			    boxes_domain_t invariant)
       { detach(); ref().backward_assign(x, e, invariant.ref()); }
       void backward_apply(operation_t op,
			   variable_t x, variable_t y, number_t k,
			   boxes_domain_t invariant)
       { detach(); ref().backward_apply(op, x, y, k, invariant.ref()); }
       void backward_apply(operation_t op,
			   variable_t x, variable_t y, variable_t z,
			   boxes_domain_t invariant)
       { detach(); ref().backward_apply(op, x, y, z, invariant.ref()); }
       
       void apply(int_conv_operation_t op, variable_t dst, variable_t src)
       { detach(); ref().apply(op, dst, src); }
       void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k)
       { detach(); ref().apply(op, x, y, k); }
       void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z)
       { detach(); ref().apply(op, x, y, z); }
       
       void assign_bool_cst(variable_t x, linear_constraint_t cst)
       { detach(); ref().assign_bool_cst(x,cst);}
       void assign_bool_var(variable_t x, variable_t y, bool is_not_y)
       { detach(); ref().assign_bool_var(x,y,is_not_y);}
       void apply_binary_bool(bool_operation_t op,variable_t x, variable_t y, variable_t z)
       { detach(); ref().apply_binary_bool(op,x,y,z);}
       void assume_bool(variable_t x, bool is_negated)
       { detach(); ref().assume_bool(x,is_negated);}
       
       void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
				    boxes_domain_t inv)
       { detach(); ref().backward_assign_bool_cst(lhs,rhs,inv.ref()); }
       void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
				     boxes_domain_t inv)
       { detach(); ref().backward_assign_bool_var(lhs,rhs,is_not_rhs,inv.ref()); }	
       void backward_apply_binary_bool(bool_operation_t op,
				       variable_t x,variable_t y,variable_t z,
				       boxes_domain_t inv)
       { detach(); ref().backward_apply_binary_bool(op,x,y,z,inv.ref()); }	
       
       /* Begin unimplemented operations */
       // array operations
       void array_init(variable_t a, linear_expression_t elem_size,
		       linear_expression_t lb_idx, linear_expression_t ub_idx, 
		       linear_expression_t val) {}      
       void array_load(variable_t lhs,
		       variable_t a, linear_expression_t elem_size,
		       linear_expression_t i) {}
       void array_store(variable_t a, linear_expression_t elem_size,
			linear_expression_t i, linear_expression_t v, 
			bool is_strong_update) {}
       void array_store(variable_t a_new, variable_t a_old,
			linear_expression_t elem_size,
			linear_expression_t i, linear_expression_t v, 
			bool is_strong_update) {}             
       void array_store_range(variable_t a, linear_expression_t elem_size,
			      linear_expression_t i, linear_expression_t j,
			      linear_expression_t v) {}
       void array_store_range(variable_t a_new, variable_t a_old,
			      linear_expression_t elem_size,
			      linear_expression_t i, linear_expression_t j,
			      linear_expression_t v) {}       
       void array_assign(variable_t lhs, variable_t rhs) {}
       // backward array operations
       void backward_array_init(variable_t a, linear_expression_t elem_size,
				linear_expression_t lb_idx, linear_expression_t ub_idx, 
				linear_expression_t val, boxes_domain_t invariant) {}
       void backward_array_load(variable_t lhs,
				variable_t a, linear_expression_t elem_size,
				linear_expression_t i, boxes_domain_t invariant) {}
       void backward_array_store(variable_t a, linear_expression_t elem_size,
				 linear_expression_t i, linear_expression_t v, 
				 bool is_strong_update, boxes_domain_t invariant) {}
       void backward_array_store(variable_t a_new, variable_t a_old,
				 linear_expression_t elem_size,
				 linear_expression_t i, linear_expression_t v, 
				 bool is_strong_update, boxes_domain_t invariant) {}       
       void backward_array_store_range(variable_t a, linear_expression_t elem_size,
				       linear_expression_t i, linear_expression_t j,
				       linear_expression_t v, boxes_domain_t invariant) {}
       void backward_array_store_range(variable_t a_new, variable_t a_old,
				       linear_expression_t elem_size,
				       linear_expression_t i, linear_expression_t j,
				       linear_expression_t v, boxes_domain_t invariant) {}       
       void backward_array_assign(variable_t lhs, variable_t rhs, boxes_domain_t invariant) {}
       // pointer operations
       void pointer_load(variable_t lhs, variable_t rhs)  {}
       void pointer_store(variable_t lhs, variable_t rhs) {} 
       void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
       void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
       void pointer_function(variable_t lhs, varname_t func) {}
       void pointer_mk_null(variable_t lhs) {}
       void pointer_assume(pointer_constraint_t cst) {}
       void pointer_assert(pointer_constraint_t cst) {}
       /* End unimplemented operations */
       
       void forget(const variable_vector_t& variables)
       { detach(); ref().forget(variables); }
       
       void project(const variable_vector_t& variables)
       { detach(); ref().project(variables); }
       
       void expand(variable_t x, variable_t new_x)
       { detach(); ref().expand(x,new_x);}

       void normalize() {}

       void minimize() {}       
       
       void write(crab_os& o) { ref().write(o); }
       
       linear_constraint_system_t to_linear_constraint_system()
       { return ref().to_linear_constraint_system(); }

       disjunctive_linear_constraint_system_t
       to_disjunctive_linear_constraint_system()
       { return ref().to_disjunctive_linear_constraint_system(); }
       
       static std::string getDomainName() { return boxes_impl_t::getDomainName(); }
       
    };     
    #endif
     
    template<typename Number, typename VariableName, int ConvexReduce, size_t LddSize>
    class checker_domain_traits<boxes_domain<Number,VariableName,ConvexReduce,LddSize>> {
    public:
      typedef boxes_domain<Number, VariableName, ConvexReduce, LddSize> this_type;
      typedef typename this_type::linear_constraint_t linear_constraint_t;
      typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;    
      
      static bool entail(this_type& lhs, const disjunctive_linear_constraint_system_t& rhs) {
	// -- trivial cases first
	if (rhs.is_false()) {
	  return false;
	} else if(rhs.is_true()) {
	  return true;
	} else if (lhs.is_bottom()) {
	  return true;
	} else if (lhs.is_top()) {
	  return false;
	}
	this_type inv = this_type::bottom();
	for (auto const& csts: rhs) {
	  this_type conj;
	  conj += csts;
	  inv  |= conj;
	}
	return (lhs & inv.complement()).is_bottom();
      }
      
      static bool entail(const disjunctive_linear_constraint_system_t& lhs, this_type& rhs) {
	// -- trivial cases first
	if (rhs.is_bottom()) {
	  return false;
	} else if (rhs.is_top()) {
	  return true;
	} else if (lhs.is_false()) {
	  return true;
	} else if (lhs.is_true()) {
	  return false;
	}
	this_type inv = this_type::bottom();
	for (auto const& csts: lhs) {
	  this_type conj;
	  conj += csts;
	  inv  |= conj;
	}
	return (inv & rhs.complement()).is_bottom();
      }
      
      static bool entail(this_type& lhs, const linear_constraint_t& rhs) {
	// -- trivial cases first
	if (lhs.is_bottom()) return true;
	if (rhs.is_tautology()) return true;
	if (rhs.is_contradiction()) return false;

	this_type inv(lhs);
	inv += rhs.negate();
	return inv.is_bottom();
      }
      
      static bool intersect(this_type& inv, const linear_constraint_t& cst) {
	// default code

	// -- trivial cases first
	if (inv.is_bottom() || cst.is_contradiction()) return false;
	if (inv.is_top() || cst.is_tautology()) return true;
	
	this_type cst_inv;
	cst_inv += cst;
	return (!(cst_inv & inv).is_bottom());
      }
      
    };
     
   } // namespace domains
}// namespace crab
#endif /* HAVE_LDD */

namespace crab {
namespace domains {
  template <typename Number, typename VariableName, int ConvexReduce, size_t LddSize> 
  struct abstract_domain_traits<boxes_domain<Number, VariableName, ConvexReduce, LddSize>> {
    typedef Number number_t;
    typedef VariableName varname_t;       
  };
}
}
  
