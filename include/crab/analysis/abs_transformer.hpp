#ifndef ABSTRACT_TRANSFORMER_HPP
#define ABSTRACT_TRANSFORMER_HPP

/* 
   Implementation of the abstract transfer functions by reducing them
   to abstract domain operations.
   
   These are the main Crab statements for which we define their abstract
   transfer functions:
   
   ARITHMETIC and BOOLEAN
     x := y bin_op z;
     x := y; 
     assume (cst) 
     assert (cst);
     x := select (cond, y, z);

   ARRAYS
     a[i] = v;    (v can be either an integer or pointer)
     v = a[i];
     a := b       (a and b arrays)

   POINTERS
     *p = q;
     p = *q;
     p := q+n 
     p := &obj;
     p := &fun
     p := null;

   FUNCTIONS
     x := foo(arg1,...,argn);
     return r;

   havoc (x);

 */

#include <boost/optional.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/noncopyable.hpp>

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/domains/linear_constraints.hpp>

namespace crab {

  namespace analyzer {

  /**
   * API abstract transformer
   **/
  template<typename Number, typename VariableName>
  class abs_transformer_api: public cfg::statement_visitor<Number, VariableName>
  {
   public:

    typedef Number number_t;    
    typedef VariableName varname_t;
    
    typedef ikos::variable <number_t, VariableName > var_t;
    typedef ikos::linear_expression<number_t, VariableName > lin_exp_t;
    typedef ikos::linear_constraint<number_t, VariableName > lin_cst_t;
    typedef ikos::linear_constraint_system<number_t, VariableName > lin_cst_sys_t;

    typedef cfg::havoc_stmt<number_t,VariableName>        havoc_t;
    typedef cfg::unreachable_stmt<number_t,VariableName>  unreach_t;
    
    typedef cfg::binary_op<number_t,VariableName>         bin_op_t;
    typedef cfg::assignment<number_t,VariableName>        assign_t;
    typedef cfg::assume_stmt<number_t,VariableName>       assume_t;
    typedef cfg::select_stmt<number_t,VariableName>       select_t;
    typedef cfg::assert_stmt<number_t,VariableName>       assert_t;

    typedef cfg::int_cast_stmt<number_t,VariableName>     int_cast_t;
    typedef cfg::callsite_stmt<number_t,VariableName>     callsite_t;
    typedef cfg::return_stmt<number_t,VariableName>       return_t;
    
    typedef cfg::array_assume_stmt<number_t,VariableName> arr_assume_t;
    typedef cfg::array_store_stmt<number_t,VariableName>  arr_store_t;
    typedef cfg::array_load_stmt<number_t,VariableName>   arr_load_t;
    typedef cfg::array_assign_stmt<number_t,VariableName> arr_assign_t;
    typedef cfg::ptr_store_stmt<number_t,VariableName>    ptr_store_t;
    typedef cfg::ptr_load_stmt<number_t,VariableName>     ptr_load_t;

    typedef cfg::ptr_assign_stmt<number_t,VariableName>   ptr_assign_t;
    typedef cfg::ptr_object_stmt<number_t,VariableName>   ptr_object_t;
    typedef cfg::ptr_function_stmt<number_t,VariableName> ptr_function_t;
    typedef cfg::ptr_null_stmt<number_t,VariableName>     ptr_null_t;
    typedef cfg::ptr_assume_stmt<number_t,VariableName>   ptr_assume_t;
    typedef cfg::ptr_assert_stmt<number_t,VariableName>   ptr_assert_t;

    typedef cfg::bool_binary_op<number_t,VariableName>    bool_bin_op_t;
    typedef cfg::bool_assign_cst<number_t,VariableName>   bool_assign_cst_t;
    typedef cfg::bool_assign_var<number_t,VariableName>   bool_assign_var_t;    
    typedef cfg::bool_assume_stmt<number_t,VariableName>  bool_assume_t;
    typedef cfg::bool_select_stmt<number_t,VariableName>  bool_select_t;
    typedef cfg::bool_assert_stmt<number_t,VariableName>  bool_assert_t;
    
   protected: 

    virtual void exec (havoc_t&) { }
    virtual void exec (unreach_t&) { }    
    virtual void exec (bin_op_t&)  { } 
    virtual void exec (assign_t&) { }
    virtual void exec (assume_t&) { }
    virtual void exec (select_t&) { }
    virtual void exec (assert_t&) { }
    virtual void exec (int_cast_t&) { }
    virtual void exec (callsite_t&) { }
    virtual void exec (return_t&) { }
    virtual void exec (arr_assume_t&) { }
    virtual void exec (arr_store_t&) { }
    virtual void exec (arr_load_t&) { }
    virtual void exec (arr_assign_t&) { }
    virtual void exec (ptr_store_t&) { }
    virtual void exec (ptr_load_t&) { }
    virtual void exec (ptr_assign_t&) { }
    virtual void exec (ptr_object_t&) { }
    virtual void exec (ptr_function_t&) { }
    virtual void exec (ptr_null_t&) { }
    virtual void exec (ptr_assume_t&) { }
    virtual void exec (ptr_assert_t&) { }
    virtual void exec (bool_bin_op_t&)  { } 
    virtual void exec (bool_assign_cst_t&) { }
    virtual void exec (bool_assign_var_t&) { }    
    virtual void exec (bool_assume_t&) { }
    virtual void exec (bool_select_t&) { }
    virtual void exec (bool_assert_t&) { }

   public: /* visitor api */

    void visit (havoc_t &s) { exec (s); }
    void visit (unreach_t &s) { exec (s); }    
    void visit (bin_op_t &s) { exec (s); }
    void visit (assign_t &s) { exec (s); }
    void visit (assume_t &s) { exec (s); }
    void visit (select_t &s) { exec (s); }
    void visit (assert_t &s) { exec (s); }
    void visit (int_cast_t &s) { exec (s); }    
    void visit (callsite_t &s) { exec (s); }
    void visit (return_t &s) { exec (s); }
    void visit (arr_assume_t &s) { exec (s); }
    void visit (arr_store_t &s) { exec (s); }
    void visit (arr_load_t &s) { exec (s); }
    void visit (arr_assign_t &s) { exec (s); }
    void visit (ptr_store_t &s) { exec (s); }
    void visit (ptr_load_t &s) { exec (s); }
    void visit (ptr_assign_t &s) { exec (s); }
    void visit (ptr_object_t &s) { exec (s); }
    void visit (ptr_function_t &s) { exec (s); }
    void visit (ptr_null_t &s) { exec (s); }
    void visit (ptr_assume_t &s) { exec (s); }
    void visit (ptr_assert_t &s) { exec (s); }
    void visit (bool_bin_op_t &s) { exec (s); }
    void visit (bool_assign_cst_t &s) { exec (s); }
    void visit (bool_assign_var_t &s) { exec (s); }    
    void visit (bool_assume_t &s) { exec (s); }
    void visit (bool_select_t &s) { exec (s); }
    void visit (bool_assert_t &s) { exec (s); }
  };


  /** 
   * Abstract forward transformer for all statements. Function calls
   * which are ingnored in a sound manner.
   **/
  template<class AbsD>
  class intra_abs_transformer: 
      public abs_transformer_api<typename AbsD::number_t, typename AbsD::varname_t>,
      public boost::noncopyable {
    
   public:
    typedef AbsD abs_dom_t;
    typedef typename abs_dom_t::number_t number_t;
    typedef typename abs_dom_t::varname_t varname_t;

   public:
    typedef abs_transformer_api <number_t, varname_t> abs_transform_api_t;
    using typename abs_transform_api_t::var_t;
    using typename abs_transform_api_t::lin_exp_t;
    using typename abs_transform_api_t::lin_cst_t;
    using typename abs_transform_api_t::lin_cst_sys_t;
    using typename abs_transform_api_t::havoc_t;
    using typename abs_transform_api_t::unreach_t;    
    using typename abs_transform_api_t::bin_op_t;
    using typename abs_transform_api_t::assign_t;
    using typename abs_transform_api_t::assume_t;
    using typename abs_transform_api_t::select_t;
    using typename abs_transform_api_t::assert_t;
    using typename abs_transform_api_t::int_cast_t;    
    using typename abs_transform_api_t::callsite_t;
    using typename abs_transform_api_t::arr_assume_t;
    using typename abs_transform_api_t::arr_load_t;
    using typename abs_transform_api_t::arr_store_t;
    using typename abs_transform_api_t::arr_assign_t;
    using typename abs_transform_api_t::ptr_load_t;
    using typename abs_transform_api_t::ptr_store_t;
    using typename abs_transform_api_t::ptr_assign_t;
    using typename abs_transform_api_t::ptr_object_t;
    using typename abs_transform_api_t::ptr_function_t;
    using typename abs_transform_api_t::ptr_null_t;
    using typename abs_transform_api_t::ptr_assume_t;
    using typename abs_transform_api_t::ptr_assert_t;
    using typename abs_transform_api_t::bool_bin_op_t;
    using typename abs_transform_api_t::bool_assign_cst_t;
    using typename abs_transform_api_t::bool_assign_var_t;    
    using typename abs_transform_api_t::bool_assume_t;
    using typename abs_transform_api_t::bool_select_t;
    using typename abs_transform_api_t::bool_assert_t;

    
   protected:
    /// XXX: the transformer does not own m_inv.
    abs_dom_t *m_inv;

   private:

    template <typename NumOrVar>
    void apply (abs_dom_t &inv, binary_operation_t op, 
		varname_t x, varname_t y, NumOrVar z)
    {
      auto op1 = conv_op<ikos::operation_t> (op);
      auto op2 = conv_op<ikos::div_operation_t> (op);
      auto op3 = conv_op<ikos::bitwise_operation_t> (op);
      if (op1)      inv.apply (*op1, x, y, z); 
      else if (op2) inv.apply (*op2, x, y, z);
      else if (op3) inv.apply (*op3, x, y, z);
      else CRAB_ERROR("unsupported binary operator", op);
    }

    abs_dom_t *get_inv () {
      if (!m_inv)
	CRAB_ERROR ("Invariant passed to transformer cannot be null!");
      return m_inv;
    }
    
   public:
    
    intra_abs_transformer (abs_dom_t* inv):
      m_inv (inv) { }
    
    virtual ~intra_abs_transformer () { }

    void set (abs_dom_t& inv) { m_inv = &inv;}
    
    abs_dom_t& inv () { return *get_inv(); }
    
    void exec (bin_op_t& stmt) {
      auto op1 = stmt.left ();
      auto op2 = stmt.right ();
      if (op1.get_variable () && op2.get_variable ()) {
        apply (*get_inv(), stmt.op (), 
               stmt.lhs ().name(), 
               (*op1.get_variable ()).name(), 
               (*op2.get_variable ()).name());
      } else {
        assert (op1.get_variable () && op2.is_constant ());
        apply (*get_inv(), stmt.op (), 
               stmt.lhs ().name (), 
               (*op1.get_variable ()).name(), 
               op2.constant ()); 
      }      
    }
    
    void exec (select_t& stmt) {
      abs_dom_t inv1 (*get_inv());
      abs_dom_t inv2 (*get_inv());
      inv1 += stmt.cond ();
      inv2 += stmt.cond ().negate ();
      if (inv2.is_bottom()) {
        inv1.assign (stmt.lhs().name (),stmt.left());
        *get_inv() = inv1;
      }
      else if (inv1.is_bottom ()) {
        inv2.assign (stmt.lhs().name (),stmt.right());
        *get_inv() = inv2;
      }
      else {
        inv1.assign (stmt.lhs().name (),stmt.left());
        inv2.assign (stmt.lhs().name (),stmt.right());
        *get_inv() = inv1 | inv2;
      }
    }
    
    void exec (assign_t& stmt) {
      get_inv()->assign (stmt.lhs().name (), lin_exp_t (stmt.rhs()));
    }
    
    void exec (assume_t& stmt) {
      *get_inv() += stmt.constraint();
    }
    
    void exec (assert_t& stmt) {
      abs_dom_t cst;
      cst += stmt.constraint();
      abs_dom_t meet = cst & *get_inv();
      if (meet.is_bottom ()) {
	// assertion does not definitely hold.
        *get_inv() = abs_dom_t::bottom (); 
      } else {
	*get_inv() += stmt.constraint ();
      }
    }

    void exec (int_cast_t &stmt){
      if (auto op = conv_op<crab::domains::int_conv_operation_t>(stmt.op())) {
	get_inv()->apply(*op,
			 stmt.dst(), stmt.dst_width(),
			 stmt.src(), stmt.src_width());
      } else {
	CRAB_ERROR("unsupported cast operator ", stmt.op());
      }
    }
    
    void exec (bool_assign_cst_t& stmt) {
      (*get_inv()).assign_bool_cst (stmt.lhs (), stmt.rhs ());
    }

    void exec (bool_assign_var_t& stmt) {
      (*get_inv()).assign_bool_var (stmt.lhs (), stmt.rhs (), stmt.is_rhs_negated());      
    }

    void exec (bool_bin_op_t& stmt) {
      if (auto op = conv_op<domains::bool_operation_t> (stmt.op()))
	(*get_inv()).apply_binary_bool (*op, stmt.lhs(),
					stmt.left(), stmt.right ());
    }
    
    
    void exec (bool_assume_t& stmt) {
      (*get_inv()).assume_bool (stmt.cond(), stmt.is_negated ());
    }

    void exec (bool_select_t& stmt) {
      abs_dom_t inv1 (*get_inv());
      abs_dom_t inv2 (*get_inv());
      const bool negate = true;
      inv1.assume_bool(stmt.cond (), !negate);
      inv2.assume_bool(stmt.cond (), negate);
      if (inv2.is_bottom()) {
        inv1.assign_bool_var(stmt.lhs(),stmt.left(), !negate);
        *get_inv() = inv1;
      }
      else if (inv1.is_bottom ()) {
        inv2.assign_bool_var(stmt.lhs(),stmt.right(), !negate);
        *get_inv() = inv2;
      }
      else {
        inv1.assign_bool_var(stmt.lhs(),stmt.left(), !negate);
        inv2.assign_bool_var(stmt.lhs(),stmt.right(), !negate);
        *get_inv() = inv1 | inv2;
      }
    }

    void exec (bool_assert_t& stmt) {
      abs_dom_t inv;
      inv.assume_bool (stmt.cond (), false);
      abs_dom_t meet = inv & *get_inv();
      if (meet.is_bottom ()) {
	// assertion does not definitely hold.	
        *get_inv() = abs_dom_t::bottom (); 
      } else {
	(*get_inv()).assume_bool (stmt.cond(), false);
      }
    }
    
    void exec (havoc_t& stmt)
    { (*get_inv()) -= stmt.variable(); }
    
    void exec (unreach_t& stmt)
    { *get_inv() = abs_dom_t::bottom (); }
    
    void exec (arr_assume_t &stmt) {
      get_inv()->array_assume (stmt.array (), stmt.array_type (), 
			       stmt.lb_index (), stmt.ub_index (), 
			       stmt.val ());
    }
    
    void exec (arr_store_t &stmt) {
      get_inv()->array_store (stmt.array(), stmt.array_type (), 
			      stmt.index (), stmt.value (),
			      stmt.elem_size(), stmt.is_singleton ());
    }
    
    void exec (arr_load_t  &stmt) {
      get_inv()->array_load (stmt.lhs (), 
			     stmt.array (), stmt.array_type (), 
			     stmt.index(), stmt.elem_size());
    }
    
    void exec (arr_assign_t  &stmt) {
      get_inv()->array_assign (stmt.lhs (), stmt.rhs (),
			       stmt.array_type ());
    }
    
    void exec (ptr_null_t & stmt)
    { get_inv()->pointer_mk_null (stmt.lhs ()); }
    
    void exec (ptr_object_t & stmt)
    { get_inv()->pointer_mk_obj (stmt.lhs (), stmt.rhs()); }
    
    void exec (ptr_assign_t & stmt)
    { get_inv()->pointer_assign (stmt.lhs (), stmt.rhs (), stmt.offset ()); }
    
    void exec (ptr_function_t & stmt)
    { get_inv()->pointer_function (stmt.lhs (), stmt.rhs ()); }
    
    void exec (ptr_load_t & stmt)
    { get_inv()->pointer_load (stmt.lhs (), stmt.rhs ()); }
    
    void exec (ptr_store_t & stmt)
    { get_inv()->pointer_store (stmt.lhs (), stmt.rhs ()); }
    
    void exec (ptr_assume_t& stmt)
    { get_inv()->pointer_assume (stmt.constraint ()); }
    
    void exec (ptr_assert_t& stmt)
    { get_inv()->pointer_assert (stmt.constraint ()); }
    
    virtual void exec (callsite_t &cs) {
      for (auto vt: cs.get_lhs())
	(*get_inv()) -= vt.first;  // havoc
    }
  }; 


  /**
   * Abstract transformer to compute necessary preconditions of error
   * states.
   **/ 
  template<class AbsD, class InvT>
  class intra_necessary_preconditions_abs_transformer: 
      public abs_transformer_api<typename AbsD::number_t,
				 typename AbsD::varname_t>,
      public boost::noncopyable {
    
   public:
    typedef AbsD abs_dom_t;
    typedef typename abs_dom_t::number_t number_t;
    typedef typename abs_dom_t::varname_t varname_t;
    typedef cfg::statement<number_t, varname_t> statement_t;    
    typedef abs_transformer_api <number_t, varname_t> abs_transform_api_t;
    using typename abs_transform_api_t::var_t;
    using typename abs_transform_api_t::lin_exp_t;
    using typename abs_transform_api_t::lin_cst_t;
    using typename abs_transform_api_t::lin_cst_sys_t;
    using typename abs_transform_api_t::havoc_t;
    using typename abs_transform_api_t::unreach_t;    
    using typename abs_transform_api_t::bin_op_t;
    using typename abs_transform_api_t::assign_t;
    using typename abs_transform_api_t::assume_t;
    using typename abs_transform_api_t::select_t;
    using typename abs_transform_api_t::assert_t;
    using typename abs_transform_api_t::int_cast_t;
    using typename abs_transform_api_t::callsite_t;
    using typename abs_transform_api_t::arr_assume_t;
    using typename abs_transform_api_t::arr_load_t;
    using typename abs_transform_api_t::arr_store_t;
    using typename abs_transform_api_t::arr_assign_t;
    using typename abs_transform_api_t::ptr_load_t;
    using typename abs_transform_api_t::ptr_store_t;
    using typename abs_transform_api_t::ptr_assign_t;
    using typename abs_transform_api_t::ptr_object_t;
    using typename abs_transform_api_t::ptr_function_t;
    using typename abs_transform_api_t::ptr_null_t;
    using typename abs_transform_api_t::ptr_assume_t;
    using typename abs_transform_api_t::ptr_assert_t;
    using typename abs_transform_api_t::bool_bin_op_t;
    using typename abs_transform_api_t::bool_assign_cst_t;
    using typename abs_transform_api_t::bool_assign_var_t;    
    using typename abs_transform_api_t::bool_assume_t;
    using typename abs_transform_api_t::bool_select_t;
    using typename abs_transform_api_t::bool_assert_t;

    
   protected:
    
    // used to compute the (necessary) preconditions
    abs_dom_t *m_pre;
    // used to refine the preconditions: map from statement_t to
    // abs_dom_t.
    InvT &m_invariants; 

   public:

    intra_necessary_preconditions_abs_transformer(abs_dom_t* post,
						  InvT &invars)
      : m_pre (post), m_invariants (invars) {
      if (!m_pre)
	CRAB_ERROR ("Postcondition cannot be null!");
    }
    
    virtual ~intra_necessary_preconditions_abs_transformer () { }

    abs_dom_t& preconditions () { return *m_pre; }
    
    void exec (bin_op_t& stmt) {

      auto op = conv_op<ikos::operation_t> (stmt.op ());
      if (!op) {
	CRAB_WARN ("backward operation ", stmt.op(), " not implemented");
	(*m_pre) -= stmt.lhs().name();
	return;
      }
      
      auto op1 = stmt.left ();
      auto op2 = stmt.right ();
      abs_dom_t invariant = m_invariants[&stmt];

      CRAB_LOG("backward",
	       crab::outs () << "** " << stmt.lhs () << " := "
	                     << op1 << " " << *op << " " << op2 << "\n"
 	                     << "FORWARD INV=" << invariant << "\n"
	                     << "POST=" << *m_pre << "\n");	       
      
      if (op1.get_variable () && op2.get_variable ()) {
        m_pre->backward_apply (*op, 
				stmt.lhs ().name(), 
				(*op1.get_variable ()).name(), 
				(*op2.get_variable ()).name(),
				invariant);
      } else {
        assert (op1.get_variable () && op2.is_constant ());
        m_pre->backward_apply (*op, 
				stmt.lhs ().name (), 
				(*op1.get_variable ()).name(), 
				op2.constant (),
				invariant); 
      }
      CRAB_LOG("backward",
	       crab::outs () << "PRE=" << *m_pre << "\n");
      
    }
    
    // select(x := cond ? e1: e2, post) can be reduced to
    //   pre: goto b_then;
    //   pre: goto b_else;
    //   b_then:
    //     assume(cond);
    //     x := e1;
    //     goto post;
    //   b_else:
    //     assume(not(cond));
    //     x := e2;
    //     goto post;
    //   post: ....
    void exec (select_t& stmt) {
      abs_dom_t old_pre = m_invariants[&stmt];

      // -- one of the two branches is false
      abs_dom_t then_inv (old_pre);                  
      then_inv += stmt.cond();      
      if (then_inv.is_bottom()) {
        m_pre->backward_assign (stmt.lhs().name(),stmt.right(),old_pre);	
	*m_pre += stmt.cond().negate();
	return;
      }
      
      abs_dom_t else_inv (old_pre);
      else_inv += stmt.cond ().negate();
      if (else_inv.is_bottom()) {
        m_pre->backward_assign (stmt.lhs().name(),stmt.left(),old_pre);	
	*m_pre += stmt.cond ();
	return;
      }

      // -- both branches can be possible so we join them
      abs_dom_t pre_then (*m_pre);
      pre_then.backward_assign(stmt.lhs().name(),stmt.left(),old_pre);      
      pre_then += stmt.cond();
      
      abs_dom_t pre_else (*m_pre);
      pre_else.backward_assign(stmt.lhs().name(),stmt.right(),old_pre);      
      pre_else += stmt.cond().negate();
      
      *m_pre = pre_then | pre_else;
    }

    // x := e
    void exec (assign_t& stmt) {
      abs_dom_t invariant = m_invariants[&stmt];
      
      CRAB_LOG("backward",
	       auto rhs = stmt.rhs();
	       crab::outs () << "** " << stmt.lhs () << " := " << rhs << "\n"
 	                     << "FORWARD INV=" << invariant << "\n"
	                     << "POST=" << *m_pre << "\n");	       
      
      m_pre->backward_assign (stmt.lhs().name (),
			      lin_exp_t (stmt.rhs()),
			      invariant);
      CRAB_LOG("backward",
	       crab::outs () << "PRE=" << *m_pre << "\n");
    }

    // assume(c)
    // the precondition must contain c so forward and backward are the same.
    void exec (assume_t& stmt) {
      CRAB_LOG("backward",
	       auto csts = stmt.constraint();
	       crab::outs () << "** assume(" << csts << ")\n"
	                     << "POST=" << *m_pre << "\n");	       
      
      *m_pre += stmt.constraint();

      CRAB_LOG("backward",
	       crab::outs () << "PRE=" << *m_pre << "\n");
      
    }

    // We are interested in computing preconditions of the error
    // states. Thus, we propagate backwards the negation of the
    // condition which represents the error states.
    void exec (assert_t& stmt) {
      CRAB_LOG("backward",
	       auto csts = stmt.constraint();
	       crab::outs () << "** assert(" << csts << ")\n"
	                     << "POST=" << *m_pre << "\n");	       

      *m_pre = abs_dom_t::top();
      *m_pre += stmt.constraint().negate ();
      
      CRAB_LOG("backward",
	       crab::outs () << "PRE=" << *m_pre << "\n");
    }
        
    // We are interested in computing preconditions of the error states.
    // Thus, an unreachable state is always safe so we propagate true.
    void exec (unreach_t& stmt) {
      *m_pre = abs_dom_t::top ();
    }

    // x := *
    // x can be anything before the assignment
    void exec (havoc_t& stmt) {
      *m_pre -= stmt.variable();
    }

    // ret := *
    virtual void exec (callsite_t &cs) {
      for (auto vt: cs.get_lhs())
	*m_pre -= vt.first;  
    }

    void exec (int_cast_t &stmt) {
      abs_dom_t invariant = m_invariants[&stmt];
      // FIXME: treats cast as an assignment ignoring bitdwith
      m_pre->backward_assign(stmt.dst(), lin_exp_t(stmt.src()), invariant);
    }
    
    void exec (bool_assign_cst_t &stmt)
    { *m_pre -= stmt.lhs(); }
    void exec (bool_assign_var_t &stmt)
    { *m_pre -= stmt.lhs(); }
    void exec (bool_bin_op_t &stmt)
    { *m_pre -= stmt.lhs(); }
    void exec (bool_select_t &stmt)
    { *m_pre -= stmt.lhs(); }
    
    void exec (bool_assume_t &stmt) {
      // same as forward
      m_pre->assume_bool (stmt.cond(), stmt.is_negated ());
    }
    
    void exec (bool_assert_t &stmt) {
      *m_pre = abs_dom_t::top();
      m_pre->assume_bool (stmt.cond(), true /*negated*/);
    }

    void exec (arr_load_t  &stmt)
    { *m_pre -= stmt.lhs(); }          
    // NOT IMPLEMENTED
    void exec (arr_assume_t &stmt) { }
    void exec (arr_store_t &stmt) { }
    void exec (arr_assign_t  &stmt) { }
    void exec (ptr_null_t &stmt) { }
    void exec (ptr_object_t &stmt) { }
    void exec (ptr_assign_t &stmt) { }
    void exec (ptr_function_t &stmt) { }
    void exec (ptr_load_t &stmt) { }
    void exec (ptr_store_t &stmt) { }
    void exec (ptr_assume_t &stmt) { }
    void exec (ptr_assert_t &stmt) { }
  }; 

    
  // Helper for function calls.
  template <typename AbsDom>
  static void unify (AbsDom &inv, variable_type ty,
		     typename AbsDom::varname_t lhs,
		     typename AbsDom::varname_t rhs)
  {
    typedef typename AbsDom::linear_expression_t linear_expression_t;
    typedef typename AbsDom::number_t number_t;
    switch (ty)
    {
      case BOOL_TYPE:
	inv.assign_bool_var(lhs, rhs, false);
	break;      
      case INT_TYPE:
      case REAL_TYPE:	
	inv.assign (lhs, linear_expression_t (rhs));
	break;
      case PTR_TYPE:
	inv.pointer_assign (lhs, rhs, number_t (0));
	break;
      case ARR_BOOL_TYPE:
      case ARR_INT_TYPE:
      case ARR_REAL_TYPE:
      case ARR_PTR_TYPE:
	inv.array_assign (lhs, rhs, ty);
	break;
      default:
	CRAB_ERROR ("unsuported type");
    }
  }
    
  /**
   * Abstract transformer specialized for computing summaries.
   * class SumTable stores the summaries
   **/
  template<class SumTable> 
  class bu_summ_abs_transformer: 
      public intra_abs_transformer <typename SumTable::abs_domain_t> {
				      
   public:

    typedef typename SumTable::abs_domain_t summ_abs_domain_t;
    typedef summ_abs_domain_t abs_dom_t;
    typedef typename abs_dom_t::number_t number_t;    

   private:

    typedef intra_abs_transformer<abs_dom_t> intra_abs_transform_t;
    typedef typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
    typedef typename abs_dom_t::varname_t varname_t;
    typedef typename abs_dom_t::linear_expression_t linear_expression_t;    
    
    SumTable* m_sum_tbl;

  public:

    using typename abs_transform_api_t::callsite_t;
        
    static void reuse_summary (abs_dom_t& caller, 
                               const callsite_t& cs,
                               const typename SumTable::Summary& summ) {
      // Before a summary is plug-in we rename it with unique variable
      // names so we avoid naming clashes in cases like for instance
      // summary variables have same names as lhs of callsites.

      
      // error if cs and function declaration associated with summ are
      // not type consistent
      summ.check_type_consistency (cs);

      CRAB_LOG("inter", 
               crab::outs () << "    Reuse summary at " << cs << "\n";
               crab::outs () << "    Summary:" << summ << "\n";); 

      std::set<varname_t> actuals, formals;
      // --- matching formal and actual parameters
      auto inputs = summ.get_renamed_inputs ();
      unsigned i=0;
      // XXX: propagating down
      for (auto p : inputs) {
        auto a = cs.get_arg_name (i);
        if (!(a == p)) {
          CRAB_LOG ("inter",
                    crab::outs () << "\t\tPropagate from caller to callee " 
                                  << p << ":=" << a << "\n");
	  unify (caller, cs.get_arg_type(i), p, a);
        }
        ++i;
        actuals.insert (a); formals.insert (p);
      }

      // --- meet caller's inv with summ
      auto sum_inv = summ.get_renamed_sum ();
      caller = caller & sum_inv;
      CRAB_LOG("inter",
	       crab::outs() << "\t\tAfter meet caller and summary: "
	                    <<  caller << "\n");

      // --- matching callsite's lhs and callee's return value 
      // XXX: propagate from the return values in the callee to the
      // lhs variables of the callsite in the caller.
      auto const &caller_vts = cs.get_lhs ();
      auto const &callee_outs = summ.get_renamed_outputs ();
      assert (caller_vts.size () == callee_outs.size ());
	
      auto caller_it = caller_vts.begin();
      auto caller_et = caller_vts.end();
      auto callee_it = callee_outs.begin();
      
      // XXX: propagating up
      for (; caller_it != caller_et; ++caller_it, ++callee_it){
        auto vt = *caller_it;
        auto r = *callee_it;
        CRAB_LOG ("inter",
                  crab::outs () << "\t\tPropagate from callee to caller " 
		                << vt.first << ":=" << r << "\n");
	unify (caller, vt.second, vt.first, r);
        actuals.insert (vt.first); formals.insert (r);
      }

      // --- remove from caller only formal parameters so we can keep
      //     as much context from the caller as possible
      std::set<varname_t> vs;
      boost::set_difference (formals, actuals, std::inserter(vs, vs.end ()));
      domains::domain_traits<abs_dom_t>::forget (caller, vs.begin (), vs.end ());

      CRAB_LOG ("inter", 
                crab::outs() << "\t\tAfter forgetting formal parameters {";
                for (auto v: vs) crab::outs () << v << ";";
                crab::outs () << "}=" <<  caller << "\n");
    }

   public:

    bu_summ_abs_transformer(abs_dom_t* inv, SumTable* sum_tbl)
        : intra_abs_transform_t(inv), m_sum_tbl (sum_tbl) { }
    
    virtual void exec (callsite_t &cs) override
    {
      if (!m_sum_tbl)
        CRAB_WARN ("The summary table is empty: ignored analysis of callsite");
      else
      {
	if (m_sum_tbl->hasSummary (cs)) {
	  auto summ = m_sum_tbl->get (cs);
	  reuse_summary (*(this->m_inv), cs, summ);
	}
	else
	{
	  CRAB_LOG("inter",
		   crab::outs() << "\tSummary not found for " << cs << "\n");
	  for (auto vt: cs.get_lhs())
	    *(this->m_inv) -= vt.first;  // havoc 
	}
      }
    }
  }; 

  /// Conversion between domains
  template<typename Domain>
  inline void convert_domains (Domain from, Domain& to) { to = from; }

  template<typename Domain1, typename Domain2>
  inline void convert_domains (Domain1 from, Domain2& to) {
    CRAB_LOG("inter",
	     crab::outs () << "Converting from "
	                   << Domain1::getDomainName () << " to "
	                   << Domain2::getDomainName () << " might lose precision if "
	                   << Domain1::getDomainName () << " is more precise than "
	                   << Domain2::getDomainName () << "\n");
	     
    for (auto cst : from.to_linear_constraint_system ())
    { to += cst; }
  }

  /**
   * Abstract transformer specialized for performing top-down forward
   * traversal while reusing summaries at the callsites.
   * class CallCtxTable stores the calling context.
   **/
  template<class SumTable, class CallCtxTable> 
  class td_summ_abs_transformer: 
      public intra_abs_transformer<typename CallCtxTable::abs_domain_t> {

   public:

    typedef typename SumTable::abs_domain_t     summ_abs_domain_t;
    typedef typename CallCtxTable::abs_domain_t call_abs_domain_t;
    typedef call_abs_domain_t abs_dom_t;
    typedef typename abs_dom_t::number_t number_t;

   private:
    
    typedef intra_abs_transformer<abs_dom_t> intra_abs_transform_t;
    
   public:
    
    typedef typename intra_abs_transform_t::abs_transform_api_t abs_transform_api_t;
    using typename abs_transform_api_t::callsite_t;

   private:
    
    SumTable* m_sum_tbl;
    CallCtxTable* m_call_tbl;

   public:
    
    td_summ_abs_transformer (abs_dom_t* inv,
			     SumTable* sum_tbl, CallCtxTable* call_tbl)
        : intra_abs_transform_t(inv), 
          m_sum_tbl (sum_tbl), m_call_tbl (call_tbl) { }
          
    virtual void exec (callsite_t &cs) override
    {      
      if (!m_sum_tbl) {
        CRAB_WARN ("The summary table is empty");
      } else if (m_sum_tbl->hasSummary (cs)) {
        auto summ = m_sum_tbl->get (cs);

        // error if cs and function declaration associated with summ
        // are not type consistent
        summ.check_type_consistency (cs);
        
        CRAB_LOG ("inter", 
                  crab::outs () << "    Pluging caller context into callee\n"
		                << "    Summary: " << summ << "\n");
	
        ///////
        /// Generate the callee context and store it.
        ///////
	typedef typename abs_dom_t::varname_t varname_t;
	
        abs_dom_t callee_ctx_inv (*(this->m_inv));
        // --- matching formal and actual parameters
        // XXX: propagating down 	
        unsigned i=0;
        auto inputs = summ.get_inputs ();
        for (varname_t p : inputs) {
          varname_t a = cs.get_arg_name (i);
          if (!(a == p))
	    unify (callee_ctx_inv, cs.get_arg_type(i), p, a);
          ++i;
        }

        // --- project only onto formal parameters
        domains::domain_traits<abs_dom_t>::project (callee_ctx_inv, 
                                                    inputs.begin (),
                                                    inputs.end ());
        // --- store the callee context
        CRAB_LOG ("inter", 
                  crab::outs() << "\t\tCallee context stored: " 
                               << callee_ctx_inv << "\n");
        m_call_tbl->insert (cs, callee_ctx_inv);          
        
        /////
        // Generate the continuation at the caller
        /////
        
        // --- convert this->m_inv to the language of summ_abs_dom_t (summ)
        summ_abs_domain_t caller_ctx_inv;
        convert_domains(*(this->m_inv), caller_ctx_inv);
	// forget the variables of the lhs of the callsite, otherwise
	// caller_ctx_inv and m_inv may be inconsistent if the lhs
	// variables are constrained before the callsite (e.g., x=-5;
	// x := abs(x);)
	for (auto vt: cs.get_lhs())
	{ *(this->m_inv) -= vt.first; }
	
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller context: " <<  caller_ctx_inv << "\n");
        
        // --- reuse summary to do the continuation
        bu_summ_abs_transformer<SumTable>::
	  reuse_summary (caller_ctx_inv, cs, summ);
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller context after plugin summary: " 
                               << caller_ctx_inv << "\n");
	
        // --- convert back inv to the language of abs_dom_t
        convert_domains(caller_ctx_inv, *(this->m_inv));
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller continuation after " 
  		               <<  cs << "=" <<  *(this->m_inv) << "\n");
        return;
      } else {
        // --- no summary found: do nothing
      }
        
      // We could not reuse a summary so we just havoc lhs of the call
      // site.      
      for (auto vt: cs.get_lhs())
      { *(this->m_inv) -= vt.first; }
    }
  }; 

  } // end namespace
} // end namespace
#endif 
