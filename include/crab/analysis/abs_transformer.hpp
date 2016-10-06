#ifndef ABSTRACT_TRANSFORMER_HPP
#define ABSTRACT_TRANSFORMER_HPP

/* 
   Implementation of the abstract transfer functions by reducing them
   to abstract domain operations.
   
   These are the main Crab statements for which we define their abstract
   transfer functions:
   
   INTEGERS 
     x := y bin_op z;
     x := y; 
     assume (cst) / assert (cst);
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

   INTEGERS/ARRAYS/POINTERS
     havoc (x);

  Comments: note that all abstract transfer functions assume z_number
  as the *only* type for numbers. Thus, changes are required to
  support other types, e.g., reals.
 */

#include <boost/optional.hpp>
#include "boost/range/algorithm/set_algorithm.hpp"

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/domains/linear_constraints.hpp>

namespace crab {

  namespace analyzer {

  using namespace cfg;
  using namespace std;

  //! API abstract transformer
  template<typename VariableName>
  class abs_transformer_api: public statement_visitor <VariableName>
  {
   public:

    typedef VariableName varname_t;
    typedef variable < z_number, VariableName > z_var_t;
    typedef linear_expression< z_number, VariableName > z_lin_exp_t;
    typedef linear_constraint< z_number, VariableName > z_lin_cst_t;
    typedef linear_constraint_system< z_number, VariableName > z_lin_cst_sys_t;

    typedef havoc_stmt<VariableName>                 havoc_t;
    typedef binary_op <z_number,VariableName>        z_bin_op_t;
    typedef assignment <z_number,VariableName>       z_assign_t;
    typedef assume_stmt <z_number,VariableName>      z_assume_t;
    typedef unreachable_stmt<VariableName>           unreach_t;
    typedef select_stmt <z_number,VariableName>      z_select_t;
    typedef assert_stmt <z_number,VariableName>      z_assert_t;
    typedef callsite_stmt<VariableName>              callsite_t;
    typedef return_stmt<VariableName>                return_t;
    typedef array_assume_stmt<z_number,VariableName> arr_assume_t;
    typedef array_store_stmt<z_number,VariableName>  arr_store_t;
    typedef array_load_stmt<z_number,VariableName>   arr_load_t;
    typedef array_assign_stmt<VariableName>          arr_assign_t;
    typedef ptr_store_stmt<VariableName>             ptr_store_t;
    typedef ptr_load_stmt<VariableName>              ptr_load_t;
    typedef ptr_assign_stmt<z_number,VariableName>   z_ptr_assign_t;
    typedef ptr_object_stmt<VariableName>            ptr_object_t;
    typedef ptr_function_stmt<VariableName>          ptr_function_t;
    typedef ptr_null_stmt<VariableName>              ptr_null_t;
    typedef ptr_assume_stmt<VariableName>            ptr_assume_t;
    typedef ptr_assert_stmt<VariableName>            ptr_assert_t;

   protected: 

    virtual void exec (havoc_t&) { }
    virtual void exec (z_bin_op_t&)  { } 
    virtual void exec (z_assign_t&) { }
    virtual void exec (z_assume_t&) { }
    virtual void exec (unreach_t&) { }
    virtual void exec (z_select_t&) { }
    virtual void exec (z_assert_t&) { }
    virtual void exec (callsite_t&) { }
    virtual void exec (return_t&) { }
    virtual void exec (arr_assume_t&) { }
    virtual void exec (arr_store_t&) { }
    virtual void exec (arr_load_t&) { }
    virtual void exec (arr_assign_t&) { }
    virtual void exec (ptr_store_t&) { }
    virtual void exec (ptr_load_t&) { }
    virtual void exec (z_ptr_assign_t&) { }
    virtual void exec (ptr_object_t&) { }
    virtual void exec (ptr_function_t&) { }
    virtual void exec (ptr_null_t&) { }
    virtual void exec (ptr_assume_t&) { }
    virtual void exec (ptr_assert_t&) { }

   public: /* visitor api */

    void visit (havoc_t &s) { exec (s); }
    void visit (z_bin_op_t &s) { exec (s); }
    void visit (z_assign_t &s) { exec (s); }
    void visit (z_assume_t &s) { exec (s); }
    void visit (unreach_t &s) { exec (s); }
    void visit (z_select_t &s) { exec (s); }
    void visit (z_assert_t &s) { exec (s); }
    void visit (callsite_t &s) { exec (s); }
    void visit (return_t &s) { exec (s); }
    void visit (arr_assume_t &s) { exec (s); }
    void visit (arr_store_t &s) { exec (s); }
    void visit (arr_load_t &s) { exec (s); }
    void visit (arr_assign_t &s) { exec (s); }
    void visit (ptr_store_t &s) { exec (s); }
    void visit (ptr_load_t &s) { exec (s); }
    void visit (z_ptr_assign_t &s) { exec (s); }
    void visit (ptr_object_t &s) { exec (s); }
    void visit (ptr_function_t &s) { exec (s); }
    void visit (ptr_null_t &s) { exec (s); }
    void visit (ptr_assume_t &s) { exec (s); }
    void visit (ptr_assert_t &s) { exec (s); }

  };


  //! Abstract transformer specialized for numerical abstract domains
  //! with arrays and pointers.
  // XXX: num_abs_transformer should be renamed to
  //      abs_transformer_impl or something like that.
  template<typename NumAbsDomain, 
           typename SumTable /*unused*/, typename CallCtxTable /*unused*/>
  class num_abs_transformer: 
        public abs_transformer_api <typename CallCtxTable::abs_domain_t::varname_t>
  {

   public:

    typedef typename SumTable::abs_domain_t summ_abs_domain_t;
    typedef typename CallCtxTable::abs_domain_t call_abs_domain_t;
    typedef NumAbsDomain abs_dom_t;

   public:

    typedef abs_transformer_api <typename CallCtxTable::abs_domain_t::varname_t> abs_transform_t;
    using typename abs_transform_t::varname_t;
    using typename abs_transform_t::z_var_t;
    using typename abs_transform_t::z_lin_exp_t;
    using typename abs_transform_t::z_lin_cst_t;
    using typename abs_transform_t::z_lin_cst_sys_t;

    using typename abs_transform_t::havoc_t;
    using typename abs_transform_t::z_bin_op_t;
    using typename abs_transform_t::z_assign_t;
    using typename abs_transform_t::z_assume_t;
    using typename abs_transform_t::z_select_t;
    using typename abs_transform_t::z_assert_t;
    using typename abs_transform_t::unreach_t;
    using typename abs_transform_t::callsite_t;
    using typename abs_transform_t::arr_assume_t;
    using typename abs_transform_t::arr_load_t;
    using typename abs_transform_t::arr_store_t;
    using typename abs_transform_t::arr_assign_t;
    using typename abs_transform_t::ptr_load_t;
    using typename abs_transform_t::ptr_store_t;
    using typename abs_transform_t::ptr_assign_t;
    using typename abs_transform_t::ptr_object_t;
    using typename abs_transform_t::ptr_function_t;
    using typename abs_transform_t::ptr_null_t;
    using typename abs_transform_t::ptr_assume_t;
    using typename abs_transform_t::ptr_assert_t;

   protected:

    abs_dom_t& m_inv;

    template<typename T>
    void apply (abs_dom_t& inv, binary_operation_t op,
                varname_t x, varname_t y, T z) {

      auto op1 = conv_op<ikos::operation_t> (op);
      auto op2 = conv_op<ikos::div_operation_t> (op);
      auto op3 = conv_op<ikos::bitwise_operation_t> (op);
      
      if (op1) 
        inv.apply (*op1, x, y, z);
      else if (op2) 
        inv.apply (*op2, x, y, z);
      else if (op3)
        inv.apply (*op3, x, y, z);
      else
        CRAB_ERROR("unsupported binary operator", op);
    }
    
   public:

    num_abs_transformer (abs_dom_t& inv): 
        m_inv (inv) { }
    
    num_abs_transformer (abs_dom_t& inv, SumTable*, CallCtxTable*): 
        m_inv (inv) { }

    abs_dom_t& inv () { return m_inv; }

    void exec (z_bin_op_t& stmt)  {
      
      z_lin_exp_t op1 = stmt.left ();
      z_lin_exp_t op2 = stmt.right ();
      
      if (op1.get_variable () && op2.get_variable ())
      {
        apply (m_inv, stmt.op (), 
               stmt.lhs ().name(), 
               (*op1.get_variable ()).name(), 
               (*op2.get_variable ()).name());
      }
      else
      {
        assert ( op1.get_variable () && op2.is_constant ());
        apply (m_inv, stmt.op (), 
               stmt.lhs ().name (), 
               (*op1.get_variable ()).name(), 
               op2.constant ()); 
      }      
    }

    void exec (z_select_t& stmt) {
    
      abs_dom_t inv1 (m_inv);
      abs_dom_t inv2 (m_inv);
      inv1 += stmt.cond ();
      inv2 += stmt.cond ().negate ();
      if (inv2.is_bottom()) {
        inv1.assign (stmt.lhs().name (),stmt.left());
        m_inv = inv1;
      }
      else if (inv1.is_bottom ()) {
        inv2.assign (stmt.lhs().name (),stmt.right());
        m_inv = inv2;
      }
      else {
        inv1.assign (stmt.lhs().name (),stmt.left());
        inv2.assign (stmt.lhs().name (),stmt.right());
        m_inv = inv1 | inv2;
      }
    }    
    
    void exec (z_assign_t& stmt) {
      m_inv.assign (stmt.lhs().name (), z_lin_exp_t (stmt.rhs()));
    }
    
    void exec (z_assume_t& stmt) {
      m_inv += stmt.constraint();
    }

    void exec (havoc_t& stmt)  {
      m_inv -= stmt.variable();
    }

    void exec (unreach_t& stmt) {
      m_inv = abs_dom_t::bottom ();
    }

    // --- default implementation is intra-procedural
    void exec (callsite_t &cs) {
      for (auto vt: cs.get_lhs())
        m_inv -= vt.first;  // havoc
    }

    void exec (z_assert_t& stmt) {
      abs_dom_t cst;
      cst += stmt.constraint();
      abs_dom_t meet = cst & m_inv;
      if (meet.is_bottom ()) {
        m_inv = abs_dom_t::bottom (); // assertion does not definitely hold.
        return;
      }
      m_inv += stmt.constraint ();
    }

    // arrays 

    void exec (arr_assume_t &stmt) {
      m_inv.array_assume (stmt.array (), stmt.array_type (), stmt.var ());
    }
    
    void exec (arr_store_t &stmt) {
      m_inv.array_store (stmt.array(), stmt.array_type (), 
                         stmt.index (), stmt.value (),
                         stmt.elem_size(), stmt.is_singleton ());
    }

    void exec (arr_load_t  &stmt) {
      m_inv.array_load (stmt.lhs (), 
                        stmt.array (), stmt.array_type (), 
                        stmt.index(), stmt.elem_size());
    }

    void exec (arr_assign_t  &stmt) {
      m_inv.array_assign (stmt.lhs (), stmt.rhs (), stmt.array_type ());
    }

    // pointers 

    void exec (ptr_null_t & stmt) { 
      m_inv.pointer_mk_null (stmt.lhs ());
    }
    
    void exec (ptr_object_t & stmt) { 
      m_inv.pointer_mk_obj (stmt.lhs (), stmt.rhs());
    }
    
    void exec (ptr_function_t & stmt) { 
      m_inv.pointer_function (stmt.lhs (), stmt.rhs ());
    }
    
    void exec (ptr_assign_t & stmt) { 
      m_inv.pointer_assign (stmt.lhs (), stmt.rhs (), stmt.offset ());
    }
    
    void exec (ptr_load_t & stmt) {
      m_inv.pointer_load (stmt.lhs (), stmt.rhs ());
    }
    
    void exec (ptr_store_t & stmt) { 
      m_inv.pointer_store (stmt.lhs (), stmt.rhs ());
    }
    
    void exec (ptr_assume_t& stmt) { 
      m_inv.pointer_assume (stmt.constraint ());
    }

    void exec (ptr_assert_t& stmt) { 
      m_inv.pointer_assert (stmt.constraint ());
    }

  }; 

  //! Transformer specialized for computing numerical summaries
  template<typename SumTable, typename CallCtxTable /*unused*/>
  class bu_summ_num_abs_transformer: 
        public num_abs_transformer <typename SumTable::abs_domain_t,
                                   SumTable, CallCtxTable> {
                                   
   public:

    typedef typename SumTable::abs_domain_t summ_abs_domain_t;
    typedef typename CallCtxTable::abs_domain_t call_abs_domain_t;
    typedef summ_abs_domain_t abs_dom_t;

   public:

    typedef num_abs_transformer <abs_dom_t, SumTable, CallCtxTable> num_abs_transform_t;
    typedef typename num_abs_transform_t::abs_transform_t abs_transform_t;

    using typename abs_transform_t::varname_t;
    using typename abs_transform_t::z_var_t;
    using typename abs_transform_t::z_lin_exp_t;

    using typename abs_transform_t::havoc_t;
    using typename abs_transform_t::z_bin_op_t;
    using typename abs_transform_t::z_assign_t;
    using typename abs_transform_t::z_assume_t;
    using typename abs_transform_t::z_select_t;
    using typename abs_transform_t::unreach_t;
    using typename abs_transform_t::z_assert_t;
    using typename abs_transform_t::arr_assume_t;
    using typename abs_transform_t::arr_load_t;
    using typename abs_transform_t::arr_store_t;
    using typename abs_transform_t::arr_assign_t;
    using typename abs_transform_t::callsite_t;
    using typename abs_transform_t::ptr_load_t;
    using typename abs_transform_t::ptr_store_t;
    using typename abs_transform_t::ptr_assign_t;
    using typename abs_transform_t::ptr_object_t;
    using typename abs_transform_t::ptr_function_t;
    using typename abs_transform_t::ptr_null_t;
    using typename abs_transform_t::ptr_assume_t;
    using typename abs_transform_t::ptr_assert_t;

   private:

    SumTable* m_sum_tbl;

   public:

    static void reuse_summary (abs_dom_t& caller, 
                               const callsite_t& cs,
                               const typename SumTable::Summary& summ) {

      // error if cs and function declaration associated with summ are
      // not type consistent
      summ.check_type_consistency (cs);

      CRAB_LOG("inter", 
               crab::outs () << "    Reuse summary at " << cs << "\n";
               crab::outs () << "    Summary:";  summ.write(crab::outs ()); 
               crab::outs () << "\n");

      std::set<varname_t> actuals, formals;
      // --- matching formal and actual parameters
      auto pars = summ.get_params ();
      unsigned i=0;
      // XXX: propagating down
      for (auto p : pars) {
        auto a = cs.get_arg_name (i);

        if (!(a == p)) {

          CRAB_LOG ("inter",
                    crab::outs () << "\t\tPropagate from caller to callee " 
                                  << p << ":=" << a << "\n");

          if (cs.get_arg_type(i) == INT_TYPE) {
            caller.assign (p, z_var_t(a));
          } else if (cs.get_arg_type(i) == PTR_TYPE) {
            caller.pointer_assign (p, a, z_number(0));
          } else {
            caller.array_assign (p, a, cs.get_arg_type(i));
          }
        }
        ++i;
        actuals.insert (a); formals.insert (p);
      }

      // --- meet caller's inv with summ
      auto sum_inv = summ.get_sum ();
      caller = caller & sum_inv;
      CRAB_LOG("inter", crab::outs() << "\t\tAfter meet caller and summary: " <<  caller << "\n");

      // --- matching callsite's lhs and callee's return value 
      // XXX: propagate from the return values in the callee to the
      // lhs variables of the callsite in the caller.
      auto const &caller_vts = cs.get_lhs ();
      auto const &callee_rs = summ.get_ret_vals ();
      assert (caller_vts.size () == callee_rs.size ());
      
      auto caller_it = caller_vts.begin();
      auto caller_et = caller_vts.end();
      auto callee_it = callee_rs.begin();
      
      // XXX: propagating up
      for (; caller_it != caller_et; ++caller_it, ++callee_it){
        auto vt = *caller_it;
        auto r = *callee_it;

        CRAB_LOG ("inter",
                  crab::outs () << "\t\tPropagate from callee to caller " 
                                << vt.first << ":=" << r << "\n");
        
        if (vt.second == INT_TYPE) {
          caller.assign(vt.first, z_var_t (r));
        } else if (vt.second == PTR_TYPE) {
          caller.pointer_assign (vt.first, r, z_number (0));
        } else { 
          caller.array_assign (vt.first, r, vt.second);
        }

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

    bu_summ_num_abs_transformer(abs_dom_t& inv, SumTable* sum_tbl)
        : num_abs_transform_t(inv), 
          m_sum_tbl (sum_tbl) { }
    
    bu_summ_num_abs_transformer(abs_dom_t& inv, SumTable* sum_tbl, CallCtxTable*)
        : num_abs_transform_t(inv), 
          m_sum_tbl (sum_tbl) { }

    abs_dom_t& inv () { return this->m_inv; }    

    void exec (havoc_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_bin_op_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_select_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assign_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assume_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (unreach_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assert_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_assume_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_store_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_load_t  &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_assign_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_load_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_store_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assign_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_object_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_function_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_null_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assume_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assert_t &stmt) { num_abs_transform_t::exec (stmt); }

    // If there is a summary for the callee we use it to compute the
    // caller's summary
    void exec (callsite_t &cs) {
      if (!m_sum_tbl) {
        CRAB_WARN ("The summary table is empty: ignored analysis of callsite");
        return;
      }

      if (m_sum_tbl->hasSummary (cs)) {
        auto &summ = m_sum_tbl->get (cs);
        reuse_summary (this->m_inv, cs, summ);
      } else {
        CRAB_LOG("inter",
                 crab::outs() << "\tSummary not found for " << cs << "\n");
        for (auto vt: cs.get_lhs())
          this->m_inv -= vt.first;  // havoc 
      }
    }
  }; 

  /// Conversion between domains
  template<typename Domain>
  inline void convert_domains (Domain from, Domain& to) { to = from; }

  template<typename Domain1, typename Domain2>
  inline void convert_domains (Domain1 from, Domain2& to) {
    CRAB_WARN("Converting from " + Domain1::getDomainName () + " to " +
              Domain2::getDomainName () + " might lose precision " +
              "if " + Domain1::getDomainName () + " is an abstraction of " +
              Domain2::getDomainName ());
    for (auto cst : from.to_linear_constraint_system ())
    { to += cst; }
  }

  // Transformer specialized for performing top-down forward
  // traversal while reusing numerical summaries at the callsites  
  template<typename SumTable, typename CallCtxTable>
  class td_summ_num_abs_transformer: 
        public num_abs_transformer<typename CallCtxTable::abs_domain_t,
                                  SumTable, CallCtxTable> {
                                   
   public:

    typedef typename SumTable::abs_domain_t summ_abs_domain_t;
    typedef typename CallCtxTable::abs_domain_t call_abs_domain_t;
    typedef call_abs_domain_t abs_dom_t;

   public:

    typedef num_abs_transformer<abs_dom_t, SumTable,CallCtxTable> num_abs_transform_t;
    typedef typename num_abs_transform_t::abs_transform_t abs_transform_t;

    using typename abs_transform_t::varname_t;
    using typename abs_transform_t::z_var_t;
    using typename abs_transform_t::z_lin_exp_t;
    using typename abs_transform_t::z_lin_cst_t;
    using typename abs_transform_t::z_lin_cst_sys_t;

    using typename abs_transform_t::havoc_t;
    using typename abs_transform_t::z_bin_op_t;
    using typename abs_transform_t::z_assign_t;
    using typename abs_transform_t::z_assume_t;
    using typename abs_transform_t::z_select_t;
    using typename abs_transform_t::unreach_t;
    using typename abs_transform_t::z_assert_t;
    using typename abs_transform_t::arr_assume_t;
    using typename abs_transform_t::arr_load_t;
    using typename abs_transform_t::arr_store_t;
    using typename abs_transform_t::arr_assign_t;
    using typename abs_transform_t::callsite_t;
    using typename abs_transform_t::ptr_load_t;
    using typename abs_transform_t::ptr_store_t;
    using typename abs_transform_t::ptr_assign_t;
    using typename abs_transform_t::ptr_object_t;
    using typename abs_transform_t::ptr_function_t;
    using typename abs_transform_t::ptr_null_t;
    using typename abs_transform_t::ptr_assume_t;
    using typename abs_transform_t::ptr_assert_t;

   private:

    typedef bu_summ_num_abs_transformer<SumTable,CallCtxTable> bu_abs_transformer_t;

    SumTable* m_sum_tbl;
    CallCtxTable* m_call_tbl;

   public:
    
    td_summ_num_abs_transformer (abs_dom_t& inv, 
                                 SumTable* sum_tbl, CallCtxTable* call_tbl)
        : num_abs_transform_t(inv), 
          m_sum_tbl (sum_tbl),
          m_call_tbl (call_tbl) { }
    
    abs_dom_t& inv () { return this->m_inv; }    

    void exec (havoc_t& stmt) { num_abs_transform_t::exec (stmt); }    
    void exec (z_bin_op_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_select_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assign_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assume_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (unreach_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (z_assert_t& stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_assume_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_store_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_load_t  &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (arr_assign_t  &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_load_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_store_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assign_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_object_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_function_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_null_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assume_t &stmt) { num_abs_transform_t::exec (stmt); }
    void exec (ptr_assert_t &stmt) { num_abs_transform_t::exec (stmt); }

    void exec (callsite_t &cs) {
      
      if (!m_sum_tbl) {
        CRAB_WARN ("The summary table is empty");
      } else if (m_sum_tbl->hasSummary (cs)) {
        auto &summ = m_sum_tbl->get (cs);

        // error if cs and function declaration associated with summ
        // are not type consistent
        summ.check_type_consistency (cs);
        
        CRAB_LOG ("inter", 
                  crab::outs () << "    Pluging caller context into callee\n");
        ///////
        /// Generate the callee context and store it.
        ///////
        abs_dom_t callee_ctx_inv (this->m_inv);
        // --- matching formal and actual parameters
        unsigned i=0;
        auto pars = summ.get_params ();
        // XXX: propagating down 
        for (auto p : pars) {
          auto a = cs.get_arg_name (i);
          if (!(a == p)) {
            if (cs.get_arg_type(i) == INT_TYPE) {
              callee_ctx_inv.assign (p, z_var_t (a));
            } else if (cs.get_arg_type(i) == PTR_TYPE) { 
              callee_ctx_inv.pointer_assign (p, a, z_number (0));
            } else {
              callee_ctx_inv.array_assign (p, a, cs.get_arg_type (i));
            } 
          }
          ++i;
        }

        // --- project only onto formal parameters
        domains::domain_traits<abs_dom_t>::project (callee_ctx_inv, 
                                                    pars.begin (),
                                                    pars.end ());
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
        convert_domains(this->m_inv, caller_ctx_inv);
        
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller context: " <<  caller_ctx_inv << "\n");
        
        // --- reuse summary to do the continuation
        bu_abs_transformer_t::reuse_summary (caller_ctx_inv, cs, summ);
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller context after plugin summary: " 
                               << caller_ctx_inv << "\n");
        
        // --- convert back inv to the language of abs_dom_t
        convert_domains(caller_ctx_inv, this->m_inv);
        CRAB_LOG ("inter",
                  crab::outs() << "\t\tCaller continuation after " 
                               <<  cs << "=" <<  this->m_inv << "\n");
        return;
      } else {
        // --- no summary found: do nothing
      }
        
      // We could not reuse a summary so we just havoc lhs of the call
      // site.
      for (auto vt: cs.get_lhs())
        this->m_inv -= vt.first; 
    }
  }; 

  } // end namespace
} // end namespace
#endif 
