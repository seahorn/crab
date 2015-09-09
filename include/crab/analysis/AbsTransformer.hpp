#ifndef ABSTRACT_TRANSFORMER_HPP
#define ABSTRACT_TRANSFORMER_HPP

#include <boost/optional.hpp>

#include <crab/cfg/Cfg.hpp>
#include <crab/domains/domain_traits.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/division_operators_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>

namespace crab {

  namespace analyzer {

  using namespace cfg;
  using namespace std;

  template<typename T>
  inline boost::optional<T> convOp (binary_operation_t op); 

  template<>
  inline boost::optional<ikos::operation_t> 
  convOp (binary_operation_t op) {     
    switch (op) {
      case BINOP_ADD: return OP_ADDITION;
      case BINOP_SUB: return OP_SUBTRACTION;
      case BINOP_MUL: return OP_MULTIPLICATION;
      case BINOP_SDIV: return OP_DIVISION;
      default: return boost::optional<ikos::operation_t> ();
    }
  }
  
  template<>
  inline boost::optional<ikos::div_operation_t> 
  convOp (binary_operation_t op) {     
    switch (op) {
      case BINOP_SDIV: return OP_SDIV;
      case BINOP_UDIV: return OP_UDIV;
      case BINOP_SREM: return OP_SREM;
      case BINOP_UREM: return OP_UREM;
      default: return boost::optional<ikos::div_operation_t> ();
    }
  }

  template<>
  inline boost::optional<ikos::bitwise_operation_t> 
  convOp (binary_operation_t op) {     
    switch (op) {
      case BINOP_AND: return OP_AND;
      case BINOP_OR: return OP_OR;
      case BINOP_XOR: return OP_XOR;
      case BINOP_SHL: return OP_SHL;
      case BINOP_LSHR: return OP_LSHR;
      case BINOP_ASHR: return OP_ASHR;
      default: return boost::optional<ikos::bitwise_operation_t> ();
    }
  }
  

  //! API abstract transformer
  template<typename VariableName, typename AbsDomain>
  class AbsTransformer: public StatementVisitor <VariableName>
  {
   public:

    typedef AbsDomain abs_dom_t;

    typedef linear_expression< z_number, VariableName > z_lin_exp_t;
    typedef linear_constraint< z_number, VariableName > z_lin_cst_t;
    typedef linear_constraint_system< z_number, VariableName > z_lin_cst_sys_t;

    typedef BinaryOp <z_number,VariableName>    z_bin_op_t;
    typedef Assignment <z_number,VariableName>  z_assign_t;
    typedef Assume <z_number,VariableName>      z_assume_t;
    typedef Havoc<VariableName>                 havoc_t;
    typedef Unreachable<VariableName>           unreach_t;
    typedef Select <z_number,VariableName>      z_select_t;
    typedef FCallSite<VariableName>             callsite_t;
    typedef Return<VariableName>                return_t;
    typedef ArrayInit<VariableName>             z_arr_init_t;
    typedef AssumeArray<z_number,VariableName>  z_assume_arr_t;
    typedef ArrayStore<z_number,VariableName>   z_arr_store_t;
    typedef ArrayLoad<z_number,VariableName>    z_arr_load_t;
    typedef PtrStore<z_number,VariableName>     z_ptr_store_t;
    typedef PtrLoad<z_number,VariableName>      z_ptr_load_t;
    typedef PtrAssign<z_number,VariableName>    z_ptr_assign_t;
    typedef PtrObject<VariableName>             ptr_object_t;
    typedef PtrFunction<VariableName>           ptr_function_t;

   protected: 

    virtual void exec (z_bin_op_t&)  { } 
    virtual void exec (z_assign_t&) { }
    virtual void exec (z_assume_t&) { }
    virtual void exec (havoc_t&) { }
    virtual void exec (unreach_t&) { }
    virtual void exec (z_select_t&) { }
    virtual void exec (callsite_t&) { }
    virtual void exec (return_t&) { }
    virtual void exec (z_arr_init_t&) { }
    virtual void exec (z_assume_arr_t&) { }
    virtual void exec (z_arr_store_t&) { }
    virtual void exec (z_arr_load_t&) { }
    virtual void exec (z_ptr_store_t&) { }
    virtual void exec (z_ptr_load_t&) { }
    virtual void exec (z_ptr_assign_t&) { }
    virtual void exec (ptr_object_t&) { }
    virtual void exec (ptr_function_t&) { }

   public: /* visitor api */

    void visit (z_bin_op_t &s) { exec (s); }
    void visit (z_assign_t &s) { exec (s); }
    void visit (z_assume_t &s) { exec (s); }
    void visit (havoc_t &s) { exec (s); }
    void visit (unreach_t &s) { exec (s); }
    void visit (z_select_t &s) { exec (s); }
    void visit (callsite_t &s) { exec (s); }
    void visit (return_t &s) { exec (s); }
    void visit (z_arr_init_t &s) { exec (s); }
    void visit (z_assume_arr_t &s) { exec (s); }
    void visit (z_arr_store_t &s) { exec (s); }
    void visit (z_arr_load_t &s) { exec (s); }
    void visit (z_ptr_store_t &s) { exec (s); }
    void visit (z_ptr_load_t &s) { exec (s); }
    void visit (z_ptr_assign_t &s) { exec (s); }
    void visit (ptr_object_t &s) { exec (s); }
    void visit (ptr_function_t &s) { exec (s); }

  };

  //! Abstract transformer specialized for a numerical abstract
  //! domain.
  template<typename NumAbsDomain>
  class NumAbsTransformer: 
        public AbsTransformer <typename NumAbsDomain::varname_t, 
                               NumAbsDomain>
  {
    typedef typename NumAbsDomain::varname_t varname_t;
    typedef AbsTransformer <varname_t, NumAbsDomain> abs_transform_t;

    using typename abs_transform_t::z_lin_exp_t;
    using typename abs_transform_t::z_bin_op_t;
    using typename abs_transform_t::z_assign_t;
    using typename abs_transform_t::z_assume_t;
    using typename abs_transform_t::z_select_t;
    using typename abs_transform_t::havoc_t;
    using typename abs_transform_t::unreach_t;
    using typename abs_transform_t::z_arr_init_t;
    using typename abs_transform_t::z_assume_arr_t;
    using typename abs_transform_t::z_arr_load_t;
    using typename abs_transform_t::z_arr_store_t;

    NumAbsDomain m_inv;
    
    template<typename T>
    void apply (NumAbsDomain& inv, binary_operation_t op,
                varname_t x, varname_t y, T z) {

      auto op1 = convOp<ikos::operation_t> (op);
      auto op2 = convOp<ikos::div_operation_t> (op);
      auto op3 = convOp<ikos::bitwise_operation_t> (op);

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

    NumAbsDomain inv() const { return m_inv; }

   public:
    
    NumAbsTransformer (NumAbsDomain inv): m_inv (inv) { }
    
    void exec (z_bin_op_t& stmt) 
    {

      
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

    void exec (z_select_t& stmt) 
    {
      NumAbsDomain inv1 (m_inv);
      NumAbsDomain inv2 (m_inv);
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
    
    void exec (z_assign_t& stmt) 
    {
      m_inv.assign (stmt.lhs().name (), z_lin_exp_t (stmt.rhs()));
    }
    
    void exec (z_assume_t& stmt) 
    {
      m_inv += stmt.constraint();
    }

    void exec (havoc_t& stmt) 
    {
      m_inv -= stmt.variable();
    }

    void exec (unreach_t& stmt) 
    {
      m_inv = NumAbsDomain::bottom ();
    }

    void exec (z_arr_init_t &stmt) 
    {
      domain_traits::array_init (m_inv, 
                                 stmt.variable (), 
                                 stmt.values ());
    }

    void exec (z_assume_arr_t &stmt) 
    {
      domain_traits::assume_array (m_inv, 
                                   stmt.variable (), 
                                   stmt.val ());
    }

    void exec (z_arr_store_t &stmt) 
    {
      if (stmt.index ().get_variable ())
      {
        auto arr = stmt.array ().name ();
        auto idx = *(stmt.index ().get_variable ());
        domain_traits::array_store (m_inv, 
                                    arr,
                                    idx.name(), 
                                    stmt.value (),
                                    stmt.n_bytes (),
                                    stmt.is_singleton ());
      }
    }

    void exec (z_arr_load_t  &stmt) 
    {
      if (stmt.index ().get_variable ())
      {
        auto idx = *(stmt.index ().get_variable ());
        domain_traits::array_load (m_inv, 
                                   stmt.lhs ().name (), 
                                   stmt.array ().name (), 
                                   idx.name (),
                                   stmt.n_bytes ());
      }
    }
  }; 

  } // end namespace
} // end namespace
#endif 
