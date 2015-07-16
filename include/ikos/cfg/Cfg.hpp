#ifndef CFG_HPP
#define CFG_HPP

/* 
 * Build a CFG to interface with IKOS
 * 
 * The CFG can support the modelling of: 
 * - integers,
 * - C-like pointers, and 
 * - arrays as any sequence of contiguous bytes either in the heap or
 *   stack.
 */

#include <stdlib.h> 

#include <boost/shared_ptr.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_set.hpp>

#include <ikos/common/types.hpp>
#include <ikos/common/bignums.hpp>
#include <ikos/algorithms/linear_constraints.hpp>
#include <ikos/domains/intervals.hpp>

namespace cfg_impl 
{
  // To convert a basic block label to a string
  template< typename T >
  inline std::string get_label_str(T e);
} 

namespace cfg 
{

  using namespace ikos;
  using namespace std;

  // The values must be such that REG <= PTR <= MEM
  enum TrackedPrecision { REG = 0, PTR = 1, MEM = 2 };

  enum VariableType { INT_TYPE, PTR_TYPE, ARR_TYPE, UNK_TYPE};

  template<typename Number, typename VariableName>
  inline std::ostream& operator<< (std::ostream &o, 
                                   const ikos::variable<Number, VariableName> &v)
  {
    auto tmp (v);
    tmp.write (o);
    return o;
  }

  template<typename Number, typename VariableName>
  inline std::ostream& operator<< (std::ostream &o, 
                                   const ikos::linear_expression<Number, VariableName> &e)
  {
    auto tmp (e);
    tmp.write (o);
    return o;
  }

  template<typename Number, typename VariableName>
  inline std::ostream& operator<< (std::ostream &o, 
                                   const ikos::linear_constraint<Number, VariableName> &c)
  {
    auto tmp (c);
    tmp.write (o);
    return o;
  }

  inline std::ostream& operator<< (std::ostream& o, VariableType t)
  {
    switch (t)
    {
      case INT_TYPE: o << "int"; break;
      case PTR_TYPE: o << "ptr"; break;
      case ARR_TYPE: o << "arr"; break;
      default: o << "unknown"; break;
    }
    return o;
  }

  template< typename VariableName>
  class Live
  {

    typedef std::vector < VariableName > live_set_t;

   public:

    typedef typename live_set_t::const_iterator  const_use_iterator;
    typedef typename live_set_t::const_iterator  const_def_iterator;

   private:

    live_set_t m_uses;
    live_set_t m_defs;

    void add (live_set_t & s, VariableName v)
    {
      auto it = std::find (s.begin (), s.end (), v);
      if (it == s.end ()) s.push_back (v);
    }

   public:

    Live() { }

    void addUse(const VariableName v){ add (m_uses,v);}
    void addDef(const VariableName v){ add (m_defs,v);}

    const_use_iterator uses_begin() const { return m_uses.begin (); }
    const_use_iterator uses_end()   const { return m_uses.end (); }
    const_use_iterator defs_begin() const { return m_defs.begin (); }
    const_use_iterator defs_end()   const { return m_defs.end (); }

    friend std::ostream& operator<<(std::ostream &o, 
                                    const Live< VariableName> &live )
    {
      o << "Use={"; 
      for (auto const& v: boost::make_iterator_range (live.uses_begin (), 
                                                      live.uses_end ()))
        o << v << ",";
      o << "} Def={"; 
      for (auto const& v: boost::make_iterator_range (live.defs_begin (), 
                                                      live.defs_end ()))
        o << v << ",";
      o << "}";
      return o;
    }

  };


  template< typename VariableName>
  struct StatementVisitor;

  template< class VariableName>
  class Statement
  {

   public:
    typedef Live<VariableName> live_t ;
    
   protected:
    live_t m_live;
    
   public:
    
    live_t getLive() const { return m_live; }
    
    virtual void accept(StatementVisitor< VariableName> *) = 0;

    virtual void write(ostream& o) const = 0 ;

    virtual boost::shared_ptr<Statement <VariableName> > clone () const = 0;

    virtual ~Statement() { }

    friend ostream& operator <<(ostream&o, const Statement<VariableName> &s)
    {
      s.write (o);
      return o;
    }

  }; 

  /*
    Basic Statements (INT_TYPE)
  */

  template< class Number, class VariableName>
  class BinaryOp: public Statement <VariableName>
  {
    
   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;

   private:

    variable_t          m_lhs;
    operation_t         m_op;
    linear_expression_t m_op1;
    linear_expression_t m_op2;

   public:
    
    BinaryOp (variable_t lhs, 
              operation_t op, 
              linear_expression_t op1, 
              linear_expression_t op2): 

        m_lhs(lhs), m_op(op), m_op1(op1), m_op2(op2) 
    { 
      this->m_live.addDef (m_lhs.name());
      for (auto v: m_op1.variables()){ this->m_live.addUse (v.name()); }         
      for (auto v: m_op2.variables()){ this->m_live.addUse (v.name()); }         
    }

    variable_t lhs () const { return m_lhs; }
    
    operation_t op () const { return m_op; }
    
    linear_expression_t left () const { return m_op1; }
    
    linear_expression_t right () const { return m_op2; }

    virtual void accept(StatementVisitor < VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef BinaryOp <Number, VariableName> BinaryOp_t;
      return boost::static_pointer_cast< Statement <VariableName>, BinaryOp_t >
          (boost::shared_ptr <BinaryOp_t> (new BinaryOp_t(m_lhs, m_op, m_op1, m_op2)));
    }

    virtual void write (ostream& o) const
    {
      o << m_lhs << " = " << m_op1;
      switch (m_op) 
      {
        case OP_ADDITION:       { o << "+"; break; }
        case OP_MULTIPLICATION: { o << "*"; break; } 
        case OP_SUBTRACTION:    { o << "-"; break; }
        case OP_DIVISION:       { o << "/"; break; }
      }
      o << m_op2 ; // << " " << this->m_live;
      return;
    }

  }; 

  template< class Number, class VariableName>
  class Assignment: public Statement<VariableName>
  {

   public:

    typedef variable< Number, VariableName >          variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    
   private:

    variable_t          m_lhs;
    linear_expression_t m_rhs;
    
   public:

    Assignment (variable_t lhs, linear_expression_t rhs): 
        m_lhs(lhs), 
        m_rhs(rhs) 
    {
      this->m_live.addDef (m_lhs.name());
      for(auto v: m_rhs.variables()) 
        this->m_live.addUse (v.name());
    }
    
    variable_t lhs () const { return m_lhs; }
    
    linear_expression_t rhs () const { return m_rhs; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Assignment <Number, VariableName> Assignment_t;
      return boost::static_pointer_cast< Statement <VariableName>, Assignment_t >
          (boost::shared_ptr <Assignment_t> (new Assignment_t(m_lhs, m_rhs)));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = " << m_rhs; // << " " << this->m_live;
      return;
    }

  }; 
    
  template<class Number, class VariableName>
  class Assume: public Statement <VariableName>
  {

   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
      
   private:

    linear_constraint_t m_cst;
    
   public:

    Assume (linear_constraint_t cst): m_cst(cst) 
    { 
      for(auto v: cst.variables())
        this->m_live.addUse (v.name()); 
    }
    
    linear_constraint_t constraint() const { return m_cst; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Assume <Number, VariableName> Assume_t;
      return boost::static_pointer_cast< Statement <VariableName>, Assume_t >
          (boost::shared_ptr <Assume_t> (new Assume_t(m_cst)));
    }
        
    virtual void write (ostream & o) const
    {
      o << "assume (" << m_cst << ")"; //  << " " << this->m_live;
      return;
    }

  }; 

  template< class VariableName>
  class Unreachable: public Statement< VariableName> 
  {
   public:
     
    Unreachable() { }
     
    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Unreachable <VariableName> Unreachable_t;
      return boost::static_pointer_cast< Statement <VariableName>, Unreachable_t >
          (boost::shared_ptr <Unreachable_t> (new Unreachable_t()));
    }
    
    virtual void write(ostream& o) const
    {
      o << "unreachable";
      return;
    }

  }; 

  template< class VariableName>
  class Havoc: public Statement< VariableName> 
  {

    VariableName m_lhs;
    
   public:

    Havoc (VariableName lhs): m_lhs(lhs) 
    { 
      this->m_live.addDef (m_lhs);
    }
     
    VariableName variable () const { return m_lhs; }
     
    virtual void accept (StatementVisitor<VariableName> *v) 
    {
      v->visit (*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Havoc <VariableName> Havoc_t;
      return boost::static_pointer_cast< Statement <VariableName>, Havoc_t >
          (boost::shared_ptr <Havoc_t> (new Havoc_t(m_lhs)));
    }
     
    void write (ostream& o) const
    {
      o << m_lhs << " =*" << " "; // << this->m_live;
      return;
    }
    
  }; 

  // select x, c, e1, e2:
  //    if c > 0 then x=e1 else x=e2
  //
  // Note that a select instruction is not strictly needed and can be
  // simulated by splitting blocks. However, frontends like LLVM can
  // generate many select instructions so we prefer to support
  // natively to avoid a blow up in the size of the CFG.
  template< class Number, class VariableName>
  class Select: public Statement <VariableName>
  {
    
   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;

   private:

    variable_t          m_lhs;
    variable_t          m_cond;
    linear_expression_t m_e1;
    linear_expression_t m_e2;

   public:
    
    Select (variable_t lhs, 
            variable_t cond, 
            linear_expression_t e1, 
            linear_expression_t e2): 

        m_lhs(lhs), m_cond(cond), m_e1(e1), m_e2(e2) 
    { 
      this->m_live.addDef (m_lhs.name());
      this->m_live.addUse (m_cond.name());
      for (auto v: m_e1.variables()){ this->m_live.addUse (v.name()); }         
      for (auto v: m_e2.variables()){ this->m_live.addUse (v.name()); }         
    }

    variable_t lhs () const { return m_lhs; }
    
    variable_t cond () const { return m_cond; }
    
    linear_expression_t left () const { return m_e1; }
    
    linear_expression_t right () const { return m_e2; }

    virtual void accept(StatementVisitor < VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Select <Number, VariableName> Select_t;
      return boost::static_pointer_cast< Statement <VariableName>, Select_t >
          (boost::shared_ptr <Select_t> (new Select_t(m_lhs, m_cond, m_e1, m_e2)));
    }

    virtual void write (ostream& o) const
    {
      o << m_lhs << " = " 
        << "(" << m_cond << " ? " << m_e1 << ":" << m_e2 << ")";
      return;
    }

  }; 

  /*
     Array statements (ARR_TYPE)
  */

  template< class VariableName>
  class ArrayInit: public Statement< VariableName> 
  {

    VariableName m_arr;
    
   public:

    ArrayInit (VariableName arr): m_arr (arr) { }
     
    VariableName variable () const { return m_arr; }
     
    virtual void accept (StatementVisitor<VariableName> *v) 
    {
      v->visit (*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef ArrayInit <VariableName> ArrayInit_t;
      return boost::static_pointer_cast< Statement <VariableName>, ArrayInit_t >
          (boost::shared_ptr <ArrayInit_t> (new ArrayInit_t(m_arr)));
    }
     
    void write (ostream& o) const
    {
      o << "array_init (" << m_arr << ")";
      return;
    }
    
  }; 

  template< class Number, class VariableName>
  class ArrayStore: public Statement<VariableName>
  {

   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    
   private:

    variable_t m_arr;
    linear_expression_t m_index;
    linear_expression_t m_value;
    bool m_is_singleton; //! whether the store writes to a singleton
                         // cell. If unknown set to false.
    
   public:

    ArrayStore (variable_t arr, linear_expression_t index, 
                linear_expression_t value, bool is_sing): 
        m_arr (arr), m_index (index),
        m_value (value), m_is_singleton (is_sing)
    {
      this->m_live.addUse (m_arr.name());
      for(auto v: m_index.variables()) 
        this->m_live.addUse (v.name());
      for(auto v: m_value.variables()) 
        this->m_live.addUse (v.name());
    }
    
    variable_t array () const { return m_arr; }
    
    linear_expression_t index () const { return m_index; }

    linear_expression_t value () const { return m_value; }

    bool is_singleton () const { return m_is_singleton;}

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef ArrayStore <Number, VariableName> array_store_t;
      return boost::static_pointer_cast< Statement <VariableName>, array_store_t>
          (boost::shared_ptr <array_store_t> (new array_store_t (m_arr, m_index,
                                                                 m_value, m_is_singleton)));
    }
    
    virtual void write(ostream& o) const
    {
      o << "array_store(" << m_arr << "," << m_index << "," << m_value << ")"; 
      return;
    }

  }; 

  template< class Number, class VariableName>
  class ArrayLoad: public Statement<VariableName>
  {

   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    
   private:

    variable_t m_lhs;
    variable_t m_array;
    linear_expression_t m_index;

   public:

    ArrayLoad (variable_t lhs, variable_t arr, linear_expression_t index): 
        m_lhs (lhs), m_array (arr), m_index (index)
    {
      this->m_live.addDef (lhs.name());
      this->m_live.addUse (m_array.name());
      for(auto v: m_index.variables()) 
        this->m_live.addUse (v.name());
    }
    
    variable_t lhs () const { return m_lhs; }

    variable_t array () const { return m_array; }

    linear_expression_t index () const { return m_index; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef ArrayLoad <Number, VariableName> array_load_t;
      return boost::static_pointer_cast< Statement <VariableName>, array_load_t>
          (boost::shared_ptr <array_load_t> (new array_load_t (m_lhs, m_array, m_index)));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = " 
        << "array_load(" << m_array << "," << m_index  << ")"; 
      return;
    }

  }; 
  
  /*
     Pointer statements (PTR_TYPE)
  */

  template< class Number, class VariableName>
  class PtrLoad: public Statement<VariableName>
  {
    // p = *q
   public:

    typedef interval <Number> interval_t;

   private:

    VariableName m_lhs;
    VariableName m_rhs;
    interval_t   m_size; //! bytes read

   public:

    PtrLoad (VariableName lhs, VariableName rhs, interval_t size): 
        m_lhs (lhs), m_rhs (rhs), m_size (size)
    {
      this->m_live.addUse (lhs);
      this->m_live.addUse (rhs);
    }
    
    VariableName lhs () const { return m_lhs; }

    VariableName rhs () const { return m_rhs; }

    interval_t size () const { return m_size; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef PtrLoad <Number,VariableName> ptr_load_t;
      return boost::static_pointer_cast< Statement <VariableName>, ptr_load_t >
          (boost::shared_ptr <ptr_load_t> (new ptr_load_t (m_lhs, m_rhs, m_size)));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = "  << "*(" << m_rhs << " + ";
      interval_t sz (m_size) ; // FIX: write method is not const in ikos domains
      sz.write (o);
      o << ")";
      return;
    }

  }; 

  template<class Number, class VariableName>
  class PtrStore: public Statement<VariableName>
  {

    // *p = q
   public:

    typedef interval <Number> interval_t;

   private:

    VariableName m_lhs;
    VariableName m_rhs;
    interval_t   m_size; //! bytes written
   
   public:

    PtrStore (VariableName lhs, VariableName rhs, interval_t size): 
        m_lhs (lhs), m_rhs (rhs), m_size (size)
    {
      this->m_live.addUse (lhs);
      this->m_live.addUse (rhs);
    }
    
    VariableName lhs () const { return m_lhs; }

    VariableName rhs () const { return m_rhs; }

    interval_t size () const { return m_size; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef PtrStore <Number,VariableName> ptr_store_t;
      return boost::static_pointer_cast< Statement <VariableName>, ptr_store_t >
          (boost::shared_ptr <ptr_store_t> (new ptr_store_t (m_lhs, m_rhs, m_size)));
    }
    
    virtual void write(ostream& o) const
    {
      o << "*(" << m_lhs << " + ";
      interval_t sz (m_size) ; // FIX: write method is not const in ikos domains
      sz.write (o);
      o << ") = "  << m_rhs;
      return;
    }

  }; 

  template< class Number, class VariableName>
  class PtrAssign: public Statement<VariableName>
  {
    //! p = q + n
   public:

    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;

   private:

    VariableName m_lhs;
    VariableName m_rhs;
    linear_expression_t m_offset;
    
   public:

    PtrAssign (VariableName lhs, VariableName rhs, linear_expression_t offset): 
        m_lhs (lhs), m_rhs (rhs), m_offset(offset)
    {
      this->m_live.addDef (lhs);
      this->m_live.addUse (rhs);
    }
    
    VariableName lhs () const { return m_lhs; }

    VariableName rhs () const { return m_rhs; }

    linear_expression_t offset () const { return m_offset; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef PtrAssign <Number, VariableName> ptr_assign_t;
      return boost::static_pointer_cast< Statement <VariableName>, ptr_assign_t >
          (boost::shared_ptr <ptr_assign_t> (new ptr_assign_t (m_lhs, m_rhs, m_offset)));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = "  << m_rhs << " + ";
      linear_expression_t off (m_offset) ; // FIX: write method is not const in ikos 
      off.write (o);
      return;
    }

  }; 

  template<class VariableName>
  class PtrObject: public Statement<VariableName>
  {
    //! lhs = &a;
    VariableName m_lhs;
    index_t m_address;
    
   public:

    PtrObject (VariableName lhs, index_t address): 
        m_lhs (lhs), m_address (address)
    {
      this->m_live.addDef (lhs);
    }
    
    VariableName lhs () const { return m_lhs; }

    index_t rhs () const { return m_address; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef PtrObject <VariableName> ptr_object_t;
      return boost::static_pointer_cast< Statement <VariableName>, ptr_object_t>
          (boost::shared_ptr <ptr_object_t> (new ptr_object_t (m_lhs, m_address)));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = "  << "&(" << m_address << ")" ;
      return;
    }

  }; 

  template<class VariableName>
  class PtrFunction: public Statement<VariableName>
  {
    // lhs = &func;
    VariableName m_lhs;
    index_t m_func;
    
   public:

    PtrFunction (VariableName lhs, VariableName func): 
        m_lhs (lhs), m_func (func)
    {
      this->m_live.addDef (lhs);
    }
    
    VariableName lhs () const { return m_lhs; }

    index_t rhs () const { return m_func; }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef PtrFunction <VariableName> ptr_function_t;
      return boost::static_pointer_cast< Statement <VariableName>, ptr_function_t>
          (boost::shared_ptr <ptr_function_t> (new ptr_function_t ( lhs (), rhs ())));
    }
    
    virtual void write(ostream& o) const
    {
      o << m_lhs << " = "  << "&(" << m_func << ")" ;
      return;
    }
  }; 

  /*
     Function calls
  */

  template< class VariableName>
  class FCallSite: public Statement<VariableName>
  {

    boost::optional<pair<VariableName,VariableType> > m_lhs;
    VariableName m_func_name;
    vector<pair<VariableName,VariableType> > m_args;

    typedef typename vector<pair<VariableName,VariableType> >::iterator arg_iterator;
    typedef typename vector<pair<VariableName,VariableType> >::const_iterator const_arg_iterator;
    
   public:

    FCallSite (VariableName func_name, 
               vector<pair<VariableName,VariableType> > args): 
        m_func_name (func_name)
    {
      std::copy (args.begin (), args.end (), std::back_inserter (m_args));
      for (auto arg:  m_args) { this->m_live.addUse (arg.first); }
    }

    FCallSite (VariableName func_name, vector<VariableName> args): 
        m_func_name (func_name)
    {
      for (auto v : args)
      {
        m_args.push_back (make_pair (v, UNK_TYPE));
        this->m_live.addUse (v);
      }
    }
    
    FCallSite (pair<VariableName,VariableType> lhs, 
               VariableName func_name, 
               vector<pair<VariableName,VariableType> > args): 
        m_lhs (boost::optional<pair<VariableName,VariableType> > (lhs)), 
        m_func_name (func_name)
    {
      std::copy (args.begin (), args.end (), std::back_inserter(m_args));
      for (auto arg:  m_args) { this->m_live.addUse (arg.first); }
      this->m_live.addDef ((*m_lhs).first);
    }

    FCallSite (VariableName lhs, 
               VariableName func_name, vector<VariableName> args): 
        m_lhs (boost::optional<pair<VariableName,VariableType> > (lhs, UNK_TYPE)), 
        m_func_name (func_name)
    {
      for (auto v : args)
      {
        m_args.push_back (make_pair (v, UNK_TYPE));
        this->m_live.addUse (v);
      }
      this->m_live.addDef ((*m_lhs).first);
    }
    

    boost::optional<VariableName> get_lhs_name () const { 
      if (m_lhs)
        return boost::optional<VariableName> ((*m_lhs).first);
      else 
        return boost::optional<VariableName> ();
    }

    VariableType get_lhs_type () const {       
      if (m_lhs) return (*m_lhs).second;
      else return UNK_TYPE;
    }
    
    VariableName get_func_name () const { 
      return m_func_name; 
    }

    unsigned get_num_args () const { return m_args.size (); }

    VariableName get_arg_name (unsigned idx) const { 
      if (idx >= m_args.size ())
        IKOS_ERROR ("Out-of-bound access to call site parameter");

      return m_args[idx].first;
    }

    VariableType get_arg_type (unsigned idx) const { 
      if (idx >= m_args.size ())
        IKOS_ERROR ("Out-of-bound access to call site parameter");
      
      return m_args[idx].second;
    }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef FCallSite <VariableName> call_site_t;
      if (m_lhs)
        return boost::static_pointer_cast< Statement <VariableName>, call_site_t>
            (boost::shared_ptr <call_site_t> (new call_site_t (*m_lhs, m_func_name, m_args)));
      else
        return boost::static_pointer_cast< Statement <VariableName>, call_site_t>
            (boost::shared_ptr <call_site_t> (new call_site_t (m_func_name, m_args)));
    }
    
    virtual void write(ostream& o) const
    {
      if (m_lhs)
        o << (*m_lhs).first << " = ";

      o << "call " << m_func_name << "(";
      for (const_arg_iterator It = m_args.begin (), Et=m_args.end (); It!=Et; )
      {
        o << It->first;
        ++It;
        if (It != Et)
          o << ",";
      }
      o << ")";

      return;
    }

  }; 

  template< class VariableName>
  class Return: public Statement<VariableName>
  {

    VariableName m_var;
    VariableType m_type;

   public:

    Return (VariableName var, VariableType type = UNK_TYPE): 
        m_var (var), m_type (type)
    {
      this->m_live.addUse (m_var); 
    }
        
    VariableName get_ret_var () const { 
      return m_var;
    }

    VariableType get_ret_type () const { 
      return m_type;
    }

    virtual void accept(StatementVisitor <VariableName> *v) 
    {
      v->visit(*this);
    }

    virtual boost::shared_ptr<Statement <VariableName> > clone () const
    {
      typedef Return <VariableName> return_t;
      return boost::static_pointer_cast< Statement <VariableName>, return_t>
          (boost::shared_ptr <return_t> (new return_t (m_var, m_type)));
    }
    
    virtual void write(ostream& o) const
    {
      o << "return " << m_var;
      return;
    }
  }; 
  
  template< class BasicBlockLabel, class VariableName>
  class Cfg;

  template< class BasicBlockLabel, class VariableName >
  class BasicBlock: public boost::noncopyable
  {

    // TODO: support for removing statements 
    friend class Cfg< BasicBlockLabel, VariableName >;
    
   public:

    // helper types to build statements
    typedef variable< z_number, VariableName > z_variable_t;
    typedef linear_expression< z_number, VariableName > z_lin_exp_t;
    typedef linear_constraint< z_number, VariableName > z_lin_cst_t;

    typedef Statement< VariableName> statement_t;
    typedef BasicBlock< BasicBlockLabel, VariableName > basic_block_t;
    
   private:

    typedef std::vector< BasicBlockLabel >   bb_id_set_t;

    typedef boost::shared_ptr< statement_t > statement_ptr;
    typedef std::vector< statement_ptr > stmt_list_t;

   public:

    typedef typename bb_id_set_t::iterator succ_iterator;
    typedef typename bb_id_set_t::const_iterator const_succ_iterator;
    typedef succ_iterator pred_iterator;
    typedef const_succ_iterator const_pred_iterator;

    typedef boost::indirect_iterator< typename stmt_list_t::iterator > iterator;
    typedef boost::indirect_iterator< typename stmt_list_t::const_iterator > const_iterator;

   private:

    typedef interval <z_number> z_interval;

    // Basic statements
    typedef BinaryOp<z_number,VariableName> z_bin_op_t;
    typedef Assignment<z_number,VariableName> z_assign_t;
    typedef Assume<z_number,VariableName> z_assume_t;
    typedef Havoc<VariableName> havoc_t;
    typedef Unreachable<VariableName> unreach_t;
    typedef Select<z_number,VariableName> z_select_t;
    // Functions
    typedef FCallSite<VariableName> callsite_t;
    typedef Return<VariableName> return_t;
    // Arrays
    typedef ArrayInit<VariableName> arr_init_t;
    typedef ArrayStore<z_number,VariableName> z_arr_store_t;
    typedef ArrayLoad<z_number,VariableName> z_arr_load_t;
    // Pointers
    typedef PtrStore<z_number,VariableName> ptr_store_t;
    typedef PtrLoad<z_number,VariableName> ptr_load_t;
    typedef PtrAssign<z_number,VariableName> ptr_assign_t;
    typedef PtrObject<VariableName> ptr_object_t;
    typedef PtrFunction<VariableName> ptr_function_t;

    typedef boost::shared_ptr<z_bin_op_t> z_bin_op_ptr;
    typedef boost::shared_ptr<z_assign_t> z_assign_ptr;
    typedef boost::shared_ptr<z_assume_t> z_assume_ptr;
    typedef boost::shared_ptr<havoc_t> havoc_ptr;      
    typedef boost::shared_ptr<unreach_t> unreach_ptr;
    typedef boost::shared_ptr<z_select_t> z_select_ptr;
    typedef boost::shared_ptr<callsite_t> callsite_ptr;      
    typedef boost::shared_ptr<return_t> return_ptr;      
    typedef boost::shared_ptr<arr_init_t> arr_init_ptr;
    typedef boost::shared_ptr<z_arr_store_t> z_arr_store_ptr;
    typedef boost::shared_ptr<z_arr_load_t> z_arr_load_ptr;    
    typedef boost::shared_ptr<ptr_store_t> ptr_store_ptr;
    typedef boost::shared_ptr<ptr_load_t> ptr_load_ptr;    
    typedef boost::shared_ptr<ptr_assign_t> ptr_assign_ptr;
    typedef boost::shared_ptr<ptr_object_t> ptr_object_ptr;    
    typedef boost::shared_ptr<ptr_function_t> ptr_function_ptr;    
   
    BasicBlockLabel m_bb_id;
    stmt_list_t m_stmts;
    bb_id_set_t m_prev, m_next;
    TrackedPrecision m_track_prec;    
    // Ideally it should be size_t to indicate any position within the
    // block. For now, we only allow to insert either at front or at
    // the back (default). Note that if insertions at the front are
    // very common we should replace stmt_list_t from a vector to a
    // deque.
    bool m_insert_point_at_front; 

    void InsertAdjacent (bb_id_set_t &c, BasicBlockLabel e)
    { 
      if (std::find(c.begin (), c.end (), e) == c.end ())
        c.push_back (e);
    }
    
    void RemoveAdjacent (bb_id_set_t &c, BasicBlockLabel e)
    {
      if (std::find(c.begin (), c.end (), e) != c.end ())
        c.erase (std::remove(c.begin (), c.end (), e), c.end ());
    }
    
    BasicBlock (BasicBlockLabel bb_id, TrackedPrecision track_prec): 
        m_bb_id (bb_id), m_track_prec (track_prec), 
        m_insert_point_at_front (false)
    { }

    static boost::shared_ptr< basic_block_t > Create (BasicBlockLabel bb_id, 
                                                     TrackedPrecision track_prec) 
    {
      return boost::shared_ptr< basic_block_t > (new basic_block_t (bb_id, track_prec));
    }
    
    void insert(statement_ptr stmt) 
    {
      if (m_insert_point_at_front)
      { 
        m_stmts.insert (m_stmts.begin(), stmt);
        m_insert_point_at_front = false;
      }
      else
        m_stmts.push_back(stmt);
    }

   public:

    //! it will be set to false after the first insertion
    void set_insert_point_front (){
      m_insert_point_at_front = true;
    }

    boost::shared_ptr <basic_block_t> clone () const
    {
      boost::shared_ptr <basic_block_t> b (new basic_block_t (label (), 
                                                            m_track_prec));

      for (auto &s : boost::make_iterator_range (begin (), end ()))
        b->m_stmts.push_back (s.clone ()); 

      for (auto id : boost::make_iterator_range (prev_blocks ())) 
        b->m_prev.push_back (id);
      
      for (auto id : boost::make_iterator_range (next_blocks ()))
        b->m_next.push_back (id);
      
      return b;
    }

    BasicBlockLabel label () const { return m_bb_id; }
    
    iterator begin()             { return boost::make_indirect_iterator (m_stmts.begin ()); }
    iterator end()               { return boost::make_indirect_iterator (m_stmts.end ()); }
    const_iterator begin() const { return boost::make_indirect_iterator (m_stmts.begin ()); }
    const_iterator end()   const { return boost::make_indirect_iterator (m_stmts.end ()); }
    
    std::size_t size() { return std::distance ( begin (), end ()); }
        
    pair<succ_iterator, succ_iterator> next_blocks ()
    { 
      return std::make_pair (m_next.begin (), m_next.end ());
    }
    
    pair<pred_iterator, pred_iterator> prev_blocks() 
    { 
      return std::make_pair (m_prev.begin (), m_prev.end ());
    }

    pair<const_succ_iterator,const_succ_iterator> next_blocks () const
    { 
      return std::make_pair (m_next.begin (), m_next.end ());
    }
    
    pair<const_pred_iterator,const_pred_iterator> prev_blocks() const
    { 
      return std::make_pair (m_prev.begin (), m_prev.end ());
    }

    void reverse()
    {
      std::swap (m_prev, m_next);
      std::reverse (m_stmts.begin (), m_stmts.end ());
    }
    
    // Add a cfg edge from *this to b
    void operator>>(basic_block_t& b) 
    {
      InsertAdjacent (m_next, b.m_bb_id);
      InsertAdjacent (b.m_prev, m_bb_id);
    }
    
    // Remove a cfg edge from *this to b
    void operator-=(basic_block_t &b)
    {
      RemoveAdjacent (m_next, b.m_bb_id);
      RemoveAdjacent (b.m_prev, m_bb_id);       
    }
    
    // insert all statements of other at the front
    void merge_front (const basic_block_t &other) 
    {
      m_stmts.insert (m_stmts.begin (), 
                      other.m_stmts.begin (), 
                      other.m_stmts.end ());
    }
    
    // insert all statements of other at the back
    void merge_back (const basic_block_t &other) 
    {
      m_stmts.insert (m_stmts.end (), 
                      other.m_stmts.begin (), 
                      other.m_stmts.end ());
    }
    
    void write(ostream& o) const
    {
      o << cfg_impl::get_label_str (m_bb_id) << ":\n";	

      for (auto const &s: *this)
        o << "  " << s << ";\n"; 
      
      o << "--> [";

      for (auto const &n : boost::make_iterator_range (next_blocks ()))
        o << cfg_impl::get_label_str (n) << ";";

      o << "]\n";
      
      return;
    }

    /// To build statements

    void add (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
    {
      insert(boost::static_pointer_cast< statement_t, z_bin_op_t >
             (z_bin_op_ptr(new z_bin_op_t(lhs, OP_ADDITION, op1, op2))));
    }
    
    void add (z_variable_t lhs, z_variable_t op1, z_number op2) 
    {
      insert(boost::static_pointer_cast< statement_t, z_bin_op_t > 
             (z_bin_op_ptr(new z_bin_op_t(lhs, OP_ADDITION, op1,  op2))));
    }
    
    void add (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
    {
      if  (op1.get_variable () && op2.get_variable ())
        add (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        add (lhs, (*op1.get_variable ()), op2.constant ());

      else if (op1.is_constant () && op2.get_variable ())
        add (lhs, (*op2.get_variable ()), op1.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        z_lin_exp_t rhs(z_number(op1.constant () + op2.constant ()));
        insert(boost::static_pointer_cast< statement_t, z_assign_t >
               (z_assign_ptr (new z_assign_t (lhs, rhs))));
      }
      else
        assert(false && "add operands unexpected");
    }
    
    void sub (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
    {
      insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
              (z_bin_op_ptr (new z_bin_op_t (lhs, OP_SUBTRACTION, op1, op2))));
    }
    
    void sub (z_variable_t lhs, z_variable_t op1, z_number op2) 
    {
      insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
              (z_bin_op_ptr (new z_bin_op_t (lhs, OP_SUBTRACTION, op1, op2))));
    }
    
    void sub (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        sub (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        sub (lhs, (*op1.get_variable ()), op2.constant ());        

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        z_lin_exp_t rhs (z_number (op1.constant () - op2.constant ()));
        insert (boost::static_pointer_cast< statement_t, z_assign_t >
                (z_assign_ptr(new z_assign_t (lhs, rhs))));
      }
      else
        assert(false && "sub operands unexpected");
    }
    
    void mul (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
    {
      insert(boost::static_pointer_cast< statement_t, z_bin_op_t >
             (z_bin_op_ptr(new z_bin_op_t (lhs, OP_MULTIPLICATION, op1, op2))));
    }
    
    void mul (z_variable_t lhs, z_variable_t op1, z_number op2) 
    {
      insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
              (z_bin_op_ptr(new z_bin_op_t (lhs, OP_MULTIPLICATION, op1, op2))));
    }
    
    void mul(z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        mul (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        mul (lhs, (*op1.get_variable()), op2.constant());

      else if (op1.is_constant () && op2.get_variable ())
        mul (lhs, (*op2.get_variable ()), op1.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        z_lin_exp_t rhs(z_number(op1.constant () * op2.constant ()));
        insert(boost::static_pointer_cast< statement_t, z_assign_t >
               (z_assign_ptr(new z_assign_t (lhs, rhs))));
      }
      else
        assert(false && "mul operands unexpected");
    }
    
    void div (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
    {
      insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
              (z_bin_op_ptr (new z_bin_op_t (lhs, OP_DIVISION, op1, op2))));
    }
    
    void div (z_variable_t lhs, z_variable_t op1, z_number op2) 
    {
      insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
              (z_bin_op_ptr (new z_bin_op_t (lhs, OP_DIVISION, op1, op2))));
    }
    
    void div (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        div (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        div (lhs, (*op1.get_variable ()), op2.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        z_lin_exp_t rhs (z_number (op1.constant() / op2.constant()));
        insert (boost::static_pointer_cast< statement_t, z_assign_t >
                (z_assign_ptr (new z_assign_t (lhs, rhs))));
      }
      else
        assert (false && "div operands unexpected");
    }
    
    void assign (z_variable_t lhs, z_lin_exp_t rhs) 
    {
      insert (boost::static_pointer_cast< statement_t, z_assign_t >
              (z_assign_ptr (new z_assign_t (lhs, rhs))));
    }
    
    void assume (z_lin_cst_t cst) 
    {
      insert (boost::static_pointer_cast< statement_t, z_assume_t >
              (z_assume_ptr (new z_assume_t (cst))));
    }
    
    void havoc(VariableName lhs) 
    {
      insert (boost::static_pointer_cast< statement_t, havoc_t > 
              (havoc_ptr (new havoc_t (lhs))));
    }
    
    void unreachable() 
    {
      insert (boost::static_pointer_cast< statement_t, unreach_t > 
              (unreach_ptr (new unreach_t ())));
    }

    void select (z_variable_t lhs, z_variable_t cond, z_number n1, z_number n2) 
    {
      insert(boost::static_pointer_cast< statement_t, z_select_t >
             (z_select_ptr(new z_select_t (lhs, cond, n1, n2))));
    }

    void select (z_variable_t lhs, z_variable_t cond, z_number n, z_variable_t v) 
    {
      insert(boost::static_pointer_cast< statement_t, z_select_t >
             (z_select_ptr(new z_select_t (lhs, cond, n, v))));
    }

    void select (z_variable_t lhs, z_variable_t cond, z_variable_t v, z_number n) 
    {
      insert(boost::static_pointer_cast< statement_t, z_select_t >
             (z_select_ptr(new z_select_t (lhs, cond, v, n))));
    }

    void select (z_variable_t lhs, z_variable_t cond, z_variable_t v1, z_variable_t v2) 
    {
      insert(boost::static_pointer_cast< statement_t, z_select_t >
             (z_select_ptr(new z_select_t (lhs, cond, v1, v2))));
    }
    
    void callsite (VariableName func, 
                   vector<pair <VariableName,VariableType> > args) 
    {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (callsite_ptr (new callsite_t(func, args))));
    }

    void callsite (VariableName func, vector<VariableName> args) 
    {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (callsite_ptr (new callsite_t(func, args))));
    }

    void callsite (pair<VariableName,VariableType> lhs, 
                   VariableName func, 
                   vector<pair <VariableName,VariableType> > args) 
    {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (callsite_ptr (new callsite_t(lhs, func, args))));
    }

    void callsite (VariableName lhs, VariableName func, vector<VariableName> args) 
    {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (callsite_ptr (new callsite_t(lhs, func, args))));
    }

    void ret (VariableName var, VariableType ty) 
    {
        insert(boost::static_pointer_cast< statement_t, return_t >
               (return_ptr (new return_t(var, ty))));
    }

    void ret (VariableName var) 
    {
        insert(boost::static_pointer_cast< statement_t, return_t >
               (return_ptr (new return_t(var))));
    }

    void array_init(VariableName arr) 
    {
      insert (boost::static_pointer_cast< statement_t, arr_init_t > 
              (arr_init_ptr (new arr_init_t (arr))));
    }

    void array_store (z_variable_t arr, z_lin_exp_t idx, 
                      z_lin_exp_t val, bool is_singleton = false) 
    {
      if (m_track_prec == MEM)
        insert(boost::static_pointer_cast< statement_t, z_arr_store_t >
               (z_arr_store_ptr (new z_arr_store_t(arr, idx, val, 
                                                   is_singleton))));
    }

    void array_load (z_variable_t lhs, z_variable_t arr, z_lin_exp_t idx) 
    {
      if (m_track_prec == MEM)
        insert(boost::static_pointer_cast< statement_t, z_arr_load_t >
               (z_arr_load_ptr (new z_arr_load_t(lhs, arr, idx))));
    }

    void ptr_store (VariableName lhs, VariableName rhs, z_interval size) 
    {
      if (m_track_prec >= PTR)
        insert(boost::static_pointer_cast< statement_t, ptr_store_t >
               (ptr_store_ptr (new ptr_store_t (lhs, rhs, size))));
    }

    void ptr_load (VariableName lhs, VariableName rhs, z_interval size) 
    {
      if (m_track_prec >= PTR)
        insert(boost::static_pointer_cast< statement_t, ptr_load_t >
               (ptr_load_ptr (new ptr_load_t (lhs, rhs, size))));
    }

    void ptr_assign (VariableName lhs, VariableName rhs, z_lin_exp_t offset) 
    {
      if (m_track_prec >= PTR)
        insert(boost::static_pointer_cast< statement_t, ptr_assign_t >
               (ptr_assign_ptr (new ptr_assign_t (lhs, rhs, offset))));
    }

    void new_object (VariableName lhs, index_t address) 
    {
      if (m_track_prec >= PTR)
        insert(boost::static_pointer_cast< statement_t, ptr_object_t >
               (ptr_object_ptr (new ptr_object_t (lhs, address))));
    }

    void new_ptr_func (VariableName lhs, index_t func) 
    {
      if (m_track_prec >= PTR)
        insert(boost::static_pointer_cast< statement_t, ptr_function_t >
               (ptr_function_ptr (new ptr_function_t (lhs, func))));
    }

    friend ostream& operator<<(ostream &o, const basic_block_t &b)
    {
      b.write (o);
      return o;
    }
    
  }; 

  template< class VariableName>
  struct StatementVisitor
  {
    typedef linear_expression< z_number, VariableName > z_lin_exp_t;
    typedef BinaryOp <z_number,VariableName> z_bin_op_t;
    typedef Assignment <z_number,VariableName> z_assign_t;
    typedef Assume <z_number,VariableName> z_assume_t;
    typedef Havoc<VariableName> havoc_t;
    typedef Unreachable<VariableName> unreach_t;
    typedef Select<z_number,VariableName> z_select_t;

    typedef FCallSite<VariableName> callsite_t;
    typedef Return<VariableName> return_t;

    typedef ArrayInit<VariableName> arr_init_t;
    typedef ArrayStore<z_number,VariableName> z_arr_store_t;
    typedef ArrayLoad<z_number,VariableName> z_arr_load_t;

    typedef PtrStore<z_number,VariableName> ptr_store_t;
    typedef PtrLoad<z_number,VariableName> ptr_load_t;
    typedef PtrAssign<z_number,VariableName> ptr_assign_t;
    typedef PtrObject<VariableName> ptr_object_t;
    typedef PtrFunction<VariableName> ptr_function_t;

    // Only implementation for basic statements is required

    virtual void visit (z_bin_op_t&) = 0;
    virtual void visit (z_assign_t&) = 0;
    virtual void visit (z_assume_t&) = 0;
    virtual void visit (havoc_t&) = 0;
    virtual void visit (unreach_t&) = 0;
    virtual void visit (z_select_t&) = 0;

    virtual void visit (callsite_t&) { };
    virtual void visit (return_t&) { };
    virtual void visit (arr_init_t&) { };
    virtual void visit (z_arr_store_t&) { };
    virtual void visit (z_arr_load_t&) { };
    virtual void visit (ptr_store_t&) { };
    virtual void visit (ptr_load_t&) { };
    virtual void visit (ptr_assign_t&) { };
    virtual void visit (ptr_object_t&) { };
    virtual void visit (ptr_function_t&) { };

    virtual ~StatementVisitor () { }
  }; 

  template< class VariableName>
  class FunctionDecl
  {

    VariableType m_lhs_type;
    VariableName m_func_name;
    vector<pair<VariableName,VariableType> > m_params;

    typedef typename vector<pair<VariableName,VariableType> >::iterator param_iterator;
    typedef typename vector<pair<VariableName,VariableType> >::const_iterator const_param_iterator;
    
   public:

    FunctionDecl (VariableName func_name, vector<VariableName> params): 
        m_lhs_type (UNK_TYPE), m_func_name (func_name)
    {
      for (auto v : params)
        m_params.push_back (make_pair (v, UNK_TYPE));
    }
    
    FunctionDecl (VariableType lhs_type, VariableName func_name, 
                  vector<pair<VariableName,VariableType> > params): 
        m_lhs_type (lhs_type), m_func_name (func_name)
    {
      std::copy (params.begin (), params.end (), std::back_inserter (m_params));
    }

    VariableType get_lhs_type () const { return m_lhs_type; }
    
    VariableName get_func_name () const { 
      return m_func_name; 
    }

    unsigned get_num_params () const { return m_params.size (); }

    VariableName get_param_name (unsigned idx) const { 
      if (idx >= m_params.size ())
        IKOS_ERROR ("Out-of-bound access to function parameter");

      return m_params[idx].first;
    }

    VariableType get_param_type (unsigned idx) const { 
      if (idx >= m_params.size ())
        IKOS_ERROR ("Out-of-bound access to function parameter");

      return m_params[idx].second;
    }
    
    void write(ostream& o) const
    {
      o << m_lhs_type << " declare " << m_func_name << "(";
      for (const_param_iterator It = m_params.begin (), Et=m_params.end (); It!=Et; )
      {
        o << It->first << ":" << It->second;
        ++It;
        if (It != Et)
          o << ",";
      }
      o << ")";

      return;
    }

    friend ostream& operator<<(ostream& o, const FunctionDecl<VariableName> &decl)
    { 
      decl.write (o);
      return o;
    }
  }; 

  template< class BasicBlockLabel, class VariableName >
  class Cfg
  {
   public:

    typedef BasicBlockLabel basic_block_label_t;
    typedef VariableName varname_t;

    typedef BasicBlock< BasicBlockLabel, VariableName > basic_block_t;   
    typedef Statement < VariableName > statement_t;

    typedef typename basic_block_t::succ_iterator succ_iterator;
    typedef typename basic_block_t::pred_iterator pred_iterator;

    typedef boost::iterator_range <succ_iterator> succ_range;
    typedef boost::iterator_range <pred_iterator> pred_range;

   private:

    typedef Cfg < BasicBlockLabel, VariableName > cfg_t;
    typedef boost::shared_ptr< basic_block_t > basic_block_ptr;
    typedef boost::unordered_map< BasicBlockLabel, basic_block_ptr > basic_block_map;
    typedef typename basic_block_map::value_type binding_t;
    typedef boost::shared_ptr< basic_block_map > basic_block_map_ptr;

    struct getRef : public std::unary_function<binding_t, basic_block_t>
    {
      getRef () { }
      basic_block_t& operator () (const binding_t &p) const { return *(p.second); }
    }; 

    struct getLabel : public std::unary_function<binding_t, BasicBlockLabel>
    {
      getLabel () { }
      BasicBlockLabel operator () (const binding_t &p) const { return p.second->label(); }
    }; 

   public:

    typedef boost::transform_iterator< getRef, 
                                       typename basic_block_map::iterator > iterator;
    typedef boost::transform_iterator< getRef, 
                                       typename basic_block_map::const_iterator > const_iterator;
    typedef boost::transform_iterator< getLabel, 
                                       typename basic_block_map::iterator > label_iterator;
    typedef boost::transform_iterator< getLabel, 
                                       typename basic_block_map::const_iterator > const_label_iterator;

   private:

    BasicBlockLabel    m_entry;
    BasicBlockLabel    m_exit;
    bool               m_has_exit;
    basic_block_map_ptr  m_blocks;
    TrackedPrecision   m_track_prec;
    //! if Cfg is associated with a function
    boost::optional<FunctionDecl<VariableName> > m_func_decl; 


    typedef boost::unordered_set< BasicBlockLabel > visited_t;
    template<typename T>
    void dfs_rec (BasicBlockLabel curId, visited_t &visited, T f) const
    {
      if (visited.find (curId) != visited.end ()) return;
      visited.insert (curId);

      const basic_block_t &cur = get_node (curId);
      f (cur);
      for (auto const n : boost::make_iterator_range (cur.next_blocks ()))
        dfs_rec (n, visited, f);
    }

    template<typename T>
    void dfs (T f) const 
    {
      visited_t visited;
      dfs_rec (m_entry, visited, f);
    }

    struct PrintBlock 
    {
      ostream &m_o;
      PrintBlock (ostream& o) : m_o (o) { }
      void operator () (const basic_block_t& B){ m_o << B ; }
    };

    Cfg (): m_blocks (basic_block_map_ptr (new basic_block_map)) { }

   public:

    Cfg (BasicBlockLabel entry, 
         TrackedPrecision track_prec = REG):  
        m_entry  (entry), 
        m_has_exit (false),
        m_blocks (basic_block_map_ptr (new basic_block_map)),
        m_track_prec (track_prec)
    {
      m_blocks->insert (binding_t (m_entry, 
                                   basic_block_t::Create (m_entry, m_track_prec)));
    }

    Cfg (BasicBlockLabel entry, 
         BasicBlockLabel exit, 
         TrackedPrecision track_prec = REG):  
        m_entry  (entry), 
        m_exit   (exit), 
        m_has_exit (true),
        m_blocks (basic_block_map_ptr (new basic_block_map)),
        m_track_prec (track_prec)
    {
      m_blocks->insert (binding_t (m_entry, 
                                   basic_block_t::Create (m_entry, m_track_prec)));
    }

    Cfg (BasicBlockLabel entry, 
         BasicBlockLabel exit, 
         FunctionDecl<VariableName> func_decl, 
         TrackedPrecision track_prec = REG):  
        m_entry (entry), 
        m_exit (exit), 
        m_has_exit (true),
        m_blocks (basic_block_map_ptr (new basic_block_map)),
        m_track_prec (track_prec),
        m_func_decl (boost::optional<FunctionDecl<VariableName> > (func_decl))
    {
      m_blocks->insert (binding_t (m_entry, 
                                   basic_block_t::Create (m_entry, m_track_prec)));
    }

    //! copy constructor will make shallow copies so use this method
    //! for deep copies.
    cfg_t clone () const
    {
      cfg_t cfg;

      cfg.m_entry = m_entry ;
      if (m_has_exit)
        cfg.m_exit = m_exit ;
      cfg.m_has_exit = m_has_exit ;
      cfg.m_func_decl = m_func_decl;
      for (auto const &BB: boost::make_iterator_range (begin (), end ()))
      {
        boost::shared_ptr <basic_block_t> copyBB = BB.clone ();
        cfg.m_blocks->insert (binding_t (copyBB->label (), copyBB));
      }
      return cfg;
    }

    boost::optional<FunctionDecl<VariableName> > get_func_decl () const { 
      return m_func_decl; 
    }

    TrackedPrecision get_track_prec () const {
      return m_track_prec;
    }

    BasicBlockLabel entry() const { return m_entry; } 
    
    bool has_exit () const { return m_has_exit; }

    BasicBlockLabel exit()  const { 
      if (has_exit ()) return m_exit; 
      
      IKOS_ERROR ("Cfg does not have an exit block");
    } 

    //! set method to mark the exit block after the Cfg has been
    //! created.
    void set_exit (BasicBlockLabel exit) { 
      m_exit = exit; 
      m_has_exit = true;
    }

    //! set method to add the function declaration after the Cfg has
    //! been created.
    void set_func_decl (FunctionDecl<VariableName> decl) { 
      m_func_decl  = boost::optional<FunctionDecl<VariableName> > (decl);
    }

    //! required by ikos fixpoint
    succ_range next_nodes (BasicBlockLabel bb_id) 
    {
      
      basic_block_t& b = get_node(bb_id);
      return boost::make_iterator_range (b.next_blocks ());
    }
            
    //! required by ikos fixpoint
    pred_range prev_nodes (BasicBlockLabel bb_id) 
    {
      basic_block_t& b = get_node(bb_id);
      return boost::make_iterator_range (b.prev_blocks ());
    }
    
    basic_block_t& insert (BasicBlockLabel bb_id) 
    {
      auto it = m_blocks->find (bb_id);
      if (it != m_blocks->end ()) return *(it->second);

      basic_block_ptr block = basic_block_t::Create (bb_id, m_track_prec);
      m_blocks->insert (binding_t (bb_id, block));
      return *block;
    }

    void remove (BasicBlockLabel bb_id)
    {
      basic_block_t& bb = get_node(bb_id) ;

      vector< pair<basic_block_t*,basic_block_t*> > dead;
      
      for (auto id : boost::make_iterator_range (bb.prev_blocks ()))
      { 
        if (bb_id != id)
        {
          basic_block_t& p = get_node(id) ;
          dead.push_back (make_pair (&p,&bb));
        }
      }
      
      for (auto id : boost::make_iterator_range (bb.next_blocks ()))
      {
        if (bb_id != id)
        {
          basic_block_t& s = get_node(id) ;
          dead.push_back (make_pair (&bb,&s));
        }
      }

      for (auto p : dead)
        (*p.first) -= (*p.second);

      m_blocks->erase (bb_id);
    }

    
    basic_block_t& get_node (BasicBlockLabel bb_id) 
    {
      auto it = m_blocks->find (bb_id);
      if (it == m_blocks->end ())
        IKOS_ERROR ("Basic block not found in the CFG");

      return *(it->second);
    }

    const basic_block_t& get_node (BasicBlockLabel bb_id) const
    {
      auto it = m_blocks->find (bb_id);
      if (it == m_blocks->end ())
        IKOS_ERROR ("Basic block not found in the CFG");

      return *(it->second);
    }
    
    //! return a begin iterator of BasicBlock's
    iterator begin() 
    {
      return boost::make_transform_iterator (m_blocks->begin (), getRef ());
    }
    
    //! return an end iterator of BasicBlock's
    iterator end() 
    {
      return boost::make_transform_iterator (m_blocks->end (), getRef ());
    }

    const_iterator begin() const
    {
      return boost::make_transform_iterator (m_blocks->begin (), getRef ());
    }
    
    const_iterator end() const
    {
      return boost::make_transform_iterator (m_blocks->end (), getRef ());
    }

    //! return a begin iterator of BasicBlockLabel's
    label_iterator label_begin() 
    {
      return boost::make_transform_iterator (m_blocks->begin (), getLabel ());
    }
    
    //! return an end iterator of BasicBlockLabel's
    label_iterator label_end() 
    {
      return boost::make_transform_iterator (m_blocks->end (), getLabel ());
    }

    const_label_iterator label_begin() const
    {
      return boost::make_transform_iterator (m_blocks->begin (), getLabel ());
    }
    
    const_label_iterator label_end() const
    {
      return boost::make_transform_iterator (m_blocks->end (), getLabel ());
    }

    size_t size () const { return std::distance (begin (), end ()); }
        
    void reverse()
    {
      if (!m_has_exit)
        IKOS_ERROR ("Cfg cannot be reversed: no exit block found");

      std::swap (m_entry, m_exit);
      for (auto &p: *this) 
        p.reverse (); 
    }

    void write (ostream& o) const
    {
      PrintBlock f (o);
      if (m_func_decl)
        o << *m_func_decl << endl;
      dfs (f);
      return;
    }

    friend ostream& operator<<(ostream &o, 
                               const Cfg< BasicBlockLabel, VariableName > &cfg)
    {
      cfg.write (o);
      return o;
    }

    void simplify ()
    {
      mergeBlocks ();        
      removeUnreachableBlocks ();
      removeUselessBlocks ();
      //after removing useless blocks there can be opportunities to
      //merge more blocks.
      mergeBlocks ();
      mergeBlocks ();
    }
    
   private:
    
    ////
    // Cfg simplifications
    ////
    
    struct HasAssertVisitor: public StatementVisitor<VariableName>
    {
      typedef typename StatementVisitor<VariableName>::z_bin_op_t z_bin_op_t;
      typedef typename StatementVisitor<VariableName>::z_assign_t z_assign_t;
      typedef typename StatementVisitor<VariableName>::z_assume_t z_assume_t;
      typedef typename StatementVisitor<VariableName>::havoc_t havoc_t;
      typedef typename StatementVisitor<VariableName>::unreach_t unreach_t;
      typedef typename StatementVisitor<VariableName>::z_select_t z_select_t;

      bool _has_assert;
      HasAssertVisitor (): _has_assert(false) { }
      void visit(z_bin_op_t&){ }  
      void visit(z_assign_t&) { }
      void visit(z_assume_t&) { _has_assert = true; }
      void visit(havoc_t&) { }
      void visit(unreach_t&){ }
      void visit(z_select_t&){ }
    };
    
    // Helpers
    bool hasOneChild (BasicBlockLabel b)
    {
      auto rng = next_nodes (b);
      return (std::distance (rng.begin (), rng.end ()) == 1);
    }
    
    bool hasOneParent (BasicBlockLabel b)
    {
      auto rng = prev_nodes (b);
      return (std::distance (rng.begin (), rng.end ()) == 1);
    }
    
    basic_block_t& getChild (BasicBlockLabel b)
    {
      assert (hasOneChild (b));
      auto rng = next_nodes (b);
      return get_node (*(rng.begin ()));
    }
    
    basic_block_t& getParent (BasicBlockLabel b)
    {
      assert (hasOneParent (b));
      auto rng = prev_nodes (b);
      return get_node (*(rng.begin ()));
    }
    
    void mergeBlocksRec (BasicBlockLabel curId, 
                         visited_t& visited)
    {
      
      if (visited.find (curId) != visited.end ()) return;
      visited.insert (curId);
      
      basic_block_t &cur = get_node (curId);
      
      if (hasOneChild (curId) && hasOneParent (curId))
      {
        basic_block_t &parent = getParent (curId);
        basic_block_t &child  = getChild (curId);
        
        HasAssertVisitor vis;
        for (auto it = cur.begin (); it != cur.end (); ++it)
          it->accept(&vis);
        
        // we don't merge two blocks with assertions
        if (!vis._has_assert)
        {
          parent.merge_back (cur);
          remove (curId);
          parent >> child;        
          mergeBlocksRec (child.label (), visited); 
          return;
        }
      }
      
      for (auto n : boost::make_iterator_range (cur.next_blocks ()))
        mergeBlocksRec (n, visited);
    }
    
    // Merges a basic block into its predecessor if there is only one
    // and the predecessor only has one successor.
    void mergeBlocks ()
    {
      visited_t visited;
      mergeBlocksRec (entry (), visited);
    }
    
    // mark reachable blocks from curId
    void markAliveBlocks (BasicBlockLabel curId, 
                          cfg_t& cfg,
                          visited_t& visited)
    {
      if (visited.count (curId) > 0) return;
      visited.insert (curId);
      for (auto child : cfg.next_nodes (curId))
        markAliveBlocks (child, cfg, visited);
    }
    
    // remove unreachable blocks
    void removeUnreachableBlocks ()
    {
      visited_t alive, dead;
      markAliveBlocks (entry (), *this, alive);
      
      for (auto const &bb : *this) 
        if (!(alive.count (bb.label ()) > 0))
          dead.insert (bb.label ());
      
      for (auto bb_id: dead)
        remove (bb_id);
    }
    
    // remove blocks that cannot reach the exit block
    void removeUselessBlocks ()
    {
      if (!has_exit ()) return;

      cfg_t cfg  = clone ();
      cfg.reverse ();
      
      visited_t useful, useless;
      markAliveBlocks (cfg.entry (), cfg, useful);
      
      for (auto const &bb : *this) 
        if (!(useful.count (bb.label ()) > 0))
          useless.insert (bb.label ());
      
      for (auto bb_id: useless)
        remove (bb_id);
    }
    
  }; 
  
} // end namespace 

#endif /* CFG_HPP */

