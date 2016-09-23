#ifndef CFG_HPP
#define CFG_HPP

/* 
 * Build a CFG to interface with the abstract domains and fixpoint
 * iterators.
 * 
 * The CFG can support the modelling of: 
 * - arithmetic operations over integers
 * - C-like pointers and arrays
 * - functions 
 * 
 * Important note: objects of the class cfg are not copyable. Instead,
 * we provide a class cfg_ref that wraps cfg references into copyable
 * and assignable objects.
 *
 * Limitations: 
 * - changes are needed to support floating point operations.
 * 
 */

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <crab/common/types.hpp>
#include <crab/common/bignums.hpp>
#include <crab/iterators/thresholds.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/discrete_domains.hpp>

#include <functional> // for wrapper_reference

namespace crab {

  namespace cfg_impl  {
     // To convert a basic block label to a string
     template< typename T >
     inline std::string get_label_str(T e);
  } 

  namespace cfg {

    using namespace ikos;
    using namespace std;

    // The values must be such that INT <= PTR <= ARR
    enum tracked_precision { INT = 0, PTR = 1, ARR = 2 };
  
    enum variable_type { INT_TYPE, PTR_TYPE, ARR_TYPE, UNK_TYPE};


    enum stmt_code {
      UNDEF = 0,
      // integers
      BIN_OP = 1, ASSIGN = 21, ASSUME = 22, UNREACH = 23, HAVOC = 24, SELECT = 25,
      ASSERT = 26,
      // arrays 
      ARR_INIT = 30, ARR_ASSUME = 31, ARR_STORE = 32, ARR_LOAD = 33,
      // pointers
      PTR_LOAD = 40, PTR_STORE = 41, PTR_ASSIGN = 42, 
      PTR_OBJECT = 43, PTR_FUNCTION = 44, PTR_NULL=45, 
      PTR_ASSUME = 46, PTR_ASSERT = 47,
      // functions calls
      CALLSITE = 50, RETURN = 51 
    }; 

    template<typename Number, typename VariableName>
    inline crab_os& operator<< (crab_os &o, 
                                const ikos::variable<Number, VariableName> &v)
    {
      auto tmp (v);
      tmp.write (o);
      return o;
    }

    template<typename Number, typename VariableName>
    inline crab_os& operator<< (crab_os &o, 
                                const ikos::linear_expression<Number, VariableName> &e)
    {
      auto tmp (e);
      tmp.write (o);
      return o;
    }

    template<typename Number, typename VariableName>
    inline crab_os& operator<< (crab_os &o, 
                                const ikos::linear_constraint<Number, VariableName> &c)
    {
      auto tmp (c);
      tmp.write (o);
      return o;
    }
  
    inline crab_os& operator<< (crab_os& o, variable_type t)
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
      typedef vector < VariableName > live_set_t;

     public:
      
      typedef typename live_set_t::const_iterator  const_use_iterator;
      typedef typename live_set_t::const_iterator  const_def_iterator;
      
     private:

      live_set_t m_uses;
      live_set_t m_defs;

      void add (live_set_t & s, VariableName v)
      {
        auto it = find (s.begin (), s.end (), v);
        if (it == s.end ()) s.push_back (v);
      }
      
     public:
      
      Live() { }
      
      void add_use(const VariableName v){ add (m_uses,v);}
      void add_def(const VariableName v){ add (m_defs,v);}
      
      const_use_iterator uses_begin() const { return m_uses.begin (); }
      const_use_iterator uses_end()   const { return m_uses.end (); }
      const_use_iterator defs_begin() const { return m_defs.begin (); }
      const_use_iterator defs_end()   const { return m_defs.end (); }

      friend crab_os& operator<<(crab_os &o, 
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

    struct debug_info {
      
      std::string m_file;
      int m_line;
      int m_col;
      
      debug_info ():
          m_file (""), m_line (-1), m_col (-1) { }

      debug_info (std::string file, unsigned line, unsigned col):
          m_file (file), m_line (line), m_col (col) { }
      
      bool operator<(const debug_info& other) const {
        return (m_file < other.m_file && 
                m_line < other.m_line && 
                m_col < other.m_col);
      }

      bool operator==(const debug_info& other) const {
        return (m_file == other.m_file && 
                m_line == other.m_line && 
                m_col == other.m_col);
      }

      bool has_debug  () const {
        return ((m_file != "") && (m_line >= 0) && (m_col >= 0));
      }

      void write (crab_os&o) const {
        o << "File  : " << m_file << "\n"
          << "Line  : " << m_line  << "\n" 
          << "Column: " << m_col << "\n";
      }
    };

    inline crab_os& operator<<(crab_os& o, const debug_info& l) {
      l.write (o);
      return o;
    }

    template< typename VariableName>
    struct statement_visitor;
  
    template< class VariableName>
    class statement
    {
      
     public:
      typedef Live<VariableName> live_t ;
      
     protected:
      live_t m_live;
      stmt_code m_stmt_code;
      debug_info m_dbg_info;

      statement (stmt_code code = UNDEF,
                 debug_info dbg_info = debug_info ()): 
          m_stmt_code (code),
          m_dbg_info (dbg_info) { }

     public:
      
      bool is_return () const { 
        return m_stmt_code == RETURN; 
      }
      bool is_bin_op () const { 
        return (m_stmt_code >= 1 && m_stmt_code <= 20); 
      }
      bool is_assign () const { 
        return (m_stmt_code == ASSIGN); 
      }
      bool is_assume () const { 
        return (m_stmt_code == ASSUME); 
      }
      bool is_select () const { 
        return (m_stmt_code == SELECT); 
      }
      bool is_assert () const { 
        return (m_stmt_code == ASSERT); 
      }
      bool is_arr_read () const { 
        return (m_stmt_code == ARR_LOAD);
      }
      bool is_ptr_read () const {
        return (m_stmt_code == PTR_LOAD); 
      }
      bool is_arr_write () const { 
        return (m_stmt_code == ARR_STORE); 
      }
      bool is_ptr_write () const { 
        return (m_stmt_code == PTR_STORE); 
      }
      bool is_ptr_null () const { 
        return (m_stmt_code == PTR_NULL); 
      }
      bool is_ptr_assume () const { 
        return (m_stmt_code == PTR_ASSUME); 
      }
      bool is_ptr_assert () const { 
        return (m_stmt_code == PTR_ASSERT); 
      }
      
     public:
      
      live_t get_live() const { return m_live; }

      debug_info get_debug_info () const { return m_dbg_info; }

      virtual void accept(statement_visitor< VariableName> *) = 0;
      
      virtual void write(crab_os& o) const = 0 ;
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const = 0;
      
      virtual ~statement() { }
      
      friend crab_os& operator <<(crab_os&o, 
                                  const statement<VariableName> &s)
      {
        s.write (o);
        return o;
      }
      
    }; 
  
    /*
      Basic statements (INT_TYPE)
    */

    template< class Number, class VariableName>
    class binary_op: public statement <VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_expression< Number, VariableName > linear_expression_t;
      
     private:
      
      variable_t          m_lhs;
      binary_operation_t  m_op;
      linear_expression_t m_op1;
      linear_expression_t m_op2;
      
     public:
      
      binary_op (variable_t lhs, 
                binary_operation_t op, 
                linear_expression_t op1, 
                linear_expression_t op2,
                debug_info dbg_info = debug_info ())
          : statement <VariableName> (BIN_OP, dbg_info),
            m_lhs(lhs), m_op(op), m_op1(op1), m_op2(op2) 
      { 
        this->m_live.add_def (m_lhs.name());
        for (auto v: m_op1.variables()){ this->m_live.add_use (v.name()); }         
        for (auto v: m_op2.variables()){ this->m_live.add_use (v.name()); }         
      }
      
      variable_t lhs () const { return m_lhs; }
      
      binary_operation_t op () const { return m_op; }
      
      linear_expression_t left () const { return m_op1; }
      
      linear_expression_t right () const { return m_op2; }
      
      virtual void accept(statement_visitor < VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef binary_op<Number, VariableName> binary_op_t;
        return boost::static_pointer_cast<statement <VariableName>, binary_op_t>
            (boost::make_shared<binary_op_t>(m_lhs, m_op, m_op1, m_op2));
      }
      
      virtual void write (crab_os& o) const
      {
        o << m_lhs << " = " << m_op1 << m_op << m_op2;
        return;
      }
    }; 

    template< class Number, class VariableName>
    class assignment: public statement<VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName >          variable_t;
      typedef linear_expression< Number, VariableName > linear_expression_t;
      
     private:
      
      variable_t          m_lhs;
      linear_expression_t m_rhs;
      
     public:
      
      assignment (variable_t lhs, linear_expression_t rhs)
          : statement <VariableName> (ASSIGN),
            m_lhs(lhs), m_rhs(rhs) 
      {
        this->m_live.add_def (m_lhs.name());
        for(auto v: m_rhs.variables()) 
          this->m_live.add_use (v.name());
      }
      
      variable_t lhs () const { return m_lhs; }
      
      linear_expression_t rhs () const { return m_rhs; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef assignment <Number, VariableName> assignment_t;
        return boost::static_pointer_cast< statement <VariableName>, assignment_t >
            (boost::make_shared<assignment_t>(m_lhs, m_rhs));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = " << m_rhs; // << " " << this->m_live;
        return;
      }
      
    }; 
    
    template<class Number, class VariableName>
    class assume_stmt: public statement <VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_constraint< Number, VariableName > linear_constraint_t;
      
     private:
      
      linear_constraint_t m_cst;
      
     public:
      
      assume_stmt (linear_constraint_t cst): 
          statement <VariableName> (ASSUME), m_cst(cst) 
      { 
        for(auto v: cst.variables())
          this->m_live.add_use (v.name()); 
      }
      
      linear_constraint_t constraint() const { return m_cst; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef assume_stmt <Number, VariableName> assume_t;
        return boost::static_pointer_cast< statement <VariableName>, assume_t >
            (boost::make_shared<assume_t>(m_cst));
      }
      
      virtual void write (crab_os & o) const
      {
        o << "assume (" << m_cst << ")"; //  << " " << this->m_live;
        return;
      }
    }; 

    template< class VariableName>
    class unreachable_stmt: public statement< VariableName> 
    {
     public:
      
      unreachable_stmt(): statement <VariableName> (UNREACH) { }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef unreachable_stmt <VariableName> unreachable_t;
        return boost::static_pointer_cast< statement <VariableName>, unreachable_t >
            (boost::make_shared<unreachable_t>());
      }
      
      virtual void write(crab_os& o) const
      {
        o << "unreachable";
        return;
      }
      
    }; 
  
    template< class VariableName>
    class havoc_stmt: public statement< VariableName> 
    {
      
      VariableName m_lhs;
      
     public:
      
      havoc_stmt (VariableName lhs): 
          statement <VariableName> (HAVOC), m_lhs(lhs)  {
        this->m_live.add_def (m_lhs);
      }
      
      VariableName variable () const { return m_lhs; }
      
      virtual void accept (statement_visitor<VariableName> *v) 
      {
        v->visit (*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef havoc_stmt <VariableName> havoc_t;
        return boost::static_pointer_cast< statement <VariableName>, havoc_t >
            (boost::make_shared<havoc_t>(m_lhs));
      }
      
      void write (crab_os& o) const
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
    class select_stmt: public statement <VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_expression< Number, VariableName > linear_expression_t;
      typedef linear_constraint< Number, VariableName > linear_constraint_t;
      
     private:
      
      variable_t          m_lhs;
      linear_constraint_t m_cond;
      linear_expression_t m_e1;
      linear_expression_t m_e2;
      
     public:
      
      select_stmt (variable_t lhs, 
              linear_constraint_t cond, 
              linear_expression_t e1, 
              linear_expression_t e2): 
          statement <VariableName> (SELECT),
          m_lhs(lhs), m_cond(cond), m_e1(e1), m_e2(e2) 
      { 
        this->m_live.add_def (m_lhs.name());
        for (auto v: m_cond.variables())
          this->m_live.add_use (v.name()); 
        for (auto v: m_e1.variables())
          this->m_live.add_use (v.name()); 
        for (auto v: m_e2.variables())
          this->m_live.add_use (v.name());
      }
      
      variable_t lhs () const { return m_lhs; }
      
      linear_constraint_t cond () const { return m_cond; }
      
      linear_expression_t left () const { return m_e1; }
      
      linear_expression_t right () const { return m_e2; }
      
      virtual void accept(statement_visitor < VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef select_stmt <Number, VariableName> select_t;
        return boost::static_pointer_cast< statement <VariableName>, select_t >
            (boost::make_shared<select_t>(m_lhs, m_cond, m_e1, m_e2));
      }
      
      virtual void write (crab_os& o) const
      {
        o << m_lhs << " = " 
          << "ite(" << m_cond << "," << m_e1 << "," << m_e2 << ")";
        return;
      }
    }; 

    template<class Number, class VariableName>
    class assert_stmt: public statement <VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_constraint< Number, VariableName > linear_constraint_t;
      
     private:
      
      linear_constraint_t m_cst;

     public:
      
      assert_stmt (linear_constraint_t cst, debug_info dbg_info = debug_info ())
          : statement <VariableName> (ASSERT, dbg_info), 
            m_cst(cst)
      { 
        for(auto v: cst.variables())
          this->m_live.add_use (v.name()); 
      }
      
      linear_constraint_t constraint() const { return m_cst; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef assert_stmt <Number, VariableName> assert_t;
        return boost::static_pointer_cast< statement <VariableName>, assert_t >
            (boost::make_shared<assert_t>(m_cst));
      }
      
      virtual void write (crab_os & o) const
      {
        o << "assert (" << m_cst << ")"; 
        return;
      }
    }; 

    /*
      Array statements (ARR_TYPE)
    */
  
    //! All the initial values of the array are statically known.
    template<class VariableName>
    class array_init_stmt: public statement< VariableName> 
    {
      VariableName m_arr; 
      vector<ikos::z_number> m_values; 
      
     public:
      array_init_stmt (VariableName arr, vector<ikos::z_number> values)
          : statement <VariableName> (ARR_INIT),
            m_arr (arr), m_values (values)  { }
      
      VariableName variable () const { return m_arr; }
      
      vector<ikos::z_number> values () const { return m_values; }
      
      virtual void accept (statement_visitor<VariableName> *v) {
        v->visit (*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef array_init_stmt <VariableName> array_init_t;
        return boost::static_pointer_cast< statement <VariableName>, array_init_t >
            (boost::make_shared<array_init_t>(m_arr, m_values));
      }
      
      void write (crab_os& o) const
      {
        o << m_arr << " = array_init ({";
        for (auto v: m_values)
          o << v << " ";
        o << "})";
        return;
      }
    }; 
  
    //! Assume all the array elements are overapproximated uniformly by
    //! some interval.
    template< class Number, class VariableName>
    class assume_array_stmt: public statement< VariableName> 
    {
      typedef bound <Number> bound_t;
      
     public:
      typedef interval <Number> interval_t;
      
     private:
      // forall i \in arr. m_val.lb () <= arr[i]  <= m_val.ub ()
      VariableName m_arr; 
      interval_t m_val;
      
     public:
      
      assume_array_stmt (VariableName arr, ikos::z_number val)
          : statement <VariableName> (ARR_ASSUME),
            m_arr (arr), m_val (bound_t (val))  { }
      
      assume_array_stmt (VariableName arr, interval_t val): 
          statement <VariableName> (ARR_ASSUME),
          m_arr (arr), m_val (val)  { }
      
      VariableName variable () const { return m_arr; }
      
      interval_t val () const { return m_val; }
      
      virtual void accept (statement_visitor<VariableName> *v) {
        v->visit (*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef assume_array_stmt <Number, VariableName> assume_array_t;
        return boost::static_pointer_cast< statement <VariableName>, assume_array_t >
            (boost::make_shared<assume_array_t>(m_arr, m_val));
      }
      
      void write (crab_os& o) const
      {
        o << "assume (forall l. " << m_arr << "[l]=";
        val ().write (o);
        o << ")";          
        return;
      }
    }; 
    
    template< class Number, class VariableName>
    class array_store_stmt: public statement<VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_expression< Number, VariableName > linear_expression_t;
      
     private:
      
      variable_t m_arr;
      linear_expression_t m_index;
      linear_expression_t m_value;
      ikos::z_number m_elem_size; // size in bytes
      
      bool m_is_singleton; //! whether the store writes to a singleton
                           //  cell. If unknown set to false.
     public:
      
      array_store_stmt (variable_t arr, 
                        linear_expression_t index, 
                        linear_expression_t value, 
                        ikos::z_number elem_size,
                        bool is_sing = false)
          : statement<VariableName>(ARR_STORE),
            m_arr (arr), m_index (index),
            m_value (value), m_elem_size (elem_size), 
            m_is_singleton (is_sing)
      {
        this->m_live.add_use (m_arr.name());
        for(auto v: m_index.variables()) 
          this->m_live.add_use (v.name());
        for(auto v: m_value.variables()) 
          this->m_live.add_use (v.name());
      }
      
      variable_t array () const { return m_arr; }
      
      linear_expression_t index () const { return m_index; }
      
      linear_expression_t value () const { return m_value; }
      
      ikos::z_number elem_size () const { return m_elem_size; }
      
      bool is_singleton () const { return m_is_singleton;}
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef array_store_stmt <Number, VariableName> array_store_t;
        return boost::static_pointer_cast< statement <VariableName>, array_store_t>
            (boost::make_shared<array_store_t>(m_arr, m_index, m_value, m_elem_size,
                                               m_is_singleton));
                                               
    }
      
      virtual void write(crab_os& o) const
      {
        o << "array_store(" << m_arr << "," << m_index << "," << m_value << ")"; 
        return;
      }
    }; 

    template< class Number, class VariableName>
    class array_load_stmt: public statement<VariableName>
    {
      
     public:
      
      typedef variable< Number, VariableName > variable_t;
      typedef linear_expression< Number, VariableName > linear_expression_t;
      
     private:
      
      variable_t m_lhs;
      variable_t m_array;
      linear_expression_t m_index;
      ikos::z_number m_elem_size; //! size in bytes
      
     public:
      
      array_load_stmt (variable_t lhs, variable_t arr, 
                 linear_expression_t index, ikos::z_number elem_size)
          : statement <VariableName> (ARR_LOAD),
            m_lhs (lhs), m_array (arr), 
            m_index (index), m_elem_size (elem_size)
      {
        this->m_live.add_def (lhs.name());
        this->m_live.add_use (m_array.name());
        for(auto v: m_index.variables()) 
          this->m_live.add_use (v.name());
      }
      
      variable_t lhs () const { return m_lhs; }
      
      variable_t array () const { return m_array; }
      
      linear_expression_t index () const { return m_index; }
      
      ikos::z_number elem_size () const { return m_elem_size; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef array_load_stmt <Number, VariableName> array_load_t;
        return boost::static_pointer_cast< statement <VariableName>, array_load_t>
            (boost::make_shared<array_load_t>(m_lhs, m_array, m_index, m_elem_size));
                                              
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = " 
          << "array_load(" << m_array << "," << m_index  << ")"; 
        return;
      }
    }; 
  
    /*
      Pointer statements (PTR_TYPE)
    */

    template<class VariableName>
    class ptr_load_stmt: public statement<VariableName>
    {
      // p = *q
      
      VariableName m_lhs;
      VariableName m_rhs;
      
     public:
      
      ptr_load_stmt (VariableName lhs, VariableName rhs, 
               debug_info dbg_info = debug_info ())
          : statement <VariableName> (PTR_LOAD, dbg_info),
            m_lhs (lhs), m_rhs (rhs)
      {
        this->m_live.add_use (lhs);
        this->m_live.add_use (rhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      VariableName rhs () const { return m_rhs; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_load_stmt <VariableName> ptr_load_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_load_t >
            (boost::make_shared<ptr_load_t>(m_lhs, m_rhs));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = "  << "*(" << m_rhs << ")";
        return;
      }
      
    }; 

    template<class VariableName>
    class ptr_store_stmt: public statement<VariableName>
    {
      // *p = q
      
      VariableName m_lhs;
      VariableName m_rhs;
      
     public:
      
      ptr_store_stmt (VariableName lhs, VariableName rhs,
                debug_info dbg_info = debug_info ())
          : statement <VariableName> (PTR_STORE, dbg_info),
            m_lhs (lhs), m_rhs (rhs)
      {
        this->m_live.add_use (lhs);
        this->m_live.add_use (rhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      VariableName rhs () const { return m_rhs; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_store_stmt <VariableName> ptr_store_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_store_t >
            (boost::make_shared<ptr_store_t>(m_lhs, m_rhs));
      }
      
      virtual void write(crab_os& o) const
      {
        o << "*(" << m_lhs << ") = " << m_rhs;
        return;
      }
    }; 

    template< class Number, class VariableName>
    class ptr_assign_stmt: public statement<VariableName>
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
      
      ptr_assign_stmt (VariableName lhs, VariableName rhs, linear_expression_t offset)
          : statement <VariableName> (PTR_ASSIGN),
            m_lhs (lhs), m_rhs (rhs), m_offset(offset)
      {
        this->m_live.add_def (lhs);
        this->m_live.add_use (rhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      VariableName rhs () const { return m_rhs; }
      
      linear_expression_t offset () const { return m_offset; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_assign_stmt <Number, VariableName> ptr_assign_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_assign_t >
            (boost::make_shared<ptr_assign_t>(m_lhs, m_rhs, m_offset));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = "  << m_rhs << " + ";
        linear_expression_t off (m_offset) ; // FIX: write method is not const in ikos 
        off.write (o);
        return;
      }
    }; 
  
    template<class VariableName>
    class ptr_object_stmt: public statement<VariableName>
    {
      //! lhs = &a;
      VariableName m_lhs;
      index_t m_address;
      
     public:
      
      ptr_object_stmt (VariableName lhs, index_t address)
          : statement <VariableName> (PTR_OBJECT),
            m_lhs (lhs), m_address (address)
      {
        this->m_live.add_def (lhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      index_t rhs () const { return m_address; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_object_stmt <VariableName> ptr_object_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_object_t>
            (boost::make_shared<ptr_object_t>(m_lhs, m_address));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = "  << "&(" << m_address << ")" ;
        return;
      }
    }; 

    template<class VariableName>
    class ptr_function_stmt: public statement<VariableName>
    {
      // lhs = &func;
      VariableName m_lhs;
      VariableName m_func; // Pre: function names are unique 
      
     public:
      
      ptr_function_stmt (VariableName lhs, VariableName func)
          : statement <VariableName> (PTR_FUNCTION),
            m_lhs (lhs), m_func (func)
      {
        this->m_live.add_def (lhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      VariableName rhs () const { return m_func; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_function_stmt <VariableName> ptr_function_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_function_t>
            (boost::make_shared<ptr_function_t>(lhs (), rhs ()));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = "  << "&(" << m_func << ")" ;
        return;
      }
    }; 

    template<class VariableName>
    class ptr_null_stmt: public statement<VariableName>
    {
      //! lhs := null;
      VariableName m_lhs;

     public:
      
      ptr_null_stmt (VariableName lhs)
          : statement <VariableName> (PTR_NULL), m_lhs (lhs) 
      {
        this->m_live.add_def (m_lhs);
      }
      
      VariableName lhs () const { return m_lhs; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_null_stmt <VariableName> ptr_null_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_null_t>
            (boost::make_shared<ptr_null_t>(m_lhs));
      }
      
      virtual void write(crab_os& o) const
      {
        o << m_lhs << " = "  << "NULL";
        return;
      }
    }; 

    template<class VariableName>
    class ptr_assume_stmt: public statement<VariableName>
    {
     public:

      typedef pointer_constraint<VariableName> ptr_cst_t;

     private:

      ptr_cst_t m_cst;
      
     public:
      
      ptr_assume_stmt (ptr_cst_t cst)
          : statement <VariableName> (PTR_ASSUME),
            m_cst (cst) 
      {
        if (!cst.is_tautology () && !cst.is_contradiction ()) {
          if (cst.is_unary ()) {
            this->m_live.add_use (cst.lhs ());
          } else {
            this->m_live.add_use (cst.lhs ());
            this->m_live.add_use (cst.rhs ());
          }
        }
      }
      
      ptr_cst_t constraint () const { return m_cst; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_assume_stmt <VariableName> ptr_assume_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_assume_t>
            (boost::make_shared<ptr_assume_t>(m_cst));
      }
      
      virtual void write(crab_os& o) const
      {
        o << "assume_ptr(" << m_cst << ")";
        return;
      }
    }; 

    template<class VariableName>
    class ptr_assert_stmt: public statement<VariableName>
    {
     public:

      typedef pointer_constraint<VariableName> ptr_cst_t;

     private:

      ptr_cst_t m_cst;
      
     public:
      
      ptr_assert_stmt (ptr_cst_t cst, debug_info dbg_info = debug_info ())
          : statement <VariableName> (PTR_ASSERT, dbg_info),
            m_cst (cst) 
      {
        if (!cst.is_tautology () && !cst.is_contradiction ()) {
          if (cst.is_unary ()) {
            this->m_live.add_use (cst.lhs ());
          } else {
            this->m_live.add_use (cst.lhs ());
            this->m_live.add_use (cst.rhs ());
          }
        }
      }
      
      ptr_cst_t constraint () const { return m_cst; }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef ptr_assert_stmt <VariableName> ptr_assert_t;
        return boost::static_pointer_cast< statement <VariableName>, ptr_assert_t>
            (boost::make_shared<ptr_assert_t>(m_cst));
      }
      
      virtual void write(crab_os& o) const
      {
        o << "assert_ptr(" << m_cst << ")";
        return;
      }
    }; 
  
  
    /*
      Function calls
    */
  
    template< class VariableName>
    class callsite_stmt: public statement<VariableName>
    {
      
      boost::optional<pair<VariableName,variable_type> > m_lhs;
      VariableName m_func_name;
      vector<pair<VariableName,variable_type> > m_args;
      
      typedef typename vector<pair<VariableName,variable_type> >::iterator arg_iterator;
      typedef typename vector<pair<VariableName,variable_type> >::const_iterator const_arg_iterator;
      
     public:
      
      callsite_stmt (VariableName func_name, 
                     vector<pair<VariableName,variable_type> > args)
          : m_func_name (func_name)
      {
        std::copy (args.begin (), args.end (), std::back_inserter (m_args));
        for (auto arg:  m_args) { this->m_live.add_use (arg.first); }
      }
      
      callsite_stmt (VariableName func_name, vector<VariableName> args)
          : statement <VariableName> (CALLSITE),
            m_func_name (func_name)
      {
        for (auto v : args)
        {
          m_args.push_back (make_pair (v, UNK_TYPE));
          this->m_live.add_use (v);
        }
      }
      
      callsite_stmt (pair<VariableName,variable_type> lhs, 
                     VariableName func_name, 
                     vector<pair<VariableName,variable_type> > args)
          : m_lhs (boost::optional<pair<VariableName,variable_type> > (lhs)), 
            m_func_name (func_name)
      {
        std::copy (args.begin (), args.end (), std::back_inserter(m_args));
        for (auto arg:  m_args) { this->m_live.add_use (arg.first); }
        this->m_live.add_def ((*m_lhs).first);
      }
      
      callsite_stmt (VariableName lhs, 
                VariableName func_name, vector<VariableName> args)
          : m_lhs (boost::optional<pair<VariableName,variable_type> > (lhs, UNK_TYPE)), 
            m_func_name (func_name)
      {
        for (auto v : args)
        {
          m_args.push_back (make_pair (v, UNK_TYPE));
          this->m_live.add_use (v);
        }
        this->m_live.add_def ((*m_lhs).first);
      }
            
      boost::optional<VariableName> get_lhs_name () const { 
        if (m_lhs)
          return boost::optional<VariableName> ((*m_lhs).first);
        else 
        return boost::optional<VariableName> ();
      }
      
      variable_type get_lhs_type () const {       
        if (m_lhs) return (*m_lhs).second;
        else return UNK_TYPE;
      }
    
      VariableName get_func_name () const { 
        return m_func_name; 
      }
      
      unsigned get_num_args () const { return m_args.size (); }

      VariableName get_arg_name (unsigned idx) const { 
        if (idx >= m_args.size ())
          CRAB_ERROR ("Out-of-bound access to call site parameter");
        
        return m_args[idx].first;
      }
      
      variable_type get_arg_type (unsigned idx) const { 
        if (idx >= m_args.size ())
        CRAB_ERROR ("Out-of-bound access to call site parameter");
        
        return m_args[idx].second;
      }
      
      virtual void accept(statement_visitor <VariableName> *v) {
        v->visit(*this);
      }
      
      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef callsite_stmt <VariableName> call_site_t;
        if (m_lhs) {
          return boost::static_pointer_cast< statement <VariableName>, 
                                             call_site_t>
              (boost::make_shared<call_site_t> (*m_lhs, m_func_name, m_args));
        }
        else {
          return boost::static_pointer_cast< statement <VariableName>, 
                                             call_site_t>
              (boost::make_shared<call_site_t>(m_func_name, m_args));
        }
      }

      virtual void write(crab_os& o) const
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
    class return_stmt: public statement<VariableName>
    {
      
      VariableName m_var;
      variable_type m_type;
      
     public:
      
      return_stmt (VariableName var, variable_type type = UNK_TYPE)
          : statement <VariableName> (RETURN),
            m_var (var), m_type (type)
      {
        this->m_live.add_use (m_var); 
      }
      
      VariableName get_ret_var () const { 
        return m_var;
      }
      
      variable_type get_ret_type () const { 
        return m_type;
      }
      
      virtual void accept(statement_visitor <VariableName> *v) 
      {
        v->visit(*this);
      }

      virtual boost::shared_ptr<statement <VariableName> > clone () const
      {
        typedef return_stmt <VariableName> return_t;
        return boost::static_pointer_cast< statement <VariableName>, return_t>
            (boost::make_shared<return_t>(m_var, m_type));
      }
      
      virtual void write(crab_os& o) const
      {
        o << "return " << m_var;
        return;
      }
    }; 
  
    template< class BasicBlockLabel, class VariableName>
    class Cfg;
  
    template< class BasicBlockLabel, class VariableName >
    class basic_block: public boost::noncopyable
    {
      
      // TODO: support for removing statements 
      friend class Cfg< BasicBlockLabel, VariableName >;
      
     public:
      
      typedef VariableName varname_t;
      typedef BasicBlockLabel basic_block_label_t;

      // helper types to build statements
      typedef variable< z_number, VariableName > z_variable_t;
      typedef linear_expression< z_number, VariableName > z_lin_exp_t;
      typedef linear_constraint< z_number, VariableName > z_lin_cst_t;
      typedef statement< VariableName> statement_t;
      typedef basic_block< BasicBlockLabel, VariableName > basic_block_t;
      typedef interval <z_number> z_interval;
      
     private:
      
      typedef vector< BasicBlockLabel > bb_id_set_t;
      typedef boost::shared_ptr< statement_t > statement_ptr;
      typedef vector< statement_ptr > stmt_list_t;
      
     public:
      
      typedef typename bb_id_set_t::iterator succ_iterator;
      typedef typename bb_id_set_t::const_iterator const_succ_iterator;
      typedef succ_iterator pred_iterator;
      typedef const_succ_iterator const_pred_iterator;
      typedef boost::indirect_iterator< typename stmt_list_t::iterator > iterator;
      typedef boost::indirect_iterator< typename stmt_list_t::const_iterator > const_iterator;
      typedef boost::indirect_iterator< typename stmt_list_t::reverse_iterator > reverse_iterator;
      typedef boost::indirect_iterator< typename stmt_list_t::const_reverse_iterator > const_reverse_iterator;
      typedef discrete_domain <VariableName> live_domain_t;
      
     public:

      // Basic statements
      typedef binary_op<z_number,VariableName> z_bin_op_t;
      typedef assignment<z_number,VariableName> z_assign_t;
      typedef assume_stmt<z_number,VariableName> z_assume_t;
      typedef havoc_stmt<VariableName> havoc_t;
      typedef unreachable_stmt<VariableName> unreach_t;
      typedef select_stmt<z_number,VariableName> z_select_t;
      typedef assert_stmt<z_number,VariableName> z_assert_t;
      // Functions
      typedef callsite_stmt<VariableName> callsite_t;
      typedef return_stmt<VariableName> return_t;
      // Arrays
      typedef array_init_stmt<VariableName> z_arr_init_t;
      typedef assume_array_stmt<z_number,VariableName> z_assume_arr_t;
      typedef array_store_stmt<z_number,VariableName> z_arr_store_t;
      typedef array_load_stmt<z_number,VariableName> z_arr_load_t;
      // Pointers
      typedef ptr_store_stmt<VariableName> ptr_store_t;
      typedef ptr_load_stmt<VariableName> ptr_load_t;
      typedef ptr_assign_stmt<z_number,VariableName> ptr_assign_t;
      typedef ptr_object_stmt<VariableName> ptr_object_t;
      typedef ptr_function_stmt<VariableName> ptr_function_t;
      typedef ptr_null_stmt<VariableName> ptr_null_t;
      typedef ptr_assume_stmt<VariableName> ptr_assume_t;
      typedef ptr_assert_stmt<VariableName> ptr_assert_t;

     private:

      typedef boost::shared_ptr<z_bin_op_t> z_bin_op_ptr;
      typedef boost::shared_ptr<z_assign_t> z_assign_ptr;
      typedef boost::shared_ptr<z_assume_t> z_assume_ptr;
      typedef boost::shared_ptr<havoc_t> havoc_ptr;      
      typedef boost::shared_ptr<unreach_t> unreach_ptr;
      typedef boost::shared_ptr<z_select_t> z_select_ptr;
      typedef boost::shared_ptr<z_assert_t> z_assert_ptr;
      typedef boost::shared_ptr<callsite_t> callsite_ptr;      
      typedef boost::shared_ptr<return_t> return_ptr;      
      typedef boost::shared_ptr<z_arr_init_t> z_arr_init_ptr;
      typedef boost::shared_ptr<z_assume_arr_t> z_assume_arr_ptr;
      typedef boost::shared_ptr<z_arr_store_t> z_arr_store_ptr;
      typedef boost::shared_ptr<z_arr_load_t> z_arr_load_ptr;    
      typedef boost::shared_ptr<ptr_store_t> ptr_store_ptr;
      typedef boost::shared_ptr<ptr_load_t> ptr_load_ptr;    
      typedef boost::shared_ptr<ptr_assign_t> ptr_assign_ptr;
      typedef boost::shared_ptr<ptr_object_t> ptr_object_ptr;    
      typedef boost::shared_ptr<ptr_function_t> ptr_function_ptr;    
      typedef boost::shared_ptr<ptr_null_t> ptr_null_ptr;    
      typedef boost::shared_ptr<ptr_assume_t> ptr_assume_ptr;    
      typedef boost::shared_ptr<ptr_assert_t> ptr_assert_ptr;    

      
      BasicBlockLabel m_bb_id;
      stmt_list_t m_stmts;
      bb_id_set_t m_prev, m_next;
      tracked_precision m_track_prec;    
      // Ideally it should be size_t to indicate any position within the
      // block. For now, we only allow to insert either at front or at
      // the back (default). Note that if insertions at the front are
      // very common we should replace stmt_list_t from a vector to a
      // deque.
      bool m_insert_point_at_front; 

      // set of used/def variables 
      live_domain_t m_live; 
      
      void insert_adjacent (bb_id_set_t &c, BasicBlockLabel e)
      { 
        if (std::find(c.begin (), c.end (), e) == c.end ())
          c.push_back (e);
      }
      
      void remove_adjacent (bb_id_set_t &c, BasicBlockLabel e)
      {
        if (std::find(c.begin (), c.end (), e) != c.end ())
          c.erase (std::remove(c.begin (), c.end (), e), c.end ());
      }
      
      basic_block (BasicBlockLabel bb_id, tracked_precision track_prec): 
          m_bb_id (bb_id), m_track_prec (track_prec), 
          m_insert_point_at_front (false), 
          m_live (live_domain_t::bottom ())
      { }
      
      static boost::shared_ptr< basic_block_t > 
      create (BasicBlockLabel bb_id, tracked_precision track_prec) 
      {
        return boost::shared_ptr<basic_block_t>(new basic_block_t(bb_id, track_prec));
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

        auto ls = stmt->get_live ();
        for (auto &v : boost::make_iterator_range (ls.uses_begin (), ls.uses_end ()))
          m_live += v;
        for (auto &v : boost::make_iterator_range (ls.defs_begin (), ls.defs_end ()))
          m_live += v;
      }
      
     public:
      
      //! it will be set to false after the first insertion
      void set_insert_point_front (){
        m_insert_point_at_front = true;
      }
      
      boost::shared_ptr<basic_block_t> clone () const
      {
        boost::shared_ptr<basic_block_t> b (new basic_block_t(label (), m_track_prec));
        
        for (auto &s : boost::make_iterator_range (begin (), end ()))
          b->m_stmts.push_back (s.clone ()); 
        
        for (auto id : boost::make_iterator_range (prev_blocks ())) 
          b->m_prev.push_back (id);
        
        for (auto id : boost::make_iterator_range (next_blocks ()))
          b->m_next.push_back (id);

        b->m_live = m_live;
        return b;
      }

      BasicBlockLabel label () const { return m_bb_id; }

      string name () const {
        return cfg_impl::get_label_str (m_bb_id); 
      }

      iterator begin() { 
        return boost::make_indirect_iterator (m_stmts.begin ()); 
      }
      iterator end() { 
        return boost::make_indirect_iterator (m_stmts.end ()); 
      }
      const_iterator begin() const { 
        return boost::make_indirect_iterator (m_stmts.begin ()); 
      }
      const_iterator end() const {
        return boost::make_indirect_iterator (m_stmts.end ()); 
      }

      reverse_iterator rbegin() {            
        return boost::make_indirect_iterator (m_stmts.rbegin ()); 
      }
      reverse_iterator rend() {              
        return boost::make_indirect_iterator (m_stmts.rend ()); 
      }
      const_reverse_iterator rbegin() const {
        return boost::make_indirect_iterator (m_stmts.rbegin ()); 
      }
      const_reverse_iterator rend() const {
        return boost::make_indirect_iterator (m_stmts.rend ()); 
      }
      
      size_t size() { return std::distance ( begin (), end ()); }

      live_domain_t& live () {
        return m_live;
      }

      live_domain_t live () const {
        return m_live;
      }

      pair<succ_iterator, succ_iterator> next_blocks ()
      { 
        return make_pair (m_next.begin (), m_next.end ());
      }
      
      pair<pred_iterator, pred_iterator> prev_blocks() 
      { 
        return make_pair (m_prev.begin (), m_prev.end ());
      }
      
      pair<const_succ_iterator,const_succ_iterator> next_blocks () const
      { 
        return make_pair (m_next.begin (), m_next.end ());
      }
      
      pair<const_pred_iterator,const_pred_iterator> prev_blocks() const
      { 
        return make_pair (m_prev.begin (), m_prev.end ());
      }

      // void reverse()
      // {
      //   std::swap (m_prev, m_next);
      //   std::reverse (m_stmts.begin (), m_stmts.end ());
      // }
      
      // Add a cfg edge from *this to b
      void operator>>(basic_block_t& b) 
      {
        insert_adjacent (m_next, b.m_bb_id);
        insert_adjacent (b.m_prev, m_bb_id);
      }
      
      // Remove a cfg edge from *this to b
      void operator-=(basic_block_t &b)
      {
        remove_adjacent (m_next, b.m_bb_id);
        remove_adjacent (b.m_prev, m_bb_id);       
      }
      
      // insert all statements of other at the front
      void merge_front (const basic_block_t &other) 
      {
        m_stmts.insert (m_stmts.begin (), 
                        other.m_stmts.begin (), 
                        other.m_stmts.end ());

        m_live = m_live | other.m_live;
      }
      
      // insert all statements of other at the back
      void merge_back (const basic_block_t &other) 
      {
        m_stmts.insert (m_stmts.end (), 
                        other.m_stmts.begin (), 
                        other.m_stmts.end ());

        m_live = m_live | other.m_live;
      }
      
      void write(crab_os& o) const
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
               (boost::make_shared<z_bin_op_t>(lhs, BINOP_ADD, op1, op2)));
      }
      
      void add (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert(boost::static_pointer_cast< statement_t, z_bin_op_t > 
               (boost::make_shared<z_bin_op_t>(lhs, BINOP_ADD, op1,  op2)));
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
                 (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("add operands unexpected");
      }
      
      void sub (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SUB, op1, op2)));
      }
      
      void sub (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SUB, op1, op2)));
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
                  (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("sub operands unexpected");
      }
      
      void mul (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert(boost::static_pointer_cast< statement_t, z_bin_op_t >
               (boost::make_shared<z_bin_op_t>(lhs, BINOP_MUL, op1, op2)));
      }
      
      void mul (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_MUL, op1, op2)));
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
                 (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("mul operands unexpected");
      }

      // signed division
      void div (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SDIV, op1, op2)));
      }
      
      void div (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SDIV, op1, op2)));
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
                  (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("div operands unexpected");
      }

      // unsigned division
      void udiv (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_UDIV, op1, op2)));
      }
      
      void udiv (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_UDIV, op1, op2)));
      }
      
      void udiv (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          udiv (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          udiv (lhs, (*op1.get_variable ()), op2.constant ());
        
        else
          CRAB_ERROR("udiv operands unexpected");
      }

      // signed rem
      void rem (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SREM, op1, op2)));
      }
      
      void rem (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_SREM, op1, op2)));
      }
      
      void rem (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          rem (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          rem (lhs, (*op1.get_variable ()), op2.constant ());
        
        else if (op1.is_constant () && op2.is_constant ()) 
        {
          z_lin_exp_t rhs (z_number (op1.constant() % op2.constant()));
          insert (boost::static_pointer_cast< statement_t, z_assign_t >
                  (boost::make_shared<z_assign_t> (lhs, rhs)));
        }
        else
          CRAB_ERROR("rem operands unexpected");
      }

      // unsigned rem
      void urem (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_UREM, op1, op2)));
      }
      
      void urem (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_UREM, op1, op2)));
      }
      
      void urem (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          urem (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          urem (lhs, (*op1.get_variable ()), op2.constant ());
        
        else
          CRAB_ERROR("urem operands unexpected");
      }

      
      void bitwise_and (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_AND, op1, op2)));
      }
      
      void bitwise_and (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t> (lhs, BINOP_AND, op1, op2)));
      }
      
      void bitwise_and (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          bitwise_and (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          bitwise_and (lhs, (*op1.get_variable ()), op2.constant ());        
        
        else if (op1.is_constant () && op2.get_variable ())
          bitwise_and (lhs, (*op2.get_variable ()), op1.constant ());

        else if (op1.is_constant () && op2.is_constant ()) 
        {
          z_lin_exp_t rhs (z_number (op1.constant () & op2.constant ()));
          insert (boost::static_pointer_cast< statement_t, z_assign_t >
                  (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("and operands unexpected");
      }

      void bitwise_or (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t> (lhs, BINOP_OR, op1, op2)));
      }
      
      void bitwise_or (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_OR, op1, op2)));
      }
      
      void bitwise_or (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          bitwise_or (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          bitwise_or (lhs, (*op1.get_variable ()), op2.constant ());        

        else if (op1.is_constant () && op2.get_variable ())
          bitwise_or (lhs, (*op2.get_variable ()), op1.constant ());
        
        else if (op1.is_constant () && op2.is_constant ()) 
        {
          z_lin_exp_t rhs (z_number (op1.constant () | op2.constant ()));
          insert (boost::static_pointer_cast< statement_t, z_assign_t >
                  (boost::make_shared<z_assign_t>(lhs, rhs)));
        }
        else
          CRAB_ERROR("or operands unexpected");
      }

      void bitwise_xor (z_variable_t lhs, z_variable_t op1, z_variable_t op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_XOR, op1, op2)));
      }
      
      void bitwise_xor (z_variable_t lhs, z_variable_t op1, z_number op2) 
      {
        insert (boost::static_pointer_cast< statement_t, z_bin_op_t >
                (boost::make_shared<z_bin_op_t>(lhs, BINOP_XOR, op1, op2)));
      }
      
      void bitwise_xor (z_variable_t lhs, z_lin_exp_t op1, z_lin_exp_t op2) 
      {
        if (op1.get_variable () && op2.get_variable ())
          bitwise_xor (lhs, (*op1.get_variable ()), (*op2.get_variable ()));
        
        else if (op1.get_variable () && op2.is_constant ())
          bitwise_xor (lhs, (*op1.get_variable ()), op2.constant ());        

        else if (op1.is_constant () && op2.get_variable ())
          bitwise_xor (lhs, (*op2.get_variable ()), op1.constant ());
        
        else if (op1.is_constant () && op2.is_constant ()) 
        {
          z_lin_exp_t rhs (z_number (op1.constant () ^ op2.constant ()));
          insert (boost::static_pointer_cast< statement_t, z_assign_t >
                  (boost::make_shared<z_assign_t> (lhs, rhs)));
        }
        else
          CRAB_ERROR("xor operands unexpected");
      }
      
      void assign (z_variable_t lhs, z_lin_exp_t rhs) 
      {
        insert (boost::static_pointer_cast< statement_t, z_assign_t >
                (boost::make_shared<z_assign_t> (lhs, rhs)));
      }
      
      void assume (z_lin_cst_t cst) 
      {
        insert (boost::static_pointer_cast< statement_t, z_assume_t >
                (boost::make_shared<z_assume_t> (cst)));
      }
      
      void havoc(VariableName lhs) 
      {
        insert (boost::static_pointer_cast< statement_t, havoc_t > 
                (boost::make_shared<havoc_t> (lhs)));
      }
      
      void unreachable() 
      {
        insert (boost::static_pointer_cast< statement_t, unreach_t > 
                (boost::make_shared<unreach_t> ()));
      }
      
      void select (z_variable_t lhs, z_variable_t v, z_lin_exp_t e1, z_lin_exp_t e2) 
      {
        z_lin_cst_t cond = (v >= z_number(1));
        insert(boost::static_pointer_cast< statement_t, z_select_t >
               (boost::make_shared<z_select_t>(lhs, cond, e1, e2)));
      }
      
      void select (z_variable_t lhs, z_lin_cst_t cond, z_lin_exp_t e1, z_lin_exp_t e2) 
      {
        insert(boost::static_pointer_cast< statement_t, z_select_t >
               (boost::make_shared<z_select_t>(lhs, cond, e1, e2)));
      }

      void assertion (z_lin_cst_t cst) 
      {
        insert (boost::static_pointer_cast< statement_t, z_assert_t >
                (boost::make_shared<z_assert_t> (cst)));
      }

      void assertion (z_lin_cst_t cst, debug_info di) 
      {
        insert (boost::static_pointer_cast< statement_t, z_assert_t >
                (boost::make_shared<z_assert_t> (cst, di)));
      }
      
      void callsite (VariableName func, 
                     vector<pair <VariableName,variable_type> > args) 
      {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (boost::make_shared<callsite_t>(func, args)));
      }
      
      void callsite (VariableName func, vector<VariableName> args) 
      {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (boost::make_shared<callsite_t>(func, args)));
      }
      
      void callsite (pair<VariableName,variable_type> lhs, 
                     VariableName func, 
                     vector<pair <VariableName,variable_type> > args) 
      {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (boost::make_shared<callsite_t>(lhs, func, args)));
      }
      
      void callsite (VariableName lhs, VariableName func, vector<VariableName> args) 
      {
        insert(boost::static_pointer_cast< statement_t, callsite_t >
               (boost::make_shared<callsite_t>(lhs, func, args)));
      }
      
      void ret (VariableName var, variable_type ty) 
      {
        insert(boost::static_pointer_cast< statement_t, return_t >
               (boost::make_shared<return_t>(var, ty)));
      }
      
      void ret (VariableName var) 
      {
        insert(boost::static_pointer_cast< statement_t, return_t >
               (boost::make_shared<return_t>(var)));
      }
      
      void array_init (VariableName a, 
                       const vector<ikos::z_number>& vals) {
        insert (boost::static_pointer_cast< statement_t, z_arr_init_t > 
                (boost::make_shared<z_arr_init_t> (a, vals)));
      }
      
      void assume_array (VariableName a, z_interval val) {
        insert (boost::static_pointer_cast< statement_t, z_assume_arr_t > 
                (boost::make_shared<z_assume_arr_t> (a, val)));
      }
      
      void assume_array (VariableName a, ikos::z_number val) {
        insert (boost::static_pointer_cast< statement_t, z_assume_arr_t > 
                (boost::make_shared<z_assume_arr_t> (a, val)));
      }
      
      void array_store (z_variable_t arr, z_lin_exp_t idx, 
                        z_lin_exp_t val, ikos::z_number elem_size, 
                        bool is_singleton = false) 
      {
        if (m_track_prec == ARR)
          insert(boost::static_pointer_cast< statement_t, z_arr_store_t >
                 (boost::make_shared<z_arr_store_t>(arr, idx, val, elem_size, is_singleton)));
      }
      
      void array_load (z_variable_t lhs, z_variable_t arr, 
                       z_lin_exp_t idx, ikos::z_number elem_size) 
      {
        if (m_track_prec == ARR)
          insert(boost::static_pointer_cast< statement_t, z_arr_load_t >
                 (boost::make_shared<z_arr_load_t>(lhs, arr, idx, elem_size)));
      }
      
      void ptr_store (VariableName lhs, VariableName rhs) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_store_t >
                 (boost::make_shared<ptr_store_t> (lhs, rhs)));
      }
      
      void ptr_load (VariableName lhs, VariableName rhs) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_load_t >
                 (boost::make_shared<ptr_load_t> (lhs, rhs)));
      }
      
      void ptr_assign (VariableName lhs, VariableName rhs, z_lin_exp_t offset) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_assign_t >
                 (boost::make_shared<ptr_assign_t> (lhs, rhs, offset)));
      }
      
      void new_object (VariableName lhs, index_t address) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_object_t >
                 (boost::make_shared<ptr_object_t> (lhs, address)));
      }
      
      void new_ptr_func (VariableName lhs, index_t func) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_function_t >
                 (boost::make_shared<ptr_function_t> (lhs, func)));
      }

      void ptr_null (VariableName lhs) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_null_t >
                 (boost::make_shared<ptr_null_t> (lhs)));
      }

      void ptr_assume (pointer_constraint<VariableName> cst) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_assume_t >
                 (boost::make_shared<ptr_assume_t> (cst)));
      }

      void ptr_assertion (pointer_constraint<VariableName> cst) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_assert_t >
                 (boost::make_shared<ptr_assert_t> (cst)));
      }

      void ptr_assertion (pointer_constraint<VariableName> cst, debug_info di) 
      {
        if (m_track_prec >= PTR)
          insert(boost::static_pointer_cast< statement_t, ptr_assert_t >
                 (boost::make_shared<ptr_assert_t> (cst, di)));
      }
      
      friend crab_os& operator<<(crab_os &o, const basic_block_t &b)
      {
        //b.write (o);
        o << cfg_impl::get_label_str (b);
        return o;
      }
      
    }; 

    // Viewing a BasicBlock with all statements reversed. Useful for
    // backward analysis.
    template<class BasicBlock> 
    class basic_block_rev {
     public:
      typedef typename BasicBlock::varname_t varname_t;
      typedef typename BasicBlock::basic_block_label_t basic_block_label_t;

      typedef basic_block_rev<BasicBlock> basic_block_rev_t;

      typedef typename BasicBlock::succ_iterator succ_iterator;
      typedef typename BasicBlock::const_succ_iterator const_succ_iterator;
      typedef succ_iterator pred_iterator;
      typedef const_succ_iterator const_pred_iterator;

      typedef typename BasicBlock::reverse_iterator iterator;
      typedef typename BasicBlock::const_reverse_iterator const_iterator;
      typedef discrete_domain<varname_t> live_domain_t;

     private:

      BasicBlock& _bb;

     public:

      basic_block_rev(BasicBlock& bb): _bb(bb) {}

      basic_block_label_t label () const { return _bb.label(); }

      std::string name () const { return _bb.name();}

      iterator begin() { return _bb.rbegin();}            

      iterator end() { return _bb.rend();}             

      const_iterator begin() const { return _bb.rbegin(); }

      const_iterator end() const { return _bb.rend();}
      
      std::size_t size() const { return std::distance ( begin (), end ()); }

      live_domain_t& live () { return _bb.live(); }

      live_domain_t live () const { return _bb.live(); }

      pair<succ_iterator, succ_iterator> next_blocks()
      { return _bb.prev_blocks(); }
      
      pair<pred_iterator, pred_iterator> prev_blocks() 
      { return _bb.next_blocks(); }
      
      pair<const_succ_iterator,const_succ_iterator> next_blocks () const
      { return _bb.prev_blocks(); }
      
      pair<const_pred_iterator,const_pred_iterator> prev_blocks() const
      { return _bb.next_blocks(); }

      void write(crab_os& o) const
      {
        o << name() << ":\n";	
        for (auto const &s: *this)
        { o << "  " << s << ";\n"; }
        o << "--> [";
        for (auto const &n : boost::make_iterator_range (next_blocks ()))
        { o << n << ";"; }
        o << "]\n";
        return;
      }     

      friend crab_os& operator<<(crab_os &o, const basic_block_rev_t &b) {
        //b.write (o);
        o << b.name();
        return o;
      }
    };
  
    template< class VariableName>
    struct statement_visitor
    {
      typedef linear_expression< z_number, VariableName > z_lin_exp_t;
      typedef binary_op <z_number,VariableName> z_bin_op_t;
      typedef assignment <z_number,VariableName> z_assign_t;
      typedef assume_stmt <z_number,VariableName> z_assume_t;
      typedef havoc_stmt<VariableName> havoc_t;
      typedef unreachable_stmt<VariableName> unreach_t;
      typedef select_stmt<z_number,VariableName> z_select_t;
      typedef assert_stmt <z_number,VariableName> z_assert_t;
      
      typedef callsite_stmt<VariableName> callsite_t;
      typedef return_stmt<VariableName> return_t;
      
      typedef array_init_stmt<VariableName> z_arr_init_t;
      typedef assume_array_stmt<z_number,VariableName> z_assume_arr_t;
      typedef array_store_stmt<z_number,VariableName> z_arr_store_t;
      typedef array_load_stmt<z_number,VariableName> z_arr_load_t;
      
      typedef ptr_store_stmt<VariableName> ptr_store_t;
      typedef ptr_load_stmt<VariableName> ptr_load_t;
      typedef ptr_assign_stmt<z_number,VariableName> ptr_assign_t;
      typedef ptr_object_stmt<VariableName> ptr_object_t;
      typedef ptr_function_stmt<VariableName> ptr_function_t;
      typedef ptr_null_stmt<VariableName> ptr_null_t;
      typedef ptr_assume_stmt<VariableName> ptr_assume_t;
      typedef ptr_assert_stmt<VariableName> ptr_assert_t;
      
      // Only implementation for basic statements is required
      
      virtual void visit (z_bin_op_t&) = 0;
      virtual void visit (z_assign_t&) = 0;
      virtual void visit (z_assume_t&) = 0;
      virtual void visit (havoc_t&) = 0;
      virtual void visit (unreach_t&) = 0;
      virtual void visit (z_select_t&) = 0;
      virtual void visit (z_assert_t&) { };
      
      virtual void visit (callsite_t&) { };
      virtual void visit (return_t&) { };
      virtual void visit (z_arr_init_t&) { };
      virtual void visit (z_assume_arr_t&) { };
      virtual void visit (z_arr_store_t&) { };
      virtual void visit (z_arr_load_t&) { };
      virtual void visit (ptr_store_t&) { };
      virtual void visit (ptr_load_t&) { };
      virtual void visit (ptr_assign_t&) { };
      virtual void visit (ptr_object_t&) { };
      virtual void visit (ptr_function_t&) { };
      virtual void visit (ptr_null_t&) { };
      virtual void visit (ptr_assume_t&) { };
      virtual void visit (ptr_assert_t&) { };
      
      virtual ~statement_visitor () { }
    }; 
    
    template< class VariableName>
    class function_decl
    {
      
      variable_type m_lhs_type;
      VariableName m_func_name;
      vector<pair<VariableName,variable_type> > m_params;
      
      typedef typename vector<pair<VariableName,variable_type> >::iterator param_iterator;
      typedef typename vector<pair<VariableName,variable_type> >::const_iterator const_param_iterator;
      
     public:
      
      function_decl (VariableName func_name, vector<VariableName> params)
          : m_lhs_type (UNK_TYPE), m_func_name (func_name)
      {
        for (auto v : params)
          m_params.push_back (make_pair (v, UNK_TYPE));
      }
      
      function_decl (variable_type lhs_type, VariableName func_name, 
                     vector<pair<VariableName,variable_type> > params)
          : m_lhs_type (lhs_type), m_func_name (func_name)
      {
        std::copy (params.begin (), params.end (), 
                   std::back_inserter (m_params));
      }
      
      variable_type get_lhs_type () const { return m_lhs_type; }
      
      VariableName get_func_name () const { return m_func_name;  }
      
      unsigned get_num_params () const { return m_params.size (); }
      
      VariableName get_param_name (unsigned idx) const { 
        if (idx >= m_params.size ())
          CRAB_ERROR ("Out-of-bound access to function parameter");
        
        return m_params[idx].first;
      }
      
      variable_type get_param_type (unsigned idx) const { 
        if (idx >= m_params.size ())
          CRAB_ERROR ("Out-of-bound access to function parameter");
        
        return m_params[idx].second;
      }

      void write(crab_os& o) const
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
      
      friend crab_os& operator<<(crab_os& o, const function_decl<VariableName> &decl)
      { 
        decl.write (o);
        return o;
      }
    }; 

    // forward declarations
    template<class Any> class cfg_rev;
    template<class Any> class cfg_ref;
     
    template< class BasicBlockLabel, class VariableName >
    class Cfg: public boost::noncopyable {
     public:

      typedef BasicBlockLabel basic_block_label_t;
      typedef basic_block_label_t node_t; // for Bgl graphs
      typedef VariableName varname_t;
      typedef function_decl<varname_t> fdecl_t;
      typedef basic_block<BasicBlockLabel, VariableName > basic_block_t;   
      typedef statement<VariableName > statement_t;
      typedef crab::iterators::thresholds<z_number> thresholds_t;

      typedef typename basic_block_t::succ_iterator succ_iterator;
      typedef typename basic_block_t::pred_iterator pred_iterator;
      typedef typename basic_block_t::const_succ_iterator const_succ_iterator;
      typedef typename basic_block_t::const_pred_iterator const_pred_iterator;
      
      typedef boost::iterator_range <succ_iterator> succ_range;
      typedef boost::iterator_range <pred_iterator> pred_range;
      typedef boost::iterator_range <const_succ_iterator> const_succ_range;
      typedef boost::iterator_range <const_pred_iterator> const_pred_range;
      
     private:
      
      typedef Cfg<BasicBlockLabel, VariableName > cfg_t;
      typedef boost::shared_ptr< basic_block_t > basic_block_ptr;
      typedef boost::unordered_map< BasicBlockLabel, basic_block_ptr > basic_block_map_t;
      typedef typename basic_block_map_t::value_type binding_t;
      typedef typename basic_block_t::live_domain_t live_domain_t;

      struct get_ref : public std::unary_function<binding_t, basic_block_t>
      {
        get_ref () { }
        basic_block_t& operator () (const binding_t &p) const { return *(p.second); }
      }; 
      
      struct getLabel : public std::unary_function<binding_t, BasicBlockLabel>
      {
        getLabel () { }
        BasicBlockLabel operator () (const binding_t &p) const { return p.second->label(); }
      }; 
      
     public:
      
      typedef boost::transform_iterator< get_ref, 
                                         typename basic_block_map_t::iterator > iterator;
      typedef boost::transform_iterator< get_ref, 
                                         typename basic_block_map_t::const_iterator > const_iterator;
      typedef boost::transform_iterator< getLabel, 
                                         typename basic_block_map_t::iterator > label_iterator;
      typedef boost::transform_iterator< getLabel, 
                                         typename basic_block_map_t::const_iterator > const_label_iterator;

      typedef typename std::vector<varname_t>::iterator var_iterator;
      typedef typename std::vector<varname_t>::const_iterator const_var_iterator;

     private:
      
      BasicBlockLabel m_entry;
      BasicBlockLabel m_exit;
      bool m_has_exit;
      basic_block_map_t m_blocks;
      tracked_precision m_track_prec;
      //! we allow to define a cfg without being associated with a
      //! function
      boost::optional<fdecl_t> m_func_decl; 
      
      
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
      
      struct print_block 
      {
        crab_os &m_o;
        print_block (crab_os& o) : m_o (o) { }
        void operator () (const basic_block_t& B){ B.write (m_o); }
      };
      
      
     public:
      
      // --- needed by crab::cg::call_graph<CFG>::cg_node
      Cfg () { }
      
      Cfg (BasicBlockLabel entry, tracked_precision track_prec = INT)
          : m_entry  (entry),  
            m_has_exit (false),
            m_track_prec (track_prec) {
        m_blocks.insert (binding_t (m_entry, 
                                    basic_block_t::create (m_entry, m_track_prec)));
      }
      
      Cfg (BasicBlockLabel entry, BasicBlockLabel exit, 
           tracked_precision track_prec = INT)
          : m_entry  (entry), 
            m_exit   (exit), 
            m_has_exit (true),
            m_track_prec (track_prec) {
        m_blocks.insert (binding_t (m_entry, 
                                    basic_block_t::create (m_entry, m_track_prec)));
      }
      
      Cfg (BasicBlockLabel entry, BasicBlockLabel exit, 
           fdecl_t func_decl, 
           tracked_precision track_prec = INT)
          : m_entry (entry), 
            m_exit (exit), 
            m_has_exit (true),
            m_track_prec (track_prec),
            m_func_decl (boost::optional<fdecl_t> (func_decl)) {
        m_blocks.insert (binding_t (m_entry, 
                                    basic_block_t::create (m_entry, m_track_prec)));
      }
      
      boost::shared_ptr<cfg_t> clone () const
      {
        boost::shared_ptr<cfg_t> _cfg (new cfg_t (m_entry, m_track_prec));
        _cfg->m_has_exit = m_has_exit ;
        if (_cfg->m_has_exit)
          _cfg->m_exit = m_exit ;
        _cfg->m_func_decl = m_func_decl;
        for (auto const &BB: boost::make_iterator_range (begin (), end ()))
        {
          boost::shared_ptr <basic_block_t> copyBB = BB.clone ();
          _cfg->m_blocks.insert (binding_t (copyBB->label (), copyBB));
        }
        return _cfg;
      }
      
      boost::optional<fdecl_t> get_func_decl () const { 
        return m_func_decl; 
      }
      
      tracked_precision get_track_prec () const {
        return m_track_prec;
      }
      
      
      bool has_exit () const { return m_has_exit; }
      
      BasicBlockLabel exit()  const { 
        if (has_exit ()) return m_exit; 
        
        CRAB_ERROR ("cfg does not have an exit block");
      } 
      
      //! set method to mark the exit block after the cfg has been
      //! created.
      void set_exit (BasicBlockLabel exit) { 
        m_exit = exit; 
        m_has_exit = true;
      }
      
      //! set method to add the function declaration after the cfg has
      //! been created.
      void set_func_decl (fdecl_t decl) { 
        m_func_decl  = boost::optional<fdecl_t> (decl);
      }

      // --- Begin ikos fixpoint API

      BasicBlockLabel entry() const { return m_entry; } 

      const_succ_range next_nodes (BasicBlockLabel bb_id) const
      {
        const basic_block_t& b = get_node(bb_id);
        return boost::make_iterator_range (b.next_blocks ());
      }
      
      const_pred_range prev_nodes (BasicBlockLabel bb_id) const
      {
        const basic_block_t& b = get_node(bb_id);
        return boost::make_iterator_range (b.prev_blocks ());
      }
      
      succ_range next_nodes (BasicBlockLabel bb_id) 
      {
        basic_block_t& b = get_node(bb_id);
        return boost::make_iterator_range (b.next_blocks ());
      }
      
      pred_range prev_nodes (BasicBlockLabel bb_id) 
      {
        basic_block_t& b = get_node(bb_id);
        return boost::make_iterator_range (b.prev_blocks ());
      }

      basic_block_t& get_node (BasicBlockLabel bb_id) 
      {
        auto it = m_blocks.find (bb_id);
        if (it == m_blocks.end ())
          CRAB_ERROR ("Basic block not found in the CFG");
        
        return *(it->second);
      }
      
      const basic_block_t& get_node (BasicBlockLabel bb_id) const
      {
        auto it = m_blocks.find (bb_id);
        if (it == m_blocks.end ())
          CRAB_ERROR ("Basic block not found in the CFG");
        
        return *(it->second);
      }

      // --- End ikos fixpoint API

      basic_block_t& insert (BasicBlockLabel bb_id) 
      {
        auto it = m_blocks.find (bb_id);
        if (it != m_blocks.end ()) return *(it->second);
        
        basic_block_ptr block = basic_block_t::create (bb_id, m_track_prec);
        m_blocks.insert (binding_t (bb_id, block));
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
        
        m_blocks.erase (bb_id);
      }
      
      // Return all variables (either used or defined) in the cfg.
      //
      // This operation is linear on the size of the cfg to still keep
      // a valid set in case a block is removed.
      std::vector<varname_t> get_vars () const {
        live_domain_t ls = live_domain_t::bottom ();
        for (auto const &b : boost::make_iterator_range (begin (), end ()))
          ls = ls | b.live ();
        // std::vector<varname_t> vars (ls.size ());
        // vars.insert (vars.end (), ls.begin (), ls.end ());
        std::vector<varname_t> vars;
        for (auto v: ls) vars.push_back (v);
        return vars;
      }
            
      //! return a begin iterator of BasicBlock's
      iterator begin() 
      {
        return boost::make_transform_iterator (m_blocks.begin (), get_ref ());
      }
      
      //! return an end iterator of BasicBlock's
      iterator end() 
      {
        return boost::make_transform_iterator (m_blocks.end (), get_ref ());
      }
      
      const_iterator begin() const
      {
        return boost::make_transform_iterator (m_blocks.begin (), get_ref ());
      }
      
      const_iterator end() const
      {
        return boost::make_transform_iterator (m_blocks.end (), get_ref ());
      }
      
      //! return a begin iterator of BasicBlockLabel's
      label_iterator label_begin() 
      {
        return boost::make_transform_iterator (m_blocks.begin (), getLabel ());
      }
      
      //! return an end iterator of BasicBlockLabel's
      label_iterator label_end() 
      {
        return boost::make_transform_iterator (m_blocks.end (), getLabel ());
      }
      
      const_label_iterator label_begin() const
      {
        return boost::make_transform_iterator (m_blocks.begin (), getLabel ());
      }
      
      const_label_iterator label_end() const
      {
        return boost::make_transform_iterator (m_blocks.end (), getLabel ());
      }
      
      size_t size () const { return std::distance (begin (), end ()); }
     
      thresholds_t initialize_thresholds_for_widening (size_t size) const {
        typedef typename thresholds_t::bound_t bound_t;

        thresholds_t thresholds (size);
        for (auto const &b : boost::make_iterator_range (begin (), end ())) {
          for (auto const&i : boost::make_iterator_range (b.begin (), b.end ())) {

            bound_t t = bound_t::plus_infinity ();
            if (i.is_assume ()) {
              auto assume_inst = static_cast<const typename basic_block_t::z_assume_t*> (&i);
              t = -(assume_inst->constraint ().expression ().constant ());
            }
            else if (i.is_select ()) {
              auto select_inst = static_cast<const typename basic_block_t::z_select_t*> (&i);
              t = -(select_inst->cond ().expression ().constant ());
            }
            
            if (t != bound_t::plus_infinity ()) {
              ////
              // For code pattern like this "if(x<k1) {x+=k2;}" note
              // that the condition (x<k1) is translated to (x<=k1-1)
              // so an useful threshold would be k1+1+k2.  Since we
              // don't keep track of how x is incremented or
              // decremented we choose arbitrarily k2=1.
              ////
              thresholds.add (t+2);
            }
          }
        }
        
        return thresholds;
      }

      // void reverse()
      // {
      //   if (!m_has_exit)
      //     CRAB_ERROR ("cfg cannot be reversed: no exit block found");
        
      //   std::swap (m_entry, m_exit);
      //   for (auto &p: *this) 
      //     p.reverse (); 
      // }
      
      void write (crab_os& o) const
      {
        print_block f (o);
        if (m_func_decl)
          o << *m_func_decl << "\n";
        dfs (f);
        return;
      }
      
      friend crab_os& operator<<(crab_os &o, 
                                 const Cfg< BasicBlockLabel, VariableName > &cfg)
      {
        cfg.write (o);
        return o;
      }
      
      void simplify ()
      {
        merge_blocks ();        
        remove_unreachable_blocks ();
        remove_useless_blocks ();
        //after removing useless blocks there can be opportunities to
        //merge more blocks.
        merge_blocks ();
        merge_blocks ();
      }
      
     private:
      
      ////
      // cfg simplifications
      ////
      
      //XXX: this is a bit adhoc. It should be probably a parameter of
      // simplify().
      struct donot_simplify_visitor: public statement_visitor<VariableName>
      {
        typedef typename statement_visitor<VariableName>::z_bin_op_t z_bin_op_t;
        typedef typename statement_visitor<VariableName>::z_assign_t z_assign_t;
        typedef typename statement_visitor<VariableName>::z_assume_t z_assume_t;
        typedef typename statement_visitor<VariableName>::havoc_t havoc_t;
        typedef typename statement_visitor<VariableName>::unreach_t unreach_t;
        typedef typename statement_visitor<VariableName>::z_select_t z_select_t;
        typedef typename statement_visitor<VariableName>::z_arr_load_t z_arr_load_t;
        
        bool _do_not_simplify;
        donot_simplify_visitor (): _do_not_simplify(false) { }
        void visit(z_bin_op_t&){ }  
        void visit(z_assign_t&) { }
        void visit(z_assume_t&) { _do_not_simplify = true; }
        void visit(z_arr_load_t&) { _do_not_simplify = true; }
        void visit(havoc_t&) { }
        void visit(unreach_t&){ }
        void visit(z_select_t&){ }
      };
      
      // Helpers
      bool has_one_child (BasicBlockLabel b)
      {
        auto rng = next_nodes (b);
        return (std::distance (rng.begin (), rng.end ()) == 1);
      }
      
      bool has_one_parent (BasicBlockLabel b)
      {
        auto rng = prev_nodes (b);
        return (std::distance (rng.begin (), rng.end ()) == 1);
      }
      
      basic_block_t& get_child (BasicBlockLabel b)
      {
        assert (has_one_child (b));
        auto rng = next_nodes (b);
        return get_node (*(rng.begin ()));
      }
      
      basic_block_t& get_parent (BasicBlockLabel b)
      {
        assert (has_one_parent (b));
        auto rng = prev_nodes (b);
        return get_node (*(rng.begin ()));
      }
      
      void merge_blocks_rec (BasicBlockLabel curId, 
                           visited_t& visited)
      {
        
        if (visited.find (curId) != visited.end ()) return;
        visited.insert (curId);
        
        basic_block_t &cur = get_node (curId);
        
        if (has_one_child (curId) && has_one_parent (curId))
        {
          basic_block_t &parent = get_parent (curId);
          basic_block_t &child  = get_child (curId);
          
          donot_simplify_visitor vis;
          for (auto it = cur.begin (); it != cur.end (); ++it)
            it->accept(&vis);
          
          if (!vis._do_not_simplify)
          {
            parent.merge_back (cur);
            remove (curId);
            parent >> child;        
            merge_blocks_rec (child.label (), visited); 
            return;
          }
        }
        
        for (auto n : boost::make_iterator_range (cur.next_blocks ()))
          merge_blocks_rec (n, visited);
      }
      
      // Merges a basic block into its predecessor if there is only one
      // and the predecessor only has one successor.
      void merge_blocks ()
      {
        visited_t visited;
        merge_blocks_rec (entry (), visited);
      }
      
      // mark reachable blocks from curId
      template<class AnyCfg>
      void mark_alive_blocks (BasicBlockLabel curId, 
                            AnyCfg& cfg,
                            visited_t& visited)
      {
        if (visited.count (curId) > 0) return;
        visited.insert (curId);
        for (auto child : cfg.next_nodes (curId))
          mark_alive_blocks (child, cfg, visited);
      }
      
      // remove unreachable blocks
      void remove_unreachable_blocks ()
      {
        visited_t alive, dead;
        mark_alive_blocks (entry (), *this, alive);
        
        for (auto const &bb : *this) 
          if (!(alive.count (bb.label ()) > 0))
            dead.insert (bb.label ());
        
        for (auto bb_id: dead)
          remove (bb_id);
      }
      
      // remove blocks that cannot reach the exit block
      void remove_useless_blocks ()
      {
        if (!has_exit ()) return;
        
        cfg_rev<cfg_ref<cfg_t> > rev_cfg (*this); 

        visited_t useful, useless;
        mark_alive_blocks (rev_cfg.entry (), rev_cfg, useful);
        
        for (auto const &bb : *this) 
          if (!(useful.count (bb.label ()) > 0))
            useless.insert (bb.label ());
        
        for (auto bb_id: useless)
          remove (bb_id);
      }
      
    }; 

    // A lightweight object that wraps a reference to a CFG into a
    // copyable, assignable object.
    template <class CFG>
    class cfg_ref {
     public:

      // CFG's typedefs
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::node_t node_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename CFG::fdecl_t fdecl_t;
      typedef typename CFG::basic_block_t basic_block_t;   
      typedef typename CFG::statement_t statement_t;
      typedef typename CFG::thresholds_t thresholds_t;

      typedef typename CFG::succ_iterator succ_iterator;
      typedef typename CFG::pred_iterator pred_iterator;
      typedef typename CFG::const_succ_iterator const_succ_iterator;
      typedef typename CFG::const_pred_iterator const_pred_iterator;
      typedef typename CFG::succ_range succ_range;
      typedef typename CFG::pred_range pred_range;
      typedef typename CFG::const_succ_range const_succ_range;
      typedef typename CFG::const_pred_range const_pred_range;
      typedef typename CFG::iterator iterator;
      typedef typename CFG::const_iterator const_iterator;
      typedef typename CFG::label_iterator label_iterator;
      typedef typename CFG::const_label_iterator const_label_iterator;
      typedef typename CFG::var_iterator var_iterator;
      typedef typename CFG::const_var_iterator const_var_iterator;

     private:

      boost::optional<reference_wrapper<CFG> > _ref;

     public:

      // --- hook needed by crab::cg::CallGraph<CFG>::CgNode
      cfg_ref () { } 

      cfg_ref (CFG &cfg)
          : _ref(reference_wrapper<CFG>(cfg)) { } 
      
      const CFG& get() const { 
        assert (_ref);
        return *_ref;
      }

      CFG& get() { 
        assert (_ref);
        return *_ref;
      }

      basic_block_label_t  entry() const {
        assert (_ref);
        return (*_ref).get().entry();
      }

      const_succ_range next_nodes (basic_block_label_t bb) const {
        assert (_ref);
        return (*_ref).get().next_nodes(bb);
      }

      const_pred_range prev_nodes (basic_block_label_t bb) const {
        assert (_ref);
        return (*_ref).get().prev_nodes(bb);
      }

      succ_range next_nodes (basic_block_label_t bb) {
        assert (_ref);
        return (*_ref).get().next_nodes(bb);
      }

      pred_range prev_nodes (basic_block_label_t bb) {
        assert (_ref);
        return (*_ref).get().prev_nodes(bb);
      }

      thresholds_t initialize_thresholds_for_widening (size_t size) const {
        assert (_ref);
        return (*_ref).get().initialize_thresholds_for_widening (size);
      }

      basic_block_t& get_node (basic_block_label_t bb) {
        assert (_ref);
        return (*_ref).get().get_node(bb);
      }
      
      const basic_block_t& get_node (basic_block_label_t bb) const {
        assert (_ref);
        return (*_ref).get().get_node(bb);
      }
 
      iterator begin() {
        assert (_ref);
        return (*_ref).get().begin();
      }
      
      iterator end() {
        assert (_ref);
        return (*_ref).get().end();
      }
      
      const_iterator begin() const {
        assert (_ref);
        return (*_ref).get().begin();
      }
      
      const_iterator end() const {
        assert (_ref);
        return (*_ref).get().end();
      }

      label_iterator label_begin() {
        assert (_ref);
        return (*_ref).get().label_begin();
      }
      
      label_iterator label_end() {
        assert (_ref);
        return (*_ref).get().label_end();
      }
      
      const_label_iterator label_begin() const {
        assert (_ref);
        return (*_ref).get().label_begin();
      }
      
      const_label_iterator label_end() const {
        assert (_ref);
        return (*_ref).get().label_end();
      }

      boost::optional<fdecl_t> get_func_decl () const { 
        assert (_ref);
        return (*_ref).get().get_func_decl();
      }

      
      bool has_exit () const {
        assert (_ref);
        return (*_ref).get().has_exit();
      }
      
      basic_block_label_t exit()  const { 
        assert (_ref);
        return (*_ref).get().exit();
      }

      
      friend crab_os& operator<<(crab_os &o, const cfg_ref<CFG> &cfg) {
        o << cfg.get();
        return o;
      }
      
      void simplify () {
        if (_ref) (*_ref).get().simplify();
      }

      // #include <boost/fusion/functional/invocation/invoke.hpp>
      // template< class... ArgTypes >
      // typename std::result_of<CFG&(ArgTypes&&...)>::type
      // operator() ( ArgTypes&&... args ) const {
      //   return boost::fusion::invoke(get(), std::forward<ArgTypes>(args)...);
      // }      
    };
  
    // Viewing a CFG with all edges and block statements
    // reversed. Useful for backward analysis.
    template<class CFGRef> // CFGRef must be copyable!
    class cfg_rev {
     public:
      typedef typename CFGRef::basic_block_label_t basic_block_label_t;
      typedef basic_block_rev<typename CFGRef::basic_block_t> basic_block_t;
      typedef basic_block_label_t node_t; // for Bgl graphs
      typedef typename CFGRef::varname_t varname_t;
      typedef typename CFGRef::fdecl_t fdecl_t;
      typedef typename CFGRef::statement_t statement_t;
      typedef typename CFGRef::thresholds_t thresholds_t;

      typedef typename CFGRef::succ_range pred_range;
      typedef typename CFGRef::pred_range succ_range;
      typedef typename CFGRef::const_succ_range const_pred_range;
      typedef typename CFGRef::const_pred_range const_succ_range;

      // For BGL
      typedef typename basic_block_t::succ_iterator succ_iterator;
      typedef typename basic_block_t::pred_iterator pred_iterator;
      typedef typename basic_block_t::const_succ_iterator const_succ_iterator;
      typedef typename basic_block_t::const_pred_iterator const_pred_iterator;

      typedef cfg_rev<CFGRef> cfg_rev_t;

     private:

      struct getRev : public std::unary_function<typename CFGRef::basic_block_t, basic_block_t> {
        const boost::unordered_map<basic_block_label_t, basic_block_t>& _rev_bbs;

        getRev(const boost::unordered_map<basic_block_label_t, basic_block_t>& rev_bbs)
            : _rev_bbs(rev_bbs) { }

        const basic_block_t& operator()(typename CFGRef::basic_block_t &bb) const {
          auto it = _rev_bbs.find(bb.label());
          if (it != _rev_bbs.end())
            return it->second;
          CRAB_ERROR ("Basic block not found in the CFG");          
        }
      }; 

     public:

      typedef boost::transform_iterator<getRev, typename CFGRef::iterator> iterator;
      typedef boost::transform_iterator<getRev, typename CFGRef::const_iterator> const_iterator;
      typedef typename CFGRef::label_iterator label_iterator;
      typedef typename CFGRef::const_label_iterator const_label_iterator;
      typedef typename CFGRef::var_iterator var_iterator;
      typedef typename CFGRef::const_var_iterator const_var_iterator;

     private:

      CFGRef _cfg;
      boost::unordered_map<basic_block_label_t, basic_block_t> _rev_bbs;
     
     public:

      // --- hook needed by crab::cg::CallGraph<CFGRef>::CgNode
      cfg_rev () { }

      cfg_rev (CFGRef cfg): _cfg(cfg) { 
        // Create basic_block_rev from BasicBlock objects
        // Note that basic_block_rev is also a view of BasicBlock so it
        // doesn't modify BasicBlock objects.
        for(auto &bb: cfg) {
          basic_block_t rev(bb);
          _rev_bbs.insert(std::make_pair(bb.label(), rev));
        }
      }

      cfg_rev(const cfg_rev_t& o)
          : _cfg(o._cfg), _rev_bbs(o._rev_bbs) { }

      cfg_rev(cfg_rev_t && o)
          : _cfg(std::move(o._cfg)), _rev_bbs(std::move(o._rev_bbs)) { }

      cfg_rev_t& operator=(const cfg_rev_t&o) {
        if (this != &o) {
          _cfg = o._cfg;
          _rev_bbs = o._rev_bbs;
        }
        return *this;
      }

      cfg_rev_t& operator=(cfg_rev_t&&o) {
        _cfg = std::move(o._cfg);
        _rev_bbs = std::move(o._rev_bbs);
        return *this;
      }

      basic_block_label_t  entry() const {
        if (!_cfg.has_exit()) CRAB_ERROR("Entry not found!");
        return _cfg.exit();
      }

      const_succ_range next_nodes (basic_block_label_t bb) const {
        return _cfg.prev_nodes(bb);
      }

      const_pred_range prev_nodes (basic_block_label_t bb) const {
        return _cfg.next_nodes(bb);
      }

      succ_range next_nodes (basic_block_label_t bb) {
        return _cfg.prev_nodes(bb);
      }

      pred_range prev_nodes (basic_block_label_t bb) {
        return _cfg.next_nodes(bb);
      }

      basic_block_t& get_node (basic_block_label_t bb_id) {
        auto it = _rev_bbs.find (bb_id);
        if (it == _rev_bbs.end()) 
          CRAB_ERROR ("Basic block not found in the CFG");
        return it->second;
      }

      const basic_block_t& get_node (basic_block_label_t bb_id) const {
        auto it = _rev_bbs.find (bb_id);
        if (it == _rev_bbs.end()) 
          CRAB_ERROR ("Basic block not found in the CFG");
        return it->second;
      }
      
      thresholds_t initialize_thresholds_for_widening (size_t size) const {
        return _cfg.initialize_thresholds_for_widening (size);
      }

      iterator begin() {
        return boost::make_transform_iterator (_cfg.begin(), getRev(_rev_bbs));
      }
      
      iterator end() {
        return boost::make_transform_iterator (_cfg.end(), getRev(_rev_bbs));
      }
      
      const_iterator begin() const {
        return boost::make_transform_iterator (_cfg.begin(), getRev(_rev_bbs));
      }
      
      const_iterator end() const {
        return boost::make_transform_iterator (_cfg.end(), getRev(_rev_bbs));
      }

      label_iterator label_begin() {
        return _cfg.label_begin();
      }
      
      label_iterator label_end() {
        return _cfg.label_end();
      }
      
      const_label_iterator label_begin() const {
        return _cfg.label_begin();
      }
      
      const_label_iterator label_end() const {
        return _cfg.label_end();
      }

      boost::optional<fdecl_t> get_func_decl () const { 
        return _cfg.get_func_decl();
      }
      
      bool has_exit () const {
        return true;
      }
      
      basic_block_label_t exit()  const { 
        return _cfg.entry();
      }

      void write(crab_os& o) const {
        if (get_func_decl())
          o << *get_func_decl() << "\n";
        for (auto &bb: *this)
        { bb.write(o); }
      }

      friend crab_os& operator<<(crab_os &o, const cfg_rev<CFGRef> &cfg) {
        cfg.write(o);
        return o;
      }
      
      void simplify () { }
      
    };

     // Helper class
    template<typename CFG>
    struct cfg_hasher {
      typedef typename CFG::basic_block_t::callsite_t callsite_t;
      typedef typename CFG::fdecl_t fdecl_t;
      
      static size_t hash (callsite_t cs) {
        size_t res = hash_value (cs.get_func_name ());
        boost::hash_combine (res, cs.get_lhs_type());
        auto num_args = cs.get_num_args ();
        for(unsigned i=0; i<num_args; i++)
          boost::hash_combine(res, cs.get_arg_type (i));
        return res;
      }
      
      static size_t hash (fdecl_t d)  {
        size_t res = hash_value (d.get_func_name ());
        boost::hash_combine (res, d.get_lhs_type ());
        for(unsigned i=0; i<d.get_num_params (); i++)
          boost::hash_combine(res, d.get_param_type (i));
        return res;
      }      
    };

    // extending boost::hash for cfg class
    template<class B, class V>
    std::size_t hash_value(Cfg<B,V> const& _cfg) {
      auto fdecl = _cfg.get_func_decl ();            
      if (!fdecl)
        CRAB_ERROR ("cannot hash a cfg because function declaration is missing");

      return cfg_hasher<Cfg<B,V> >::hash(*fdecl);
    }

    template<class CFG>
    std::size_t hash_value(cfg_ref<CFG> const& _cfg) {
      auto fdecl = _cfg.get().get_func_decl ();            
      if (!fdecl)
        CRAB_ERROR ("cannot hash a cfg because function declaration is missing");

      return cfg_hasher<cfg_ref<CFG> >::hash(*fdecl);
    }

    template<class CFGRef>
    std::size_t hash_value(cfg_rev<CFGRef> const& _cfg) {
      auto fdecl = _cfg.get().get_func_decl ();            
      if (!fdecl)
        CRAB_ERROR ("cannot hash a cfg because function declaration is missing");

      return cfg_hasher<cfg_rev<CFGRef> >::hash(*fdecl);
    }

    template<class B, class V>
    bool operator==(Cfg<B,V> const& a, Cfg<B,V> const& b) {
      return hash_value (a) == hash_value (b);
    }

    template<class CFG>
    bool operator==(cfg_ref<CFG> const& a, cfg_ref<CFG> const& b) {
      return hash_value (a) == hash_value (b);
    }

    template<class CFGRef>
    bool operator==(cfg_rev<CFGRef> const& a, cfg_rev<CFGRef> const& b) {
      return hash_value (a) == hash_value (b);
    }

  } // end namespace cfg
}  // end namespace crab

#endif /* CFG_HPP */

