#ifndef CFG_HPP
#define CFG_HPP

/* 
 * Build a CFG to interface with IKOS
 */


#include <ikos_domains/common.hpp>
#include <ikos_domains/bignums.hpp>
#include <ikos_domains/linear_constraints.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/noncopyable.hpp>
#include <boost/unordered_set.hpp>

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

    friend std::ostream& operator<<(std::ostream &o, const Live< VariableName> &live )
    {
      o << "Use={"; 
      for (auto const& v: boost::make_iterator_range (live.uses_begin (), live.uses_end ()))
        o << v << ",";
      o << "} Def={"; 
      for (auto const& v: boost::make_iterator_range (live.defs_begin (), live.defs_end ()))
        o << v << ",";
      o << "}";
      return o;
    }

  };


  template< typename VariableName>
  class StatementVisitor;

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

    virtual ostream& write(ostream& o) const = 0 ;

    virtual boost::shared_ptr<Statement <VariableName> > clone () const = 0;

    virtual ~Statement() { }

    friend ostream& operator <<(ostream&o, const Statement<VariableName> &s)
    {
      s.write (o);
      return o;
    }

  }; 

  template< class Number, class VariableName>
  class BinaryOp: public Statement <VariableName>
  {
    
   public:

    typedef variable< Number, VariableName >          variable_t;
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

    variable_t lhs () { return m_lhs; }
    
    operation_t op () { return m_op; }
    
    linear_expression_t left () { return m_op1; }
    
    linear_expression_t right () { return m_op2; }

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

    virtual ostream& write (ostream& o) const
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
      return o;
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
    
    variable_t lhs () { return m_lhs; }
    
    linear_expression_t rhs () { return m_rhs; }

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
    
    virtual ostream& write(ostream& o) const
    {
      o << m_lhs << " = " << m_rhs; // << " " << this->m_live;
      return o;
    }

  }; 
    
  template<class Number, class VariableName>
  class Assume: public Statement <VariableName>
  {

   public:

    typedef variable< Number, VariableName >          variable_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
      
   private:

    linear_constraint_t m_cst;
    
   public:

    Assume (linear_constraint_t cst): m_cst(cst) 
    { 
      for(auto v: cst.variables())
        this->m_live.addUse (v.name()); 
    }
    
    linear_constraint_t constraint() { return m_cst; }

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
        
    virtual ostream& write (ostream & o) const
    {
      o << "assume (" << m_cst << ")"; //  << " " << this->m_live;
      return o;
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
    
    virtual ostream& write(ostream& o) const
    {
      o << "unreachable";
      return o;
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
     
    VariableName variable () { return m_lhs; }
     
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
     
    ostream& write (ostream& o) const
    {
      o << m_lhs << " =*" << " "; // << this->m_live;
      return o;
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

    // To build statements
    typedef z_number                                    ZNumber;
    typedef variable< ZNumber, VariableName >           ZVariable;
    typedef linear_expression< ZNumber, VariableName >  ZLinearExpression;
    typedef linear_constraint< ZNumber, VariableName >  ZLinearConstraint;

    typedef Statement< VariableName>                    Statement_t;
    typedef BasicBlock< BasicBlockLabel, VariableName > BasicBlock_t;
    
   private:

    typedef std::vector< BasicBlockLabel >    bb_id_set_t;
    typedef boost::shared_ptr< bb_id_set_t >  bb_id_set_ptr;
    typedef boost::shared_ptr< BasicBlock_t > basic_block_ptr;
    typedef boost::shared_ptr< Statement_t >  statement_ptr;
    typedef std::vector< statement_ptr >      stmt_list_t;
    typedef boost::shared_ptr< stmt_list_t >  stmt_list_ptr;

   public:

    typedef boost::iterator_range< typename bb_id_set_t::iterator > succ_iterator;
    typedef boost::iterator_range< typename bb_id_set_t::iterator > pred_iterator;
    typedef boost::indirect_iterator< typename stmt_list_t::iterator > iterator;
    typedef boost::iterator_range< typename bb_id_set_t::const_iterator > const_succ_iterator;
    typedef boost::iterator_range< typename bb_id_set_t::const_iterator > const_pred_iterator;
    typedef boost::indirect_iterator< typename stmt_list_t::const_iterator > const_iterator;

   private:

    typedef BinaryOp< ZNumber, VariableName >  ZBinaryOp;
    typedef Assignment< ZNumber, VariableName> ZAssignment;
    typedef Assume< ZNumber, VariableName >    ZAssume;
    typedef Havoc< VariableName >              Havoc_t;
    typedef Unreachable< VariableName >        Unreachable_t;

    typedef boost::shared_ptr< ZBinaryOp >     ZBinaryOp_ptr;
    typedef boost::shared_ptr< ZAssignment >   ZAssignment_ptr;
    typedef boost::shared_ptr< ZAssume >       ZAssume_ptr;
    typedef boost::shared_ptr< Havoc_t >       Havoc_ptr;      
    typedef boost::shared_ptr< Unreachable_t > Unreachable_ptr;      
    
    
    BasicBlockLabel m_bb_id;
    stmt_list_ptr   m_stmts;
    bb_id_set_ptr   m_prev, m_next;
    
    void InsertAdjacent (bb_id_set_ptr c, BasicBlockLabel e)
    { 
      if (std::find(c->begin (), c->end (), e) == c->end ())
        c->push_back (e);
    }
    
    void RemoveAdjacent (bb_id_set_ptr c, BasicBlockLabel e)
    {
      if (std::find(c->begin (), c->end (), e) != c->end ())
        c->erase (std::remove(c->begin (), c->end (), e), c->end ());
    }
    
    BasicBlock (BasicBlockLabel bb_id): 
      m_bb_id (bb_id), 
      m_stmts (stmt_list_ptr (new stmt_list_t)), 
      m_prev (bb_id_set_ptr (new bb_id_set_t)), 
      m_next (bb_id_set_ptr (new bb_id_set_t))  
    { }

    static basic_block_ptr Create (BasicBlockLabel bb_id) 
    {
      return basic_block_ptr (new BasicBlock_t (bb_id));
    }
    
    void insert(statement_ptr stmt) 
    {
      this->m_stmts->push_back(stmt);
    }

   public:

    boost::shared_ptr <BasicBlock_t> clone () const
    {
      boost::shared_ptr <BasicBlock_t> b (new BasicBlock_t (label ()));

      for (auto &s : boost::make_iterator_range (begin (), end ()))
        b->m_stmts->push_back (s.clone ()); 

      for (auto id : prev_blocks ()) 
        b->m_prev->push_back (id);
      
      for (auto id : next_blocks ())
        b->m_next->push_back (id);
      
      return b;
    }

    BasicBlockLabel label () const { return m_bb_id; }
    
    iterator begin()             { return boost::make_indirect_iterator (m_stmts->begin ()); }
    iterator end()               { return boost::make_indirect_iterator (m_stmts->end ()); }
    const_iterator begin() const { return boost::make_indirect_iterator (m_stmts->begin ()); }
    const_iterator end()   const { return boost::make_indirect_iterator (m_stmts->end ()); }
    
    std::size_t size() { return std::distance ( begin (), end ()); }
        
    succ_iterator next_blocks ()
    { 
      return boost::make_iterator_range ( m_next->begin (), m_next->end ());
    }
    
    pred_iterator prev_blocks() 
    { 
      return boost::make_iterator_range ( m_prev->begin (), m_prev->end ());
    }

    const_succ_iterator next_blocks () const
    { 
      return boost::make_iterator_range ( m_next->begin (), m_next->end ());
    }
    
    const_pred_iterator prev_blocks() const
    { 
      return boost::make_iterator_range ( m_prev->begin (), m_prev->end ());
    }
    
    void reverse()
    {
      std::swap (m_prev, m_next);
      std::reverse (m_stmts->begin (), m_stmts->end ());
    }
    
    // Add a cfg edge from *this to b
    void operator>>(BasicBlock_t& b) 
    {
      InsertAdjacent (m_next, b.m_bb_id);
      InsertAdjacent (b.m_prev, m_bb_id);
    }
    
    // Remove a cfg edge from *this to b
    void operator-=(BasicBlock_t &b)
    {
      RemoveAdjacent (m_next, b.m_bb_id);
      RemoveAdjacent (b.m_prev, m_bb_id);       
    }
    
    // insert all statements of other at the front
    void merge_front (const BasicBlock_t &other) 
    {
      m_stmts->insert (m_stmts->begin (), 
                       other.m_stmts->begin (), 
                       other.m_stmts->end ());
    }
    
    // insert all statements of other at the back
    void merge_back (const BasicBlock_t &other) 
    {
      m_stmts->insert (m_stmts->end (), 
                       other.m_stmts->begin (), 
                       other.m_stmts->end ());
    }
    
    ostream& write(ostream& o) const
    {
      o << cfg_impl::get_label_str (m_bb_id) << ":\n";	

      for (auto const &s: *this)
        o << "  " << s << ";\n"; 
      
      o << "--> [";

      for (auto const &n : next_blocks ())
        o << cfg_impl::get_label_str (n) << ";";

      o << "]\n";
      
      return o;
    }

    /// To build statements

    void add (ZVariable lhs, ZVariable op1, ZVariable op2) 
    {
      insert(boost::static_pointer_cast< Statement_t, ZBinaryOp >
             (ZBinaryOp_ptr(new ZBinaryOp(lhs, OP_ADDITION, op1, op2))));
    }
    
    void add (ZVariable lhs, ZVariable op1, ZNumber op2) 
    {
      insert(boost::static_pointer_cast< Statement_t, ZBinaryOp > 
             (ZBinaryOp_ptr(new ZBinaryOp(lhs, OP_ADDITION, op1,  op2))));
    }
    
    void add (ZVariable lhs, ZLinearExpression op1, ZLinearExpression op2) 
    {
      if  (op1.get_variable () && op2.get_variable ())
        add (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        add (lhs, (*op1.get_variable ()), op2.constant ());

      else if (op1.is_constant () && op2.get_variable ())
        add (lhs, (*op2.get_variable ()), op1.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        ZLinearExpression rhs(ZNumber(op1.constant () + op2.constant ()));
        this->insert(boost::static_pointer_cast< Statement_t, ZAssignment >
                     (ZAssignment_ptr (new ZAssignment (lhs, rhs))));
      }
      else
        assert(false && "add operands unexpected");
    }
    
    void sub (ZVariable lhs, ZVariable op1, ZVariable op2) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZBinaryOp >
              (ZBinaryOp_ptr (new ZBinaryOp (lhs, OP_SUBTRACTION, op1, op2))));
    }
    
    void sub (ZVariable lhs, ZVariable op1, ZNumber op2) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZBinaryOp >
              (ZBinaryOp_ptr (new ZBinaryOp (lhs, OP_SUBTRACTION, op1, op2))));
    }
    
    void sub (ZVariable lhs, ZLinearExpression op1, ZLinearExpression op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        sub (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        sub (lhs, (*op1.get_variable ()), op2.constant ());        

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        ZLinearExpression rhs (ZNumber (op1.constant () - op2.constant ()));
        insert (boost::static_pointer_cast< Statement_t, ZAssignment >
                (ZAssignment_ptr(new ZAssignment (lhs, rhs))));
      }
      else
        assert(false && "sub operands unexpected");
    }
    
    void mul (ZVariable lhs, ZVariable op1, ZVariable op2) 
    {
      insert(boost::static_pointer_cast< Statement_t, ZBinaryOp >
             (ZBinaryOp_ptr(new ZBinaryOp (lhs, OP_MULTIPLICATION, op1, op2))));
    }
    
    void mul (ZVariable lhs, ZVariable op1, ZNumber op2) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZBinaryOp >
              (ZBinaryOp_ptr(new ZBinaryOp (lhs, OP_MULTIPLICATION, op1, op2))));
    }
    
    void mul(ZVariable lhs, ZLinearExpression op1, ZLinearExpression op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        mul (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        mul (lhs, (*op1.get_variable()), op2.constant());

      else if (op1.is_constant () && op2.get_variable ())
        mul (lhs, (*op2.get_variable ()), op1.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        ZLinearExpression rhs(ZNumber(op1.constant () * op2.constant ()));
        insert(boost::static_pointer_cast< Statement_t, ZAssignment >
               (ZAssignment_ptr(new ZAssignment (lhs, rhs))));
      }
      else
        assert(false && "mul operands unexpected");
    }
    
    void div (ZVariable lhs, ZVariable op1, ZVariable op2) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZBinaryOp >
              (ZBinaryOp_ptr(new ZBinaryOp (lhs, OP_DIVISION, op1, op2))));
    }
    
    void div (ZVariable lhs, ZVariable op1, ZNumber op2) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZBinaryOp >
              (ZBinaryOp_ptr (new ZBinaryOp (lhs, OP_DIVISION, op1, op2))));
    }
    
    void div (ZVariable lhs, ZLinearExpression op1, ZLinearExpression op2) 
    {
      if (op1.get_variable () && op2.get_variable ())
        div (lhs, (*op1.get_variable ()), (*op2.get_variable ()));

      else if (op1.get_variable () && op2.is_constant ())
        div (lhs, (*op1.get_variable ()), op2.constant ());

      else if (op1.is_constant () && op2.is_constant ()) 
      {
        ZLinearExpression rhs (ZNumber (op1.constant() / op2.constant()));
        insert (boost::static_pointer_cast< Statement_t, ZAssignment >
                (ZAssignment_ptr (new ZAssignment (lhs, rhs))));
      }
      assert (false && "div operands unexpected");
    }
    
    void assign (ZVariable lhs, ZLinearExpression rhs) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZAssignment >
              (ZAssignment_ptr (new ZAssignment (lhs, rhs))));
    }
    
    void assume (ZLinearConstraint cst) 
    {
      insert (boost::static_pointer_cast< Statement_t, ZAssume >
              (ZAssume_ptr (new ZAssume (cst))));
    }
    
    void havoc(VariableName lhs) 
    {
      insert (boost::static_pointer_cast< Statement_t, Havoc_t > 
              (Havoc_ptr (new Havoc_t (lhs))));
    }
    
    void unreachable() 
    {
      insert (boost::static_pointer_cast< Statement_t, Unreachable_t > 
              (Unreachable_ptr (new Unreachable_t ())));
    }

    friend ostream& operator<<(ostream &o, const BasicBlock_t &b)
    {
      b.write (o);
      return o;
    }
    
  }; 

  template< class VariableName>
  struct StatementVisitor
  {
    typedef z_number                             ZNumber;
    typedef linear_expression< ZNumber, VariableName > ZLinearExpression;
    typedef BinaryOp < ZNumber, VariableName >   ZBinaryOp;
    typedef Assignment < ZNumber, VariableName > ZAssignment;
    typedef Assume < ZNumber, VariableName >     ZAssume;
    typedef Havoc< VariableName >                Havoc_t;
    typedef Unreachable< VariableName >          Unreachable_t;
    
    virtual void visit (ZBinaryOp&) = 0;
    virtual void visit (ZAssignment&) = 0;
    virtual void visit (ZAssume&) = 0;
    virtual void visit (Havoc_t&) = 0;
    virtual void visit (Unreachable_t&) = 0;

    virtual ~StatementVisitor () { }
  }; 

  template< class BasicBlockLabel, class VariableName >
  class Cfg
  {
   public:

    typedef BasicBlock< BasicBlockLabel, VariableName > BasicBlock_t;   
    typedef Statement < VariableName >                  Statement_t;

    typedef typename BasicBlock_t::succ_iterator succ_iterator;
    typedef typename BasicBlock_t::pred_iterator pred_iterator;

   private:

    typedef Cfg < BasicBlockLabel, VariableName > cfg_t;
    typedef boost::shared_ptr< BasicBlock_t > BasicBlock_ptr;
    typedef boost::unordered_map< BasicBlockLabel, BasicBlock_ptr > BasicBlockMap;
    typedef typename BasicBlockMap::value_type Binding;
    typedef boost::shared_ptr< BasicBlockMap > BasicBlockMap_ptr;

    template<typename P, typename R> 
    struct getSecond : public std::unary_function<P, R>
    {
      getSecond () { }
      R& operator () (const P &p) const { return *(p.second); }
    }; 

   public:

    typedef boost::transform_iterator< getSecond< Binding, BasicBlock_t >, 
                                       typename BasicBlockMap::iterator > iterator;
    typedef boost::transform_iterator< getSecond< Binding, BasicBlock_t >, 
                                       typename BasicBlockMap::const_iterator > const_iterator;

   private:

    BasicBlockLabel    m_entry;
    BasicBlockLabel    m_exit;
    bool               m_has_exit;
    BasicBlockMap_ptr  m_blocks;

    typedef boost::unordered_set< BasicBlockLabel > visited_t;

    template<typename T>
    void dfs_rec (BasicBlockLabel curId, visited_t &visited, T f) 
    {
      if (visited.find (curId) != visited.end ()) return;
      visited.insert (curId);

      BasicBlock_t &cur = get_node (curId);
      f (cur);
      for (auto n : cur.next_blocks ())
        dfs_rec (n, visited, f);
    }

    template<typename T>
    void dfs (T f) 
    {
      visited_t visited;
      dfs_rec (m_entry, visited, f);
    }

    struct PrintBlock 
    {
      ostream &m_o;
      PrintBlock (ostream& o) : m_o (o) { }
      void operator () (const BasicBlock_t& B){ m_o << B ; }
    };

    Cfg (): m_blocks (BasicBlockMap_ptr (new BasicBlockMap)) { }

   public:

    Cfg (BasicBlockLabel entry):  
        m_entry  (entry), 
        m_has_exit (false),
        m_blocks (BasicBlockMap_ptr (new BasicBlockMap))
    {
      m_blocks->insert (Binding (m_entry, 
                                 BasicBlock_t::Create (m_entry)));
    }

    Cfg (BasicBlockLabel entry, BasicBlockLabel exit):  
        m_entry  (entry), 
        m_exit   (exit), 
        m_has_exit (true),
        m_blocks (BasicBlockMap_ptr (new BasicBlockMap))
    {
      m_blocks->insert (Binding (m_entry, 
                                 BasicBlock_t::Create (m_entry)));
    }

    //! copy constructor will make shallow copies so use this method
    //! for deep copies.
    cfg_t clone () const
    {
      cfg_t cfg;
      cfg.m_entry     = m_entry ;
      if (m_has_exit)
        cfg.m_exit      = m_exit ;
      cfg.m_has_exit  = m_has_exit ;

      for (auto const &BB: boost::make_iterator_range (begin (), end ()))
      {
        boost::shared_ptr <BasicBlock_t> copyBB = BB.clone ();
        cfg.m_blocks->insert (Binding (copyBB->label (), copyBB));
      }
      return cfg;
    }
        
    BasicBlockLabel entry() const { return m_entry; } 
    
    bool has_exit () const { return m_has_exit; }
    BasicBlockLabel exit()  const { assert (has_exit ()); return m_exit; } 
    void set_exit (BasicBlockLabel exit) { m_exit = exit; m_has_exit = true;}
    
    succ_iterator next_nodes (BasicBlockLabel bb_id) 
    {
      
      BasicBlock_t& b = get_node(bb_id);
      return b.next_blocks ();
    }
    
    pred_iterator prev_nodes (BasicBlockLabel bb_id) 
    {
      BasicBlock_t& b = get_node(bb_id);
      return b.prev_blocks ();
    }
    
    BasicBlock_t& insert (BasicBlockLabel bb_id) 
    {
      auto it = m_blocks->find (bb_id);
      if (it != m_blocks->end ()) return *(it->second);

      BasicBlock_ptr block = BasicBlock_t::Create (bb_id);
      m_blocks->insert (Binding (bb_id, block));
      return *block;
    }

    void remove (BasicBlockLabel bb_id)
    {
      BasicBlock_t& bb = get_node(bb_id) ;
      
      for (auto id : bb.prev_blocks ())
      { 
        if (bb_id != id)
        {
          BasicBlock_t& p = get_node(id) ;
          p -= bb;
        }
      }
      
      for (auto id : bb.next_blocks ())
      {
        if (bb_id != id)
        {
          BasicBlock_t& s = get_node(id) ;
          bb -= s;
        }
      }

      m_blocks->erase (bb_id);
    }

    
    BasicBlock_t& get_node (BasicBlockLabel bb_id) 
    {
      auto it = m_blocks->find (bb_id);
      assert (it != m_blocks->end () && 
              "basic block not found in the cfg");
      return *(it->second);
    }
    
    
    iterator begin() 
    {
      return boost::make_transform_iterator (m_blocks->begin (), 
                                             getSecond <Binding, BasicBlock_t> ());
    }
    
    iterator end() 
    {
      return boost::make_transform_iterator (m_blocks->end (), 
                                             getSecond <Binding, BasicBlock_t> ());
    }

    const_iterator begin() const
    {
      return boost::make_transform_iterator (m_blocks->begin (), 
                                             getSecond <Binding, BasicBlock_t> ());
    }
    
    const_iterator end() const
    {
      return boost::make_transform_iterator (m_blocks->end (), 
                                             getSecond <Binding, BasicBlock_t> ());
    }

    size_t size () const { return std::distance (begin (), end ()); }
        
    void reverse()
    {
      assert (m_has_exit);     

      std::swap (m_entry, m_exit);
      for (auto &p: *this) 
        p.reverse (); 
    }

    ostream& write (ostream& o) 
    {
      PrintBlock f (o);
      o << "CFG blocks= " << size () << endl;
      dfs (f);
      return o;
    }

    friend ostream& operator<<(ostream &o, 
                               Cfg< BasicBlockLabel, VariableName > &cfg)
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
      typedef typename StatementVisitor<VariableName>::ZBinaryOp     ZBinaryOp;
      typedef typename StatementVisitor<VariableName>::ZAssignment   ZAssignment;
      typedef typename StatementVisitor<VariableName>::ZAssume       ZAssume;
      typedef typename StatementVisitor<VariableName>::Havoc_t       Havoc;
      typedef typename StatementVisitor<VariableName>::Unreachable_t Unreachable;
      bool _has_assert;
      HasAssertVisitor (): _has_assert(false) { }
      void visit(ZBinaryOp&){ }  
      void visit(ZAssignment&) { }
      void visit(ZAssume&) { _has_assert = true; }
      void visit(Havoc&) { }
      void visit(Unreachable&){ }
    };
    
    // Helpers
    bool hasOneChild (BasicBlockLabel b)
    {
      auto It = this->next_nodes (b);
      return (std::distance (It.begin (), It.end ()) == 1);
    }
    
    bool hasOneParent (BasicBlockLabel b)
    {
      auto It = this->prev_nodes (b);
      return (std::distance (It.begin (), It.end ()) == 1);
    }
    
    BasicBlock_t& getChild (BasicBlockLabel b)
    {
      assert (hasOneChild (b));
      auto It = this->next_nodes (b);
      return this->get_node (*(It.begin()));
    }
    
    BasicBlock_t& getParent (BasicBlockLabel b)
    {
      assert (hasOneParent (b));
      auto It = this->prev_nodes (b);
      return this->get_node (*(It.begin()));
    }
    
    void mergeBlocksRec (BasicBlockLabel curId, 
                         visited_t& visited)
    {
      
      if (visited.find (curId) != visited.end ()) return;
      visited.insert (curId);
      
      BasicBlock_t &cur = this->get_node (curId);
      
      if (hasOneChild (curId) && hasOneParent (curId))
      {
        BasicBlock_t &parent = getParent (curId);
        BasicBlock_t &child  = getChild (curId);
        
        HasAssertVisitor vis;
        for (auto it = cur.begin (); it != cur.end (); ++it)
          it->accept(&vis);
        
        // we don't merge two blocks with assertions
        if (!vis._has_assert)
        {
          parent.merge_back (cur);
          this->remove (curId);
          parent >> child;        
          mergeBlocksRec (child.label (), visited); 
          return;
        }
      }
      
      for (auto n : cur.next_blocks ())
        mergeBlocksRec (n, visited);
    }
    
    // Merges a basic block into its predecessor if there is only one
    // and the predecessor only has one successor.
    void mergeBlocks ()
    {
      visited_t visited;
      mergeBlocksRec (this->entry (), visited);
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
      markAliveBlocks (this->entry (), *this, alive);
      
      for (auto const &bb : *this) 
        if (!(alive.count (bb.label ()) > 0))
          dead.insert (bb.label ());
      
      for (auto bb_id: dead)
        this->remove (bb_id);
    }
    
    // remove blocks that cannot reach the exit block
    void removeUselessBlocks ()
    {
      if (!this->has_exit ()) return;

      cfg_t cfg  = this->clone ();
      cfg.reverse ();
      
      visited_t useful, useless;
      markAliveBlocks (cfg.entry (), cfg, useful);
      
      for (auto const &bb : *this) 
        if (!(useful.count (bb.label ()) > 0))
          useless.insert (bb.label ());
      
      for (auto bb_id: useless)
        this->remove (bb_id);
    }
    
  }; 
  
} // end namespace 

#endif /* CFG_HPP */
