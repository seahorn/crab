#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/Liveness.hpp>
#include <ikos/domains/domain_traits.hpp>

namespace analyzer
{
  using namespace cfg;
  using namespace std;

  template<typename VariableName, typename NumAbsDomain>
  class StatementAnalyzer: public StatementVisitor <VariableName>
  {
    NumAbsDomain m_inv;

    using typename StatementVisitor<VariableName>::ZLinearExpression;

    // Statements
    using typename StatementVisitor<VariableName>::ZBinaryOp;
    using typename StatementVisitor<VariableName>::ZAssignment;
    using typename StatementVisitor<VariableName>::ZAssume;
    using typename StatementVisitor<VariableName>::Havoc_t;
    using typename StatementVisitor<VariableName>::Unreachable_t;
    using typename StatementVisitor<VariableName>::ArrayInit_t;
    using typename StatementVisitor<VariableName>::ZArrayStore;
    using typename StatementVisitor<VariableName>::ZArrayLoad;

   public:
    
    StatementAnalyzer (NumAbsDomain inv): 
        StatementVisitor<VariableName> (), m_inv (inv)
    { }
    
    NumAbsDomain inv() { return m_inv; }

    void visit(ZBinaryOp& stmt) 
    {
      ZLinearExpression op1 = stmt.left ();
      ZLinearExpression op2 = stmt.right ();
      
      if (op1.get_variable () && op2.get_variable ())
      {
        m_inv.apply (stmt.op (), 
                     stmt.lhs ().name(), 
                     (*op1.get_variable ()).name(), 
                     (*op2.get_variable ()).name());
      }
      else
      {
        assert ( op1.get_variable () && op2.is_constant ());
        m_inv.apply (stmt.op (), 
                     stmt.lhs ().name (), 
                     (*op1.get_variable ()).name(), 
                     op2.constant ()); 
      }      
    }
    
    void visit(ZAssignment& stmt) 
    {
      m_inv.assign (stmt.lhs().name(), 
                    ZLinearExpression(stmt.rhs()));
    }
    
    void visit(ZAssume& stmt) 
    {
      m_inv += stmt.constraint();
    }

    void visit(Havoc_t& stmt) 
    {
      m_inv -= stmt.variable();
    }

    void visit(Unreachable_t& stmt) 
    {
      m_inv = NumAbsDomain::bottom ();
    }

    void visit(ArrayInit_t & stmt) 
    {
      domain_traits::array_init (m_inv, stmt.variable ());
    }

    void visit(ZArrayStore & stmt) 
    {
      if (stmt.index ().get_variable ())
      {
        auto arr = stmt.array ().name ();
        auto idx = *(stmt.index ().get_variable ());
        domain_traits::array_store (m_inv, 
                                    arr,
                                    idx.name(), 
                                    stmt.value (),
                                    stmt.is_singleton ());
      }
    }

    void visit(ZArrayLoad & stmt) 
    {
      if (stmt.index ().get_variable ())
      {
        auto idx = *(stmt.index ().get_variable ());
        domain_traits::array_load (m_inv, 
                                   stmt.lhs ().name (), 
                                   stmt.array ().name (), 
                                   idx.name ());
      }
    }

  }; 


   template< typename BasicBlockLabel, typename VariableName, typename CFG, 
             typename VariableFactory, typename AbsDomain>
   class FwdAnalyzer: 
      public interleaved_fwd_fixpoint_iterator< BasicBlockLabel, CFG, AbsDomain > 
   {
     
    public:
     
     typedef interleaved_fwd_fixpoint_iterator< BasicBlockLabel, CFG, AbsDomain > fwd_iterator_t;

    private:

     typedef boost::unordered_map< BasicBlockLabel, AbsDomain> invariant_map_t;    

     typedef Liveness< BasicBlockLabel, CFG, VariableName> liveness_t;     
     typedef typename liveness_t::live_set_t live_set_t;     

     CFG              m_cfg;
     VariableFactory& m_vfac;
     liveness_t       m_live;
     invariant_map_t  m_pre_map;
     
     //! Given a basic block and the invariant at the entry it produces
     //! the invariant at the exit of the block.
     AbsDomain analyze (BasicBlockLabel node, AbsDomain pre) 
     { 
       BasicBlock<BasicBlockLabel, VariableName> &b = m_cfg.get_node (node);
       StatementAnalyzer<VariableName, AbsDomain> vis (pre);
       for (auto &s : b) { s.accept (&vis); }
       AbsDomain post = vis.inv ();

       // prune dead variables 
       if (post.is_bottom() || post.is_top()) return post;

       auto dead = m_live.dead_exit (node);
       domain_traits::forget (post, dead.begin (), dead.end ());
       return post;
     } 
     
     void process_pre (BasicBlockLabel node, AbsDomain inv) 
     {//cout << "Pre at " << node << ": " << inv << endl; 
       auto it = m_pre_map.find (node);
       if (it == m_pre_map.end())
         m_pre_map.insert(typename invariant_map_t::value_type (node, inv));
     }
     
     void process_post (BasicBlockLabel node, AbsDomain inv) 
     {//cout << "Post at " << node << ": " << inv << endl; 
     }
     
    public:
     
     FwdAnalyzer (CFG cfg,  VariableFactory &vfac, bool runLive=false): 
         fwd_iterator_t (cfg), 
         m_cfg (cfg), m_vfac(vfac), m_live (m_cfg)
     { 
       if (runLive)
         m_live.exec ();
     }

     //! Trigger the fixpoint computation 
     void Run (AbsDomain inv) 
     { this->run (inv); }      
     

     //! Propagate invariants at the statement level
     template < typename Statement >
     AbsDomain AnalyzeStmt (Statement stmt, AbsDomain pre) {
       StatementAnalyzer<VariableName, AbsDomain> vis (pre);
       vis.visit (stmt);
       return vis.inv ();
     }

     //! Return the invariants that hold at the entry of bb
     AbsDomain operator[] (BasicBlockLabel b) const
     {
       auto it = m_pre_map.find (b);
       if (it == m_pre_map.end ())
         return AbsDomain::top ();
       else
         return it->second;
     }

   }; 
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
