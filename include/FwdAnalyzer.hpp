#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <include/Cfg.hpp>
#include <include/VarFactory.hpp>
#include <include/Liveness.hpp>

namespace analyzer
{
  using namespace cfg;
  using namespace std;

  template<typename VariableName, typename NumAbsDomain>
  class StatementAnalyzer: public StatementVisitor <VariableName>
  {
    NumAbsDomain  m_inv;

    using typename StatementVisitor<VariableName>::ZBinaryOp;
    using typename StatementVisitor<VariableName>::ZAssignment;
    using typename StatementVisitor<VariableName>::ZAssume;
    using typename StatementVisitor<VariableName>::Havoc_t;
    using typename StatementVisitor<VariableName>::Unreachable_t;
    using typename StatementVisitor<VariableName>::ZLinearExpression;

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
    
  }; 


   template< typename BasicBlockLabel, typename VariableName, typename CFG, 
             typename VariableFactory, typename AbsDomain>
   class FwdAnalyzer: 
      public interleaved_fwd_fixpoint_iterator< BasicBlockLabel, CFG, AbsDomain > 
   {
     
    public:
     
     typedef interleaved_fwd_fixpoint_iterator< BasicBlockLabel, CFG, AbsDomain > fwd_iterator_t;

    private:

     typedef Liveness< BasicBlockLabel, CFG, VariableName> liveness_t;     
     typedef typename liveness_t::live_set_t live_set_t;     

     CFG              m_cfg;
     VariableFactory& m_vfac;
     liveness_t       m_live;

     //! Using liveness information delete all dead variables from inv
     void prune_dead_variables (BasicBlockLabel node_name, AbsDomain &inv)
     {
       if (inv.is_bottom() || inv.is_top()) return;
       for (auto v: m_live.dead (node_name)) { inv -= v;  }
     }

     //! Given a basic block and the invariant at the entry it produces
     //! the invariant at the exit of the block.
     AbsDomain analyze (BasicBlockLabel node, AbsDomain pre) 
     { 
       prune_dead_variables(node, pre);

       BasicBlock<BasicBlockLabel, VariableName> &b = m_cfg.get_node (node);
       StatementAnalyzer<VariableName, AbsDomain> vis (pre);
       for (auto &s : b) { s.accept (&vis); }
       return vis.inv();
     } 
     
     void process_pre (BasicBlockLabel node, AbsDomain inv) 
     { cout << "Pre at " << node << ": " << inv << endl; }
     
     void process_post (BasicBlockLabel node, AbsDomain inv) 
     { cout << "Post at " << node << ": " << inv << endl; }
     
    public:
     
     FwdAnalyzer (CFG cfg,  VariableFactory &vfac, bool runLive=true): 
         fwd_iterator_t (cfg), 
         m_cfg (cfg), m_vfac(vfac), m_live (m_cfg)
     { 
       if (runLive)
       {
         //cout << "Liveness=\n";
         m_live.exec ();
         //cout << m_live << "\n";
       }
     }
     
     void Run (AbsDomain inv) 
     {       
       cout << "Running " << inv.getDomainName () << "... \n";
       this->run (inv); 
     }

   }; 
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
