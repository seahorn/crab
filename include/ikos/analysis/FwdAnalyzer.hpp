#ifndef FWD_ANALYZER_HPP
#define FWD_ANALYZER_HPP

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/Liveness.hpp>

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

     typedef boost::unordered_map< BasicBlockLabel, AbsDomain> invariant_map_t;    

     typedef Liveness< BasicBlockLabel, CFG, VariableName> liveness_t;     
     typedef typename liveness_t::live_set_t live_set_t;     

     CFG              m_cfg;
     VariableFactory& m_vfac;
     liveness_t       m_live;
     invariant_map_t  m_pre_map;
     

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
     {
       //cout << "Pre at " << node << ": " << inv << endl; 
       auto it = m_pre_map.find (node);
       if (it == m_pre_map.end())
       {
         prune_dead_variables (node, inv);
         m_pre_map.insert(typename invariant_map_t::value_type (node, inv));
       }
     }
     
     void process_post (BasicBlockLabel node, AbsDomain inv) 
     { //cout << "Post at " << node << ": " << inv << endl; 
     }
     
    public:
     
     FwdAnalyzer (CFG cfg,  VariableFactory &vfac, bool runLive=false): 
         fwd_iterator_t (cfg), 
         m_cfg (cfg), m_vfac(vfac), m_live (m_cfg)
     { 
       if (runLive)
         m_live.exec ();
     }
     
     void Run (AbsDomain inv) 
     { this->run (inv); }      
     

    //! return the invariants that hold at the entry of bb
    AbsDomain operator[] (BasicBlockLabel b) const
    {
      auto it = m_pre_map.find (b);
      if (it == m_pre_map.end ())
        return AbsDomain::bottom ();
      else
        return it->second;
    }

   }; 
} // end namespace

#endif /* FWD_ANALYZER_HPP*/
