#ifndef MUAZ_ANALYZER_HPP
#define MUAZ_ANALYZER_HPP

#include <ikos/muaz.hpp>  
#include <ikos/fwd_fixpoint_iterators.hpp>         
#include <ikos/linear_constraints.hpp> 
#include <ikos/intervals.hpp>                      
#include <ikos/intervals_congruences.hpp>                      
#include <ikos/octagons.hpp>                      
#include <ikos/dbm.hpp>                      

using namespace std;

namespace ikos_impl
{
  using namespace ikos::muaz;

  // Visitor for analyzing statements (i.e., apply transfer functions)
  typedef statement_visitor< varname_t, string >   statement_visitor_t;

  // CFG
  typedef variable< z_number, varname_t >          z_var;
  typedef muaz_cfg< varname_t, string >            cfg_t;
  typedef cfg_t::basic_block_id_t                  basic_block_id_t;
  typedef cfg_t::basic_block_t                     basic_block_t;
  typedef cfg_t::statement_t                       statement_t;

  // Manipulate linear expressions
  typedef linear_expression< z_number, varname_t > linear_expression_t;

  // Numerical domains
  typedef interval_domain< z_number, varname_t >              interval_domain_t;
  typedef interval_congruence_domain< z_number, varname_t >   interval_congruences_domain_t;
  typedef DBM< z_number, varname_t >                          dbm_domain_t;
  typedef octagon< z_number, varname_t >                      octagon_domain_t;
}


namespace ikos
{

  using namespace ikos_impl;
  using namespace ikos::muaz;

  template<typename NumAbsDomain>
  class StatementAnalyzer: public statement_visitor_t
  {

    NumAbsDomain  _inv;

   public:
    
    StatementAnalyzer (NumAbsDomain inv): 
      statement_visitor_t (),
      _inv (inv)
    { }
    
    NumAbsDomain inv() { return this->_inv; }

    void visit(z_binary_operation_t& stmt) 
    {
      _inv.apply (stmt.operation(), 
                  stmt.lhs().name(), 
                  stmt.left_operand().name(), stmt.right_operand().name());
    }
    
    void visit(z_linear_assignment_t& stmt) 
    {
      _inv.assign (stmt.lhs().name(), linear_expression_t(stmt.rhs()));
    }
    
    void visit(z_linear_assertion_t& stmt) 
    {
      _inv += stmt.constraint();
    }
    
    void visit(z_array_write_t &stmt)  { }
    
    void visit(z_array_read_t  &stmt)  { }
    
    void visit(checkpoint_t&) { }
    
  }; 


   template<typename AbsDomain>
   class FwdAnalyzer: 
      public interleaved_fwd_fixpoint_iterator< basic_block_id_t, cfg_t, AbsDomain > 
   {
     
    public:
     
     typedef interleaved_fwd_fixpoint_iterator< basic_block_id_t, cfg_t, AbsDomain > fwd_iterator_t;
     
    private:

     VariableFactory & _vfac;
     
     //! Given a basic block and the invariant at the entry it produces
     //! the invariant at the exit of the block.
     AbsDomain analyze(basic_block_id_t node, AbsDomain pre) 
     { 
       basic_block_t b = this->get_cfg().get_node(node);
       StatementAnalyzer<AbsDomain> vis (pre);
       for (auto &s : b) { s.accept (&vis); }
       return vis.inv();
     } 
     
     void process_pre(basic_block_id_t node, AbsDomain inv) 
     { 
       cout << "Pre at " << node << ": " << inv << endl;
     }
     
     void process_post(basic_block_id_t node, AbsDomain inv) 
     { 
       cout << "Post at " << node << ": " << inv << endl;
     }
     
    public:
     
     FwdAnalyzer (cfg_t cfg,  VariableFactory &vfac): 
         fwd_iterator_t (cfg), _vfac(vfac) { }
     
     void Run (AbsDomain inv) 
     {
       try { this->run (inv); }
       catch (error &e)
       {
         cerr << "ikos error: " << e << endl;
         exit (EXIT_FAILURE);
       }
     }

   }; 
} // end namespace

#endif 
