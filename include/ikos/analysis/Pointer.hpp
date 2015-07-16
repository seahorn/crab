#ifndef POINTER_ANALYSIS_HPP
#define POINTER_ANALYSIS_HPP

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

#include <ikos/common/types.hpp>
#include <ikos/domains/pta.hpp>
#include <ikos/domains/intervals.hpp>
#include <ikos/domains/discrete_domains.hpp>

namespace analyzer
{

  using namespace std;
  using namespace cfg;

  class PtrAddr {
    index_t m_addr;

   public:

    PtrAddr (index_t addr): m_addr (addr) { }

    bool operator==(PtrAddr other) const { return m_addr == other.m_addr;}

    bool operator<(PtrAddr other) const  { return m_addr < other.m_addr;}

    size_t index () const { return m_addr; }

    size_t hash_value () const { return index (); }

    void write (ostream &o){
      o << "&(" << m_addr << ")";
    }
  };

  //! Run a flow-insensitive pointer analysis for a set of CFGs.
  template<typename CFG, typename VariableFactory>
  class Pointer
  {
    typedef typename CFG::varname_t varname_t;

    typedef interval_domain <z_number, varname_t> num_domain_t;
    typedef interval <z_number> interval_t;

    //! Generate points-to constraints from statements using
    //! invariants to resolve offsets.
    template < typename NumInvGen>
    struct GenBasicBlockCons: public StatementVisitor <varname_t> 
    {
      using typename StatementVisitor<varname_t>::z_bin_op_t;
      using typename StatementVisitor<varname_t>::z_assign_t;
      using typename StatementVisitor<varname_t>::z_assume_t;
      using typename StatementVisitor<varname_t>::havoc_t;
      using typename StatementVisitor<varname_t>::unreach_t;
      using typename StatementVisitor<varname_t>::z_select_t;
      
      using typename StatementVisitor<varname_t>::callsite_t;
      using typename StatementVisitor<varname_t>::return_t;
      using typename StatementVisitor<varname_t>::ptr_load_t;
      using typename StatementVisitor<varname_t>::ptr_store_t;
      using typename StatementVisitor<varname_t>::ptr_assign_t;
      using typename StatementVisitor<varname_t>::ptr_object_t;
      using typename StatementVisitor<varname_t>::ptr_function_t;
      
      typedef FunctionDecl<varname_t> FunctionDecl_t;
      typedef boost::unordered_map< varname_t, pointer_var > pt_var_map_t;
      
      VariableFactory& m_vfac;
      pta_system* m_cs; 
      num_domain_t m_inv;
      NumInvGen* m_inv_gen; 
      pt_var_map_t*m_pt_map; 
      boost::optional<varname_t> m_func_name; 

      //! Return the interval [0,0]
      interval <z_number> zero() const { return interval <z_number> (0,0);}
      
      //! Create a pointer variable
      pointer_var new_pointer_var (varname_t v) {
        auto it = m_pt_map->find(v);
        if (it != m_pt_map->end())
          return it->second;
        
        pointer_var pt = ikos::mk_pointer_var (v.index());
        m_pt_map->insert (make_pair (v, pt));
        return pt;
      }
      
      //! Create a new parameter
      pointer_var new_param_ref (varname_t fname, unsigned param) {
        auto par_ref = ikos::mk_param_ref(new_pointer_var(fname), param);
        std::ostringstream buf;
        par_ref->print (buf);
        pointer_var p = new_pointer_var (m_vfac[buf.str()]);
        return p;
      }
      
      //! Create a new return
      pointer_var new_return_ref (varname_t fname) {
        auto ret_ref = ikos::mk_return_ref(new_pointer_var(fname));
        std::ostringstream buf;
        ret_ref->print (buf);
        pointer_var p = new_pointer_var (m_vfac[buf.str()]);
        return p;
      }
      
      GenBasicBlockCons (VariableFactory& vfac,
                         pta_system* cs, 
                         num_domain_t inv, 
                         NumInvGen* inv_gen, 
                         pt_var_map_t* pt_map,
                         boost::optional<varname_t> func_name): 
          m_vfac (vfac), m_cs (cs), m_inv (inv), m_inv_gen (inv_gen), 
          m_pt_map (pt_map), m_func_name (func_name)
      { }


      void gen_func_decl_cons  (FunctionDecl_t & decl) { 
        unsigned num_params = decl.get_num_params ();
        for (unsigned i=0; i<num_params; i++)
        {
          if (decl.get_param_type (i) == PTR_TYPE)
          {
            pointer_var fp = new_pointer_var (decl.get_param_name (i));
            (*m_cs) += fp == (new_param_ref (decl.get_func_name (), i) + zero());              
          }
        }
      }
      
      void visit (callsite_t & stmt) { 
        unsigned num_args = stmt.get_num_args ();
        for (unsigned i=0; i<num_args; i++)
        {
          if (stmt.get_arg_type (i) == PTR_TYPE)
          {
            pointer_var ap = new_pointer_var (stmt.get_arg_name (i));
            pointer_var fp = new_param_ref (stmt.get_func_name (), i);
            (*m_cs) += fp == (ap + zero());              
          }
        }
        
        auto lhs_name = stmt.get_lhs_name ();
        if (lhs_name && (stmt.get_lhs_type () == PTR_TYPE))
        {
          pointer_var p = new_pointer_var (*lhs_name);
          (*m_cs) += (p == new_return_ref (stmt.get_func_name ()) + zero ());
        }
        
      }
      
      void visit (return_t & stmt) { 
        if (stmt.get_ret_type () == PTR_TYPE && (m_func_name))
        {
          pointer_var p = new_pointer_var (stmt.get_ret_var ());
          (*m_cs) += (new_return_ref (*m_func_name) == p + zero ());
        }

      }
      
      void visit (ptr_object_t & stmt) { 
        pointer_var lhs = new_pointer_var (stmt.lhs ());
        auto ref = ikos::mk_object_ref (stmt.rhs (), zero());
        (*m_cs) += lhs == ref;
      }
      
      void visit (ptr_function_t & stmt) { 
        pointer_var lhs = new_pointer_var (stmt.lhs ());
        auto ref = ikos::mk_function_ref (stmt.rhs ());
        (*m_cs) += lhs == ref;
        }
      
      void visit (ptr_assign_t & stmt) { 
        pointer_var lhs = new_pointer_var (stmt.lhs ());
        pointer_var rhs = new_pointer_var (stmt.rhs ());
        interval_t  offset = interval_t::top ();
        if (stmt.offset().get_variable ())
          offset = m_inv [(*stmt.offset().get_variable ()).name ()];
        else if (stmt.offset ().is_constant ())
          offset = interval_t (stmt.offset ().constant ());
        else
        {
          cerr << "Pointer arithmetic does not allow arbitrary index expressions.\n";
          std::exit (EXIT_FAILURE);
        }
        (*m_cs) += lhs == (rhs + offset );
      }
      
      void visit (ptr_load_t & stmt) {
        pointer_var lhs = new_pointer_var (stmt.lhs ());
        pointer_var rhs = new_pointer_var (stmt.rhs ());
        interval_t size = stmt.size ();
        (*m_cs) += lhs *= (rhs + size);
      }
      
      void visit (ptr_store_t & stmt) { 
        pointer_var lhs = new_pointer_var (stmt.lhs ());
        pointer_var rhs = new_pointer_var (stmt.rhs ());
        interval_t size = stmt.size ();
        (*m_cs) += (lhs + size) << rhs;          
        }
      
      void visit (z_bin_op_t& stmt) 
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }

      void visit (z_assign_t& stmt)
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }
      
      void visit (z_assume_t& stmt)
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }
      
      void visit (havoc_t& stmt)
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }
      
      void visit (unreach_t& stmt)
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }

      void visit (z_select_t& stmt) 
      { m_inv = m_inv_gen->AnalyzeStmt (stmt, m_inv); }
      
    };

    typedef NumAbsTransformer <varname_t, num_domain_t> num_abs_tr_t;
    typedef FwdAnalyzer<CFG, num_abs_tr_t, VariableFactory> num_inv_gen_t;
    typedef typename GenBasicBlockCons < num_inv_gen_t>::pt_var_map_t pt_var_map_t;

    //! for external queries
    typedef boost::unordered_map< varname_t,
                                  pair< discrete_domain< PtrAddr >, interval_t > >
    ptr_map_t;

    VariableFactory& m_vfac;
    pta_system m_cs;
    pt_var_map_t m_pt_var_map; 

    ptr_map_t m_ptr_map; //! to store results

    typedef GenBasicBlockCons<num_inv_gen_t> gen_bb_cons_t;

    public:
    
      Pointer (VariableFactory &vfac): m_vfac (vfac)
      { }
    
    
      //////
      // Generation of constraints for functions
      //////
      // for each function definition "func(p1,...,pn){...}" do
      //   p1 >= param(func,1) 
      //   ...
      //   pn >= param(func,n) 
      //
      // for each call site "x = func(a1,,...an);"
      //   param (func,1) >= a1
      //   ...
      //   param (func,1) >= an
      //   x >= ret (func)
      //
      // for each given return statement "ret p"
      //   ret (func) >= p
      ///////
    
      void gen_constraints (CFG cfg) 
      {
        // 1) Run a numerical analysis to infer numerical invariants
        //    about offsets.
        const bool run_live = true;
        num_inv_gen_t It (cfg, m_vfac, run_live);
        It.Run (num_domain_t::top ());
        
        // 2) Gen points-to constraints
        auto func_decl = cfg.get_func_decl ();
        boost::optional<varname_t> func_name;
        if (func_decl)
        {
          gen_bb_cons_t vis (m_vfac, &m_cs, num_domain_t::top (), 
                             &It, &m_pt_var_map, func_name);
          vis.gen_func_decl_cons (*func_decl);
          func_name = (*func_decl).get_func_name ();
        }

        for (auto &b : cfg)
        {
          gen_bb_cons_t vis (m_vfac, &m_cs, It [b.label ()], 
                             &It, &m_pt_var_map, func_name);
          for (auto &s: b) { s.accept (&vis); }
        }
      }

      void solve () {
        
        // 3) Solve points-to set constraints
        // cout << "Generated constraints:\n" << m_cs << endl;

        m_cs.solve ();
        
        // 4) Store results for future queries
        for (auto pt_to_var : m_pt_var_map)
        {
          address_set pt_set = m_cs.get (pt_to_var.second).first;
          auto pt_varset = ((pt_set.begin() == pt_set.end())
                            ? discrete_domain<PtrAddr>::top()
                            : discrete_domain<PtrAddr>::bottom());
          
          interval_t offset = m_cs.get (pt_to_var.second).second;
          if (offset.is_bottom())
            offset = interval_t::top();
          
          for (auto pi : pt_set) { pt_varset = pt_varset | PtrAddr (pi); }
          auto v = make_pair (pt_to_var.first, make_pair(pt_varset, offset));
          m_ptr_map.insert(v);
        }
      }

      pair< discrete_domain< varname_t >, 
            interval_t > operator[] (varname_t p) const
      {
        return m_ptr_map [p];
      }
    
      void write (ostream& o) const
      {
        o << "{";
        for (auto e : m_ptr_map)
        {
          o << e.first << "->" << "(" 
            << e.second.first << "," << e.second.second << ")" << ";";
        }
        o << "}";
      }
   }; 

   template <typename CFG, typename VariableFactory >
   ostream& operator << (ostream& o, const Pointer<CFG, VariableFactory> &pta)
   {
     pta.write (o);
     return o;
   }

} // end namespace 

#endif 
