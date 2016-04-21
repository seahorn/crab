#ifndef PROPERTY_CHECKER_HPP
#define PROPERTY_CHECKER_HPP

/* 
   Base class for a property checker
 */

#include <crab/cfg/Cfg.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/nullity.hpp>

namespace crab {

  namespace checker {

  using namespace cfg;

  typedef enum { _SAFE, _ERR, _WARN, _UNREACH } check_kind_t;

  // Toy database to store invariants. We may want to replace it with
  // a permanent external database.
  class ChecksDB {
    typedef std::pair<DebugInfo, check_kind_t> check_t;
    typedef std::set<check_t> checks_db_t;

    checks_db_t m_db;
    unsigned m_total_safe;
    unsigned m_total_err;
    unsigned m_total_unreach;
    unsigned m_total_warn;

   public:

    ChecksDB ():
        m_total_safe(0), m_total_err(0), 
        m_total_unreach(0), m_total_warn (0) { }
    
    // add an entry in the database
    void add (check_kind_t status, DebugInfo dbg = DebugInfo () ) {
      switch (status) {
        case _SAFE: m_total_safe++;break;
        case _ERR : m_total_err++;break;
        case _UNREACH: m_total_unreach++;break;
        default: m_total_warn++;
      }
      if (dbg.has_debug ())
        m_db.insert (check_t (dbg, status));
    }
    
    // merge two databases
    void operator+=(const ChecksDB& other) {
      m_db.insert (other.m_db.begin (), other.m_db.end ());
      m_total_safe += other.m_total_safe;
      m_total_err += other.m_total_err;
      m_total_warn += other.m_total_warn;
      m_total_unreach += other.m_total_unreach;
    }
    
    void write (std::ostream& o) const {
      std::vector<unsigned> cnts = { m_total_safe, m_total_err, 
                                     m_total_warn, m_total_unreach };
      unsigned MaxValLen = 0;
      for (auto c: cnts) {
        MaxValLen = std::max(MaxValLen,
                             (unsigned)std::to_string(c).size());
      }
      
      o << std::string ((int) MaxValLen - std::to_string(m_total_safe).size(), ' ') 
        << m_total_safe << std::string (2, ' ') << "Number of total safe checks\n";
      o << std::string ((int) MaxValLen - std::to_string(m_total_err).size(), ' ') 
        << m_total_err << std::string (2, ' ') << "Number of total error checks\n";
      o << std::string ((int) MaxValLen - std::to_string(m_total_warn).size(), ' ') 
        << m_total_warn << std::string (2, ' ') << "Number of total warning checks\n";
      o << std::string ((int) MaxValLen - std::to_string(m_total_unreach).size(), ' ') 
        << m_total_unreach << std::string (2, ' ') << "Number of total unreachable checks\n";

      unsigned MaxFileLen = 0;
      for (auto const &p: m_db) {
        MaxFileLen = std::max(MaxFileLen,
                              (unsigned)p.first.m_file.size());
      }
      
      for (auto const &p: m_db) {
        switch (p.second) {
          case _SAFE:    o << "safe: "; break;
          case _ERR:     o << "error: "; break;
          case _UNREACH: o << "unreachable: "; break;
          default:       o << "warning: "; break;
        }
        o << p.first.m_file << std::string ((int) MaxFileLen - p.first.m_file.size(), ' ') 
          << std::string (2, ' ')
          << " line " << p.first.m_line 
          << " col " << p.first.m_col << "\n";
      }
    }
    
  };

  #define LOG_SAFE(VERBOSE,INV,PROP,DEBUG_LOC)       \
  this->m_db.add (_SAFE);                            \
  if (VERBOSE >=3) {                                 \
    std::cout << " --- SAFE -----------------\n";    \
    if (DEBUG_LOC.has_debug ())                      \
      std::cout << DEBUG_LOC << "\n";                \
    std::cout << "Property : " << PROP << "\n";      \
    std::cout << "Invariant: " << INV << "\n";       \
    std::cout << " -----------------------------\n"; \
  }

  #define LOG_WARN(VERBOSE,INV,PROP,DEBUG_LOC)       \
  this->m_db.add (_WARN, DEBUG_LOC);                 \
  if (VERBOSE >=2) {                                 \
    std::cout << " --- WARNING -----------------\n"; \
    if (DEBUG_LOC.has_debug ())                      \
      std::cout << DEBUG_LOC << "\n";                \
    std::cout << "Property : " << PROP << "\n";      \
    std::cout << "Invariant: " << INV << "\n";       \
    std::cout << " -----------------------------\n"; \
  }

  #define LOG_ERR(VERBOSE,INV,PROP,DEBUG_LOC)        \
    this->m_db.add (_ERR, DEBUG_LOC);                \
  if (VERBOSE >=1) {                                 \
    std::cout << " --- ERROR -----------------\n";   \
    if (DEBUG_LOC.has_debug ())                      \
      std::cout << DEBUG_LOC << "\n";                \
    std::cout << "Property : " << PROP << "\n";      \
    std::cout << "Invariant: " << INV << "\n";       \
    std::cout << " -----------------------------\n"; \
  }


  namespace num_dom_detail {
     // To avoid compilation errors if the abstract domain is not a
     // numerical domain
     template<typename Domain, typename VariableName>
     struct GetAs {
       typedef crab::domains::interval<z_number> interval_t;
       typedef linear_constraint< z_number, VariableName> z_lin_cst_t;

       Domain& m_inv;
       GetAs (Domain& inv): m_inv (inv) { }

       interval_t operator[](VariableName v){ return m_inv [v]; }

       bool entails(z_lin_cst_t cst) { 
         if (cst.is_tautology ()) return true;

         if (cst.is_contradiction ()) return false;

         Domain inv;
         inv += cst;
         if (inv.is_top ()) return false;
         return m_inv <= inv; 
       }

       bool intersect(z_lin_cst_t cst) {
         if (cst.is_tautology ()) return true;

         if (cst.is_contradiction ()) return false;
         
         Domain inv;
         inv += cst;
         Domain meet = m_inv & inv;
         return (!meet.is_bottom ());
       }
     };

     // TO BE COMPLETED: add here any non-numerical domain
     template<typename VariableName>
     struct GetAs <crab::domains::nullity_domain<VariableName>, VariableName> {
       typedef crab::domains::interval<z_number> interval_t;
       typedef linear_constraint<z_number,VariableName> z_lin_cst_t;
       typedef crab::domains::nullity_domain<VariableName> nullity_domain_t;
       
       GetAs (nullity_domain_t&) { }
       interval_t operator[](VariableName){ return interval_t::top (); }
       bool entails(z_lin_cst_t) { return false; }
       bool intersect(z_lin_cst_t) { return true; }  
       
     };
  } // end num_dom_detail namespace

  template<typename Analyzer>
  class PropertyChecker: public StatementVisitor <typename Analyzer::varname_t> {
   public:
    typedef typename Analyzer::abs_tr_ptr abs_tr_ptr;
    typedef typename Analyzer::varname_t varname_t;
    typedef typename Analyzer::abs_dom_t abs_dom_t;

    typedef variable < z_number, varname_t > z_var_t;
    typedef linear_expression< z_number, varname_t > z_lin_exp_t;
    typedef linear_constraint< z_number, varname_t > z_lin_cst_t;
    typedef linear_constraint_system< z_number, varname_t > z_lin_cst_sys_t;
    typedef BinaryOp <z_number,varname_t> z_bin_op_t;
    typedef Assignment <z_number,varname_t> z_assign_t;
    typedef Assume <z_number,varname_t> z_assume_t;
    typedef Assert <z_number,varname_t> z_assert_t;
    typedef Havoc<varname_t> havoc_t;
    typedef Unreachable<varname_t> unreach_t;
    typedef Select <z_number,varname_t> z_select_t;
    typedef FCallSite<varname_t> callsite_t;
    typedef Return<varname_t> return_t;
    typedef ArrayInit<varname_t> z_arr_init_t;
    typedef AssumeArray<z_number,varname_t> z_assume_arr_t;
    typedef ArrayStore<z_number,varname_t> z_arr_store_t;
    typedef ArrayLoad<z_number,varname_t> z_arr_load_t;
    typedef PtrStore<z_number,varname_t> z_ptr_store_t;
    typedef PtrLoad<z_number,varname_t> z_ptr_load_t;
    typedef PtrAssign<z_number,varname_t> z_ptr_assign_t;
    typedef PtrObject<varname_t> ptr_object_t;
    typedef PtrFunction<varname_t> ptr_function_t;
    typedef PtrNull<varname_t> ptr_null_t;
    typedef PtrAssume<varname_t> ptr_assume_t;
    typedef PtrAssert<varname_t> ptr_assert_t;

    typedef std::set<std::pair<DebugInfo, check_kind_t> > CheckResultsDB;
 
   protected: 

    abs_tr_ptr m_abs_tr; // it can be null
    int m_verbose;
    ChecksDB m_db; // Store debug information about the checks
   
    virtual void check (z_assert_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }

    virtual void check (z_bin_op_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    } 
      
    virtual void check (z_assign_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_assume_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (havoc_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (unreach_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_select_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (callsite_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (return_t& s) {
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt 
    }
      
    virtual void check (z_arr_init_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_assume_arr_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_arr_store_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_arr_load_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_ptr_store_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_ptr_load_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (z_ptr_assign_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_object_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_function_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_null_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_assume_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_assert_t&s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
   public: 

    /* Visitor API */
    void visit (z_bin_op_t &s) { check (s); }
    void visit (z_assign_t &s) { check (s); }
    void visit (z_assume_t &s) { check (s); }
    void visit (havoc_t &s) { check (s); }
    void visit (unreach_t &s) { check (s); }
    void visit (z_select_t &s) { check (s); }
    void visit (z_assert_t &s) { check (s); }
    void visit (callsite_t &s) { check (s); }
    void visit (return_t &s) { check (s); }
    void visit (z_arr_init_t &s) { check (s); }
    void visit (z_assume_arr_t &s) { check (s); }
    void visit (z_arr_store_t &s) { check (s); }
    void visit (z_arr_load_t &s) { check (s); }
    void visit (z_ptr_store_t &s) { check (s); }
    void visit (z_ptr_load_t &s) { check (s); }
    void visit (z_ptr_assign_t &s) { check (s); }
    void visit (ptr_object_t &s) { check (s); }
    void visit (ptr_function_t &s) { check (s); }
    void visit (ptr_null_t &s) { check (s); }
    void visit (ptr_assert_t &s) { check (s); }
    
    PropertyChecker (int verbose): 
        m_abs_tr (nullptr), m_verbose (verbose) { }
    
    void set (abs_tr_ptr abs_tr) {
      m_abs_tr = abs_tr;
    }
    
    const ChecksDB& get_db () const { return m_db; }
    
    ChecksDB get_db () { return m_db; }

    virtual std::string get_property_name () const {
       return "dummy property";
    }

    void write (std::ostream& o) const {
      o << get_property_name () << "\n";
      m_db.write (o);
    }

  };  
  } // end namespace
} // end namespace
#endif 
