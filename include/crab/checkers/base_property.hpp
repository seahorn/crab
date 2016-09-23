#ifndef PROPERTY_CHECKER_HPP
#define PROPERTY_CHECKER_HPP

/* 
   Base class for a property checker
 */

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/nullity.hpp>

namespace crab {

  namespace checker {

  using namespace cfg;

  typedef enum { _SAFE, _ERR, _WARN, _UNREACH } check_kind_t;

  // Toy database to store invariants. We may want to replace it with
  // a permanent external database.
  class checks_db {
    typedef std::pair<debug_info, check_kind_t> check_t;
    typedef std::set<check_t> checks_db_t;

    checks_db_t m_db;
    unsigned m_total_safe;
    unsigned m_total_err;
    unsigned m_total_unreach;
    unsigned m_total_warn;

   public:

    checks_db ()
        : m_total_safe(0), m_total_err(0), 
          m_total_unreach(0), m_total_warn (0) { }

    unsigned get_total_safe () const { return m_total_safe + m_total_unreach; }

    unsigned get_total_warning () const { return m_total_warn; }

    unsigned get_total_error () const { return m_total_err; }

    // add an entry in the database
    void add (check_kind_t status, debug_info dbg = debug_info () ) {
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
    void operator+=(const checks_db& other) {
      m_db.insert (other.m_db.begin (), other.m_db.end ());
      m_total_safe += other.m_total_safe;
      m_total_err += other.m_total_err;
      m_total_warn += other.m_total_warn;
      m_total_unreach += other.m_total_unreach;
    }
    
    void write (crab_os& o) const {
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
        // print all checks here
        // o << p.first.m_file << std::string ((int) MaxFileLen - p.first.m_file.size(), ' ') 
        //   << std::string (2, ' ')
        //   << " line " << p.first.m_line 
        //   << " col " << p.first.m_col << "\n";
      }
    }
    
  };

  #define LOG_SAFE(VERBOSE,INV,PROP,DEBUG_LOC)           \
  do {                                                   \
    this->m_db.add (_SAFE);                              \
    if (VERBOSE >=3) {                                   \
      crab::outs() << " --- SAFE -----------------\n";   \
      if (DEBUG_LOC.has_debug ())                        \
         crab::outs() << DEBUG_LOC << "\n";              \
      crab::outs() << "Property : " << PROP << "\n";     \
      crab::outs() << "Invariant: " << INV << "\n";      \
      crab::outs() << " -----------------------------\n";\
    }                                                    \
  } while (false);

  #define LOG_WARN(VERBOSE,INV,PROP,DEBUG_LOC)            \
  do {                                                    \
    this->m_db.add (_WARN, DEBUG_LOC);                    \
    if (VERBOSE >=2) {                                    \
      crab::outs() << " --- WARNING -----------------\n"; \
      if (DEBUG_LOC.has_debug ())                         \
         crab::outs() << DEBUG_LOC << "\n";               \
      crab::outs() << "Property : " << PROP << "\n";      \
      crab::outs() << "Invariant: " << INV << "\n";       \
      crab::outs() << " -----------------------------\n"; \
    }                                                     \
  } while (false);  

  #define LOG_ERR(VERBOSE,INV,PROP,DEBUG_LOC)             \
  do {                                                    \
    this->m_db.add (_ERR, DEBUG_LOC);                     \
    if (VERBOSE >=1) {                                    \
       crab::outs() << " --- ERROR -----------------\n";  \
       if (DEBUG_LOC.has_debug ())                        \
           crab::outs() << DEBUG_LOC << "\n";             \
       crab::outs() << "Property : " << PROP << "\n";     \
       crab::outs() << "Invariant: " << INV << "\n";      \
       crab::outs() << " -----------------------------\n";\
    }                                                     \
  } while (false); 


  namespace num_dom_detail {

     template<typename Domain>
     struct checker_ops {
       typedef typename Domain::varname_t varname_t;
       typedef crab::domains::interval<z_number> interval_t;
       typedef linear_constraint< z_number, varname_t> z_lin_cst_t;
       Domain& m_inv;
       checker_ops (Domain& inv): m_inv (inv) { }

       interval_t operator[](varname_t v){ return m_inv [v]; }

       // if the domain is not numerical then return always false
       bool entails(z_lin_cst_t cst) { 
         // special cases first
         if (m_inv.is_bottom()) return true;

         if (cst.is_tautology ()) return true;

         if (cst.is_contradiction ()) return false;

         Domain inv;
         inv += cst;
         if (inv.is_top ()) return false;
         return m_inv <= inv; 
       }

       // if the domain is not numerical then return always true
       bool intersect(z_lin_cst_t cst) {
         // special cases first
         if (m_inv.is_bottom () || cst.is_contradiction ()) return false;

         if (m_inv.is_top () || cst.is_tautology ()) return true;

         Domain inv;
         inv += cst;
         Domain meet = m_inv & inv;
         return (!meet.is_bottom ());
       }
     };

  } // end num_dom_detail namespace

  template<typename Analyzer>
  class property_checker: public statement_visitor <typename Analyzer::varname_t> {
   public:
    typedef typename Analyzer::abs_tr_ptr abs_tr_ptr;
    typedef typename Analyzer::varname_t varname_t;
    typedef typename Analyzer::abs_dom_t abs_dom_t;

    typedef variable < z_number, varname_t > z_var_t;
    typedef linear_expression< z_number, varname_t > z_lin_exp_t;
    typedef linear_constraint< z_number, varname_t > z_lin_cst_t;
    typedef linear_constraint_system< z_number, varname_t > z_lin_cst_sys_t;

    typedef binary_op<z_number,varname_t>         z_bin_op_t;
    typedef assignment<z_number,varname_t>        z_assign_t;
    typedef assume_stmt<z_number,varname_t>       z_assume_t;
    typedef assert_stmt<z_number,varname_t>       z_assert_t;
    typedef havoc_stmt<varname_t>                 havoc_t;
    typedef unreachable_stmt<varname_t>           unreach_t;
    typedef select_stmt <z_number,varname_t>      z_select_t;
    typedef callsite_stmt<varname_t>              callsite_t;
    typedef return_stmt<varname_t>                return_t;
    typedef array_init_stmt<varname_t>            z_arr_init_t;
    typedef assume_array_stmt<z_number,varname_t> z_assume_arr_t;
    typedef array_store_stmt<z_number,varname_t>  z_arr_store_t;
    typedef array_load_stmt<z_number,varname_t>   z_arr_load_t;
    typedef ptr_store_stmt<varname_t>             ptr_store_t;
    typedef ptr_load_stmt<varname_t>              ptr_load_t;
    typedef ptr_assign_stmt<z_number,varname_t>   z_ptr_assign_t;
    typedef ptr_object_stmt<varname_t>            ptr_object_t;
    typedef ptr_function_stmt<varname_t>          ptr_function_t;
    typedef ptr_null_stmt<varname_t>              ptr_null_t;
    typedef ptr_assume_stmt<varname_t>            ptr_assume_t;
    typedef ptr_assert_stmt<varname_t>            ptr_assert_t;

    typedef std::set<std::pair<debug_info, check_kind_t> > check_results_db;
 
   protected: 

    abs_tr_ptr m_abs_tr; // it can be null
    int m_verbose;
    checks_db m_db; // Store debug information about the checks
   
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
      
    virtual void check (ptr_store_t& s) { 
      if (!this->m_abs_tr) return;        
        s.accept (&*this->m_abs_tr); // propagate m_inv to the next stmt
    }
      
    virtual void check (ptr_load_t& s) { 
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
    void visit (ptr_store_t &s) { check (s); }
    void visit (ptr_load_t &s) { check (s); }
    void visit (z_ptr_assign_t &s) { check (s); }
    void visit (ptr_object_t &s) { check (s); }
    void visit (ptr_function_t &s) { check (s); }
    void visit (ptr_null_t &s) { check (s); }
    void visit (ptr_assert_t &s) { check (s); }
    
    property_checker (int verbose): 
        m_abs_tr (nullptr), m_verbose (verbose) { }
    
    void set (abs_tr_ptr abs_tr) {
      m_abs_tr = abs_tr;
    }
    
    const checks_db& get_db () const { return m_db; }
    
    checks_db get_db () { return m_db; }

    virtual std::string get_property_name () const {
       return "dummy property";
    }

    void write (crab_os& o) const {
      o << get_property_name () << "\n";
      m_db.write (o);
    }

  };  
  } // end namespace
} // end namespace
#endif 
