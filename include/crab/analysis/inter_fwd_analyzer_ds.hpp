#ifndef INTER_FWD_ANALYSIS_DATASTRUCTURES_HPP
#define INTER_FWD_ANALYSIS_DATASTRUCTURES_HPP

#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <crab/common/types.hpp>
#include <crab/cfg/cfg.hpp>

/* 
 * Datastructures needed during inter-procedural analysis
 */

namespace crab {

  namespace analyzer {

    // XXX: Currently we assume a _contex-insensitive_ approach so we
    // store a single summary per function. For implementing a
    // context-sensitive approach we would need to allow multiple
    // summaries per function.

    /* Store the calling contexts of each function */  
    template <typename CFG, typename AbsDomain>
    class call_ctx_table: boost::noncopyable {
     public:

      typedef typename CFG::basic_block_t::callsite_t callsite_t;
      typedef typename CFG::fdecl_t fdecl_t;
      typedef typename CFG::varname_t varname_t;
      typedef AbsDomain abs_domain_t;

     private:
      typedef boost::unordered_map <std::size_t, AbsDomain> call_table_t;
      
      call_table_t m_call_table;

      // XXX: assume context-insensitive analysis so it will merge all
      // calling contexts using abstract domain's join keeping a
      // single calling context per function.
      void insert_helper (std::size_t func_key, AbsDomain inv) {
        auto it = m_call_table.find (func_key);
        if (it != m_call_table.end ()) {
          //crab::CrabStats::count (AbsDomain::getDomainName() + ".count.join");
          //crab::ScopedCrabStats __st__(AbsDomain::getDomainName() + ".join");
          it->second = it->second | inv;
        }
        else
          m_call_table.insert (std::make_pair (func_key, inv));
      }

     public:
      
      call_ctx_table() { }

      void insert (callsite_t cs, AbsDomain inv) {
        insert_helper (cfg::cfg_hasher<CFG>::hash (cs), inv);
      }

      void insert (fdecl_t d, AbsDomain inv) {
        insert_helper (cfg::cfg_hasher<CFG>::hash (d), inv);
      }

      AbsDomain get_call_ctx (fdecl_t d) const {
        auto it = m_call_table.find (cfg::cfg_hasher<CFG>::hash (d));
        if (it != m_call_table.end ())
          return it->second;
        else 
          return AbsDomain::top ();
      }

    };

    /* Store the summaries for each function*/
    template<typename CFG, typename AbsDomain>
    class summary_table: boost::noncopyable {

     public:
      
      typedef typename CFG::basic_block_t::callsite_t callsite_t;
      typedef typename CFG::fdecl_t fdecl_t;
      typedef AbsDomain abs_domain_t;
      typedef typename CFG::varname_t varname_t;

      // A summary is an input-output relationship between function
      // parameters. The relationship can be as expressive as
      // AbsDomain.
      class Summary {

        // --- function info
        fdecl_t m_fdecl;
        // --- Summary involving only m_params + m_ret_vals variables
        abs_domain_t m_sum;
        // --- Keep all the input parameters of the function
        std::vector <varname_t> m_inputs;
        // --- Keep all the output parameters of the function
        std::vector <varname_t> m_outputs;
        
       public:
        
        Summary (fdecl_t fdecl,
                 abs_domain_t sum, 
                 const std::vector<varname_t> &inputs,
                 const std::vector<varname_t> &outputs):
            m_fdecl (fdecl), 
            m_sum (sum), 
            m_inputs (inputs), m_outputs (outputs) { 
            
          if (m_fdecl.get_num_inputs() != m_inputs.size ())
            CRAB_ERROR ("Mismatch between function declaration and summary parameters");
          if (m_fdecl.get_num_outputs() != m_outputs.size ())
            CRAB_ERROR ("Mismatch between function declaration and summary return vals");
        }
        
        fdecl_t get_fdecl () const { return m_fdecl; }
        
        abs_domain_t get_sum () const { return m_sum;}
        
        const std::vector<varname_t>& get_inputs () const { return m_inputs;}

        const std::vector<varname_t>& get_outputs () const { return m_outputs;}

        // Check type consistency between function declaration and callsite
        // XXXX: this is needed because we don't type check crab
        // programs.
        void check_type_consistency (callsite_t cs) const {

          if (m_fdecl.get_num_inputs() != cs.get_num_args()) 
            CRAB_ERROR ("Mismatch between number of callsite and summary parameters");
          if (m_fdecl.get_num_outputs() != cs.get_lhs().size())
            CRAB_ERROR ("Mismatch between number of callsite and summary return values");
         
          for (unsigned i=0; i < cs.get_num_args (); i++){
            if (m_fdecl.get_input_type(i) != cs.get_arg_type (i))
              CRAB_ERROR ("Mismatch between type of callsite and summary parameter");
          }

          for (unsigned i=0; i < cs.get_lhs().size(); i++){
            if (m_fdecl.get_output_type(i) != cs.get_lhs()[i].second)
              CRAB_ERROR ("Mismatch between type of callsite and summary return value");
          }
        }

        void write(crab_os &o) const {
          o << m_fdecl << " --> " << " variables = {";
          for (auto const &in: m_inputs) 
            o << in << ";";
          for (auto const &out: m_outputs) 
            o << out << ";";
          o << "}";
          abs_domain_t tmp (m_sum);
          o << " summary = " << tmp;
        }

      };

     private:

      typedef boost::shared_ptr <Summary> summary_ptr;
      typedef boost::unordered_map <std::size_t, summary_ptr> summary_table_t;
      
      summary_table_t m_sum_table;
      
     public:

      summary_table () { }
      
      // insert summary information
      void insert (fdecl_t d, 
                   AbsDomain sum,
                   const std::vector<varname_t>& inputs,
                   const std::vector<varname_t>& outputs) {

        std::vector<varname_t> ins(inputs.begin(), inputs.end ());
        std::vector<varname_t> outs(outputs.begin(), outputs.end ());
        summary_ptr sum_tuple (new Summary (d, sum, ins, outs));
        m_sum_table.insert (std::make_pair (cfg::cfg_hasher<CFG>::hash (d), sum_tuple));
      }

      // return true if there is a summary
      bool hasSummary (callsite_t cs) const {
        auto it = m_sum_table.find (cfg::cfg_hasher<CFG>::hash (cs));
        return (it != m_sum_table.end ());
      }

      bool hasSummary (fdecl_t d) const {
        auto it = m_sum_table.find (cfg::cfg_hasher<CFG>::hash (d));
        return (it != m_sum_table.end ());
      }

      // get the summary
      Summary& get (callsite_t cs) const {
        auto it = m_sum_table.find (cfg::cfg_hasher<CFG>::hash (cs));
        assert (it != m_sum_table.end ());
        
        return *(it->second);
      }

      Summary& get (fdecl_t d) const {
        auto it = m_sum_table.find (cfg::cfg_hasher<CFG>::hash (d));
        assert (it != m_sum_table.end ());
        
        return *(it->second);
      }

      void write (crab_os &o) const {
        o << "--- Begin summary table: \n";
        for (auto const &p: m_sum_table) {
          p.second->write (o);
          o << "\n";
        }
        o << "--- End summary table\n";
      }
      
      friend crab_os& operator<<(crab_os& o, const summary_table<CFG, AbsDomain>&t) {
        t.write (o);
        return o;
      }

    };

  } // end namespace
} // end namespace

#endif 
