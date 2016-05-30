#ifndef INTER_ANALYSIS_DATASTRUCTURES_HPP
#define INTER_ANALYSIS_DATASTRUCTURES_HPP

#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <crab/common/types.hpp>
#include <crab/cfg/cfg.hpp>

/* Datastructures needed during inter-procedural analysis */

namespace crab {

  namespace analyzer {

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

      // Assume context-insensitive analysis so it will merge all
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
          m_call_table.insert (make_pair (func_key, inv));
      }

     public:
      
      call_ctx_table() { }

      void insert (callsite_t cs, AbsDomain inv) {
        insert_helper (cfg_hasher<CFG>::hash (cs), inv);
      }

      void insert (fdecl_t d, AbsDomain inv) {
        insert_helper (cfg_hasher<CFG>::hash (d), inv);
      }

      AbsDomain get_call_ctx (fdecl_t d) const {
        auto it = m_call_table.find (cfg_hasher<CFG>::hash (d));
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
        // --- Summary involving only m_params + m_out variables
        abs_domain_t m_sum;
        // --- Keep the returned value of the function if any
        boost::optional <varname_t> m_ret;
        // --- Keep all the formal parameters of the function. They can
        //     be INT_TYPE or PTR_TYPE inputs as well as input/output
        //     ARR_TYPE parameters.
        std::vector <varname_t> m_params;
        
       public:
        
        Summary (fdecl_t fdecl,
                 abs_domain_t sum, 
                 boost::optional<varname_t> ret,
                 const vector<varname_t> &params): 
            m_fdecl (fdecl), 
            m_sum (sum), 
            m_ret (ret), 
            m_params(params)  { }
        
        fdecl_t get_fdecl () const { return m_fdecl; }
        
        abs_domain_t get_sum () const { return m_sum;}
        
        const vector<varname_t>& get_params () const { return m_params;}

        boost::optional<varname_t> get_ret_val () const { return m_ret;}

        void write(crab_os &o)  {
          o << m_fdecl << " --> " << " variables = {";
          for (auto const &p: m_params) {
            o << p << ";";
          }
          if (m_ret) o << *m_ret;
          o << "}" << " summary = " << m_sum;
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
                   boost::optional<varname_t> ret,
                   const std::vector<varname_t>& params) {

        std::vector<varname_t> ps (params.begin(), params.end ());
        summary_ptr sum_tuple (new Summary (d, sum, ret, ps));
        m_sum_table.insert (std::make_pair (cfg_hasher<CFG>::hash (d), sum_tuple));
      }

      // return true if there is a summary
      bool hasSummary (callsite_t cs) const {
        auto it = m_sum_table.find (cfg_hasher<CFG>::hash (cs));
        return (it != m_sum_table.end ());
      }

      bool hasSummary (fdecl_t d) const {
        auto it = m_sum_table.find (cfg_hasher<CFG>::hash (d));
        return (it != m_sum_table.end ());
      }

      // get the summary
      Summary& get (callsite_t cs) const {
        auto it = m_sum_table.find (cfg_hasher<CFG>::hash (cs));
        assert (it != m_sum_table.end ());
        
        return *(it->second);
      }

      Summary& get (fdecl_t d) const {
        auto it = m_sum_table.find (cfg_hasher<CFG>::hash (d));
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
