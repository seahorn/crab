#ifndef KILLGEN_FIXPOINT_ITERATOR_HPP
#define KILLGEN_FIXPOINT_ITERATOR_HPP

/* A specialized fixpoint iterator for kill-gen problems */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/cfg/CfgBgl.hpp>
#include <crab/analysis/graphs/Sccg.hpp>
#include <crab/analysis/graphs/TopoOrder.hpp>

namespace crab {

  namespace iterators {

    // A quick wrapper for a generic kill-gen domain
    template<class Element>
    class killgen_domain: public ikos::writeable {
         
     private:

      typedef killgen_domain<Element> killgen_domain_t;
      typedef discrete_domain<Element> discrete_domain_t;
      
     public:

      typedef typename discrete_domain_t::iterator iterator;
      typedef Element element_t;

     private:

      discrete_domain_t _inv;
      
     public:
      
      killgen_domain(discrete_domain_t inv)
          : ikos::writeable(), _inv(inv){ } 
      
     public:
      
      static killgen_domain_t top() {
        return killgen_domain(discrete_domain_t::top());
      }
      
      static killgen_domain_t bottom() {
        return killgen_domain(discrete_domain_t::bottom());
      }
      
     public:
      
      killgen_domain()
          : ikos::writeable(), _inv(discrete_domain_t::bottom()){ }
      
      killgen_domain(Element e)
          : ikos::writeable(), _inv(e) { }
      
      killgen_domain(const killgen_domain_t &o)
          : ikos::writeable(), _inv(o._inv) { } 
          
      killgen_domain(killgen_domain_t &&o)
          : ikos::writeable(), _inv(std::move(o._inv)) { } 
      
      killgen_domain_t& operator=(const killgen_domain_t &other) {
        if (this != &other) 
          _inv = other._inv;
        return *this;
      }
      
      killgen_domain_t& operator=(killgen_domain_t &&other) {
        _inv = std::move(other._inv);
        return *this;
      }
      
      iterator begin() { return _inv.begin(); }
      
      iterator end() { return _inv.end(); }
      
      unsigned size() { return _inv.size(); }
      
      bool is_bottom() { return _inv.is_bottom(); }
         
      bool is_top() { return _inv.is_top(); }
         
      bool operator<=(killgen_domain_t other) {
        if (is_bottom ()) 
          return true;
        else if (other.is_top ())
          return true;
        else
          return (_inv <= other._inv);
      }
         
      void operator-=(Element x) {
        if (is_bottom ()) 
          return;
        _inv -= x;
      }
         
      void operator-=(killgen_domain_t other) {
        if (is_bottom () || other.is_bottom ()) 
          return;
        
        if (!other._inv.is_top()) {
          for (auto v : other) 
            _inv -= v; 
        }
      }
      
      void operator+=(Element x) {
        if (is_top ()) 
          return;
        _inv += x;
      }
      
      void operator+=(killgen_domain_t other) {
        if (is_top () || other.is_bottom ()) {
          return;
        } else if (other.is_top ()) {
          _inv = discrete_domain_t::top();
        } else {
          _inv = (_inv | other._inv);
        }
      }
      
      killgen_domain_t operator|(killgen_domain_t other) {
        return (_inv | other._inv);
      }
      
      killgen_domain_t operator&(killgen_domain_t other) {
        return (_inv & other._inv);
      }
         
      void write(crab_os& o) { _inv.write(o); }
      
    }; 

    // API for a kill-gen analysis
    template<class CFG, class Element>
    class killgen_analysis {

     public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;    
      typedef killgen_domain<Element> killgen_domain_t;

     protected:

      CFG _cfg;

     public:

      killgen_analysis (CFG cfg): _cfg(cfg) { }

      virtual ~killgen_analysis() { }

      // whether forward or backward analysis
      virtual bool is_forward () = 0;

      // initial state
      virtual killgen_domain_t entry() = 0;
 
      // (optional) initialization for the fixpoint
      virtual void init_fixpoint () = 0;

      // confluence operator
      virtual killgen_domain_t merge(killgen_domain_t, killgen_domain_t) = 0;

      // analyze a basic block
      virtual killgen_domain_t analyze (basic_block_label_t, killgen_domain_t) = 0;

      // analysis name
      virtual std::string name () = 0;
    };

    // A simple fixpoint for a killgen analysis
    template<class CFG, class KillGenAnalysis>
    class killgen_fixpoint_iterator {

     public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename KillGenAnalysis::killgen_domain_t killgen_domain_t;
      typedef typename boost::unordered_map<basic_block_label_t,killgen_domain_t>::iterator iterator;
      typedef typename boost::unordered_map<basic_block_label_t,killgen_domain_t>::const_iterator const_iterator;
      
     protected:

      CFG _cfg;
      boost::unordered_map<basic_block_label_t, killgen_domain_t> _in_map;
      boost::unordered_map<basic_block_label_t, killgen_domain_t> _out_map;

     private:

      KillGenAnalysis _analysis;

      /// XXX: run_bwd_fixpo(G) is equivalent to run_fwd_fixpo(reverse(G)).
      ///      However, weak_rev_topo_sort(G) != weak_topo_sort(reverse(G))
      /// For instance, for a G=(V,E) where
      ///   V= {v1,v2, v3, v4, v5}, 
      ///   E= {(v1,v2), (v1,v3), (v2,v4), (v4,v1), (v3,v5)}
      /// (1) weak_rev_topo_sort(cfg)=[v5,v3,v4,v2,v1] 
      /// (2) weak_topo_sort(reverse(cfg))=[v5,v3,v2,v4,v1] or even
      ///     worse [v5,v3,v1,v4,v2] if vertices in the same scc are
      ///     traversed in preorder.
      /// For a backward analysis, (1) will converge faster.
      /// For all of this, we decide not to reverse graphs and have
      /// two dual versions for the forward and backward analyses.

      void run_fwd_fixpo (std::vector<typename CFG::node_t> &order,
                          unsigned &iterations){

        order = crab::analyzer::graph_algo::weak_topo_sort(_cfg);
        assert (order.size () == std::distance(_cfg.begin(), _cfg.end()));
        bool change = true;
        iterations = 0;
        while (change) {
          change = false;
          ++iterations;
          for (auto &n: order) {
            auto in = _analysis.entry();
            for (auto p: _cfg.prev_nodes (n))
              in = _analysis.merge(in, _out_map[p]); 
            auto old_out = _out_map[n];
            auto out = _analysis.analyze(n, in);
            if (!(out <= old_out)) {
              _out_map[n] = _analysis.merge(out, old_out);
              change = true;
            } else 
              _in_map[n] = in;
          }
        }
      }
      
      void run_bwd_fixpo (std::vector<typename CFG::node_t> &order,
                          unsigned &iterations){

        order = crab::analyzer::graph_algo::weak_rev_topo_sort(_cfg);
        assert (order.size () == std::distance(_cfg.begin(), _cfg.end()));
        bool change = true;
        iterations = 0;
        while (change) {
          change = false;
          ++iterations;
          for (auto &n: order) {
            auto out = _analysis.entry();
            for (auto p: _cfg.next_nodes (n))
              out = _analysis.merge(out, _in_map[p]); 
            auto old_in = _in_map[n];
            auto in = _analysis.analyze(n, out);
            if (!(in <= old_in)) {
              _in_map[n] = _analysis.merge(in, old_in);
              change = true;
            } else 
              _out_map[n] = out;
          }
        }
      }

     public:

      killgen_fixpoint_iterator (CFG cfg) : _cfg (cfg), _analysis (_cfg) { }

      void release_memory () {
        _in_map.clear();
        _out_map.clear();
      }

      void run() { 
        crab::ScopedCrabStats __st__(_analysis.name());

        _analysis.init_fixpoint(); 

        std::vector<typename CFG::node_t> order;
        unsigned iterations = 0;
        if (_analysis.is_forward())
          run_fwd_fixpo(order, iterations);
        else
          run_bwd_fixpo(order, iterations);

        CRAB_LOG(_analysis.name(), 
                 crab::outs()  << "fixpoint ordering={"; 
                 bool first=true;
                 for (auto &v : order) {
                   if (!first) crab::outs() << ",";
                   first=false;
                   crab::outs() << cfg_impl::get_label_str(v); 
                 }
                 crab::outs() << "}\n";); 
        

        CRAB_LOG(_analysis.name(), 
                 crab::outs() << _analysis.name() << ": " 
                              << "fixpoint reached in " << iterations << " iterations.\n"); 
        
        CRAB_LOG(_analysis.name(), 
                 crab::outs() << _analysis.name() << " sets:\n";
                 for (auto n: boost::make_iterator_range (_cfg.label_begin (),
                                                          _cfg.label_end ())) {
                   crab::outs() << cfg_impl::get_label_str(n) << " "
                                << "IN="  << _in_map[n]  << " "
                                << "OUT=" << _out_map[n] << "\n"; 
                 }
                 crab::outs() << "\n";);
      }      

      iterator in_begin() { return _in_map.begin(); }
      iterator in_end() { return _in_map.end(); } 
      const_iterator in_begin() const { return _in_map.begin(); }
      const_iterator in_end() const { return _in_map.end(); } 

      iterator out_begin() { return _out_map.begin(); } 
      iterator out_end() { return _out_map.end(); }
      const_iterator out_begin() const { return _out_map.begin(); } 
      const_iterator out_end() const { return _out_map.end(); }

   }; 

  } // end namespace iterators
} // end namespace crab

#endif 
