#pragma once 

/* Liveness analysis */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/types.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>

#include <boost/noncopyable.hpp>

namespace crab {

  namespace analyzer {

   template<typename V>
   using varset_domain = crab::domains::flat_killgen_domain<V>;       

    /** 
     * Define the main operations for the liveness variable analysis:
     * compute for each basic block the set of live variables, i.e.,
     * variables that might be used in the future 
     **/    
   template<class CFG>
   class liveness_analysis_operations:
      public crab::iterators::
      killgen_operations_api<CFG, varset_domain<typename CFG::variable_t>> {

   public:
     
     typedef varset_domain<typename CFG::variable_t> varset_domain_t;
     typedef typename CFG::basic_block_label_t basic_block_label_t;

   private:
     
     typedef crab::iterators::killgen_operations_api<CFG, varset_domain_t> parent_type;
     typedef std::pair<varset_domain_t, varset_domain_t> binding_t;
     typedef boost::unordered_map<basic_block_label_t, binding_t> liveness_map_t;
     
     liveness_map_t m_liveness_map;

    public:

     liveness_analysis_operations(CFG cfg): parent_type(cfg) { }

     virtual bool is_forward() { return false; }

     virtual varset_domain_t entry() {
      return varset_domain_t::bottom();
     }

     virtual varset_domain_t merge(varset_domain_t d1,  varset_domain_t d2) {
       return d1 | d2;
     }

     virtual void init_fixpoint() {
       for (auto &b: boost::make_iterator_range(this->_cfg.begin(),this->_cfg.end())){
	 varset_domain_t kill, gen;
	 for (auto &s: boost::make_iterator_range(b.rbegin(),b.rend())) { 
	   auto live = s.get_live();
	   for (auto d: boost::make_iterator_range(live.defs_begin(), 
						   live.defs_end())) {
	     kill += d; 
	     gen -= d;
	   }
	   for (auto u: boost::make_iterator_range(live.uses_begin(), 
						   live.uses_end())) {
	     gen  += u; 
	   }
	 }
	 m_liveness_map.insert(std::make_pair(b.label(), binding_t(kill, gen)));
       }
     }
     
      virtual varset_domain_t analyze(basic_block_label_t bb_id,
				      varset_domain_t in) {
        auto it = m_liveness_map.find(bb_id);
        assert(it != m_liveness_map.end());
        in -= it->second.first;
        in += it->second.second;
        return in;
      }

      virtual std::string name() { return "liveness";}
   };

    /** Live variable analysis **/
    template<typename CFG>
    class liveness_analysis:
      boost::noncopyable, 
      public crab::iterators::
      killgen_fixpoint_iterator<CFG,liveness_analysis_operations<CFG>>{

      typedef liveness_analysis_operations<CFG> liveness_analysis_operations_t;
      typedef crab::iterators::
      killgen_fixpoint_iterator<CFG,liveness_analysis_operations_t> killgen_fixpoint_iterator_t; 
      
    public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::statement_t statement_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename liveness_analysis_operations_t::varset_domain_t varset_domain_t;
      
    private:
      
      // output of the analysis: map basic blocks to set of live
      // variables at the end of the blocks
      boost::unordered_map<basic_block_label_t, varset_domain_t> _live_map;
      
    public:
      
      liveness_analysis(CFG cfg)
	: killgen_fixpoint_iterator_t(cfg)
      { }
      
      void exec() { 
	this->run();
	for (auto p : boost::make_iterator_range(this->out_begin(), this->out_end())) {
	  _live_map.insert({p.first, p.second});
	  CRAB_LOG("liveness-live",
		   crab::outs() << cfg_impl::get_label_str(p.first) 
		                << " live variables=" << p.second <<"\n";);
	  
	} 
	this->release_memory();
      }
      
      varset_domain_t get(basic_block_label_t bb) const {
	auto it = _live_map.find(bb);
	if (it != _live_map.end()) {
	  return it->second;
	} else {
	  return varset_domain_t::bottom();
	}
      }
      
      void write (crab_os &o) const {
	o << "TODO: print liveness analysis results\n";
      }
      
    }; 
  
    template<typename CFG>
    crab_os& operator<<(crab_os& o, const liveness_analysis<CFG> &l) {
      l.write (o);
      return o;
    }

    /** 
     * Dead variable analysis.  
     * 
     * FIXME: the name "liveness" is not that great so we should
     *        change it at some point. For that we need to change crab
     *        clients.
     **/
    template<typename CFG>
    class liveness: boost::noncopyable {
    public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::basic_block_t basic_block_t;    
      typedef typename CFG::statement_t statement_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename CFG::variable_t variable_t;    
      typedef varset_domain<variable_t> varset_domain_t;
      
    private:
      
      typedef liveness_analysis<CFG> liveness_analysis_t;
      
      // the cfg
      CFG _cfg;
      
      // output of the analysis: map basic blocks to set of dead
      // variables at the end of the blocks
      boost::unordered_map<basic_block_label_t, varset_domain_t> _dead_map;
      
      // statistics 
      unsigned _max_live;
      unsigned _total_live;
      unsigned _total_blks;
      
    public:
      
      //for backward compatibility
      //XXX: maybe unused already
      typedef varset_domain_t set_t;
      
      liveness(CFG cfg)
	: _cfg(cfg)
	, _max_live(0)
	, _total_live(0)
	, _total_blks (0) { }
    
    
      void exec() {
	/** Remove dead variables locally **/
	
	liveness_analysis_t live(_cfg);
	live.exec();
	for (auto &bb: boost::make_iterator_range(_cfg.begin(), _cfg.end())) {
	  varset_domain_t live_set = live.get(bb.label());
	  if (live_set.is_bottom()) continue;

	  varset_domain_t dead_set = this->_cfg.get_node(bb.label()).live();
	  // dead variables = (USE(bb) U DEF(bb)) \ live_out(bb)
	  dead_set -= live_set;
	  CRAB_LOG("liveness",
		   crab::outs() << cfg_impl::get_label_str(bb.label()) 
		                << " dead variables=" << dead_set <<"\n";);
	  _dead_map.insert (std::make_pair(bb.label(), dead_set));
	  // update statistics
	  _total_live += live_set.size ();
	  _max_live = std::max (_max_live, live_set.size ());
	  _total_blks ++;
	}      
      }
    
      // Return the set of dead variables at the exit of block bb
      varset_domain_t dead_exit (basic_block_label_t bb) const {
	auto it = _dead_map.find(bb);
	if (it == _dead_map.end()) {
	  return varset_domain_t();
	} else {
	  return it->second;
	}
      }
    
      void get_stats(unsigned& total_live, 
		     unsigned& max_live_per_blk,
		     unsigned& avg_live_per_blk)  const {
	total_live = _total_live;
	max_live_per_blk = _max_live;
	avg_live_per_blk = (_total_blks == 0 ? 0 : (int) _total_live / _total_blks);
      }
      
      void write (crab_os &o) const {
	o << "TODO: printing dead variable analysis results\n";
      }         
    }; 
    
   template<typename CFG>
   crab_os& operator<<(crab_os& o, const liveness<CFG> &l) {
     l.write(o);
     return o;
   }

  } // end namespace analyzer
} // end namespace crab
