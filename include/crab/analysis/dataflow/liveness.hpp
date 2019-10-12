#pragma once 

/* Liveness analysis */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/killgen_domain.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>

#include <unordered_map>
#include <boost/range/iterator_range.hpp>

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
     typedef std::unordered_map<basic_block_label_t, binding_t> liveness_map_t;
     
     liveness_map_t m_liveness_map;

    public:

     liveness_analysis_operations(CFG cfg): parent_type(cfg) { }

     virtual bool is_forward() { return false; }

     virtual varset_domain_t entry() {
       varset_domain_t res = varset_domain_t::bottom();
       if (this->m_cfg.has_func_decl()) {
	 auto fdecl = this->m_cfg.get_func_decl();
	 for(unsigned i=0, e= fdecl.get_num_outputs(); i<e; ++i) {
	   res += fdecl.get_output_name(i);
	 }
       }
      return res;
     }

     virtual varset_domain_t merge(varset_domain_t d1,  varset_domain_t d2) {
       return d1 | d2;
     }

     virtual void init_fixpoint() {
       for (auto &b: boost::make_iterator_range(this->m_cfg.begin(),this->m_cfg.end())){
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
      public crab::iterators::
      killgen_fixpoint_iterator<CFG,liveness_analysis_operations<CFG>>{

      typedef liveness_analysis_operations<CFG> liveness_analysis_operations_t;
      typedef crab::iterators::
      killgen_fixpoint_iterator<CFG,liveness_analysis_operations_t> killgen_fixpoint_iterator_t; 

      liveness_analysis(const liveness_analysis<CFG>& other) = delete;
      liveness_analysis<CFG>& operator=(const liveness_analysis<CFG>& other) = delete;
      
    public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::statement_t statement_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename liveness_analysis_operations_t::varset_domain_t varset_domain_t;
      
    private:
      
      // output of the analysis: map basic blocks to set of live
      // variables at the end of the blocks
      std::unordered_map<basic_block_label_t, varset_domain_t> _live_map;
      
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
    inline crab_os& operator<<(crab_os& o, const liveness_analysis<CFG> &l) {
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
    class liveness {
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
      CFG m_cfg;
      
      // output of the analysis: map basic blocks to set of dead
      // variables at the end of the blocks
      std::unordered_map<basic_block_label_t, varset_domain_t> _dead_map;
      
      // statistics 
      unsigned m_max_live;
      unsigned m_total_live;
      unsigned m_total_blocks;
      
    public:
      
      //for backward compatibility
      //XXX: maybe unused already
      typedef varset_domain_t set_t;
      
      liveness(CFG cfg)
	: m_cfg(cfg)
	, m_max_live(0)
	, m_total_live(0)
	, m_total_blocks (0) { }

      liveness(const liveness<CFG>& other) = delete;
      
      liveness<CFG>& operator=(const liveness<CFG>& other) = delete;
          
      void exec() {
	/** Remove dead variables locally **/
	
	liveness_analysis_t live(m_cfg);
	live.exec();
	for (auto &bb: boost::make_iterator_range(m_cfg.begin(), m_cfg.end())) {
	  varset_domain_t live_set = live.get(bb.label());
	  if (live_set.is_bottom()) continue;

	  varset_domain_t dead_set = m_cfg.get_node(bb.label()).live();
	  // dead variables = (USE(bb) U DEF(bb)) \ live_out(bb)
	  dead_set -= live_set;
	  CRAB_LOG("liveness",
		   crab::outs() << cfg_impl::get_label_str(bb.label()) 
		                << " dead variables=" << dead_set <<"\n";);
	  _dead_map.insert (std::make_pair(bb.label(), dead_set));
	  // update statistics
	  m_total_live += live_set.size ();
	  m_max_live = std::max (m_max_live, live_set.size ());
	  m_total_blocks ++;
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
	total_live = m_total_live;
	max_live_per_blk = m_max_live;
	avg_live_per_blk = (m_total_blocks == 0 ? 0 : (int) m_total_live / m_total_blocks);
      }
      
      void write (crab_os &o) const {
	o << "TODO: printing dead variable analysis results\n";
      }         
    }; 
    
   template<typename CFG>
   inline crab_os& operator<<(crab_os& o, const liveness<CFG> &l) {
     l.write(o);
     return o;
   }

  } // end namespace analyzer
} // end namespace crab
