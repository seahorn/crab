#ifndef LIVENESS_ANALYSIS_HPP
#define LIVENESS_ANALYSIS_HPP

/* Liveness analysis */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/types.hpp>
#include <crab/iterators/killgen_fixpoint_iterator.hpp>

#include <boost/noncopyable.hpp>

namespace crab {

  namespace analyzer {

   template<typename V>
   using liveness_domain = crab::domains::flat_killgen_domain<V>;       

    
   template<class CFG>
   class liveness_ops:
      public crab::iterators::
      killgen_operations_api<CFG, liveness_domain<typename CFG::varname_t> > {

   public:
     
     typedef liveness_domain<typename CFG::varname_t> liveness_domain_t;
     typedef typename CFG::basic_block_label_t basic_block_label_t;

   private:
     
     typedef crab::iterators::killgen_operations_api<CFG, liveness_domain_t> this_type;
     typedef std::pair<liveness_domain_t, liveness_domain_t> binding_t;
     typedef boost::unordered_map<basic_block_label_t, binding_t> liveness_map_t;
     
     liveness_map_t _liveness_map;

    public:

     liveness_ops (CFG cfg): this_type (cfg) { }

     virtual bool is_forward () { return false; }

     virtual liveness_domain_t entry() 
     { return liveness_domain_t::bottom();}

     virtual liveness_domain_t merge(liveness_domain_t d1,  liveness_domain_t d2)
     { return d1 | d2; }

     virtual void init_fixpoint() {
        for (auto &b: boost::make_iterator_range(this->_cfg.begin(),this->_cfg.end())) {
          liveness_domain_t kill, gen;
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
          _liveness_map.insert(std::make_pair(b.label(), binding_t(kill, gen)));
        }
      }

      virtual liveness_domain_t analyze (basic_block_label_t bb_id, liveness_domain_t in) {
        auto it = _liveness_map.find (bb_id);
        assert(it != _liveness_map.end ());
        in -= it->second.first;
        in += it->second.second;
        return in;
      }

      virtual std::string name () { return "liveness";}
   };

   //! Live variable analysis
   template<typename CFG>
   class liveness:
      boost::noncopyable, 
      public crab::iterators::killgen_fixpoint_iterator <CFG, liveness_ops<CFG> >{
     
    public:

     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::statement_t statement_t;
     typedef typename CFG::varname_t varname_t;
     
    private:

     typedef liveness_ops<CFG> liveness_ops_t;
     typedef typename liveness_ops_t::liveness_domain_t liveness_domain_t;
     typedef crab::iterators::killgen_fixpoint_iterator<CFG,liveness_ops_t> 
             killgen_fixpoint_iterator_t;
     
    public:

     //for backward compatibility
     //XXX: maybe unused already
     typedef liveness_domain_t set_t;

    private:

     // output of the analysis: map basic blocks to set of dead
     // variables at the end of the blocks
     boost::unordered_map<basic_block_label_t, liveness_domain_t> _dead_map;

     // statistics 
     unsigned _max_live;
     unsigned _total_live;
     unsigned _total_blks;

     void process_post (basic_block_label_t bb, 
                        liveness_domain_t live_out) {
       // --- Collect dead variables at the exit of bb
       if (!live_out.is_bottom ()) {
         liveness_domain_t dead_set = this->_cfg.get_node(bb).live();
         dead_set -= live_out;
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << dead_set <<"\n";);
         _dead_map.insert (std::make_pair(bb, dead_set));
         // update statistics
         _total_live += live_out.size ();
         _max_live = std::max (_max_live, live_out.size ());
         _total_blks ++;
       } else {
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << live_out <<"\n";);
       }
     }
          
    public:

     liveness (CFG cfg)
         : killgen_fixpoint_iterator_t(cfg), 
           _max_live (0), _total_live (0), _total_blks (0) { }

     void exec() { 
       this->run();
       for (auto p : boost::make_iterator_range(this->out_begin(), this->out_end()))
       { process_post (p.first, p.second); } 
       this->release_memory();
     }

     //! return the set of dead variables at the exit of block bb
     liveness_domain_t dead_exit (basic_block_label_t bb) const {
       auto it = _dead_map.find(bb);
       if (it == _dead_map.end()) 
         return liveness_domain_t();
       else 
         return it->second; 
     }
       
     void get_stats (unsigned& total_live, 
                     unsigned& max_live_per_blk,
                     unsigned& avg_live_per_blk)  const {
       total_live = _total_live;
       max_live_per_blk = _max_live;
       avg_live_per_blk = (_total_blks == 0 ? 0 : (int) _total_live / _total_blks);
     }

     // TODO
     void write (crab_os &o) const { }         
   }; 

   template <typename CFG>
   crab_os& operator << (crab_os& o, const liveness<CFG> &l) {
     l.write (o);
     return o;
   }

  } // end namespace analyzer
} // end namespace crab
#endif 
