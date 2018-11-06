#pragma once 

#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/iterators/fwd_fixpoint_iterators.hpp>
#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/inter_fwd_analyzer_ds.hpp>
#include <crab/analysis/dataflow/liveness.hpp>
#include <crab/domains/domain_traits.hpp>

#include "boost/range/algorithm/set_algorithm.hpp"

namespace crab {

  namespace analyzer {

    // template<typename CFG>
    // inline std::vector<typename CFG::varname_t> find_return_vars (const CFG& cfg)
    // {
    //   typedef typename CFG::varname_t varname_t;
    //   typedef typename CFG::number_t number_t;      
    //   std::vector<varname_t> res;

    //   if (cfg.has_exit ()) {
    //     auto const &bb = cfg.get_node (cfg.exit ());
    //     for (auto const &s : boost::make_iterator_range (bb.begin(), bb.end())) {
    //       if (s.is_return ()) {                
    //         auto ret_stmt =
    // 	      static_cast<const crab::cfg::return_stmt<number_t,varname_t> *> (&s);
    //         auto const &ret_typed_vars = ret_stmt->get_ret_vals ();
    //         res.reserve (ret_typed_vars.size ());
    //         for (auto vt: ret_typed_vars)
    //           res.push_back (vt.first);
    //         return res;
    //       }
    //     }
    //   }
    //   return res;
    // }

    /**
     * Only for internal use.
     * Implementation of an intra-procedural forward analysis.
     * 
     * Perform a standard forward flow-sensitive analysis. AbsTr
     * defines the abstract transfer functions as well as which
     * operations are modelled.
     **/
    template< typename CFG, typename AbsTr>
    class fwd_analyzer: 
        public ikos::interleaved_fwd_fixpoint_iterator<typename CFG::basic_block_label_t, 
						       CFG, 
						       typename AbsTr::abs_dom_t>,
        public boost::noncopyable {

     public:

      typedef CFG cfg_t;
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename AbsTr::abs_dom_t abs_dom_t;
      typedef AbsTr* abs_tr_ptr;
      
     private:

      typedef ikos::interleaved_fwd_fixpoint_iterator<basic_block_label_t, CFG, abs_dom_t>
      fwd_iterator_t;
      
     public:

      typedef typename fwd_iterator_t::assumption_map_t assumption_map_t;
      typedef liveness<CFG> liveness_t;     
      typedef typename fwd_iterator_t::wto_t wto_t;
      typedef typename fwd_iterator_t::iterator iterator;
      typedef typename fwd_iterator_t::const_iterator const_iterator;

     private:
      
      typedef typename liveness_t::set_t live_set_t;     

      abs_tr_ptr m_abs_tr; // the abstract transformer
      const liveness_t* m_live;
      live_set_t m_formals;
      
      void prune_dead_variables (abs_dom_t &inv, basic_block_label_t node) {
        if (!m_live) return;

	crab::ScopedCrabStats __st__("Pruning dead variables");
	
        if (inv.is_bottom() || inv.is_top()) return;
        auto dead = m_live->dead_exit (node);       

        dead -= m_formals;
        domains::domain_traits<abs_dom_t>::forget(inv, dead.begin(), dead.end()); 
      }

      //! Given a basic block and the invariant at the entry it produces
      //! the invariant at the exit of the block.
      void analyze (basic_block_label_t node, abs_dom_t &inv) {
        auto &b = this->get_cfg().get_node (node);
	// XXX: set takes a reference to inv so no copies here
	m_abs_tr->set (inv);
        for (auto &s : b) { s.accept (m_abs_tr); }
        prune_dead_variables (inv, node);
      } 
      
      void process_pre (basic_block_label_t node, abs_dom_t inv) {}
      void process_post (basic_block_label_t node, abs_dom_t inv) {}
      
     public:

      fwd_analyzer (CFG cfg, const wto_t *wto, abs_tr_ptr abs_tr,
		    // fixpoint parameters
                    unsigned int widening_delay,
                    unsigned int descending_iters,
                    size_t jump_set_size,
		    // live can be nullptr if no live info is available
		    const liveness_t* live)
	: fwd_iterator_t(cfg, wto,
			 widening_delay, descending_iters, jump_set_size,
			 false /*disable processor*/), 
	  m_abs_tr(abs_tr),
	  m_live(live) {
	CRAB_VERBOSE_IF(1, crab::outs() << "Type checking CFG ... ";);
	crab::CrabStats::resume("CFG type checking");
	crab::cfg::type_checker<CFG> tc(this->_cfg);
	tc.run();
	crab::CrabStats::stop("CFG type checking");	
	CRAB_VERBOSE_IF(1, crab::outs() << "OK\n";);
	
        if (live) {
          // --- collect input and output parameters 
          if (auto fdecl = this->get_cfg ().get_func_decl ()) {
	    for (unsigned i=0; i < (*fdecl).get_num_inputs();i++)
	      m_formals += (*fdecl).get_input_name (i);
	    for (unsigned i=0; i < (*fdecl).get_num_outputs();i++)
	      m_formals += (*fdecl).get_output_name (i);
	  }
        }

      }
      
      iterator       pre_begin ()       { return this->_pre.begin(); } 
      iterator       pre_end ()         { return this->_pre.end();   }
      const_iterator pre_begin () const { return this->_pre.begin(); }
      const_iterator pre_end ()   const { return this->_pre.end();   }
      
      iterator       post_begin ()       { return this->_post.begin(); } 
      iterator       post_end ()         { return this->_post.end();   }
      const_iterator post_begin () const { return this->_post.begin(); }
      const_iterator post_end ()   const { return this->_post.end();   }
      
      //! Trigger the fixpoint computation 
      void Run ()  {
        // initialization of static data
        domains::domain_traits<abs_dom_t>::do_initialization (this->get_cfg());
        // XXX: inv was created before the static data is initialized
        //      so it won't contain that data.
        this->run(m_abs_tr->inv());         
      }      

      void Run (basic_block_label_t entry, assumption_map_t &assumptions)  {
        // initialization of static data
        domains::domain_traits<abs_dom_t>::do_initialization (this->get_cfg());
        // XXX: inv was created before the static data is initialized
        //      so it won't contain that data.	
        this->run(entry, m_abs_tr->inv(), assumptions);         
      }      
      
      //! Propagate inv through statements
      abs_tr_ptr get_abs_transformer (abs_dom_t &inv) {
	/// XXX: set takes a reference to inv so no copies here
	m_abs_tr->set (inv);
	return m_abs_tr;
      }

      //! Return the invariants that hold at the entry of b
      inline abs_dom_t operator[](basic_block_label_t b) const {
        return get_pre(b);
      }
      
      //! Return the invariants that hold at the entry of b
      abs_dom_t get_pre(basic_block_label_t b) const {
        auto it = this->_pre.find (b);
        if (it == this->_pre.end ()) {
          return abs_dom_t::bottom();
	  // if the basic block is not in the invariant table it must
	  // be because it was not reached by the analysis. We
	  // returned top but it never had real effect because
	  // process_pre made sure that all unreachable blocks were in
	  // the invariant table with a bottom invariant. This was
	  // just a waste of space.
 	  // 
	  // return abs_dom_t::top ();
	} else {
          return it->second;
	}
      }
      
      //! Return the invariants that hold at the exit of b
      abs_dom_t get_post(basic_block_label_t b) const {
        auto it = this->_post.find (b);
        if (it == this->_post.end ()) {
          return abs_dom_t::bottom();	  
          //return abs_dom_t::top ();
	} else {
          return it->second;
	}
      }

      //! Return the WTO of the CFG. The WTO contains also how many
      //! times each head was visited by the fixpoint iterator.
      const wto_t& get_WTO () const {
	return this->get_wto ();
      }

      // clear all invariants (pre and post)
      void clear() {
	this->_pre.clear();
	this->_post.clear();
      }
      
    }; 

    /**
     * Wrapper for fwd_analyzer_class. The main difference with
     * fwd_analyzer class is that here we create an abstract
     * transformer instance while fwd_analyzer does not.
     **/
    template<typename CFG, typename AbsDomain, typename AbsTr>
    class intra_fwd_analyzer_wrapper
    {
      typedef fwd_analyzer<CFG, AbsTr> fwd_analyzer_t;

      AbsDomain m_init;      
      AbsTr m_abs_tr;
      fwd_analyzer_t m_analyzer;
      
    public:

      typedef AbsDomain abs_dom_t;
      typedef liveness<CFG> liveness_t;
      
      typedef CFG cfg_t;
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename CFG::number_t number_t;
      typedef typename fwd_analyzer_t::abs_tr_ptr abs_tr_ptr;
      typedef typename fwd_analyzer_t::wto_t wto_t;      
      
    public:

      typedef typename fwd_analyzer_t::assumption_map_t assumption_map_t;
      typedef typename fwd_analyzer_t::iterator iterator;
      typedef typename fwd_analyzer_t::const_iterator const_iterator;

    public:

      // simplest API to trigger the intra-procedural forward analysis
      intra_fwd_analyzer_wrapper(CFG cfg, 
				 // fixpoint parameters
				 unsigned int widening_delay=1,
				 unsigned int descending_iters=UINT_MAX,
				 size_t jump_set_size=0):
	m_init(AbsDomain::top()),
	m_abs_tr(&m_init),
	m_analyzer(cfg, nullptr, &m_abs_tr, 
		   widening_delay, descending_iters, jump_set_size, nullptr) { }
      
      intra_fwd_analyzer_wrapper (CFG cfg, AbsDomain init,
				  // liveness info
				  const liveness_t* live = nullptr,
				  // fixpoint parameters
				  unsigned int widening_delay=1,
				  unsigned int descending_iters=UINT_MAX,
				  size_t jump_set_size=0):
	m_init(init),
	m_abs_tr(&m_init),
	m_analyzer(cfg, nullptr, &m_abs_tr, 
		   widening_delay, descending_iters, jump_set_size,
		   live) { }

      intra_fwd_analyzer_wrapper (CFG cfg,
				  // avoid precompute wto if already available
				  // it can be null				  
				  const wto_t *wto, 
				  AbsDomain init,
				  // liveness info
				  const liveness_t* live = nullptr,
				  // fixpoint parameters
				  unsigned int widening_delay=1,
				  unsigned int descending_iters=UINT_MAX,
				  size_t jump_set_size=0):
	m_init(init),
	m_abs_tr(&m_init),
	m_analyzer(cfg, wto, &m_abs_tr, 
		   widening_delay, descending_iters, jump_set_size,
		   live) { }
      
      
      
      iterator       pre_begin ()       { return m_analyzer.pre_begin();} 
      iterator       pre_end ()         { return m_analyzer.pre_end();}
      const_iterator pre_begin () const { return m_analyzer.pre_begin();}
      const_iterator pre_end ()   const { return m_analyzer.pre_end();}
      
      iterator       post_begin ()       { return m_analyzer.post_begin();}
      iterator       post_end ()         { return m_analyzer.post_end();}
      const_iterator post_begin () const { return m_analyzer.post_begin();}
      const_iterator post_end ()   const { return m_analyzer.post_end();}

      void run() { m_analyzer.Run();}

      void run(basic_block_label_t entry, assumption_map_t &assumptions)
      { m_analyzer.Run(entry, assumptions);}
      
      abs_dom_t operator[] (basic_block_label_t b) const
      { return m_analyzer[b]; }
      
      abs_dom_t get_pre (basic_block_label_t b) const
      { return m_analyzer.get_pre (b); }
      
      abs_dom_t get_post (basic_block_label_t b) const
      { return m_analyzer.get_post (b); }

      void clear() { m_analyzer.clear(); }
      
      CFG get_cfg ()
      { return m_analyzer.get_cfg (); }
      
      abs_tr_ptr get_abs_transformer (abs_dom_t &inv)
      { return m_analyzer.get_abs_transformer (inv); }

      const wto_t& get_wto () const {
	return m_analyzer.get_WTO ();
      }

    };


    /**
     * External api
     **/
    template<typename CFG, typename AbsDomain>
    using intra_fwd_analyzer =
      intra_fwd_analyzer_wrapper<CFG, AbsDomain,
				 intra_abs_transformer<AbsDomain> >;

    
  } // end namespace
} // end namespace

