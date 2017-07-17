#ifndef CRAB_BACKWARD_ANALYZER_HPP
#define CRAB_BACKWARD_ANALYZER_HPP

#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/iterators/fwd_fixpoint_iterators.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/abs_transformer.hpp>

namespace crab {
  
  namespace analyzer {

    /** 
     * Compute necessary preconditions by computing the abstract
     * least fixpoint of the reversed CFG. At each step the
     * precondition is refined with an invariant computed by a
     * forward analysis.
     **/      
    template <typename NodeName, typename CFG, typename AbsDom>
    class necessary_preconditions_fixpoint_iterator:
      public ikos::interleaved_fwd_fixpoint_iterator
                   <NodeName, cfg::cfg_rev<CFG>, AbsDom> {
      
      typedef ikos::interleaved_fwd_fixpoint_iterator
      <NodeName, cfg::cfg_rev<CFG>, AbsDom> fixpoint_iterator_t;
      
      typedef intra_abs_transformer<AbsDom> abs_fwd_tr_t;	
      typedef typename CFG::basic_block_label_t bb_label_t;
      typedef typename CFG::statement_t stmt_t; 
      typedef boost::unordered_map<bb_label_t, AbsDom> bb_abstract_map_t;
      typedef boost::unordered_map<stmt_t*, AbsDom> pp_abstract_map_t;      

      CFG m_cfg;
      // Postcondition (i.e, final states) that we want to propagate backwards
      AbsDom m_postcond;
      // necessary preconditions
      bb_abstract_map_t m_preconditions;
      // invariants
      bb_abstract_map_t m_invariants;	
      
      /**
       * Compute necessary preconditions for a basic block
       **/
      virtual void analyze(bb_label_t node, AbsDom &precond) override {
	typedef intra_necessary_preconditions_abs_transformer
	  <AbsDom, pp_abstract_map_t> abs_bwd_tr_t;
	
	auto &bb = m_cfg.get_node (node);

	CRAB_LOG("backward",
		 crab::outs () << "Post at " << node << ": "
		               << precond << "\n");
		 
	// invariants that hold at the entry of the block
	AbsDom invariant = m_invariants [node];
	// rebuild local invariants that hold at each program point.
	abs_fwd_tr_t F(&invariant);
	pp_abstract_map_t pp_invariants;
	for(auto &s: boost::make_iterator_range(bb.begin(),bb.end())) {
	  pp_invariants.insert (std::make_pair(&s, F.inv()));
	  s.accept (&F);
	}
	
	CRAB_LOG("backward",
		 crab::outs ()
		 << "Done forward propagation at each program point \n"
		 << "Starting backward propagation ... \n");

	// compute precondition at the entry of the block
	abs_bwd_tr_t B(&precond, pp_invariants); 
	for(auto &s: boost::make_iterator_range(bb.rbegin(),bb.rend()))
	  s.accept (&B);
	CRAB_LOG("backward",
		 crab::outs () << "Pre at " << node << ": "
		               << precond << "\n");
	
      }
      
      virtual void process_pre(bb_label_t /*node*/,
			       AbsDom /*postcond*/) override { }
      
      /**
       *  Store necessary preconditions 
       **/
      virtual void process_post (bb_label_t node,
				 AbsDom precond) override {
	m_preconditions.insert(make_pair (node, precond));
      }
      
    public:

      typedef bb_abstract_map_t precond_map_t;     
      typedef typename precond_map_t::iterator iterator;
      typedef typename precond_map_t::const_iterator const_iterator;
      
      necessary_preconditions_fixpoint_iterator
      (CFG cfg, AbsDom postcond,
       /* fixpoint parameters */
       unsigned int widening_delay,
       unsigned int descending_iterations,
       size_t jump_set_size)
	: fixpoint_iterator_t(cfg::cfg_rev<CFG> (cfg),
			      widening_delay, descending_iterations,
			      jump_set_size),
	  m_cfg (cfg),
	  m_postcond (postcond) { }
      
      template <typename Range>
      void Run (Range invariants) {
	m_invariants = bb_abstract_map_t(invariants.begin(), invariants.end());
	this->run (m_postcond);
      }
      
      iterator begin () { return m_preconditions.begin(); } 
      
      iterator end () { return m_preconditions.end();   }
      
      const_iterator begin () const { return m_preconditions.begin(); }
      
      const_iterator end () const { return m_preconditions.end();   }
      
      AbsDom operator[](bb_label_t node) const {
	auto it = m_preconditions.find (node);
	if (it != m_preconditions.end())
	  return it->second;
	else
	  return AbsDom::top();
      }
    };

    /**
     *  A forward-backward analyzer to compute necessary
     *  preconditions based on Cousot&Cousot's JLP'92
     **/
    template<typename CFG, typename AbsDom>
    class intra_forward_backward_analyzer: public boost::noncopyable {
    public:

      typedef CFG cfg_t;
      typedef typename CFG::basic_block_label_t bb_label_t;
      typedef typename CFG::varname_t varname_t;
      typedef typename CFG::number_t number_t;
      // used for checkers
      typedef AbsDom abs_dom_t;
      typedef boost::shared_ptr<intra_abs_transformer<AbsDom> > abs_tr_ptr;
      
    private:
      
      typedef intra_fwd_analyzer<CFG, AbsDom> fwd_analyzer_t;
      typedef necessary_preconditions_fixpoint_iterator<bb_label_t, CFG, AbsDom>
      bwd_fixpoint_iterator_t;
      typedef typename bwd_fixpoint_iterator_t::precond_map_t precond_map_t;
      
      CFG m_cfg;
      precond_map_t m_preconditions;
      // XXX: if a checker uses this analyzer to prove a property it
      // will need an abstract transformer to rebuild local invariants.
      abs_tr_ptr m_checker_abs_tr;
      
    public:
      
      typedef typename bwd_fixpoint_iterator_t::iterator iterator;
      typedef typename bwd_fixpoint_iterator_t::const_iterator const_iterator;

      intra_forward_backward_analyzer (CFG cfg)
	: m_cfg(cfg),
	  m_checker_abs_tr(new intra_abs_transformer<AbsDom> (nullptr))
      {}
	
      /**
       * Perform the refining forward-backward loop.
       **/
      void run (AbsDom init_states, AbsDom final_states,
		// parameters for each forward or backward analysis
		unsigned int widening_delay=1,
		unsigned int descending_iters=UINT_MAX,
		size_t jump_set_size=0) {

	CRAB_LOG ("backward",
		  crab::outs() << "Initial states=" << init_states << "\n";
		  crab::outs() << "Final states=" << final_states << "\n");

	unsigned iters = 0;
	while (true) {
	  iters++;

	  CRAB_LOG ("backward",
		    crab::outs() << "Iteration " << iters << "\n"
		                 << "Starting forward pass ...\n");
	  
	  // run forward analysis computing invariants
	  fwd_analyzer_t F (m_cfg, init_states, nullptr,
			    widening_delay, descending_iters, jump_set_size);
	  F.run ();

	  CRAB_LOG("backward",
		   crab::outs () << "Forward analysis done:\n";
		   for (auto &kv: boost::make_iterator_range (F.pre_begin(),
							      F.pre_end ())) {
		     crab::outs () << cfg_impl::get_label_str (kv.first)
				   << ":\n" << kv.second << "\n";
		   }
		   crab::outs () << "\nStarting backward pass ...\n";);
	  
	  // run backward analysis computing necessary preconditions
	  // refined with invariants
	  bwd_fixpoint_iterator_t B (m_cfg, final_states,
				     widening_delay, descending_iters,
				     jump_set_size);
	  B.Run (boost::make_iterator_range(F.pre_begin(), F.pre_end()));

	  CRAB_LOG("backward",
		   crab::outs () << "Backward analysis done:\n";
		   for (auto &kv: boost::make_iterator_range (B.begin(),
							      B.end ())) {
		     crab::outs () << cfg_impl::get_label_str (kv.first)
				   << ":\n" << kv.second << "\n";
		   }
		   crab::outs () << "\n");
	  
	  AbsDom new_init_states = B[m_cfg.entry ()];
	  if (init_states <= new_init_states) {
            CRAB_LOG ("backward",
		      crab::outs() << "Cannot refine more initial states.\n");
	    // store results for future queries
	    m_preconditions = precond_map_t (B.begin(), B.end());
	    break;
	  } else {
            CRAB_LOG ("backward",
		      crab::outs() << "Initial state refined: \n"
		                   << init_states << " ==> ");
	    
	    if (iters > descending_iters) {
	      // ensure termination by using narrowing
	      init_states = init_states && new_init_states;
	    } else {
	      init_states = init_states & new_init_states;
	    }

            CRAB_LOG ("backward",
		      crab::outs() << init_states << "\n");
	    
	  }
	} // end while true

	CRAB_LOG ("backward",
		  crab::outs() << "Done after " << iters
		               << " iterations.\n");
      }
      
      AbsDom operator[](bb_label_t b) const {
	auto it = m_preconditions.find (b);
	if (it != m_preconditions.end())
	  return it->second;
	else
	  return AbsDom::top();
      }

      // used for checkers
      CFG get_cfg (void) { return m_cfg; }

      // used for checkers
      abs_tr_ptr get_abs_transformer (AbsDom &inv) {
	m_checker_abs_tr->set (inv);
	return m_checker_abs_tr;
      }
      
	      
    };

  } // end namespace
} // end namespace

#endif
