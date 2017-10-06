#ifndef CRAB_BACKWARD_ANALYZER_HPP
#define CRAB_BACKWARD_ANALYZER_HPP

#include <crab/cfg/cfg.hpp>
#include <crab/cfg/var_factory.hpp>
#include <crab/iterators/fwd_fixpoint_iterators.hpp>
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/analysis/abs_transformer.hpp>
#include <crab/analysis/dataflow/liveness.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/noncopyable.hpp>
#include <boost/range/iterator_range.hpp>

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

	CRAB_LOG("backward-fixpoint",
		 crab::outs () << "Post at " << cfg_impl::get_label_str(node) << ": "
		               << precond << "\n");
		 
	// invariants that hold at the entry of the block
	AbsDom invariant = m_invariants [node];
	// rebuild local invariants that hold at each program point.
	abs_fwd_tr_t F(&invariant);
	pp_abstract_map_t pp_invariants;
	for(auto &s: boost::make_iterator_range(bb.begin(),bb.end())) {
	  pp_invariants.insert (std::make_pair(&s, F.inv()));
	  CRAB_LOG("backward-fixpoint",
		   crab::outs () << "\tRebuilding at statement " << s << " inv="
		                 << F.inv() << "\n");
	  s.accept (&F);	  	  
	}
	
	CRAB_LOG("backward-fixpoint",
		 crab::outs ()
		 << "Done forward propagation at each program point \n"
		 << "Starting backward propagation ... \n");

	// compute precondition at the entry of the block
	abs_bwd_tr_t B(&precond, pp_invariants); 
	for(auto &s: boost::make_iterator_range(bb.rbegin(),bb.rend()))
	  s.accept (&B);
	CRAB_LOG("backward-fixpoint",
		 crab::outs () << "Pre at " << cfg_impl::get_label_str(node) << ": "
		               << precond << "\n");
	
      }
      
      virtual void process_pre(bb_label_t /*node*/,
			       AbsDom /*postcond*/) override { }
      
      /**
       *  Store necessary preconditions 
       **/
      virtual void process_post (bb_label_t node,
				 AbsDom precond) override {
	m_preconditions.insert(std::make_pair (node, precond));
      }
      
    public:

      typedef typename fixpoint_iterator_t::wto_t wto_t;
      typedef bb_abstract_map_t precond_map_t;     
      typedef typename precond_map_t::iterator iterator;
      typedef typename precond_map_t::const_iterator const_iterator;
      
      necessary_preconditions_fixpoint_iterator
      (CFG cfg, const wto_t *wto, AbsDom postcond,
       /* fixpoint parameters */
       unsigned int widening_delay,
       unsigned int descending_iterations,
       size_t jump_set_size)
	: fixpoint_iterator_t(cfg::cfg_rev<CFG> (cfg), wto,
			      widening_delay, descending_iterations, jump_set_size),
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
      
      const wto_t& get_WTO() const {
	return this->get_wto();
      }
      
    };

    /**
     *  A forward-backward analyzer to compute necessary
     *  preconditions based on Cousot&Cousot's JLP'92
     *
     * The API of this class is such that from outside it looks pretty
     * much like a forward analyzer.
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
      typedef boost::unordered_map<bb_label_t, AbsDom> invariant_map_t; 
      typedef typename bwd_fixpoint_iterator_t::precond_map_t precond_map_t;
      typedef typename bwd_fixpoint_iterator_t::wto_t bwd_wto_t;
      typedef liveness<CFG> liveness_t;     
      
    public:
      
      typedef typename fwd_analyzer_t::assumption_map_t assumption_map_t;
      // bwd_wto_t and wto_t are different types because bwd_wto_t is
      // over the reversed CFG.
      typedef typename fwd_analyzer_t::wto_t wto_t;

    private:
      
      CFG m_cfg;
      // We keep the two wto's (from forward and reversed CFGs) to
      // avoid recompute them during the below iterative process.
      // Only the forward wto is exposed to outside clients.
      const wto_t* m_wto;
      const bwd_wto_t* m_b_wto;      
      
      // to keep the results of the last iteration
      invariant_map_t m_pre_invariants;
      invariant_map_t m_post_invariants;
      precond_map_t m_preconditions;
      
      // if a checker uses this analyzer to prove a property it will
      // need an abstract transformer to rebuild local invariants.
      abs_tr_ptr m_checker_abs_tr;

      void store_forward_analysis_results (fwd_analyzer_t &f) {
	m_pre_invariants = invariant_map_t(f.pre_begin(), f.pre_end());
	m_post_invariants = invariant_map_t(f.post_begin(), f.post_end());	    
      }
      
      void store_analysis_results (fwd_analyzer_t &f, bwd_fixpoint_iterator_t &b) {
	store_forward_analysis_results(f);
	m_preconditions = precond_map_t(b.begin(), b.end());
      }
      
    public:

      typedef typename bwd_fixpoint_iterator_t::iterator iterator;
      typedef typename bwd_fixpoint_iterator_t::const_iterator const_iterator;

      intra_forward_backward_analyzer (CFG cfg)
	: m_cfg(cfg),
	  m_wto(nullptr),
	  m_b_wto(nullptr),
	  m_checker_abs_tr(new intra_abs_transformer<AbsDom> (nullptr))
      {}

      ~intra_forward_backward_analyzer () {
	if (m_wto) delete m_wto;
	if (m_b_wto) delete m_b_wto;
      }
      
      /**
       * Perform the refining forward-backward loop.
       **/
      void run (AbsDom init_states, AbsDom final_states,
		// behaves as a standard forward analysis
		bool only_forward,
		// assumptions
		assumption_map_t &assumptions,
		// liveness information
		const liveness_t* live, 
		// parameters for each forward or backward analysis
		unsigned int widening_delay=1,
		unsigned int descending_iters=UINT_MAX,
		size_t jump_set_size=0) {

	CRAB_LOG ("backward",
		  crab::outs() << "Initial states=" << init_states << "\n";
		  crab::outs() << "Final states=" << final_states << "\n");

	crab::CrabStats::count ("CombinedForwardBackward.invocations");	
	unsigned iters = 0;
	while (true) {
	  iters++;
          crab::CrabStats::count ("CombinedForwardBackward.iterations");
	  CRAB_VERBOSE_IF(1, crab::outs() << "Iteration " << iters << "\n" 	  
			                  << "Started forward analysis.\n";);

	  crab::CrabStats::resume ("CombinedForwardBackward.ForwardPass");
	  // run forward analysis computing invariants
	  fwd_analyzer_t F (m_cfg, m_wto, init_states, live,
			    widening_delay, descending_iters, jump_set_size);
	  F.run (assumptions);
	  crab::CrabStats::stop ("CombinedForwardBackward.ForwardPass");	  
	  
	  // reuse wto for next iteration
	  if (iters == 1) m_wto = new wto_t(F.get_wto ());

	  CRAB_VERBOSE_IF(1, crab::outs () << "Finished forward analysis.\n";);
	  
	  CRAB_LOG("backward",
		   for (auto &kv: boost::make_iterator_range (F.pre_begin(),
							      F.pre_end ())) {
		     crab::outs () << cfg_impl::get_label_str (kv.first)
				   << ":\n" << kv.second << "\n";
		   });

	  if (only_forward) {
	    store_forward_analysis_results(F);
	    CRAB_VERBOSE_IF(1, crab::outs () << "\nSkipped backward pass.\n";);
	    break;
	  }
	  
	  CRAB_VERBOSE_IF(1, crab::outs () << "Started backward analysis.\n";);		   

	  crab::CrabStats::resume ("CombinedForwardBackward.BackwardPass");	  
	  // run backward analysis computing necessary preconditions
	  // refined with invariants
	  bwd_fixpoint_iterator_t B (m_cfg, m_b_wto, final_states,
				     widening_delay, descending_iters, jump_set_size);
	  B.Run (boost::make_iterator_range(F.pre_begin(), F.pre_end()));
	  crab::CrabStats::stop ("CombinedForwardBackward.BackwardPass");	  

	  CRAB_VERBOSE_IF(1, crab::outs () << "Finished backward analysis.\n";);
	  
	  CRAB_LOG("backward",
		   for (auto &kv: boost::make_iterator_range (B.begin(),
							      B.end ())) {
		     crab::outs () << cfg_impl::get_label_str (kv.first)
				   << ":\n" << kv.second << "\n";
		   }
		   crab::outs () << "\n");
	  
	  AbsDom new_init_states = B[m_cfg.entry ()];
	  if (init_states <= new_init_states) {
            //CRAB_LOG ("backward",
	    //           crab::outs() << "Cannot refine more initial states.\n");
	    store_analysis_results(F,B);
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

	  // reuse wto for next iteration
	  if (iters == 1) m_b_wto = new bwd_wto_t(B.get_WTO ());
	  
	} // end while true

	CRAB_VERBOSE_IF(1, crab::outs() << "Combined forward+backward analysis done after "
			                << iters << " iterations.\n";); 	
      }

      // Return the invariants that hold at the entry of b
      AbsDom operator[] (bb_label_t b) const {
        return get_pre (b);
      }
      
      // Return the invariants that hold at the entry of b
      AbsDom get_pre (bb_label_t b) const { 
        auto it = m_pre_invariants.find (b);
        if (it == m_pre_invariants.end ())
          return abs_dom_t::top ();
        else
          return it->second;
      }
      
      // Return the invariants that hold at the exit of b
      AbsDom get_post (bb_label_t b) const {
        auto it = m_post_invariants.find (b);
        if (it == m_post_invariants.end ())
          return abs_dom_t::top ();
        else 
          return it->second;      
      }

      // Return the necessary preconditions for basic block b
      AbsDom get_preconditions(bb_label_t b) const {
	auto it = m_preconditions.find (b);
	if (it != m_preconditions.end())
	  return it->second;
	else
	  return AbsDom::top();
      }

      // Return the wto of the cfg
      const wto_t& get_wto () const {
	assert (m_wto);
	return *m_wto;
      }
      
      // API for checkers
      CFG get_cfg (void) { return m_cfg; }

      // API for checkers
      abs_tr_ptr get_abs_transformer (AbsDom &inv) {
	m_checker_abs_tr->set (inv);
	return m_checker_abs_tr;
      }
      
	      
    };

  } // end namespace
} // end namespace

#endif
