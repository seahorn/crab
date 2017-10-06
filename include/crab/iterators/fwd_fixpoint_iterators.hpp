/*******************************************************************************
 *
 * Forward fixpoint iterators of varying complexity and precision.
 *
 * The interleaved fixpoint iterator is described in G. Amato and F. Scozzari's
 * paper: Localizing widening and narrowing. In Proceedings of SAS 2013,
 * pages 25-42. LNCS 7935, 2013.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#ifndef IKOS_FWD_FIXPOINT_ITERATORS_HPP
#define IKOS_FWD_FIXPOINT_ITERATORS_HPP

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/iterators/wto.hpp>
#include <crab/iterators/fixpoint_iterators_api.hpp>
#include <crab/iterators/thresholds.hpp>

namespace ikos {

  namespace interleaved_fwd_fixpoint_iterator_impl {
    
    template< typename NodeName, typename CFG, typename AbstractValue >
    class wto_iterator;

    template< typename NodeName, typename CFG, typename AbstractValue >
    class wto_processor;
    
  } // namespace interleaved_fwd_fixpoint_iterator_impl
  
  template< typename NodeName, typename CFG, typename AbstractValue >
  class interleaved_fwd_fixpoint_iterator: 
      public forward_fixpoint_iterator< NodeName, CFG, AbstractValue > {

    friend class interleaved_fwd_fixpoint_iterator_impl::wto_iterator< NodeName, CFG, AbstractValue >;

  public:
    typedef wto< NodeName, CFG> wto_t;

  private:
    typedef boost::unordered_map<NodeName,AbstractValue> invariant_table_t;
    typedef boost::shared_ptr<invariant_table_t> invariant_table_ptr;
    typedef interleaved_fwd_fixpoint_iterator_impl::wto_iterator< NodeName, CFG, AbstractValue > wto_iterator_t;
    typedef interleaved_fwd_fixpoint_iterator_impl::wto_processor< NodeName, CFG, AbstractValue > wto_processor_t;
    typedef crab::iterators::thresholds_t thresholds_t;
    
    CFG _cfg;
    wto_t _wto;
    invariant_table_ptr _pre, _post;
    // number of iterations until triggering widening
    unsigned int _widening_delay;
    // number of narrowing iterations. If the narrowing operator is
    // indeed a narrowing operator this parameter is not
    // needed. However, there are abstract domains for which an actual
    // narrowing operation is not available so we must enforce
    // termination.
    unsigned int _descending_iterations;
    // whether jump set is used for widening
    bool _use_widening_jump_set;    
    // set of thresholds to jump during widening
    thresholds_t _jump_set;

    void set(invariant_table_ptr table, NodeName node, const AbstractValue& v) {
      std::pair< typename invariant_table_t::iterator, bool > res = 
          table->insert(std::make_pair(node, v));
      if (!res.second) {
        (res.first)->second = v;
      }
    }
    
    void set_pre(NodeName node, const AbstractValue& v) {
      this->set(this->_pre, node, v);
    }

    void set_post(NodeName node, const AbstractValue& v) {
      this->set(this->_post, node, v);
    }

    AbstractValue get(invariant_table_ptr table, NodeName n) {
      typename invariant_table_t::iterator it = table->find(n);
      if (it != table->end()) {
        return it->second;
      } else {
        return AbstractValue::bottom();
      }
    }

    AbstractValue extrapolate(NodeName node, unsigned int iteration, 
                              AbstractValue before, AbstractValue after) {

      CRAB_VERBOSE_IF(2, crab::outs() << "Increasing iteration=" << iteration << "\n"
		      << "Widening at " << crab::cfg_impl::get_label_str(node) << "\n";);      

      if (iteration <= _widening_delay) {
        CRAB_LOG("fixpo",
                 auto widen_res = before | after;
                 crab::outs() << "Prev   : " << before << "\n"
                           << "Current: " << after << "\n"
                           << "Res    : " << widen_res << "\n");
        return before | after; 
      } else {
        CRAB_LOG("fixpo",
                 crab::outs() << "Prev   : " << before << "\n"
                           << "Current: " << after << "\n");
        
        if (_use_widening_jump_set) {
          CRAB_LOG("fixpo",
                   auto widen_res = before.widening_thresholds (after, _jump_set);
                   crab::outs() << "Res    : " << widen_res << "\n");
          return before.widening_thresholds (after, _jump_set);
        } else {
          CRAB_LOG("fixpo",
                   auto widen_res = before || after;
                   crab::outs() << "Res    : " << widen_res << "\n");
          return before || after;
        }
      }
    }

    AbstractValue refine(NodeName node, unsigned int iteration, 
                         AbstractValue before, AbstractValue after) {

      CRAB_VERBOSE_IF(2, 
		      crab::outs() << "Decreasing iteration=" << iteration << "\n"
		      << "Narrowing at " << crab::cfg_impl::get_label_str(node) << "\n";);

      if (iteration == 1) {
        CRAB_LOG("fixpo",
                 auto narrow_res = before && after;
                 crab::outs() << "Prev   : " << before << "\n"
                           << "Current: " << after << "\n"
                           << "Res    : " << narrow_res << "\n");
        return before & after; 
      } else {
        CRAB_LOG("fixpo",
                 auto narrow_res = before && after;
                 crab::outs() << "Prev   : " << before << "\n"
                 << "Current: " << after << "\n"
                           << "Res    : " << narrow_res << "\n");
        return before && after; 
      }
    }

    void initialize_thresholds(size_t jump_set_size) {
      if (_use_widening_jump_set) {
        crab::CrabStats::resume ("Fixpo");
        // select statically some widening points to jump to.
        _jump_set = _cfg.initialize_thresholds_for_widening(jump_set_size);
        crab::CrabStats::stop ("Fixpo");
      }      
    }
    
  public:

    typedef invariant_table_t assumption_map_t;
    
    interleaved_fwd_fixpoint_iterator(CFG cfg, const wto_t *wto,
                                      unsigned int widening_delay,
                                      unsigned int descending_iterations,
                                      size_t jump_set_size): 
        _cfg(cfg),
        _wto(!wto ? cfg: *wto),
        _pre(boost::make_shared<invariant_table_t>()),
        _post(boost::make_shared<invariant_table_t>()),
        _widening_delay(widening_delay),
        _descending_iterations(descending_iterations),
        _use_widening_jump_set (jump_set_size > 0) {
      initialize_thresholds(jump_set_size);
    }
    
    virtual ~interleaved_fwd_fixpoint_iterator() { }
    
    CFG get_cfg() const {
      return this->_cfg;
    }

    const wto_t& get_wto() const {
      return this->_wto;
    }

    AbstractValue get_pre(NodeName node) {
      return this->get(this->_pre, node);
    }
    
    AbstractValue get_post(NodeName node) {
      return this->get(this->_post, node);
    }

    void run(AbstractValue init) {
      crab::ScopedCrabStats __st__("Fixpo");
      CRAB_VERBOSE_IF(2, crab::outs() << "== Started fixpoint\n");
      this->set_pre(this->_cfg.entry(), init);
      wto_iterator_t iterator(this);
      this->_wto.accept(&iterator);
      wto_processor_t processor(this);
      this->_wto.accept(&processor);
      CRAB_VERBOSE_IF(2, crab::outs() << "== Fixpoint reached.\n");
      reset ();
    }

    void run(AbstractValue init, assumption_map_t &assumptions) {
      crab::ScopedCrabStats __st__("Fixpo");
      CRAB_VERBOSE_IF(2, crab::outs() << "== Started fixpoint\n");      
      this->set_pre(this->_cfg.entry(), init);
      wto_iterator_t iterator(this, &assumptions);
      this->_wto.accept(&iterator);
      wto_processor_t processor(this);
      this->_wto.accept(&processor);
      CRAB_VERBOSE_IF(2, crab::outs() << "== Fixpoint reached.\n");      
    }

    void reset () {
      this->_pre.reset();
      this->_post.reset();      
    }
        
  }; // class interleaved_fwd_fixpoint_iterator

  namespace interleaved_fwd_fixpoint_iterator_impl {
    
    template< typename NodeName, typename CFG, typename AbstractValue >
    class wto_iterator: public wto_component_visitor<NodeName, CFG> {
      
    public:
      typedef interleaved_fwd_fixpoint_iterator<NodeName, CFG, AbstractValue> interleaved_iterator_t;
      typedef wto_vertex<NodeName, CFG> wto_vertex_t;
      typedef wto_cycle<NodeName, CFG> wto_cycle_t;
      typedef wto<NodeName, CFG> wto_t;
      typedef typename wto_t::wto_nesting_t wto_nesting_t;
      typedef typename interleaved_iterator_t::assumption_map_t assumption_map_t;
      
    private:
      interleaved_iterator_t *_iterator;
      assumption_map_t *_assumptions;

      inline AbstractValue strengthen (NodeName n, AbstractValue inv) {
	if (_assumptions) {
	  auto it = _assumptions->find(n);
	  if (it != _assumptions->end ()) {
	    CRAB_LOG ("fixpo",
	    	      crab::outs() << "Before assumption at " << n << ":"
		                   << inv << "\n");
		                   
	    inv = inv & it->second;
	    CRAB_LOG ("fixpo",
	    	      crab::outs() << "After assumption at " << n << ":"
	    	                   << inv << "\n");
	  }
	}
	return inv;
      }
		  
    public:
      wto_iterator(interleaved_iterator_t *iterator):
	_iterator(iterator), _assumptions (nullptr) { }

      wto_iterator(interleaved_iterator_t *iterator,
		   assumption_map_t *assumptions):
	_iterator(iterator), _assumptions (assumptions) { }      
      
      void visit(wto_vertex_t& vertex) {
        AbstractValue pre;
        NodeName node = vertex.node();
        if (node == this->_iterator->get_cfg().entry()) {
          pre = this->_iterator->get_pre(node);
        } else {
          auto prev_nodes = this->_iterator->_cfg.prev_nodes(node);
          pre = AbstractValue::bottom();
          CRAB_VERBOSE_IF (2,
		           crab::outs() << "Joining predecessors of "
			   << crab::cfg_impl::get_label_str(node) << "\n");
          for (NodeName prev : prev_nodes) {
            pre |= this->_iterator->get_post(prev);  
          }
	  pre = strengthen (node, pre);	  
          this->_iterator->set_pre(node, pre);
        }
	
        CRAB_VERBOSE_IF (2, crab::outs() << "Analyzing node "
			 << crab::cfg_impl::get_label_str(node) << "\n");
        this->_iterator->analyze(node, pre);
        this->_iterator->set_post(node, pre);
      }
      
      void visit(wto_cycle_t& cycle) {
        NodeName head = cycle.head();
        CRAB_VERBOSE_IF (2,
			 crab::outs() << "** Analyzing loop with head "
			 << crab::cfg_impl::get_label_str(head) << "\n");
	
        wto_nesting_t cycle_nesting = this->_iterator->_wto.nesting(head);
        auto prev_nodes = this->_iterator->_cfg.prev_nodes(head);
        AbstractValue pre = AbstractValue::bottom();
	CRAB_VERBOSE_IF (2, crab::outs() << "Joining predecessors of "
			 << crab::cfg_impl::get_label_str(head) << "\n");
        for (NodeName prev : prev_nodes) {
          if (!(this->_iterator->_wto.nesting(prev) > cycle_nesting)) {
            pre |= this->_iterator->get_post(prev); 
          }
        }
	pre = strengthen (head, pre);
	
        for(unsigned int iteration = 1; ; ++iteration) {
	  // keep track of how many times the cycle is visited by the fixpoint
	  cycle.increment_fixpo_visits ();
	  
          // Increasing iteration sequence with widening
          this->_iterator->set_pre(head, pre);
          AbstractValue post(pre); 
          CRAB_VERBOSE_IF(2, crab::outs() << "Analyzing node "
			  << crab::cfg_impl::get_label_str(head) << "\n");
          this->_iterator->analyze(head, post);
          this->_iterator->set_post(head, post);
          for (typename wto_cycle_t::iterator it = cycle.begin();
	       it != cycle.end(); ++it) {
            it->accept(this);
          }
          AbstractValue new_pre = AbstractValue::bottom();
          for (NodeName prev : prev_nodes) {
            new_pre |= this->_iterator->get_post(prev); 
          }
          if (new_pre <= pre) {
            // Post-fixpoint reached
            CRAB_VERBOSE_IF(2, crab::outs() << "post-fixpoint reached\n");
            this->_iterator->set_pre(head, new_pre);
            pre = new_pre;
            break;
          } else {
            pre = this->_iterator->extrapolate(head, iteration, pre, new_pre);
          }
        }

	
	if (this->_iterator->_descending_iterations == 0) {
	  // no narrowing
	  return;
	}

	CRAB_VERBOSE_IF(2, crab::outs () << "Started narrowing phase\n";);
	
        for(unsigned int iteration = 1; ; ++iteration) {
          // Decreasing iteration sequence with narrowing
          AbstractValue post(pre); 
          CRAB_VERBOSE_IF(2,crab::outs() << "Analyzing node "
			  << crab::cfg_impl::get_label_str(head) << "\n");
          this->_iterator->analyze(head, post);
          this->_iterator->set_post(head, post);
          for (typename wto_cycle_t::iterator it = cycle.begin();
	       it != cycle.end(); ++it) {
            it->accept(this);
          }
          AbstractValue new_pre = AbstractValue::bottom();
          for (NodeName prev : prev_nodes) {
            new_pre |= this->_iterator->get_post(prev); 
          }
          if (pre <= new_pre) {
            CRAB_VERBOSE_IF(2, crab::outs() << "No more refinement possible.\n");
            // No more refinement possible (pre == new_pre)
            break;
          } else {
            if (iteration > this->_iterator->_descending_iterations) break; 
            pre = this->_iterator->refine(head, iteration, pre, new_pre);
            this->_iterator->set_pre(head, pre);
          }
        }
        CRAB_VERBOSE_IF (2,
			 crab::outs() << "** Finished loop with head "
			 << crab::cfg_impl::get_label_str(head) << "\n");
	
      }
      
    }; // class wto_iterator
  
    template< typename NodeName, typename CFG, typename AbstractValue >
    class wto_processor: public wto_component_visitor< NodeName, CFG > {

    public:
      typedef interleaved_fwd_fixpoint_iterator< NodeName, CFG, AbstractValue > interleaved_iterator_t;
      typedef wto_vertex< NodeName, CFG > wto_vertex_t;
      typedef wto_cycle< NodeName, CFG > wto_cycle_t;

    private:
      interleaved_iterator_t *_iterator;
      
    public:
      wto_processor(interleaved_iterator_t *iterator): _iterator(iterator) { }
      
      void visit(wto_vertex_t& vertex) {
        NodeName node = vertex.node();
        this->_iterator->process_pre(node, this->_iterator->get_pre(node));
        this->_iterator->process_post(node, this->_iterator->get_post(node));
      }
      
      void visit(wto_cycle_t& cycle) {
        NodeName head = cycle.head();
        this->_iterator->process_pre(head, this->_iterator->get_pre(head));
        this->_iterator->process_post(head, this->_iterator->get_post(head));
        for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end(); ++it) {
          it->accept(this);
        }	
      }
      
    }; // class wto_processor
  
  } // interleaved_fwd_fixpoint_iterator_impl  
} // namespace ikos
#endif // IKOS_FWD_FIXPOINT_ITERATORS
