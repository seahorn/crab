/*******************************************************************************
 *
 * Construction and management of weak topological orderings (WTOs).
 *
 * The construction of weak topological orderings is based on F. Bourdoncle's
 * paper: "Efficient chaotic iteration strategies with widenings", Formal
 * Methods in Programming and Their Applications, 1993, pages 128-141.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 * Contributors: Jorge A. Navas (jorge.navas@sri.com)
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

#pragma once 

#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/container/slist.hpp>
#include <boost/unordered_map.hpp>

#include <crab/common/types.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/domains/interval.hpp>


// Define RECURSIVE_WTO to use older, recursive version.  It is
// retained for a while for performance comparison purposes.
// #define RECURSIVE_WTO

namespace ikos {

  template<typename NodeName, typename CFG>
  class wto;

  template<typename NodeName, typename CFG>
  class wto_vertex;

  template<typename NodeName, typename CFG>
  class wto_cycle;

  template<typename NodeName, typename CFG>
  class wto_component_visitor;

  template<typename NodeName, typename CFG>
  class wto_nesting {
    
    friend class wto<NodeName, CFG>;
    friend class wto_vertex<NodeName, CFG>;
    friend class wto_cycle<NodeName, CFG>;
    
  public:
    typedef wto_nesting<NodeName, CFG> wto_nesting_t;
    
  private:
    typedef std::vector<NodeName> node_list_t;
    typedef boost::shared_ptr<node_list_t> node_list_ptr;

    node_list_ptr _nodes;
    
  public:
    typedef typename node_list_t::iterator iterator;
    typedef typename node_list_t::const_iterator const_iterator;
    
  private:
    wto_nesting(node_list_ptr l): _nodes(boost::make_shared<node_list_t>(*l)) { }

    int compare(wto_nesting_t& other) const {
      const_iterator this_it = this->begin(), other_it = other.begin();
      while (this_it != this->end()) {
	if (other_it == other.end()) {
	  return 1;
	} else if (*this_it == *other_it) {
	  ++this_it;
	  ++other_it;
	} else {
	  return 2; // Nestings are not comparable
	}
      }
      if (other_it == other.end()) {
	return 0;
      } else {
	return -1;
      }
    }
    
  public:
    wto_nesting(): _nodes(boost::make_shared<node_list_t>()) { }

    void operator+=(NodeName n) {
      this->_nodes = boost::make_shared<node_list_t>(*(this->_nodes));
      this->_nodes->push_back(n);
    }
    
    wto_nesting_t operator+(NodeName n) {
      wto_nesting_t res(this->_nodes);
      res._nodes->push_back(n);
      return res;
    }

    iterator begin() {
      return _nodes->begin();
    }

    iterator end() {
      return _nodes->end();
    }

    const_iterator begin() const {
      return _nodes->begin();
    }

    const_iterator end() const {
      return _nodes->end();
    }
    
    wto_nesting_t operator^(wto_nesting_t other) const {
      wto_nesting_t res;
      for (const_iterator this_it = this->begin(), other_it = other.begin(); 
           this_it != this->end() && other_it != other.end(); 
           ++this_it, ++other_it) {
        if (*this_it == *other_it) {
          res._nodes->push_back(*this_it);
        } else {
          break;
        }
      }
      return res;
    }
    
    bool operator<=(wto_nesting_t other) const {
      return this->compare(other) <= 0;
    }
    
    bool operator==(wto_nesting_t other) const {
      return this->compare(other) == 0;
    }
    
    bool operator>=(wto_nesting_t other) const {
      return this->operator<=(other, *this);
    }
    
    bool operator>(wto_nesting_t other) const {
      return this->compare(other) == 1;
    }
    
    void write(crab::crab_os& o) const {
      o << "[";
      for (const_iterator it = this->begin(); it != this->end(); ) {
	NodeName n = *it;
	o << n;
	++it;
	if (it != this->end()) {
	  o << ", ";
	}
      }
      o << "]";
    }
  }; // class nesting

  template<typename NodeName, typename CFG>
  inline crab::crab_os& operator<<(crab::crab_os& o, const wto_nesting<NodeName, CFG>& n) {
    n.write(o);
    return o;
  }

  template<typename NodeName, typename CFG>
  class wto_component {

  public:
    typedef wto_nesting<NodeName, CFG> wto_nesting_t;
    
    virtual void accept(wto_component_visitor<NodeName, CFG> *) = 0;

    virtual ~wto_component() { }

    virtual void write(crab::crab_os& os) const = 0;

  }; // class wto_component

  template<typename NodeName, typename CFG>
  inline crab::crab_os& operator<<(crab::crab_os& o, const wto_component<NodeName, CFG>& c) {
    c.write(o);
    return o;
  }

  template<typename NodeName, typename CFG>
  class wto_vertex: public wto_component<NodeName, CFG> {

    friend class wto<NodeName, CFG>;

    NodeName _node;

    wto_vertex(NodeName node): _node(node) { }

  public:
    NodeName node() {
      return this->_node;
    }

    void accept(wto_component_visitor<NodeName, CFG> *v) {
      v->visit(*this);
    }
    
    void write(crab::crab_os& o) const {
      o <<this->_node;
    }   

  }; // class wto_vertex

  template<typename NodeName, typename CFG>
  class wto_cycle: public wto_component<NodeName, CFG> {

    friend class wto<NodeName, CFG>;
    
  public:
    typedef wto_component<NodeName, CFG> wto_component_t;
    
  private:
    typedef boost::shared_ptr<wto_component_t> wto_component_ptr;
    typedef boost::container::slist<wto_component_ptr> wto_component_list_t;
    typedef boost::shared_ptr<wto_component_list_t> wto_component_list_ptr;

    NodeName _head;
    wto_component_list_ptr _wto_components;
    // number of times the wto cycle is analyzed by the fixpoint iterator
    unsigned _num_fixpo;
    
    wto_cycle(NodeName head, 
              wto_component_list_ptr wto_components): 
      _head(head), _wto_components(wto_components), _num_fixpo(0) { }
    
  public:

    typedef boost::indirect_iterator<typename wto_component_list_t::iterator> iterator;
    typedef boost::indirect_iterator<typename wto_component_list_t::const_iterator> const_iterator;
    
    NodeName head() {
      return this->_head;
    }
    
    void accept(wto_component_visitor<NodeName, CFG> *v) {
      v->visit(*this);
    }
    
    iterator begin() {
      return boost::make_indirect_iterator(_wto_components->begin());
    }
    
    iterator end() {
      return boost::make_indirect_iterator(_wto_components->end());      
    }

    const_iterator begin() const {
      return boost::make_indirect_iterator(_wto_components->begin());
    }
    
    const_iterator end() const {
      return boost::make_indirect_iterator(_wto_components->end());      
    }
    
    void increment_fixpo_visits () {
      _num_fixpo++;
    }

    unsigned get_fixpo_visits () const {
      return _num_fixpo;
    }
    
    void write(crab::crab_os& o) const {
      o << "(" << this->_head;
      if (!this->_wto_components->empty()) {
	o << " ";
	for (const_iterator it = this->begin(); it != this->end(); ) {
	  const wto_component_t& c = *it;
	  o << c;
	  ++it;
	  if (it != this->end()) {
	    o << " ";
	  }
	}
      }
      o << ")";
      if (this->_num_fixpo > 0)
	o << "^{" << this->_num_fixpo << "}";
    }
    
  }; // class wto_cycle
  
  template<typename NodeName, typename CFG>
  class wto_component_visitor {

  public:
    typedef wto_vertex<NodeName, CFG> wto_vertex_t;
    typedef wto_cycle<NodeName, CFG> wto_cycle_t;

    virtual void visit(wto_vertex_t&) = 0;
    virtual void visit(wto_cycle_t&) = 0;
    virtual ~wto_component_visitor() {}

  }; // class wto_component_visitor

  template<typename NodeName, typename CFG>
  class wto {
    
  public:
    typedef wto_nesting<NodeName, CFG> wto_nesting_t;
    typedef wto_component<NodeName, CFG> wto_component_t;
    typedef wto_vertex<NodeName, CFG> wto_vertex_t;
    typedef wto_cycle<NodeName, CFG> wto_cycle_t;
    typedef wto<NodeName, CFG> wto_t;
  
  private:
    typedef boost::shared_ptr<wto_component_t> wto_component_ptr;
    typedef boost::shared_ptr<wto_vertex_t> wto_vertex_ptr;
    typedef boost::shared_ptr<wto_cycle_t> wto_cycle_ptr;
    typedef boost::container::slist<wto_component_ptr> wto_component_list_t;
    typedef boost::shared_ptr<wto_component_list_t> wto_component_list_ptr;
    typedef bound<z_number> dfn_t;
    typedef boost::unordered_map<NodeName, dfn_t> dfn_table_t;
    typedef boost::shared_ptr<dfn_table_t> dfn_table_ptr;
    typedef std::vector<NodeName> stack_t;
    typedef boost::shared_ptr<stack_t> stack_ptr;
    typedef boost::unordered_map<NodeName, wto_nesting_t> nesting_table_t;
    typedef boost::shared_ptr<nesting_table_t> nesting_table_ptr;
    
    wto_component_list_ptr _wto_components;
    dfn_table_ptr _dfn_table;
    dfn_t _num;
    stack_ptr _stack;
    nesting_table_ptr _nesting_table;

    class nesting_builder: public wto_component_visitor<NodeName, CFG> {
      
    public:
      typedef wto_vertex<NodeName, CFG> wto_vertex_t;
      typedef wto_cycle<NodeName, CFG> wto_cycle_t;
      
    private:
      wto_nesting_t _nesting;
      nesting_table_ptr _nesting_table;
      
    public:
      nesting_builder(nesting_table_ptr nesting_table): 
	_nesting_table(nesting_table) { }

      void visit(wto_cycle_t& cycle) {
        NodeName head = cycle.head();
        wto_nesting_t previous_nesting = this->_nesting;
        this->_nesting_table->insert(std::make_pair(head, this->_nesting));
        this->_nesting += head;
        for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end(); ++it) {
          it->accept(this);
        }
        this->_nesting = previous_nesting;
      }
      
      void visit(wto_vertex_t& vertex) {
        this->_nesting_table->insert(std::make_pair(vertex.node(), this->_nesting));
      }
      
    }; // class nesting_builder

    dfn_t get_dfn(NodeName n) {
      typename dfn_table_t::iterator it = this->_dfn_table->find(n);
      if (it == this->_dfn_table->end()) {
        return 0;
      } else {
        return it->second;
      }
    }
    
    void set_dfn(NodeName n, dfn_t dfn) {
      std::pair<typename dfn_table_t::iterator, bool> res = 
	this->_dfn_table->insert(std::make_pair(n, dfn));
      if (!res.second) {
        (res.first)->second = dfn;
      }
    }
    
    NodeName pop() {
      if (this->_stack->empty()) {
        CRAB_ERROR("WTO computation: empty stack");
      } else {
        NodeName top = this->_stack->back();
        this->_stack->pop_back();
        return top;
      }
    }

    void push(NodeName n) {
      this->_stack->push_back(n);
    }

    wto_cycle_ptr component(CFG cfg, NodeName vertex) {
      auto partition = boost::make_shared<wto_component_list_t>();
      auto next_nodes = cfg.next_nodes(vertex);
      for (NodeName succ : next_nodes) {
        if (this->get_dfn(succ) == 0) {
          this->visit(cfg, succ, partition);
        }
      }
      return wto_cycle_ptr(new wto_cycle_t(vertex, partition));
    }
    
    #ifndef RECURSIVE_WTO
    struct visit_stack_elem {
      typedef typename CFG::succ_range succ_range;
      typedef typename CFG::succ_iterator succ_iterator;      
      NodeName _node;
      succ_iterator _it; // begin iterator for node's successors
      succ_iterator _et; // end iterator for node's successors
      dfn_t _min;        // smallest dfn number of any (direct or
			 // indirect) node's successor through node's
			 // DFS subtree, included node.
      
      visit_stack_elem(NodeName node, succ_range succs, dfn_t min)
	: _node(node)
	, _it(succs.begin())
	, _et(succs.end())
	, _min(min) {}
    };
    
    void visit(CFG cfg, NodeName vertex, wto_component_list_ptr partition) {
      typedef std::vector<visit_stack_elem> visit_stack_t;      
      visit_stack_t visit_stack;
      std::set<NodeName> loop_nodes;
      
      /* discover vertex */
      push(vertex);
      _num += 1;
      set_dfn(vertex, _num);
      
      visit_stack.push_back(visit_stack_elem(vertex, cfg.next_nodes(vertex), _num));
      CRAB_LOG("wto-nonrec",
	       crab::outs() << "WTO: Node " << vertex << ": dfs num=" << _num << "\n";);      
      while (!visit_stack.empty()) {
	/*
	 * Perform dfs.
	 * 
	 * When this loop terminates, visit_stack.back()_node's children
	 * have been processed.  For each loop iteration we push in
	 * visit_stack one more descendant.
	 */
	while (visit_stack.back()._it != visit_stack.back()._et) {
	  NodeName child = *visit_stack.back()._it++;
	  dfn_t child_dfn = get_dfn(child);
	  if (child_dfn == 0) {
	    /* discover new vertex */
	    push(child);
	    _num += 1;
	    set_dfn(child, _num);
	    visit_stack.push_back(visit_stack_elem(child,cfg.next_nodes(child), _num));
	    CRAB_LOG("wto-nonrec",
		     crab::outs() << "WTO: Node " << child << ": dfs num=" << _num << "\n";);
	  } else {
	    if (child_dfn <= visit_stack.back()._min) {
	      visit_stack.back()._min = child_dfn;
	      CRAB_LOG("wto-nonrec",
		       crab::outs() << "WTO: loop found " << child << "\n";);
	      loop_nodes.insert(child);
	    }
	  }
	}
	

	// propagate min from child to parent
	NodeName visiting_node = visit_stack.back()._node;
	dfn_t min_visiting_node = visit_stack.back()._min;
	bool is_loop = loop_nodes.count(visiting_node)> 0;
	visit_stack.pop_back();
	if (!visit_stack.empty() && visit_stack.back()._min > min_visiting_node) {
	  visit_stack.back()._min = min_visiting_node;
	}

	auto dfn_visiting_node = get_dfn(visiting_node);
	CRAB_LOG("wto-nonrec",
	    crab::outs() << "WTO: popped node " << visiting_node
	                 << " dfs num= " <<  dfn_visiting_node 
	                 << ": min=" << min_visiting_node << "\n";);
	
	if (min_visiting_node == get_dfn(visiting_node)) {
	  CRAB_LOG("wto-nonrec",
		   crab::outs() << "WTO: BEGIN building partition for node "
	                        << visiting_node << "\n";);
	  set_dfn(visiting_node, dfn_t::plus_infinity());
	  NodeName element = pop();
	  if (is_loop) {
	    while (!(element == visiting_node)) {
	      set_dfn(element, 0);
	      CRAB_LOG("wto-nonrec",
		       crab::outs () << "\tWTO: node " << element << ": dfn num=0\n";);
	      element = pop();
	    }
	    CRAB_LOG("wto-nonrec",
		     crab::outs() << "\tWTO: adding component starting from "
		                 << visiting_node << "\n";);
	    partition->push_front(boost::static_pointer_cast<wto_component_t,wto_cycle_t> 
				  (component(cfg, visiting_node)));
	  } else {
	    CRAB_LOG("wto-nonrec",
		     crab::outs() << "\tWTO: adding vertex " << visiting_node << "\n";);
	    partition->push_front(boost::static_pointer_cast< wto_component_t, wto_vertex_t>
				  (wto_vertex_ptr(new wto_vertex_t(visiting_node))));
	  }
	  CRAB_LOG("wto-nonrec", crab::outs() << "WTO: END building partition\n";);
	}
      } // end while (!visit_stack.empty())
    }
    
    #else
    dfn_t visit(CFG cfg, NodeName vertex, wto_component_list_ptr partition) {
      dfn_t head = 0, min = 0;
      bool loop;
      NodeName element;

      this->push(vertex);
      this->_num += 1;
      head = this->_num;
      this->set_dfn(vertex, head);
      loop = false;
      auto next_nodes = cfg.next_nodes(vertex);
      for (NodeName succ: next_nodes) {
        dfn_t succ_dfn = this->get_dfn(succ);
        if (succ_dfn == 0) {
          min = this->visit(cfg, succ, partition);
        } else {
          min = succ_dfn;
        }
        if (min <= head) {
          head = min;
	  loop = true;
        }
      }
      if (head == this->get_dfn(vertex)) {
        this->set_dfn(vertex, dfn_t::plus_infinity());
        element = this->pop();
        if (loop) {
          while (!(element == vertex)) {
            this->set_dfn(element, 0);
            element = this->pop();
          }
          partition->push_front(boost::static_pointer_cast<wto_component_t, 
                                wto_cycle_t>(this->component(cfg, vertex)));
        } else {
          partition->push_front(boost::static_pointer_cast<wto_component_t, 
                                wto_vertex_t>(wto_vertex_ptr(new wto_vertex_t(vertex))));
	}
      }
      return head;
    }
    #endif
    
    void build_nesting() {
      nesting_builder builder(this->_nesting_table);
      for (iterator it = this->begin(); it != this->end(); ++it) {
        it->accept(&builder);
      }
    }

  public:

    typedef boost::indirect_iterator<typename wto_component_list_t::iterator> iterator;
    typedef boost::indirect_iterator<typename wto_component_list_t::const_iterator> const_iterator;
    
    wto(CFG cfg): 
      _wto_components(boost::make_shared<wto_component_list_t>()), 
      _dfn_table(boost::make_shared<dfn_table_t>()), 
      _num(0), _stack(boost::make_shared<stack_t>()), 
      _nesting_table(boost::make_shared<nesting_table_t>()) {
      crab::ScopedCrabStats __st__("Fixpo.WTO");

      this->visit(cfg, cfg.entry(), this->_wto_components);
      this->_dfn_table.reset();
      this->_stack.reset();
      this->build_nesting();
    }

    // deep copy
    wto(const wto_t &other):
      _wto_components(boost::make_shared<wto_component_list_t>(*other._wto_components)),
      _dfn_table(other._dfn_table ?
		 boost::make_shared<dfn_table_t>(*other._dfn_table):
		 nullptr),
      _num(other._num),
      _stack(other._stack ?
	     boost::make_shared<stack_t>(*other._stack):
	     nullptr),
      _nesting_table(boost::make_shared<nesting_table_t>(*other._nesting_table)) { }

    wto(const wto_t &&other):
      _wto_components(boost::move(other._wto_components)),
      _dfn_table(boost::move(other._dfn_table)),
      _num(other._num),
      _stack(boost::move(other._stack)),
      _nesting_table(boost::move(other._nesting_table)) { }      
      
    wto_t& operator=(const wto_t &other) {
      if (this != &other) {
	this->_wto_components = other._wto_components;
	this->_dfn_table = other._dfn_table;
	this->_num = other._num;
	this->_stack = other._stack;
	this->_nesting_table = other._nesting_table;
      }
      return *this;
    }

    iterator begin() {
      return boost::make_indirect_iterator(_wto_components->begin());
    }
    
    iterator end() {
      return boost::make_indirect_iterator(_wto_components->end());      
    }

    const_iterator begin() const {
      return boost::make_indirect_iterator(_wto_components->begin());
    }
    
    const_iterator end() const {
      return boost::make_indirect_iterator(_wto_components->end());      
    }


    wto_nesting_t nesting(NodeName n) {
      typename nesting_table_t::iterator it = this->_nesting_table->find(n);
      if (it == this->_nesting_table->end()) {
        CRAB_ERROR("WTO nesting: node ", n," not found");
      } else {
        return it->second;
      }
    }

    void accept(wto_component_visitor<NodeName, CFG> *v) {
      for (iterator it = this->begin(); it != this->end(); ++it) {
        it->accept(v);
      }
    }
    
    void write(crab::crab_os& o) const {
      for (const_iterator it = this->begin(); it != this->end(); ) {
        const wto_component_t& c = *it;
        o << c;
        ++it;
        if (it != this->end()) {
          o << " ";
        }
      }      
    }
    
    friend crab::crab_os& operator<<(crab::crab_os &o, const wto_t &wto) {
      wto.write(o);
      return o;
    }    
    
  }; // class wto

} // namespace ikos

