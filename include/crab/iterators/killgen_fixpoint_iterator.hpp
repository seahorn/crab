#ifndef KILLGEN_FIXPOINT_ITERATOR_HPP
#define KILLGEN_FIXPOINT_ITERATOR_HPP

/**
  * Specialized fixpoint iterators and domains for kill-gen problems.
  */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg_bgl.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>

#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/patricia_trees.hpp>

#include <boost/optional.hpp>

namespace crab {

  namespace domains {

    // A wrapper for discrete_domain (i.e,. set of Element)
    template<class Element>
    class flat_killgen_domain: public ikos::writeable {
         
     private:

      typedef flat_killgen_domain<Element> flat_killgen_domain_t;
      typedef discrete_domain<Element> discrete_domain_t;
      
     public:

      typedef typename discrete_domain_t::iterator iterator;
      typedef Element element_t;

     private:

      discrete_domain_t _inv;
      
     public:
      
      flat_killgen_domain(discrete_domain_t inv)
          : ikos::writeable(), _inv(inv){ } 
      
      static flat_killgen_domain_t top() {
        return flat_killgen_domain(discrete_domain_t::top());
      }
      
      static flat_killgen_domain_t bottom() {
        return flat_killgen_domain(discrete_domain_t::bottom());
      }
      
      flat_killgen_domain()
          : ikos::writeable(), _inv(discrete_domain_t::bottom()){ }
      
      flat_killgen_domain(Element e)
          : ikos::writeable(), _inv(e) { }
      
      flat_killgen_domain(const flat_killgen_domain_t &o)
          : ikos::writeable(), _inv(o._inv) { } 
          
      flat_killgen_domain(flat_killgen_domain_t &&o)
          : ikos::writeable(), _inv(std::move(o._inv)) { } 
      
      flat_killgen_domain_t& operator=(const flat_killgen_domain_t &other) {
        if (this != &other) 
          _inv = other._inv;
        return *this;
      }
      
      flat_killgen_domain_t& operator=(flat_killgen_domain_t &&other) {
        _inv = std::move(other._inv);
        return *this;
      }
      
      iterator begin() { return _inv.begin(); }
      
      iterator end() { return _inv.end(); }
      
      unsigned size() { return _inv.size(); }
      
      bool is_bottom() { return _inv.is_bottom(); }
         
      bool is_top() { return _inv.is_top(); }

      bool operator==(flat_killgen_domain_t other) {
	return *this <= other && other <= *this;
      }
      
      bool operator<=(flat_killgen_domain_t other) {
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
         
      void operator-=(flat_killgen_domain_t other) {
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
      
      void operator+=(flat_killgen_domain_t other) {
        if (is_top () || other.is_bottom ()) {
          return;
        } else if (other.is_top ()) {
          _inv = discrete_domain_t::top();
        } else {
          _inv = (_inv | other._inv);
        }
      }
      
      flat_killgen_domain_t operator|(flat_killgen_domain_t other) {
        return (_inv | other._inv);
      }
      
      flat_killgen_domain_t operator&(flat_killgen_domain_t other) {
        return (_inv & other._inv);
      }
         
      void write(crab_os& o) { _inv.write(o); }
      
    };

    // To represent sets of pairs (Key,Value). 
    // Bottom means empty set rather than failure.
    template <typename Key, typename Value>
    class separate_killgen_domain: public ikos::writeable {
      
    private:
      typedef ikos::patricia_tree<Key,Value> patricia_tree_t;
      typedef typename patricia_tree_t::unary_op_t unary_op_t;
      typedef typename patricia_tree_t::binary_op_t binary_op_t;
      typedef typename patricia_tree_t::partial_order_t partial_order_t;

    public:
      typedef separate_killgen_domain<Key,Value > separate_killgen_domain_t;
      typedef typename patricia_tree_t::iterator iterator;
      typedef Key key_type;
      typedef Value value_type;
      
    private:
      bool _is_top;
      patricia_tree_t _tree;
    
    public: 
      class bottom_found { };
      
      class join_op: public binary_op_t {
	boost::optional< Value > apply(Value x, Value y) {
	  Value z = x.operator|(y);
	  if (z.is_top()) {
	    return boost::optional<Value>();
	  } else {
	    return boost::optional<Value>(z);
	  }
	}
	bool default_is_absorbing() { return false; }
      }; // class join_op

      class meet_op: public binary_op_t {
	boost::optional< Value > apply(Value x, Value y) {
	  Value z = x.operator&(y);
	  if (z.is_bottom()) {
	    throw bottom_found();
	  } else {
	    return boost::optional<Value>(z);
	  }
	};
	bool default_is_absorbing() { return true; }
      }; // class meet_op
    
      class domain_po: public partial_order_t {
	bool leq(Value x, Value y) { return x.operator<=(y); }
	bool default_is_top() { return false; }
      }; // class domain_po
    
   public:
      
      static separate_killgen_domain_t top() {
	return separate_killgen_domain_t (true);
      }
      
      static separate_killgen_domain_t bottom() {
	return separate_killgen_domain_t (false);
      }
    
   private:
      
      static patricia_tree_t apply_operation(binary_op_t& o, 
					     patricia_tree_t t1, 
					     patricia_tree_t t2) {
	t1.merge_with(t2, o);
	return t1;
      }
    
      separate_killgen_domain(patricia_tree_t t)
	: _is_top(false), _tree(t) { }
      
      separate_killgen_domain(bool b)
	: _is_top(b) { }
    
    public:
      
      separate_killgen_domain()
	: _is_top(false), _tree (patricia_tree_t()) { }

      separate_killgen_domain(const separate_killgen_domain_t& o)
	: _is_top(o._is_top), _tree(o._tree) { }
    
      separate_killgen_domain_t& operator=(separate_killgen_domain_t o) {
	this->_is_top = o._is_top;
	this->_tree = o._tree;
	return *this;
      }

      iterator begin() const {
	if (this->is_top()) {
	  CRAB_ERROR("Separate killgen domain: trying to invoke iterator on top");
	} else {
	  return this->_tree.begin();
	}
      }
    
      iterator end() const {
	if (this->is_top()) {
	  CRAB_ERROR("Separate killgen domain: trying to invoke iterator on top");
	} else {
	  return this->_tree.end();
	}
      }

      bool is_top() const {
	return _is_top;
      }
      
      bool is_bottom() const {
	return (!is_top () && _tree.empty ());
      }
    
    
      bool operator<=(separate_killgen_domain_t o) {
	domain_po po; 
	return (o.is_top() || (!is_top() && (_tree.leq (o._tree, po))));
      }
    
      separate_killgen_domain_t operator|(separate_killgen_domain_t o) {
	if (is_top() || o.is_top ()) {
	  return separate_killgen_domain_t::top();
	} else {
	  join_op op;
	  return separate_killgen_domain_t(apply_operation(op, _tree, o._tree));
	}
      }
    
      separate_killgen_domain_t operator&(separate_killgen_domain_t o) {
	if (is_top ()) {
	  return o;
	} else if (o.is_top()) {
	  return *this;
	} else {
	  try {
	    meet_op op;
	    return separate_killgen_domain_t(apply_operation(op, _tree, o._tree));
	  }
	  catch (bottom_found& exc) {
	    return separate_killgen_domain_t::bottom ();
	  }
	}
      }

      void set(Key k, Value v) {
	if (!is_top ()) {
	  // if (v.is_bottom()) {
	  //   this->_tree.remove(k);
	  // } else {
	  //   this->_tree.insert(k, v);
	  // }
	  this->_tree.insert(k, v);	  
	}
      }
    
      separate_killgen_domain_t& operator-=(Key k) {
	if (!is_top ()) {
	  _tree.remove(k);
	}
	return *this;
      }
    
      Value operator[](Key k) {
	if (is_top ())
	  return Value::top ();
	else {
	  boost::optional< Value > v = _tree.lookup(k);
	  if (v) {
	    return *v;
	  } else {
	    return Value::bottom();
	  }
	}
      }
    
      void write(crab::crab_os& o) {
	if (this->is_top()) {
	  o << "{...}";
	} if (_tree.empty ()) {
	  o << "_|_";
	}
	else {
	  o << "{";
	  for (typename patricia_tree_t::iterator it = this->_tree.begin(); 
	       it != this->_tree.end(); ) {
	    Key k = it->first;
	    k.write(o);
	    o << " -> ";
          Value v = it->second;
          v.write(o);
          ++it;
          if (it != this->_tree.end()) {
            o << "; ";
	  }
	  }
	  o << "}";
	}
      }
    }; // class separate_killgen_domain
    
  } // end namespace domains

  
  namespace iterators {
    
    // API for a kill-gen analysis operations
    template<class CFG, class Dom>
    class killgen_operations_api {

     public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;    
      typedef Dom killgen_domain_t;

     protected:

      CFG _cfg;

     public:

      killgen_operations_api (CFG cfg): _cfg(cfg) { }

      virtual ~killgen_operations_api() { }

      // whether forward or backward analysis
      virtual bool is_forward () = 0;

      // initial state
      virtual Dom entry() = 0;
 
      // (optional) initialization for the fixpoint
      virtual void init_fixpoint () = 0;

      // confluence operator
      virtual Dom merge(Dom, Dom) = 0;

      // analyze a basic block
      virtual Dom analyze (basic_block_label_t, Dom) = 0;

      // analysis name
      virtual std::string name () = 0;
    };

    // A simple fixpoint for a killgen analysis
    template<class CFG, class KgAnalysisOps>
    class killgen_fixpoint_iterator {

     public:
      
      typedef typename CFG::basic_block_label_t basic_block_label_t;
      typedef typename KgAnalysisOps::killgen_domain_t killgen_domain_t;
      typedef boost::unordered_map<basic_block_label_t,killgen_domain_t> inv_map_t;
      typedef typename inv_map_t::iterator iterator;
      typedef typename inv_map_t::const_iterator const_iterator;
      
     protected:

      CFG _cfg;
      inv_map_t _in_map;
      inv_map_t _out_map;

     private:

      KgAnalysisOps _analysis;

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

      killgen_fixpoint_iterator (CFG cfg)
	: _cfg (cfg), _analysis (_cfg) { }

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
