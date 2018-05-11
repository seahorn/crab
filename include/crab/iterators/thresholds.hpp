#pragma once

#include <crab/common/bignums.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/common/debug.hpp>
#include <crab/cfg/cfg.hpp>
#include <crab/iterators/wto.hpp>

#include <boost/unordered_map.hpp>
#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <climits>

namespace crab {

   namespace iterators {

     /** 
	 Class that represents a thresholds used by the widening operator
     **/
     template<typename Number>
     class thresholds {
       
      private:

       /// XXX: internal representation of a threshold
       typedef ikos::bound <Number> bound_t;
       
       std::vector<bound_t> m_thresholds;
       unsigned int m_size;

       template<class B1, class B2>
       static B2 convert_bounds_impl (B1 b1) {
	 B2 b2 (0); // some initial value it doesn't matter which one
	 ikos::bounds_impl::convert_bounds (b1, b2);
	 return b2;
       }
       
      public:
       
       thresholds (int size = UINT_MAX) : m_size (size) { 
         m_thresholds.push_back (bound_t::minus_infinity ());
         m_thresholds.push_back (0);
	 // useful thresholds for wrapped domains
	 #if 0
	 m_thresholds.push_back (bound_t("-2147483648"));	 
	 m_thresholds.push_back (bound_t("-32768"));
	 m_thresholds.push_back (bound_t("-128"));
	 m_thresholds.push_back (bound_t("127"));
	 m_thresholds.push_back (bound_t("255"));
	 m_thresholds.push_back (bound_t("32767"));
	 m_thresholds.push_back (bound_t("65535"));	 
	 m_thresholds.push_back (bound_t("2147483647"));
	 #endif 
         m_thresholds.push_back (bound_t::plus_infinity ());	 
       }

       unsigned size (void) const { return m_thresholds.size(); }
       
       template<typename N>
       void add (ikos::bound<N> v1) { 
         if (m_thresholds.size () < m_size) {
       	   bound_t v = convert_bounds_impl<ikos::bound<N>, bound_t> (v1);
           if (std::find
	       (m_thresholds.begin (), m_thresholds.end (), v) == m_thresholds.end ()) {
             auto ub = std::upper_bound (m_thresholds.begin (), m_thresholds.end (), v);
             m_thresholds.insert (ub, v);
           }
         }
       }
              
       template<typename N>
       ikos::bound<N> get_next (ikos::bound<N> v1) const {
	 if (v1.is_plus_infinity()) return v1;
	 bound_t v = convert_bounds_impl<ikos::bound<N>, bound_t> (v1);
	 bound_t t = m_thresholds [m_thresholds.size () - 1];	 
         auto ub = std::upper_bound (m_thresholds.begin (), m_thresholds.end (), v);         
         if (ub != m_thresholds.end ())
	   t = *ub;
	 return convert_bounds_impl<bound_t, ikos::bound<N> > (t);
       }

       template<typename N>
       ikos::bound<N> get_prev (ikos::bound<N> v1) const {
	 if (v1.is_minus_infinity()) return v1;
	 bound_t v = convert_bounds_impl<ikos::bound<N>, bound_t> (v1);	 
         auto lb = std::lower_bound (m_thresholds.begin (), m_thresholds.end (), v);         
         if (lb != m_thresholds.end ()) {
           --lb;
           if (lb != m_thresholds.end ()) {
	     return convert_bounds_impl<bound_t, ikos::bound<N> > (*lb);	 	     
	   }
         }
	 return convert_bounds_impl<bound_t, ikos::bound<N> > (m_thresholds [0]);
       }

       void write (crab_os &o) const { 
         o << "{";
         for (typename std::vector<bound_t>::const_iterator it = m_thresholds.begin (), 
                  et= m_thresholds.end (); it != et; ) {
           bound_t b (*it);
           b.write (o);
           ++it;
           if (it != m_thresholds.end ())
             o << ",";
         }
         o << "}";
       }
       
     };

     template<typename Number>
     crab_os& operator<< (crab_os& o, const thresholds<Number>& t) {
       t.write (o);
       return o;
     }

     // explicit instantiation
     typedef thresholds<ikos::q_number> thresholds_t;     
     
     /**
	Collect thresholds per wto cycle  (i.e. loop)
     **/
     template< typename NodeName, typename CFG>
     class wto_thresholds: public ikos::wto_component_visitor< NodeName, CFG > {
       
     public:
       
       typedef ikos::wto_vertex< NodeName, CFG > wto_vertex_t;
       typedef ikos::wto_cycle< NodeName, CFG > wto_cycle_t;
       typedef crab::iterators::thresholds_t thresholds_t;      
       typedef boost::unordered_map<NodeName, thresholds_t> thresholds_map_t;
       
     private:
       
       // the cfg
       CFG m_cfg;
       // maximum number of thresholds
       size_t m_max_size;
       // keep a set of thresholds per wto head
       thresholds_map_t m_head_to_thresholds;
       // the top of the stack is the current wto head
       std::vector<NodeName> m_stack;

       typedef typename CFG::basic_block_t basic_block_t;
       typedef typename CFG::number_t number_t;
       typedef typename CFG::varname_t varname_t;       
       typedef ikos::bound<number_t> bound_t;

       typedef crab::cfg::assume_stmt<number_t,varname_t> assume_t;
       typedef crab::cfg::select_stmt<number_t,varname_t> select_t;
       typedef crab::cfg::assignment<number_t,varname_t>  assign_t;
       
       void get_thresholds(const basic_block_t &bb, thresholds_t &thresholds) {
	 for (auto const&i : boost::make_iterator_range (bb.begin (), bb.end ())) {
	   bound_t t = bound_t::plus_infinity ();
	   if (i.is_assume ()) {
	     auto a = static_cast<const assume_t*> (&i);
	     t = bound_t(-(a->constraint ().expression ().constant ()));
	   }
	   else if (i.is_select ()) {
	     auto s = static_cast<const select_t*> (&i);
	     t = bound_t(-(s->cond ().expression ().constant ()));
	   }
	   else if (i.is_assign ()) {
	     auto a = static_cast<const assign_t*> (&i);
	     if (a->rhs ().is_constant ())
	       t = bound_t(-(a->rhs().constant ()));
	   }	    
	   
	   if (t != bound_t::plus_infinity ()) {
	     // XXX: for code pattern like this "while(x<t) {x+=k;}"
	     // note that the condition (x<t) is translated to
	     // (x<=t-1) so an useful threshold would be t+1+k.
	     // Since we don't keep track of how x is incremented or
	     // decremented we choose arbitrarily k=1.
	     thresholds.add(t+2);
	   }
	 }
       }
       
     public:
       
       wto_thresholds(CFG cfg, size_t max_size)
	 : m_cfg(cfg), m_max_size (max_size) { }
       
       void visit(wto_vertex_t& vertex) {
	 if (m_stack.empty()) return;
	 
	 NodeName head = m_stack.back();
	 auto it = m_head_to_thresholds.find(head);
	 if (it != m_head_to_thresholds.end()) {
	   thresholds_t & thresholds = it->second;
	   typename CFG::basic_block_t & bb = m_cfg.get_node(vertex.node());
	   get_thresholds(bb, thresholds);
	 } else {
	   CRAB_ERROR("No head found while gathering thresholds");
	 }
       }
       
       void visit(wto_cycle_t& cycle) {
	 thresholds_t thresholds(m_max_size);
	 typename CFG::basic_block_t & bb = m_cfg.get_node(cycle.head());
	 get_thresholds(bb, thresholds);

	 #if 1
	 // XXX: if we want to consider constants from loop
	 // initializations
	 for (auto pre : boost::make_iterator_range(bb.prev_blocks())) {
	   if (pre != cycle.head()) {
	     typename CFG::basic_block_t & pred_bb = m_cfg.get_node(pre);
	     get_thresholds(pred_bb, thresholds);
	   }
	 }
	 #endif 
	 
	 m_head_to_thresholds.insert(std::make_pair(cycle.head(), thresholds));
	 m_stack.push_back(cycle.head());
	 for (typename wto_cycle_t::iterator it = cycle.begin(); it != cycle.end(); ++it) {
	   it->accept(this);
	 }
	 m_stack.pop_back();
       }
       
       const thresholds_map_t& get_thresholds_map() const {
	 return m_head_to_thresholds;
       }
       
       void write (crab::crab_os& o) const {
	 for (auto &kv: m_head_to_thresholds) {
	   o << crab::cfg_impl::get_label_str(kv.first) << "=" << kv.second << "\n";
	 }
       }
       
     }; // class wto_thresholds
     
     template<typename NodeName, typename CFG>
     crab::crab_os& operator<<(crab::crab_os& o, const wto_thresholds<NodeName, CFG>& t) {
       t.write (o);
       return o;
     }
     
   } // end namespace iterators
} // end namespace crab
