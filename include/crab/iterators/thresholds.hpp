#ifndef THRESHOLDS_HPP
#define THRESHOLDS_HPP

#include <crab/common/bignums.hpp>
#include <crab/domains/intervals.hpp>
#include <algorithm>
#include <climits>

namespace crab {

   namespace iterators {

     
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
         m_thresholds.push_back (bound_t::plus_infinity ());
       }

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
   
   }
}


#endif 
