#ifndef THRESHOLDS_HPP
#define THRESHOLDS_HPP

#include <crab/domains/intervals.hpp>
#include <algorithm>
#include <climits>

using namespace std;

namespace crab {

   namespace iterators {
 
     template<typename Number>
     class thresholds {
       
      public:
       
       typedef bound <Number> bound_t;
       
      private:
       
       vector<bound_t> m_thresholds;
       int m_size;

      public:
       
       thresholds (int size = INT_MAX) : m_size (size) { 
         m_thresholds.push_back (bound_t::minus_infinity ());
         m_thresholds.push_back (0);
         m_thresholds.push_back (bound_t::plus_infinity ());
       }

       void add (bound_t v) { 
         if (m_thresholds.size () < m_size) {
           if (std::find (m_thresholds.begin (), m_thresholds.end (), v) == m_thresholds.end ()) {
             auto ub = std::upper_bound (m_thresholds.begin (), m_thresholds.end (), v);
             m_thresholds.insert (ub, v);
           }
         }
       }

       bound_t get_next (bound_t v) const { 
         auto ub = std::upper_bound (m_thresholds.begin (), m_thresholds.end (), v);         
         if (ub != m_thresholds.end ())
           return *ub;
         
         return m_thresholds [m_thresholds.size () - 1];
       }
       
       bound_t get_prev (bound_t v) const { 
         auto lb = std::lower_bound (m_thresholds.begin (), m_thresholds.end (), v);         
         if (lb != m_thresholds.end ()) {
           --lb;
           if (lb != m_thresholds.end ())
             return *lb;
         }
         return m_thresholds [0];
       }
       
       void write (crab_os &o) const { 
         o << "{";
         for (typename vector<bound_t>::const_iterator it = m_thresholds.begin (), 
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
   
   }
}


#endif 
