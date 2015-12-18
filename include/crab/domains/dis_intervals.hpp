#ifndef DISJUNCTIVE_INTERVALS_HPP
#define DISJUNCTIVE_INTERVALS_HPP

/* 
  DisIntervals: a domain of disjunction of intervals described in the
  paper "Clousot: Static Contract Checking with Abstract
  Interpretation" by Fahndrich and Logozzo. It is based on the
  implementation available in CodeContracts.
*/

/// Uncomment for enabling debug information
// #include <crab/common/dbg.hpp>

#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>
#include <crab/domains/division_operators_api.hpp>

using namespace boost;
using namespace ikos;
using namespace std;

namespace crab {
   namespace domains {

   template< typename Number >
   class dis_interval: public writeable {
    
    public:

     typedef bound< Number > bound_t;
     typedef interval< Number > interval_t;
     typedef dis_interval <Number> dis_interval_t;

    private:

     typedef enum { BOT, FINITE, TOP } state_t;
     typedef vector<interval_t> list_intervals_t;
    
    public:

     typedef typename list_intervals_t::iterator iterator;
     typedef typename list_intervals_t::const_iterator const_iterator;

    private:

     state_t _state;
     list_intervals_t _list;
     
    public:
     
     static dis_interval_t top () { 
       return dis_interval_t (TOP);
     } 
     
     static dis_interval_t bottom() {
       return dis_interval_t (BOT);
     }
     
    private:

     // FIXME: this is assuming Number=z_number
     static bool are_consecutive (const interval_t& i1, const interval_t& i2) {
       assert (i1.lb().is_finite () && i1.ub().is_finite () && 
               i2.lb().is_finite () && i2.ub().is_finite ());

       return ( (i1.lb () <= i2.lb () && i1.ub () <= i2.ub ()) &&
                (i1.ub () + Number(1) == i2.lb ())) ||
              ( (i2.lb () <= i1.lb () && i2.ub () <= i1.ub ()) &&
                (i2.ub () + Number(1) == i1.lb ()));
     }
     
     static bool overlap (const interval_t& i1, const interval_t& i2) {
       assert (i1.lb().is_finite () && i1.ub().is_finite () && 
               i2.lb().is_finite () && i2.ub().is_finite ());

       return (i2.lb () <= i1.ub () && i1.lb () <= i2.ub ());
     }

     struct IsOnTheLeft{
       bool operator() (const interval_t& i1, const interval_t& i2) const {
         assert (i1.lb().is_finite () && i1.ub().is_finite () && 
                 i2.lb().is_finite () && i2.ub().is_finite ());

         return ((i1.ub () <= i2.lb ()) && (i1.ub () != i2.lb ()));
       }
     };
     
     static bool check_well_formed (const dis_interval_t &x) {

       if (x.is_top () || x.is_bottom ()) {
         return true;
       }

       if (!x.is_finite ())  {
         CRAB_ERROR ("sanity check -- state should be finite\n",
                     x, " not well formed");
         return false;
       }
       
       if (x._list.empty ()) {
         CRAB_ERROR ("sanity check -- list cannot be empty\n",
                     x, " not well formed");
         return false;
       }

       if (x._list.size () == 1) {
         if (x._list[0].is_top () || x._list [0].is_bottom ()) {
           CRAB_ERROR ("sanity check -- cannot be top or bottom\n",
                       x, " not well formed");
           return false;
         } else {
           return true;
         }
       }

       // -- check the list of intervals is strictly sorted. This also
       // -- checks for duplicates but it won't complain with two
       // -- intervals like [0,2], [3,4]
       IsOnTheLeft comp;
       list_intervals_t tmp (x._list);
       unsigned prev = 0;
       for (unsigned cur=1; cur < tmp.size (); cur++,prev++) {
         if (comp (tmp [prev], tmp [cur]))
           continue;
         
         CRAB_ERROR ("sanity check -- list is not strictly sorted: ",
                     tmp [prev], " not leq ", tmp[cur], "\n",
                     x, " not well formed");
         return false;
       }

       return true;
     }
     
     list_intervals_t normalize (list_intervals_t l, bool& is_bottom) {
       
       if (l.size () <= 1) {
         CRAB_DEBUG ("-- Normalize: singleton ", l[0]);
         is_bottom = false;
         return l;
       }
       
       IsOnTheLeft comp;
       std::sort (l.begin (), l.end(), comp);
       
       list_intervals_t res;
       interval_t prev = interval_t::top (); 
       unsigned int bottoms = 0;

       for (unsigned int i=0; i < l.size (); ++i){
         interval_t intv = l[i];
         
         if (prev == intv) { 
           CRAB_DEBUG ("-- Normalize: duplicate");
           continue;
         }
         
         if (intv.is_bottom ()) {
           CRAB_DEBUG ("-- Normalize: bottom interval");
           bottoms++;
           continue;
         }

         if (intv.is_top ()) {
           CRAB_DEBUG ("-- Normalize: top interval");
           is_bottom = false;
           return list_intervals_t ();
         }
         
         if (!prev.is_top ()) { 
           
           bool refined = true;
           while (refined && res.size () > 0) {
             prev = res [res.size () -1];
             if (overlap (prev, intv) || are_consecutive (prev,intv)) {
               CRAB_DEBUG ("-- Normalize: overlapping or consecutive intervals", 
                           prev,  " and ", intv);               
               res.pop_back ();
               intv = prev | intv;
             }
             else if (intv <= prev) {
               CRAB_DEBUG ("-- Normalize: skipping subsumed interval",
                           prev,  " and ", intv);               
               goto next_iter;
             }
             else {
               refined = false;
             }
           }
           
         }
         
         prev = intv;
         CRAB_DEBUG ("-- Normalize: adding ",  intv);
         if (intv.is_top ()) {
           CRAB_DEBUG ("-- Normalize: top interval");
           is_bottom = false;
           return list_intervals_t ();
         }

         res.push_back (intv);
         next_iter: ;;
         
       }
       
       is_bottom = (bottoms == l.size ());
       CRAB_DEBUG ("-- Normalize: number of bottoms = ", bottoms);

       //cout << "-- Normalize result="; for (auto i: res) { cout << i << "|";} cout << "\n";
       return res;
     }
     
     
     dis_interval(state_t state): _state (state) { }

     dis_interval(list_intervals_t l, bool Normalize = true): 
         _state (FINITE), _list (l) { 

       // TODO: make it a template parameter
       const std::size_t max_num_disjunctions = 50;

       if (Normalize) {
         bool is_bottom = false;
         list_intervals_t res = normalize (_list, is_bottom);
         if (is_bottom) {
           _state = BOT;
           _list.clear ();
         } else if (res.empty ()) {
           _state = TOP;
           _list.clear ();
         }
         else {
           std::swap (_list, res);
         }
       }


       if (is_finite () && _list.size () >= max_num_disjunctions) {
         // TODO: rather than merging all intervals do a more graceful
         // degradation e.g., start by merging the nearest intervals.
         CRAB_WARN (" reached maximum allowed number of disjunctions. ",
                    "Merging all intervals ... ");
         interval_t alljoined = approx (_list);
         _list.clear ();
         _list.push_back (alljoined);
       }

       assert (check_well_formed (*this));
     }
          
     // pre: x is normalized
     interval_t approx (list_intervals_t x) const {
       if (x.empty ()) 
         CRAB_ERROR ("list should not be empty");
       
       if (x.size () <= 1) 
         return x[0];
       else 
         return (x [0] | x [x.size () - 1]);
     }
     
    public:

     dis_interval (): writeable (), _state (TOP) { }

     dis_interval(interval_t i): 
         writeable (), _state (FINITE) { 
       if (i.is_top ()) 
         _state = TOP;
       else if (i.is_bottom ())
         _state = BOT;
       else 
         _list.push_back (i);
     }

     dis_interval(const dis_interval_t& i): 
         writeable(), _state (i._state), _list (i._list) { }

     dis_interval(dis_interval_t&& i): 
         writeable(), 
         _state (std::move (i._state)),
         _list (std::move (i._list)) { }
     
     dis_interval_t& operator=(const dis_interval_t& i){
       if (this != &i){
         _state = i._state;
         _list = i._list;
       }
       return *this;
     }

     dis_interval_t& operator=(dis_interval_t&& i){
       if (this != &i){
         _state = std::move (i._state);
         _list = std::move (i._list);
       }
       return *this;
     }
     
     bool is_bottom() const { 
       return (_state == BOT); 
     }
     
     bool is_top() const {
       return (_state == TOP);
     }
     
     bool is_finite () const {
       return (_state == FINITE);
     }

     // for the interval solver
     dis_interval_t lower_half_line () const {
       if (is_bottom()) {
         return bottom();
       } else {
         auto f = [](interval_t x){ return x.lower_half_line ();};
         return apply_unary_op (*this, f);
       }
     }

     // for the interval solver
     dis_interval_t upper_half_line () const {
       if (is_bottom()) {
         return bottom();
       } else {
         auto f = [](interval_t x){ return x.upper_half_line ();};
         return apply_unary_op (*this, f);
       }
     }

     // for the interval solver
     boost::optional< Number > singleton() const {
       interval_t i = approx ();
       return i.singleton ();
     }

     // iterate over the disjunctive intervals
     iterator begin () { return _list.begin (); }

     iterator end () { return _list.end (); }

     const_iterator begin () const { return _list.begin (); }

     const_iterator end () const { return _list.end (); }
       
     interval_t approx () const {
       if (is_bottom ()) return interval_t::bottom ();
       if (is_top ()) return interval_t::top ();
       return approx (_list);
     }

     // required for separate_domains
     // pre: *this and o are normalized
     bool operator==(dis_interval_t o) const {
#if 0
       // -- semantic check
       return (*this <= o && o <= *this);
#else
       // -- syntactic check
       if (is_bottom () && o.is_bottom ()) 
         return true;
       else if (is_top () && o.is_top ())
         return true;
       else if (is_bottom () || o.is_bottom ())
         return false;
       else if (is_top () || o.is_top ())
         return false;
       else {
         if (_list.size () != o._list.size ())
           return false;
         else {
           unsigned int j=0;
           for (unsigned int i=0; i < _list.size (); ++i,++j){
             if (! (_list [i] == o._list [j]))
               return false;
           }
           return true;
         }
       }
#endif 
     }

     bool operator<=(dis_interval_t o) const {
       if (this->is_bottom()) {
         return true;
       } else if (o.is_bottom()) {
         return false;
       } else {

         unsigned j=0;
         for (unsigned int i=0; i < _list.size (); i++) {
           for (; j < o._list.size (); j++){
             if (_list[i] <= o._list[j])
               goto next_iter;
           }
           return false; // _list[i] is not included in any interval of o._list
           next_iter: ;;
         }
         return true;
       }
     }
     
     // pre: *this and o are normalized
     dis_interval_t operator|(dis_interval_t o) const {

       CRAB_DEBUG ("Join of ", *this, " and ", o);

       if (this->is_bottom()) {
         return o;
       } else if (o.is_bottom()) {
         return *this;
       } else if (this->is_top ()) {
         return *this;
       } else if (o.is_top ()) {
         return o;
       } else {
         
         unsigned int i=0;
         unsigned int j=0;
         list_intervals_t res;
         res.reserve (_list.size () + o._list.size ());
         
         while (i < _list.size () && j < o._list.size ()) {
           CRAB_DEBUG ("Join -- left operand =", _list[i], " right operand=", o._list[j]);

           if (_list[i].is_top () || o._list[j].is_top ()) {
             CRAB_DEBUG ("Join -- One of the operands is top");
             return dis_interval_t ();
           }
           else if (_list[i].is_bottom ()) {
             CRAB_DEBUG ("Join -- Left operand is bottom");
             i++;
           }
           else if (o._list[j].is_bottom ()) {
             CRAB_DEBUG ("Join -- Right operand is bottom");
             j++;
           }
           else if (_list[i] == o._list[j]) {
             CRAB_DEBUG ("Join -- Left operand is equal to right");
             res.push_back(_list[i]);
             i++; j++;
           }
           else if (_list[i] <= o._list[j]) {
             CRAB_DEBUG ("Join -- Left operand is included in the right");
             res.push_back(o._list[j]);
             i++;j++; 
           }
           else if (o._list[j] <= _list[i]) {
             CRAB_DEBUG ("Join -- Right operand is included in the left");
             res.push_back (_list[i]);
             i++;j++;
           }
           else if (overlap (_list[i], o._list[j]) || 
                    are_consecutive (_list[i], o._list[j])) {
             CRAB_DEBUG ("Join -- Left ",_list[i], " and right ", o._list[j],
                         " operands overlap or are consecutive");
             res.push_back (_list [i] | o._list[j]);
             i++; j++;
           }
           else if (IsOnTheLeft () (_list[i], o._list [j])) {
             CRAB_DEBUG ("Join -- Left operand ", _list[i] , 
                         " is on the left of the right operand ", o._list[j]);
             res.push_back (_list [i]);
             i++;
           }
           else {
             assert (IsOnTheLeft () (o._list [j], _list [i]));
                     
             CRAB_DEBUG ("Join -- Right operand ", o._list[j] ,
                         " is on the left of the right operand ", _list[i]);
             res.push_back (o._list [j]);
             j++;
           }
         }
         
         // consume the rest of left operand
         while (i < _list.size ()){
           auto intv = _list [i];
           
           CRAB_DEBUG ("Join -- Adding the rest of left operand ", intv);
           
           bool refined = true;
           while (refined && res.size () > 0) {
             auto prev = res [res.size () -1];
             if (overlap (prev, intv) || are_consecutive (prev,intv)) {
               CRAB_DEBUG ("\t-- Join : overlapping or consecutive intervals", 
                           prev,  " and ", intv);               
               res.pop_back ();
               intv = prev | intv;
             }
             else if (intv <= prev) {
               CRAB_DEBUG ("\t-- Join: skipping subsumed interval",
                           prev,  " and ", intv);               
               goto next_iter_1;
             }
             else {
               refined = false;
             }
           }

           res.push_back (intv);
           next_iter_1:
           i++;
         }

         // consume the rest of right operand
         while (j < o._list.size ()){
           auto intv = o._list [j];
           CRAB_DEBUG ("Join -- Adding the rest of right operands ", intv);
           
           bool refined = true;
           while (refined && res.size () > 0) {
             auto prev = res [res.size () -1];
             if (overlap (prev, intv) || are_consecutive (prev,intv)) {
               CRAB_DEBUG ("\t-- Join : overlapping or consecutive intervals", 
                           prev,  " and ", intv);               
               res.pop_back ();
               intv = prev | intv;
             }
             else if (intv <= prev) {
               CRAB_DEBUG ("\t-- Join: skipping subsumed interval",
                           prev,  " and ", intv);               
               goto next_iter_2;
             }
             else {
               refined = false;
             }
           }

           res.push_back (intv);
           next_iter_2:
           j++;
         }

         if (res.empty ()) {
           CRAB_DEBUG ("Join result=_|_");
           return dis_interval_t (BOT);
         }
         else if (res.size () <= 1 && res[0].is_top ()) {
           CRAB_DEBUG ("Joing result=[-oo,+oo]");
           return dis_interval_t (TOP);
         }
         else {
           // It needs normalization. E.g., the join of {[0, 7] | [9, 11]}
           // and {[0, 6] | [8, 11]} returns {[0, 7] | [8, 11]} which
           // should be further simplified to [0,11].

           dis_interval_t join (res);
           CRAB_DEBUG ("Joing result=", join);

           return join;
         }
       }
     }
         
     // pre: *this and o are normalized     
     dis_interval_t operator&(dis_interval_t o) const {
       if (this->is_bottom() || o.is_bottom()) {
         return this->bottom();
       } else if (this->is_top ()) { 
         return o;
       } else if (o.is_top ()) {
         return *this;
       } else {

         list_intervals_t res;
         res.reserve (_list.size () + o._list.size ());

         bool is_bot = true;
         for (unsigned int i=0; i< _list.size (); ++i){
           for (unsigned int j=0; j< o._list.size (); ++j){
             auto meet = _list[i] & o._list[j];
             if (!meet.is_bottom ()) {
               res.push_back (meet);
               is_bot = false;
             }
           }
         }

         if (is_bot) {
           return (this->bottom ());
         } else if (res.empty ()) {
           return (this->top ());
         } else {
           return dis_interval_t (res); 
         }
       }
     }


    private:

     struct WidenOp {
       interval_t apply (interval_t before, interval_t after){ 
         return before || after;
       }
     };
     
     template<typename Thresholds>
     struct WidenWithThresholdsOp {
       const Thresholds & m_ts;
       
       WidenWithThresholdsOp (const Thresholds &ts): m_ts (ts) { }
       
       interval_t apply (interval_t before, interval_t after) { 
         return before.widening_thresholds (after, m_ts);
       }
     };
     
     template <typename WidenOp>
     dis_interval_t widening (dis_interval_t o, WidenOp widen_op) const {
       if (this->is_bottom()) {
         return o;
       } else if (o.is_bottom()) {
         return *this;
       } else if (this->is_top ()) {
         return *this;
       } else if (o.is_top ()) {
         return o;
       } else {
         
         // --- trivial cases first 
         if (_list.size () == 1 && o._list.size () == 1) {  
           return dis_interval_t (widen_op.apply (_list[0], o._list[0]));
         }
         
         if (_list.size () == 1 && o._list.size () > 1) {
           interval_t x = approx (o._list);
           return dis_interval_t (widen_op.apply (_list [0], x));
         }
         
         if (_list.size () > 1 && o._list.size () == 1) {
           interval_t x = approx (_list);
           return dis_interval_t (widen_op.apply (x, o._list [0]));
         }
         
         assert (_list.size () >= 2 && o._list.size () >= 2);
         
         // The widening implemented in CodeContracts widens the
         // extremes and keep only stable intervals. For this query:
         // widening( [1,1] | [4, 4], [1,1] | [3,3] | [4, 4]) the
         // result would to be [1,1] | [4,4]. But this is not even an
         // upper bound of the right argument so it cannot be a
         // widening.

         // -- widen the extremes
         interval_t lb_widen = widen_op.apply (_list [0], o._list [0]);
         interval_t ub_widen = widen_op.apply (_list [_list.size () - 1],
                                               o._list [o._list.size () - 1]);
         
         list_intervals_t res;
         res.reserve (_list.size () + o._list.size ());
         
         res.push_back (lb_widen);

         // for (unsigned int i=1; i < _list.size () - 1; i++) {
         //   for (unsigned int j=1; j < o._list.size () - 1; j++) {
         //     if (o._list [j] <= _list [i]) { 
         //       res.push_back (_list [i]);
         //       break;
         //     }
         //   }
         // }

         // keep all the intervals, normalize will do the rest
         res.insert (res.end (), _list.begin () + 1, _list.end () - 1);
         res.insert (res.end (), o._list.begin () + 1, o._list.end () - 1);

         res.push_back (ub_widen);

         return dis_interval_t (res); 
       }
     }

    public:

     // pre: *this and o are normalized
     dis_interval_t operator|| (dis_interval_t o) const {
       WidenOp op;
       return widening (o, op);
     }
     
     // pre: *this and o are normalized
     template<typename Thresholds>
     dis_interval_t widening_thresholds (dis_interval_t o, const Thresholds &ts) {
       WidenWithThresholdsOp <Thresholds> op (ts);
       return widening (o, op);
     }
     
     // pre: *this and o are normalized
     dis_interval_t operator&&(dis_interval_t o) const {
       CRAB_WARN (" DisIntervals narrowing operator replaced with meet");

       return (*this & o);
     }

    private:

     // pre: x and y are normalized
     template<typename BinOp>
     dis_interval_t apply_bin_op (dis_interval_t x, dis_interval_t y, 
                                  BinOp op, 
                                  bool shortcut_top) const {
       
       // if shortcut_top is true then the result is top if one of the
       // operands is top

       if (x.is_bottom () || y.is_bottom ())
         return this->bottom ();

       if (x.is_top () && y.is_top ())
         return this->top ();

       if (shortcut_top && (x.is_top () || y.is_top ()))
         return this->top ();

       list_intervals_t res;
       res.reserve (x._list.size () + y._list.size ());
       bool is_bot = true;

       if (shortcut_top || (x.is_finite () && y.is_finite ())) {
         for (unsigned int i=0; i < x._list.size (); ++i){ 
           for (unsigned int j=0; j < y._list.size (); ++j){ 
             interval_t intv = op (x._list[i], y._list[j]);
             if (intv.is_bottom ())
               continue;

             if (intv.is_top ())
               return this->top ();
             
             is_bot = false;
             res.push_back (intv);
             
           }
         }
       }
       else {
         
         assert (x.is_top () || y.is_top ());

         dis_interval_t non_top_arg = (!x.is_top () ? x : y);
         bool is_non_top_left = (!x.is_top ());
         for (unsigned int i=0; i< non_top_arg._list.size (); ++i) {
           
           interval_t intv = (is_non_top_left ? 
                              op (non_top_arg._list[i], interval_t::top ()) :
                              op (interval_t::top () , non_top_arg._list[i]));

           if (intv.is_bottom ())
             continue;

           if (intv.is_top ())
             return this->top ();

           is_bot = false;
           res.push_back (intv);
             
         }
       }
       
       if (is_bot) 
         return this->bottom ();
       else
         return dis_interval_t (res);
     }
     

     // pre: x is normalized
     template<typename UnaryOp>
     dis_interval_t apply_unary_op (dis_interval_t x, UnaryOp op) const {
                                    
       if (x.is_bottom ())
         return this->bottom ();

       if (x.is_top ())
         return this->top ();

       assert (x.is_finite ());

       if (x._list.empty ())
         CRAB_ERROR ("list should not be empty");

       list_intervals_t res;
       res.reserve (x._list.size ());
       bool is_bot = true;
       
       for (unsigned int i=0; i < x._list.size (); ++i){ 
         interval_t intv = op (x._list[i]);
         if (intv.is_bottom ())
           continue;
         if (intv.is_top ())
           return this->top ();
         is_bot = false;
         res.push_back (intv);
       }
       
       if (is_bot) 
         return this->bottom ();
       else
         return dis_interval_t (res);
     }
     
    public:

     dis_interval_t operator+(dis_interval_t x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a + b;};
         return apply_bin_op (*this, x, f, true);
       }
     }
     
     dis_interval_t& operator+=(dis_interval_t x) {
       return this->operator=(this->operator+(x));
     }
     
     dis_interval_t operator-() const {
       if (this->is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a){ return -a;};
         return apply_unary_op (*this, f);
       }
     }
    
     dis_interval_t operator-(dis_interval_t x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a - b;};
         return apply_bin_op (*this, x, f, true); 
       }
     }

     dis_interval_t& operator-=(dis_interval_t x) {
       return this->operator=(this->operator-(x));
     }
    
     dis_interval_t operator*(dis_interval_t x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a * b;};
         return apply_bin_op (*this, x, f, true); 
       }
     }
    
     dis_interval_t& operator*=(dis_interval_t x) {
       return this->operator=(this->operator*(x));
     }
    
     dis_interval_t operator/(dis_interval_t x) {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a / b;};
         return apply_bin_op (*this, x, f, false); 
       }
     }
     
     dis_interval_t& operator/=(dis_interval_t x) {
       return this->operator=(this->operator/(x));
     }   

     void normalize () {
       if (is_bottom () || is_top ()) return;
       bool is_bot = false;
       list_intervals_t res = normalize (_list, is_bot);

       if (is_bot) {
         _state = BOT;
       } else if (res.empty ()) {
         _state = TOP;
       } else {
         _state = FINITE;
         std::swap (_list, res);
       }
     }

     
    void write(std::ostream& o) {
      if (is_bottom()) {
        o << "_|_";
      } else if (is_top ()) {
        o << "[-oo,+oo]";
      }
      else {
        assert (_state == FINITE);

        for (typename list_intervals_t::iterator it = _list.begin (), 
                 et= _list.end (); it!=et; ){
          o << *it;
          ++it;
          if (it != et)
            o << " | ";
        } 
      }
    }    
    
    // division and remainder operations

     dis_interval_t UDiv (dis_interval_t x ) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a / b;};
         return apply_bin_op (*this, x, f, false); 
       }
     }
     
     dis_interval_t SRem(dis_interval_t x)  const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.SRem(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }

     dis_interval_t URem(dis_interval_t x)  const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.URem(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }

     // bitwise operations

     dis_interval_t Trunc(unsigned /* width */) const {
       return *this;
     }
     
     dis_interval_t ZExt(unsigned /* width */) const {
       return *this;
     }
     
     dis_interval_t SExt(unsigned /* width */) const {
       return *this;
    }

     dis_interval_t And(dis_interval_t x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.And(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }
    
     dis_interval_t Or(dis_interval_t x)  const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.Or(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }
    
     dis_interval_t Xor(dis_interval_t x) const { 
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.Xor(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }
     
     dis_interval_t Shl(dis_interval_t x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.Shl(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }
     
     dis_interval_t LShr(dis_interval_t  x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.LShr(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }

     dis_interval_t AShr(dis_interval_t  x) const {
       if (this->is_bottom() || x.is_bottom()) {
         return this->bottom();
       } else {
         auto f = [](interval_t a, interval_t b){ return a.AShr(b);};
         return apply_bin_op (*this, x, f, false); 
       }
     }
     
   };//  class dis_interval

   template<typename N>
   std::ostream& operator<<(std::ostream& o, dis_interval<N> i) {
     i.write (o);
     return o;
   }

 } // end namespace domains
} // end namespace crab 
 
namespace ikos {
  namespace intervals_impl {
     /// --- for interval solver of disequalities

     typedef crab::domains::dis_interval <z_number> dis_z_interval_t;
     typedef crab::domains::dis_interval <q_number> dis_q_interval_t;

     template<>
     inline dis_z_interval_t trim_bound(dis_z_interval_t  x, z_number c) {

       if (x.is_bottom ())
         return x;

       dis_z_interval_t res = dis_z_interval_t::bottom ();

       if (x.is_top ()) {
         res = res | dis_z_interval_t (z_interval(c-1).lower_half_line ());
         res = res | dis_z_interval_t (z_interval(c+1).upper_half_line ());
       }
       else {
         for (auto i: boost::make_iterator_range (x.begin (), x.end ())) {
           if (i.lb() == c) {
             res = res | dis_z_interval_t (z_interval(c + 1, i.ub()));
           } else if (i.ub() == c) {
             res = res | dis_z_interval_t (z_interval(i.lb(), c - 1));
           } else {
             res = res | dis_z_interval_t (z_interval(i.lb(), c - 1));
             res = res | dis_z_interval_t (z_interval(c + 1, i.ub()));
           }
         }
       }

       return res;
     }

    template<>
    inline dis_q_interval_t trim_bound(dis_q_interval_t i, q_number /* c */) { 
      // No refinement possible for disequations over rational numbers
      return i;
    }

  } // end namespace intervals_impl
} // end namespace ikos

namespace crab {
  namespace domains {

   template<typename Number, typename VariableName>
   class dis_interval_domain: 
         public ikos::writeable, 
         public numerical_domain< Number, VariableName>,
         public bitwise_operators< Number, VariableName >, 
         public division_operators< Number, VariableName > {
     
    public:
     using typename numerical_domain< Number, VariableName>::linear_expression_t;
     using typename numerical_domain< Number, VariableName>::linear_constraint_t;
     using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
     using typename numerical_domain< Number, VariableName>::variable_t;
     using typename numerical_domain< Number, VariableName>::number_t;
     using typename numerical_domain< Number, VariableName>::varname_t;

     typedef interval <Number> interval_t;
     typedef interval_domain <Number, VariableName> interval_domain_t;

     typedef dis_interval <Number> dis_interval_t;
     typedef dis_interval_domain <Number, VariableName> dis_interval_domain_t;

    private:

     typedef separate_domain< VariableName, dis_interval_t > separate_domain_t;
     typedef linear_interval_solver< Number, VariableName, separate_domain_t > solver_t;

    public:

     typedef typename separate_domain_t::iterator iterator;

    private:

     separate_domain_t _env;

     dis_interval_domain (separate_domain_t env): _env(env) { }

    public:
          
     static dis_interval_domain_t top() {
       return dis_interval_domain (separate_domain_t::top ());
     }

     static dis_interval_domain_t bottom() {
       return dis_interval_domain (separate_domain_t::bottom ());
     }

     dis_interval_domain(): _env(separate_domain_t::top()) { }    

     dis_interval_domain (const dis_interval_domain_t& o): 
         writeable(), 
         numerical_domain< Number, VariableName >(), 
         bitwise_operators< Number, VariableName >(),
         division_operators< Number, VariableName >(),
         _env(o._env) { }
     
     dis_interval_domain_t& operator=(const dis_interval_domain_t& o) {
       if (this != &o)
         this->_env = o._env;
       return *this;
     }

     bool is_bottom() const {
       return this->_env.is_bottom();
     }
     
     bool is_top() const {
       return this->_env.is_top();
     }

     iterator begin() {
       return this->_env.begin();
     }
     
     iterator end() {
       return this->_env.end();
     }

     bool operator<=(dis_interval_domain_t e) {
       //cout << "*** Leq " << *this << " and " << e << "\n";
       bool res = this->_env <= e._env;
       //cout << "result=" << res << "\n";
       return res;
     }
     
     void operator|=(dis_interval_domain_t e) {
       this->_env = this->_env | e._env;
     }
     
     dis_interval_domain_t operator|(dis_interval_domain_t e) {
       //cout << "*** Join " << *this << " and " << e << "\n";
       dis_interval_domain_t res (this->_env | e._env);
       //cout << "result=" << res << "\n";
       return res;
     }
     
     dis_interval_domain_t operator&(dis_interval_domain_t e) {
       //cout << "*** Meet " << *this << " and " << e << "\n";
       dis_interval_domain_t res (this->_env & e._env);
       //cout << "result=" << res << "\n";
       return res;
     }

     dis_interval_domain_t operator||(dis_interval_domain_t e) {
       //cout << "*** Widening " << *this << " and " << e << "\n";
       dis_interval_domain_t res (this->_env || e._env);
       //cout << "result=" << res << "\n";
       return res;
     }
     
     template<typename Thresholds>
     dis_interval_domain_t widening_thresholds (dis_interval_domain_t e, const Thresholds &ts) {
       //cout << "*** Widening w/ thresholds " << *this << " and " << e << "\n";
       dis_interval_domain_t res = this->_env.widening_thresholds (e._env, ts);
       //cout << "result=" << res << "\n";
       return res;
     }
     
     dis_interval_domain_t operator&&(dis_interval_domain_t e) {
       //cout << "*** Narrowing " << *this << " and " << e << "\n";
       dis_interval_domain_t res (this->_env && e._env);
       //cout << "result=" << res << "\n";
       return res;
     }
     
     void operator-=(VariableName v) {
       this->_env -= v;
     }

     interval_t operator[](VariableName v)  {
       dis_interval_t x = this->_env [v];
       return x.approx ();
     }
     
     void set(VariableName v, interval_t intv) {
       this->_env.set (v, dis_interval_t (intv));
     }
     
     void operator+=(linear_constraint_system_t csts) {
       if (!this->is_bottom()) {
         //cout << "*** add constraints " << csts << " in " << *this << "\n";
         const size_t threshold = 10; // make this template parameter
         solver_t solver(csts, threshold);
         solver.run(this->_env);
         //cout << "result=" << *this << "\n";
       }
     }

     void assign (VariableName x, linear_expression_t e)  {
       //cout << "*** " <<  x << ":=" << e << " in " << *this << "\n";
       if (boost::optional<variable_t> v = e.get_variable ()) {
         this->_env.set(x, this->_env [(*v).name ()]);
       }
       else {
         dis_interval_t r (e.constant());
         for (auto t: e)
           r += dis_interval_t (t.first) * this->_env[t.second.name()];
         this->_env.set(x, r);
       }
       //cout << "result=" << *this << "\n";
     }
     
     void apply (operation_t op, VariableName x, VariableName y, Number z) {
       //cout << "*** " << x << ":=" << y << op << z << " in " << *this << "\n";
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi (z);
       dis_interval_t xi = dis_interval_t::top ();
       switch (op) {
         case OP_ADDITION:  xi = yi + zi; break; 
         case OP_SUBTRACTION: xi = yi - zi; break; 
         case OP_MULTIPLICATION:  xi = yi * zi; break; 
         case OP_DIVISION: xi = yi / zi; break;
       }
       this->_env.set(x, xi);
       //cout << "result=" << *this << "\n";
     }

     void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
       //cout << "*** " << x << ":=" << y << op << z << " in " << *this << "\n";
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi = this->_env[z];
       dis_interval_t xi = dis_interval_t::top();
       switch (op) {
         case OP_ADDITION: xi = yi + zi; break;
         case OP_SUBTRACTION: xi = yi - zi; break; 
         case OP_MULTIPLICATION:  xi = yi * zi; break; 
         case OP_DIVISION: xi = yi / zi; break; 
       }
       this->_env.set(x, xi);
       //cout << "result=" << *this << "\n";
     }
     
     
     void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width)  {
       dis_interval_t yi = this->_env[y];
       dis_interval_t xi = dis_interval_t::top ();
       switch(op){
         case OP_TRUNC: xi = yi.Trunc(width); break;
         case OP_ZEXT: xi = yi.ZExt(width); break;
         case OP_SEXT: xi = yi.SExt(width); break;
       }      
       this->_env.set(x, xi);
     }
     
     void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
       dis_interval_t yi (k);
       dis_interval_t xi = dis_interval_t::top ();
       switch(op){
         case OP_TRUNC: xi = yi.Trunc(width); break;
         case OP_ZEXT: xi = yi.ZExt(width); break;
         case OP_SEXT: xi = yi.SExt(width); break;
       }      
       this->_env.set(x, xi);
     }
     
     void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi = this->_env[z];
       dis_interval_t xi = dis_interval_t::top ();
       switch (op) {
         case OP_AND: xi = yi.And(zi); break;
         case OP_OR: xi = yi.Or(zi); break;
         case OP_XOR: xi = yi.Xor(zi); break;
         case OP_SHL: xi = yi.Shl(zi); break;
         case OP_LSHR: xi = yi.LShr(zi); break;
         case OP_ASHR: xi = yi.AShr(zi); break;
       }
       this->_env.set(x, xi);
     }
        
     void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi (k);
       dis_interval_t xi = dis_interval_t::top ();
       switch (op) {
         case OP_AND: xi = yi.And(zi); break;
         case OP_OR: xi = yi.Or(zi); break;
         case OP_XOR: xi = yi.Xor(zi); break;
         case OP_SHL: xi = yi.Shl(zi); break;
         case OP_LSHR: xi = yi.LShr(zi); break;
         case OP_ASHR: xi = yi.AShr(zi); break;
       }
       this->_env.set(x, xi);
     }
     
     void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi = this->_env[z];
       dis_interval_t xi = dis_interval_t::top ();
      
       switch (op) {
         case OP_SDIV: xi = yi / zi; break;
         case OP_UDIV: xi = yi.UDiv(zi); break;
         case OP_SREM: xi = yi.SRem(zi); break;
         case OP_UREM: xi = yi.URem(zi); break;
       }
       this->_env.set(x, xi);
     }
     
     void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
       dis_interval_t yi = this->_env[y];
       dis_interval_t zi (k);
       dis_interval_t xi = dis_interval_t::top ();
       
       switch (op) {
         case OP_SDIV: xi = yi / zi; break;
         case OP_UDIV: xi = yi.UDiv(zi); break;
         case OP_SREM: xi = yi.SRem(zi); break;
         case OP_UREM: xi = yi.URem(zi); break;
       }
       this->_env.set(x, xi);

     }

     void expand (VariableName x, VariableName new_x) {
       this->_env.set (new_x, this->_env [x]);
     }

     template<typename Range>
     void project (Range vs) {
       separate_domain_t env;
       for (auto v : vs) {
         env.set (v, this->_env [v]);
       }
       std::swap (_env, env);
     }

     void normalize () {
       if (is_bottom () || is_top ()) return;

       separate_domain_t env;
       for (auto p : boost::make_iterator_range (_env.begin (), _env.end ())) {
         env.set (p.first, p.second.normalize ());
       }
       std::swap (_env, env);
     }

     interval_domain_t approx () const {
       if (is_bottom ()) 
         return interval_domain_t::bottom ();
       else if (is_top ())
         return interval_domain_t::top ();
       else {
         interval_domain_t res;
         for (auto const p: boost::make_iterator_range (_env.begin (), _env.end ())){
           res.set (p.first, p.second.approx ());
         }
         return res;
       }
     }
     
     linear_constraint_system_t to_linear_constraint_system () {
       interval_domain_t intervals = this->approx ();
       return intervals.to_linear_constraint_system ();
     }
     
     void write(ostream& o) {
      this->_env.write(o);
     }
     
     static const char* getDomainName () {
       return "DisjunctiveIntervals";
     }  
   }; 
   
   } // namespace domains

   namespace domain_traits {
      using namespace domains;

      template <typename Number, typename VariableName>
      void normalize (dis_interval_domain<Number,VariableName>& inv) {
        inv.normalize ();
      }

      template <typename Number, typename VariableName, typename Iterator>
      void project (dis_interval_domain <Number, VariableName>& inv, 
                    Iterator begin, Iterator end) {
        inv.project (boost::make_iterator_range (begin, end));
      }
     
      template <typename Number, typename VariableName>
      void expand (dis_interval_domain <Number, VariableName>& inv, 
                   VariableName x, VariableName new_x) {
        inv.expand (x, new_x);
      }

   } // namespace domain_traits

}// namespace crab
#endif 
