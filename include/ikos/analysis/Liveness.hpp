#ifndef LIVENESS_ANALYSIS_HPP
#define LIVENESS_ANALYSIS_HPP

/* Liveness analysis */

#include <ikos/analysis/BwdAnalyzer.hpp>
#include <boost/range/algorithm/copy.hpp>

// set
#include <set>
// bitset
#include <boost/dynamic_bitset.hpp>
#include <ikos/common/mergeable_map.hpp>
// sorted vector
#include <algorithm>
#include <vector>
// // discrete domain
// #include <ikos/domains/discrete_domains.hpp>

namespace analyzer 
{

  using namespace std;
  using namespace cfg;

  namespace liveness_set_impl
  {

   template< typename Element>
   class liveness_domain
   {
  
     typedef liveness_domain< Element > liveness_domain_t;
     typedef std::set<Element> ElemVector;
     
    public:

     typedef typename ElemVector::iterator       iterator;
     typedef typename ElemVector::const_iterator const_iterator;
     
    private:

     bool       m_is_top;
     ElemVector m_inv;
     
     liveness_domain (ElemVector inv, bool is_top): 
         m_is_top(is_top), m_inv(inv)
     { 
       if (m_is_top)
         m_inv.clear ();
     }
     
    public:
     
     static liveness_domain_t top() 
     {
       return liveness_domain (ElemVector(), true);
     }
     
     static liveness_domain_t bottom() 
     {
       return liveness_domain (ElemVector(), false);
     }
     
     liveness_domain(): m_is_top (false) { }
     
     liveness_domain(Element e): m_is_top(false) 
     {
       m_inv.insert(e);
     }
     
     ~liveness_domain() { m_inv.clear(); }
     
     liveness_domain (const liveness_domain_t &other): 
         m_is_top(other.m_is_top) 
     {
       if (!is_top ())
         boost::copy (other.m_inv, std::inserter(m_inv, m_inv.end ()));
     }
     
     liveness_domain_t& operator= (liveness_domain_t other) 
     {
       m_is_top = other.m_is_top;
       if (m_is_top) m_inv.clear();
       else m_inv.swap(other.m_inv); 
       return *this;
     }
     
     iterator begin()             { return m_inv.begin(); }
     iterator end()               { return m_inv.end();   }
     const_iterator begin() const { return m_inv.begin(); }
     const_iterator end()   const { return m_inv.end();   }
     
     // boost::iterator_range<iterator> get_iterator_range()
     // {
     //   return boost::make_iterator_range (begin (), end ());
     // }
     
     unsigned size() const { return std::distance (begin (), end ()); }
     
     bool is_top() { return m_is_top; }
     
     bool is_bottom () { return (!is_top() && (size() == 0)); }
     
     bool operator<= (liveness_domain_t other) 
     {
       if (is_bottom() || other.is_top())    return true;
       if (is_top()    || other.is_bottom()) return false;
       return std::includes(m_inv.begin (), m_inv.end (), 
                            other.m_inv.begin (), other.m_inv.end ());
     }
     
     void operator-= (liveness_domain_t other) 
     {
       if (size () == 0 || other.size () == 0) return;
       if (is_top () || other.is_top ()) return;
       for (auto &v : other){ m_inv.erase (v); }
       //this->_inv.erase(other._inv.begin (), other._inv.end ());
     }

     void operator+= (liveness_domain_t other) 
     {
       if (is_top()) return;
       if (other.is_top())
       {
         m_is_top = true;
         m_inv.clear();
       }
       else
         m_inv.insert( other.m_inv.begin (), other.m_inv.end ());
     }

     liveness_domain_t operator| (liveness_domain_t other) 
     {
       if (is_top() || other.is_top())
         return liveness_domain_t::top();
       
       if (is_bottom())
         return other;
       else if (other.is_bottom())
         return *this;
       else
       {
         ElemVector s(m_inv);
         s.insert(other.m_inv.begin (), other.m_inv.end ());
         return liveness_domain_t (s, false);
       }
     }
  
     liveness_domain_t operator& (liveness_domain_t other) 
     {
       if (is_bottom() || other.is_bottom())
         return liveness_domain_t::bottom();
       
       if (is_top()) return other;
       else if (other.is_top()) return *this;
       else
       {
         ElemVector s;
         std::set_intersection(m_inv.begin (), m_inv.end (), 
                               other.m_inv.begin (), other.m_inv.end (), 
                               std::inserter (s, s.end()));
         return liveness_domain_t (s, false);
       }
     }
     
     liveness_domain_t operator|| (liveness_domain_t other) 
     {
       return this->operator|(other);
     }
  
     liveness_domain_t operator&& (liveness_domain_t other) 
     {
       return this->operator&(other);
     }
     
     ostream& write(ostream& o) {
       if (is_top ())  o << "{...}";
       else if (is_bottom ())   o << "_|_";
       else 
       {
         o << "{";
         for (iterator it = begin(); it != end(); ) 
         {
           o << *it;  ++it;
           if (it != end()) 
             o << "; ";
         }
         o << "}";
       }
       return o;
     }  
   }; 
  } // end namespace

  namespace liveness_vector_impl
  {
     template< typename Element>
     class liveness_domain
     {

       template <class T, class Compare = std::less<T> >
       class sorted_vector 
       {
         
         vector<T> V;
         Compare cmp; 
         
         typedef sorted_vector <T, Compare> sorted_vector_t;
    
        public:
    
         typedef typename vector<T>::iterator iterator;
         typedef typename vector<T>::const_iterator const_iterator;
         typedef typename vector<T>::value_type value_type;
         
        public:
         
         sorted_vector(const Compare& c = Compare()) : V(), cmp(c) {}
         
         template <class InputIterator>
         sorted_vector(InputIterator first, InputIterator last, 
                       const Compare& c = Compare())
             : V(first, last), cmp(c)
         {
           std::sort(begin(), end(), cmp);
         }
         
         sorted_vector_t& operator=(sorted_vector_t other)
         {
           this->V   = other.V;
           this->cmp = other.cmp;
           return *this;
         }
         
         
         iterator begin()             { return V.begin(); }
         iterator end()               { return V.end();   }
         const_iterator begin() const { return V.begin(); }
         const_iterator end() const   { return V.end();   }
         
         void clear() { V.clear();} 
    
         unsigned size() const {return V.size();}
    
         iterator insert(const T& t) 
         {
           iterator i = lower_bound(begin(), end(), t, cmp);
           if (i == end() || cmp(t, *i))
             V.insert(i, t);
           return i;
         }
         
         // iterator insert(iterator i, const T& t)
         // {
         //   if (i == end() || cmp(t, *i))
         //     V.insert(i, t);
         //   return i;
         // }
         
         bool includes(sorted_vector_t other) const
         {
           for (auto e : *this)
             if (other.find(e) == other.end()) return false;
           return true;
         }

         sorted_vector_t intersection(sorted_vector_t other) const
         {
           sorted_vector_t res;
           for (auto e : *this)
             if (other.find(e) != other.end()){ res.insert(e); }
           return res;
         }
         
         const_iterator find(const T& t) const 
         {
           const_iterator i = lower_bound(begin(), end(), t, cmp);
           return i == end() || cmp(t, *i) ? end() : i;
         }
         
         void erase(const T& t)
         {
           if (V.empty()) return;
           V.erase(std::remove(V.begin() , V.end(), t));
         }
    
         template<class Iterator>
         void erase(Iterator first, Iterator last){
           if (first == last) return;
           Iterator new_last = last;
           for ( auto t : boost::make_iterator_range( first, last))      
             new_last = std::remove(V.begin (), V.end (), t);
           
           if (new_last != last)
             V.erase(new_last);
         }
       };
  
      private:
       typedef liveness_domain< Element > liveness_domain_t;
       typedef sorted_vector<Element> ElemVector;
       
      public:
       typedef typename ElemVector::iterator iterator;
       typedef typename ElemVector::const_iterator const_iterator;
       
      private:
       bool       m_is_top;
       ElemVector m_inv;
       
       liveness_domain(ElemVector inv, bool is_top): 
           m_is_top(is_top), m_inv(inv)
       { 
         if (m_is_top) m_inv.clear ();
       }

      public:
       static liveness_domain_t top() 
       {
         return liveness_domain(ElemVector(), true);
       }
       
       static liveness_domain_t bottom() 
       {
         return liveness_domain(ElemVector(), false);
       }
       
       liveness_domain(): m_is_top(false) { }

       liveness_domain(Element e): m_is_top(false) 
       {
         m_inv.insert(e);
       }
       
       ~liveness_domain() { m_inv.clear(); }

       liveness_domain(const liveness_domain_t &other): 
         m_is_top(other.m_is_top) 
       { 
         if (!is_top ())
           m_inv = ElemVector(other.m_inv.begin (), other.m_inv.end ());
       }
       
       liveness_domain_t& operator=(liveness_domain_t other) 
       {
         m_is_top = other.m_is_top;
         if (m_is_top) 
           m_inv.clear();
         else
           m_inv = ElemVector(other.m_inv.begin (), other.m_inv.end ());
         return *this;
       }
       
       iterator begin()             { return m_inv.begin(); }
       iterator end()               { return m_inv.end();   }
       const_iterator begin() const { return m_inv.begin(); }
       const_iterator end()   const { return m_inv.end();   }
       
       // boost::iterator_range<iterator> get_iterator_range()
       // {
       //   return boost::make_iterator_range(begin(), end());
       // }

       unsigned size() const { return m_inv.size();} 

       bool is_top() { return m_is_top; }

       bool is_bottom() { return (!is_top() && (size() == 0)); }
  
       bool operator<=(liveness_domain_t other) 
       {
         if (is_bottom() || other.is_top())    return true;
         if (is_top()    || other.is_bottom()) return false;
         return m_inv.includes(other.m_inv);
       }
       
       void operator-=(liveness_domain_t other) 
       {
         if (size () == 0 || other.size () == 0)
           return;
         if (is_top () || other.is_top ())
           return;
         m_inv.erase(other.m_inv.begin (), other.m_inv.end ());
       }

       void operator+=(liveness_domain_t other) 
       {
         if (is_top()) return;
         if (other.is_top())
         {
           m_is_top = true;
           m_inv.clear();
         }
         else
         {
           for (auto e: boost::make_iterator_range( other.m_inv.begin (), 
                                                    other.m_inv.end ()))        
             m_inv.insert (e);
         }
       }
       
       liveness_domain_t operator|(liveness_domain_t other) 
       {
         if (is_top() || other.is_top())
           return liveness_domain_t::top();
         
         if (is_bottom())
           return other;
         else if (other.is_bottom())
           return *this;
         else
         {
           ElemVector s(m_inv);
           for (auto e: boost::make_iterator_range( other.m_inv.begin (), 
                                                    other.m_inv.end ()))       
             s.insert (e);
           return liveness_domain_t (s, false);
         }
       }
  
       liveness_domain_t operator&(liveness_domain_t other) 
       {
         if (is_bottom() || other.is_bottom())
           return liveness_domain_t::bottom();
         
         if (is_top())
           return other;
         else if (other.is_top())
           return *this;
         else
         {
           ElemVector s = m_inv.intersection(other.m_inv);
           return liveness_domain_t (s, false);
         }
       }

       liveness_domain_t operator||(liveness_domain_t other) 
       {
         return this->operator|(other);
       }
       
       liveness_domain_t operator&&(liveness_domain_t other) 
       {
         return this->operator&(other);
       }
       
       ostream& write(ostream& o) {
         if (is_top ()) 
           o << "{...}";
         else if (is_bottom ()) 
           o << "_|_";
         else 
         {
           o << "{";
           for (iterator it = begin(); it != end(); ) 
           {
             o << *it; 
             ++it;
             if (it != end()) 
               o << "; ";
           }
           o << "}";
         }
         return o;
       }  
     }; 
  } // end namespace

  namespace liveness_bitset_impl{

  template< typename Element, unsigned N = 50>
  class liveness_domain
  {

    template< typename T>
    class BitSet
    {
      
      // Wrapper for keys in patricia trees
      template< typename X>
      class Index  {    
        index_t   _idx;  
       public:
        Index(X n): _idx(n){ }
        Index(const Index &o): _idx(o._idx) { }
        Index& operator=(Index o){
          this->_idx = o._idx;
          return *this;
        }
        index_t index() { return _idx; }
        index_t index() const{ return _idx; }
        bool operator<(const Index&o) const { return (_idx < o._idx); }
        std::ostream& write(std::ostream& o) {
          o << _idx;
          return o;
        }
      }; // end class index
      
      typedef BitSet< T > BitSet_t;
      typedef unsigned bit_t;
      typedef mergeable_map< Index<bit_t> , T> bit_to_elem_map_t;
      boost::dynamic_bitset<> _bitset;
      bit_to_elem_map_t _map;
      
      void resize(boost::dynamic_bitset<> &b, 
                  boost::dynamic_bitset<>::size_type out_idx)
      {
        boost::dynamic_bitset<>::size_type curr_sz = b.size ();    
        assert(out_idx >= curr_sz);
        boost::dynamic_bitset<>::size_type delta = out_idx - curr_sz;
        b.resize( (curr_sz + delta) * 2);
      }
      
      void unify_sizes(boost::dynamic_bitset<> &b1, 
                       boost::dynamic_bitset<> &b2) 
      {
        boost::dynamic_bitset<>::size_type s1 = b1.size ();    
        boost::dynamic_bitset<>::size_type s2 = b2.size ();    
        if (s1 == s2) 
          return;
        
        if (s1 > s2)
          b2.resize(s1);
        else
          b1.resize(s2);
        
        assert(b1.size () == b2.size ());
      }
      
     private:
      
      BitSet(boost::dynamic_bitset<> bitset, bit_to_elem_map_t map) : 
          _bitset(bitset), _map(map) { }
      
     public:
      
      BitSet(unsigned sz): _bitset(sz) { }
    
      BitSet(const BitSet_t &other): 
        _bitset(other._bitset), _map(other._map) { }
      
      BitSet_t& operator= (BitSet_t other)
      { 
        this->_bitset = other._bitset;
        this->_map    = other._map;
        return *this;
      }
      
      void clear()
      {
        this->_bitset.clear();
        this->_map.clear();
      }
      
      void set (const T &x, bool val)
      {
        if (x.index () >= this->_bitset.size())
          resize(this->_bitset, x.index ());
        
        this->_bitset.set ( x.index (), val); 
        this->_map.set ( x.index (), x);
      }
      
      bool operator[](const T &x) const
      {
        return this->_bitset[x.index ()];
      }
      
      unsigned size() const { return this->_bitset.count();}

      // set inclusion
      bool operator<=(const BitSet_t &other) 
      {
        boost::dynamic_bitset<> other_bitset(other._bitset);
        boost::dynamic_bitset<> bitset(this->_bitset);
        unify_sizes(bitset, other_bitset);
        return ( (bitset |= other_bitset) == other_bitset);
      }
      
      // set union
      BitSet_t operator|(BitSet_t other)
      {
        boost::dynamic_bitset<> other_bitset = other._bitset;
        boost::dynamic_bitset<> bitset(this->_bitset);
        unify_sizes(bitset, other_bitset);
        bitset |= other_bitset;
        return BitSet_t(bitset, this->_map | other._map);
      }
      
      // set intersection
      BitSet_t operator&(BitSet_t other)
      {
        boost::dynamic_bitset<> other_bitset = other._bitset;
        boost::dynamic_bitset<> bitset(this->_bitset);
        unify_sizes(bitset, other_bitset);
        bitset &= other_bitset;
        return BitSet_t(bitset, this->_map | other._map);
      }
      
      // set difference
      BitSet_t& operator-=(const BitSet_t &other)
      {
        boost::dynamic_bitset<> other_bitset(other._bitset);
        unify_sizes(this->_bitset, other_bitset);
        this->_bitset -= other_bitset;
        return *this;
      }

      std::set<T> to_set() 
      {
        std::set<T> s;
        for (boost::dynamic_bitset<>::size_type i=0; i < this->_bitset.size(); i++)
        {
          if (!this->_bitset[i]) continue;
          boost::optional<T> v = this->_map[i];
          assert (v && "bit cannot mapped to a variable name");
          s.insert(*v);
        }
        return s;
      }
      
      ostream& write(ostream &o)
      {
        if (size())
          o << "_|_" ;
        else
          o << "{"; for (auto &x : to_set()) {  o << x << ";" ; } o << "}";
        return o;
      }
      
    }; //end class BitSet
  
   public:
    typedef typename std::set<Element>::iterator       iterator;
    typedef typename std::set<Element>::const_iterator const_iterator;
    
   private:
    typedef liveness_domain< Element, N > liveness_domain_t;
    typedef BitSet<Element> BitSet_t;
    
    bool     m_is_top;
    BitSet_t m_inv;

    liveness_domain(BitSet_t inv, bool is_top): 
        m_is_top(is_top), m_inv(inv)
    { 
      if (m_is_top) m_inv.clear ();
    }

   public:
    static liveness_domain_t top() 
    {
      return liveness_domain(BitSet_t(N), true);
    }
  
    static liveness_domain_t bottom() 
    {
      return liveness_domain(BitSet_t(N), false);
    }
  
    liveness_domain(): m_is_top(false), m_inv(BitSet_t(N)) { }

    liveness_domain(Element e): m_is_top(false), m_inv(BitSet_t(N)) 
    {
      m_inv.set(e, true);
    }

    ~liveness_domain() {  m_inv.clear(); }

    // boost::iterator_range<iterator> get_iterator_range()
    // {
    //   auto s = this->_inv.to_set();
    //   return boost::make_iterator_range(s.begin(), s.end());
    // }
    
    // expensive operation
    std::set<Element> to_set (){ return m_inv.to_set(); }
    // expensive operations
    iterator begin ()       { return to_set ().begin (); }
    iterator end ()         { return to_set ().end (); }
    iterator begin () const { return to_set ().begin (); }
    iterator end ()   const { return to_set ().end (); }
    
    liveness_domain(const liveness_domain_t &other): 
      m_is_top(other.m_is_top), m_inv(other.m_inv) { } 

    liveness_domain_t& operator=(liveness_domain_t other) 
    {
      m_is_top = other.m_is_top;
      m_inv    = other.m_inv;
      if (m_is_top) m_inv.clear();
      return *this;
    }
  
    unsigned size() const { return m_inv.size(); }
    
    bool is_top() { return m_is_top; }
    
    bool is_bottom() { return (!is_top() && (size() == 0)); }
    
    bool operator<=(liveness_domain_t other) 
    {
      if (is_bottom() || other.is_top())    return true;
      if (is_top()    || other.is_bottom()) return false;
      return m_inv <= other.m_inv;
    }

    void operator-=(liveness_domain_t other) 
    {
      if (size () == 0 || other.size () == 0)
        return;
      if (is_top () || other.is_top ())
        return;
      m_inv -= other.m_inv;
    }

    void operator+=(liveness_domain_t other) 
    {
      if (is_top()) return;
      if (other.is_top())
      {
        m_is_top = true;
        m_inv.clear();
      }
      else
        m_inv = m_inv | other.m_inv;
    }

    liveness_domain_t operator|(liveness_domain_t other) 
    {
      if (is_top() || other.is_top())
        return liveness_domain_t::top();
      
      if (is_bottom())
        return other;
      else if (other.is_bottom())
        return *this;
      else
        return liveness_domain_t (m_inv | other.m_inv , false);
    }
  
    liveness_domain_t operator&(liveness_domain_t other) 
    {
      if (is_bottom() || other.is_bottom())
        return liveness_domain_t::bottom();
      
      if (is_top())
        return other;
      else if (other.is_top())
        return *this;
      else
        return liveness_domain_t (m_inv & other.m_inv , false);
    }

    liveness_domain_t operator||(liveness_domain_t other) 
    {
      return this->operator|(other);
    }
  
    liveness_domain_t operator&&(liveness_domain_t other) 
    {
      return this->operator&(other);
    }
    
    ostream& write(ostream& o) 
    {
      if (is_top ()) 
        o << "{...}";
      else if (is_bottom ()) 
        o << "_|_";
      else 
        m_inv.write(o);
      return o;
    }  
  }; 
  } // end namespace

// namespace liveness_discrete_impl{

// template< typename Element>
// class liveness_domain: public ikos::writeable {
  
//  private:
//   typedef liveness_domain< Element > liveness_domain_t;

//  private:
//   typedef discrete_domain< Element > discrete_domain_t;
  
//  public:
//   typedef typename discrete_domain_t::iterator iterator;
  
//  private:
//   discrete_domain_t _inv;

//   liveness_domain(discrete_domain_t inv): 
//      ikos::writeable(), _inv(inv){ }

//  public:
//   static liveness_domain_t top() {
//     return liveness_domain(discrete_domain_t::top());
//   }
  
//   static liveness_domain_t bottom() {
//     return liveness_domain(discrete_domain_t::bottom());
//   }
  
//  public:

//   liveness_domain(): _inv(discrete_domain_t::bottom()){ }

//   liveness_domain(Element e): ikos::writeable(), _inv(e) { }

//   liveness_domain(const liveness_domain_t &other): 
//       ikos::writeable(), _inv(other._inv) { }

//   liveness_domain_t& operator=(liveness_domain_t other) {
//     this->_inv = other._inv;
//     return *this;
//   }
  
//   iterator begin() {
//     return this->_inv.begin();
//   }

//   iterator end() {
//     return this->_inv.end();
//   }

//   boost::iterator_range<iterator> get_iterator_range()
//   {
//     return boost::make_iterator_range(begin(), end());
//   }

//   unsigned size() {
//     return this->_inv.size();
//   }

//   bool is_bottom() {
//     return this->_inv.is_bottom();
//   }

//   bool is_top() {
//     return this->_inv.is_top();
//   }
  
//   bool operator<=(liveness_domain_t other) {
//     return (this->_inv <= other._inv);
//   }

//   void operator-=(Element x) {
//     this->_inv -= x;
//   }

//   void operator-=(liveness_domain_t other) {
//     if (!other._inv.is_top())
//     {
//       for ( auto v : other) { this->_inv -= v; }
//     }
//   }

//   void operator+=(Element x) {
//     this->_inv += x;
//   }

//   void operator+=(liveness_domain_t other) {
//     this->_inv = (this->_inv | other._inv);
//   }

//   liveness_domain_t operator|(liveness_domain_t other) {
//     return (this->_inv | other._inv);
//   }
  
//   liveness_domain_t operator&(liveness_domain_t other) {
//     return (this->_inv & other._inv);
//   }

//   liveness_domain_t operator||(liveness_domain_t other) {
//     return this->operator|(other);
//   }
  
//   liveness_domain_t operator&&(liveness_domain_t other) {
//     return this->operator&(other);
//   }

//   ostream& write(ostream& o) {
//     return this->_inv.write(o);
//   }  

// }; 
// } // end namespace

   //! Live variable analysis
   template<typename CFG>
   class Liveness: 
      public backward_fp_iterator< typename CFG::basic_block_label_t, 
                                   CFG, 
                                   liveness_set_impl::liveness_domain <typename CFG::varname_t> >
   {
    public:
     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::varname_t varname_t;
     typedef std::set<varname_t> live_set_t;
     
    private:
     typedef liveness_set_impl::liveness_domain<varname_t> liveness_domain_t;
     typedef boost::unordered_map< basic_block_label_t, live_set_t > liveness_map_t;
     typedef std::pair< liveness_domain_t, liveness_domain_t >   kill_gen_t;
     typedef boost::unordered_map< basic_block_label_t, kill_gen_t>  kill_gen_map_t;
     typedef typename liveness_map_t::value_type l_binding_t;
     typedef typename kill_gen_map_t::value_type kg_binding_t;

     // for external queries
     liveness_map_t m_live_map;
     liveness_map_t m_dead_map;
     // for internal use
     kill_gen_map_t m_kill_gen_map;
     live_set_t m_all_vars;

    private:

     // precompute use/def sets
     void init()
     {
       for (auto &b: boost::make_iterator_range ( this->get_cfg().begin (), 
                                                  this->get_cfg().end ()))
       {
         liveness_domain_t kill, gen;
         for (auto &s: b)
         {
           auto live = s.getLive();
           for (auto d: boost::make_iterator_range (live.defs_begin (), live.defs_end ()))
           {
             kill += d; 
             gen -= d;
             m_all_vars.insert (d);
           }
           for (auto u: boost::make_iterator_range (live.uses_begin (), live.uses_end ()))
           {
             gen  += u; 
             m_all_vars.insert (u);
           }
         }
         m_kill_gen_map.insert ( kg_binding_t (b.label (), kill_gen_t (kill, gen)));
       }
     }
     
    public:

     Liveness (CFG cfg): 
        backward_fp_iterator<basic_block_label_t, CFG, liveness_domain_t> (cfg) 
     { init(); }

     void exec() { this->run (liveness_domain_t::bottom());  }

     //! return the set of dead variables at the exit of block bb
     live_set_t dead_exit (basic_block_label_t bb) const
     {
       auto it = m_dead_map.find(bb);
       if (it == m_dead_map.end()) return live_set_t ();
       else return it->second; 
     }

     //! return the set of live variables at the entry of block bb
     boost::optional <live_set_t> live_entry (basic_block_label_t bb) const
     {
       auto it = m_live_map.find(bb);
       if (it == m_live_map.end())
         return boost::optional <live_set_t> ();
       else
         return boost::optional <live_set_t> (it->second); 
     }

     ostream& write (ostream &o) const
     {
       for (auto p: m_live_map)
       {
         o << p.first << " :{";
         for (auto v : p.second){ o << v << ";" ;  }
         o << "}\n";
       }
       return o;
     }

    private:

     liveness_domain_t analyze (basic_block_label_t bb_id, liveness_domain_t inv)
     {
       auto it = m_kill_gen_map.find (bb_id);
       assert(it != m_kill_gen_map.end ());
       kill_gen_t p = it->second;
       inv -= p.first;
       inv += p.second;
       return inv;
     }

     void check_pre (basic_block_label_t bb, liveness_domain_t pre)
     {
       // Collect live variables at the entry of bb
       live_set_t live_set;
       if (!pre.is_bottom())
       {
         for (auto v: boost::make_iterator_range (pre.begin (), pre.end ()))
         { live_set.insert (v); }
       }
       m_live_map.insert (l_binding_t (bb, live_set));
     }
     
     void check_post (basic_block_label_t bb, liveness_domain_t post)
     { 
       // Collect dead variables at the exit of bb
       live_set_t dead_set;
       if (!post.is_bottom ())
       { 
         live_set_t live_set;
         for (auto v: boost::make_iterator_range (post.begin (), post.end ()))
         { live_set.insert (v); }
         
         std::set_difference( m_all_vars.begin (), m_all_vars.end (), 
                              live_set.begin   (), live_set.end (), 
                              std::inserter (dead_set, dead_set.end()));
       }
       m_dead_map.insert (l_binding_t (bb, dead_set));
     }
     
   }; 

   template <typename CFG>
   ostream& operator << (ostream& o, const Liveness<CFG> &l)
   {
     l.write (o);
     return o;
   }

} // end namespace 

#endif 
