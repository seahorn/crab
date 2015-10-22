#ifndef LIVENESS_ANALYSIS_HPP
#define LIVENESS_ANALYSIS_HPP

/* Liveness analysis */

#include <crab/analysis/BwdAnalyzer.hpp>
#include <crab/domains/discrete_domains.hpp>

namespace crab {

  namespace analyzer {

    using namespace cfg;

    namespace liveness_discrete_impl{

       template< typename Element>
       class liveness_domain: public ikos::writeable {
         
        private:
         typedef liveness_domain< Element > liveness_domain_t;
         typedef discrete_domain< Element > discrete_domain_t;
         
        public:
         typedef typename discrete_domain_t::iterator iterator;
         
        private:
         discrete_domain_t _inv;

        public:
         // hook to speedup check_post in Liveness class
         liveness_domain(discrete_domain_t inv): 
             ikos::writeable(), _inv(inv){ }

        public:
         static liveness_domain_t top() {
           return liveness_domain(discrete_domain_t::top());
         }
         
         static liveness_domain_t bottom() {
           return liveness_domain(discrete_domain_t::bottom());
         }
         
        public:
         
         liveness_domain(): _inv(discrete_domain_t::bottom()){ }
         
         liveness_domain(Element e): ikos::writeable(), _inv(e) { }
         
         liveness_domain(const liveness_domain_t &other): 
             ikos::writeable(), _inv(other._inv) { }
         
         liveness_domain_t& operator=(liveness_domain_t other) {
           this->_inv = other._inv;
           return *this;
         }
         
         iterator begin() {
           return this->_inv.begin();
         }
         
         iterator end() {
           return this->_inv.end();
         }
         
         unsigned size() {
           return this->_inv.size();
         }
         
         bool is_bottom() {
           return this->_inv.is_bottom();
         }
         
         bool is_top() {
           return this->_inv.is_top();
         }
         
         bool operator<=(liveness_domain_t other) {
           return (this->_inv <= other._inv);
         }
         
         void operator-=(Element x) {
           this->_inv -= x;
         }
         
         void operator-=(liveness_domain_t other) {
           if (!other._inv.is_top()) {
             for ( auto v : other) 
               this->_inv -= v; 
           }
         }
         
         void operator+=(Element x) {
           this->_inv += x;
         }
         
         void operator+=(liveness_domain_t other) {
           this->_inv = (this->_inv | other._inv);
         }
         
         liveness_domain_t operator|(liveness_domain_t other) {
           return (this->_inv | other._inv);
         }
         
         liveness_domain_t operator&(liveness_domain_t other) {
           return (this->_inv & other._inv);
         }
         
         liveness_domain_t operator||(liveness_domain_t other) {
           return this->operator|(other);
         }
      
         liveness_domain_t operator&&(liveness_domain_t other) {
           return this->operator&(other);
         }
         
         void write(ostream& o) {
           this->_inv.write(o);
         }  
         
       }; 

    } // end namespace liveness_discrete_impl

#if 0
    // The implementation using discrete_domain is significantly faster
    // than the one using std::set
    namespace liveness_set_impl {
   
      template< typename Element>
      class liveness_domain {
        
        typedef liveness_domain< Element > liveness_domain_t;
        typedef std::set<Element> ElemSet;
        
       public:
        
        typedef typename ElemSet::iterator       iterator;
        typedef typename ElemSet::const_iterator const_iterator;
        
       private:
        
        bool       m_is_top;
        ElemSet m_inv;
        
        liveness_domain (ElemSet inv, bool is_top): 
            m_is_top(is_top), m_inv(inv) {
          if (m_is_top)
            m_inv.clear ();
        }
        
       public:
        
        static liveness_domain_t top() {
          return liveness_domain (ElemSet(), true);
        }
        
        static liveness_domain_t bottom() {
          return liveness_domain (ElemSet(), false);
        }
        
        liveness_domain(): m_is_top (false) { }
        
        liveness_domain(Element e): m_is_top(false)  {
          m_inv.insert(e);
        }
        
        ~liveness_domain() { 
          //m_inv.clear(); 
        }
        
        liveness_domain (const liveness_domain_t &other): 
            m_is_top(other.m_is_top), m_inv (other.m_inv) {
          if (is_top ())
            m_inv.clear ();
        }
        
        liveness_domain_t& operator=(const liveness_domain_t &other) {
          if (this != &other) {
            m_is_top = other.m_is_top;
            m_inv = other.m_inv;
          }
          return *this;
        }
        
        iterator begin()             { return m_inv.begin(); }
        iterator end()               { return m_inv.end();   }
        const_iterator begin() const { return m_inv.begin(); }
        const_iterator end()   const { return m_inv.end();   }
        
        unsigned size() const { return m_inv.size (); }
        
        bool is_top() const { return m_is_top; }
        
        bool is_bottom () const { return (!is_top() && (size() == 0)); }
        
        bool operator<=(const liveness_domain_t &other) {
          if (is_bottom() || other.is_top())    return true;
          if (is_top()    || other.is_bottom()) return false;
          
          return std::includes(m_inv.begin (), m_inv.end (), 
                               other.m_inv.begin (), other.m_inv.end ());
        }
        
        void operator-=(const liveness_domain_t &other) {
          if (is_bottom () || is_bottom ()) 
            return;
          else if (is_top ())
            CRAB_ERROR ("not defined set difference if first operand is top");
          else if (other.is_top ())  {
            *this = liveness_domain_t::bottom ();
          }
          else {
            for (auto v : other) 
              m_inv.erase (v); 
            //m_inv.erase(other.m_inv.begin (), other.m_inv.end ());
          }
        }
        
        void operator+=(const liveness_domain_t &other) {
          if (is_top()) 
            return;
          else if (is_bottom () || other.is_top())
            *this = other;
          else if (other.is_bottom ())
            return;
          else {
          m_inv.insert( other.m_inv.begin (), other.m_inv.end ());
          }
        }
        
        // Apply kill-gen dataflow equation by composing (inv -= kill ; inv += get)
        void apply (const liveness_domain_t &kill, const liveness_domain_t &gen){
          if (is_bottom () || kill.is_bottom ()) {
            this->operator+= (gen);
          }
          else if (is_top ())
            CRAB_ERROR ("not defined set difference if first operand is top");
          else if (kill.is_top ()) {
            *this = gen;
          }
          else if (gen.is_top ()) { 
            *this = liveness_domain_t::top ();
          }
          else if (gen.is_bottom ()) {
            this->operator-= (kill);
          }
          else {
            // minimize the number of deletions
            ElemSet s;
            std::set_difference(kill.m_inv.begin (), kill.m_inv.end (), 
                                gen.m_inv.begin (), gen.m_inv.end (), 
                                std::inserter (s, s.end()));
            for (auto &v : s) 
              m_inv.erase (v); 
            
            m_inv.insert( gen.m_inv.begin (), gen.m_inv.end ());
          }         
        }
        
        
        liveness_domain_t operator|(const liveness_domain_t &other) {
          
          if (is_top() || other.is_top())
            return liveness_domain_t::top();
          else if (is_bottom())
            return other;
          else if (other.is_bottom())
            return *this;
          else {
            ElemSet s(m_inv);
            s.insert(other.m_inv.begin (), other.m_inv.end ());
            return liveness_domain_t (s, false);
          }
        }
        
        liveness_domain_t operator&(const liveness_domain_t &other) {
          if (is_bottom() || other.is_bottom())
            return liveness_domain_t::bottom();
          else if (is_top()) return other;
          else if (other.is_top()) return *this;
          else {
            ElemSet s;
            std::set_intersection(m_inv.begin (), m_inv.end (), 
                                  other.m_inv.begin (), other.m_inv.end (), 
                                  std::inserter (s, s.end()));
            return liveness_domain_t (s, false);
          }
        }
        
        liveness_domain_t operator|| (const liveness_domain_t &other) {
          return this->operator|(other);
        }
        
        liveness_domain_t operator&& (const liveness_domain_t &other) {
          return this->operator&(other);
        }
        
        void write(ostream& o) {
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
        }  
      }; 
    } // end namespace liveness_set_impl
#endif 

   //! Live variable analysis
   template<typename CFG>
   class Liveness: 
        public backward_fp_iterator< typename CFG::basic_block_label_t, 
                                     CFG, 
                                     liveness_discrete_impl::liveness_domain <typename CFG::varname_t> > {

    public:
     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::varname_t varname_t;
     typedef liveness_discrete_impl::liveness_domain<varname_t> liveness_domain_t;
     typedef liveness_domain_t set_t;
     
    private:
     typedef boost::shared_ptr <set_t> set_ptr;
     typedef boost::unordered_map< basic_block_label_t, set_ptr > liveness_map_t;
     typedef std::pair< liveness_domain_t, liveness_domain_t >   kill_gen_t;
     typedef boost::unordered_map< basic_block_label_t, kill_gen_t>  kill_gen_map_t;
     typedef typename liveness_map_t::value_type l_binding_t;
     typedef typename kill_gen_map_t::value_type kg_binding_t;

     // for external queries
     liveness_map_t m_live_map;
     liveness_map_t m_dead_map;
     // for internal use
     kill_gen_map_t m_kill_gen_map;

     bool m_has_exec;

     // statistics 
     unsigned m_max_live;
     unsigned m_total_live;
     unsigned m_total_blks;

    private:

     // precompute use/def sets
     void init()
     {
       for (auto &b: boost::make_iterator_range ( this->get_cfg().begin (), 
                                                  this->get_cfg().end ())) {
         liveness_domain_t kill, gen;
         for (auto &s: b) {
           auto live = s.getLive();
           for (auto d: boost::make_iterator_range (live.defs_begin (), 
                                                    live.defs_end ())) {
             kill += d; 
             gen -= d;
           }
           for (auto u: boost::make_iterator_range (live.uses_begin (), 
                                                    live.uses_end ())) {
             gen  += u; 
           }
         }
         m_kill_gen_map.insert ( kg_binding_t (b.label (), 
                                               kill_gen_t (kill, gen)));
       }
     }
     
    public:

     Liveness (CFG cfg): 
         backward_fp_iterator
         <basic_block_label_t, CFG, liveness_domain_t> (cfg),
         m_has_exec (false),
         m_max_live (0), m_total_live (0), m_total_blks (0) {
       init(); 
     }

     void exec() { 
       if (m_has_exec) {
         CRAB_WARN ("Trying to execute liveness twice!");
         return;
       }

       this->run (liveness_domain_t::bottom()); 
       m_has_exec = true;
     }

     //! return the set of dead variables at the exit of block bb
     set_t dead_exit (basic_block_label_t bb) const {
       auto it = m_dead_map.find(bb);
       if (it == m_dead_map.end()) 
         return set_t ();
       else 
         return *(it->second); 
     }

     void get_stats (unsigned& total_live, 
                     unsigned& max_live_per_blk,
                     unsigned& avg_live_per_blk)  const {
       total_live = m_total_live;
       max_live_per_blk = m_max_live;
       avg_live_per_blk = (m_total_blks == 0 ? 0 : (int) m_total_live / m_total_blks);
     }

     // //! return the set of live variables at the entry of block bb
     // set_t live_entry (basic_block_label_t bb) const {
     //   auto it = m_live_map.find(bb);
     //   if (it == m_live_map.end())
     //     return set_t ();
     //   else
     //     return *(it->second);
     // }

     void write (ostream &o) const {
       for (auto p: m_live_map) {
         o << p.first << " :{";
         for (auto v : p.second)
           o << v << ";" ;
         o << "}\n";
       }
     }

    private:

     liveness_domain_t analyze (basic_block_label_t bb_id, 
                                liveness_domain_t inv) {
       auto it = m_kill_gen_map.find (bb_id);
       assert(it != m_kill_gen_map.end ());
       kill_gen_t p = it->second;

       inv -= p.first;
       inv += p.second;

       return inv;
     }

     void check_pre (basic_block_label_t bb, 
                     liveness_domain_t live_in) {

       // // --- Collect live variables at the entry of bb
       // if (!live_in.is_bottom()) {
       //   set_ptr live_in_set = set_ptr (new set_t (live_in));         
       //   m_live_map.insert (l_binding_t (bb, live_in_set));
       // }
     }
     
     void check_post (basic_block_label_t bb, 
                      liveness_domain_t live_out) {

       // --- Collect dead variables at the exit of bb

       if (!live_out.is_bottom ()) {
         set_ptr dead_set (new set_t (this->get_cfg().get_node (bb).live ()));
         *dead_set -= live_out;
         m_dead_map.insert (l_binding_t (bb, dead_set));


         // update statistics
         m_total_live += live_out.size ();
         m_max_live = std::max (m_max_live, live_out.size ());
         m_total_blks ++;
       }
     }
          
   }; 

   template <typename CFG>
   ostream& operator << (ostream& o, const Liveness<CFG> &l) {
     l.write (o);
     return o;
   }

  } // end namespace analyzer

} // end namespace crab

#endif 
