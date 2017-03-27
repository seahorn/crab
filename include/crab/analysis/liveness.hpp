#ifndef LIVENESS_ANALYSIS_HPP
#define LIVENESS_ANALYSIS_HPP

/* Liveness analysis */

#include <crab/common/stats.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/types.hpp>
#include <boost/noncopyable.hpp>

#if 1

#include <crab/iterators/killgen_fixpoint_iterator.hpp>

namespace crab {

  namespace analyzer {

   template<typename V>
   using liveness_domain = crab::domains::flat_killgen_domain<V>;       

    
   template<class CFG>
   class liveness_ops:
      public crab::iterators::
      killgen_operations_api<CFG, liveness_domain<typename CFG::varname_t> > {

   public:
     
     typedef liveness_domain<typename CFG::varname_t> liveness_domain_t;
     typedef typename CFG::basic_block_label_t basic_block_label_t;

   private:
     
     typedef crab::iterators::killgen_operations_api<CFG, liveness_domain_t> this_type;
     typedef std::pair<liveness_domain_t, liveness_domain_t> binding_t;
     typedef boost::unordered_map<basic_block_label_t, binding_t> liveness_map_t;
     
     liveness_map_t _liveness_map;

    public:

     liveness_ops (CFG cfg): this_type (cfg) { }

     virtual bool is_forward () { return false; }

     virtual liveness_domain_t entry() 
     { return liveness_domain_t::bottom();}

     virtual liveness_domain_t merge(liveness_domain_t d1,  liveness_domain_t d2)
     { return d1 | d2; }

     virtual void init_fixpoint() {
        for (auto &b: boost::make_iterator_range(this->_cfg.begin(),this->_cfg.end())) {
          liveness_domain_t kill, gen;
          for (auto &s: boost::make_iterator_range(b.rbegin(),b.rend())) { 
             auto live = s.get_live();
             for (auto d: boost::make_iterator_range(live.defs_begin(), 
                                                     live.defs_end())) {
               kill += d; 
               gen -= d;
             }
             for (auto u: boost::make_iterator_range(live.uses_begin(), 
                                                     live.uses_end())) {
               gen  += u; 
             }
          }
          _liveness_map.insert(std::make_pair(b.label(), binding_t(kill, gen)));
        }
      }

      virtual liveness_domain_t analyze (basic_block_label_t bb_id, liveness_domain_t in) {
        auto it = _liveness_map.find (bb_id);
        assert(it != _liveness_map.end ());
        in -= it->second.first;
        in += it->second.second;
        return in;
      }

      virtual std::string name () { return "liveness";}
   };

   //! Live variable analysis
   template<typename CFG>
   class liveness:
      boost::noncopyable, 
      public crab::iterators::killgen_fixpoint_iterator <CFG, liveness_ops<CFG> >{
     
    public:

     typedef typename CFG::basic_block_label_t basic_block_label_t;
     typedef typename CFG::statement_t statement_t;
     typedef typename CFG::varname_t varname_t;
     
    private:

     typedef liveness_ops<CFG> liveness_ops_t;
     typedef typename liveness_ops_t::liveness_domain_t liveness_domain_t;
     typedef crab::iterators::killgen_fixpoint_iterator<CFG,liveness_ops_t> 
             killgen_fixpoint_iterator_t;
     
    public:

     //for backward compatibility
     //XXX: maybe unused already
     typedef liveness_domain_t set_t;

    private:

     // output of the analysis: map basic blocks to set of dead
     // variables at the end of the blocks
     boost::unordered_map<basic_block_label_t, liveness_domain_t> _dead_map;

     // statistics 
     unsigned _max_live;
     unsigned _total_live;
     unsigned _total_blks;

     void process_post (basic_block_label_t bb, 
                        liveness_domain_t live_out) {
       // --- Collect dead variables at the exit of bb
       if (!live_out.is_bottom ()) {
         liveness_domain_t dead_set = this->_cfg.get_node(bb).live();
         dead_set -= live_out;
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << dead_set <<"\n";);
         _dead_map.insert (std::make_pair(bb, dead_set));
         // update statistics
         _total_live += live_out.size ();
         _max_live = std::max (_max_live, live_out.size ());
         _total_blks ++;
       } else {
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << live_out <<"\n";);
       }
     }
          
    public:

     liveness (CFG cfg)
         : killgen_fixpoint_iterator_t(cfg), 
           _max_live (0), _total_live (0), _total_blks (0) { }

     void exec() { 
       this->run();
       for (auto p : boost::make_iterator_range(this->out_begin(), this->out_end()))
       { process_post (p.first, p.second); } 
       this->release_memory();
     }

     //! return the set of dead variables at the exit of block bb
     liveness_domain_t dead_exit (basic_block_label_t bb) const {
       auto it = _dead_map.find(bb);
       if (it == _dead_map.end()) 
         return liveness_domain_t();
       else 
         return it->second; 
     }
       
     void get_stats (unsigned& total_live, 
                     unsigned& max_live_per_blk,
                     unsigned& avg_live_per_blk)  const {
       total_live = _total_live;
       max_live_per_blk = _max_live;
       avg_live_per_blk = (_total_blks == 0 ? 0 : (int) _total_live / _total_blks);
     }

     // TODO
     void write (crab_os &o) const { }         
   }; 

   template <typename CFG>
   crab_os& operator << (crab_os& o, const liveness<CFG> &l) {
     l.write (o);
     return o;
   }

  } // end namespace analyzer
} // end namespace crab
#else

#include <crab/domains/discrete_domains.hpp>
#include <crab/cfg/cfg_bgl.hpp>
#include <crab/analysis/graphs/sccg.hpp>
#include <crab/analysis/graphs/topo_order.hpp>

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
         // hook to speedup check_post in liveness class
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
         
         liveness_domain_t& operator=(const liveness_domain_t &other) {
           if (this != &other) 
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
           if (is_bottom ()) 
             return true;
           else if (other.is_top ())
             return true;
           else
             return (this->_inv <= other._inv);
         }
         
         void operator-=(Element x) {
           if (is_bottom ()) 
             return;
           
           this->_inv -= x;
         }
         
         void operator-=(liveness_domain_t other) {
           if (is_bottom () || other.is_bottom ()) 
             return;

           if (!other._inv.is_top()) {
             for ( auto v : other) 
               this->_inv -= v; 
           }
         }
         
         void operator+=(Element x) {
           if (is_top ()) 
             return;

           this->_inv += x;
         }
         
         void operator+=(liveness_domain_t other) {
           if (is_top () || other.is_bottom ()) {
             return;
           }
           else if (other.is_top ()) {
             _inv = discrete_domain_t::top();
           }
           else {
             this->_inv = (this->_inv | other._inv);
           }
         }
         
         liveness_domain_t operator|(liveness_domain_t other) {
           return (this->_inv | other._inv);
         }
         
         liveness_domain_t operator&(liveness_domain_t other) {
           return (this->_inv & other._inv);
         }
         
         liveness_domain_t operator||(liveness_domain_t other) {
           // liveness domain is finite so no infinite AC
           return other;
         }
      
         liveness_domain_t operator&&(liveness_domain_t other) {
           // liveness domain is finite so no infinite DC
           return other;
         }
         
         void write(crab_os& o) {
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
        
        bool m_is_top;
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
        
        void write(crab_os& o) {
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
   class liveness: boost::noncopyable {

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

     CFG m_cfg;
     boost::unordered_map< basic_block_label_t, set_t> m_in_map;
     boost::unordered_map< basic_block_label_t, set_t> m_out_map;

     // for internal use
     kill_gen_map_t m_kill_gen_map;
     // for external queries
     liveness_map_t m_dead_map;

     // avoid run twice
     bool m_has_exec;

     // statistics 
     unsigned m_max_live;
     unsigned m_total_live;
     unsigned m_total_blks;

    private:

     // precompute use/def sets
     void init() {
       for (auto &b: boost::make_iterator_range(m_cfg.begin(),m_cfg.end())) {
         
         liveness_domain_t kill, gen;
         for (auto &s: boost::make_iterator_range(b.rbegin(),b.rend())) { 
                                                   
           auto live = s.get_live();
           for (auto d: boost::make_iterator_range(live.defs_begin(), 
                                                   live.defs_end())) {
             kill += d; 
             gen -= d;
           }
           for (auto u: boost::make_iterator_range(live.uses_begin(), 
                                                   live.uses_end())) {
             gen  += u; 
           }
         }
         m_kill_gen_map.insert ( kg_binding_t (b.label (), 
                                               kill_gen_t (kill, gen)));
       }
     }
     
    public:

     liveness (CFG cfg)
         : m_cfg (cfg), m_has_exec (false),
           m_max_live (0), m_total_live (0), m_total_blks (0) {
       crab::ScopedCrabStats __st__("liveness");
       init(); 
     }

     void exec() { 
        crab::ScopedCrabStats __st__("liveness");

       if (m_has_exec) {
         CRAB_WARN ("Trying to execute liveness twice!");
         return;
       }

        std::vector<typename CFG::node_t> order = 
            crab::analyzer::graph_algo::weak_rev_topo_sort(m_cfg);

       assert (order.size () == std::distance(m_cfg.begin(), m_cfg.end()));

       CRAB_LOG("liveness", 
                crab::outs()  << "liveness fixpoint ordering={"; 
                for (auto &v : order)
                  crab::outs() << cfg_impl::get_label_str(v) << " -- "; 
                crab::outs() << "}\n";); 

       bool change = true;
       unsigned iterations = 0;
       while (change) {
         change = false;
         ++iterations;
         for (auto &n : order) {
           liveness_domain_t out = liveness_domain_t::bottom();
           for (auto p: m_cfg.next_nodes (n))
             out = out | m_in_map [p]; 
           liveness_domain_t in = analyze (n, out);
           liveness_domain_t old_in = m_in_map [n];
           if (!(in <= old_in)) {
             m_in_map [n] = in | old_in;
             change = true;
           }
           else
             m_out_map [n] = out;
         }
       }
       
       for (auto &n : order) {
         process_post (n, m_out_map [n]);
       }

      CRAB_LOG("liveness", 
               crab::outs() << "liveness fixpoint reached in " 
                         << iterations << " iterations \n"); 
       m_has_exec = true;
       
      CRAB_LOG("liveness", 
               crab::outs() << "liveness sets: \n";
               for (auto n: boost::make_iterator_range (m_cfg.label_begin (),
                                                        m_cfg.label_end ())) {
                 crab::outs() << cfg_impl::get_label_str(n) << " "
                              << "OUT=" << m_out_map [n]  << " "
                              << "IN=" << m_in_map [n] << "\n";
               }
               crab::outs() << "\n";);
       
       // --- Keep a small memory footprint in client analyses
       m_in_map.clear ();
       m_out_map.clear ();
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

     void write (crab_os &o) const {
     }
     
    private:

      
     liveness_domain_t analyze (basic_block_label_t bb_id, 
                                liveness_domain_t live_out) {
       auto it = m_kill_gen_map.find (bb_id);
       assert(it != m_kill_gen_map.end ());

       // --- compute live_in
       live_out -= it->second.first;
       live_out += it->second.second;
       return live_out;
     }

     void process_post (basic_block_label_t bb, 
                        liveness_domain_t live_out) {

       // --- Collect dead variables at the exit of bb
       if (!live_out.is_bottom ()) {
         auto dead_set = boost::make_shared<set_t>(m_cfg.get_node (bb).live ());
         *dead_set -= live_out;
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << *dead_set <<"\n";);
         m_dead_map.insert (l_binding_t (bb, dead_set));
         // update statistics
         m_total_live += live_out.size ();
         m_max_live = std::max (m_max_live, live_out.size ());
         m_total_blks ++;
       } else {
         CRAB_LOG("liveness",
                  crab::outs() << cfg_impl::get_label_str(bb) 
                               << " dead variables=" << live_out <<"\n";);
       }
     }
          
   }; 

   template <typename CFG>
   crab_os& operator << (crab_os& o, const liveness<CFG> &l) {
     l.write (o);
     return o;
   }
  } // end namespace analyzer
} // end namespace crab
#endif 

#endif 
