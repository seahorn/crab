#ifndef BACKWARD_ANALYZER_HPP
#define BACKWARD_ANALYZER_HPP

#include <ikos/fwd_fixpoint_iterators.hpp>
#include <boost/noncopyable.hpp>

namespace analyzer 
{ 
  using namespace ikos;

  template< typename BasicBlockLabel, typename CFG, typename AbsDomain>
  class backward_fp_iterator : 
        private interleaved_fwd_fixpoint_iterator < BasicBlockLabel, CFG, AbsDomain >,
        boost::noncopyable 
  {

    typedef interleaved_fwd_fixpoint_iterator < BasicBlockLabel, CFG, AbsDomain > fwd_fixpoint_t;

    inline CFG reverse(const CFG &cfg) const
    {
      if (cfg.has_exit ())
      {
        CFG rev = cfg.clone ();
        rev.reverse();
        return rev;
      }
      else return cfg;
    }

   public:

    backward_fp_iterator(CFG cfg): fwd_fixpoint_t (reverse (cfg)) {  }  
    
    virtual ~backward_fp_iterator(){ }
    
    void run (AbsDomain inv) 
    { if (get_cfg ().has_exit ()) fwd_fixpoint_t::run(inv); }
    
    CFG get_cfg() 
    { return fwd_fixpoint_t::get_cfg(); }
    
    virtual AbsDomain analyze (BasicBlockLabel bb_id, AbsDomain post) = 0;
    
    // check of a cfg node starting from the pre
    virtual void check_pre (BasicBlockLabel bb_id, AbsDomain pre) = 0;
    
    // check of a cfg node starting from the post
    virtual void check_post (BasicBlockLabel bb_id, AbsDomain post) = 0;
 
   private:
    
    void process_pre (BasicBlockLabel  bb_id, AbsDomain inv) 
    { if (get_cfg ().has_exit ()) check_post (bb_id, inv); }

    void process_post (BasicBlockLabel  bb_id, AbsDomain inv) 
    { if (get_cfg ().has_exit ()) check_pre (bb_id, inv);  }
    
  }; // end backward_fp_iterator class

} // end namespace 

#endif /*BACKWARD_ANALYZER_HPP*/
