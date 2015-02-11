/*******************************************************************************
 * Generic API for fixpoint iterators.
 ******************************************************************************/

#ifndef IKOS_FIXPOINT_ITERATORS_API_HPP
#define IKOS_FIXPOINT_ITERATORS_API_HPP

namespace ikos {

  template< typename NodeName, typename CFG, typename AbstractValue >
  class forward_fixpoint_iterator {

  public:
    virtual AbstractValue analyze(NodeName, AbstractValue) = 0;
    
    virtual void process_pre(NodeName, AbstractValue) = 0;
    
    virtual void process_post(NodeName, AbstractValue) = 0;
    
    virtual ~forward_fixpoint_iterator() { }
    
  }; // class fixpoint_iterator
  
} // namespace ikos

#endif // IKOS_FIXPOINT_ITERATORS_API_HPP
