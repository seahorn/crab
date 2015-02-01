/*******************************************************************************
 * Generic API for control-flow graphs.
 ******************************************************************************/


#ifndef IKOS_CFG_API_HPP
#define IKOS_CFG_API_HPP

#include <ikos/collections.hpp>

namespace ikos {

  template< typename NodeName, typename Node >
  class cfg {
    
  public:
    typedef collection< NodeName > node_collection_t;

  public:
    virtual NodeName entry() = 0;

    virtual Node get_node(NodeName) = 0;

    virtual node_collection_t next_nodes(NodeName) = 0;

    virtual node_collection_t prev_nodes(NodeName) = 0;

    virtual ~cfg() { }

  }; // class cfg
  
} // namespace ikos

#endif // IKOS_CFG_API_HPP
