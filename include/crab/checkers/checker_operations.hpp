#pragma once

namespace crab {
namespace checker {  
namespace checker_operations {
    
  // Return true if inv entails cst.
  // if cst cannot be represented by Domain then return false.
  template<typename Domain>
  bool entail(Domain inv, typename Domain::linear_constraint_t cst) { 
    if (inv.is_bottom()) return true;
    if (cst.is_tautology ()) return true;
    if (cst.is_contradiction ()) return false;
    
    Domain cst_inv;
    cst_inv += cst;
    // cst cannot be represented by the domain.
    if (cst_inv.is_top ()) return false;
    
    return inv <= cst_inv; 
  }
  
  // Return true if cst intersects with inv
  template<typename Domain>
  bool intersect(typename Domain::linear_constraint_t cst, Domain inv) {
    if (inv.is_bottom () || cst.is_contradiction ()) return false;
    if (inv.is_top () || cst.is_tautology ()) return true;
    
    Domain cst_inv;
    cst_inv += cst;
    return (!(cst_inv & inv).is_bottom ());
  }
  
} 
}
}
