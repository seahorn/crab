#define LDD_DEBUG 1

/**
   C++ interface for LDD library kindly provided by Ufo
 */

#pragma once

#include "cuddInt.h"
#include "ldd.h"
#undef TRUE
#undef FALSE
#include "tvpi.h"

// ugly but it disables cudd macro that might conflict with boost::lexical_cast
#ifdef fail
#undef fail
#endif

#include <functional>
#include <memory>

namespace crab {
namespace domains {
namespace ldd {

/** An LDD node */
typedef std::shared_ptr<LddNode> LddNodePtr;
typedef std::pair<LddNodePtr, LddNodePtr> LddNodePtrPair;
typedef std::vector<LddNodePtr> LddNodePtrVector;

/** Custom deleter for LddNodePtr*/
class ldd_node_deleter {
private:
  LddManager *ldd;

public:
  ldd_node_deleter(LddManager *m) : ldd(m) {}
  void operator()(LddNode *p) {
    if (p)
      Ldd_RecursiveDeref(ldd, p);
  }
  LddManager *getManager() const { return ldd; }
};

/** Constructs a LddNodePtr from an LddManager and LddNode */
inline LddNodePtr lddPtr(LddManager *ldd, LddNode *n) {
  // -- if n is not NULL, inc the reference
  if (n) {
    Ldd_Ref(n);
    // -- wrap in a shared pointer interface
    LddNodePtr p(n, ldd_node_deleter(ldd));
    return p;
  }
#ifdef LDD_DEBUG
  Cudd_ReduceHeap(Ldd_GetCudd(ldd), CUDD_REORDER_SAME, 0);
  assert(Cudd_DebugCheck(Ldd_GetCudd(ldd)) == 0);
  assert(Cudd_CheckKeys(Ldd_GetCudd(ldd)) == 0);
  Ldd_SanityCheck(ldd);
#endif
  return LddNodePtr();
}

// XXX:: get_deleter requires RTTI which we try to avoid.
// We comment getLddManager since it's dead code.
/** extract a manager from an LddNodePtr*/
// inline LddManager* getLddManager (LddNodePtr n) {
//   if (ldd_node_deleter const* pd = std::get_deleter<ldd_node_deleter>(n))
//     return pd->getManager ();
//   return NULL;
// }

inline size_t hash_value(LddNodePtr n) {
  return (size_t)(((size_t)(&*n)) * DD_P1);
}

template <typename Number> theory_t *create_box_theory(size_t sz);

template <> inline theory_t *create_box_theory<ikos::z_number>(size_t sz) {
  return tvpi_create_boxz_theory(sz);
}

template <> inline theory_t *create_box_theory<ikos::q_number>(size_t sz) {
  return tvpi_create_box_theory(sz);
}

} // namespace ldd
} // namespace domains
} // namespace crab

namespace std {
template <> struct hash<crab::domains::ldd::LddNodePtr> {
  std::size_t operator()(const crab::domains::ldd::LddNodePtr &v) const {
    return crab::domains::ldd::hash_value(v);
  }
};
} // namespace std
