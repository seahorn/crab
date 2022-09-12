#pragma once

#include <type_traits>

#include <crab/domains/graphs/adapt_sgraph.hpp>
#include <crab/domains/graphs/ht_graph.hpp>
#include <crab/domains/graphs/pt_graph.hpp>
#include <crab/domains/graphs/sparse_graph.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/numbers/safeint.hpp>

namespace crab {
namespace domains {
namespace DBM_impl {

// All of these representations are implementations of a
// sparse weighted graph. They differ on the datastructures
// used to store successors and predecessors
enum GraphRep {
  // sparse-map and sparse-sets
  ss = 1,
  // adaptive sparse-map and sparse-sets
  adapt_ss = 2,
  // patricia tree-maps and patricia tree-sets
  pt = 3,
  // hash table and hash sets
  ht = 4
};

/** DBM weights (Wt) can be represented using one of the following
 * types:
 *
 * 1) basic integer type: e.g., int64_t
 * 2) safei64
 * 3) z_number
 *
 * 1) is the fastest but things can go wrong if some DBM
 * operation overflows.
 * 
 * 2) is slower than 1) but it checks for overflow before any DBM
 * operation.
 * 
 * 3) is the slowest and it represents weights using unbounded
 * mathematical integers so overflow is not a concern but it might not
 * be what you need when reasoning about programs with wraparound
 * semantics.
 **/

template <typename Number, GraphRep Graph = GraphRep::adapt_ss>
class DefaultParams {
public:
  using Wt = int64_t;

  using graph_t = typename std::conditional<
      (Graph == ss), SparseWtGraph<Wt>,
      typename std::conditional<
          (Graph == adapt_ss), AdaptGraph<Wt>,
          typename std::conditional<(Graph == pt), PtGraph<Wt>,
                                    HtGraph<Wt>>::type>::type>::type;
};

// We don't use GraphRep::adapt_ss because having problems
// realloc'ed Number objects.
template <typename Number, GraphRep Graph = GraphRep::ss>
class BigNumDefaultParams {
public:
  using Wt = Number;

  using graph_t = typename std::conditional<
      (Graph == ss), SparseWtGraph<Wt>,
      typename std::conditional<
          (Graph == adapt_ss), AdaptGraph<Wt>,
          typename std::conditional<(Graph == pt), PtGraph<Wt>,
                                    HtGraph<Wt>>::type>::type>::type;
};

template <typename Number, GraphRep Graph = GraphRep::adapt_ss>
class SafeInt64DefaultParams {
public:
  using Wt = safe_i64;

  using graph_t = typename std::conditional<
      (Graph == ss), SparseWtGraph<Wt>,
      typename std::conditional<
          (Graph == adapt_ss), AdaptGraph<Wt>,
          typename std::conditional<(Graph == pt), PtGraph<Wt>,
                                    HtGraph<Wt>>::type>::type>::type;
};

/**
 * Helper to translate from Number to DBM Wt (graph weights).  Number
 * is the template parameter of the DBM-based abstract domain to
 * represent a number. Number might not fit into Wt type.
 **/
template <typename Number, typename Wt> struct NtoW {
  static Wt convert(const Number &n, bool &overflow) {
    overflow = false;
    return static_cast<Wt>(n);
  }
};

template <> struct NtoW<ikos::z_number, int64_t> {
  static int64_t convert(const ikos::z_number &n, bool &overflow) {
    overflow = false;
    if (!n.fits_int64()) {
      overflow = true;
      return 0;
    }
    return static_cast<int64_t>(n);
  }
};

template <> struct NtoW<ikos::z_number, safe_i64> {
  static safe_i64 convert(const ikos::z_number &n, bool &overflow) {
    overflow = false;
    if (!n.fits_int64()) {
      overflow = true;
      return 0;
    }
    return safe_i64(n);
  }
};

} // namespace DBM_impl
} // namespace domains
} // namespace crab
