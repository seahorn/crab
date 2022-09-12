#pragma once

#include <crab/support/os.hpp>
#include <memory>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

/**
 * An adaptive sparse-map.
 *
 * Starts off as an unsorted vector, switching to a
 * sparse-set when |S| >= sparse_threshold
 * WARNING: Assumes Val is a basic type (so doesn't need a ctor/dtor call)
 **/

namespace crab {
namespace adapt_sgraph_impl {
inline void *malloc_fail(size_t size) {
  void *m = malloc(size);
  if (!m) {
    CRAB_ERROR("Allocation failure");
  }
  return m;
}

inline void *realloc_fail(void *ptr, size_t size) {
  void *m = realloc(ptr, size);
  if (!m) {
    CRAB_ERROR("Allocation failure");
  }
  return m;
}
} // namespace adapt_sgraph_impl

template <class Val> class AdaptSMap {
  enum { sparse_threshold = 8 };

public:
  using key_t = uint16_t;
  using val_t = Val;
  class elt_t {
  public:
    elt_t(key_t _k, const val_t &_v) : key(_k), val(_v) {}
    elt_t(const elt_t &o) = default;
    elt_t &operator=(const elt_t &o) = default;

    key_t key;
    val_t val;
  };

private:
  size_t sz;
  size_t dense_maxsz;
  size_t sparse_ub;
  elt_t *dense;
  key_t *sparse;

public:
  AdaptSMap(void)
      : sz(0), dense_maxsz(sparse_threshold), sparse_ub(10),
        dense((elt_t *)adapt_sgraph_impl::malloc_fail(sizeof(elt_t) *
                                                      sparse_threshold)),
        sparse(nullptr) {}

  AdaptSMap(AdaptSMap &&o)
      : sz(o.sz), dense_maxsz(o.dense_maxsz), sparse_ub(o.sparse_ub),
        dense(o.dense), sparse(o.sparse) {
    o.dense = nullptr;
    o.sparse = nullptr;
    o.sz = 0;
    o.dense_maxsz = 0;
    o.sparse_ub = 0;
  }

  AdaptSMap(const AdaptSMap &o)
      : sz(o.sz), dense_maxsz(o.dense_maxsz), sparse_ub(o.sparse_ub),
        dense((elt_t *)adapt_sgraph_impl::malloc_fail(sizeof(elt_t) *
                                                      dense_maxsz)),
        sparse(nullptr) {

    memcpy(static_cast<void *>(dense), o.dense, sizeof(elt_t) * sz);
    if (o.sparse) {
      sparse =
          (key_t *)adapt_sgraph_impl::malloc_fail(sizeof(key_t) * sparse_ub);
      for (key_t idx = 0; idx < sz; idx++)
        sparse[dense[idx].key] = idx;
    }
  }

  AdaptSMap &operator=(const AdaptSMap &o) {
    if (this != &o) {
      if (dense_maxsz < o.dense_maxsz) {
        dense_maxsz = o.dense_maxsz;
        dense = (elt_t *)adapt_sgraph_impl::realloc_fail(
            static_cast<void *>(dense), sizeof(elt_t) * dense_maxsz);
      }
      sz = o.sz;
      memcpy(static_cast<void *>(dense), o.dense, sizeof(elt_t) * sz);

      if (o.sparse) {
        if (!sparse || sparse_ub < o.sparse_ub) {
          sparse_ub = o.sparse_ub;
          sparse = (key_t *)adapt_sgraph_impl::realloc_fail(
              static_cast<void *>(sparse), sizeof(key_t) * sparse_ub);
        }
      }

      if (sparse) {
        for (key_t idx = 0; idx < sz; idx++)
          sparse[dense[idx].key] = idx;
      }
    }

    return *this;
  }

  AdaptSMap &operator=(AdaptSMap &&o) {
    if (this != &o) {
      if (dense)
        free(dense);
      if (sparse)
        free(sparse);

      dense = o.dense;
      o.dense = nullptr;
      sparse = o.sparse;
      o.sparse = nullptr;

      sz = o.sz;
      o.sz = 0;
      dense_maxsz = o.dense_maxsz;
      o.dense_maxsz = 0;
      sparse_ub = o.sparse_ub;
      o.sparse_ub = 0;
    }
    return *this;
  }

  ~AdaptSMap(void) {
    if (dense)
      free(dense);
    if (sparse)
      free(sparse);
  }

  size_t size(void) const { return sz; }

  class key_iter_t {
  public:
    key_iter_t(void) : e(nullptr) {}
    key_iter_t(elt_t *_e) : e(_e) {}

    // XXX: to make sure that we always return the same address
    // for the "empty" iterator, otherwise we can trigger
    // undefined behavior.
    static key_iter_t empty_iterator() {
      static std::unique_ptr<key_iter_t> it = nullptr;
      if (!it)
        it = std::unique_ptr<key_iter_t>(new key_iter_t());
      return *it;
    }

    key_t operator*(void) const { return (*e).key; }
    bool operator!=(const key_iter_t &o) const { return e < o.e; }
    key_iter_t &operator++(void) {
      ++e;
      return *this;
    }

    elt_t *e;
  };
  using elt_iter_t = elt_t *;

  class key_range_t {
  public:
    using iterator = key_iter_t;

    key_range_t(elt_t *_e, size_t _sz) : e(_e), sz(_sz) {}
    size_t size(void) const { return sz; }

    key_iter_t begin(void) const { return key_iter_t(e); }
    key_iter_t end(void) const { return key_iter_t(e + sz); }

    elt_t *e;
    size_t sz;
  };

  class elt_range_t {
  public:
    using iterator = elt_iter_t;

    elt_range_t(elt_t *_e, size_t _sz) : e(_e), sz(_sz) {}
    elt_range_t(const elt_range_t &o) : e(o.e), sz(o.sz) {}
    size_t size(void) const { return sz; }
    elt_iter_t begin(void) const { return elt_iter_t(e); }
    elt_iter_t end(void) const { return elt_iter_t(e + sz); }

    elt_t *e;
    size_t sz;
  };

  elt_range_t elts(void) const { return elt_range_t(dense, sz); }
  key_range_t keys(void) const { return key_range_t(dense, sz); }

  bool elem(key_t k) const {
    if (sparse) {
      if (k >= sparse_ub) {
        return false;
      }
      int idx = sparse[k];
      // WARNING: Valgrind will complain about uninitialized memory
      // being accessed. AddressSanitizer does not care about
      // uninitialized memory but MemorySanitizer might also
      // complain. This code is not portable. We are triggering
      // undefined behavior (dense[idx].key might be uninitialized)
      // and relying on the fact that compilers such as gcc and Clang
      // will not optimized code.
      //
      // A possible solution is to allocate sparse with calloc and
      // hoping that the OS does lazily the initialization.
      return (idx < sz) && dense[idx].key == k;
    } else {
      for (key_t ke : keys()) {
        if (ke == k)
          return true;
      }
      return false;
    }
  }

  bool lookup(key_t k, val_t *v_out) const {
    if (sparse) {
      if (k >= sparse_ub) {
        return false;
      }
      int idx = sparse[k];
      // SEE ABOVE WARNING
      if (idx < sz && dense[idx].key == k) {
        (*v_out) = dense[idx].val;
        return true;
      }
      return false;
    } else {
      for (elt_t elt : elts()) {
        if (elt.key == k) {
          (*v_out) = elt.val;
          return true;
        }
      }
      return false;
    }
  }

  // precondition: k \in S
  void remove(key_t k) {
    --sz;
    elt_t repl = dense[sz];
    if (sparse) {
      if (k >= sparse_ub) {
        return;
      }
      int idx = sparse[k];
      dense[idx] = repl;
      sparse[repl.key] = idx;
    } else {
      elt_t *e = dense;
      while (e->key != k)
        ++e;
      *e = repl;
    }
  }

  // precondition: k \notin S
  void add(key_t k, const val_t &v) {
    if (dense_maxsz <= sz)
      growDense(sz + 1);

    dense[sz] = elt_t(k, v);
    if (sparse) {
      if (sparse_ub <= k)
        growSparse(k + 1);
      sparse[k] = sz;
    }
    sz++;
  }

  void growDense(size_t new_max) {
    assert(dense_maxsz < new_max);

    while (dense_maxsz < new_max)
      dense_maxsz *= 2;
    elt_t *new_dense = (elt_t *)adapt_sgraph_impl::realloc_fail(
        static_cast<void *>(dense), sizeof(elt_t) * dense_maxsz);
    dense = new_dense;

    if (!sparse) {
      // After resizing the first time, we switch to an sset.
      key_t key_max = 0;
      for (key_t k : keys())
        key_max = std::max(key_max, k);

      sparse_ub = key_max + 1;
      sparse =
          (key_t *)adapt_sgraph_impl::malloc_fail(sizeof(key_t) * sparse_ub);
      key_t idx = 0;
      for (key_t k : keys())
        sparse[k] = idx++;
    }
  }

  void growSparse(size_t new_ub) {
    while (sparse_ub < new_ub)
      sparse_ub *= 2;
    key_t *new_sparse =
        (key_t *)adapt_sgraph_impl::malloc_fail(sizeof(key_t) * (sparse_ub));
    free(sparse);
    sparse = new_sparse;

    key_t idx = 0;
    for (key_t k : keys())
      sparse[k] = idx++;
  }

  void clear(void) { sz = 0; }
};
} // namespace crab
#pragma GCC diagnostic pop
