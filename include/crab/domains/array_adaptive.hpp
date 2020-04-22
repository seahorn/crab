/*******************************************************************************
 * Adaptive array domain
 *
 * Initially, an array is modeled by mapping sequences of consecutive
 * bytes (segments) to cells. A cell is pair <offset, size> where:
 *
 * - offset is an unsigned number
 * - size is an unsigned number
 *
 * A cell, when associated to an array A, is mapped to a scalar
 * variable in the base domain representing the byte contents of the
 * array segment A[offset,...,offset+size-1]
 *
 * The domain is general enough to represent any possible sequence of
 * consecutive bytes including sequences of bytes starting at the same
 * offsets but different sizes, overlapping sequences starting at
 * different offsets, etc.
 *
 * However, the domain only keeps track precisely of all bytes
 * contents as long as array writes use constant indexes. If an array
 * write with a non-constant index occurs then the array is _smashed_
 * into a single summarized variable. After that, all array writes are
 * modeled as weak updates but we still ensure that the number of
 * accessed bytes and the offset are consistent with previous array
 * accesses.
 ******************************************************************************/

#pragma once

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

#include <crab/domains/array_smashing.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/patricia_trees.hpp>

#include <algorithm>
#include <functional>
#include <set>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

// forward declaration
template <typename Domain, class Params> class array_adaptive_domain;

namespace array_adaptive_impl {

class DefaultParams {
public:
  /* options for array smashing */  
  enum { is_smashable = 1 };
  enum { smash_at_nonzero_offset = 1};
  enum { max_smashable_cells = 512};
  /* options for array expansion */
  // only used for now by array_store_range operation
  enum { max_array_size = 512 };
};

class NoSmashableParams {
public:
  /* options for array smashing */
  enum { is_smashable = 0 };
  enum { smash_at_nonzero_offset = 0};
  enum { max_smashable_cells = 512};
  /* options for array expansion */
  // only used for now by array_store_range operation  
  enum { max_array_size = 512 };
};

// Trivial constant propagation lattice
class cp_domain_t {
  typedef bound<ikos::z_number> bound_t;

  bool m_is_bottom;
  bound_t m_val;

  void set_to_top() {
    m_is_bottom = false;
    m_val = bound_t::plus_infinity();
  }

  void set_to_bot() {
    m_is_bottom = true;
    m_val = bound_t::plus_infinity();
  }

  cp_domain_t(bool is_bottom)
      : m_is_bottom(is_bottom), m_val(bound_t::plus_infinity()) {}

public:
  cp_domain_t() : m_is_bottom(false), m_val(bound_t::plus_infinity()) {}

  cp_domain_t(int64_t sz) : m_is_bottom(false), m_val(sz) {}

  cp_domain_t(const interval<ikos::z_number> &sz)
      : m_is_bottom(false), m_val(bound_t::plus_infinity()) {

    if (auto v_opt = sz.singleton()) {
      m_val = bound_t(*v_opt);
    }
  }

  bool is_top() const {
    return (!is_bottom() && (m_val == bound_t::plus_infinity()));
  }

  bool is_bottom() const { return m_is_bottom; }

  bool is_zero() const { return !m_is_bottom && m_val == bound_t(0); }

  bool is_negative() const { return !m_is_bottom && m_val < bound_t(0); }

  bound_t val() const {
    if (is_bottom()) {
      CRAB_ERROR("cp_domain_t::val cannot be called on bottom");
    }
    return m_val;
  }

  boost::optional<uint64_t> get_uint64_val() const {
    if (is_bottom()) {
      CRAB_ERROR("cp_domain_t::get_uint64_val cannot be called on bottom");
    }

    if (m_val.number()) {
      ikos::z_number n = *(m_val.number());
      if (n.fits_int64()) {
        if (static_cast<int64_t>(n) > 0) {
          return (uint64_t) static_cast<int64_t>(n);
        }
      }
    }
    return boost::optional<uint64_t>();
  }

  static cp_domain_t bottom() { return cp_domain_t(true); }

  static cp_domain_t top() { return cp_domain_t(false); }

  bool operator==(const cp_domain_t &o) const {
    return *this <= o && o <= *this;
  }

  bool operator<=(const cp_domain_t &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    } else {
      return (m_val == o.m_val);
    }
  }

  void operator|=(const cp_domain_t &o) {
    if (is_bottom() || o.is_top()) {
      *this = o;
    } else if (o.is_bottom() || is_top()) {
      // do nothing
    } else if (m_val != o.m_val) {
      set_to_top();
    } else {
      // do nothing
    }
  }

  cp_domain_t operator|(const cp_domain_t &o) const {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() || o.is_top()) {
      return cp_domain_t();
    } else {
      if (m_val == o.m_val) {
        return *this;
      } else {
        return cp_domain_t();
      }
    }
  }

  cp_domain_t operator&(const cp_domain_t &o) const {
    if (is_bottom() || o.is_bottom()) {
      return bottom();
    } else if (is_top()) {
      return o;
    } else if (o.is_top()) {
      return *this;
    } else {
      if (m_val == o.m_val) {
        return *this;
      } else {
        return bottom();
      }
    }
  }

  cp_domain_t operator||(const cp_domain_t &o) const { return operator|(o); }

  cp_domain_t operator&&(const cp_domain_t &o) const { return operator&(o); }

  void write(crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else {
      o << m_val;
    }
  }
};

// forward declaration
class offset_map;

/*
 * Wrapper for using ikos::index_t as patricia_tree keys
 */
class offset_t {
  ikos::index_t m_val;

public:
  explicit offset_t(ikos::index_t v) : m_val(v) {}

  ikos::index_t index() const { return m_val; }

  size_t hash() const {
    // casting to size_t may overflow but it shouldn't affect correctness
    return std::hash<size_t>{}(static_cast<size_t>(m_val));
  }

  bool operator<(const offset_t &o) const { return m_val < o.m_val; }

  bool operator==(const offset_t &o) const { return m_val == o.m_val; }

  bool operator!=(const offset_t &o) const { return !(*this == o); }

  offset_t operator%(const offset_t &o) const {
    return offset_t(m_val % o.m_val);
  }

  offset_t operator-(const offset_t &o) const {
    return offset_t(m_val - o.m_val);
  }

  bool is_negative() const { return m_val < 0; }

  bool is_zero() const { return m_val == 0; }

  void write(crab::crab_os &o) const { o << m_val; }

  friend crab::crab_os &operator<<(crab::crab_os &o, const offset_t &v) {
    v.write(o);
    return o;
  }
};

/*
 *  A synthetic cell is used to give a symbolic name to the byte
 *  contents of some array segment represented by
 *
 *     [m_offset, m_offset+1,...,m_offset+m_size-1]
 *
 */
class cell_t {
private:
  friend class offset_map_t;
  typedef ikos::interval<ikos::z_number> interval_t;
  offset_t m_offset;
  uint64_t m_size;
  //// Boolean flag to indicate if the cell has been removed.  When
  //// smashing is enabled we need to know which cells have been
  //// created and join them regardless whether they have been removed
  //// or not. If this flag is enabled then overlap and
  //// symbolic_overlap pretend the cell does not exist and return
  //// always false but the cell is not destroyed.
  bool m_removed;

  // Only offset_map class can create cells
  cell_t() : m_offset(0), m_size(0), m_removed(false) {}

  cell_t(offset_t offset, uint64_t size)
      : m_offset(offset), m_size(size), m_removed(false) {}

  static interval_t to_interval(const offset_t o, uint64_t size) {
    interval_t i(o.index(), o.index() + size - 1);
    return i;
  }

  interval_t to_interval() const {
    return to_interval(get_offset(), get_size());
  }

public:
  bool is_null() const { return (m_offset.index() == 0 && m_size == 0); }

  offset_t get_offset() const { return m_offset; }

  size_t get_size() const { return m_size; }

  cell_t clone(void) const {
    cell_t new_cell;
    new_cell.m_offset = m_offset;
    new_cell.m_size = m_size;
    new_cell.m_removed = m_removed;
    return new_cell;
  }

  void mark_as_removed(bool v) { m_removed = v; }

  bool is_removed(void) const { return m_removed; }

  size_t hash() const {
    size_t h1 = m_offset.hash();
    size_t h2 = std::hash<size_t>{}(m_size);
    return h1 ^ (h2 << 1);
  }

  // inclusion test
  bool operator<=(const cell_t &o) const {
    interval_t x = to_interval();
    interval_t y = o.to_interval();
    return x <= y;
  }

  bool operator==(const cell_t &o) const {
    return (get_offset() == o.get_offset() && get_size() == o.get_size());
  }

  bool operator<(const cell_t &o) const {
    if (get_offset() < o.get_offset()) {
      return true;
    } else if (get_offset() == o.get_offset()) {
      return get_size() < o.get_size();
    } else {
      return false;
    }
  }

  // Return true if [o, o+size) definitely overlaps with the cell,
  // where o is a constant expression.
  bool overlap(const offset_t &o, uint64_t size) const {
    if (m_removed)
      return false;

    interval_t x = to_interval();
    interval_t y = to_interval(o, size);
    bool res = (!(x & y).is_bottom());
    CRAB_LOG("array-adaptive-overlap", crab::outs() << "**Checking if " << x
                                                    << " overlaps with " << y
                                                    << "=" << res << "\n";);
    return res;
  }

  // Return true if [symb_lb, symb_ub] may overlap with the cell,
  // where symb_lb and symb_ub are not constant expressions.
  template <typename Dom>
  bool symbolic_overlap(const typename Dom::linear_expression_t &symb_lb,
                        const typename Dom::linear_expression_t &symb_ub,
                        const Dom &dom) const {
    if (m_removed)
      return false;

    typedef typename Dom::linear_expression_t linear_expression_t;

    interval_t x = to_interval();
    assert(x.lb().is_finite());
    assert(x.ub().is_finite());
    linear_expression_t lb(*(x.lb().number()));
    linear_expression_t ub(*(x.ub().number()));

    CRAB_LOG("array-adaptive-overlap", Dom tmp(dom);
             linear_expression_t tmp_symb_lb(symb_lb);
             linear_expression_t tmp_symb_ub(symb_ub);
             crab::outs() << "**Checking if " << *this
                          << " overlaps with symbolic "
                          << "[" << tmp_symb_lb << "," << tmp_symb_ub << "]"
                          << " with abstract state=" << tmp << "\n";);

    Dom tmp1(dom);
    tmp1 += (lb >= symb_lb);
    tmp1 += (lb <= symb_ub);
    if (!tmp1.is_bottom()) {
      CRAB_LOG("array-adaptive-overlap", crab::outs() << "\tyes.\n";);
      return true;
    }

    Dom tmp2(dom);
    tmp2 += (ub >= symb_lb);
    tmp2 += (ub <= symb_ub);
    if (!tmp2.is_bottom()) {
      CRAB_LOG("array-adaptive-overlap", crab::outs() << "\tyes.\n";);
      return true;
    }

    CRAB_LOG("array-adaptive-overlap", crab::outs() << "\tno.\n";);
    return false;
  }

  void write(crab::crab_os &o) const {
    if (is_null()) {
      o << "NULL";
    } else {
      o << to_interval() << " ";
      if (m_removed) {
        o << "R";
      }
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, const cell_t &c) {
    c.write(o);
    return o;
  }
};

namespace cell_set_impl {
template <typename Set> inline Set set_intersection(Set &s1, Set &s2) {
  Set s3;
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::inserter(s3, s3.end()));
  return s3;
}

template <typename Set> inline Set set_union(Set &s1, Set &s2) {
  Set s3;
  std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 std::inserter(s3, s3.end()));
  return s3;
}

template <typename Set> inline bool set_inclusion(Set &s1, Set &s2) {
  Set s3;
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(s3, s3.end()));
  return s3.empty();
}

template <typename Set> inline Set set_difference(Set &s1, Set &s2) {
  Set s3;
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(s3, s3.end()));
  return s3;
}

template <typename Set>
inline Set set_difference(Set &s1, const typename Set::key_type &e) {
  Set s3;
  Set s2;
  s2.insert(e);
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(s3, s3.end()));
  return s3;
}
} // end namespace cell_set_impl

/*
 * A Patricia tree that maps numerical offsets to synthetic cells.
 */
class offset_map_t {
private:
  template <typename Dom, class Params>
  friend class crab::domains::array_adaptive_domain;

  typedef std::set<cell_t> cell_set_t;
  typedef crab::variable_type type_t;

  /*
    The keys in the patricia tree are processing in big-endian
    order. This means that the keys are sorted. Sortedeness is
    very important to perform efficiently operations such as
    checking for overlap cells. Since keys are treated as bit
    patterns, negative offsets can be used but they are treated
    as large unsigned numbers.
  */
  typedef patricia_tree<offset_t, cell_set_t> patricia_tree_t;
  typedef typename patricia_tree_t::binary_op_t binary_op_t;
  typedef typename patricia_tree_t::partial_order_t partial_order_t;

  patricia_tree_t m_map;

  // for algorithm::lower_bound and algorithm::upper_bound
  struct compare_binding_t {
    bool operator()(const typename patricia_tree_t::binding_t &kv,
                    const offset_t &o) const {
      return kv.first < o;
    }
    bool operator()(const offset_t &o,
                    const typename patricia_tree_t::binding_t &kv) const {
      return o < kv.first;
    }
    bool operator()(const typename patricia_tree_t::binding_t &kv1,
                    const typename patricia_tree_t::binding_t &kv2) const {
      return kv1.first < kv2.first;
    }
  };

  patricia_tree_t apply_operation(binary_op_t &o, patricia_tree_t t1,
                                  patricia_tree_t t2) {
    bool res = t1.merge_with(t2, o);
    if (res) {
      CRAB_ERROR("array_adaptive::offset_map should not return bottom");
    }
    return t1;
  }

  class join_op : public binary_op_t {
    // apply is called when two bindings (one each from a
    // different map) have the same key(i.e., offset).
    std::pair<bool, boost::optional<cell_set_t>> apply(cell_set_t x,
                                                       cell_set_t y) {
      return {false, cell_set_impl::set_union(x, y)};
    }
    // if one map does not have a key in the other map we add it.
    bool default_is_absorbing() { return false; }
  };

  class meet_op : public binary_op_t {
    std::pair<bool, boost::optional<cell_set_t>> apply(cell_set_t x,
                                                       cell_set_t y) {
      return {false, cell_set_impl::set_union(x, y)};
    }
    // if one map does not have a key in the other map we ignore
    // it.
    bool default_is_absorbing() { return true; }
  };

  class domain_po : public partial_order_t {
    bool leq(cell_set_t x, cell_set_t y) {
      return cell_set_impl::set_inclusion(x, y);
    }
    // default value is bottom (i.e., empty map)
    bool default_is_top() { return false; }
  }; // class domain_po

  // Delete completely the cell
  void erase_cell(const cell_t &c) {
    if (boost::optional<cell_set_t> cells = m_map.lookup(c.get_offset())) {
      if ((*cells).erase(c) > 0) {
        m_map.remove(c.get_offset());
        if (!(*cells).empty()) {
          // a bit of a waste ...
          m_map.insert(c.get_offset(), *cells);
        }
      }
    }
  }

  // Pretend the cell is removed by marking it as "removed"
  void remove_cell(const cell_t &c) {
    if (boost::optional<cell_set_t> cells = m_map.lookup(c.get_offset())) {
      if ((*cells).count(c) > 0) {
        m_map.remove(c.get_offset());
        cell_set_t ncells = cell_set_impl::set_difference(*cells, c);
        cell_t nc = c.clone();
        nc.mark_as_removed(true);
        ncells.insert(nc);
        m_map.insert(c.get_offset(), ncells);
      }
    }
  }

  void insert_cell(const cell_t &c) {
    if (boost::optional<cell_set_t> cells = m_map.lookup(c.get_offset())) {
      if ((*cells).insert(c).second) {
        m_map.remove(c.get_offset());
        m_map.insert(c.get_offset(), *cells);
      }
    } else {
      cell_set_t new_cells;
      new_cells.insert(c);
      m_map.insert(c.get_offset(), new_cells);
    }
  }

  cell_t get_cell(offset_t o, uint64_t size) const {
    if (boost::optional<cell_set_t> cells = m_map.lookup(o)) {
      cell_t tmp(o, size);
      auto it = (*cells).find(tmp);
      if (it != (*cells).end()) {
        return *it;
      }
    }
    // not found
    return cell_t();
  }

  // create a fresh _unamed_ cell
  cell_t mk_cell(offset_t o, uint64_t size) {
    cell_t c = get_cell(o, size);
    if (c.is_null()) {
      cell_t nc(o, size);
      insert_cell(nc);
      c = nc;
      CRAB_LOG("array-adaptive", crab::outs()
                                     << "**Created cell " << c << "\n";);
    }
    if (c.is_null()) {
      CRAB_ERROR("cannot create a null cell from offset ", o, " and size ",
                 size);
    }
    return c;
  }

  offset_map_t(patricia_tree_t &&m) : m_map(std::move(m)) {}

public:
  offset_map_t() {}

  offset_map_t(const offset_map_t &o) : m_map(o.m_map) {}

  offset_map_t(const offset_map_t &&o) : m_map(std::move(o.m_map)) {}

  offset_map_t &operator=(const offset_map_t &o) {
    if (this != &o) {
      m_map = o.m_map;
    }
    return *this;
  }

  offset_map_t &operator=(const offset_map_t &&o) {
    if (this != &o) {
      m_map = std::move(o.m_map);
    }
    return *this;
  }

  bool empty() const { return m_map.empty(); }

  std::size_t size() const { return m_map.size(); }

  // leq operator
  bool operator<=(const offset_map_t &o) const {
    domain_po po;
    return m_map.leq(o.m_map, po);
  }

  // set union: if two cells with same offset do not agree on
  // size then they are ignored.
  offset_map_t operator|(const offset_map_t &o) {
    join_op op;
    offset_map_t res = offset_map_t(apply_operation(op, m_map, o.m_map));
    CRAB_LOG("array-adaptive-offset-map",
             crab::outs() << "offset_map join\n"
                          << *this << "\nand\n"
                          << o << "\nResult=" << res << "\n";);
    return res;
  }

  // set intersection: if two cells with same offset do not agree
  // on size then they are ignored.
  offset_map_t operator&(const offset_map_t &o) {
    meet_op op;
    offset_map_t res = offset_map_t(apply_operation(op, m_map, o.m_map));
    CRAB_LOG("array-adaptive-offset-map",
             crab::outs() << "offset_map meet\n"
                          << *this << "\nand\n"
                          << o << "\nResult=" << res << "\n";);
    return res;
  }

  // Completely delete the cell from the offset map
  void erase(const cell_t &c) { erase_cell(c); }

  void erase(const std::vector<cell_t> &cells) {
    for (unsigned i = 0, e = cells.size(); i < e; ++i) {
      erase(cells[i]);
    }
  }

  // Pretend the cell is removed so no other cells overlap with it but
  // the cell is not actually deleted from the offset map. We need to
  // know all created cells when smashing occurs.
  void remove(const cell_t &c) { remove_cell(c); }

  void remove(const std::vector<cell_t> &cells) {
    for (unsigned i = 0, e = cells.size(); i < e; ++i) {
      remove(cells[i]);
    }
  }

  // cells are sorted by offset
  std::vector<cell_t> get_all_cells() const {
    std::vector<cell_t> res;
    for (auto it = m_map.begin(), et = m_map.end(); it != et; ++it) {
      auto const &o_cells = it->second;
      for (auto &c : o_cells) {
        res.push_back(c);
      }
    }
    return res;
  }

  // Return in out all cells that might overlap with (o, size).
  //
  // It is not marked as const because we insert temporary a cell.
  // However, upon completion this method leaves unmodified the object.
  void get_overlap_cells(offset_t o, uint64_t size, std::vector<cell_t> &out) {
    compare_binding_t comp;

    bool added = false;
    cell_t c = get_cell(o, size);
    if (c.is_null()) {
      // we need to add a temporary cell for (o, size)
      c = cell_t(o, size);
      insert_cell(c);
      added = true;
    }

    auto lb_it = std::lower_bound(m_map.begin(), m_map.end(), o, comp);
    if (lb_it != m_map.end()) {
      // Store m_map[begin,...,lb_it] into a vector so that we can
      // go backwards from lb_it.
      //
      // TODO: give support for reverse iterator in patricia_tree.
      std::vector<cell_set_t> upto_lb;
      upto_lb.reserve(std::distance(m_map.begin(), lb_it));
      for (auto it = m_map.begin(), et = lb_it; it != et; ++it) {
        upto_lb.push_back(it->second);
      }
      upto_lb.push_back(lb_it->second);

      for (int i = upto_lb.size() - 1; i >= 0; --i) {
        ///////
        // All the cells in upto_lb[i] have the same offset. They
        // just differ in the size.
        //
        // If none of the cells in upto_lb[i] overlap with (o, size)
        // we can stop.
        ////////
        bool continue_outer_loop = false;
        for (const cell_t &x : upto_lb[i]) {
          if (x.overlap(o, size)) {
            if (!(x == c)) {
              // FIXME: we might have some duplicates. this is a very drastic
              // solution.
              if (std::find(out.begin(), out.end(), x) == out.end()) {
                out.push_back(x);
              }
            }
            continue_outer_loop = true;
          }
        }
        if (!continue_outer_loop) {
          break;
        }
      }
    }

    // search for overlapping cells > o
    auto ub_it = std::upper_bound(m_map.begin(), m_map.end(), o, comp);
    for (; ub_it != m_map.end(); ++ub_it) {
      bool continue_outer_loop = false;
      for (const cell_t &x : ub_it->second) {
        if (x.overlap(o, size)) {
          // FIXME: we might have some duplicates. this is a very drastic
          // solution.
          if (std::find(out.begin(), out.end(), x) == out.end()) {
            out.push_back(x);
          }
          continue_outer_loop = true;
        }
      }
      if (!continue_outer_loop) {
        break;
      }
    }

    // do not forget the rest of overlapping cells == o
    for (auto it = ++lb_it, et = ub_it; it != et; ++it) {
      bool continue_outer_loop = false;
      for (const cell_t &x : it->second) {
        if (x == c) { // we dont put it in out
          continue;
        }
        if (x.overlap(o, size)) {
          if (!(x == c)) {
            if (std::find(out.begin(), out.end(), x) == out.end()) {
              out.push_back(x);
            }
          }
          continue_outer_loop = true;
        }
      }
      if (!continue_outer_loop) {
        break;
      }
    }

    if (added) {
      // remove the temporary cell for (o, size)
      assert(!c.is_null());
      erase_cell(c);
    }

    CRAB_LOG(
        "array-adaptive-overlap", crab::outs()
                                      << "**Overlap set between \n"
                                      << *this << "\nand "
                                      << "(" << o << "," << size << ")={";
        for (unsigned i = 0, e = out.size(); i < e;) {
          crab::outs() << out[i];
          ++i;
          if (i < e) {
            crab::outs() << ",";
          }
        } crab::outs()
        << "}\n";);
  }

  template <typename Dom>
  void get_overlap_cells_symbolic_offset(
      const Dom &dom, const typename Dom::linear_expression_t &symb_lb,
      const typename Dom::linear_expression_t &symb_ub,
      std::vector<cell_t> &out) const {

    for (auto it = m_map.begin(), et = m_map.end(); it != et; ++it) {
      const cell_set_t &o_cells = it->second;
      // All cells in o_cells have the same offset. They only differ
      // in the size. If the largest cell overlaps with [offset,
      // offset + size) then the rest of cells are considered to
      // overlap. This is an over-approximation because [offset,
      // offset+size) can overlap with the largest cell but it
      // doesn't necessarily overlap with smaller cells. For
      // efficiency, we assume it overlaps with all.
      cell_t largest_cell;
      for (auto &c : o_cells) {
        if (largest_cell.is_null()) {
          largest_cell = c;
        } else {
          assert(c.get_offset() == largest_cell.get_offset());
          if (largest_cell < c) {
            largest_cell = c;
          }
        }
      }
      if (!largest_cell.is_null()) {
        if (largest_cell.symbolic_overlap(symb_lb, symb_ub, dom)) {
          for (auto &c : o_cells) {
            out.push_back(c);
          }
        }
      }
    }
  }

  void clear(void) { m_map.clear(); }

  void write(crab::crab_os &o) const {
    if (m_map.empty()) {
      o << "empty";
    } else {
      for (auto it = m_map.begin(), et = m_map.end(); it != et; ++it) {
        const cell_set_t &cells = it->second;
        o << "{";
        for (auto cit = cells.begin(), cet = cells.end(); cit != cet;) {
          if ((*cit).is_removed()) {
            ++cit;
            continue;
          }
          o << *cit;
          ++cit;
          if (cit != cet) {
            o << ",";
          }
        }
        o << "}";
      }
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, const offset_map_t &m) {
    m.write(o);
    return o;
  }

  /* Operations needed if used as value in a patricia tree */
  bool operator==(const offset_map_t &o) const {
    return *this <= o && o <= *this;
  }
  bool is_top() const { return empty(); }
  bool is_bottom() const { return false; }
  /*
     a patricia tree only calls bottom if operator[] is called over
     a bottom state. Thus, we will make sure that we don't call
     operator[] in that case.
  */
  static offset_map_t bottom() {
    CRAB_ERROR("offset_map::bottom() cannot be called");
  }
  /* Top is called when a key is not found in a patricia tree */
  static offset_map_t top() { return offset_map_t(); }
};
} // end namespace array_adaptive_impl

template <typename NumDomain, class Params = array_adaptive_impl::DefaultParams>
class array_adaptive_domain final
    : public abstract_domain<array_adaptive_domain<NumDomain, Params>> {

public:
  typedef typename NumDomain::number_t number_t;
  typedef typename NumDomain::varname_t varname_t;

private:
  typedef array_adaptive_domain<NumDomain, Params> array_adaptive_domain_t;
  typedef abstract_domain<array_adaptive_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  typedef typename NumDomain::variable_t variable_t;
  typedef typename NumDomain::variable_vector_t variable_vector_t;
  typedef crab::pointer_constraint<variable_t> ptr_cst_t;
  typedef interval<number_t> interval_t;
  typedef array_smashing<NumDomain> base_domain_t;

private:
  typedef crab::variable_type type_t;
  typedef array_adaptive_impl::offset_t offset_t;
  typedef array_adaptive_impl::offset_map_t offset_map_t;
  typedef array_adaptive_impl::cell_t cell_t;
  typedef array_adaptive_impl::cp_domain_t cp_domain_t;
  struct cell_varmap_hasher {
    std::size_t operator()(const std::pair<variable_t, cell_t> &p) const {
      size_t h1 = p.first.hash();
      size_t h2 = p.second.hash();
      return h1 ^ (h2 << 1);
    }
  };
  struct cell_varmap_equal {
    std::size_t operator()(const std::pair<variable_t, cell_t> &p1,
                           const std::pair<variable_t, cell_t> &p2) const {
      return ((p1.first == p2.first) ? p1.second == p2.second : false);
    }
  };
  typedef std::unordered_map<std::pair<variable_t, cell_t>, variable_t,
                             cell_varmap_hasher, cell_varmap_equal>
      cell_varmap_t;
  typedef std::unordered_map<variable_t, variable_t> smashed_varmap_t;

  class array_state {
    // whether the array has been smashed
    bool m_is_smashed;
    // element size of the array if smashed
    cp_domain_t m_element_sz;
    // precise array contents if no smashed
    offset_map_t m_offset_map;

    static bool consistent_offset(offset_t o, size_t elem_size) {
      if (o.is_negative()) {
        CRAB_LOG("array-adaptive-smash",
                 crab::outs() << "cannot smash because negative offset\n";);
        return false;
      }
      size_t n = static_cast<size_t>(o.index());
      if (n % elem_size != 0) {
        CRAB_LOG("array-adaptive-smash", crab::outs()
                                             << "cannot smash because " << n
                                             << "%" << elem_size << "!= 0\n");
        return false;
      }
      return true;
    }

    // Smash the array if it's safe to do so.
    void smash_array(const variable_t &a, const cp_domain_t &elem_sz,
                     cell_varmap_t &cvm, smashed_varmap_t &svm,
                     base_domain_t &base_dom) {

      // we only smash the array if elem_sz is a constant value and
      // all array elements are consistent wrt elem_sz. Leave the
      // array without smashing is always sound.

      std::vector<cell_t> cells = get_offset_map().get_all_cells();
      if (cells.empty()) {
        return;
      }

      if (cells.size() > Params::max_smashable_cells) {
	// Smashing is expensive because it will go over all cells
        // performing one join per weak update even if the smashed
        // array is already unconstrained. We don't smash if the
        // number of cells to be smashed is too large.
	CRAB_WARN("array adaptive did not smash array because its size is greater than ",
		  Params::max_smashable_cells);
	return;
      }
      
      if (!Params::smash_at_nonzero_offset && !cells[0].get_offset().is_zero()) {
        return;
      }

      if (boost::optional<uint64_t> sz_opt = elem_sz.get_uint64_val()) {
        bool can_be_smashed = true;
        variable_t smashed_a =
            array_adaptive_domain_t::get_smashed_variable(a, svm);
        for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
          const cell_t &c = cells[k];
          if (!consistent_offset(c.get_offset(), *sz_opt)) {
            can_be_smashed = false;
            break;
          }
          const bool is_strong_update = (k == 0);
          linear_expression_t idx(number_t(c.get_offset().index()));
          auto it = cvm.find({a, c});
          if (it != cvm.end()) {
            const variable_t &c_scalar_var = it->second;
            base_dom.array_store(smashed_a, c.get_size(), idx, c_scalar_var,
                                 is_strong_update);
          } else {
            can_be_smashed = false;
            break;
          }
        }

        if (can_be_smashed) {
          for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
            const cell_t &c = cells[k];
            auto it = cvm.find({a, c});
            if (it != cvm.end()) {
              const variable_t &c_scalar_var = it->second;
              // remove the synthethic cell from the base domain and
              // from the offset map
              base_dom -= c_scalar_var;
            }
            get_offset_map().erase(c);
            cvm.erase({a, c});
          }
          m_is_smashed = true;
          m_element_sz = elem_sz;
        } else {
          base_dom -= smashed_a;
        }
      }
    }

    void do_sanity_checks() const {
      if (m_element_sz.is_bottom()) {
        CRAB_ERROR("array_state::m_element_sz cannot be bottom");
      }
      if (m_offset_map.is_bottom()) {
        CRAB_ERROR("array_state::m_offset_map cannot be bottom");
      }
    }

  public:
    array_state() : m_is_smashed(false), m_element_sz((int64_t)0) {}

    array_state(bool &&is_smashed, cp_domain_t &&sz, offset_map_t &&om)
        : m_is_smashed(std::move(is_smashed)), m_element_sz(std::move(sz)),
          m_offset_map(std::move(om)) {
      do_sanity_checks();
    }

    array_state(const array_state &o)
        : m_is_smashed(o.m_is_smashed), m_element_sz(o.m_element_sz),
          m_offset_map(o.m_offset_map) {
      do_sanity_checks();
    }

    array_state(const array_state &&o)
        : m_is_smashed(std::move(o.m_is_smashed)),
          m_element_sz(std::move(o.m_element_sz)),
          m_offset_map(std::move(o.m_offset_map)) {
      do_sanity_checks();
    }

    array_state &operator=(const array_state &o) {
      if (this != &o) {
        m_is_smashed = o.m_is_smashed;
        m_element_sz = o.m_element_sz;
        m_offset_map = o.m_offset_map;
      }
      do_sanity_checks();
      return *this;
    }

    array_state &operator=(const array_state &&o) {
      if (this != &o) {
        m_is_smashed = std::move(o.m_is_smashed);
        m_element_sz = std::move(o.m_element_sz);
        m_offset_map = std::move(o.m_offset_map);
      }
      do_sanity_checks();
      return *this;
    }

    /*** begin mergeable_map API ***/
    array_state join(const variable_t &v, array_state &o,
                     cell_varmap_t &cvm_left, smashed_varmap_t &svm_left,
                     base_domain_t &dom_left, cell_varmap_t &cvm_right,
                     smashed_varmap_t &svm_right, base_domain_t &dom_right) {
      if (m_is_smashed && !o.m_is_smashed) {
        o.smash_array(v, get_element_sz(), cvm_right, svm_right, dom_right);
      } else if (!m_is_smashed && o.m_is_smashed) {
        smash_array(v, o.get_element_sz(), cvm_left, svm_left, dom_left);
      }
      return array_state(m_is_smashed | o.m_is_smashed,
                         m_element_sz | o.m_element_sz,
                         m_offset_map | o.m_offset_map);
    }

    array_state meet(const variable_t &v, array_state &o,
                     cell_varmap_t &cvm_left, smashed_varmap_t &svm_left,
                     base_domain_t &dom_left, cell_varmap_t &cvm_right,
                     smashed_varmap_t &svm_right, base_domain_t &dom_right) {
      if (m_is_smashed && !o.m_is_smashed) {
        o.smash_array(v, get_element_sz(), cvm_right, svm_right, dom_right);
      } else if (!m_is_smashed && o.m_is_smashed) {
        smash_array(v, o.get_element_sz(), cvm_left, svm_left, dom_left);
      }
      return array_state(m_is_smashed & o.m_is_smashed,
                         m_element_sz & o.m_element_sz,
                         m_offset_map & o.m_offset_map);
    }

    bool operator==(const array_state &o) {
      if (m_is_smashed != o.m_is_smashed) {
        return false;
      }
      if (m_is_smashed) {
        return (m_element_sz == o.m_element_sz);
      } else {
        return m_offset_map == o.m_offset_map;
      }
    }

    bool is_top() const { return false; }
    bool is_bottom() const { return false; }
    /*
     a patricia tree only calls bottom if operator[] is called over
     a bottom state. Thus, we will make sure that we don't call
     operator[] in that case.
    */
    static array_state bottom() {
      CRAB_ERROR("array_state::bottom() cannot be called");
    }
    /* Top is called when a key is not found in a patricia tree */
    static array_state top() {
      CRAB_ERROR("array_state::top() cannot be called");
      // return array_state();
    }

    /*** end mergeable_map API ***/

    bool is_smashed() const { return m_is_smashed; }

    void set_smashed(bool v) { m_is_smashed = v; }

    offset_map_t &get_offset_map() { return m_offset_map; }

    const offset_map_t &get_offset_map() const { return m_offset_map; }

    cp_domain_t &get_element_sz() { return m_element_sz; }

    const cp_domain_t &get_element_sz() const { return m_element_sz; }

    static bool can_be_smashed(const std::vector<cell_t> &cells,
                               uint64_t elem_sz,
                               bool allow_start_at_nonzero_offset) {
      if (cells.empty()) {
        return false;
      }

      offset_t start_offset = cells[0].get_offset();

      if (!allow_start_at_nonzero_offset) {
        if (!start_offset.is_zero()) {
          CRAB_LOG("array-adaptive-smash",
                   crab::outs() << "cannot smash because array does not start "
                                   "at offset 0\n";);
          return false;
        }
      }

      for (unsigned i = 0, e = cells.size(); i < e; ++i) {
        const cell_t &c = cells[i];
        if (c.get_size() != elem_sz) {
          CRAB_LOG("array-adaptive-smash",
                   crab::outs() << "cannot smash because array elements have "
                                   "different sizes\n");
          return false;
        }
        // note that we adjust the cell's offset
        if (!consistent_offset(c.get_offset() - start_offset, elem_sz)) {
          return false;
        }
      }
      return true;
    }

    bool can_be_smashed(uint64_t elem_sz) const {
      if (m_is_smashed) {
        // already smashed, bail out ...
        return false;
      }
      std::vector<cell_t> cells = m_offset_map.get_all_cells();
      return can_be_smashed(cells, elem_sz,
			    Params::smash_at_nonzero_offset);
    }

    void write(crab_os &o) const {
      if (m_is_smashed) {
        o << "smashed with element size=";
        m_element_sz.write(o);
      } else {
        m_offset_map.write(o);
      }
    }
  };

  class array_state_map_t {
  private:
    typedef patricia_tree<variable_t, array_state> patricia_tree_t;
    typedef typename patricia_tree_t::key_binary_op_t key_binary_op_t;

  public:
    typedef typename patricia_tree_t::iterator iterator;

  private:
    patricia_tree_t m_tree;

    class join_op : public key_binary_op_t {
      cell_varmap_t &m_cvm_left;
      smashed_varmap_t &m_svm_left;
      base_domain_t &m_dom_left;
      cell_varmap_t &m_cvm_right;
      smashed_varmap_t &m_svm_right;
      base_domain_t &m_dom_right;

    public:
      join_op(cell_varmap_t &cvm_left, smashed_varmap_t &svm_left,
              base_domain_t &dom_left, cell_varmap_t &cvm_right,
              smashed_varmap_t &svm_right, base_domain_t &dom_right)
          : m_cvm_left(cvm_left), m_svm_left(svm_left), m_dom_left(dom_left),
            m_cvm_right(cvm_right), m_svm_right(svm_right),
            m_dom_right(dom_right) {}
      std::pair<bool, boost::optional<array_state>>
      apply(const variable_t &k, array_state x, array_state y) {
        array_state z = x.join(k, y, m_cvm_left, m_svm_left, m_dom_left,
                               m_cvm_right, m_svm_right, m_dom_right);
        return {false, boost::optional<array_state>(z)};
      }

      bool default_is_absorbing() { return true; }
    }; // class join_op

    class meet_op : public key_binary_op_t {
      cell_varmap_t &m_cvm_left;
      smashed_varmap_t &m_svm_left;
      base_domain_t &m_dom_left;
      cell_varmap_t &m_cvm_right;
      smashed_varmap_t &m_svm_right;
      base_domain_t &m_dom_right;

    public:
      meet_op(cell_varmap_t &cvm_left, smashed_varmap_t &svm_left,
              base_domain_t &dom_left, cell_varmap_t &cvm_right,
              smashed_varmap_t &svm_right, base_domain_t &dom_right)
          : m_cvm_left(cvm_left), m_svm_left(svm_left), m_dom_left(dom_left),
            m_cvm_right(cvm_right), m_svm_right(svm_right),
            m_dom_right(dom_right) {}
      std::pair<bool, boost::optional<array_state>>
      apply(const variable_t &k, array_state x, array_state y) {
        array_state z = x.meet(k, y, m_cvm_left, m_svm_left, m_dom_left,
                               m_cvm_right, m_svm_right, m_dom_right);
        return {false, boost::optional<array_state>(z)};
      }
      bool default_is_absorbing() { return false; }
    }; // class meet_op

    patricia_tree_t apply_operation(key_binary_op_t &o, patricia_tree_t t1,
                                    patricia_tree_t t2) {
      bool res = t1.merge_with(t2, o);
      if (res) {
        CRAB_ERROR(
            "array_adaptive::array_state_map_t should not return bottom");
      }
      return t1;
    }

    array_state_map_t(patricia_tree_t &&t) : m_tree(std::move(t)) {}

  public:
    array_state_map_t() {}

    array_state_map_t(const array_state_map_t &o) : m_tree(o.m_tree) {}

    array_state_map_t(const array_state_map_t &&o)
        : m_tree(std::move(o.m_tree)) {}

    array_state_map_t &operator=(const array_state_map_t &o) {
      if (this != &o) {
        m_tree = o.m_tree;
      }
      return *this;
    }

    array_state_map_t &operator=(const array_state_map_t &&o) {
      if (this != &o) {
        m_tree = std::move(o.m_tree);
      }
      return *this;
    }

    iterator begin() const { return m_tree.begin(); }

    iterator end() const { return m_tree.end(); }

    size_t size() const { return m_tree.size(); }

    // Join
    array_state_map_t join(array_state_map_t &o, cell_varmap_t &cvm_left,
                           smashed_varmap_t &svm_left, base_domain_t &dom_left,
                           cell_varmap_t &cvm_right,
                           smashed_varmap_t &svm_right,
                           base_domain_t &dom_right) {
      join_op op(cvm_left, svm_left, dom_left, cvm_right, svm_right, dom_right);
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree);
      return array_state_map_t(std::move(res));
    }

    // Meet
    array_state_map_t meet(array_state_map_t &o, cell_varmap_t &cvm_left,
                           smashed_varmap_t &svm_left, base_domain_t &dom_left,
                           cell_varmap_t &cvm_right,
                           smashed_varmap_t &svm_right,
                           base_domain_t &dom_right) {
      meet_op op(cvm_left, svm_left, dom_left, cvm_right, svm_right, dom_right);
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree);
      return array_state_map_t(std::move(res));
    }

    void set(variable_t k, array_state v) { m_tree.insert(k, v); }

    array_state_map_t &operator-=(const variable_t &k) {
      m_tree.remove(k);
      return *this;
    }

    const array_state *find(const variable_t &k) const {
      return m_tree.find(k);
    }

    void write(crab::crab_os &o) const {
      o << "{";
      for (auto it = m_tree.begin(); it != m_tree.end();) {
        variable_t k = it->first;
        k.write(o);
        o << " -> ";
        array_state v = it->second;
        v.write(o);
        ++it;
        if (it != this->_tree.end()) {
          o << "; ";
        }
      }
      o << "}";
    }

    friend crab::crab_os &operator<<(crab::crab_os &o,
                                     const array_state_map_t &m) {
      m.write(o);
      return o;
    }
  }; // class array_state_map_t

  // -- scalar domain containing scalar variables, synthetic scalar
  // -- variables from cells and synthetic summarized smashed
  // -- variables.
  base_domain_t m_inv;

  // -- map an array variable to its synthetic cells
  array_state_map_t m_array_map;

  // -- map a synthetic cell to a scalar variable used in m_inv
  cell_varmap_t m_cell_varmap;

  // -- map an array variable to its smashed version used in m_inv
  smashed_varmap_t m_smashed_varmap;

private:
  const array_state &lookup_array_state(const variable_t &v) {
    if (is_bottom()) {
      CRAB_ERROR("cannot call lookup_array_state on bottom");
    }
    const array_state *as = m_array_map.find(v);
    if (as) {
      return *as;
    }
    array_state s;
    m_array_map.set(v, s);
    as = m_array_map.find(v);
    if (!as) {
      CRAB_ERROR("array_state::lookup_array_state returned null");
    }
    return *as;
  }

  static std::string mk_scalar_name(varname_t a, offset_t o, uint64_t size) {
    crab::crab_string_os os;
    os << a << "[";
    if (size == 1) {
      os << o;
    } else {
      os << o << "..." << o.index() + size - 1;
    }
    os << "]";
    return os.str();
  }

  static type_t get_array_element_type(type_t array_type) {
    if (array_type == ARR_BOOL_TYPE) {
      return BOOL_TYPE;
    } else if (array_type == ARR_INT_TYPE) {
      return INT_TYPE;
    } else if (array_type == ARR_REAL_TYPE) {
      return REAL_TYPE;
    } else {
      assert(array_type == ARR_PTR_TYPE);
      return PTR_TYPE;
    }
  }

  cell_t mk_cell(variable_t a, offset_t o, uint64_t sz, offset_map_t &om) {
    // create first an unnamed cell
    cell_t c = om.mk_cell(o, sz);
    assert(!c.is_null());
    // assign a scalar variable to the cell
    auto &vfac = a.name().get_var_factory();
    std::string vname = mk_scalar_name(a.name(), o, sz);
    type_t vtype = get_array_element_type(a.get_type());
    variable_t scalar_var(vfac.get(vname), vtype, sz);
    m_cell_varmap.insert({{a, c}, scalar_var});
    // return the cell
    return c;
  }

  typedef boost::optional<variable_t> variable_opt_t;
  variable_opt_t
  get_scalar(const std::pair<const variable_t &, const cell_t &> &p) {
    if (!p.first.is_array_type()) {
      CRAB_ERROR("array_adaptive::get_scalar only if array variable");
    }
    auto it = m_cell_varmap.find(p);
    if (it != m_cell_varmap.end()) {
      return it->second;
    } else {
      return variable_opt_t();
    }
  }

  static variable_t get_smashed_variable(variable_t a, smashed_varmap_t &svm) {
    if (!a.is_array_type()) {
      CRAB_ERROR(
          "array_adaptive::get_smashed_variable only takes array variables");
    }

    auto mk_smashed_variable = [](variable_t v) {
      assert(v.is_array_type());
      auto &vfac = v.name().get_var_factory();
      crab::crab_string_os os;
      os << "smashed(" << v << ")";
      return variable_t(vfac.get(os.str()), v.get_type());
    };

    auto it = svm.find(a);
    if (it != svm.end()) {
      return it->second;
    } else {
      variable_t smashed_var = mk_smashed_variable(a);
      svm.insert({a, smashed_var});
      return smashed_var;
    }
  }

  std::vector<variable_t> get_array_variables() const {
    std::vector<variable_t> res;
    res.reserve(m_array_map.size());
    for (auto it = m_array_map.begin(), et = m_array_map.end(); it != et;
         ++it) {
      res.push_back(it->first);
    }
    return res;
  }

  void forget_array(const variable_t &v) {
    if (!v.is_array_type()) {
      CRAB_ERROR("cannot call forget_array on a non-array variable");
    }

    const array_state &as = lookup_array_state(v);
    if (!as.is_smashed()) {
      /// We extract all the synthetic cells from the array and forget
      /// them from the underlying abstract domain.
      const offset_map_t &om = as.get_offset_map();
      std::vector<cell_t> cells = om.get_all_cells();
      for (auto &c : cells) {
        variable_opt_t v_opt = get_scalar({v, c});
        if (v_opt) {
          m_inv -= *v_opt;
        }
      }
    } else {
      m_inv -= get_smashed_variable(v, m_smashed_varmap);
    }
    m_array_map -= v;
  }

  interval_t to_interval(linear_expression_t expr, base_domain_t inv) {
    interval_t r(expr.constant());
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      interval_t c(it->first);
      r += c * inv[it->second];
    }
    return r;
  }

  interval_t to_interval(linear_expression_t expr) {
    return to_interval(expr, m_inv);
  }

  void kill_cells(const variable_t &a, const std::vector<cell_t> &cells,
                  offset_map_t &offset_map) {

    assert(a.is_array_type());

    if (!cells.empty()) {
      // Forget the scalars from the numerical domain
      for (unsigned i = 0, e = cells.size(); i < e; ++i) {
        const cell_t &c = cells[i];
        if (variable_opt_t c_scalar_opt = get_scalar({a, c})) {
          m_inv -= *c_scalar_opt;
        }
      }
      if (!Params::is_smashable) {
        // Delete completely the cells. If needed again they they will
        // be re-created.
        for (unsigned i = 0, e = cells.size(); i < e; ++i) {
          const cell_t &c = cells[i];
          offset_map.erase(c);
          m_cell_varmap.erase({a, c});
        }
      } else {
        // if an array is smashable then we don't delete the cells
        // from the offset map. Otherwise, smashing might be unsound.
        // Instead, we mark them as "removed" so they cannot overlap
        // with other cells.
        offset_map.remove(cells);
      }
    }
  }

  // Helper that assign rhs to lhs by switching to the version with
  // the right type.
  void do_assign(variable_t lhs, variable_t rhs) {
    if (lhs.get_type() != rhs.get_type()) {
      CRAB_ERROR("array_adaptive assignment with different types");
    }
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      m_inv.assign_bool_var(lhs, rhs, false);
      break;
    case INT_TYPE:
    case REAL_TYPE:
      m_inv.assign(lhs, rhs);
      break;
    case PTR_TYPE:
      m_inv.pointer_assign(lhs, rhs, number_t(0));
      break;
    default:;
      CRAB_ERROR("array_adaptive assignment with unexpected type");
    }
  }

  // helper to assign a cell into a variable
  void do_assign(variable_t lhs, const variable_t &a, cell_t rhs_c) {
    if (!a.is_array_type()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }

    if (variable_opt_t rhs_v_opt = get_scalar({a, rhs_c})) {
      do_assign(lhs, *rhs_v_opt);
    } else {
      CRAB_LOG("array-adaptive",
               CRAB_WARN("array_adaptive cell without scalar in do_assign"););
      m_inv -= lhs;
    }
  }

  // helper to assign a linear expression into a cell
  void do_assign(const variable_t &a, cell_t lhs_c, linear_expression_t v) {
    if (!a.is_array_type()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }
    variable_opt_t lhs_v_opt = get_scalar({a, lhs_c});
    if (!lhs_v_opt) {
      CRAB_LOG("array-adaptive",
               CRAB_WARN("array_adaptive cell without scalar in do_assign"););
      return;
    }

    variable_t lhs = *lhs_v_opt;
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      if (v.is_constant()) {
        if (v.constant() >= number_t(1)) {
          m_inv.assign_bool_cst(lhs, linear_constraint_t::get_true());
        } else {
          m_inv.assign_bool_cst(lhs, linear_constraint_t::get_false());
        }
      } else if (auto var = v.get_variable()) {
        m_inv.assign_bool_var(lhs, (*var), false);
      }
      break;
    case INT_TYPE:
    case REAL_TYPE:
      m_inv.assign(lhs, v);
      break;
    case PTR_TYPE:
      if (v.is_constant() && v.constant() == number_t(0)) {
        m_inv.pointer_mk_null(lhs);
      } else if (auto var = v.get_variable()) {
        m_inv.pointer_assign(lhs, (*var), number_t(0));
      }
      break;
    default:;
      CRAB_ERROR("array_adaptive assignment with unexpected type");
    }
  }

  // Helper that assign backward rhs to lhs by switching to the
  // version with the right type.
  void do_backward_assign(variable_t lhs, variable_t rhs, base_domain_t &dom) {
    if (lhs.get_type() != rhs.get_type()) {
      CRAB_ERROR("array_adaptive backward assignment with different types");
    }
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      m_inv.backward_assign_bool_var(lhs, rhs, false, dom);
      break;
    case INT_TYPE:
    case REAL_TYPE:
      m_inv.backward_assign(lhs, rhs, dom);
      break;
    case PTR_TYPE:
      CRAB_WARN("array_adaptive backward pointer assignment not implemented");
      break;
    default:;
      CRAB_ERROR("array_adaptive backward_assignment with unexpected type");
    }
  }

  // helper to assign backward a cell into a variable
  void do_backward_assign(variable_t lhs, const variable_t &a, cell_t rhs_c,
                          base_domain_t &dom) {
    if (!a.is_array_type()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }
    variable_opt_t rhs_v_opt = get_scalar({a, rhs_c});
    if (!rhs_v_opt) {
      CRAB_LOG(
          "array-adaptive",
          CRAB_WARN(
              "array_adaptive cell without scalar in do_backward_assign"););
      return;
    }
    do_backward_assign(lhs, *rhs_v_opt, dom);
  }

  // helper to assign backward a linear expression into a cell
  void do_backward_assign(const variable_t &a, cell_t lhs_c,
                          linear_expression_t v, base_domain_t &dom) {
    if (!a.is_array_type()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }
    variable_opt_t lhs_v_opt = get_scalar({a, lhs_c});
    if (!lhs_v_opt) {
      CRAB_LOG(
          "array-adaptive",
          CRAB_WARN(
              "array_adaptive cell without scalar in do_backward_assign"););
      return;
    }
    variable_t lhs = *lhs_v_opt;
    switch (lhs.get_type()) {
    case BOOL_TYPE:
      if (v.is_constant()) {
        if (v.constant() >= number_t(1)) {
          m_inv.backward_assign_bool_cst(lhs, linear_constraint_t::get_true(),
                                         dom);
        } else {
          m_inv.backward_assign_bool_cst(lhs, linear_constraint_t::get_false(),
                                         dom);
        }
      } else if (auto var = v.get_variable()) {
        m_inv.backward_assign_bool_var(lhs, (*var), false, dom);
      }
      break;
    case INT_TYPE:
    case REAL_TYPE:
      m_inv.backward_assign(lhs, v, dom);
      break;
    case PTR_TYPE:
      CRAB_WARN("array_adaptive backward pointer assignment not implemented");
      break;
    default:;
      CRAB_ERROR("array_adaptive backward assignment with unexpected type");
    }
  }

  // The internal representation contains summarized variables of
  // array type and add them as dimensions in the underlying numerical
  // domain. This is OK but it shouldn't be exposed outside via linear
  // constraints.
  //
  // XXX: we should also probably filter out scalar variables
  // originated from cells.
  linear_constraint_system_t
  filter_nonscalar_vars(linear_constraint_system_t &&csts) {
    linear_constraint_system_t res;
    for (auto &cst : csts) {
      auto vars = cst.variables();
      if (std::all_of(vars.begin(), vars.end(), [](const variable_t &v) {
            return v.is_int_type() || v.is_bool_type();
          })) {
        res += cst;
      }
    }
    return res;
  }

  uint64_t check_and_get_elem_size(linear_expression_t elem_size) {
    interval_t i_elem_size = to_interval(elem_size);
    if (boost::optional<number_t> n_bytes = i_elem_size.singleton()) {
      if (static_cast<int64_t>(*n_bytes) > 0 &&
          static_cast<int64_t>(*n_bytes) <=
              std::numeric_limits<uint64_t>::max()) {
        return (uint64_t) static_cast<int64_t>(*n_bytes);
      }
    }
    CRAB_ERROR("array adaptive domain expects constant array element sizes ",
               "between 1 and ", std::numeric_limits<uint64_t>::max(),
               ". Found ", elem_size);
  }

  void do_renaming_for_join(base_domain_t &left_dom, base_domain_t &right_dom,
                            smashed_varmap_t &left_svm,
                            smashed_varmap_t &right_svm,
                            cell_varmap_t &left_cvm, cell_varmap_t &right_cvm,
                            smashed_varmap_t &out_svm, cell_varmap_t &out_cvm) {

    std::vector<variable_t> old_vars_left, old_vars_right, new_vars;

    // Common renaming for smashed scalars
    for (auto &kv : left_svm) {
      variable_t &v1 = kv.second;
      auto &vfac = v1.name().get_var_factory();
      auto it = right_svm.find(kv.first);
      if (it != right_svm.end()) {
        variable_t &v2 = it->second;
        if (v1 != v2) {
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                          v1.get_bitwidth());
          old_vars_left.push_back(v1);
          old_vars_right.push_back(v2);
          new_vars.push_back(outv);
          out_svm.insert({kv.first, outv});
        } else {
          out_svm.insert(kv);
        }
      }
    }

    left_dom.rename(old_vars_left, new_vars);
    right_dom.rename(old_vars_right, new_vars);
    
    old_vars_left.clear();
    old_vars_right.clear();
    new_vars.clear();

    // Common renaming for cell scalars
    for (auto &kv : left_cvm) {
      variable_t &v1 = kv.second;
      auto &vfac = v1.name().get_var_factory();
      auto it = right_cvm.find(kv.first);
      if (it != right_cvm.end()) {
        variable_t &v2 = it->second;
        if (v1 != v2) {
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                          v1.get_bitwidth());
          old_vars_left.push_back(v1);
          old_vars_right.push_back(v2);
          new_vars.push_back(outv);
          out_cvm.insert({kv.first, outv});
        } else {
          out_cvm.insert({kv.first, v1});
        }
      }
    }

    /// Rename the base domains
    left_dom.rename(old_vars_left, new_vars);
    right_dom.rename(old_vars_right, new_vars);
  }

  void do_renaming_for_meet(base_domain_t &left_dom, base_domain_t &right_dom,
                            smashed_varmap_t &left_svm,
                            smashed_varmap_t &right_svm,
                            cell_varmap_t &left_cvm, cell_varmap_t &right_cvm,
                            smashed_varmap_t &out_svm, cell_varmap_t &out_cvm) {

    std::vector<variable_t> old_vars_left, old_vars_right, new_vars;

    /* Figure out common renaming for smashed scalars */

    // Add all mappings from the left operand
    for (auto &kv : left_svm) {
      variable_t &v1 = kv.second;
      auto &vfac = v1.name().get_var_factory();
      auto it = right_svm.find(kv.first);
      if (it != right_svm.end()) {
        variable_t &v2 = it->second;
        if (v1 != v2) {
          // same key but different scalar -> create a fresh common scalar
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                          v1.get_bitwidth());
          old_vars_left.push_back(v1);
          old_vars_right.push_back(v2);
          new_vars.push_back(outv);
          out_svm.insert({kv.first, outv});
          continue;
        }
      }
      out_svm.insert(kv);
    }
    // Add the rest of mappings from the right operand
    for (auto &kv : right_svm) {
      auto it = left_svm.find(kv.first);
      if (it == left_svm.end()) {
        out_svm.insert(kv);
      }
    }
    left_dom.rename(old_vars_left, new_vars);
    right_dom.rename(old_vars_right, new_vars);
 
    old_vars_left.clear();
    old_vars_right.clear();
    new_vars.clear();

    /* Figure out common renaming for cell scalars */

    // Add all mappings from the left operand
    for (auto &kv : left_cvm) {
      variable_t &v1 = kv.second;
      auto &vfac = v1.name().get_var_factory();
      auto it = right_cvm.find(kv.first);
      if (it != right_cvm.end()) {
        variable_t &v2 = it->second;
        if (v1 != v2) {
          // same key but different scalar -> create a fresh common scalar
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                          v1.get_bitwidth());
          old_vars_left.push_back(v1);
          old_vars_right.push_back(v2);
          new_vars.push_back(outv);
          out_cvm.insert({kv.first, outv});
          continue;
        }
      }
      out_cvm.insert(kv);
    }

    // Add the rest of mappings from the right operand
    for (auto &kv : right_cvm) {
      auto it = left_cvm.find(kv.first);
      if (it == left_cvm.end()) {
        out_cvm.insert(kv);
      }
    }

    left_dom.rename(old_vars_left, new_vars);
    right_dom.rename(old_vars_right, new_vars);
  }

  array_adaptive_domain(base_domain_t &&inv, array_state_map_t &&amap,
                        cell_varmap_t &&cvarmap, smashed_varmap_t &&svarmap)
      : m_inv(std::move(inv)), m_array_map(std::move(amap)),
        m_cell_varmap(std::move(cvarmap)),
        m_smashed_varmap(std::move(svarmap)) {}

public:
  array_adaptive_domain(bool is_bottom = false)
      : m_inv(is_bottom ? base_domain_t::bottom() : base_domain_t::top()) {}

  void set_to_top() {
    array_adaptive_domain abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    array_adaptive_domain abs(true);
    std::swap(*this, abs);
  }

  array_adaptive_domain(const array_adaptive_domain_t &other)
      : m_inv(other.m_inv), m_array_map(other.m_array_map),
        m_cell_varmap(other.m_cell_varmap),
        m_smashed_varmap(other.m_smashed_varmap) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  array_adaptive_domain(const array_adaptive_domain_t &&other)
      : m_inv(std::move(other.m_inv)),
        m_array_map(std::move(other.m_array_map)),
        m_cell_varmap(std::move(other.m_cell_varmap)),
        m_smashed_varmap(std::move(other.m_smashed_varmap)) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  array_adaptive_domain_t &operator=(const array_adaptive_domain_t &other) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &other) {
      m_inv = other.m_inv;
      m_array_map = other.m_array_map;
      m_cell_varmap = other.m_cell_varmap;
      m_smashed_varmap = other.m_smashed_varmap;
    }
    return *this;
  }

  array_adaptive_domain_t &operator=(const array_adaptive_domain_t &&other) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &other) {
      m_inv = std::move(other.m_inv);
      m_array_map = std::move(other.m_array_map);
      m_cell_varmap = std::move(other.m_cell_varmap);
      m_smashed_varmap = std::move(other.m_smashed_varmap);
    }
    return *this;
  }

  bool is_bottom() { return (m_inv.is_bottom()); }

  bool is_top() { return (m_inv.is_top()); }

  bool operator<=(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");

    if (is_bottom()) {
      return true;
    } else if (other.is_top()) {
      return true;
    } else {

      CRAB_LOG("array-adaptive",
               crab::outs() << "Check if " << *this << " <= " << other << "\n");

      base_domain_t left_dom(m_inv);
      std::vector<variable_t> old_vars_left, old_vars_right, new_vars;

      /* Figure out common renaming */

      for (auto &kv : m_smashed_varmap) {
        variable_t &v1 = kv.second;
        auto &vfac = v1.name().get_var_factory();
        auto it = other.m_smashed_varmap.find(kv.first);
        // cell exists in both
        if (it != other.m_smashed_varmap.end()) {
          variable_t &v2 = it->second;
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          // same name and type but different variable id
          if (v1 != v2) {
            variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                            v1.get_bitwidth());
            old_vars_left.push_back(v1);
            old_vars_right.push_back(v2);
            new_vars.push_back(outv);
          }
        }
      }
      left_dom.rename(old_vars_left, new_vars);
      other.m_inv.rename(old_vars_right, new_vars);
      
      old_vars_left.clear();
      old_vars_right.clear();
      new_vars.clear();

      for (auto &kv : m_cell_varmap) {
        variable_t &v1 = kv.second;
        auto &vfac = v1.name().get_var_factory();
        auto it = other.m_cell_varmap.find(kv.first);
        // cell exists in both
        if (it != other.m_cell_varmap.end()) {
          variable_t &v2 = it->second;
          assert(v1.name().str() == v2.name().str());
          assert(v1.get_type() == v2.get_type());
          assert(v1.get_bitwidth() == v2.get_bitwidth());
          // same name and type but different variable id
          if (v1 != v2) {
            variable_t outv(vfac.get(v1.name().str()), v1.get_type(),
                            v1.get_bitwidth());
            old_vars_left.push_back(v1);
            old_vars_right.push_back(v2);
            new_vars.push_back(outv);
          }
        }
      }
      left_dom.rename(old_vars_left, new_vars);
      other.m_inv.rename(old_vars_right, new_vars);
      
      // We need to be careful if one array state is smashed and the
      // other is not.
      bool res = (left_dom <= other.m_inv);
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  bool operator==(array_adaptive_domain_t other) {
    return (m_inv <= other.m_inv && other.m_inv <= m_inv);
  }

  void operator|=(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    CRAB_LOG("array-adaptive",
	     crab::outs() << "Join " << *this << " and " << other << "\n";);
    
    if (other.is_bottom() || is_top()) {
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);      
      return;
    } else if (is_bottom() || other.is_top()) {
      *this = other;
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);      
    } else {

      // this must be done before the renaming
      m_array_map = std::move(m_array_map.join(
          other.m_array_map, m_cell_varmap, m_smashed_varmap, m_inv,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv));

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;
      do_renaming_for_join(m_inv, other.m_inv, m_smashed_varmap,
                           other.m_smashed_varmap, m_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);
      m_inv |= other.m_inv;
      std::swap(m_cell_varmap, out_cell_varmap);
      std::swap(m_smashed_varmap, out_smashed_varmap);
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
    }
  }

  array_adaptive_domain_t operator|(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    if (other.is_bottom() || is_top()) {
      return *this;
    } else if (is_bottom() || other.is_top()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive",
               crab::outs() << "Join " << *this << " and " << other << "\n";);

      base_domain_t left_dom(m_inv);
      cell_varmap_t left_cell_varmap(m_cell_varmap);
      smashed_varmap_t left_smashed_varmap(m_smashed_varmap);

      // Must be done before the renaming.
      auto out_array_map = std::move(m_array_map.join(
          other.m_array_map, left_cell_varmap, left_smashed_varmap, left_dom,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv));

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;
      do_renaming_for_join(left_dom, other.m_inv, left_smashed_varmap,
                           other.m_smashed_varmap, left_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);

      array_adaptive_domain_t res(left_dom | other.m_inv,
				  std::move(out_array_map),
                                  std::move(out_cell_varmap),
                                  std::move(out_smashed_varmap));

      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t operator&(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive",
               crab::outs() << "Meet " << *this << " and " << other << "\n";);

      base_domain_t left_dom(m_inv);
      cell_varmap_t left_cell_varmap(m_cell_varmap);
      smashed_varmap_t left_smashed_varmap(m_smashed_varmap);

      // Must be done before the renaming.
      auto out_array_map = m_array_map.meet(
          other.m_array_map, left_cell_varmap, left_smashed_varmap, left_dom,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv);

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;
      do_renaming_for_meet(left_dom, other.m_inv, left_smashed_varmap,
                           other.m_smashed_varmap, left_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);

      array_adaptive_domain_t res(left_dom & other.m_inv, std::move(out_array_map),
                                  std::move(out_cell_varmap),
                                  std::move(out_smashed_varmap));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t operator||(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    if (other.is_bottom()) {
      return *this;
    } else if (is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Widening " << *this << " and "
                                              << other << "\n";);

      base_domain_t left_dom(m_inv);
      cell_varmap_t left_cell_varmap(m_cell_varmap);
      smashed_varmap_t left_smashed_varmap(m_smashed_varmap);

      // Must be done before the renaming.
      auto out_array_map = m_array_map.join(
          other.m_array_map, left_cell_varmap, left_smashed_varmap, left_dom,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv);

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;

      do_renaming_for_join(left_dom, other.m_inv, left_smashed_varmap,
                           other.m_smashed_varmap, left_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);

      array_adaptive_domain_t res(left_dom || other.m_inv, std::move(out_array_map),
                                  std::move(out_cell_varmap),
                                  std::move(out_smashed_varmap));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t
  widening_thresholds(array_adaptive_domain_t other,
                      const iterators::thresholds<number_t> &ts) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    if (other.is_bottom()) {
      return *this;
    } else if (is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Widening " << *this << " and "
                                              << other << "\n";);

      base_domain_t left_dom(m_inv);
      cell_varmap_t left_cell_varmap(m_cell_varmap);
      smashed_varmap_t left_smashed_varmap(m_smashed_varmap);

      // Must be done before the renaming.
      auto out_array_map = m_array_map.join(
          other.m_array_map, left_cell_varmap, left_smashed_varmap, left_dom,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv);

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;
      do_renaming_for_join(left_dom, other.m_inv, left_smashed_varmap,
                           other.m_smashed_varmap, left_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);

      array_adaptive_domain_t res(
          left_dom.widening_thresholds(other.m_inv, ts), std::move(out_array_map),
          std::move(out_cell_varmap), std::move(out_smashed_varmap));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t operator&&(array_adaptive_domain_t other) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    if (is_bottom()) {
      return *this;
    } else if (other.is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Narrowing " << *this
                                              << " and " << other << "\n";);

      base_domain_t left_dom(m_inv);
      cell_varmap_t left_cell_varmap(m_cell_varmap);
      smashed_varmap_t left_smashed_varmap(m_smashed_varmap);

      // Must be done before the renaming.
      auto out_array_map = m_array_map.join(
          other.m_array_map, left_cell_varmap, left_smashed_varmap, left_dom,
          other.m_cell_varmap, other.m_smashed_varmap, other.m_inv);

      smashed_varmap_t out_smashed_varmap;
      cell_varmap_t out_cell_varmap;

      do_renaming_for_meet(left_dom, other.m_inv, left_smashed_varmap,
                           other.m_smashed_varmap, left_cell_varmap,
                           other.m_cell_varmap, out_smashed_varmap,
                           out_cell_varmap);

      array_adaptive_domain_t res(left_dom && other.m_inv, std::move(out_array_map),
                                  std::move(out_cell_varmap),
                                  std::move(out_smashed_varmap));

      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  void forget(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    variable_vector_t scalar_variables;
    scalar_variables.reserve(variables.size());
    for (variable_t v : variables) {
      if (v.is_array_type()) {
	CRAB_LOG("array-adaptive",
		 crab::outs() << "Forget array variable " << v << "\n";);
        forget_array(v);
      } else {
	CRAB_LOG("array-adaptive",
		 crab::outs() << "Forget scalar variable " << v << "\n";);	
        scalar_variables.push_back(v);
      }
    }
    m_inv.forget(scalar_variables);
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    // if we must keep array variable v then we need to keep all its
    // synthetic cells in m_inv.

    variable_vector_t keep_vars;
    std::set<variable_t> keep_array_varset;
    for (variable_t v : variables) {
      if (v.is_array_type()) {
        keep_array_varset.insert(v);
      } else {
        keep_vars.push_back(v);
      }
    }

    std::vector<variable_t> array_variables = get_array_variables();
    for (variable_t v : array_variables) {
      if (keep_array_varset.count(v) > 0) {
        // keep all cells of v
        const array_state &as = lookup_array_state(v);
        if (!as.is_smashed()) {
          const offset_map_t &om = as.get_offset_map();
          std::vector<cell_t> cells = om.get_all_cells();
          for (auto &c : cells) {
            variable_opt_t v_opt = get_scalar({v, c});
            if (v_opt) {
              keep_vars.push_back(*v_opt);
            }
          }
        } else {
          keep_vars.push_back(get_smashed_variable(v, m_smashed_varmap));
        }
      } else {
        // no need to remove any cell from m_inv because m_inv.project will do
        m_array_map -= v;
      }
    }

    // Finally we project
    m_inv.project(keep_vars);
  }

  void expand(variable_t var, variable_t new_var) {
    CRAB_WARN("array adaptive expand not implemented");
  }

  void normalize() { CRAB_WARN("array adaptive normalize not implemented"); }

  void minimize() { m_inv.minimize(); }

  interval_t operator[](variable_t v) {
    if (!v.is_array_type()) {
      return m_inv[v];
    } else {
      return interval_t::top();
    }
  }
  
  /* begin intrinsics operations */  
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    if ((std::all_of(inputs.begin(), inputs.end(),
		    [](const variable_t &v) {
		      return v.is_int_type() || v.is_bool_type();
		    })) &&
	(std::all_of(outputs.begin(), outputs.end(),
		     [](const variable_t &v) {
		       return v.is_int_type() || v.is_bool_type();
		     }))) {
      m_inv.intrinsic(name, inputs, outputs);
    } else { 
      CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
    }
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  array_adaptive_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

    m_inv += csts;

    CRAB_LOG("array-adaptive",
             crab::outs() << "assume(" << csts << ")  " << *this << "\n";);
  }

  void operator-=(variable_t var) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom()) {
      return;
    }
    
    if (var.is_array_type()) {
      forget_array(var);
    } else {
      m_inv -= var;
    }
  }

  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    m_inv.assign(x, e);

    CRAB_LOG("array-adaptive", crab::outs() << "apply " << x << " := " << e
                                            << " " << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    m_inv.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    m_inv.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void backward_assign(variable_t x, linear_expression_t e,
                       array_adaptive_domain_t inv) {
    m_inv.backward_assign(x, e, inv.get_content_domain());
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      array_adaptive_domain_t inv) {
    m_inv.backward_apply(op, x, y, z, inv.get_content_domain());
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      array_adaptive_domain_t inv) {
    m_inv.backward_apply(op, x, y, z, inv.get_content_domain());
  }

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    m_inv.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    m_inv.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    m_inv.apply(op, x, y, k);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << k << " " << *this << "\n";);
  }

  // boolean operators
  virtual void assign_bool_cst(variable_t lhs,
                               linear_constraint_t rhs) override {
    m_inv.assign_bool_cst(lhs, rhs);
  }

  virtual void assign_bool_var(variable_t lhs, variable_t rhs,
                               bool is_not_rhs) override {
    m_inv.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  virtual void apply_binary_bool(bool_operation_t op, variable_t x,
                                 variable_t y, variable_t z) override {
    m_inv.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(variable_t v, bool is_negated) override {
    m_inv.assume_bool(v, is_negated);
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                        array_adaptive_domain_t inv) {
    m_inv.backward_assign_bool_cst(lhs, rhs, inv.get_content_domain());
  }

  virtual void backward_assign_bool_var(variable_t lhs, variable_t rhs,
                                        bool is_not_rhs,
                                        array_adaptive_domain_t inv) {
    m_inv.backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                   inv.get_content_domain());
  }

  virtual void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                          variable_t y, variable_t z,
                                          array_adaptive_domain_t inv) {
    m_inv.backward_apply_binary_bool(op, x, y, z, inv.get_content_domain());
  }

  // pointer_operators_api
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    m_inv.pointer_load(lhs, rhs);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    m_inv.pointer_store(lhs, rhs);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    m_inv.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    m_inv.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    m_inv.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    m_inv.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(ptr_cst_t cst) override {
    m_inv.pointer_assume(cst);
  }

  virtual void pointer_assert(ptr_cst_t cst) override {
    m_inv.pointer_assert(cst);
  }

  // array_operators_api

  // array_init returns a fresh array where all elements between
  // lb_idx and ub_idx are initialized to val. Thus, the first thing
  // we need to do is to kill existing cells.
  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.array_init");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_init");

    if (is_bottom())
      return;

    const array_state &as = lookup_array_state(a);
    // The array shouldn't be smashed yet
    if (!as.is_smashed()) {
      array_state next_as(as); // important to make the copy
      offset_map_t &om = next_as.get_offset_map();
      std::vector<cell_t> old_cells = om.get_all_cells();
      if (!old_cells.empty()) {
        kill_cells(a, old_cells, om);
        m_array_map.set(a, next_as);
      }
    }

    array_store_range(a, elem_size, lb_idx, ub_idx, val);
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {
    crab::CrabStats::count(getDomainName() + ".count.array_load");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_load");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);
    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      // Check smashed array is consistent with elem_size
      const cp_domain_t &a_elem_size = as.get_element_sz();
      cp_domain_t cp_e_sz((int64_t)e_sz);
      cp_e_sz |= a_elem_size;
      if (!cp_e_sz.is_top()) {
        variable_t smashed_a = get_smashed_variable(a, m_smashed_varmap);
        m_inv.array_load(lhs, smashed_a, elem_size, i);
        goto array_load_end;
      } else {
        // lhs will be forgotten
      }
    } else {
      interval_t ii = to_interval(i);
      if (boost::optional<number_t> n = ii.singleton()) {
        array_state next_as(as); // important to make the copy
        offset_map_t &offset_map = next_as.get_offset_map();
        offset_t o(static_cast<int64_t>(*n));
        std::vector<cell_t> cells;
        offset_map.get_overlap_cells(o, e_sz, cells);
	CRAB_LOG("array-adaptive",	
		 crab::outs() << "Number of overlapping cells=" << cells.size() << "\n";);
	
        if (!cells.empty()) {
          CRAB_LOG("array-adaptive",
                   CRAB_WARN("Ignored read from cell ", a, "[", o, "...",
                             o.index() + e_sz - 1, "]",
                             " because it overlaps with ", cells.size(),
                             " cells"););
          /*
            TODO: we can apply here "Value Recomposition" 'a la'
            Mine'06 to construct values of some type from a sequence
            of bytes. It can be endian-independent but it would more
            precise if we choose between little- and big-endian.
          */
        } else {
          cell_t c = mk_cell(a, o, e_sz, offset_map);
          // Here it's ok to do assignment (instead of expand)
          // because c is not a summarized variable. Otherwise, it
          // would be unsound.
          do_assign(lhs, a, c);
          m_array_map.set(a, next_as);
          goto array_load_end;
        }
      } else {
        linear_expression_t symb_lb(i);
        linear_expression_t symb_ub(i + number_t(e_sz - 1));
        std::vector<cell_t> cells;
        const offset_map_t &offset_map = as.get_offset_map();
        offset_map.get_overlap_cells_symbolic_offset(m_inv, symb_lb, symb_ub,
                                                     cells);
        // XXX: if we have a large array that is never smashed but we
        // do many reads with symbolic offsets then it might be better
        // to smash the array so that each read is cheaper.
        if (Params::is_smashable) {
          if (array_state::can_be_smashed(cells, e_sz, true)) {
            // we smash all overlapping cells into a fresh array
            // (summarized) variable
            auto &vfac = a.name().get_var_factory();
            variable_t fresh_var(vfac.get(), a.get_type());
            bool found_cell_without_scalar = false;
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              auto c_scalar_opt = get_scalar({a, c});
              if (!c_scalar_opt) {
		CRAB_LOG("array-adaptive",
			 CRAB_WARN("array adaptive: ignored array load from ", a,
				   " because non-constant array index ", i, "=", ii,
				   " because found unnamed cell ", c););
		
                found_cell_without_scalar = true;
                break;
              }
              const bool is_strong_update = (k == 0);
              m_inv.array_store(fresh_var, elem_size, i, *c_scalar_opt,
                                is_strong_update);
            }
            if (found_cell_without_scalar) {
              m_inv -= lhs;
            } else {
              // we read from the temporary summarized variable
              m_inv.array_load(lhs, fresh_var, elem_size, i);
            }
            // we forget the temporary summarized variable
            m_inv -= fresh_var;
            goto array_load_end;
          } else {
            CRAB_LOG("array-adaptive",
                     CRAB_WARN("array adaptive: ignored array load from ", a,
                               " because non-constant array index ", i, "=", ii,
                               " and cannot smash the array"););
          }
        } else {
          CRAB_LOG("array-adaptive",
                   CRAB_WARN("array adaptive: ignored array load from ", a,
                             " because non-constant array index ", i, "=",
                             ii););
        }
      }
    }
    m_inv -= lhs;

  array_load_end:
    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << lhs << ":=" << a << "[" << i << "..." << ub
                          << "]  -- " << *this << "\n";);
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool is_strong_update) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);
    const array_state &as = lookup_array_state(a);

    if (as.is_smashed()) {
      variable_t smashed_a = get_smashed_variable(a, m_smashed_varmap);
      const cp_domain_t &a_elem_size = as.get_element_sz();
      cp_domain_t cp_e_sz((int64_t)e_sz);
      cp_e_sz |= a_elem_size;
      if (!cp_e_sz.is_top()) {
        m_inv.array_store(smashed_a, elem_size, i, val, is_strong_update);
      } else {
        m_inv -= smashed_a;
      }
    } else {
      interval_t ii = to_interval(i);
      array_state next_as(as);
      offset_map_t &offset_map = next_as.get_offset_map();
      if (boost::optional<number_t> n = ii.singleton()) {
        // -- Constant index: kill overlapping cells + perform strong update
        std::vector<cell_t> cells;
        offset_t o(static_cast<int64_t>(*n));
        offset_map.get_overlap_cells(o, e_sz, cells);
        if (cells.size() > 0) {
          CRAB_LOG("array-adaptive",
                   CRAB_WARN("Killed ", cells.size(),
                             " overlapping cells with ", "[", o, "...",
                             o.index() + e_sz - 1, "]", " before writing."));
          kill_cells(a, cells, offset_map);
        }
        // Perform scalar update
        // -- create a new cell it there is no one already
        cell_t c = mk_cell(a, o, e_sz, offset_map);
        // -- strong update
        do_assign(a, c, val);
      } else {
        // -- Non-constant index: kill overlapping cells

        CRAB_LOG("array-adaptive", crab::outs() << "array write to " << a
                                                << " with non-constant index "
                                                << i << "=" << ii << "\n";);


        bool smashed = false; // whether smashing took place
        if (Params::is_smashable) {
	  std::vector<cell_t> cells = offset_map.get_all_cells();
          if (next_as.can_be_smashed(e_sz) &&
	      // Smashing is expensive because it will go over all cells
	      // performing one join per weak update even if the smashed
	      // array is already unconstrained. We don't smash if the
	      // number of cells to be smashed is too large.
	      (cells.size() <= Params::max_smashable_cells)) {
            smashed = true;
            CRAB_LOG("array-adaptive-smash",
                     crab::outs() << "Array " << a << " will be smashed\n";);
            bool found_cell_without_scalar = false;
            variable_t smashed_a = get_smashed_variable(a, m_smashed_varmap);
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              auto c_scalar_opt = get_scalar({a, c});
              if (!c_scalar_opt) {
                found_cell_without_scalar = true;
                break;
              }
              const bool is_strong_update = (k == 0);
              m_inv.array_store(smashed_a, elem_size, i, *c_scalar_opt,
                                is_strong_update);
              CRAB_LOG("array-adaptive-smash",
                       crab::outs() << "\tAfter smashing " << *c_scalar_opt
                                    << "=" << m_inv << "\n";);
            }

            if (found_cell_without_scalar) {
              m_inv -= smashed_a;
            } else {
              // Finally the array store
              m_inv.array_store(smashed_a, elem_size, i, val, is_strong_update);
            }

            // The removal of cells from offset_map must be done after
            // the array has been fully smashed.
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              if (auto c_scalar_opt = get_scalar({a, c})) {
                // destroy the synthethic cell from the base domain and
                // from the offset map
                m_inv -= *c_scalar_opt;
              }
              offset_map.erase(c);
              m_cell_varmap.erase({a, c});
            }
            next_as.set_smashed(true);
            next_as.get_element_sz() = cp_domain_t(number_t(e_sz));

            CRAB_LOG("array-adaptive-smash",
                     crab::outs() << "Array " << a
                                  << " has been smashed:" << m_inv << "\n";);
          } else {
            CRAB_LOG("array-adaptive",
		     if (cells.size() > Params::max_smashable_cells) {
		       CRAB_WARN("Array ", a,
				 " cannot be smashed because too many cells ",
				 cells.size(),
				 ". Array write at index ", i, "=", ii,
				 " is safely ignored");
		     } else {
		       CRAB_WARN("Array ", a,
				 " cannot be smashed so array write at index ", i,
				 "=", ii, " is safely ignored");
		     });
          }
        }

        if (!smashed) {
          linear_expression_t symb_lb(i);
          linear_expression_t symb_ub(i + number_t(e_sz - 1));
          std::vector<cell_t> cells;
          offset_map.get_overlap_cells_symbolic_offset(m_inv, symb_lb, symb_ub,
                                                       cells);
          CRAB_LOG(
              "array-adaptive", crab::outs() << "Killed cells: {";
              for (unsigned j = 0; j < cells.size();) {
                crab::outs() << cells[j];
                ++j;
                if (j < cells.size()) {
                  crab::outs() << ",";
                }
              } crab::outs()
              << "}\n";);

          kill_cells(a, cells, offset_map);
        }
      }
      m_array_map.set(a, next_as);
    }
    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << a << "[" << i << "..." << ub << "]:=" << val
                          << " -- " << *this << "\n";);
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t elem_size, linear_expression_t i,
                           linear_expression_t val,
                           bool /*is_strong_update*/) override {
    CRAB_WARN("array_store in the array adaptive domain not implemented");
  }

  // Perform array stores over an array segment [lb_idx, ub_idx]
  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t lb_idx,
                                 linear_expression_t ub_idx,
                                 linear_expression_t val) override {
    crab::CrabStats::count(getDomainName() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() + ".array_store_range");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);
    interval_t lb_i = to_interval(lb_idx);
    auto lb = lb_i.singleton();
    if (!lb) {
      CRAB_WARN("array adaptive store range ignored because ", "lower bound",
                lb_idx, " is not constant");
      return;
    }

    interval_t ub_i = to_interval(ub_idx);
    auto ub = ub_i.singleton();
    if (!ub) {
      CRAB_WARN("array adaptive store range ignored because ", "upper bound ",
                ub_idx, " is not constant");
      return;
    }

    if (!(*lb <= *ub)) {
      CRAB_WARN("array adaptive store range ignored because lower bound ",
		*lb, " is not less or equal than upper bound ", *ub);
      return;
    }
    
    number_t num_elems = (*ub - *lb) / e_sz;
    number_t e = *ub;
    if (num_elems > Params::max_array_size) {
      e = *lb + ((number_t(Params::max_array_size) -1) * e_sz);
      CRAB_WARN("array adaptive store range will ignore indexes greater than ", e);
    }

    for (number_t i = *lb; i <= e;) {
      array_store(a, elem_size, i, val, false);
      i = i + e_sz;
    }
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t lb_idx,
                                 linear_expression_t ub_idx,
                                 linear_expression_t val) override {
    CRAB_WARN("array_store_range in the array adaptive domain not implemented");
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    CRAB_LOG("array-adaptive",
	     crab::outs() << "Array assign " << lhs << " := " << rhs << "\n";);

    if (is_bottom()) {
      return;
    }
    
    const array_state &as = lookup_array_state(rhs);
    if (!as.is_smashed()) {
      offset_map_t lhs_om;
      const offset_map_t &rhs_om = as.get_offset_map();
      std::map<variable_t, variable_t> renmap;
      CRAB_LOG("array-adaptive-array-assign", crab::outs() << "Not smashed\n";);       
      std::vector<cell_t> cells = rhs_om.get_all_cells();
      for (auto &c : cells) {
        variable_opt_t c_scalar_opt = get_scalar({rhs, c});
        if (!c_scalar_opt) {
          continue;
        }
        // Create a new cell for lhs from rhs's cell.
        cell_t new_c = mk_cell(lhs, c.get_offset(), c.get_size(), lhs_om);
        variable_opt_t new_c_scalar_opt = get_scalar({lhs, new_c});
        assert(new_c_scalar_opt);
        renmap.insert({*new_c_scalar_opt, *c_scalar_opt});
      }

      cp_domain_t elem_sz = as.get_element_sz();
      m_array_map.set(
          lhs, array_state(false, std::move(elem_sz), std::move(lhs_om)));

      CRAB_LOG("array-adaptive-array-assign",      
	       crab::outs() << "array variables={";
	       std::vector<variable_t> array_variables = get_array_variables();
	       for(unsigned i=0,e=array_variables.size();i<e;++i) {
		 crab::outs() << array_variables[i] <<";";
	       }
	       crab::outs() << "}\n";);
      
      for (auto &kv : renmap) {
	CRAB_LOG("array-adaptive-array-assign",      	
		 crab::outs() << "Base domain assign " << kv.first << " := " << kv.second << "\n";);
        do_assign(kv.first, kv.second);
      }
    } else if (Params::is_smashable) {
      CRAB_LOG("array-adaptive-array-assign", crab::outs() << "Smashed\n";); 
      variable_t smashed_lhs = get_smashed_variable(lhs, m_smashed_varmap);
      variable_t smashed_rhs = get_smashed_variable(rhs, m_smashed_varmap);
      m_inv.array_assign(smashed_lhs, smashed_rhs);
      m_array_map.set(lhs, as);
    }
    CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
  }

  // backward array operations

  virtual void backward_array_init(variable_t a, linear_expression_t elem_size,
                                   linear_expression_t lb_idx,
                                   linear_expression_t ub_idx,
                                   linear_expression_t val,
                                   array_adaptive_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_init");

    if (is_bottom())
      return;

    // make all array cells uninitialized
    const array_state &as = lookup_array_state(a);
    if (!as.is_smashed()) {
      array_state next_as(as);
      offset_map_t &om = next_as.get_offset_map();
      std::vector<cell_t> old_cells = om.get_all_cells();
      if (!old_cells.empty()) {
        kill_cells(a, old_cells, om);
        m_array_map.set(a, next_as);
      }
    } else {
      CRAB_WARN("array_adaptive::backward_array_init not implemented if array "
                "smashed");
    }

    // meet with forward invariant
    *this = *this & invariant;
  }

  virtual void backward_array_load(variable_t lhs, variable_t a,
                                   linear_expression_t elem_size,
                                   linear_expression_t i,
                                   array_adaptive_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_load");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);

    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      CRAB_WARN("array_adaptive::backward_array_load not implemented if array "
                "smashed");
    } else {
      // XXX: we use the forward invariant to extract the array index
      interval_t ii = to_interval(i, invariant.get_content_domain());
      if (boost::optional<number_t> n = ii.singleton()) {
        array_state next_as(as);
        offset_map_t &om = next_as.get_offset_map();
        offset_t o(static_cast<int64_t>(*n));
        cell_t c = mk_cell(a, o, e_sz, om);
        do_backward_assign(lhs, a, c, invariant.get_content_domain());
        m_array_map.set(a, next_as);
      } else {
        CRAB_LOG("array-adaptive",
                 CRAB_WARN("array index is not a constant value"););
        // -- Forget lhs
        m_inv -= lhs;
        // -- Meet with forward invariant
        *this = *this & invariant;
      }
    }

    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << "BACKWARD " << lhs << ":=" << a << "[" << i
                          << "..." << ub << "]  -- " << *this << "\n";);
  }

  virtual void
  backward_array_store(variable_t a, linear_expression_t elem_size,
                       linear_expression_t i, linear_expression_t val,
                       bool /*is_strong_update*/,
                       array_adaptive_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_array_store");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);

    // XXX: we use the forward invariant to extract the array index
    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      CRAB_WARN("array_adaptive::backward_array_store not implemented if array "
                "smashed");
    } else {
      array_state next_as(as);
      offset_map_t &om = next_as.get_offset_map();
      // XXX: we use the forward invariant to extract the array index
      interval_t ii = to_interval(i, invariant.get_content_domain());
      if (boost::optional<number_t> n = ii.singleton()) {
        // -- Constant index and the store updated one single cell:
        // -- backward assign in the base domain.
        offset_t o(static_cast<int64_t>(*n));
        std::vector<cell_t> cells;
        om.get_overlap_cells(o, e_sz, cells);
        // post: forall c \in cells:: c != [o,e_sz)
        // that is, get_overlap_cells returns cells different from [o, e_sz)
        if (cells.size() >= 1) {
          kill_cells(a, cells, om);
          *this = *this & invariant;
        } else {
          // c might be in m_inv or not.
          cell_t c = mk_cell(a, o, e_sz, om);
          do_backward_assign(a, c, val, invariant.get_content_domain());
        }
      } else {
        // TODOX: smash the array if needed

        // -- Non-constant index or multiple overlapping cells: kill
        // -- overlapping cells and meet with forward invariant.
        linear_expression_t symb_lb(i);
        linear_expression_t symb_ub(i + number_t(e_sz - 1));
        std::vector<cell_t> cells;
        om.get_overlap_cells_symbolic_offset(invariant.get_content_domain(),
                                             symb_lb, symb_ub, cells);
        kill_cells(a, cells, om);

        *this = *this & invariant;
      }
      m_array_map.set(a, next_as);
    }

    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << "BACKWARD " << a << "[" << i << "..." << ub
                          << "]:=" << val << " -- " << *this << "\n";);
  }

  virtual void
  backward_array_store(variable_t a_new, variable_t a_old,
                       linear_expression_t elem_size, linear_expression_t i,
                       linear_expression_t val, bool /*is_strong_update*/,
                       array_adaptive_domain_t invariant) override {
    CRAB_WARN("backward_array_store in array_adaptive domain not implemented");
  }

  virtual void backward_array_store_range(
      variable_t a, linear_expression_t elem_size, linear_expression_t lb_idx,
      linear_expression_t ub_idx, linear_expression_t val,
      array_adaptive_domain_t invariant) override {
    crab::CrabStats::count(getDomainName() +
                           ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".count.backward_array_store_range");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);

    interval_t lb_i = to_interval(lb_idx, invariant.get_content_domain());
    auto lb = lb_i.singleton();
    if (!lb) {
      return;
    }

    interval_t ub_i = to_interval(ub_idx, invariant.get_content_domain());
    auto ub = ub_i.singleton();
    if (!ub) {
      return;
    }

    if (!(*lb <= *ub)) {
      return;
    }
    
    number_t num_elems = (*ub - *lb) / e_sz;
    number_t e = *ub;
    if (num_elems > Params::max_array_size) {
      e = *lb + ((number_t(Params::max_array_size) -1) * e_sz);
    }

    for (number_t i = *lb; i <= e;) {
      backward_array_store(a, elem_size, i, val, false, invariant);
      i = i + e_sz;
    }
  }

  virtual void backward_array_store_range(
      variable_t a_new, variable_t a_old, linear_expression_t elem_size,
      linear_expression_t lb_idx, linear_expression_t ub_idx,
      linear_expression_t val, array_adaptive_domain_t invariant) override {
    CRAB_WARN(
        "backward_array_store_range in array_adaptive domain not implemented");
  }

  virtual void
  backward_array_assign(variable_t lhs, variable_t rhs,
                        array_adaptive_domain_t invariant) override {
    CRAB_WARN("backward_array_assign in array_adaptive domain not implemented");
  }

  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    return filter_nonscalar_vars(
        std::move(m_inv.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    disjunctive_linear_constraint_system_t res;
    auto disj_csts = m_inv.to_disjunctive_linear_constraint_system();
    for (auto &csts : disj_csts) {
      auto filtered_csts = filter_nonscalar_vars(std::move(csts));
      if (!filtered_csts.is_true()) {
        res += filtered_csts;
      }
    }
    return res;
  }

  base_domain_t get_content_domain() const { return m_inv; }

  base_domain_t &get_content_domain() { return m_inv; }

  void write(crab_os &o) {
    o << m_inv;
    CRAB_LOG("array-adaptive-print-details",
	     crab::outs() << "\n";
	     crab::outs () << "=== CELLS PER ARRAY === \n";
	     for (auto it = m_array_map.begin(), et = m_array_map.end(); it!=et; ++it) {
	       const variable_t  &v = it->first;
	       const array_state &as = it->second;
	       crab::outs() << "ARRAY VAR " <<  v << ":\n";
	       as.write(crab::outs());
	       crab::outs() << "\n";
	     }
	     crab::outs() << "=== MAP FROM CELLS TO SCALARS === \n";
	     for (auto const& kv: m_cell_varmap) {
	       crab::outs() << "\t" << kv.first.first << "#" << kv.first.second << " -> "
			    << kv.second << "\n";
	     });
  }

  static std::string getDomainName() {
    std::string name("ArrayAdaptive(" + base_domain_t::getDomainName() + ")");
    return name;
  }

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    CRAB_WARN("array_adaptive::rename not implemented");
    
    // m_inv.rename(from, to);
    // for (auto &v : from) {
    //   if (v.is_array_type()) {
    //     CRAB_WARN("TODO: array_adaptive::rename array variable");
    //   }
    // }
  }

}; // end array_adaptive_domain

template <typename Dom, class Params>
struct abstract_domain_traits<array_adaptive_domain<Dom, Params>> {
  typedef typename Dom::number_t number_t;
  typedef typename Dom::varname_t varname_t;
};

template <typename Dom, class Params>
class checker_domain_traits<array_adaptive_domain<Dom, Params>> {
public:
  typedef array_adaptive_domain<Dom, Params> this_type;
  typedef typename this_type::base_domain_t base_domain_t;
  typedef typename this_type::linear_constraint_t linear_constraint_t;
  typedef typename this_type::disjunctive_linear_constraint_system_t
      disjunctive_linear_constraint_system_t;

  static bool entail(this_type &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    base_domain_t &lhs_dom = lhs.get_content_domain();
    return checker_domain_traits<base_domain_t>::entail(lhs_dom, rhs);
  }

  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     this_type &rhs) {
    base_domain_t &rhs_dom = rhs.get_content_domain();
    return checker_domain_traits<base_domain_t>::entail(lhs, rhs_dom);
  }

  static bool entail(this_type &lhs, const linear_constraint_t &rhs) {
    base_domain_t &lhs_dom = lhs.get_content_domain();
    return checker_domain_traits<base_domain_t>::entail(lhs_dom, rhs);
  }

  static bool intersect(this_type &inv, const linear_constraint_t &cst) {
    base_domain_t &dom = inv.get_content_domain();
    return checker_domain_traits<base_domain_t>::intersect(dom, cst);
  }
};

} // namespace domains
} // namespace crab
