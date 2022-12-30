/*******************************************************************************
 * Adaptive array domain
 *
 * Initially, an array is modeled by mapping sequences of consecutive
 * bytes (segments) to cells. A cell is pair <offset, size> where:
 *
 * - offset is an unsigned number
 * - size is an unsigned number
 *
 * A cell, when associated to an array A, is mapped to a scalar ghost
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
 *
 * Ghost variables:
 *
 * - Each concrete array cell is modeled by a scalar ghost
 *   variable. This is managed by cell_ghost_man.
 *
 * - Ghost variables are created by the same variable factory used to
 *   generate the CrabIR. We don't reuse ghost variables so there is a
 *   small risk of running out of ghost variables. If that happens, a
 *   runtime error will be reported.
 *
 * - POSSIBLE BREAK OF MODULARITY of the base domain: A more important
 *   issue is that currently we pass ghost scalar variables to the
 *   base domain that can change during the analysis. If the base
 *   domain uses these ghost variables to create new ghost variables
 *   we will get into trouble. Currently, these ghost scalars are
 *   passed directly to boolean/numerical abstract domains so we
 *   should be fine with current Crab domains because they don't
 *   create ghost variables but this needs to be revisited anytime a
 *   new domain is plugin.
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/array_smashing.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/indexable.hpp>

#include <algorithm>
#include <functional>
#include <set>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

// forward declaration
template <typename Domain> class array_adaptive_domain;

namespace array_adaptive_impl {

// forward declaration
class offset_map;

/*
 * Wrapper for using ikos::index_t as patricia_tree keys
 */
class offset_t : public indexable {
  ikos::index_t m_val;

public:
  explicit offset_t(ikos::index_t v);

  virtual ikos::index_t index() const override;

  size_t hash() const;

  bool operator<(const offset_t &o) const;

  bool operator==(const offset_t &o) const;

  bool operator!=(const offset_t &o) const;

  offset_t operator%(const offset_t &o) const;

  offset_t operator-(const offset_t &o) const;

  bool is_negative() const;

  bool is_zero() const;

  virtual void write(crab::crab_os &o) const override;

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
  using interval_t = ikos::interval<ikos::z_number>;
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
  cell_t();

  cell_t(offset_t offset, uint64_t size);

  static interval_t to_interval(const offset_t o, uint64_t size);

  interval_t to_interval() const;

public:
  bool is_null() const;

  offset_t get_offset() const;

  size_t get_size() const;

  cell_t clone(void) const;

  void mark_as_removed(bool v);

  bool is_removed(void) const;

  size_t hash() const;

  // inclusion test
  bool operator<=(const cell_t &o) const;

  bool operator==(const cell_t &o) const;

  bool operator<(const cell_t &o) const;

  // Return true if [o, o+size) definitely overlaps with the cell,
  // where o is a constant expression.
  bool overlap(const offset_t &o, uint64_t size) const;

  // Return true if [symb_lb, symb_ub] may overlap with the cell,
  // where symb_lb and symb_ub are not constant expressions.
  template <typename Dom>
  bool symbolic_overlap(const typename Dom::linear_expression_t &symb_lb,
                        const typename Dom::linear_expression_t &symb_ub,
                        const Dom &dom) const {
    if (m_removed)
      return false;

    using linear_expression_t = typename Dom::linear_expression_t;

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

  void write(crab::crab_os &o) const;

  friend crab::crab_os &operator<<(crab::crab_os &o, const cell_t &c) {
    c.write(o);
    return o;
  }
};

namespace cell_set_impl {
template <typename Set>
inline Set set_intersection(const Set &s1, const Set &s2) {
  Set s3;
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::inserter(s3, s3.end()));
  return s3;
}

template <typename Set> inline Set set_union(const Set &s1, const Set &s2) {
  Set s3;
  std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
                 std::inserter(s3, s3.end()));
  return s3;
}

template <typename Set>
inline bool set_inclusion(const Set &s1, const Set &s2) {
  Set s3;
  std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::inserter(s3, s3.end()));
  return s3.empty();
}

template <typename Set>
inline Set set_difference(const Set &s1, const Set &s2) {
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
  template <typename Dom> friend class crab::domains::array_adaptive_domain;

  using cell_set_t = std::set<cell_t>;

  /*
    The keys in the patricia tree are processing in big-endian
    order. This means that the keys are sorted. Sortedeness is
    very important to perform efficiently operations such as
    checking for overlap cells. Since keys are treated as bit
    patterns, negative offsets can be used but they are treated
    as large unsigned numbers.
  */
  using patricia_tree_t = ikos::patricia_tree<offset_t, cell_set_t>;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

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
                                  const patricia_tree_t &t2) const {
    bool res = t1.merge_with(t2, o);
    if (res) {
      CRAB_ERROR("array_adaptive::offset_map should not return bottom");
    }
    return t1;
  }

  class join_op : public binary_op_t {
    // apply is called when two bindings (one each from a
    // different map) have the same key(i.e., offset).
    std::pair<bool, boost::optional<cell_set_t>>
    apply(const offset_t &, const cell_set_t &x, const cell_set_t &y) override {
      return {false, cell_set_impl::set_union(x, y)};
    }
    // if one map does not have a key in the other map we add it.
    bool default_is_absorbing() override { return false; }
  };

  class meet_op : public binary_op_t {
    std::pair<bool, boost::optional<cell_set_t>>
    apply(const offset_t &, const cell_set_t &x, const cell_set_t &y) override {
      return {false, cell_set_impl::set_union(x, y)};
    }
    // if one map does not have a key in the other map we ignore
    // it.
    bool default_is_absorbing() override { return true; }
  };

  class domain_po : public partial_order_t {
    bool leq(const cell_set_t &x, const cell_set_t &y) override {
      return cell_set_impl::set_inclusion(x, y);
    }
    // default value is bottom (i.e., empty map)
    bool default_is_top() override { return false; }
  }; // class domain_po

  // Delete completely the cell
  void erase_cell(const cell_t &c);

  // Pretend the cell is removed by marking it as "removed"
  void remove_cell(const cell_t &c);

  void insert_cell(const cell_t &c);

  cell_t get_cell(const offset_t &o, uint64_t size) const;

  // create a fresh _unamed_ cell
  cell_t mk_cell(const offset_t &o, uint64_t size /*bytes*/);

  offset_map_t(patricia_tree_t &&m);

public:
  offset_map_t();

  offset_map_t(const offset_map_t &o);

  offset_map_t(const offset_map_t &&o);

  offset_map_t &operator=(const offset_map_t &o);

  offset_map_t &operator=(const offset_map_t &&o);

  bool empty() const;

  std::size_t size() const;

  // leq operator
  bool operator<=(const offset_map_t &o) const;

  // set union: if two cells with same offset do not agree on
  // size then they are ignored.
  offset_map_t operator|(const offset_map_t &o) const;

  // set intersection: if two cells with same offset do not agree
  // on size then they are ignored.
  offset_map_t operator&(const offset_map_t &o) const;

  // Completely delete the cell from the offset map
  void erase(const cell_t &c);

  void erase(const std::vector<cell_t> &cells);

  // Pretend the cell is removed so no other cells overlap with it but
  // the cell is not actually deleted from the offset map. We need to
  // know all created cells when smashing occurs.
  void remove(const cell_t &c);

  void remove(const std::vector<cell_t> &cells);

  // cells are sorted by offset
  std::vector<cell_t> get_all_cells() const;

  unsigned get_number_cells() const;

  // Return in out all cells that might overlap with (o, size).
  //
  // It is not marked as const because we insert temporary a cell.
  // However, upon completion this method leaves unmodified the object.
  void get_overlap_cells(const offset_t &o, uint64_t size,
                         std::vector<cell_t> &out);

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

  void clear(void);

  void write(crab::crab_os &o) const;

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

/* A map from an array variable to a vector of pairs of cell and
 * scalar variable.
 *
 * The vector is sorted so that finding/inserting/erasing cells is
 * faster.
 */
template <typename Variable> class cell_varmap {
  using variable_t = Variable;

public:
  using sorted_cell_var_vector = std::vector<std::pair<cell_t, variable_t>>;

private:
  struct less_cell_var {
    bool operator()(const std::pair<cell_t, variable_t> &p,
                    const cell_t &c) const {
      return p.first < c;
    }
  };
  std::unordered_map<variable_t, sorted_cell_var_vector> m_map;

public:
  using iterator =
      typename std::unordered_map<variable_t, sorted_cell_var_vector>::iterator;
  using const_iterator =
      typename std::unordered_map<variable_t,
                                  sorted_cell_var_vector>::const_iterator;
  using cell_var_iterator = typename sorted_cell_var_vector::iterator;

  cell_varmap() {}

  iterator begin() { return m_map.begin(); }
  iterator end() { return m_map.end(); }
  const_iterator begin() const { return m_map.begin(); }
  const_iterator end() const { return m_map.end(); }

  cell_var_iterator begin_cells(const variable_t &array) {
    return m_map[array].begin();
  }
  cell_var_iterator end_cells(const variable_t &array) {
    return m_map[array].end();
  }

  boost::optional<variable_t> find(const variable_t &array,
                                   const cell_t &c) const {
    auto it = m_map.find(array);
    if (it == m_map.end()) {
      return boost::none;
    }

    const sorted_cell_var_vector &cell_vars = it->second;
    less_cell_var cmp;
    auto lb = std::lower_bound(cell_vars.begin(), cell_vars.end(), c, cmp);
    if (lb != cell_vars.end() && !(c < (*lb).first)) {
      return (*lb).second;
    } else {
      return boost::none;
    }
  }

  void erase(const variable_t &array, const cell_t &c) {
    auto it = m_map.find(array);
    if (it == m_map.end()) {
      return;
    }

    sorted_cell_var_vector &cell_vars = it->second;
    less_cell_var cmp;
    auto lb = std::lower_bound(cell_vars.begin(), cell_vars.end(), c, cmp);
    if (lb != cell_vars.end() && !(c < (*lb).first)) {
      cell_vars.erase(lb);
    }
  }

  void erase(const variable_t &array) { m_map.erase(array); }

  void insert(const variable_t &array, const cell_t &c,
              const variable_t &base_v) {
    std::vector<std::pair<cell_t, variable_t>> vs;
    auto it = m_map.insert({array, vs}).first;
    sorted_cell_var_vector &cell_vars = it->second;
    less_cell_var cmp;
    auto lb = std::lower_bound(cell_vars.begin(), cell_vars.end(), c, cmp);
    if (lb != cell_vars.end() && !(c < (*lb).first)) {
      // already exists
      return;
    }
    cell_vars.insert(lb, {c, base_v});
  }

  void insert(const variable_t &array) {
    sorted_cell_var_vector vec;
    auto res = m_map.insert({array, vec});
    if (!res.second) {
      CRAB_ERROR("cell_varmap::insert already exists");
    }
  }

  void write(crab_os &o) const {
    for (auto &kv : m_map) {
      o << kv.first << ":\n";
      for (auto &cv : kv.second) {
        o << "\t" << cv.first << " --> " << cv.second << "\n";
      }
    }
  }
};

/// Model each array cell with a ghost scalar variable.
///
/// TOIMPROVE: we can avoid expensive join and meet operations because
/// each cell can be indentified by a concrete offset so we can use
/// strings to name concrete offsets. This would also solve the risk
/// of breaking modularity of the base domain as mentioned above.
template <typename Domain> class cell_ghost_man {
  using variable_t = typename Domain::variable_t;
  using variable_vector_t = typename Domain::variable_vector_t;
  using varname_t = typename variable_t::varname_t;
  using type_t = typename variable_t::type_t;
  using offset_t = array_adaptive_impl::offset_t;
  using cell_t = array_adaptive_impl::cell_t;
  using ghost_map_t = array_adaptive_impl::cell_varmap<variable_t>;
  using this_type = cell_ghost_man<Domain>;

  ghost_map_t m_map;

  static variable_type_kind get_array_element_type(type_t array_type) {
    if (array_type.is_bool_array()) {
      return BOOL_TYPE;
    } else if (array_type.is_integer_array()) {
      return INT_TYPE;
    } else {
      assert(array_type.is_real_array());
      return REAL_TYPE;
    }
  }

  static std::string mk_scalar_name(const varname_t &a, const offset_t &o,
                                    uint64_t size) {
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

  variable_t make_ghost(const variable_t &a, const cell_t &c) {
    assert(!c.is_null());
    // assign a scalar variable to the cell
    auto &vfac = const_cast<varname_t *>(&(a.name()))->get_var_factory();
    std::string vname = mk_scalar_name(a.name(), c.get_offset(), c.get_size());
    variable_type_kind vtype_kind = get_array_element_type(a.get_type());
    variable_t scalar_var(
        vfac.get(vname), vtype_kind,
        (vtype_kind == BOOL_TYPE
             ? 1
             : (vtype_kind == INT_TYPE ? 8 * c.get_size() : 0)));
    return scalar_var;
  }

  variable_t insert_ghost(const variable_t &a, const cell_t &c) {
    variable_t ghost_v = make_ghost(a, c);
    m_map.insert(a, c, ghost_v);
    return ghost_v;
  }

  cell_ghost_man(ghost_map_t &&m) : m_map(std::move(m)) {}

public:
  cell_ghost_man() {}

  variable_t get_or_insert_ghost(const variable_t &a, const cell_t &c) {
    if (boost::optional<variable_t> scalar_v = m_map.find(a, c)) {
      return *scalar_v;
    } else {
      return insert_ghost(a, c);
    }
  }

  boost::optional<variable_t> get_ghost(const variable_t &a,
                                        const cell_t &c) const {
    return m_map.find(a, c);
  }

  void erase_all(const variable_t &a) { m_map.erase(a); }

  void erase(const variable_t &a, const cell_t &c) { m_map.erase(a, c); }

  this_type join(const this_type &right_man, Domain &left_val,
                 Domain &right_val) const {

    std::vector<variable_t> old_ghost_vars_left, old_ghost_vars_right,
        new_ghost_vars;
    ghost_map_t out_map;
    const ghost_map_t &left_map = m_map;
    const ghost_map_t &right_map = right_man.m_map;

    for (auto const &kv : left_map) {
      for (auto const &cv : kv.second) {
        const variable_t &array_var = kv.first;
        const cell_t &c = cv.first;
        const variable_t &v1 = cv.second;
        auto &vfac = const_cast<varname_t *>(&(v1.name()))->get_var_factory();
        if (boost::optional<variable_t> v2 = right_map.find(array_var, c)) {
          if (v1 != *v2) {
            assert(v1.name().str() == (*v2).name().str());
            assert(v1.get_type() == (*v2).get_type());
            variable_t outv(vfac.get(v1.name().str()), v1.get_type());
            old_ghost_vars_left.push_back(v1);
            old_ghost_vars_right.push_back(*v2);
            new_ghost_vars.push_back(outv);
            out_map.insert(array_var, c, outv);
          } else {
            out_map.insert(array_var, c, v1);
          }
        }
      }
    }
    /// project might be necessary to avoid keeping variables that
    /// exist in the base domain but they doen't exist on m_map.
    ///
    /// If such a variable exists only on either left_val or right_val
    /// the join removes it.  However, if we have the same variable in
    /// both left_val and right_val then the join will preserve it.

    // left_val.project(old_ghost_vars_left);
    left_val.rename(old_ghost_vars_left, new_ghost_vars);
    // right_val.project(old_ghost_vars_right);
    right_val.rename(old_ghost_vars_right, new_ghost_vars);

    return this_type(std::move(out_map));
  }

  this_type meet(const this_type &right_man, Domain &left_val,
                 Domain &right_val) const {

    std::vector<variable_t> old_ghost_vars_left, old_ghost_vars_right,
        new_ghost_vars;
    ghost_map_t out_map;
    const ghost_map_t &left_map = m_map;
    const ghost_map_t &right_map = right_man.m_map;

    // Add all mappings from the left operand
    for (auto const &kv : left_map) {
      for (auto const &cv : kv.second) {
        const variable_t &array_var = kv.first;
        const cell_t &c = cv.first;
        const variable_t &v1 = cv.second;
        auto &vfac = const_cast<varname_t *>(&(v1.name()))->get_var_factory();
        if (boost::optional<variable_t> v2 = right_map.find(array_var, c)) {
          if (v1 != *v2) {
            // same key but different ghost -> create a fresh common ghost
            assert(v1.name().str() == (*v2).name().str());
            assert(v1.get_type() == (*v2).get_type());
            variable_t outv(vfac.get(v1.name().str()), v1.get_type());
            old_ghost_vars_left.push_back(v1);
            old_ghost_vars_right.push_back(*v2);
            new_ghost_vars.push_back(outv);
            out_map.insert(array_var, c, outv);
            continue;
          }
        }
        out_map.insert(array_var, c, v1);
      }
    }

    // Add the rest of mappings from the right operand
    for (auto const &kv : right_map) {
      for (auto const &cv : kv.second) {
        const variable_t &array_var = kv.first;
        const cell_t &c = cv.first;
        const variable_t &v = cv.second;
        if (!left_map.find(array_var, c)) {
          out_map.insert(array_var, c, v);
        }
      }
    }
    /// See comments in join about project.
    // left_val.project(old_ghost_vars_left);
    left_val.rename(old_ghost_vars_left, new_ghost_vars);
    // right_val.project(old_ghost_vars_right);
    right_val.rename(old_ghost_vars_right, new_ghost_vars);

    return this_type(std::move(out_map));
  }

  void rename(const variable_t &old_array, const variable_t &new_array,
              std::function<cell_t(const offset_t &, uint64_t)> mk_unnamed_cell,
              variable_vector_t &old_ghosts, variable_vector_t &new_ghosts) {

    // We insert here an empty vector so that no resizing can happen
    // while we iterate over m_map.
    m_map.insert(new_array);
    for (auto it = m_map.begin_cells(old_array),
              et = m_map.end_cells(old_array);
         it != et; ++it) {
      const cell_t &old_c = (*it).first;
      const variable_t &old_ghost = (*it).second;
      /// Modify m_map but it doesn't invalidate iterators
      /// because new_array is already in m_map.
      cell_t new_c = mk_unnamed_cell(old_c.get_offset(), old_c.get_size());
      variable_t new_ghost = insert_ghost(new_array, new_c);
      old_ghosts.push_back(old_ghost);
      new_ghosts.push_back(new_ghost);
    }
  }

  void write(crab_os &o) const {
    for (auto const &kv : m_map) {
      for (auto const &cv : kv.second) {
        crab::outs() << kv.first << "#" << cv.first << " -> " << cv.second
                     << "\n";
      }
    }
  }
};

// Constant value extended with bottom and top.
class constant_value {
  using bound_t = ikos::bound<ikos::z_number>;

  bool m_is_bottom;
  bound_t m_val;

  void set_to_top();

  void set_to_bot();

  constant_value(bool is_bottom);

public:
  constant_value();

  constant_value(int64_t sz);

  constant_value(const ikos::interval<ikos::z_number> &sz);

  bool is_top() const;

  bool is_bottom() const;

  bool is_zero() const;

  bool is_negative() const;

  bound_t val() const;

  boost::optional<uint64_t> get_uint64_val() const;

  static constant_value bottom();

  static constant_value top();

  bool operator==(const constant_value &o) const;

  bool operator<=(const constant_value &o) const;

  void operator|=(const constant_value &o);

  constant_value operator|(const constant_value &o) const;

  constant_value operator&(const constant_value &o) const;

  constant_value operator||(const constant_value &o) const;

  constant_value operator&&(const constant_value &o) const;

  void write(crab_os &o) const;
};

} // end namespace array_adaptive_impl

template <typename NumDomain>
class array_adaptive_domain final
    : public abstract_domain_api<array_adaptive_domain<NumDomain>> {

public:
  using number_t = typename NumDomain::number_t;
  using varname_t = typename NumDomain::varname_t;

private:
  using array_adaptive_domain_t = array_adaptive_domain<NumDomain>;
  using abstract_domain_t = abstract_domain_api<array_adaptive_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using base_domain_t = array_smashing<NumDomain>;

private:
  using type_t = typename variable_t::type_t;
  using offset_t = array_adaptive_impl::offset_t;
  using offset_map_t = array_adaptive_impl::offset_map_t;
  using cell_t = array_adaptive_impl::cell_t;
  using constant_value = array_adaptive_impl::constant_value;
  using cell_ghost_man_t = array_adaptive_impl::cell_ghost_man<base_domain_t>;

  class array_state {
    // whether the array has been smashed
    bool m_is_smashed;
    // element size of the array if smashed
    constant_value m_element_sz;
    // precise array contents if no smashed
    offset_map_t m_offset_map;

    static bool consistent_offset(const offset_t &o, size_t elem_size) {
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
    void smash_array(const variable_t &a, const constant_value &elem_sz,
                     cell_ghost_man_t &cm, base_domain_t &base_dom) {

      // we only smash the array if elem_sz is a constant value and
      // all array elements are consistent wrt elem_sz. Leave the
      // array without smashing is always sound.

      std::vector<cell_t> cells = get_offset_map().get_all_cells();
      if (cells.empty()) {
        return;
      }

      if (cells.size() >
          crab_domain_params_man::get().array_adaptive_max_smashable_cells()) {
        // Smashing is expensive because it will go over all cells
        // performing one join per weak update even if the smashed
        // array is already unconstrained. We don't smash if the
        // number of cells to be smashed is too large.
        CRAB_WARN(
            "array adaptive did not smash array because its size is "
            "greater than ",
            crab_domain_params_man::get().array_adaptive_max_smashable_cells());
        return;
      }

      if (!crab_domain_params_man::get()
               .array_adaptive_smash_at_nonzero_offset() &&
          !cells[0].get_offset().is_zero()) {
        return;
      }

      if (boost::optional<uint64_t> sz_opt = elem_sz.get_uint64_val()) {
        bool can_be_smashed = true;
        for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
          const cell_t &c = cells[k];
          if (!consistent_offset(c.get_offset(), *sz_opt)) {
            CRAB_LOG(
                "array-adaptive",
                CRAB_WARN("cannot smashing because of inconsistent offsets"););
            can_be_smashed = false;
            break;
          }
          const bool is_strong_update = (k == 0);
          linear_expression_t idx(number_t(c.get_offset().index()));
          if (boost::optional<variable_t> c_scalar_var = cm.get_ghost(a, c)) {
            base_dom.array_store(a, c.get_size(), idx, *c_scalar_var,
                                 is_strong_update);
          } else {
            CRAB_LOG(
                "array-adaptive",
                CRAB_WARN(
                    "cannot smashing because of scalar variable for array ", a,
                    " and cell ", c, " not found in ");
                cm.write(crab::outs()); crab::outs() << "\n";);
            can_be_smashed = false;
            break;
          }
        }

        if (can_be_smashed) {
          for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
            const cell_t &c = cells[k];
            if (boost::optional<variable_t> c_scalar_var = cm.get_ghost(a, c)) {
              // remove the synthethic cell from the base domain and
              // from the offset map
              base_dom -= *c_scalar_var;
            }
            get_offset_map().erase(c);
            cm.erase(a, c);
          }
          m_is_smashed = true;
          m_element_sz = elem_sz;
        } else {
          base_dom -= a;
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

    array_state(bool &&is_smashed, constant_value &&sz, offset_map_t &&om)
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
    array_state join(const variable_t &v, const array_state &o,
                     cell_ghost_man_t &cm_left, base_domain_t &dom_left,
                     cell_ghost_man_t &cm_right,
                     base_domain_t &dom_right) const {
      if (m_is_smashed && !o.m_is_smashed) {
        array_state right(o);
        right.smash_array(v, get_element_sz(), cm_right, dom_right);
        return array_state(m_is_smashed | right.m_is_smashed,
                           m_element_sz | right.m_element_sz,
                           m_offset_map | right.m_offset_map);
      } else if (!m_is_smashed && o.m_is_smashed) {
        array_state left(*this);
        left.smash_array(v, o.get_element_sz(), cm_left, dom_left);
        return array_state(left.m_is_smashed | o.m_is_smashed,
                           left.m_element_sz | o.m_element_sz,
                           left.m_offset_map | o.m_offset_map);
      } else {
        return array_state(m_is_smashed | o.m_is_smashed,
                           m_element_sz | o.m_element_sz,
                           m_offset_map | o.m_offset_map);
      }
    }

    array_state meet(const variable_t &v, const array_state &o,
                     cell_ghost_man_t &cm_left, base_domain_t &dom_left,
                     cell_ghost_man_t &cm_right,
                     base_domain_t &dom_right) const {
      if (m_is_smashed && !o.m_is_smashed) {
        array_state right(o);
        right.smash_array(v, get_element_sz(), cm_right, dom_right);
        return array_state(m_is_smashed & right.m_is_smashed,
                           m_element_sz & right.m_element_sz,
                           m_offset_map & right.m_offset_map);
      } else if (!m_is_smashed && o.m_is_smashed) {
        array_state left(*this);
        left.smash_array(v, o.get_element_sz(), cm_left, dom_left);
        return array_state(left.m_is_smashed & o.m_is_smashed,
                           left.m_element_sz & o.m_element_sz,
                           left.m_offset_map & o.m_offset_map);
      } else {
        return array_state(m_is_smashed & o.m_is_smashed,
                           m_element_sz & o.m_element_sz,
                           m_offset_map & o.m_offset_map);
      }
    }

    bool operator==(const array_state &o) const {
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

    constant_value &get_element_sz() { return m_element_sz; }

    const constant_value &get_element_sz() const { return m_element_sz; }

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
                            crab_domain_params_man::get()
                                .array_adaptive_smash_at_nonzero_offset());
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
    using patricia_tree_t = ikos::patricia_tree<variable_t, array_state>;
    using binary_op_t = typename patricia_tree_t::binary_op_t;

  public:
    using iterator = typename patricia_tree_t::iterator;

  private:
    patricia_tree_t m_tree;

    class join_op : public binary_op_t {
      cell_ghost_man_t &m_cm_left;
      base_domain_t &m_dom_left;
      cell_ghost_man_t &m_cm_right;
      base_domain_t &m_dom_right;

    public:
      join_op(cell_ghost_man_t &cm_left, base_domain_t &dom_left,
              cell_ghost_man_t &cm_right, base_domain_t &dom_right)
          : m_cm_left(cm_left), m_dom_left(dom_left), m_cm_right(cm_right),
            m_dom_right(dom_right) {}

      std::pair<bool, boost::optional<array_state>>
      apply(const variable_t &k, const array_state &x,
            const array_state &y) override {
        array_state z =
            x.join(k, y, m_cm_left, m_dom_left, m_cm_right, m_dom_right);
        return {false, boost::optional<array_state>(z)};
      }

      bool default_is_absorbing() override { return true; }
    }; // class join_op

    class meet_op : public binary_op_t {
      cell_ghost_man_t &m_cm_left;
      base_domain_t &m_dom_left;
      cell_ghost_man_t &m_cm_right;
      base_domain_t &m_dom_right;

    public:
      meet_op(cell_ghost_man_t &cm_left, base_domain_t &dom_left,
              cell_ghost_man_t &cm_right, base_domain_t &dom_right)
          : m_cm_left(cm_left), m_dom_left(dom_left), m_cm_right(cm_right),
            m_dom_right(dom_right) {}

      std::pair<bool, boost::optional<array_state>>
      apply(const variable_t &k, const array_state &x,
            const array_state &y) override {
        array_state z =
            x.meet(k, y, m_cm_left, m_dom_left, m_cm_right, m_dom_right);
        return {false, boost::optional<array_state>(z)};
      }
      bool default_is_absorbing() override { return false; }
    }; // class meet_op

    patricia_tree_t apply_operation(binary_op_t &o, patricia_tree_t t1,
                                    const patricia_tree_t &t2) const {
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
    array_state_map_t join(const array_state_map_t &o,
                           cell_ghost_man_t &cm_left, base_domain_t &dom_left,
                           cell_ghost_man_t &cm_right,
                           base_domain_t &dom_right) const {
      join_op op(cm_left, dom_left, cm_right, dom_right);
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree);
      return array_state_map_t(std::move(res));
    }

    // Meet
    array_state_map_t meet(const array_state_map_t &o,
                           cell_ghost_man_t &cm_left, base_domain_t &dom_left,
                           cell_ghost_man_t &cm_right,
                           base_domain_t &dom_right) const {
      meet_op op(cm_left, dom_left, cm_right, dom_right);
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

    // Assume that from does not have duplicates.
    void rename(const std::vector<variable_t> &from,
                const std::vector<variable_t> &to) {
      if (from.size() != to.size()) {
        CRAB_ERROR("array_adaptive::array_state_t::rename received input "
                   "vectors of different sizes");
      }

      if (m_tree.size() == 0) {
        return;
      }

      for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
        variable_t k = from[i];
        variable_t new_k = to[i];
        if (k == new_k) { // nothing to rename
          continue;
        }

        if (::crab::CrabSanityCheckFlag) {
          if (m_tree.lookup(new_k)) {
            CRAB_ERROR("array_adaptive::array_state_t:rename assumes that  ",
                       new_k, " does not exist");
          }
        }

        if (boost::optional<array_state> k_val_opt = m_tree.lookup(k)) {
          if (!(*k_val_opt).is_top()) {
            m_tree.insert(new_k, *k_val_opt);
          }
          m_tree.remove(k);
        }
      }
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
        if (it != m_tree.end()) {
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

  // -- scalar domain containing scalar variables, ghost scalar
  //    variables from cells and ghost summarized smashed variables.
  base_domain_t m_base_dom;
  // -- map an array variable to its synthetic cells
  array_state_map_t m_array_map;
  // -- manage ghost variables to model concrete arrays
  cell_ghost_man_t m_cell_ghost_man;

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

  // Return a named cell, i.e., a cell with a scalar variable
  // associated to it.
  std::pair<cell_t, variable_t> mk_named_cell(const variable_t &a,
                                              const offset_t &o,
                                              uint64_t sz /*bytes*/,
                                              offset_map_t &om) {
    // create first an unnamed cell
    cell_t c = om.mk_cell(o, sz);
    variable_t scalar_v = m_cell_ghost_man.get_or_insert_ghost(a, c);
    return {c, scalar_v};
  }

  using variable_opt_t = boost::optional<variable_t>;
  variable_opt_t get_scalar(const variable_t &array_v, const cell_t &c) {
    if (!array_v.get_type().is_array()) {
      CRAB_ERROR("array_adaptive::get_scalar only if array variable");
    }
    return m_cell_ghost_man.get_ghost(array_v, c);
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
    if (!v.get_type().is_array()) {
      CRAB_ERROR("cannot call forget_array on a non-array variable");
    }

    std::vector<variable_t> base_vars;
    const array_state &as = lookup_array_state(v);
    if (!as.is_smashed()) {
      /// We extract all the synthetic cells from the array and forget
      /// them from the underlying abstract domain.
      const offset_map_t &om = as.get_offset_map();
      std::vector<cell_t> cells = om.get_all_cells();
      for (auto &c : cells) {
        if (variable_opt_t v_opt = get_scalar(v, c)) {
          base_vars.push_back(*v_opt);
        }
      }
    } else {
      base_vars.push_back(v);
    }

    m_array_map -= v;
    m_cell_ghost_man.erase_all(v);
    m_base_dom.forget(base_vars);
  }

  interval_t to_interval(const linear_expression_t &expr, base_domain_t inv) {
    interval_t r(expr.constant());
    for (auto kv : expr) {
      interval_t c(kv.first);
      r += c * inv[kv.second];
    }
    return r;
  }

  interval_t to_interval(const linear_expression_t &expr) {
    return to_interval(expr, m_base_dom);
  }

  void kill_cells(const variable_t &a, const std::vector<cell_t> &cells,
                  offset_map_t &offset_map) {

    assert(a.get_type().is_array());

    if (!cells.empty()) {
      // Forget the scalars from the numerical domain
      for (unsigned i = 0, e = cells.size(); i < e; ++i) {
        const cell_t &c = cells[i];
        if (variable_opt_t c_scalar_opt = get_scalar(a, c)) {
          m_base_dom -= *c_scalar_opt;
        }
      }
      if (!crab_domain_params_man::get().array_adaptive_is_smashable()) {
        // Delete completely the cells. If needed again they they will
        // be re-created.
        for (unsigned i = 0, e = cells.size(); i < e; ++i) {
          const cell_t &c = cells[i];
          offset_map.erase(c);
          m_cell_ghost_man.erase(a, c);
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
  void do_assign(const variable_t &lhs, const variable_t &rhs) {
    if (lhs.get_type() != rhs.get_type()) {
      CRAB_ERROR("array_adaptive assignment ", lhs, ":=", rhs,
                 " with different types");
    }
    auto lhs_ty = lhs.get_type();
    if (lhs_ty.is_bool()) {
      m_base_dom.assign_bool_var(lhs, rhs, false);
    } else if (lhs_ty.is_integer() || lhs_ty.is_real()) {
      m_base_dom.assign(lhs, rhs);
    } else {
      CRAB_ERROR(
          "array_adaptive assignment with unexpected array element type");
    }
  }

  // helper to assign an array store's value
  void do_assign(const variable_t &lhs, const linear_expression_t &v) {
    auto lhs_ty = lhs.get_type();
    if (lhs_ty.is_bool()) {
      if (v.is_constant()) {
        if (v.constant() >= number_t(1)) {
          m_base_dom.assign_bool_cst(lhs, linear_constraint_t::get_true());
        } else {
          m_base_dom.assign_bool_cst(lhs, linear_constraint_t::get_false());
        }
      } else if (auto var = v.get_variable()) {
        m_base_dom.assign_bool_var(lhs, (*var), false);
      }
    } else if (lhs_ty.is_integer() || lhs_ty.is_real()) {
      m_base_dom.assign(lhs, v);
    } else {
      CRAB_ERROR(
          "array_adaptive assignment with unexpected array element type");
    }
  }

  // Helper that assign backward rhs to lhs by switching to the
  // version with the right type.
  void do_backward_assign(const variable_t &lhs, const variable_t &rhs,
                          const base_domain_t &dom) {
    if (lhs.get_type() != rhs.get_type()) {
      CRAB_ERROR("array_adaptive backward assignment with different types");
    }
    auto lhs_ty = lhs.get_type();
    if (lhs_ty.is_bool()) {
      m_base_dom.backward_assign_bool_var(lhs, rhs, false, dom);
    } else if (lhs_ty.is_integer() || lhs_ty.is_real()) {
      m_base_dom.backward_assign(lhs, rhs, dom);
    } else {
      CRAB_ERROR("array_adaptive backward_assignment with unexpected array "
                 "element type");
    }
  }

  // helper to assign backward a cell into a variable
  void do_backward_assign(const variable_t &lhs, const variable_t &a,
                          const cell_t &rhs_c, const base_domain_t &dom) {
    if (!a.get_type().is_array()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }
    variable_opt_t rhs_v_opt = get_scalar(a, rhs_c);
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
  void do_backward_assign(const variable_t &a, const cell_t &lhs_c,
                          const linear_expression_t &v,
                          const base_domain_t &dom) {
    if (!a.get_type().is_array()) {
      CRAB_ERROR("array_adaptive assignment 1st argument must be array type");
    }
    variable_opt_t lhs_v_opt = get_scalar(a, lhs_c);
    if (!lhs_v_opt) {
      CRAB_LOG(
          "array-adaptive",
          CRAB_WARN(
              "array_adaptive cell without scalar in do_backward_assign"););
      return;
    }
    variable_t lhs = *lhs_v_opt;
    auto lhs_ty = lhs.get_type();
    if (lhs_ty.is_bool()) {
      if (v.is_constant()) {
        if (v.constant() >= number_t(1)) {
          m_base_dom.backward_assign_bool_cst(
              lhs, linear_constraint_t::get_true(), dom);
        } else {
          m_base_dom.backward_assign_bool_cst(
              lhs, linear_constraint_t::get_false(), dom);
        }
      } else if (auto var = v.get_variable()) {
        m_base_dom.backward_assign_bool_var(lhs, (*var), false, dom);
      }
    } else if (lhs_ty.is_integer() || lhs_ty.is_real()) {
      m_base_dom.backward_assign(lhs, v, dom);
    } else {
      CRAB_ERROR("array_adaptive backward assignment with unexpected array "
                 "element type");
    }
  }

  // The internal representation contains summarized variables of
  // array type and add them as dimensions in the underlying numerical
  // domain. They shouldn't be exposed outside via linear constraints.
  //
  // XXX: we should also probably filter out scalar variables
  // originated from cells.
  linear_constraint_system_t
  filter_nonscalar_vars(linear_constraint_system_t &&csts) const {
    linear_constraint_system_t res;
    for (auto const &cst : csts) {
      if (std::all_of(
              cst.expression().variables_begin(),
              cst.expression().variables_end(), [](const variable_t &v) {
                return v.get_type().is_integer() || v.get_type().is_bool();
              })) {
        res += cst;
      }
    }
    return res;
  }

  uint64_t check_and_get_elem_size(const linear_expression_t &elem_size) {
    interval_t i_elem_size = to_interval(elem_size);
    if (boost::optional<number_t> n_bytes = i_elem_size.singleton()) {
      if (static_cast<int64_t>(*n_bytes) > 0) {
        return (uint64_t) static_cast<int64_t>(*n_bytes);
      }
    }
    CRAB_ERROR("array adaptive domain expects constant array element sizes ",
               "between 1 and ", std::numeric_limits<uint64_t>::max(),
               ". Found ", elem_size);
  }

  struct renamed_meet_state {
    base_domain_t left_dom;
    base_domain_t right_dom;
    array_state_map_t array_map;
    cell_ghost_man_t cell_ghost_man;
    renamed_meet_state(base_domain_t &&left, base_domain_t &&right,
                       array_state_map_t &&_array_map,
                       cell_ghost_man_t &&_cell_ghost_man)
        : left_dom(std::move(left)), right_dom(std::move(right)),
          array_map(std::move(_array_map)),
          cell_ghost_man(std::move(_cell_ghost_man)) {}
  };

  // It performs all the renaming needed for meet.  The only operation
  // is not done is the actual meet of the based domains.
  renamed_meet_state
  rename_with_meet_semantics(const array_adaptive_domain_t &other) const {
    // Trivial cases are covered elsewhere(i.e., this and other are
    // bottom and/or top)
    base_domain_t left_dom(m_base_dom);
    cell_ghost_man_t left_cell_ghost_man(m_cell_ghost_man);

    base_domain_t right_dom(other.m_base_dom);
    cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

    // Must be done before the renaming.
    auto out_array_map =
        m_array_map.meet(other.m_array_map, left_cell_ghost_man, left_dom,
                         right_cell_ghost_man, right_dom);

    cell_ghost_man_t out_cell_ghost_man =
        left_cell_ghost_man.meet(right_cell_ghost_man, left_dom, right_dom);

    renamed_meet_state res(std::move(left_dom), std::move(right_dom),
                           std::move(out_array_map),
                           std::move(out_cell_ghost_man));
    return res;
  }
  
  array_adaptive_domain(base_domain_t &&inv, array_state_map_t &&amap,
                        cell_ghost_man_t &&cgman)
      : m_base_dom(std::move(inv)), m_array_map(std::move(amap)),
        m_cell_ghost_man(std::move(cgman)) {}

public:
  array_adaptive_domain(bool is_bottom = false) {
    if (is_bottom) {
      m_base_dom.set_to_bottom();
    } else {
      m_base_dom.set_to_top();
    }
  }

  array_adaptive_domain make_top() const override {
    array_adaptive_domain out(false);
    return out;
  }

  array_adaptive_domain make_bottom() const override {
    array_adaptive_domain out(true);
    return out;
  }

  void set_to_top() override {
    array_adaptive_domain abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    array_adaptive_domain abs(true);
    std::swap(*this, abs);
  }

  array_adaptive_domain(const array_adaptive_domain_t &other)
      : m_base_dom(other.m_base_dom), m_array_map(other.m_array_map),
        m_cell_ghost_man(other.m_cell_ghost_man) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  array_adaptive_domain(const array_adaptive_domain_t &&other)
      : m_base_dom(std::move(other.m_base_dom)),
        m_array_map(std::move(other.m_array_map)),
        m_cell_ghost_man(std::move(other.m_cell_ghost_man)) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  array_adaptive_domain_t &operator=(const array_adaptive_domain_t &other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      m_base_dom = other.m_base_dom;
      m_array_map = other.m_array_map;
      m_cell_ghost_man = other.m_cell_ghost_man;
    }
    return *this;
  }

  array_adaptive_domain_t &operator=(const array_adaptive_domain_t &&other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      m_base_dom = std::move(other.m_base_dom);
      m_array_map = std::move(other.m_array_map);
      m_cell_ghost_man = std::move(other.m_cell_ghost_man);
    }
    return *this;
  }

  bool is_bottom() const override { return (m_base_dom.is_bottom()); }

  bool is_top() const override { return (m_base_dom.is_top()); }

  bool operator<=(const array_adaptive_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom()) {
      return true;
    } else if (other.is_top()) {
      return true;
    } else {

      CRAB_LOG("array-adaptive", array_adaptive_domain_t left(*this);
               array_adaptive_domain_t right(other);
               crab::outs() << "Check if " << left << " <= " << right << "\n");

      base_domain_t left_dom(m_base_dom);
      base_domain_t right_dom(other.m_base_dom);

      m_cell_ghost_man.join(other.m_cell_ghost_man, left_dom, right_dom);

      // We need to be careful if one array state is smashed and the
      // other is not.
      bool res = (left_dom <= right_dom);
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  bool operator==(array_adaptive_domain_t other) {
    return (m_base_dom <= other.m_base_dom && other.m_base_dom <= m_base_dom);
  }

  void operator|=(const array_adaptive_domain_t &other) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    CRAB_LOG("array-adaptive",
             crab::outs() << "Join " << *this << " and " << other << "\n";);

    if (other.is_bottom() || is_top()) {
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
      return;
    } else if (is_bottom() || other.is_top()) {
      *this = other;
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
    } else {

      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

      // this must be done before the renaming
      m_array_map = std::move(
          m_array_map.join(other.m_array_map, m_cell_ghost_man, m_base_dom,
                           right_cell_ghost_man, right_dom));

      cell_ghost_man_t out_cell_ghost_man =
          m_cell_ghost_man.join(right_cell_ghost_man, m_base_dom, right_dom);

      m_base_dom |= right_dom;
      std::swap(m_cell_ghost_man, out_cell_ghost_man);
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
    }
  }

  array_adaptive_domain_t
  operator|(const array_adaptive_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    if (other.is_bottom() || is_top()) {
      return *this;
    } else if (is_bottom() || other.is_top()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive",
               crab::outs() << "Join " << *this << " and " << other << "\n";);

      base_domain_t left_dom(m_base_dom);
      cell_ghost_man_t left_cell_ghost_man(m_cell_ghost_man);

      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

      // Must be done before the renaming.
      auto out_array_map = std::move(
          m_array_map.join(other.m_array_map, left_cell_ghost_man, left_dom,
                           right_cell_ghost_man, right_dom));

      cell_ghost_man_t out_cell_ghost_man =
          left_cell_ghost_man.join(right_cell_ghost_man, left_dom, right_dom);

      array_adaptive_domain_t res(left_dom | right_dom,
                                  std::move(out_array_map),
                                  std::move(out_cell_ghost_man));

      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t
  operator&(const array_adaptive_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive",
               crab::outs() << "Meet " << *this << " and " << other << "\n";);
      auto s = rename_with_meet_semantics(other);
      array_adaptive_domain_t res(s.left_dom & s.right_dom,
                                  std::move(s.array_map),
                                  std::move(s.cell_ghost_man));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  void operator&=(const array_adaptive_domain_t &other) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (is_top() || other.is_bottom()) {
      *this = other;
    } else {
      CRAB_LOG("array-adaptive",
	       crab::outs() << "Meet " << *this << " and " << other << "\n";);


      base_domain_t &left_dom = m_base_dom;
      cell_ghost_man_t &left_cell_ghost_man = m_cell_ghost_man;
      
      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);
      
      // Must be done before the renaming.
      auto out_array_map =
        m_array_map.meet(other.m_array_map, left_cell_ghost_man, left_dom,
                         right_cell_ghost_man, right_dom);

      cell_ghost_man_t out_cell_ghost_man =
        left_cell_ghost_man.meet(right_cell_ghost_man, left_dom, right_dom);

      left_dom &= right_dom;
      m_array_map = std::move(out_array_map);
      m_cell_ghost_man = std::move(out_cell_ghost_man);
					   
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
    }
  }
  
  array_adaptive_domain_t
  operator||(const array_adaptive_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    if (other.is_bottom()) {
      return *this;
    } else if (is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Widening " << *this << " and "
                                              << other << "\n";);

      base_domain_t left_dom(m_base_dom);
      cell_ghost_man_t left_cell_ghost_man(m_cell_ghost_man);

      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

      // Must be done before the renaming.
      auto out_array_map =
          m_array_map.join(other.m_array_map, left_cell_ghost_man, left_dom,
                           right_cell_ghost_man, right_dom);

      cell_ghost_man_t out_cell_ghost_man =
          left_cell_ghost_man.join(right_cell_ghost_man, left_dom, right_dom);

      array_adaptive_domain_t res(left_dom || right_dom,
                                  std::move(out_array_map),
                                  std::move(out_cell_ghost_man));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t
  widening_thresholds(const array_adaptive_domain_t &other,
                      const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    if (other.is_bottom()) {
      return *this;
    } else if (is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Widening " << *this << " and "
                                              << other << "\n";);

      base_domain_t left_dom(m_base_dom);
      cell_ghost_man_t left_cell_ghost_man(m_cell_ghost_man);

      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

      // Must be done before the renaming.
      auto out_array_map =
          m_array_map.join(other.m_array_map, left_cell_ghost_man, left_dom,
                           right_cell_ghost_man, right_dom);

      cell_ghost_man_t out_cell_ghost_man =
          left_cell_ghost_man.join(right_cell_ghost_man, left_dom, right_dom);

      array_adaptive_domain_t res(left_dom.widening_thresholds(right_dom, ts),
                                  std::move(out_array_map),
                                  std::move(out_cell_ghost_man));
      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  array_adaptive_domain_t
  operator&&(const array_adaptive_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    if (is_bottom()) {
      return *this;
    } else if (other.is_bottom()) {
      return other;
    } else {
      CRAB_LOG("array-adaptive", crab::outs() << "Narrowing " << *this
                                              << " and " << other << "\n";);

      base_domain_t left_dom(m_base_dom);
      cell_ghost_man_t left_cell_ghost_man(m_cell_ghost_man);

      base_domain_t right_dom(other.m_base_dom);
      cell_ghost_man_t right_cell_ghost_man(other.m_cell_ghost_man);

      // Must be done before the renaming.
      auto out_array_map =
          m_array_map.join(other.m_array_map, left_cell_ghost_man, left_dom,
                           right_cell_ghost_man, right_dom);

      cell_ghost_man_t out_cell_ghost_man =
          left_cell_ghost_man.meet(right_cell_ghost_man, left_dom, right_dom);

      array_adaptive_domain_t res(left_dom && right_dom,
                                  std::move(out_array_map),
                                  std::move(out_cell_ghost_man));

      CRAB_LOG("array-adaptive", crab::outs() << "Res=" << res << "\n";);
      return res;
    }
  }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    variable_vector_t scalar_variables;
    scalar_variables.reserve(variables.size());
    for (variable_t v : variables) {
      if (v.get_type().is_array()) {
        CRAB_LOG("array-adaptive",
                 crab::outs() << "Forget array variable " << v << "\n";);
        forget_array(v);
      } else {
        CRAB_LOG("array-adaptive",
                 crab::outs() << "Forget scalar variable " << v << "\n";);
        scalar_variables.push_back(v);
      }
    }
    m_base_dom.forget(scalar_variables);
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    // if we must keep array variable v then we need to keep all its
    // synthetic cells in m_base_dom.

    variable_vector_t keep_vars;
    std::set<variable_t> keep_arrays;
    for (variable_t v : variables) {
      if (v.get_type().is_array()) {
        keep_arrays.insert(v);
      } else {
        keep_vars.push_back(v);
      }
    }

    std::vector<variable_t> array_variables = get_array_variables();
    for (variable_t v : array_variables) {
      if (keep_arrays.count(v) > 0) {
        // keep all cells of v
        const array_state &as = lookup_array_state(v);
        if (!as.is_smashed()) {
          const offset_map_t &om = as.get_offset_map();
          std::vector<cell_t> cells = om.get_all_cells();
          for (auto &c : cells) {
            if (variable_opt_t v_opt = get_scalar(v, c)) {
              keep_vars.push_back(*v_opt);
            }
          }
        } else {
          keep_vars.push_back(v);
        }
      } else {
        m_array_map -= v;
        m_cell_ghost_man.erase_all(v);
      }
    }

    // Finally we project
    m_base_dom.project(keep_vars);
  }

  void normalize() override {
    CRAB_WARN("array adaptive normalize not implemented");
  }

  void minimize() override { m_base_dom.minimize(); }

  virtual interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      if (!v.get_type().is_array()) {
        return m_base_dom.at(v);
      } else {
        return interval_t::top();
      }
    }
  }

  virtual interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      if (!v.get_type().is_array()) {
        return m_base_dom.at(v);
      } else {
        return interval_t::top();
      }
    }
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    if ((std::all_of(inputs.begin(), inputs.end(),
                     [](const variable_or_constant_t &v) {
                       return v.get_type().is_integer() ||
                              v.get_type().is_bool();
                     })) &&
        (std::all_of(outputs.begin(), outputs.end(), [](const variable_t &v) {
          return v.get_type().is_integer() || v.get_type().is_bool();
        }))) {
      m_base_dom.intrinsic(name, inputs, outputs);
    } else {
      CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
    }
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const array_adaptive_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    m_base_dom += csts;

    CRAB_LOG("array-adaptive",
             crab::outs() << "assume(" << csts << ")  " << *this << "\n";);
  }

  bool entails(const linear_constraint_t &cst) const override {
    return m_base_dom.entails(cst);
  }
  
  void operator-=(const variable_t &var) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom()) {
      return;
    }

    if (var.get_type().is_array()) {
      forget_array(var);
    } else {
      m_base_dom -= var;
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    m_base_dom.assign(x, e);

    CRAB_LOG("array-adaptive", crab::outs() << "apply " << x << " := " << e
                                            << " " << *this << "\n";);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    m_base_dom.weak_assign(x, e);

    CRAB_LOG("array-adaptive", crab::outs() << "weak_assign(" << x << "," << e
	                                    << ")=" << *this << "\n";);
  }
    
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    m_base_dom.select(lhs, cond, e1, e2);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_assign(x, e, s.right_dom);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_apply(op, x, y, z, s.right_dom);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_apply(op, x, y, z, s.right_dom);    
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    m_base_dom.apply(op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    m_base_dom.apply(op, x, y, z);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << z << " " << *this << "\n";);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    m_base_dom.apply(op, x, y, k);

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "apply " << x << " := " << y << " " << op
                                   << " " << k << " " << *this << "\n";);
  }

  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    m_base_dom.assign_bool_cst(lhs, rhs);
  }

  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) override {
    m_base_dom.assign_bool_ref_cst(lhs, rhs);
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    m_base_dom.assign_bool_var(lhs, rhs, is_not_rhs);
  }

  virtual void weak_assign_bool_cst(const variable_t &lhs,
				    const linear_constraint_t &rhs) override {
    m_base_dom.weak_assign_bool_cst(lhs, rhs);
  }

  virtual void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
				    bool is_not_rhs) override {
    m_base_dom.weak_assign_bool_var(lhs, rhs, is_not_rhs);
  }
  
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    m_base_dom.apply_binary_bool(op, x, y, z);
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    m_base_dom.assume_bool(v, is_negated);
  }

  virtual void select_bool(const variable_t &lhs, const variable_t &cond,
                           const variable_t &b1,
                           const variable_t &b2) override {
    m_base_dom.select_bool(lhs, cond, b1, b2);
  }

  // backward boolean operators
  virtual void
  backward_assign_bool_cst(const variable_t &lhs,
                           const linear_constraint_t &rhs,
                           const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_assign_bool_cst(lhs, rhs, s.right_dom);
  }

  virtual void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_assign_bool_ref_cst(lhs, rhs, s.right_dom);
  }

  virtual void
  backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                           bool is_not_rhs,
                           const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_assign_bool_var(lhs, rhs, is_not_rhs, s.right_dom);
  }

  virtual void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const array_adaptive_domain_t &inv) override {
    auto s = rename_with_meet_semantics(inv);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    m_base_dom.backward_apply_binary_bool(op, x, y, z, s.right_dom);
  }

  /// array_adaptive is a functor domain that implements all
  /// operations except region/reference operations.
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(array_adaptive_domain_t)

  // array_operators_api

  // array_init returns a fresh array where all elements between
  // lb_idx and ub_idx are initialized to val. Thus, the first thing
  // we need to do is to kill existing cells.
  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".array_init");

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

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".array_load");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);
    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      // Check smashed array is consistent with elem_size
      const constant_value &a_elem_size = as.get_element_sz();
      constant_value cp_e_sz((int64_t)e_sz);
      cp_e_sz |= a_elem_size;
      if (!cp_e_sz.is_top()) {
        m_base_dom.array_load(lhs, a, elem_size, i);
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
                 crab::outs() << "Number of overlapping cells=" << cells.size()
                              << "\n";);

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
          variable_t rhs = mk_named_cell(a, o, e_sz, offset_map).second;
          // Here it's ok to do assignment (instead of expand)
          // because c is not a summarized variable. Otherwise, it
          // would be unsound.
          do_assign(lhs, rhs);
          m_array_map.set(a, next_as);
          goto array_load_end;
        }
      } else {
        linear_expression_t symb_lb(i);
        linear_expression_t symb_ub(i + number_t(e_sz - 1));
        std::vector<cell_t> cells;
        const offset_map_t &offset_map = as.get_offset_map();
        offset_map.get_overlap_cells_symbolic_offset(m_base_dom, symb_lb,
                                                     symb_ub, cells);
        // XXX: if we have a large array that is never smashed but we
        // do many reads with symbolic offsets then it might be better
        // to smash the array so that each read is cheaper.
        if (crab_domain_params_man::get().array_adaptive_is_smashable()) {
          if (array_state::can_be_smashed(cells, e_sz, true)) {
            // we smash all overlapping cells into a temporary array
            // (summarized) variable
            auto &vfac =
                const_cast<varname_t *>(&(a.name()))->get_var_factory();
            variable_t tmp_var(vfac.get(), a.get_type());
            bool found_cell_without_scalar = false;
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              auto c_scalar_opt = get_scalar(a, c);
              if (!c_scalar_opt) {
                CRAB_LOG("array-adaptive",
                         CRAB_WARN("array adaptive: ignored array load from ",
                                   a, " because non-constant array index ", i,
                                   "=", ii, " because found unnamed cell ",
                                   c););

                found_cell_without_scalar = true;
                break;
              }
              const bool is_strong_update = (k == 0);
              m_base_dom.array_store(tmp_var, elem_size, i, *c_scalar_opt,
                                     is_strong_update);
            }
            if (found_cell_without_scalar) {
              m_base_dom -= lhs;
            } else {
              // we read from the temporary summarized variable
              m_base_dom.array_load(lhs, tmp_var, elem_size, i);
            }
            // we forget the temporary summarized variable
            m_base_dom -= tmp_var;
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
    m_base_dom -= lhs;

  array_load_end:
    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << lhs << ":=" << a << "[" << i << "..." << ub
                          << "]  -- " << *this << "\n";);
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);
    const array_state &as = lookup_array_state(a);

    if (as.is_smashed()) {
      const constant_value &a_elem_size = as.get_element_sz();
      constant_value cp_e_sz((int64_t)e_sz);
      cp_e_sz |= a_elem_size;
      if (!cp_e_sz.is_top()) {
        m_base_dom.array_store(a, elem_size, i, val, is_strong_update);
      } else {
        m_base_dom -= a;
      }
    } else {
      interval_t ii = to_interval(i);
      array_state next_as(as);
      offset_map_t &offset_map = next_as.get_offset_map();
      boost::optional<number_t> n_opt = ii.singleton();
      if (n_opt &&
          (offset_map.get_number_cells() <
           crab_domain_params_man::get().array_adaptive_max_array_size())) {

        // -- Constant index: kill overlapping cells + perform strong update
        std::vector<cell_t> cells;
        offset_t o(static_cast<int64_t>(*n_opt));
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
        variable_t scalar_v = mk_named_cell(a, o, e_sz, offset_map).second;
        // -- strong update
        do_assign(scalar_v, val);
      } else {
        // -- Non-constant index: kill overlapping cells

        if (n_opt) {
          CRAB_LOG("array-adaptive",
                   crab::outs()
                       << "array write to " << a << " with constant index " << i
                       << "=" << ii << " but array size exceeded threshold of "
                       << crab_domain_params_man::get()
                              .array_adaptive_max_array_size()
                       << " so smashing is happening.\n";);
        } else {
          CRAB_LOG("array-adaptive", crab::outs() << "array write to " << a
                                                  << " with non-constant index "
                                                  << i << "=" << ii << "\n";);
        }

        bool smashed = false; // whether smashing took place
        if (crab_domain_params_man::get().array_adaptive_is_smashable()) {
          std::vector<cell_t> cells = offset_map.get_all_cells();
          if (next_as.can_be_smashed(e_sz) &&
              // Smashing is expensive because it will go over all cells
              // performing one join per weak update even if the smashed
              // array is already unconstrained. We don't smash if the
              // number of cells to be smashed is too large.
              (cells.size() <= crab_domain_params_man::get()
                                   .array_adaptive_max_smashable_cells())) {
            smashed = true;
            CRAB_LOG("array-adaptive-smash",
                     crab::outs() << "Array " << a << " will be smashed\n";);
            bool found_cell_without_scalar = false;
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              auto c_scalar_opt = get_scalar(a, c);
              if (!c_scalar_opt) {
                found_cell_without_scalar = true;
                break;
              }
              const bool is_strong_update = (k == 0);
              m_base_dom.array_store(a, elem_size, i, *c_scalar_opt,
                                     is_strong_update);
              CRAB_LOG("array-adaptive-smash",
                       crab::outs() << "\tAfter smashing " << *c_scalar_opt
                                    << "=" << m_base_dom << "\n";);
            }

            if (found_cell_without_scalar) {
              m_base_dom -= a;
            } else {
              // Finally the array store
              m_base_dom.array_store(a, elem_size, i, val, is_strong_update);
            }

            // The removal of cells from offset_map must be done after
            // the array has been fully smashed.
            for (unsigned k = 0, num_cells = cells.size(); k < num_cells; ++k) {
              const cell_t &c = cells[k];
              if (auto c_scalar_opt = get_scalar(a, c)) {
                // destroy the synthethic cell from the base domain and
                // from the offset map
                m_base_dom -= *c_scalar_opt;
              }
              offset_map.erase(c);
              m_cell_ghost_man.erase(a, c);
            }
            next_as.set_smashed(true);
            next_as.get_element_sz() = constant_value(number_t(e_sz));

            CRAB_LOG("array-adaptive-smash",
                     crab::outs() << "Array " << a << " has been smashed:"
                                  << m_base_dom << "\n";);
          } else {
            CRAB_LOG(
                "array-adaptive",
                if (cells.size() > crab_domain_params_man::get()
                                       .array_adaptive_max_smashable_cells()) {
                  CRAB_WARN("Array ", a,
                            " cannot be smashed because too many cells ",
                            cells.size(), ". Array write at index ", i, "=", ii,
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
          offset_map.get_overlap_cells_symbolic_offset(m_base_dom, symb_lb,
                                                       symb_ub, cells);
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

  // Perform array stores over an array segment [lb_idx, ub_idx]
  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &lb_idx,
                                 const linear_expression_t &ub_idx,
                                 const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store_range");

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
      CRAB_WARN("array adaptive store range ignored because lower bound ", *lb,
                " is not less or equal than upper bound ", *ub);
      return;
    }

    number_t num_elems = (*ub - *lb) / e_sz;
    number_t e = *ub;
    if (num_elems >
        crab_domain_params_man::get().array_adaptive_max_array_size()) {
      e = *lb +
          ((number_t(
                crab_domain_params_man::get().array_adaptive_max_array_size()) -
            1) *
           e_sz);
      CRAB_WARN("array adaptive store range will ignore indexes greater than ",
                e);
    }

    for (number_t i = *lb; i <= e;) {
      array_store(a, elem_size, i, val, false);
      i = i + e_sz;
    }
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
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
        variable_opt_t c_scalar_opt = get_scalar(rhs, c);
        if (!c_scalar_opt) {
          continue;
        }
        // Create a new cell for lhs from rhs's cell.
        auto named_cell =
            mk_named_cell(lhs, c.get_offset(), c.get_size(), lhs_om);
        variable_t new_c_scalar = named_cell.second;
        renmap.insert({new_c_scalar, *c_scalar_opt});
      }

      constant_value elem_sz = as.get_element_sz();
      m_array_map.set(
          lhs, array_state(false, std::move(elem_sz), std::move(lhs_om)));

      CRAB_LOG(
          "array-adaptive-array-assign", crab::outs() << "array variables={";
          std::vector<variable_t> array_variables = get_array_variables();
          for (unsigned i = 0, e = array_variables.size(); i < e;
               ++i) { crab::outs() << array_variables[i] << ";"; } crab::outs()
          << "}\n";);

      for (auto &kv : renmap) {
        CRAB_LOG("array-adaptive-array-assign",
                 crab::outs() << "Base domain assign " << kv.first
                              << " := " << kv.second << "\n";);
        do_assign(kv.first, kv.second);
      }
    } else if (crab_domain_params_man::get().array_adaptive_is_smashable()) {
      CRAB_LOG("array-adaptive-array-assign", crab::outs() << "Smashed\n";);
      m_base_dom.array_assign(lhs, rhs);
      m_array_map.set(lhs, as);
    }
    CRAB_LOG("array-adaptive", crab::outs() << "Res=" << *this << "\n";);
  }

  // backward array operations

  virtual void
  backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                      const linear_expression_t &lb_idx,
                      const linear_expression_t &ub_idx,
                      const linear_expression_t &val,
                      const array_adaptive_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_init");

    if (is_bottom()) {
      return;
    }

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

  virtual void
  backward_array_load(const variable_t &lhs, const variable_t &a,
                      const linear_expression_t &elem_size,
                      const linear_expression_t &i,
                      const array_adaptive_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_load");

    if (is_bottom()) {
      return;
    }

    auto s = rename_with_meet_semantics(invariant);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
       
    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      CRAB_WARN("array_adaptive::backward_array_load not implemented if array "
                "smashed");
    } else {
      // We use the forward invariant to extract the array index.
      // it's ok that invariant is not renamed here.
      interval_t ii = to_interval(i, invariant.get_content_domain());
      if (boost::optional<number_t> n = ii.singleton()) {
        array_state next_as(as);
        offset_map_t &om = next_as.get_offset_map();
        offset_t o(static_cast<int64_t>(*n));
        uint64_t e_sz = check_and_get_elem_size(elem_size);
        // We need to ensure when do_backward_assign is called that
        // m_base_dom and s.right_dom (renamed version of invariant.m_base_dom)
        // have been properly renamed.
        cell_t c = mk_named_cell(a, o, e_sz, om).first;
        do_backward_assign(lhs, a, c, s.right_dom /*invariant.m_base_dom*/);
        m_array_map.set(a, next_as);
      } else {
        CRAB_LOG("array-adaptive",
                 CRAB_WARN("array index is not a constant value"););
        // -- Forget lhs
        m_base_dom -= lhs;
        // -- Meet with forward invariant
        //
        // m_base_dom and s.right_dom have been properly renamed so we
        // don't need need the next line which would do again the
        // renaming.
        //*this = *this & invariant;
        m_base_dom = m_base_dom & s.right_dom;
      }
    }

    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << "BACKWARD " << lhs << ":=" << a << "[" << i
                          << "..." << ub << "]  -- " << *this << "\n";);
  }

  virtual void backward_array_store(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &val,
      bool /*is_strong_update*/,
      const array_adaptive_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store");

    if (is_bottom()) {
      return;
    }

    auto s = rename_with_meet_semantics(invariant);
    m_base_dom = std::move(s.left_dom);
    m_array_map = std::move(s.array_map);
    m_cell_ghost_man = std::move(s.cell_ghost_man);
    
    uint64_t e_sz = check_and_get_elem_size(elem_size);
    const array_state &as = lookup_array_state(a);
    if (as.is_smashed()) {
      CRAB_WARN("array_adaptive::backward_array_store not implemented if array "
                "smashed");
    } else {
      array_state next_as(as);
      offset_map_t &om = next_as.get_offset_map();
      // We use the forward invariant to extract the array index.
      // it's ok to use "invariant" without being renamed.
      interval_t ii = to_interval(i, invariant.m_base_dom);
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
          // m_base_dom and s.right_dom have been properly renamed so we
          // don't need the next line which would do again the
          // renaming.
          //
          //*this = *this & invariant;
          m_base_dom = m_base_dom & s.right_dom;
        } else {
          // We need to ensure when do_backward_assign is called that
          // m_base_dom and s.right_dom (renamed version of
          // invariant.m_base_dom) have been properly renamed.
          cell_t c = mk_named_cell(a, o, e_sz, om).first;
          do_backward_assign(a, c, val, s.right_dom /*invariant.m_base_dom*/);
        }
      } else {
        // TODOX: smash the array if needed
        // -- Non-constant index or multiple overlapping cells: kill
        // -- overlapping cells and meet with forward invariant.
        linear_expression_t symb_lb(i);
        linear_expression_t symb_ub(i + number_t(e_sz - 1));
        std::vector<cell_t> cells;
        om.get_overlap_cells_symbolic_offset(
            s.right_dom /*invariant.m_base_dom*/, symb_lb, symb_ub, cells);
        kill_cells(a, cells, om);
        // -- Meet with forward invariant
        //
        // m_base_dom and s.right_dom have been properly renamed so we
        // don't need need the next line which would do again the
        // renaming.
        //
        //*this = *this & invariant;
        m_base_dom = m_base_dom & s.right_dom;
      }
      m_array_map.set(a, next_as);
    }

    CRAB_LOG("array-adaptive", linear_expression_t ub = i + elem_size - 1;
             crab::outs() << "BACKWARD " << a << "[" << i << "..." << ub
                          << "]:=" << val << " -- " << *this << "\n";);
  }

  virtual void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
      const linear_expression_t &val,
      const array_adaptive_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".count.backward_array_store_range");

    if (is_bottom())
      return;

    uint64_t e_sz = check_and_get_elem_size(elem_size);

    // make copy to avoid one extra copy
    base_domain_t base_dom(invariant.m_base_dom);
    interval_t lb_i = to_interval(lb_idx, base_dom);
    auto lb = lb_i.singleton();
    if (!lb) {
      return;
    }

    interval_t ub_i = to_interval(ub_idx, base_dom);
    auto ub = ub_i.singleton();
    if (!ub) {
      return;
    }

    if (!(*lb <= *ub)) {
      return;
    }

    number_t num_elems = (*ub - *lb) / e_sz;
    number_t e = *ub;
    if (num_elems >
        crab_domain_params_man::get().array_adaptive_max_array_size()) {
      e = *lb +
          ((number_t(
                crab_domain_params_man::get().array_adaptive_max_array_size()) -
            1) *
           e_sz);
    }

    for (number_t i = *lb; i <= e;) {
      backward_array_store(a, elem_size, i, val, false, invariant);
      i = i + e_sz;
    }
  }

  virtual void
  backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                        const array_adaptive_domain_t &invariant) override {
    CRAB_WARN("backward_array_assign in array_adaptive domain not implemented");
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    return filter_nonscalar_vars(
        std::move(m_base_dom.to_linear_constraint_system()));
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    disjunctive_linear_constraint_system_t res;
    auto disj_csts = m_base_dom.to_disjunctive_linear_constraint_system();
    for (auto &csts : disj_csts) {
      auto filtered_csts = filter_nonscalar_vars(std::move(csts));
      if (!filtered_csts.is_true()) {
        res += filtered_csts;
      }
    }
    return res;
  }

  base_domain_t get_content_domain() const { return m_base_dom; }

  base_domain_t &get_content_domain() { return m_base_dom; }

  void write(crab_os &o) const override {
    o << m_base_dom;

    CRAB_LOG(
        "array-adaptive-print-details", crab::outs() << "\n";
        crab::outs() << "=== CELLS PER ARRAY === \n";
        for (auto it = m_array_map.begin(), et = m_array_map.end(); it != et;
             ++it) {
          const variable_t &v = it->first;
          const array_state &as = it->second;
          crab::outs() << "ARRAY VAR " << v << ":\n";
          as.write(crab::outs());
          crab::outs() << "\n";
        } crab::outs()
        << "=== CELL GHOST VARIABLES === \n";
        m_cell_ghost_man.write(crab::outs()););
  }

  std::string domain_name() const override {
    std::string name("ArrayAdaptive(" + m_base_dom.domain_name() + ")");
    return name;
  }

  void expand(const variable_t &v, const variable_t &new_v) override {
    if (is_bottom() || is_top()) {
      return;
    }

    if (v.get_type() != new_v.get_type()) {
      CRAB_ERROR(domain_name(), "::expand must preserve same type");
    }

    if (v.get_type().is_array()) {
      CRAB_WARN(domain_name(), "::expand not implemented for array variable");
    } else {
      m_base_dom.expand(v, new_v);
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (is_bottom() || is_top()) {
      return;
    }

    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename expects vectors same sizes");
    }

    CRAB_LOG("array-adaptive", crab::outs()
                                   << "Before renaming " << *this << "\n";);

    // Split into array and scalar variables
    variable_vector_t old_array_vars, new_array_vars, old_base_vars,
        new_base_vars;
    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      variable_t old_v = from[i];
      variable_t new_v = to[i];
      if (old_v.get_type() != new_v.get_type()) {
        CRAB_ERROR(domain_name(), "::rename must preserve same type");
      }
      if (new_v.get_type().is_array()) {
        old_array_vars.push_back(old_v);
        new_array_vars.push_back(new_v);
      } else {
        old_base_vars.push_back(old_v);
        new_base_vars.push_back(new_v);
      }
    }

    unsigned num_arr_vars = old_array_vars.size();
    for (unsigned i = 0; i < num_arr_vars; ++i) {
      variable_t old_v = old_array_vars[i];
      variable_t new_v = new_array_vars[i];

      // Rename m_cell_ghost_man and m_array_map
      if (const array_state *old_as = m_array_map.find(old_v)) {
        array_state new_as;
        new_as.set_smashed(old_as->is_smashed());
        constant_value &new_cp_dom = new_as.get_element_sz();
        new_cp_dom = old_as->get_element_sz();
        offset_map_t &offset_map = new_as.get_offset_map();

        if (!old_as->is_smashed()) {
          m_cell_ghost_man.rename(
              old_v, new_v,
              [&offset_map](const offset_t &o, uint64_t sz) -> cell_t {
                return offset_map.mk_cell(o, sz);
              },
              old_base_vars, new_base_vars);
        } else {
          old_base_vars.push_back(old_v);
          new_base_vars.push_back(new_v);
        }

        m_array_map -= old_v;
        m_array_map.set(new_v, new_as);
      }
      m_cell_ghost_man.erase_all(old_v);
    } // end for

    m_base_dom.rename(old_base_vars, new_base_vars);
    CRAB_LOG("array-adaptive", crab::outs()
                                   << "After renaming " << *this << "\n";);
  }

}; // end array_adaptive_domain

template <typename Dom>
struct abstract_domain_traits<array_adaptive_domain<Dom>> {
  using number_t = typename Dom::number_t;
  using varname_t = typename Dom::varname_t;
};

// template <typename Dom>
// class checker_domain_traits<array_adaptive_domain<Dom>> {
// public:
//   using this_type = array_adaptive_domain<Dom>;
//   using base_domain_t = typename this_type::base_domain_t;
//   using linear_constraint_t = typename this_type::linear_constraint_t;
//   using disjunctive_linear_constraint_system_t =
//       typename this_type::disjunctive_linear_constraint_system_t;
//   static bool entail(this_type &lhs,
//                      const disjunctive_linear_constraint_system_t &rhs) {
//     base_domain_t &lhs_dom = lhs.get_content_domain();
//     return checker_domain_traits<base_domain_t>::entail(lhs_dom, rhs);
//   }
//   static bool entail(const disjunctive_linear_constraint_system_t &lhs,
//                      this_type &rhs) {
//     base_domain_t &rhs_dom = rhs.get_content_domain();
//     return checker_domain_traits<base_domain_t>::entail(lhs, rhs_dom);
//   }
//   static bool intersect(this_type &inv, const linear_constraint_t &cst) {
//     base_domain_t &dom = inv.get_content_domain();
//     return checker_domain_traits<base_domain_t>::intersect(dom, cst);
//   }
// };

} // namespace domains
} // namespace crab
