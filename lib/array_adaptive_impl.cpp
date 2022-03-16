#include <crab/domains/array_adaptive.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace domains {
namespace array_adaptive_impl {

using bound_t = ikos::bound<ikos::z_number>;

void constant_value::set_to_top() {
  m_is_bottom = false;
  m_val = bound_t::plus_infinity();
}

void constant_value::set_to_bot() {
  m_is_bottom = true;
  m_val = bound_t::plus_infinity();
}

constant_value::constant_value(bool is_bottom)
    : m_is_bottom(is_bottom), m_val(bound_t::plus_infinity()) {}

constant_value::constant_value()
    : m_is_bottom(false), m_val(bound_t::plus_infinity()) {}

constant_value::constant_value(int64_t sz) : m_is_bottom(false), m_val(sz) {}

constant_value::constant_value(const ikos::interval<ikos::z_number> &sz)
    : m_is_bottom(false), m_val(bound_t::plus_infinity()) {

  if (auto v_opt = sz.singleton()) {
    m_val = bound_t(*v_opt);
  }
}

bool constant_value::is_top() const {
  return (!is_bottom() && (m_val == bound_t::plus_infinity()));
}

bool constant_value::is_bottom() const { return m_is_bottom; }

bool constant_value::is_zero() const {
  return !m_is_bottom && m_val == bound_t(0);
}

bool constant_value::is_negative() const {
  return !m_is_bottom && m_val < bound_t(0);
}

bound_t constant_value::val() const {
  if (is_bottom()) {
    CRAB_ERROR("constant_value::val cannot be called on bottom");
  }
  return m_val;
}

boost::optional<uint64_t> constant_value::get_uint64_val() const {
  if (is_bottom()) {
    CRAB_ERROR("constant_value::get_uint64_val cannot be called on bottom");
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

constant_value constant_value::bottom() { return constant_value(true); }

constant_value constant_value::top() { return constant_value(false); }

bool constant_value::operator==(const constant_value &o) const {
  return *this <= o && o <= *this;
}

bool constant_value::operator<=(const constant_value &o) const {
  if (is_bottom() || o.is_top()) {
    return true;
  } else if (is_top() || o.is_bottom()) {
    return false;
  } else {
    return (m_val == o.m_val);
  }
}

void constant_value::operator|=(const constant_value &o) {
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

constant_value constant_value::operator|(const constant_value &o) const {
  if (is_bottom()) {
    return o;
  } else if (o.is_bottom()) {
    return *this;
  } else if (is_top() || o.is_top()) {
    return constant_value();
  } else {
    if (m_val == o.m_val) {
      return *this;
    } else {
      return constant_value();
    }
  }
}

constant_value constant_value::operator&(const constant_value &o) const {
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

constant_value constant_value::operator||(const constant_value &o) const {
  return operator|(o);
}

constant_value constant_value::operator&&(const constant_value &o) const {
  return operator&(o);
}

void constant_value::write(crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else {
    o << m_val;
  }
}

offset_t::offset_t(ikos::index_t v) : m_val(v) {}

ikos::index_t offset_t::index() const { return m_val; }

size_t offset_t::hash() const {
  // casting to size_t may overflow but it shouldn't affect correctness
  return std::hash<size_t>{}(static_cast<size_t>(m_val));
}

bool offset_t::operator<(const offset_t &o) const { return m_val < o.m_val; }

bool offset_t::operator==(const offset_t &o) const { return m_val == o.m_val; }

bool offset_t::operator!=(const offset_t &o) const { return !(*this == o); }

offset_t offset_t::operator%(const offset_t &o) const {
  return offset_t(m_val % o.m_val);
}

offset_t offset_t::operator-(const offset_t &o) const {
  return offset_t(m_val - o.m_val);
}

bool offset_t::is_negative() const { return m_val < 0; }

bool offset_t::is_zero() const { return m_val == 0; }

void offset_t::write(crab::crab_os &o) const { o << m_val; }

cell_t::cell_t() : m_offset(0), m_size(0), m_removed(false) {}

cell_t::cell_t(offset_t offset, uint64_t size)
    : m_offset(offset), m_size(size), m_removed(false) {}

cell_t::interval_t cell_t::to_interval(const offset_t o, uint64_t size) {
  cell_t::interval_t i(o.index(), o.index() + size - 1);
  return i;
}

cell_t::interval_t cell_t::to_interval() const {
  return to_interval(get_offset(), get_size());
}

bool cell_t::is_null() const { return (m_offset.index() == 0 && m_size == 0); }

offset_t cell_t::get_offset() const { return m_offset; }

size_t cell_t::get_size() const { return m_size; }

cell_t cell_t::clone(void) const {
  cell_t new_cell;
  new_cell.m_offset = m_offset;
  new_cell.m_size = m_size;
  new_cell.m_removed = m_removed;
  return new_cell;
}

void cell_t::mark_as_removed(bool v) { m_removed = v; }

bool cell_t::is_removed(void) const { return m_removed; }

size_t cell_t::hash() const {
  size_t h1 = m_offset.hash();
  size_t h2 = std::hash<size_t>{}(m_size);
  return h1 ^ (h2 << 1);
}

// inclusion test
bool cell_t::operator<=(const cell_t &o) const {
  cell_t::interval_t x = to_interval();
  cell_t::interval_t y = o.to_interval();
  return x <= y;
}

bool cell_t::operator==(const cell_t &o) const {
  return (get_offset() == o.get_offset() && get_size() == o.get_size());
}

bool cell_t::operator<(const cell_t &o) const {
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
bool cell_t::overlap(const offset_t &o, uint64_t size) const {
  if (m_removed)
    return false;

  cell_t::interval_t x = to_interval();
  cell_t::interval_t y = to_interval(o, size);
  bool res = (!(x & y).is_bottom());
  CRAB_LOG("array-adaptive-overlap", crab::outs() << "**Checking if " << x
                                                  << " overlaps with " << y
                                                  << "=" << res << "\n";);
  return res;
}

void cell_t::write(crab::crab_os &o) const {
  if (is_null()) {
    o << "NULL";
  } else {
    o << to_interval() << " ";
    if (m_removed) {
      o << "R";
    }
  }
}

// Delete completely the cell
void offset_map_t::erase_cell(const cell_t &c) {
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
void offset_map_t::remove_cell(const cell_t &c) {
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

void offset_map_t::insert_cell(const cell_t &c) {
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

cell_t offset_map_t::get_cell(const offset_t &o, uint64_t size) const {
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
cell_t offset_map_t::mk_cell(const offset_t &o, uint64_t size /*bytes*/) {
  cell_t c = get_cell(o, size);
  if (c.is_null()) {
    cell_t nc(o, size);
    insert_cell(nc);
    c = nc;
    CRAB_LOG("array-adaptive", crab::outs() << "**Created cell " << c << "\n";);
  }
  if (c.is_null()) {
    CRAB_ERROR("cannot create a null cell from offset ", o, " and size ", size);
  }
  return c;
}

offset_map_t::offset_map_t(patricia_tree_t &&m) : m_map(std::move(m)) {}

offset_map_t::offset_map_t() {}

offset_map_t::offset_map_t(const offset_map_t &o) : m_map(o.m_map) {}

offset_map_t::offset_map_t(const offset_map_t &&o)
    : m_map(std::move(o.m_map)) {}

offset_map_t &offset_map_t::operator=(const offset_map_t &o) {
  if (this != &o) {
    m_map = o.m_map;
  }
  return *this;
}

offset_map_t &offset_map_t::operator=(const offset_map_t &&o) {
  if (this != &o) {
    m_map = std::move(o.m_map);
  }
  return *this;
}

bool offset_map_t::empty() const { return m_map.empty(); }

std::size_t offset_map_t::size() const { return m_map.size(); }

// leq operator
bool offset_map_t::operator<=(const offset_map_t &o) const {
  domain_po po;
  return m_map.leq(o.m_map, po);
}

// set union: if two cells with same offset do not agree on
// size then they are ignored.
offset_map_t offset_map_t::operator|(const offset_map_t &o) const {
  join_op op;
  offset_map_t res = offset_map_t(apply_operation(op, m_map, o.m_map));
  CRAB_LOG("array-adaptive-offset-map", crab::outs() << "offset_map join\n"
                                                     << *this << "\nand\n"
                                                     << o << "\nResult=" << res
                                                     << "\n";);
  return res;
}

// set intersection: if two cells with same offset do not agree
// on size then they are ignored.
offset_map_t offset_map_t::operator&(const offset_map_t &o) const {
  meet_op op;
  offset_map_t res = offset_map_t(apply_operation(op, m_map, o.m_map));
  CRAB_LOG("array-adaptive-offset-map", crab::outs() << "offset_map meet\n"
                                                     << *this << "\nand\n"
                                                     << o << "\nResult=" << res
                                                     << "\n";);
  return res;
}

// Completely delete the cell from the offset map
void offset_map_t::erase(const cell_t &c) { erase_cell(c); }

void offset_map_t::erase(const std::vector<cell_t> &cells) {
  for (unsigned i = 0, e = cells.size(); i < e; ++i) {
    erase(cells[i]);
  }
}

// Pretend the cell is removed so no other cells overlap with it but
// the cell is not actually deleted from the offset map. We need to
// know all created cells when smashing occurs.
void offset_map_t::remove(const cell_t &c) { remove_cell(c); }

void offset_map_t::remove(const std::vector<cell_t> &cells) {
  for (unsigned i = 0, e = cells.size(); i < e; ++i) {
    remove(cells[i]);
  }
}

// cells are sorted by offset
std::vector<cell_t> offset_map_t::get_all_cells() const {
  std::vector<cell_t> res;
  for (auto it = m_map.begin(), et = m_map.end(); it != et; ++it) {
    auto const &o_cells = it->second;
    for (auto &c : o_cells) {
      res.push_back(c);
    }
  }
  return res;
}

unsigned offset_map_t::get_number_cells() const {
  unsigned count = 0;
  for (auto it = m_map.begin(), et = m_map.end(); it != et; ++it) {
    auto const &o_cells = it->second;
    count += o_cells.size();
  }
  return count;
}

// Return in out all cells that might overlap with (o, size).
//
// It is not marked as const because we insert temporary a cell.
// However, upon completion this method leaves unmodified the object.
void offset_map_t::get_overlap_cells(const offset_t &o, uint64_t size,
                                     std::vector<cell_t> &out) {
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

void offset_map_t::clear(void) { m_map.clear(); }

void offset_map_t::write(crab::crab_os &o) const {
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

} // end namespace array_adaptive_impl
} // end namespace domains
} // namespace crab
