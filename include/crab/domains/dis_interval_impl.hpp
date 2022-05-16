#pragma once

#include <boost/range/iterator_range.hpp>
#include <crab/domains/dis_interval.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace domains {

template <typename Number> dis_interval<Number> dis_interval<Number>::top() {
  return dis_interval<Number>(dis_interval<Number>::TOP);
}

template <typename Number> dis_interval<Number> dis_interval<Number>::bottom() {
  return dis_interval<Number>(dis_interval<Number>::BOT);
}

// REVISIT: this is correct only if Number=z_number
template <typename Number>
bool dis_interval<Number>::are_consecutive(
    const ikos::interval<Number> &i1, const ikos::interval<Number> &i2) const {
  return ((i1.lb() <= i2.lb() && i1.ub() <= i2.ub()) &&
          (i1.ub() + Number(1) == i2.lb())) ||
         ((i2.lb() <= i1.lb() && i2.ub() <= i1.ub()) &&
          (i2.ub() + Number(1) == i1.lb()));
}

template <typename Number>
bool dis_interval<Number>::overlap(const ikos::interval<Number> &i1,
                                   const ikos::interval<Number> &i2) const {
  return !((i1 & i2).is_bottom());
}

template <typename Number> struct IsOnTheLeft {
  bool operator()(const ikos::interval<Number> &i1,
                  const ikos::interval<Number> &i2) const {
    return ((i1.ub() <= i2.lb()) && (i1.ub() != i2.lb()));
  }
};

template <typename Number>
bool dis_interval<Number>::check_well_formed(
    const dis_interval<Number> &x) const {

  if (x.is_top() || x.is_bottom()) {
    return true;
  }

  if (!x.is_finite()) {
    CRAB_ERROR("sanity check -- state should be finite\n", x,
               " not well formed");
    return false;
  }

  if (x.m_list.empty()) {
    CRAB_ERROR("sanity check -- list cannot be empty\n", x, " not well formed");
    return false;
  }

  if (x.m_list.size() == 1) {
    if (x.m_list[0].is_top() || x.m_list[0].is_bottom()) {
      CRAB_ERROR("sanity check -- cannot be top or bottom\n", x,
                 " not well formed");
      return false;
    } else {
      return true;
    }
  }

  // -- check the list of intervals is strictly sorted. This also
  // -- checks for duplicates but it won't complain with two
  // -- intervals like [0,2], [3,4]
  IsOnTheLeft<Number> comp;
  list_intervals_t tmp(x.m_list);
  unsigned prev = 0;
  for (unsigned cur = 1; cur < tmp.size(); cur++, prev++) {
    if (comp(tmp[prev], tmp[cur]))
      continue;

    CRAB_ERROR("sanity check -- list is not strictly sorted: ", tmp[prev],
               " not leq ", tmp[cur], "\n", x, " not well formed");
    return false;
  }

  return true;
}

template <typename Number>
typename dis_interval<Number>::list_intervals_t dis_interval<Number>::normalize(
    typename dis_interval<Number>::list_intervals_t l, bool &is_bottom) const {

  if (l.size() <= 1) {
    CRAB_LOG("disint", crab::outs()
                           << "-- Normalize: singleton " << l[0] << "\n");
    is_bottom = false;
    return l;
  }

  std::sort(
      l.begin(), l.end(),
      [](const ikos::interval<Number> &a, const ikos::interval<Number> &b) {
        return (a.lb() <= b.lb() && a.lb() != b.lb());
      });

  typename dis_interval<Number>::list_intervals_t res;
  ikos::interval<Number> prev = ikos::interval<Number>::top();
  unsigned int bottoms = 0;

  for (unsigned int i = 0; i < l.size(); ++i) {
    ikos::interval<Number> intv = l[i];

    if (prev == intv) {
      CRAB_LOG("disint", crab::outs() << "-- Normalize: duplicate"
                                      << "\n");
      continue;
    }

    if (intv.is_bottom()) {
      CRAB_LOG("disint", crab::outs() << "-- Normalize: bottom interval"
                                      << "\n");
      bottoms++;
      continue;
    }

    if (intv.is_top()) {
      CRAB_LOG("disint", crab::outs() << "-- Normalize: top interval"
                                      << "\n");
      is_bottom = false;
      return typename dis_interval<Number>::list_intervals_t();
    }

    if (!prev.is_top()) {

      bool refined = true;
      while (refined && res.size() > 0) {
        prev = res[res.size() - 1];
        if (overlap(prev, intv) || are_consecutive(prev, intv)) {
          CRAB_LOG("disint",
                   crab::outs()
                       << "-- Normalize: overlapping or consecutive intervals"
                       << prev << " and " << intv << "\n");
          res.pop_back();
          intv = prev | intv;
        } else if (intv <= prev) {
          CRAB_LOG("disint", crab::outs()
                                 << "-- Normalize: skipping subsumed interval"
                                 << prev << " and " << intv << "\n");
          goto next_iter;
        } else {
          refined = false;
        }
      }
    }

    prev = intv;
    CRAB_LOG("disint", crab::outs() << "-- Normalize: adding " << intv << "\n");
    if (intv.is_top()) {
      CRAB_LOG("disint", crab::outs() << "-- Normalize: top interval"
                                      << "\n");
      is_bottom = false;
      return typename dis_interval<Number>::list_intervals_t();
    }

    res.push_back(intv);
  next_iter:;
    ;
  }

  is_bottom = (bottoms == l.size());
  CRAB_LOG("disint", crab::outs() << "-- Normalize: number of bottoms = "
                                  << bottoms << "\n");
  // crab::outs() << "-- Normalize result="; for (auto i: res) { crab::outs()
  // << i << "|";} crab::outs() << "\n";
  return res;
}

template <typename Number>
dis_interval<Number>::dis_interval(state_t state) : m_state(state) {}

template <typename Number>
dis_interval<Number>::dis_interval(
    typename dis_interval<Number>::list_intervals_t l, bool apply_normalize)
    : m_state(FINITE), m_list(l) {

  // TODO: make it a template parameter
  const std::size_t max_num_disjunctions = 50;

  if (apply_normalize) {
    bool is_bottom = false;
    typename dis_interval<Number>::list_intervals_t res =
        normalize(m_list, is_bottom);
    if (is_bottom) {
      m_state = BOT;
      m_list.clear();
    } else if (res.empty()) {
      m_state = TOP;
      m_list.clear();
    } else {
      m_list = std::move(res);
    }
  }

  if (is_finite() && m_list.size() >= max_num_disjunctions) {
    // TODO: rather than merging all intervals do a more graceful
    // degradation e.g., start by merging the nearest intervals.
    CRAB_WARN(" reached maximum allowed number of disjunctions. ",
              "Merging all intervals ... ");
    ikos::interval<Number> alljoined = approx(m_list);
    m_list.clear();
    m_list.push_back(alljoined);
  }

  assert(check_well_formed(*this));
}

// pre: x is normalized
template <typename Number>
ikos::interval<Number> dis_interval<Number>::approx(
    typename dis_interval<Number>::list_intervals_t x) const {
  if (x.empty())
    CRAB_ERROR("list should not be empty");

  if (x.size() <= 1)
    return x[0];
  else
    return (x[0] | x[x.size() - 1]);
}

template <typename Number>
dis_interval<Number>::dis_interval() : m_state(TOP) {}

template <typename Number>
dis_interval<Number>::dis_interval(ikos::interval<Number> i) : m_state(FINITE) {
  if (i.is_top())
    m_state = TOP;
  else if (i.is_bottom())
    m_state = BOT;
  else
    m_list.push_back(i);
}

template <typename Number> bool dis_interval<Number>::is_bottom() const {
  return (m_state == BOT);
}

template <typename Number> bool dis_interval<Number>::is_top() const {
  return (m_state == TOP);
}

template <typename Number> bool dis_interval<Number>::is_finite() const {
  return (m_state == FINITE);
}

template <typename Number>
dis_interval<Number> dis_interval<Number>::lower_half_line() const {
  if (is_bottom()) {
    return bottom();
  } else {
    auto f = [](ikos::interval<Number> x) -> ikos::interval<Number> {
      return x.lower_half_line();
    };
    return apply_unary_op(*this, f);
  }
}

template <typename Number>
dis_interval<Number> dis_interval<Number>::upper_half_line() const {
  if (is_bottom()) {
    return bottom();
  } else {
    auto f = [](ikos::interval<Number> x) -> ikos::interval<Number> {
      return x.upper_half_line();
    };
    return apply_unary_op(*this, f);
  }
}

template <typename Number>
boost::optional<Number> dis_interval<Number>::singleton() const {
  ikos::interval<Number> i = approx();
  return i.singleton();
}

template <typename Number>
typename dis_interval<Number>::iterator dis_interval<Number>::begin() {
  return m_list.begin();
}

template <typename Number>
typename dis_interval<Number>::iterator dis_interval<Number>::end() {
  return m_list.end();
}

template <typename Number>
typename dis_interval<Number>::const_iterator
dis_interval<Number>::begin() const {
  return m_list.begin();
}

template <typename Number>
typename dis_interval<Number>::const_iterator
dis_interval<Number>::end() const {
  return m_list.end();
}

template <typename Number>
ikos::interval<Number> dis_interval<Number>::approx() const {
  if (is_bottom())
    return ikos::interval<Number>::bottom();
  if (is_top())
    return ikos::interval<Number>::top();
  return approx(m_list);
}

// required for separate_domains
// pre: *this and o are normalized
template <typename Number>
bool dis_interval<Number>::operator==(const dis_interval<Number> &o) const {
#if 0
       // -- semantic check
       return (*this <= o && o <= *this);
#else
  // -- syntactic check
  if (is_bottom() && o.is_bottom())
    return true;
  else if (is_top() && o.is_top())
    return true;
  else if (is_bottom() || o.is_bottom())
    return false;
  else if (is_top() || o.is_top())
    return false;
  else {
    if (m_list.size() != o.m_list.size())
      return false;
    else {
      unsigned int j = 0;
      for (unsigned int i = 0; i < m_list.size(); ++i, ++j) {
        if (!(m_list[i] == o.m_list[j]))
          return false;
      }
      return true;
    }
  }
#endif
}

template <typename Number>
bool dis_interval<Number>::operator<=(const dis_interval<Number> &o) const {
  if (this->is_bottom()) {
    return true;
  } else if (o.is_bottom()) {
    return false;
  } else {

    unsigned j = 0;
    for (unsigned int i = 0; i < m_list.size(); i++) {
      for (; j < o.m_list.size(); j++) {
        if (m_list[i] <= o.m_list[j])
          goto next_iter;
      }
      return false; // m_list[i] is not included in any interval of o.m_list
    next_iter:;
      ;
    }
    return true;
  }
}

// pre: *this and o are normalized
template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator|(const dis_interval<Number> &o) const {

  CRAB_LOG("disint", crab::outs()
                         << "Join of " << *this << " and " << o << "\n");

  if (this->is_bottom()) {
    return o;
  } else if (o.is_bottom()) {
    return *this;
  } else if (this->is_top()) {
    return *this;
  } else if (o.is_top()) {
    return o;
  } else {
    IsOnTheLeft<Number> isLeftCmp;
    unsigned int i = 0;
    unsigned int j = 0;
    typename dis_interval<Number>::list_intervals_t res;
    res.reserve(m_list.size() + o.m_list.size());

    while (i < m_list.size() && j < o.m_list.size()) {
      CRAB_LOG("disint", crab::outs()
                             << "Join -- left operand =" << m_list[i]
                             << " right operand=" << o.m_list[j] << "\n");

      if (m_list[i].is_top() || o.m_list[j].is_top()) {
        CRAB_LOG("disint", crab::outs() << "Join -- One of the operands is top"
                                        << "\n");
        return dis_interval<Number>();
      } else if (m_list[i].is_bottom()) {
        CRAB_LOG("disint", crab::outs() << "Join -- Left operand is bottom"
                                        << "\n");
        i++;
      } else if (o.m_list[j].is_bottom()) {
        CRAB_LOG("disint", crab::outs() << "Join -- Right operand is bottom"
                                        << "\n");
        j++;
      } else if (m_list[i] == o.m_list[j]) {
        CRAB_LOG("disint", crab::outs()
                               << "Join -- Left operand is equal to right"
                               << "\n");
        res.push_back(m_list[i]);
        i++;
        j++;
      } else if (m_list[i] <= o.m_list[j]) {
        CRAB_LOG("disint",
                 crab::outs() << "Join -- Left operand is included in the right"
                              << "\n");
        res.push_back(o.m_list[j]);
        i++;
        j++;
      } else if (o.m_list[j] <= m_list[i]) {
        CRAB_LOG("disint",
                 crab::outs() << "Join -- Right operand is included in the left"
                              << "\n");
        res.push_back(m_list[i]);
        i++;
        j++;
      } else if (overlap(m_list[i], o.m_list[j]) ||
                 are_consecutive(m_list[i], o.m_list[j])) {
        CRAB_LOG("disint", crab::outs()
                               << "Join -- Left " << m_list[i] << " and right "
                               << o.m_list[j]
                               << " operands overlap or are consecutive"
                               << "\n");
        res.push_back(m_list[i] | o.m_list[j]);
        i++;
        j++;
      } else if (isLeftCmp(m_list[i], o.m_list[j])) {
        CRAB_LOG("disint", crab::outs()
                               << "Join -- Left operand " << m_list[i]
                               << " is on the left of the right operand "
                               << o.m_list[j] << "\n");
        res.push_back(m_list[i]);
        i++;
      } else {
        assert(isLeftCmp(o.m_list[j], m_list[i]));

        CRAB_LOG("disint", crab::outs()
                               << "Join -- Right operand " << o.m_list[j]
                               << " is on the left of the right operand "
                               << m_list[i] << "\n");

        res.push_back(o.m_list[j]);
        j++;
      }
    }

    // consume the rest of left operand
    while (i < m_list.size()) {
      auto intv = m_list[i];

      CRAB_LOG("disint", crab::outs()
                             << "Join -- Adding the rest of left operand "
                             << intv << "\n");

      bool refined = true;
      while (refined && res.size() > 0) {
        auto prev = res[res.size() - 1];
        if (overlap(prev, intv) || are_consecutive(prev, intv)) {
          CRAB_LOG("disint",
                   crab::outs()
                       << "\t-- Join : overlapping or consecutive intervals"
                       << prev << " and " << intv << "\n");
          res.pop_back();
          intv = prev | intv;
        } else if (intv <= prev) {
          CRAB_LOG("disint", crab::outs()
                                 << "\t-- Join: skipping subsumed interval"
                                 << prev << " and " << intv << "\n");
          goto next_iter_1;
        } else {
          refined = false;
        }
      }

      res.push_back(intv);
    next_iter_1:
      i++;
    }

    // consume the rest of right operand
    while (j < o.m_list.size()) {
      auto intv = o.m_list[j];
      CRAB_LOG("disint", crab::outs()
                             << "Join -- Adding the rest of right operands "
                             << intv << "\n");
      bool refined = true;
      while (refined && res.size() > 0) {
        auto prev = res[res.size() - 1];
        if (overlap(prev, intv) || are_consecutive(prev, intv)) {
          CRAB_LOG("disint",
                   crab::outs()
                       << "\t-- Join : overlapping or consecutive intervals"
                       << prev << " and " << intv << "\n");
          res.pop_back();
          intv = prev | intv;
        } else if (intv <= prev) {
          CRAB_LOG("disint", crab::outs()
                                 << "\t-- Join: skipping subsumed interval"
                                 << prev << " and " << intv << "\n");
          goto next_iter_2;
        } else {
          refined = false;
        }
      }

      res.push_back(intv);
    next_iter_2:
      j++;
    }

    if (res.empty()) {
      CRAB_LOG("disint", crab::outs() << "Join result=_|_"
                                      << "\n");
      return dis_interval<Number>(BOT);
    } else if (res.size() <= 1 && res[0].is_top()) {
      CRAB_LOG("disint", crab::outs() << "Joing result=[-oo,+oo]"
                                      << "\n");
      return dis_interval<Number>(TOP);
    } else {
      // It needs normalization. E.g., the join of {[0, 7] | [9, 11]}
      // and {[0, 6] | [8, 11]} returns {[0, 7] | [8, 11]} which
      // should be further simplified to [0,11].

      dis_interval<Number> join(res);
      CRAB_LOG("disint", crab::outs() << "Joing result=" << join << "\n");
      return join;
    }
  }
}

// pre: *this and o are normalized
template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator&(const dis_interval<Number> &o) const {
  if (this->is_bottom() || o.is_bottom()) {
    return this->bottom();
  } else if (this->is_top()) {
    return o;
  } else if (o.is_top()) {
    return *this;
  } else {

    typename dis_interval<Number>::list_intervals_t res;
    res.reserve(m_list.size() + o.m_list.size());

    bool is_bot = true;
    for (unsigned int i = 0; i < m_list.size(); ++i) {
      for (unsigned int j = 0; j < o.m_list.size(); ++j) {
        auto meet = m_list[i] & o.m_list[j];
        if (!meet.is_bottom()) {
          res.push_back(meet);
          is_bot = false;
        }
      }
    }

    if (is_bot) {
      return (this->bottom());
    } else if (res.empty()) {
      return (this->top());
    } else {
      return dis_interval<Number>(res);
    }
  }
}

template <typename Number>
ikos::interval<Number>
dis_interval<Number>::BasicWidenOp::apply(const ikos::interval<Number> &before,
                                          const ikos::interval<Number> &after) {
  return before || after;
}

template <typename Number>
dis_interval<Number>::WidenWithThresholdsOp::WidenWithThresholdsOp(
    const thresholds<Number> &ts)
    : m_ts(ts) {}

template <typename Number>
ikos::interval<Number> dis_interval<Number>::WidenWithThresholdsOp::apply(
    const ikos::interval<Number> &before, const ikos::interval<Number> &after) {
  return before.widening_thresholds(after, m_ts);
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::widening(const dis_interval<Number> &o,
                               WidenOp &widen_op) const {
  if (this->is_bottom()) {
    return o;
  } else if (o.is_bottom()) {
    return *this;
  } else if (this->is_top()) {
    return *this;
  } else if (o.is_top()) {
    return o;
  } else {

    // --- trivial cases first
    if (m_list.size() == 1 && o.m_list.size() == 1) {
      return dis_interval<Number>(widen_op.apply(m_list[0], o.m_list[0]));
    }

    if (m_list.size() == 1 && o.m_list.size() > 1) {
      ikos::interval<Number> x = approx(o.m_list);
      return dis_interval<Number>(widen_op.apply(m_list[0], x));
    }

    if (m_list.size() > 1 && o.m_list.size() == 1) {
      ikos::interval<Number> x = approx(m_list);
      return dis_interval<Number>(widen_op.apply(x, o.m_list[0]));
    }

    assert(m_list.size() >= 2 && o.m_list.size() >= 2);

    // The widening implemented in CodeContracts widens the
    // extremes and keep only stable intervals. For this query:
    // widening( [1,1] | [4, 4], [1,1] | [3,3] | [4, 4]) the
    // result would to be [1,1] | [4,4]. But this is not even an
    // upper bound of the right argument so it cannot be a
    // widening.

    // -- widen the extremes
    ikos::interval<Number> lb_widen = widen_op.apply(m_list[0], o.m_list[0]);
    ikos::interval<Number> ub_widen = widen_op.apply(
        m_list[m_list.size() - 1], o.m_list[o.m_list.size() - 1]);

    typename dis_interval<Number>::list_intervals_t res;
    res.reserve(m_list.size() + o.m_list.size());

    res.push_back(lb_widen);

    // for (unsigned int i=1; i < m_list.size() - 1; i++) {
    //   for (unsigned int j=1; j < o.m_list.size() - 1; j++) {
    //     if (o.m_list [j] <= m_list [i]) {
    //       res.push_back(m_list [i]);
    //       break;
    //     }
    //   }
    // }

    // keep all the intervals, normalize will do the rest
    res.insert(res.end(), m_list.begin() + 1, m_list.end() - 1);
    res.insert(res.end(), o.m_list.begin() + 1, o.m_list.end() - 1);

    res.push_back(ub_widen);

    return dis_interval<Number>(res);
  }
}

// pre: *this and o are normalized
template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator||(const dis_interval<Number> &o) const {
  BasicWidenOp op;
  return widening(o, op);
}

// pre: *this and o are normalized
template <typename Number>
dis_interval<Number>
dis_interval<Number>::widening_thresholds(const dis_interval<Number> &o,
                                          const thresholds<Number> &ts) const {
  WidenWithThresholdsOp op(ts);
  return widening(o, op);
}

// pre: *this and o are normalized
template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator&&(const dis_interval<Number> &o) const {
  // CRAB_WARN(" DisIntervals narrowing operator replaced with meet");
  return (*this & o);
}

// pre: x and y are normalized
template <typename Number>
dis_interval<Number> dis_interval<Number>::apply_bin_op(
    const dis_interval<Number> &x, const dis_interval<Number> &y,
    std::function<ikos::interval<Number>(ikos::interval<Number>,
                                         ikos::interval<Number>)>
        op,
    bool shortcut_top) const {

  // if shortcut_top is true then the result is top if one of the
  // operands is top

  if (x.is_bottom() || y.is_bottom())
    return this->bottom();

  if (x.is_top() && y.is_top())
    return this->top();

  if (shortcut_top && (x.is_top() || y.is_top()))
    return this->top();

  typename dis_interval<Number>::list_intervals_t res;
  res.reserve(x.m_list.size() + y.m_list.size());
  bool is_bot = true;

  if (shortcut_top || (x.is_finite() && y.is_finite())) {
    for (unsigned int i = 0; i < x.m_list.size(); ++i) {
      for (unsigned int j = 0; j < y.m_list.size(); ++j) {
        ikos::interval<Number> intv = op(x.m_list[i], y.m_list[j]);
        if (intv.is_bottom())
          continue;

        if (intv.is_top())
          return this->top();

        is_bot = false;
        res.push_back(intv);
      }
    }
  } else {

    assert(x.is_top() || y.is_top());

    dis_interval<Number> non_top_arg = (!x.is_top() ? x : y);
    bool is_non_top_left = (!x.is_top());
    for (unsigned int i = 0; i < non_top_arg.m_list.size(); ++i) {

      ikos::interval<Number> intv =
          (is_non_top_left
               ? op(non_top_arg.m_list[i], ikos::interval<Number>::top())
               : op(ikos::interval<Number>::top(), non_top_arg.m_list[i]));

      if (intv.is_bottom())
        continue;

      if (intv.is_top())
        return this->top();

      is_bot = false;
      res.push_back(intv);
    }
  }

  if (is_bot)
    return this->bottom();
  else
    return dis_interval<Number>(res);
}

// pre: x is normalized
template <typename Number>
dis_interval<Number> dis_interval<Number>::apply_unary_op(
    const dis_interval<Number> &x,
    std::function<ikos::interval<Number>(ikos::interval<Number>)> op) const {

  if (x.is_bottom())
    return this->bottom();

  if (x.is_top())
    return this->top();

  assert(x.is_finite());

  if (x.m_list.empty())
    CRAB_ERROR("list should not be empty");

  typename dis_interval<Number>::list_intervals_t res;
  res.reserve(x.m_list.size());
  bool is_bot = true;

  for (unsigned int i = 0; i < x.m_list.size(); ++i) {
    ikos::interval<Number> intv = op(x.m_list[i]);
    if (intv.is_bottom())
      continue;
    if (intv.is_top())
      return this->top();
    is_bot = false;
    res.push_back(intv);
  }

  if (is_bot)
    return this->bottom();
  else
    return dis_interval<Number>(res);
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator+(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a + b;
    };
    return apply_bin_op(*this, x, f, true);
  }
}

template <typename Number>
dis_interval<Number> &
dis_interval<Number>::operator+=(const dis_interval<Number> &x) {
  return this->operator=(this->operator+(x));
}

template <typename Number>
dis_interval<Number> dis_interval<Number>::operator-() const {
  if (this->is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a) -> ikos::interval<Number> {
      return -a;
    };
    return apply_unary_op(*this, f);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator-(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a - b;
    };
    return apply_bin_op(*this, x, f, true);
  }
}

template <typename Number>
dis_interval<Number> &
dis_interval<Number>::operator-=(const dis_interval<Number> &x) {
  return this->operator=(this->operator-(x));
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator*(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a * b;
    };
    return apply_bin_op(*this, x, f, true);
  }
}

template <typename Number>
dis_interval<Number> &
dis_interval<Number>::operator*=(const dis_interval<Number> &x) {
  return this->operator=(this->operator*(x));
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::operator/(const dis_interval<Number> &x) {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a / b;
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number> &
dis_interval<Number>::operator/=(const dis_interval<Number> &x) {
  return this->operator=(this->operator/(x));
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::UDiv(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a / b;
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::SRem(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.SRem(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::URem(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.URem(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::And(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.And(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::Or(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.Or(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::Xor(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.Xor(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::Shl(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.Shl(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::LShr(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.LShr(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number>
dis_interval<Number>
dis_interval<Number>::AShr(const dis_interval<Number> &x) const {
  if (this->is_bottom() || x.is_bottom()) {
    return this->bottom();
  } else {
    auto f = [](ikos::interval<Number> a, ikos::interval<Number> b) {
      return a.AShr(b);
    };
    return apply_bin_op(*this, x, f, false);
  }
}

template <typename Number> void dis_interval<Number>::normalize() {
  if (is_bottom() || is_top())
    return;
  bool is_bot = false;
  typename dis_interval<Number>::list_intervals_t res =
      normalize(m_list, is_bot);

  if (is_bot) {
    m_state = dis_interval<Number>::BOT;
  } else if (res.empty()) {
    m_state = dis_interval<Number>::TOP;
  } else {
    m_state = dis_interval<Number>::FINITE;
    std::swap(m_list, res);
  }
}

template <typename Number> void dis_interval<Number>::write(crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
  } else if (is_top()) {
    o << "[-oo,+oo]";
  } else {
    assert(m_state == dis_interval<Number>::FINITE);

    for (typename dis_interval<Number>::list_intervals_t::const_iterator
             it = m_list.begin(),
             et = m_list.end();
         it != et;) {
      o << *it;
      ++it;
      if (it != et)
        o << " | ";
    }
  }
}

} // end namespace domains
} // end namespace crab
