#pragma once

/*******************************************************************************
 *  Non-overlapping sequence of intervals
 *******************************************************************************/

#include <crab/domains/interval.hpp>
#include <crab/domains/linear_interval_solver.hpp>
#include <crab/fixpoint/thresholds.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {

template <typename Number> class dis_interval {
public:
  using bound_t = ikos::bound<Number>;
  using interval_t = ikos::interval<Number>;
  using dis_interval_t = dis_interval<Number>;

private:
  using state_t = enum { BOT, FINITE, TOP };
  using list_intervals_t = std::vector<interval_t>;

public:
  using iterator = typename list_intervals_t::iterator;
  using const_iterator = typename list_intervals_t::const_iterator;

private:
  state_t m_state;
  list_intervals_t m_list;

  bool are_consecutive(const interval_t &i1, const interval_t &i2) const;

  bool overlap(const interval_t &i1, const interval_t &i2) const;

  bool check_well_formed(const dis_interval_t &x) const;

  list_intervals_t normalize(list_intervals_t l, bool &is_bottom) const;

  dis_interval(state_t state);

  dis_interval(list_intervals_t l, bool Normalize = true);

  // pre: x is normalized
  interval_t approx(list_intervals_t x) const;

  struct WidenOp {
    virtual ~WidenOp(){}
    virtual interval_t apply(const interval_t &before,
                             const interval_t &after) = 0;
  };

  struct BasicWidenOp : public WidenOp {
    virtual ~BasicWidenOp(){}
    virtual interval_t apply(const interval_t &before,
                             const interval_t &after) override;
  };

  struct WidenWithThresholdsOp : public WidenOp {
    const thresholds<Number> &m_ts;
    WidenWithThresholdsOp(const thresholds<Number> &ts);
    virtual interval_t apply(const interval_t &before,
                             const interval_t &after) override;
  };

  dis_interval_t widening(const dis_interval_t &o, WidenOp &widen_op) const;

  dis_interval_t
  apply_bin_op(const dis_interval_t &x, const dis_interval_t &y,
               std::function<ikos::interval<Number>(ikos::interval<Number>,
                                                    ikos::interval<Number>)>
                   op,
               bool shortcut_top) const;

  dis_interval_t apply_unary_op(
      const dis_interval_t &x,
      std::function<ikos::interval<Number>(ikos::interval<Number>)> op) const;

public:
  dis_interval();

  dis_interval(interval_t i);

  dis_interval(const dis_interval_t &i) = default;

  dis_interval(dis_interval_t &&i) = default;

  dis_interval_t &operator=(const dis_interval_t &i) = default;

  dis_interval_t &operator=(dis_interval_t &&i) = default;

  static dis_interval_t top();

  static dis_interval_t bottom();

  bool is_bottom() const;

  bool is_top() const;

  bool is_finite() const;

  dis_interval_t lower_half_line() const;

  dis_interval_t upper_half_line() const;

  boost::optional<Number> singleton() const;

  iterator begin();

  iterator end();

  const_iterator begin() const;

  const_iterator end() const;

  interval_t approx() const;

  bool operator==(const dis_interval_t &o) const;

  bool operator<=(const dis_interval_t &o) const;

  dis_interval_t operator|(const dis_interval_t &o) const;

  dis_interval_t operator&(const dis_interval_t &o) const;

  dis_interval_t operator||(const dis_interval_t &o) const;

  dis_interval_t widening_thresholds(const dis_interval_t &o,
                                     const crab::thresholds<Number> &ts) const;

  dis_interval_t operator&&(const dis_interval_t &o) const;

  dis_interval_t operator+(const dis_interval_t &x) const;

  dis_interval_t &operator+=(const dis_interval_t &x);

  dis_interval_t operator-() const;

  dis_interval_t operator-(const dis_interval_t &x) const;

  dis_interval_t &operator-=(const dis_interval_t &x);

  dis_interval_t operator*(const dis_interval_t &x) const;

  dis_interval_t &operator*=(const dis_interval_t &x);

  dis_interval_t operator/(const dis_interval_t &x);

  dis_interval_t &operator/=(const dis_interval_t &x);

  dis_interval_t UDiv(const dis_interval_t &x) const;

  dis_interval_t SRem(const dis_interval_t &x) const;

  dis_interval_t URem(const dis_interval_t &x) const;

  dis_interval_t And(const dis_interval_t &x) const;

  dis_interval_t Or(const dis_interval_t &x) const;

  dis_interval_t Xor(const dis_interval_t &x) const;

  dis_interval_t Shl(const dis_interval_t &x) const;

  dis_interval_t LShr(const dis_interval_t &x) const;

  dis_interval_t AShr(const dis_interval_t &x) const;

  void normalize();

  void write(crab_os &o) const;

  friend crab_os &operator<<(crab_os &o, const dis_interval<Number> &i) {
    i.write(o);
    return o;
  }
};
} // end namespace domains
} // end namespace crab

namespace ikos {
namespace linear_interval_solver_impl {
template <>
crab::domains::dis_interval<z_number>
trim_interval(const crab::domains::dis_interval<z_number> &x,
              const crab::domains::dis_interval<z_number> &y);
template <>
crab::domains::dis_interval<q_number>
trim_interval(const crab::domains::dis_interval<q_number> &i,
              const crab::domains::dis_interval<q_number> &j);
template <>
crab::domains::dis_interval<z_number>
lower_half_line(const crab::domains::dis_interval<z_number> &i, bool is_signed);
template <>
crab::domains::dis_interval<q_number>
lower_half_line(const crab::domains::dis_interval<q_number> &i, bool is_signed);
template <>
crab::domains::dis_interval<z_number>
upper_half_line(const crab::domains::dis_interval<z_number> &i, bool is_signed);
template <>
crab::domains::dis_interval<q_number>
upper_half_line(const crab::domains::dis_interval<q_number> &i, bool is_signed);
} // namespace linear_interval_solver_impl
} // end namespace ikos
