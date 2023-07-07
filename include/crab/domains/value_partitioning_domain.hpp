#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
namespace value_partitioning_domain_impl {
template <typename D> class partition;
} // end namespace value_partitioning_domain_impl

#define VALUE_PARTITION_START "value_partition_start"
#define VALUE_PARTITION_END "value_partition_end"
/*
 * State partitioning based on the values (over-approximated as an
 * interval) of a single *numerical* variable (see Tutorial on Static
 * Inference of Numeric Invariants by Abstract Interpretation by
 * A. Mine, section 6.3.3). The partitioning variable is chosen
 * dynamically via a special intrinsics.  Unlike Mine's description in
 * his tutorial, the partitions are created and maintained
 * dynamically. To ensure termination, widening is applied
 * element-wise only if both operands have the same
 * partitions. Otherwise, all partitions are merged before applying
 * widening in the base domain.  The concretization is the union of
 * all partitions.
 */
template <typename NumDomain> class product_value_partitioning_domain;

template <typename NumDomain>
class value_partitioning_domain final
    : public abstract_domain_api<value_partitioning_domain<NumDomain>> {
  friend class product_value_partitioning_domain<NumDomain>;

public:
  using value_partitioning_domain_t = value_partitioning_domain<NumDomain>;
  using abstract_domain_api_t =
      abstract_domain_api<value_partitioning_domain_t>;
  using typename abstract_domain_api_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_t;
  using typename abstract_domain_api_t::linear_expression_t;
  using typename abstract_domain_api_t::number_t;
  using typename abstract_domain_api_t::reference_constraint_t;
  using typename abstract_domain_api_t::variable_or_constant_t;
  using typename abstract_domain_api_t::variable_or_constant_vector_t;
  using typename abstract_domain_api_t::variable_t;
  using typename abstract_domain_api_t::variable_vector_t;
  using typename abstract_domain_api_t::varname_t;

  using partition_t = value_partitioning_domain_impl::partition<NumDomain>;
  using interval_t = ikos::interval<number_t>;
  using partition_vector = std::vector<partition_t>;
  using partition_iterator = typename partition_vector::iterator;
  using const_partition_iterator = typename partition_vector::const_iterator;

private:
  // the partitioning variable
  boost::optional<variable_t> m_variable;
  // partitions over m_variable's values: are sorted and they don't
  // overlap
  partition_vector m_partitions;

  value_partitioning_domain(boost::optional<variable_t> variable,
                            partition_vector &&partitions)
      : m_variable(variable), m_partitions(std::move(partitions)) {}

  /* update all partitions after the partitioning variable has been modified */
  void update_partitions(void) {
    if (!m_variable) {
      return;
    }
    CRAB_LOG("partition",
             crab::outs() << "Before updating partitions: " << *this << "\n";);

    // Update intervals and remove bottom partitions.
    // We traverse in reverse order so that removing partitions is constant
    for (auto it = m_partitions.end(); it != m_partitions.begin();) {
      --it;
      interval_t interval = it->get_dom()[*m_variable];
      if (interval.is_bottom()) {
        if (m_partitions.size() > 1) {
          it = m_partitions.erase(it);
        } else {
          it->get_interval() = interval_t::top();
          return;
        }
      } else {
        it->get_interval() = std::move(interval);
      }
    }

    // Sort partitions
    std::sort(m_partitions.begin(), m_partitions.end(),
              [](const partition_t &a, const partition_t &b) {
                return a.get_interval().lb() < b.get_interval().lb();
              });

    // Merge partitions if necessary.
    // We traverse in reverse order so that removing partitions is constant
    auto it = m_partitions.end();
    --it;
    for (; it != m_partitions.begin() && m_partitions.size() > 1;) {
      --it;
      auto next_it = it;
      ++next_it;
      if (it->get_interval().ub() >= next_it->get_interval().lb()) {
        // merge two partitions
        it->join_interval(*next_it);
	auto &v1 = it->get_dom();
	auto const&v2 = next_it->get_dom();
        v1 |= v2;
        it = (m_partitions.erase(next_it));
	--it;
      }
    }
    CRAB_LOG("partition",
             crab::outs() << "After updating partitions: " << *this << "\n";);
  }

  // For meet, widening, and narrowing. Joins and inclusion are
  // handled separately.
  value_partitioning_domain_t apply_binary_op(
      const value_partitioning_domain_t &other, const std::string &opName,
      std::function<interval_t(const interval_t &, const interval_t &)>
          intervalOp,
      std::function<NumDomain(const NumDomain &, const NumDomain &)> domainOp)
      const {

    // trivial cases (bot/top) handled in the caller
    assert(!is_bottom());
    assert(!is_top());
    assert(!other.is_bottom());
    assert(!other.is_top());

    if (!m_variable && !other.m_variable) {
      CRAB_LOG("partition", crab::outs() << opName << " (no partitions) "
                                         << *this << " and " << other << "=";);
      // no partitions
      value_partitioning_domain_t out(*this);
      out.m_partitions[0].get_dom() = domainOp(out.m_partitions[0].get_dom(),
                                               other.m_partitions[0].get_dom());
      CRAB_LOG("partition", crab::outs() << out << "\n";);
      return out;
    } else if (!(m_variable == other.m_variable)) {
      CRAB_LOG("partition",
               crab::outs()
                   << opName << " (operands with different partition variable) "
                   << *this << " and " << other << "=";);
      value_partitioning_domain_t out(*this);
      out.remove_partitions();
      out.m_partitions[0].get_dom() = domainOp(
          out.m_partitions[0].get_dom(), other.merge_partitions().get_dom());
      CRAB_LOG("partition", crab::outs() << out << "\n";);
      return out;
    } else if (has_same_partitions(other)) {
      CRAB_LOG("partition", crab::outs()
                                << opName << " (operands with same partitions) "
                                << *this << " and " << other << "=";);
      // element-wise op if both operands have the same partitions
      partition_vector out_partitions;
      out_partitions.reserve(m_partitions.size());
      for (unsigned i = 0, sz = m_partitions.size(); i < sz; i++) {
        auto range = m_partitions[i].get_interval();
        partition_t out_partition(
            std::move(range),
            std::move(domainOp(m_partitions[i].get_dom(),
                               other.m_partitions[i].get_dom())));
        out_partitions.emplace_back(std::move(out_partition));
      }
      value_partitioning_domain_t out(m_variable, std::move(out_partitions));
      CRAB_LOG("partition", crab::outs() << out << "\n";);
      return out;
    } else {
      CRAB_LOG("partition", crab::outs()
                                << opName
                                << " (operands with different partitions) "
                                << *this << " and " << other << "=";);
      // otherwise, merge all partitions in both operands and apply op
      value_partitioning_domain_t out(*this);
      partition_t smashed_out = out.merge_partitions();
      out.m_partitions.clear();
      out.m_partitions.emplace_back(std::move(smashed_out));
      partition_t smashed_other = other.merge_partitions();
      out.m_partitions[0].get_interval() = intervalOp(
          out.m_partitions[0].get_interval(), smashed_other.get_interval());
      out.m_partitions[0].get_dom() =
          domainOp(out.m_partitions[0].get_dom(), smashed_other.get_dom());
      CRAB_LOG("partition", crab::outs() << out << "\n";);
      return out;
    }
  }

  bool contains(const linear_constraint_t &cst, const variable_t &v) const {
    auto sorted_variables = cst.variables();
    auto lower = std::lower_bound(sorted_variables.begin(),
				  sorted_variables.end(), v);
    return (!(lower == sorted_variables.end() || v < *lower));
  }
  
  bool contains(const linear_constraint_system_t &csts, const variable_t &v) const {
    return std::any_of(csts.begin(), csts.end(),
		       [this,&v](const linear_constraint_t &cst) {
			 return contains(cst, v);
		       });
  }

  void add_constraints(const linear_constraint_system_t &csts) {
    assert(!is_bottom());
    assert(!csts.is_true());
    assert(!csts.is_false());
    
    CRAB_LOG("partition", crab::outs() << "Adding " << csts << "\n";);
    partition_vector vec;
    vec.reserve(m_partitions.size());
    for (auto &partition : m_partitions) {
      partition.get_dom() += csts;
      // remove bottom partitions
      if (partition.get_dom().is_bottom()) {
        continue;
      }
      vec.emplace_back(partition);
    }
    std::swap(m_partitions, vec);

    // update partitions
    if (m_variable) {
      
      if (contains(csts, *m_variable)) {
        update_partitions();
      }
    }
    CRAB_LOG("partition", crab::outs() << "Res=" << *this << "\n";);
  }

  void split_relevant_disequalities(const linear_constraint_system_t &csts,
				    const variable_t &v,
				    linear_constraint_system_t &diseqs,
				    linear_constraint_system_t &non_diseqs) {
    for (auto it=csts.begin(), et=csts.end();it!=et;++it) {
      const linear_constraint_t &cst = *it;
      if (!cst.is_disequation()) {
	non_diseqs += cst;
      } else {
	if (contains(cst, v)) {
	  diseqs += cst;
	} else {
	  non_diseqs += cst;
	}
      }
    }
  }
  
public:
  value_partitioning_domain() : m_variable(boost::none) {
    m_partitions.reserve(1);
    m_partitions.emplace_back(partition_t::top());
  }
  value_partitioning_domain(const value_partitioning_domain_t &other) = default;
  value_partitioning_domain(value_partitioning_domain_t &&other) = default;
  value_partitioning_domain_t &
  operator=(const value_partitioning_domain_t &other) = default;
  value_partitioning_domain_t &
  operator=(value_partitioning_domain_t &&other) = default;

  bool is_bottom() const override {
    for (auto &partition : m_partitions) {
      if (!partition.get_dom().is_bottom()) {
        return false;
      }
    }
    return true;
  }

  bool is_top() const override {
    for (auto &partition : m_partitions) {
      if (!partition.get_dom().is_top()) {
        return false;
      }
    }
    return true;
  }

  value_partitioning_domain_t make_bottom() const override {
    value_partitioning_domain_t res;
    res.set_to_bottom();
    return res;
  }

  value_partitioning_domain_t make_top() const override {
    value_partitioning_domain_t res;
    return res;
  }

  void set_to_top() override {
    m_partitions.clear();
    m_partitions.emplace_back(partition_t::top());
  }

  void set_to_bottom() override {
    m_partitions.clear();
    m_partitions.emplace_back(partition_t::bottom());
  }

  /* Begin specialized API of the value_partitioning_domain */

  boost::optional<variable_t> get_variable() const { return m_variable; }

  // passed by reference so that we can set m_variable if none
  boost::optional<variable_t> &get_variable() { return m_variable; }

  bool has_same_partitions(const value_partitioning_domain_t &other) const {
    if (m_variable != other.m_variable ||
        m_partitions.size() != other.m_partitions.size()) {
      return false;
    }
    for (unsigned i = 0, sz = m_partitions.size(); i < sz; ++i) {
      if (!(m_partitions[i].get_interval() ==
            other.m_partitions[i].get_interval())) {
        return false;
      }
    }
    return true;
  }

  void remove_partitions(void) {
    if (!m_variable) {
      if (m_partitions.size() != 1) {
        CRAB_ERROR("remove_partitions expects partition of size 1");
      }
      return;
    }
    partition_t partition = merge_partitions();
    m_partitions.clear();
    m_partitions.emplace_back(std::move(partition));
    m_partitions[0].get_interval() = interval_t::top();
    m_variable = boost::none;
  }

  partition_t merge_partitions() const {
    if (m_partitions.empty()) {
      CRAB_ERROR("merge_partitions() expects at least one partition");
    }

    partition_t smashed_p = m_partitions[0];
    for (unsigned i = 1, sz = m_partitions.size(); i < sz; ++i) {
      smashed_p.join(m_partitions[i]);
    }
    return smashed_p;
  }

  partition_iterator partitions_begin() {
    if (is_bottom()) {
      CRAB_ERROR("partitions_begin() on bottom");
    }
    return m_partitions.begin();
  }

  partition_iterator partitions_end() {
    if (is_bottom()) {
      CRAB_ERROR("partitions_end() on bottom");
    }
    return m_partitions.end();
  }

  const_partition_iterator partitions_begin() const {
    if (is_bottom()) {
      CRAB_ERROR("partitions_begin() on bottom");
    }
    return m_partitions.begin();
  }

  const_partition_iterator partitions_end() const {
    if (is_bottom()) {
      CRAB_ERROR("partitions_end() on bottom");
    }
    return m_partitions.end();
  }
  /* End specialized API of the value_partitioning_domain */

  bool operator<=(const value_partitioning_domain_t &other) const override {
    if (is_bottom()) {
      return true;
    } else if (other.is_top()) {
      return true;
    } else if (!m_variable && !other.m_variable) {
      return (m_partitions[0].get_dom() <= other.m_partitions[0].get_dom());
    } else if (!(m_variable == other.m_variable)) {
      partition_t smashed_this = merge_partitions();
      partition_t smashed_other = other.merge_partitions();
      return smashed_this.get_dom() <= smashed_other.get_dom();
    } else {
      for (auto this_it = m_partitions.begin(), this_et = m_partitions.end(),
                other_it = other.m_partitions.begin(),
                other_et = other.m_partitions.end();
           this_it != this_et;) {
        if (other_it == other_et ||
            this_it->get_interval().ub() < other_it->get_interval().lb()) {
          if (!this_it->get_dom().is_bottom()) {
            return false;
          }
          ++this_it;
        } else if (other_it->get_interval().ub() <
                   this_it->get_interval().lb()) {
          ++other_it;
        } else if (this_it->get_interval().ub() <=
                   other_it->get_interval().ub()) {
          if (!(this_it->get_dom() <= other_it->get_dom())) {
            return false;
          }
          ++this_it;
        } else {
          // The partition on the left overlaps one or more partitions on the
          // right
          NumDomain other_dom = other_it->get_dom();
          for (auto it = ++other_it;
               it != other_et &&
               this_it->get_interval().ub() >= it->get_interval().lb();
               ++it) {
            other_dom |= it->get_dom();
          }
          if (!(this_it->get_dom() <= other_dom)) {
            return false;
          }
          ++this_it;
        }
      }
      return true;
    }
  }

  void operator|=(const value_partitioning_domain_t &other) override {
    CRAB_LOG("partition",
             crab::outs() << "Join " << *this << " and " << other << "=\n";);
    if (is_bottom()) {
      *this = other;
    } else if (other.is_bottom()) {
      // do nothing
    } else if (!m_variable && !other.m_variable) {
      CRAB_LOG("partition", crab::outs() << "\tNo partitions\n");
      m_partitions[0].get_dom() |= other.m_partitions[0].get_dom();
    } else if (!(m_variable == other.m_variable)) {
      CRAB_LOG("partition",
               crab::outs() << "\tPartitioning on two different variables\n";);
      remove_partitions(); // ensure only one partition from now on and
                           // m_partitions[0].get_interval() is top
      for (auto &partition : other.m_partitions) {
        m_partitions[0].get_dom() |= partition.get_dom();
      }
    } else {
      CRAB_LOG("partition", crab::outs()
                                << "\tPartitioning on the same variable\n";);
      auto left_it = m_partitions.begin();
      auto right_it = other.m_partitions.begin();
      auto right_et = other.m_partitions.end();
      while (right_it != right_et) {
        if (left_it == m_partitions.end()) {
          // No more on the left so we just append the rest on the right
          m_partitions.insert(m_partitions.end(), right_it, right_et);
          break;
        } else if (left_it->get_interval().ub() <
                   right_it->get_interval().lb()) {
          // The interval on the left appears before than the one on the right
          ++left_it;
        } else if (right_it->get_interval().ub() <
                   left_it->get_interval().lb()) {
          // The interval on the right appears before than the one on the left
          left_it = m_partitions.insert(left_it, *right_it);
	  ++left_it; // skip inserted element
          ++right_it;
        } else {
          // The interval on the left and the one on the right overlap
	  left_it->join_interval(*right_it);
          auto it = left_it;
          ++it;
          for (; it != m_partitions.end() &&
                 left_it->get_interval().ub() >= it->get_interval().lb();) {
	    left_it->join(*it);
            it = m_partitions.erase(it);
            left_it = it;
            --left_it;
          }
          // Join with the partition on the right
          left_it->get_dom() |= right_it->get_dom();
          ++right_it;
        }
      }
    }
    CRAB_LOG("partition", crab::outs() << *this << "\n";);
  }

  value_partitioning_domain_t
  operator|(const value_partitioning_domain_t &other) const override {
    value_partitioning_domain_t tmp(*this);
    tmp |= other;
    return tmp;
  }

  void operator&=(const value_partitioning_domain_t &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (is_top() || other.is_bottom()) {
      *this = other;
    } else {
      // similar to apply_binary_op but avoiding the copy of the left
      // operand.
      if (!m_variable && !other.m_variable) {
        CRAB_LOG("partition", crab::outs() << " Meet (no partitions) " << *this
                                           << " and " << other << "=";);
        // no partitions
        m_partitions[0].get_dom() &= other.m_partitions[0].get_dom();
        CRAB_LOG("partition", crab::outs() << *this << "\n";);
      } else if (!(m_variable == other.m_variable)) {
        CRAB_LOG("partition",
                 crab::outs()
                     << " Meet (operands with different partition variable) "
                     << *this << " and " << other << "=";);
        remove_partitions();
        m_partitions[0].get_dom() &= other.merge_partitions().get_dom();
        CRAB_LOG("partition", crab::outs() << *this << "\n";);
      } else if (has_same_partitions(other)) {
        CRAB_LOG("partition", crab::outs()
                                  << " Meet (operands with same partitions) "
                                  << *this << " and " << other << "=";);
        // element-wise op if both operands have the same partitions
        for (unsigned i = 0, sz = m_partitions.size(); i < sz; i++) {
          m_partitions[i].get_dom() &= other.m_partitions[i].get_dom();
        }
        CRAB_LOG("partition", crab::outs() << *this << "\n";);
      } else {
        CRAB_LOG("partition",
                 crab::outs() << " Meet (operands with different partitions) "
                              << *this << " and " << other << "=";);
        // otherwise, merge all partitions in both operands and apply op
        partition_t smashed_left = merge_partitions();
        m_partitions.clear();
        m_partitions.emplace_back(std::move(smashed_left));
        partition_t smashed_right = other.merge_partitions();
        m_partitions[0].meet(smashed_right);
        CRAB_LOG("partition", crab::outs() << *this << "\n";);
      }
    }
  }

  value_partitioning_domain_t
  operator&(const value_partitioning_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      // The solution here is sub-optimal
      return apply_binary_op(
          other, "Meet",
          [](const interval_t &i1, const interval_t &i2) { return i1 & i2; },
          [](const NumDomain &v1, const NumDomain &v2) { return v1 & v2; });
    }
  }

  value_partitioning_domain_t
  operator||(const value_partitioning_domain_t &other) const override {
    CRAB_LOG("partition",
             crab::outs() << "Widening " << *this << " and " << other << "\n";);
    if (is_bottom() || other.is_top()) {
      CRAB_LOG("partition", crab::outs() << "Result=" << other << "\n";);
      return other;
    } else if (other.is_bottom() || is_top()) {
      CRAB_LOG("partition", crab::outs() << "Result=" << *this << "\n";);
      return *this;
    } else {
      return apply_binary_op(
          other, "Widening",
          [](const interval_t &i1, const interval_t &i2) { return i1 || i2; },
          [](const NumDomain &v1, const NumDomain &v2) { return v1 || v2; });
    }
  }

  value_partitioning_domain_t
  widening_thresholds(const value_partitioning_domain_t &other,
                      const thresholds<number_t> &ts) const override {
    CRAB_LOG("partition",
             crab::outs() << "Widening " << *this << " and " << other << "\n";);
    if (is_bottom() || other.is_top()) {
      CRAB_LOG("partition", crab::outs() << "Result=" << other << "\n";);
      return other;
    } else if (other.is_bottom() || is_top()) {
      CRAB_LOG("partition", crab::outs() << "Result=" << *this << "\n";);
      return *this;
    } else {
      return apply_binary_op(
          other, "Widening w/ thresholds",
          [&ts](const interval_t &i1, const interval_t &i2) {
            return i1.widening_thresholds(i2, ts);
          },
          [&ts](const NumDomain &v1, const NumDomain &v2) {
            return v1.widening_thresholds(v2, ts);
          });
    }
  }

  value_partitioning_domain_t
  operator&&(const value_partitioning_domain_t &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      // the solution here is sub-optimal
      return apply_binary_op(
          other, "Narrowing",
          [](const interval_t &i1, const interval_t &i2) { return i1 && i2; },
          [](const NumDomain &v1, const NumDomain &v2) { return v1 && v2; });
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().assign(x, e);
      }
      if (m_variable && *m_variable == x) {
	update_partitions();
      } 
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().weak_assign(x, e);
      }
      if (m_variable && *m_variable == x) {
	update_partitions();
      } 
    }
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().apply(op, x, y, z);
      }
      if (m_variable && *m_variable == x) {
        update_partitions();
      }
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().apply(op, x, y, k);
      }
      if (m_variable && *m_variable == x) {
        update_partitions();
      }
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().apply(op, dst, src);
      }
      if (m_variable && *m_variable == dst) {
        update_partitions();
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().apply(op, x, y, z);
      }
      if (m_variable && *m_variable == x) {
        update_partitions();
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().apply(op, x, y, k);
      }
      if (m_variable && *m_variable == x) {
        update_partitions();
      }
    }
  }

  
  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom() || csts.is_true()) {
      return;
    }
    if (csts.is_false()) {
      set_to_bottom();
      return;
    }

    if (m_variable) {
      linear_constraint_system_t diseqs, non_diseqs;
      split_relevant_disequalities(csts, *m_variable, diseqs, non_diseqs);

      add_constraints(non_diseqs);
      // Disequalities that use the partition variable will add new
      // partitions
      for (auto it=diseqs.begin(), et=diseqs.end();it!=et;++it) {
	auto &cst = *it;
	linear_constraint_t ltC(cst.expression(),
				linear_constraint_t::STRICT_INEQUALITY);
	linear_constraint_t gtC(-cst.expression(),
			      linear_constraint_t::STRICT_INEQUALITY);
	value_partitioning_domain_t copy(*this);
	add_constraints(ltC);
	copy.add_constraints(gtC);
	operator|=(copy);
      }
    } else {
      add_constraints(csts);
    }
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        if (!partition.get_dom().entails(cst)) {
	  return false;
	}
      }
    }
    return true;
  }
  
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const value_partitioning_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const value_partitioning_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const value_partitioning_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  DEFAULT_SELECT(value_partitioning_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(value_partitioning_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(value_partitioning_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(value_partitioning_domain_t)

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    if (is_bottom()) {
      return;
    }

    if (name == VALUE_PARTITION_START && inputs.size() == 1 &&
        inputs[0].is_variable() && outputs.empty()) {
      variable_t partition_variable = inputs[0].get_variable();
      if (m_variable) {
        // do nothing
      } else {
        if (!partition_variable.get_type().is_integer()) {
          CRAB_ERROR("Partitioning only on integer variables ",
                     partition_variable);
        }
        m_variable = partition_variable;
        update_partitions();
        CRAB_LOG("partition", crab::outs() << "Adding partition variable: "
                                           << partition_variable << "\n"
                                           << *this << "\n";);
      }
    } else if (name == VALUE_PARTITION_END && inputs.size() == 1 &&
               inputs[0].is_variable() && outputs.empty()) {
      variable_t partition_variable = inputs[0].get_variable();
      if (m_variable && *m_variable == partition_variable) {
        remove_partitions();
        CRAB_LOG("partition", crab::outs() << "Removed partition variable: "
                                           << partition_variable << "\n"
                                           << *this << "\n";);
      }
    } else {
      for (auto &partition : m_partitions) {
        partition.get_dom().intrinsic(name, inputs, outputs);
      }
    }
  }

  void
  backward_intrinsic(std::string name,
                     const variable_or_constant_vector_t &inputs,
                     const variable_vector_t &outputs,
                     const value_partitioning_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  void operator-=(const variable_t &v) override {
    if (!is_bottom()) {
      if (m_variable && *m_variable == v) {
        remove_partitions();
        m_partitions[0].get_dom() -= v;
      } else {
        for (auto &partition : m_partitions) {
          partition.get_dom() -= v;
        }
      }
    }
  }

  interval_t operator[](const variable_t &v) override {
    auto smashed = merge_partitions().get_dom();
    return smashed[v];
  }

  interval_t at(const variable_t &v) const override {
    auto smashed = merge_partitions().get_dom();
    return smashed.at(v);
  }

  void normalize() override {}

  void minimize() override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      if (m_variable) {
        for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
          if (m_variable == from[i]) {
            m_variable = to[i];
            break;
          }
        }
      }
      for (auto &partition : m_partitions) {
        partition.get_dom().rename(from, to);
      }
    }
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    if (!is_bottom()) {
      for (auto &partition : m_partitions) {
        partition.get_dom().expand(x, new_x);
      }
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      if (m_variable) {
        if (std::any_of(variables.begin(), variables.end(),
                        [this](const variable_t &v) {
                          return (*(this->m_variable) == v);
                        })) {
          remove_partitions();
        }
      }
      for (auto &partition : m_partitions) {
        partition.get_dom().forget(variables);
      }
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      if (m_variable) {
        if (!std::any_of(variables.begin(), variables.end(),
                         [this](const variable_t &v) {
                           return ((*this->m_variable) == v);
                         })) {
          remove_partitions();
        }
      }
      for (auto &partition : m_partitions) {
        partition.get_dom().project(variables);
      }
    }
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    auto smashed = merge_partitions().get_dom();
    return smashed.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    if (is_bottom()) {
      disjunctive_linear_constraint_system_t res(true);
      return res;
    } else if (is_top()) {
      disjunctive_linear_constraint_system_t res(false);
      return res;
    } else {
      disjunctive_linear_constraint_system_t res(true);
      for (auto &partition : m_partitions) {
        res += partition.get_dom().to_linear_constraint_system();
      }
      return res;
    }
  }

  void write(crab::crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      if (m_variable) {
        o << "{";
        for (unsigned i = 0, sz = m_partitions.size(); i < sz;) {
          o << *m_variable << "=" << m_partitions[i].get_interval() << " => ";
          o << m_partitions[i].get_dom();
          ++i;
          if (i < sz) {
            o << ", ";
          }
        }
        o << "}";
      } else {
        o << m_partitions[0].get_dom();
      }
    }
  }

  std::string domain_name() const override {
    return std::string("ValuePartitioning(") +
           m_partitions[0].get_dom().domain_name() + ")";
  }
};

namespace value_partitioning_domain_impl {
template <typename NumDomain> class partition {
  using number_t = typename NumDomain::number_t;
  using interval_t = ikos::interval<number_t>;

  interval_t m_interval;
  NumDomain m_dom;

public:
  partition(interval_t &&i, NumDomain &&d)
      : m_interval(std::move(i)), m_dom(std::move(d)) {}

  partition(const partition &other) = default;
  partition(partition &&other) = default;  
  partition &operator=(const partition &other) = default;
  partition &operator=(partition &&other) = default;  
    
  static partition top() {    
    NumDomain top_base;
    partition res(interval_t::top(), std::move(top_base));
    return res;
  }

  static partition bottom() {
    NumDomain bot_base;
    bot_base.set_to_bottom();
    partition res(interval_t::bottom(), std::move(bot_base));
    return res;
  }


  void join(const partition &other) {
    m_interval = m_interval | other.m_interval;
    m_dom |= other.m_dom;
  }

  partition join(const partition &other) const {
    return partition(m_interval | other.m_interval, m_dom | other.m_dom);
  }

  void meet(const partition &other) {
    m_interval = m_interval & other.m_interval;
    m_dom &= other.m_dom;
  }

  // join only the interval part of the partition
  void join_interval(const partition &other) {
    m_interval = m_interval | other.get_interval();
  }
  
  NumDomain &get_dom() { return m_dom; }

  const NumDomain &get_dom() const { return m_dom; }

  interval_t &get_interval() { return m_interval; }

  const interval_t &get_interval() const { return m_interval; }
};
} // namespace value_partitioning_domain_impl

template <typename Domain>
struct abstract_domain_traits<value_partitioning_domain<Domain>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

/*
 * Lift value_partitioning_domain to a Cartesian product.  This allows
 * to perform simultaneously value partitioning on multiple
 * variables. The concretization is the intersection of the
 * concretization of each element of the Cartesian product.
 *
 * Note that this domain is equivalent to run a separate
 * value_partitioning_domain analysis on each variable. Thus, this
 * domain is less precise than other disjunctive domain constructions
 * such as the reduced cardinal power.
 */
template <typename NumDomain>
class product_value_partitioning_domain final
    : public abstract_domain_api<product_value_partitioning_domain<NumDomain>> {

public:
  using this_type = product_value_partitioning_domain<NumDomain>;
  using abstract_domain_api_t = abstract_domain_api<this_type>;
  using typename abstract_domain_api_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_t;
  using typename abstract_domain_api_t::linear_expression_t;
  using typename abstract_domain_api_t::number_t;
  using typename abstract_domain_api_t::reference_constraint_t;
  using typename abstract_domain_api_t::variable_or_constant_t;
  using typename abstract_domain_api_t::variable_or_constant_vector_t;
  using typename abstract_domain_api_t::variable_t;
  using typename abstract_domain_api_t::variable_vector_t;
  using typename abstract_domain_api_t::varname_t;
  using interval_t = ikos::interval<number_t>;
  using value_partitioning_domain_t = value_partitioning_domain<NumDomain>;

private:
  NumDomain m_absval;
  std::vector<value_partitioning_domain_t> m_product;

  bool has_partitions() const { return !m_product.empty(); }

  void remove_all_partitions() { m_product.clear(); }

  /** A state is normalized if:
   *
   *  let N be m_product.size()
   *  1. (N == 0)  OR
   *  2. (N >= 1  AND
   *      for all i \in [0,N-1] :: m_product[i].get_variable()  AND
   *      for all i,j such that i<j :: m_product[i].get_variable() <
   *                                   m_product[j].get_variable())
   **/
  void normalize(product_value_partitioning_domain &val) const {
    crab::CrabStats::count(domain_name() + ".count.normalize");
    crab::ScopedCrabStats __st__(domain_name() + ".normalize");

    if (!has_partitions()) {
      return;
    }

    // Remove any element without partitioning variable. This is okay
    // because we have intersection semantics.
    size_t size = val.m_product.size();
    for (auto it = val.m_product.begin(); it != val.m_product.end();) {
      if ((!it->get_variable()) && (size > 1)) {
        it = val.m_product.erase(it);
        --size;
      } else {
        ++it;
      }
    }

    // Sort all components so then binary operations are linear
    if (val.m_product.size() >= 2) {
      std::sort(val.m_product.begin(), val.m_product.end(),
                [](const value_partitioning_domain_t &v1,
                   const value_partitioning_domain_t &v2) {
                  assert(v1.get_variable());
                  assert(v2.get_variable());
                  return *(v1.get_variable()) < *(v2.get_variable());
                });
    } else if (val.m_product.size() == 1) {
      if (!val.m_product[0].get_variable()) {
        val.m_absval = val.m_product[0].merge_partitions().get_dom();
        val.remove_all_partitions();
      }
    }

    #if 0
    // Sanity check
    for (auto const&e: m_product) {
      if (!e.get_variable()) {
	CRAB_ERROR(domain_name(), " normalization did not normalize");
      }
    }
    #endif 
    
  }

  // Sanity check
  void ERROR_IF_NOT_PARTITION_VAR(const value_partitioning_domain_t &e, unsigned line) const {
    if (!e.get_variable()) {
      CRAB_ERROR(domain_name(), " expected partitioning variable at line ", line);
    }
  }

  // Merge the whole Cartesian product into a single
  // value_partitioning_domain_t with partition_variable.
  NumDomain merge_product() const {
    if (is_bottom()) {
      CRAB_ERROR("called merge_partitions() on bottom");
    }

    if (!has_partitions()) {
      return m_absval;
    }

    NumDomain val = m_product[0].merge_partitions().get_dom();
    for (unsigned i = 1, sz = m_product.size(); i < sz; ++i) {
      // We meet here because the concretization of the Cartesian
      // product is the intersection of the concretization of each
      // product element.
      val = val & m_product[i].merge_partitions().get_dom();
    }
    return val;
  }

  value_partitioning_domain_t
  merge_product(const variable_t &partition_variable) const {
    using partition_t = typename value_partitioning_domain_t::partition_t;
    NumDomain val = merge_product();
    partition_t partition(val.at(partition_variable), std::move(val));
    value_partitioning_domain_t p(partition_variable, {partition});
    return p;
  }

  product_value_partitioning_domain(NumDomain &&absval)
      : m_absval(std::move(absval)) {}

  bool contains(const variable_t &var) const {
    if (is_bottom()) {
      CRAB_ERROR("called contains on bottom");
    }

    if (!has_partitions()) {
      return false;
    }

    for (auto const &elem : m_product) {
      boost::optional<variable_t> var_opt = elem.get_variable();
      if (var_opt && (*var_opt) == var) {
        return true;
      }
    }
    return false;
  }

  // return null if var has no partitions
  value_partitioning_domain_t *get_partition(const variable_t &var) {
    if (is_bottom()) {
      CRAB_ERROR("called get_partition(...) on bottom");
    }

    if (!has_partitions()) {
      return nullptr;
    }

    for (auto &elem : m_product) {
      boost::optional<variable_t> var_opt = elem.get_variable();
      if (var_opt && (*var_opt) == var) {
        return &elem;
      }
    }
    return nullptr;
  }
  
  void add_partition(value_partitioning_domain_t &&elem) {
    if (is_bottom()) {
      return;
    }

    if (elem.is_bottom()) {
      set_to_bottom();
      return;
    }

    CRAB_LOG("product-partition",
             crab::outs() << "Adding partition " << elem << " -- ";);

    // it's safe to throw away m_absval
    m_absval.set_to_top();

    size_t size = m_product.size();
    if (size == 0) {
      m_product.emplace_back(std::move(elem));
    } else if (size == 1) {
      ERROR_IF_NOT_PARTITION_VAR(m_product[0], __LINE__);
      if (elem.get_variable()) {
        m_product.emplace_back(std::move(elem));
        if (!(m_product[0].get_variable() < m_product[1].get_variable())) {
          std::swap(m_product[0], m_product[1]);
        }
      } else {
        // do nothing
      }
    } else {
      // At this point, all elements in m_product have partitioning
      // variables and m_product is sorted.
      if (elem.get_variable()) {
        auto lb = std::lower_bound(
            m_product.begin(), m_product.end(), elem,
            [this](const value_partitioning_domain_t &v1,
                   const value_partitioning_domain_t &v2) {
              ERROR_IF_NOT_PARTITION_VAR(v1, __LINE__);
              ERROR_IF_NOT_PARTITION_VAR(v2, __LINE__);
              return *(v1.get_variable()) < *(v2.get_variable());
            });
        if (lb != m_product.end() &&
            !(*(elem.get_variable()) < *((*lb).get_variable()))) {
          // already exists: do nothing
        } else {
          m_product.insert(lb, std::move(elem));
        }
      }
    }
    CRAB_LOG("product-partition", crab::outs() << *this << "\n";);
  }

  void remove_partition(const variable_t &var) {
    if (is_bottom()) {
      return;
    }

    if (!has_partitions()) {
      return;
    }

    for (auto &e : m_product) {
      if (e.get_variable() && (*(e.get_variable()) == var)) {
        e.remove_partitions();
      }
    }

    // REVISIT: we could avoid fully normalization
    normalize();
    CRAB_LOG("product-partition", crab::outs()
                                      << "After removing partitions on " << var
                                      << " -- " << *this << "\n";);
  }

  template <typename T> struct binary_op {
    virtual void apply(T &left, const T &right) const = 0;
  };

  template <typename T> struct join_op : public binary_op<T> {
    void apply(T &left, const T &right) const override { left |= right; }
  };

  template <typename T> struct meet_op : public binary_op<T> {
    void apply(T &left, const T &right) const override { left &= right; }
  };

  template <typename T> struct widen_op : public binary_op<T> {
    void apply(T &left, const T &right) const override { left = left || right; }
  };

  template <typename T> struct widen_with_thresholds_op : public binary_op<T> {
    const thresholds<number_t> &m_ts;
    widen_with_thresholds_op(const thresholds<number_t> &ts) : m_ts(ts) {}
    void apply(T &left, const T &right) const override {
      left = left.widening_thresholds(right, m_ts);
    }
  };

  template <typename T> struct narrow_op : public binary_op<T> {
    void apply(T &left, const T &right) const override { left = left && right; }
  };

  using binary_op_t = binary_op<value_partitioning_domain_t>;
  using binary_base_op_t = binary_op<NumDomain>;
  using join_op_t = join_op<value_partitioning_domain_t>;
  using join_base_op_t = join_op<NumDomain>;
  using meet_op_t = meet_op<value_partitioning_domain_t>;
  using meet_base_op_t = meet_op<NumDomain>;
  using widen_op_t = widen_op<value_partitioning_domain_t>;
  using widen_base_op_t = widen_op<NumDomain>;
  using widen_with_thresholds_op_t =
      widen_with_thresholds_op<value_partitioning_domain_t>;
  using widen_with_thresholds_base_op_t = widen_with_thresholds_op<NumDomain>;
  using narrow_op_t = narrow_op<value_partitioning_domain_t>;
  using narrow_base_op_t = narrow_op<NumDomain>;

  bool is_leq(const this_type &left, const this_type &right) const {
    assert(!left.is_bottom());
    assert(!left.is_top());
    assert(!right.is_bottom());
    assert(!right.is_top());

    if (!left.has_partitions() && !right.has_partitions()) {
      return (left.m_absval <= right.m_absval);
    } else if (!left.has_partitions() && right.has_partitions()) {
      return (left.m_absval <= right.merge_product());
    } else if (left.has_partitions() && !right.has_partitions()) {
      return (left.merge_product() <= right.m_absval);
    }

    // Since the more partitions the more precise, return true
    // if left's partitions are a superset of the right's
    // partitions
    for (auto left_it = left.m_product.begin(), left_et = left.m_product.end(),
              right_it = right.m_product.begin(),
              right_et = right.m_product.end();
         right_it != right_et;) {
      if (left_it == left_et) {
        // we don't have more left partitions but some right
        // partitions
        return false;
      } else {
        const value_partitioning_domain_t &left_val = *left_it;
        const value_partitioning_domain_t &right_val = *right_it;
        ERROR_IF_NOT_PARTITION_VAR(left_val, __LINE__);
        ERROR_IF_NOT_PARTITION_VAR(right_val, __LINE__);

        auto left_var = *(left_val.get_variable());
        auto right_var = *(right_val.get_variable());
        if (left_var < right_var) {
          ++left_it;
        } else if (left_var == right_var) {
          if (!(left_val <= right_val)) {
            return false;
          }
          ++left_it;
          ++right_it;
        } else {
          // some right's partition won't be on the left so we can
          // bail out.
          return false;
        }
      }
    }
    return true;
  }

  // Recall that the more partitions, the more precise is an abstract
  // state. Thus, for the join or widening we do the intersection of
  // the two Cartesian products and apply operation on each compatible
  // (same variable partitioning) pair of elements.
  void join_or_widening(this_type &left, const this_type &right,
                        binary_op_t &op, binary_base_op_t &base_op) const {
    assert(!left.is_bottom());
    assert(!left.is_top());
    assert(!right.is_bottom());
    assert(!right.is_top());

    if (!left.has_partitions() && !right.has_partitions()) {
      base_op.apply(left.m_absval, right.m_absval);
    } else if (!left.has_partitions() && right.has_partitions()) {
      auto smashed_right = right.merge_product();
      base_op.apply(left.m_absval, smashed_right);
    } else if (left.has_partitions() && !right.has_partitions()) {
      auto smashed_left = left.merge_product();
      base_op.apply(smashed_left, right.m_absval);
      left.m_absval = std::move(smashed_left);
      left.remove_all_partitions();
    } else {
      // FIXME: merged_left might not be used and this is an expensive op.
      NumDomain merged_left = left.merge_product();

      auto right_it = right.m_product.begin();
      auto right_et = right.m_product.end();
      for (auto left_it = left.m_product.begin(); left_it != left.m_product.end();) {
        if (right_it == right_et) {
          left.m_product.erase(left_it, left.m_product.end());
          break;
        } else {
          value_partitioning_domain_t &left_val = *left_it;
          const value_partitioning_domain_t &right_val = *right_it;
          ERROR_IF_NOT_PARTITION_VAR(left_val, __LINE__);
          ERROR_IF_NOT_PARTITION_VAR(right_val, __LINE__);

          auto left_var = *(left_val.get_variable());
          auto right_var = *(right_val.get_variable());

          if (left_var < right_var) {
            left_it = left.m_product.erase(left_it);
          } else if (left_var == right_var) {
            op.apply(left_val, right_val);
            ++left_it;
            ++right_it;
          } else {
            ++right_it;
          }
        }
      }
      // At this point left can be empty if the left and right do not
      // have common partitions
      if (left.m_product.empty()) {
        NumDomain merged_right = right.merge_product();
        base_op.apply(merged_left, merged_right);
        left.m_absval = std::move(merged_left);
      }
    }
  }

  // Recall that the more partitions, the more precise is an abstract
  // state.  Thus, for the meet or narrowing we do the union of all
  // elements of the two Cartesian products and apply operation on
  // each compatible (same variable partitioning) pair of elements.
  void meet_or_narrowing(this_type &left, const this_type &right,
                         binary_op_t &op, binary_base_op_t &base_op) const {
    assert(!left.is_bottom());
    assert(!left.is_top());
    assert(!right.is_bottom());
    assert(!right.is_top());

    if (!left.has_partitions() && !right.has_partitions()) {
      base_op.apply(left.m_absval, right.m_absval);
    } else if (!left.has_partitions() && right.has_partitions()) {
      auto smashed_right = right.merge_product();
      base_op.apply(left.m_absval, smashed_right);
    } else if (left.has_partitions() && !right.has_partitions()) {
      auto smashed_left = left.merge_product();
      base_op.apply(smashed_left, right.m_absval);
      left.m_absval = std::move(smashed_left);
      left.remove_all_partitions();
    } else {
      auto right_it = right.m_product.begin();
      auto right_et = right.m_product.end();
      for (auto left_it = left.m_product.begin(); left_it != left.m_product.end();) {
        if (right_it == right_et) {
          return;
        } else {
          value_partitioning_domain_t &left_val = *left_it;
          const value_partitioning_domain_t &right_val = *right_it;
          ERROR_IF_NOT_PARTITION_VAR(left_val, __LINE__);
          ERROR_IF_NOT_PARTITION_VAR(right_val, __LINE__);

          auto left_var = *(left_val.get_variable());
          auto right_var = *(right_val.get_variable());

          if (left_var < right_var) {
            ++left_it;
          } else if (left_var == right_var) {
            op.apply(left_val, right_val);
            if (left_val.is_bottom()) {
              left.set_to_bottom();
              return;
            }
            ++left_it;
            ++right_it;
          } else {
            left_it = (left.m_product.insert(left_it, right_val))++;
            ++right_it;
          }
        }
      }
      // append the rest of right's partitions on the left.
      left.m_product.insert(left.m_product.end(), right_it, right_et);
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e, bool weak) {
    if (!is_bottom()) {
      if (!has_partitions()) {
	if (!weak) {
	  m_absval.assign(x, e);
	} else {
	  m_absval.weak_assign(x, e);
	} 
      } else {
	
	if (boost::optional<variable_t> y = e.get_variable()) {
	  // x := y
	  // 
	  if (value_partitioning_domain_t *partition_y = get_partition(*y)) {
	    // We create a new partition for x if we have one for y
	    if (!contains(x)) {
	      // The user didn't select x as a partitioning variable.
	      // This step is a bit problematic. If we don't add a new
	      // partition for x then we might lose precision if x
	      // flows later to y which can happen in loops. If we add
	      // a new partition then the analysis might blow up if
	      // too many partitions are added.
	      CRAB_WARN(domain_name(),
			" assigning a partitioning variable ", *y,
			" to a non-partitioning one ", x);
	    }  
	    value_partitioning_domain_t partition_x(*partition_y);
	    partition_x.m_variable = x;
	    remove_partition(x); 
	    add_partition(std::move(partition_x));
	    // TODO(IMPORTANT): if x is not partitioning variable
	    // chosen by the user then we should remove the partition
	    // for x when we remove the partition of y. This is not
	    // done at the moment.
	  }
	}
	
        for (auto &elem : m_product) {
	  if (!weak) {
	    elem.assign(x, e);
	  } else {
	    elem.weak_assign(x, e);
	  } 
        }
      }
    }
  }
  
public:
  product_value_partitioning_domain() {
    NumDomain top;
    m_absval = std::move(top);
  }

  product_value_partitioning_domain(const this_type &other) = default;

  // product_value_partitioning_domain(const this_type &other):
  //   m_absval(other.m_absval),
  //   m_product(other.m_product) {
  //   crab::CrabStats::count(domain_name() + ".count.copy");
  //   crab::ScopedCrabStats __st__(domain_name() + ".copy");
  // }

  this_type &operator=(const this_type &other) = default;

  // this_type &operator=(const this_type &other) {
  //   crab::CrabStats::count(domain_name() + ".count.copy");
  //   crab::ScopedCrabStats __st__(domain_name() + ".copy");
  //   if (this != &other) {
  //     m_absval = other.m_absval;
  //     m_product = other.m_product;
  //   }
  //   return *this;
  // }

  product_value_partitioning_domain(this_type &&other) = default;
  this_type &operator=(this_type &&other) = default;

  bool is_bottom() const override {
    if (!has_partitions()) {
      return m_absval.is_bottom();
    }

    for (auto const &elem : m_product) {
      if (elem.is_bottom()) {
        return true;
      }
    }
    return false;
  }

  bool is_top() const override {
    if (!has_partitions()) {
      return m_absval.is_top();
    }

    for (auto const &elem : m_product) {
      if (!elem.is_top()) {
        return false;
      }
    }
    return true;
  }

  this_type make_bottom() const override {
    this_type res;
    res.set_to_bottom();
    return res;
  }

  this_type make_top() const override {
    this_type res;
    return res;
  }

  void set_to_top() override {
    if (!has_partitions()) {
      m_absval.set_to_top();
    } else {
      remove_all_partitions();
      NumDomain top;
      m_absval = std::move(top);
    }
  }

  void set_to_bottom() override {
    if (!has_partitions()) {
      m_absval.set_to_bottom();
    } else {
      remove_all_partitions();
      NumDomain bot;
      bot.set_to_bottom();
      m_absval = std::move(bot);
    }
  }

  bool operator<=(const this_type &other) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom() || other.is_top()) {
      return true;
    } else if (other.is_bottom() || is_top()) {
      return false;
    } else {
      return is_leq(*this, other);
    }
  }

  void operator|=(const this_type &other) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        m_absval |= other.m_absval;
      } else {
        join_op_t op;
        join_base_op_t base_op;
        join_or_widening(*this, other, op, base_op);
      }
    }
  }

  this_type operator|(const this_type &other) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return this_type(m_absval | other.m_absval);
      } else {
        join_op_t op;
        join_base_op_t base_op;
        this_type left(*this);
        join_or_widening(left, other, op, base_op);
        return left;
      }
    }
  }

  this_type operator&(const this_type &other) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return this_type(m_absval & other.m_absval);
      } else {
        meet_op_t op;
        meet_base_op_t base_op;
        this_type left(*this);
        meet_or_narrowing(left, other, op, base_op);
        return left;
      }
    }
  }

  void operator&=(const this_type &other) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (is_top() || other.is_bottom()) {
      *this = other;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return m_absval &= other.m_absval;
      } else {
        meet_op_t op;
        meet_base_op_t base_op;
        meet_or_narrowing(*this, other, op, base_op);
      }
    }
  }

  this_type operator||(const this_type &other) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return this_type(m_absval || other.m_absval);
      } else {
        widen_op_t op;
        widen_base_op_t base_op;
        this_type left(*this);
        join_or_widening(left, other, op, base_op);
        return left;
      }
    }
  }

  this_type widening_thresholds(const this_type &other,
                                const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return this_type(m_absval.widening_thresholds(other.m_absval, ts));
      } else {
        widen_with_thresholds_op_t op(ts);
        widen_with_thresholds_base_op_t base_op(ts);
        this_type left(*this);
        join_or_widening(left, other, op, base_op);
        return left;
      }
    }
  }

  this_type operator&&(const this_type &other) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (is_top() || other.is_bottom()) {
      return other;
    } else {
      if (!has_partitions() && !other.has_partitions()) {
        return this_type(m_absval && other.m_absval);
      } else {
        narrow_op_t op;
        narrow_base_op_t base_op;
        this_type left(*this);
        meet_or_narrowing(left, other, op, base_op);
        return left;
      }
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    assign(x, e, false /*!weak*/);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    assign(x, e, true /*weak*/);
  }
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.apply(op, x, y, z);
      } else {
        for (auto &elem : m_product) {
          elem.apply(op, x, y, z);
        }
      }
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.apply(op, x, y, k);
      } else {
        for (auto &elem : m_product) {
          elem.apply(op, x, y, k);
        }
      }
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.apply(op, dst, src);
      } else {
        for (auto &elem : m_product) {
          elem.apply(op, dst, src);
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.apply(op, x, y, z);
      } else {
        for (auto &elem : m_product) {
          elem.apply(op, x, y, z);
        }
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.apply(op, x, y, k);
      } else {
        for (auto &elem : m_product) {
          elem.apply(op, x, y, k);
        }
      }
    }
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom() || csts.is_true()) {
      return;
    }

    if (csts.is_false()) {
      set_to_bottom();
      return;
    }

    if (!has_partitions()) {
      m_absval += csts;
    } else {
      for (auto &elem : m_product) {
        elem += csts;
        if (elem.is_bottom()) {
          set_to_bottom();
          return;
        }
      }
    }
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (!is_bottom()) {
      if (!has_partitions()) {
	return m_absval.entails(cst);
      } else {
	for (auto &elem: m_product) {
	  if (!elem.entails(cst)) {
	    return false;
	  }
	}
      }
    }
    return true;
  }
  
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const this_type &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  DEFAULT_SELECT(this_type)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(this_type)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(this_type)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(this_type)

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    if (is_bottom()) {
      return;
    }

    if (name == VALUE_PARTITION_START && outputs.empty()) {
      for (auto const &input : inputs) {
        if (input.is_variable()) {
          variable_t partition_var = input.get_variable();
          if (!contains(partition_var)) {
            add_partition(merge_product(partition_var));
          }
        } else {
          CRAB_ERROR(VALUE_PARTITION_START, " only expects variables");
        }
      }
    } else if (name == VALUE_PARTITION_END && outputs.empty()) {
      for (auto const &input : inputs) {
        if (input.is_variable()) {
          variable_t partition_var = input.get_variable();
          remove_partition(partition_var);
        } else {
          CRAB_ERROR(VALUE_PARTITION_END, " only expects variables");
        }
      }
    } else {
      for (auto &elem : m_product) {
        elem.intrinsic(name, inputs, outputs);
      }
    }
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_type &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  void operator-=(const variable_t &v) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval -= v;
      } else {
        for (auto &elem : m_product) {
          elem -= v;
        }
        normalize();
      }
    }
  }

  interval_t operator[](const variable_t &v) override {
    normalize();
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      if (!has_partitions()) {
        return m_absval[v];
      } else {
        for (auto &e : m_product) {
          if (!e.get_variable()) {
            return e[v];
          }
        }
        CRAB_ERROR("operator[] unreachable");
      }
    }
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      if (!has_partitions()) {
        return m_absval.at(v);
      } else {
        for (auto const &e : m_product) {
          if (!e.get_variable()) {
            return e.at(v);
          }
        }
        CRAB_ERROR("at unreachable");
      }
    }
  }

  void normalize() override { normalize(*this); }
  void minimize() override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.rename(from, to);
      } else {
        for (auto &elem : m_product) {
          elem.rename(from, to);
        }
	normalize();	
      }
    }
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.expand(x, new_x);
      } else {
        for (auto &elem : m_product) {
          elem.expand(x, new_x);
        }
        normalize();
      }
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.forget(variables);
      } else {
        for (auto &elem : m_product) {
          elem.forget(variables);
        }
        normalize();
      }
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      if (!has_partitions()) {
        m_absval.project(variables);
      } else {
        for (auto &elem : m_product) {
          elem.project(variables);
        }
        normalize();
      }
    }
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      linear_constraint_system_t out;
      out += linear_constraint_t::get_false();
      return out;
    } else if (is_top()) {
      linear_constraint_system_t out;
      out += linear_constraint_t::get_true();
      return out;
    } else {
      if (!has_partitions()) {
        return m_absval.to_linear_constraint_system();
      } else {
        for (auto const &e : m_product) {
          if (!e.get_variable()) {
            return e.to_linear_constraint_system();
          }
        }
        CRAB_ERROR("to_linear_constraint_system unreachable");
      }
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    if (is_bottom()) {
      disjunctive_linear_constraint_system_t res(true);
      return res;
    } else if (is_top()) {
      disjunctive_linear_constraint_system_t res(false);
      return res;
    } else {
      if (!has_partitions()) {
        return m_absval.to_disjunctive_linear_constraint_system();
      } else {
        CRAB_WARN("TODO ", domain_name(),
                  "::to_disjunctive_linear_constraint_system");
        disjunctive_linear_constraint_system_t res(false);
        return res;
      }
    }
  }

  void write(crab::crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      if (!has_partitions()) {
        m_absval.write(o);
      } else {
        o << "(";
        for (unsigned i = 0, sz = m_product.size(); i < sz;) {
          o << m_product[i];
          ++i;
          if (i < sz) {
            o << ", ";
          }
        }
        o << ")";
      }
    }
  }

  std::string domain_name() const override {
    value_partitioning_domain_t top;
    return top.domain_name();
  }
};

template <typename Domain>
struct abstract_domain_traits<product_value_partitioning_domain<Domain>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

} // namespace domains
} // namespace crab
