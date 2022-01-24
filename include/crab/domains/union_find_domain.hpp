#pragma once

#include <crab/domains/discrete_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>

/**
 * A standard union-find equipped with lattice operations (inclusion,
 * lub, meet, etc).
 **/

namespace crab {
namespace domains {

template <class Variable, class Domain> class equivalence_class {
public:
  using variable_t = Variable;
  using domain_t = Domain;

private:
  std::size_t m_rank;
  Domain m_val;

public:
  explicit equivalence_class(Domain val) : m_rank(0), m_val(val) {}

  std::size_t &get_rank() { return m_rank; }

  std::size_t get_rank() const { return m_rank; }

  Domain &get_domain() { return m_val; }

  const Domain &get_domain() const { return m_val; }

}; // end class equivalence_class

/*
  Partition a set of variables into equivalence classes. Apart from
  lattice operations, the domain provides standard union-find
  operations:

  - make(v)
  - union(v1, v2)
  - find(v)

  In addition, each equivalence class has associated an abstract state
  of type Domain.

  The method "top()" returns an union-find in which all variables are
  in the same equivalence class. Since the domain does not know the
  universe of variables, the method "is_top()" only succeeds if the
  union_find_domain state was created by calling directly "top()".
  The method "bottom()" represents an unreachable/failure abstract
  state. Do not use "top()" if you intend to generate an empty
  union-find. Instead, use the constructor without parameters to
  produce an empty union-find.
*/
template <class Variable, class Domain> class union_find_domain {
public:
  using variable_t = Variable;
  using domain_t = Domain;
  enum class lattice_val { bottom, top, neither_top_nor_bot };

private:
  using union_find_domain_t = union_find_domain<Variable, Domain>;
  using variable_set_t = ikos::discrete_domain<variable_t>;
  using parents_map_t = std::unordered_map<variable_t, variable_t>;
  using equivalence_class_t = equivalence_class<variable_t, Domain>;
  using equiv_class_vars_t = std::unordered_map<variable_t, variable_set_t>;
  using classes_map_t = std::unordered_map<variable_t, equivalence_class_t>;

  // immediate link to the parent
  parents_map_t m_parents;
  // each root has its equivalence class
  classes_map_t m_classes;
  lattice_val m_val;

  void clear() {
    m_parents.clear();
    m_classes.clear();
  }

  // Build a map from representative to a set with all the variables
  // in the equivalence class.
  equiv_class_vars_t equiv_classes_vars() {
    equiv_class_vars_t res;
    for (auto &kv : m_parents) {
      variable_t rep = find(kv.second);
      auto it = res.find(rep);
      if (it == res.end()) {
        variable_set_t varset = variable_set_t::bottom();
        varset += kv.first;
        res.insert({rep, varset});
      } else {
        variable_set_t &varset = it->second;
        varset += kv.first;
      }
    }
    return res;
  }

  // Return the representative after joining all variables in vars.
  // Pre-condition: if dom != none then forall v \in vars:: contains(v)
  boost::optional<variable_t>
  merge_vars(const variable_set_t &vars,
             boost::optional<domain_t> dom = boost::none) {
    auto it = vars.begin();
    auto et = vars.end();
    if (it == et) {
      return boost::none;
    }
    variable_t v = *it;
    if (dom && !contains(v)) {
      make(v, *dom);
    }
    ++it;

    if (it == et) {
      return (contains(v) ? boost::optional<variable_t>(find(v)) : boost::none);
    }

    boost::optional<variable_t> repr;
    for (; it != et; ++it) {
      if (dom && !contains(*it)) {
        make(*it, *dom);
      }
      repr = join(v, *it);
    }
    return repr;
  }

  void print(crab_os &o) const {
    o << "({";
    for (auto it = m_parents.begin(), et = m_parents.end(); it != et;) {
      variable_t key = it->first, value = it->second;
      o << key << " -> " << value;
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "}, {";
    for (auto it = m_classes.begin(), et = m_classes.end(); it != et;) {
      variable_t rep = it->first;
      equivalence_class_t ec = it->second;
      o << rep << " -> " << ec.get_domain();
      ++it;
      if (it != et) {
        o << ", ";
      }
    }
    o << "})";
  }

  void get_all_variables(variable_set_t &out) const {
    for (auto &kv: m_parents) {
      out += kv.first;
    }
  }

  union_find_domain_t join_or_widening(const union_find_domain_t &o, bool is_join) const {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() | o.is_top()) {
      union_find_domain_t res;
      return res;
    } else {
      CRAB_LOG("union-find",
               crab::outs() << (is_join ? "Join" : "Widen")
	                    << "\n\t" << *this << "\n\t" << o << "\n";);
      union_find_domain_t left(*this);
      union_find_domain_t right(o);

      // Keep only common variables      
      variable_set_t left_vars = variable_set_t::bottom();
      left.get_all_variables(left_vars);
      for (auto v: left_vars) {
	if (!right.contains(v)) {
	  left.forget(v);
	}
      }
      variable_set_t right_vars = variable_set_t::bottom();      
      right.get_all_variables(right_vars);      
      for (auto v: right_vars) {
	if (!left.contains(v)) {
	  right.forget(v);
	}
      }
      
      equiv_class_vars_t right_equiv_classes = right.equiv_classes_vars();
      equiv_class_vars_t left_equiv_classes = left.equiv_classes_vars();
      
      // Merge equivalence classes from the left to the right while
      // joining/widening the domains associated to the equivalence
      // classes.
      //
      // The merging on the right is needed so that right_dom is updated.
      for (auto &kv : left_equiv_classes) {
        boost::optional<variable_t> right_repr = right.merge_vars(kv.second);
        if (!right_repr)
          continue; // this shouldn't happen
        domain_t &left_dom = left.get_equiv_class(kv.first).get_domain();
        const domain_t &right_dom =
            right.get_equiv_class(*right_repr).get_domain();
        left_dom = (is_join ? (left_dom | right_dom): (left_dom || right_dom));
      }

      // Merge equivalence classes from the right to the left while
      // joining/widening the domains associated to the equivalence
      // classes.
      for (auto &kv : right_equiv_classes) {
        boost::optional<variable_t> left_repr = left.merge_vars(kv.second);
        if (!left_repr)
          continue; // this shouldn't happen
        domain_t &left_dom = left.get_equiv_class(*left_repr).get_domain();
        const domain_t &right_dom =
            right.get_equiv_class(kv.first).get_domain();
	left_dom = (is_join ? (left_dom | right_dom): (left_dom || right_dom));	
      }
      CRAB_LOG("union-find", crab::outs() << "Res=" << left << "\n";);
      return left;
    }
  }

  union_find_domain_t meet_or_narrowing(const union_find_domain_t &o, bool is_meet) const {
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    } else {
      /* TOIMPROVE: this operation is imprecise because we still merge
       * the equivalence classes.
       */
      union_find_domain_t left(*this);
      union_find_domain_t right(o);
      equiv_class_vars_t right_equiv_classes = right.equiv_classes_vars();
      equiv_class_vars_t left_equiv_classes = left.equiv_classes_vars();

      // Merge equivalence classes from the right to the left while
      // applying meet/narrowing
      //
      // The merging on the right is needed so that right_dom is updated.
      for (auto &kv : left_equiv_classes) {
        domain_t &left_dom = left.get_equiv_class(kv.first).get_domain();
        // add variables on the right if they do not exist
        boost::optional<variable_t> right_repr =
            right.merge_vars(kv.second, left_dom);
        if (!right_repr)
          continue; // this shouldn't happen
        const domain_t &right_dom =
            right.get_equiv_class(*right_repr).get_domain();
        left_dom = (is_meet ? left_dom & right_dom: left_dom && right_dom);
        if (left_dom.is_bottom()) {
          left.set_to_bottom();
          return left;
        }
      }

      // Merge equivalence classes from the right to the left while
      // applying meet/narrowing
      for (auto &kv : right_equiv_classes) {
        const domain_t &right_dom =
            right.get_equiv_class(kv.first).get_domain();
        boost::optional<variable_t> left_repr =
            left.merge_vars(kv.second, right_dom);
        if (!left_repr)
          continue; // this shouldn't happen
        domain_t &left_dom = left.get_equiv_class(*left_repr).get_domain();
	left_dom = (is_meet ? left_dom & right_dom: left_dom && right_dom);
        if (left_dom.is_bottom()) {
          left.set_to_bottom();
          return left;
        }
      }
      return left;
    }
  }
  
public:
  union_find_domain(lattice_val val = lattice_val::neither_top_nor_bot)
      : m_val(val) {}
  union_find_domain(const union_find_domain_t &o) = default;
  union_find_domain(union_find_domain_t &&o) = default;
  union_find_domain_t &operator=(const union_find_domain_t &o) = default;
  union_find_domain_t &operator=(union_find_domain_t &&o) = default;

  /* =============================================================
   *    Standard union-find operations: make, find, and union
   * =============================================================
   */

  // Pre-condition: !contains(v)
  void make(const variable_t &v, Domain val) {
    if (is_bottom()) {
      CRAB_ERROR("calling union_find_domain::make on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("calling union_find_domain::make on top");
    }
    if (contains(v)) {
      CRAB_ERROR("variable already exists when called union_find_domain::make");
    }
    m_parents.insert({v, v});
    m_classes.insert({v, equivalence_class_t(val)});
  }

  // Pre-condition: contains(v)
  // NOTE: it is not "const" because it does path-compression
  variable_t &find(const variable_t &v) {
    if (!contains(v)) {
      assert(false);
      CRAB_ERROR("called union_find_domain::find on a non-existing variable ",
                 v, " in ", *this);
    }
    variable_t &parent = m_parents.at(v);
    if (parent == v) {
      return parent;
    } else {
      return parent = find(parent);
    }
  }

  // Pre-condition: contains(v)
  // find operation without path-compression
  variable_t &find(const variable_t &v) const {
    if (!contains(v)) {
      assert(false);
      CRAB_ERROR("called union_find_domain::find on a non-existing variable ",
                 v, " in ", *this);
    }
    variable_t &parent = m_parents.at(v);
    if (parent == v) {
      return parent;
    } else {
      return find(parent);
    }
  }  

  // Return the representative after merging x and y in the same
  // equivalence class.
  // Pre-condition: contains(x) && contains(y)
  variable_t join(const variable_t &x, const variable_t &y) {
    if (!contains(x)) {
      CRAB_ERROR("called union_find_domain::join on a non-existing variable ",
                 x, " in ", *this);
    }
    if (!contains(y)) {
      CRAB_ERROR("called union_find_domain::join on a non-existing variable ",
                 y, " in ", *this);
    }

    variable_t rep_x = find(x);
    variable_t rep_y = find(y);
    if (rep_x != rep_y) {
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      equivalence_class_t &ec_y = m_classes.at(rep_y);
      domain_t &dom_x = ec_x.get_domain();
      domain_t &dom_y = ec_y.get_domain();
      if (ec_x.get_rank() > ec_y.get_rank()) {
        dom_x = dom_x | dom_y;
        m_parents.at(rep_y) = rep_x;
        m_classes.erase(rep_y);
        return rep_x;
      } else {
        dom_y = dom_y | dom_x;
        m_parents.at(rep_x) = rep_y;
        if (ec_x.get_rank() == ec_y.get_rank()) {
          ec_y.get_rank()++;
        }
        m_classes.erase(rep_x);
        return rep_y;
      }
    } else {
      return rep_x;
    }
  }

  bool contains(const variable_t &v) const {
    return m_parents.find(v) != m_parents.end();
  }

  // Pre-condition: contains(v)
  equivalence_class_t &get_equiv_class(const variable_t &v) {
    return m_classes.at(find(v));
  }

  void remove_equiv_class(const variable_t &v) {
    if (!contains(v)) {
      return;
    }
    variable_t rep_v = find(v);
    m_classes.erase(rep_v);

    equiv_class_vars_t map = equiv_classes_vars();
    auto it = map.find(rep_v);
    if (it != map.end()) {
      for (auto v : it->second) {
        m_parents.erase(v);
      }
    }
  }

  /* =============================================================
   *                     Abstract domain API
   * =============================================================
   */

  static union_find_domain_t top() {
    return union_find_domain_t(lattice_val::top);
  }

  static union_find_domain_t bottom() {
    return union_find_domain_t(lattice_val::bottom);
  }

  bool is_bottom() const { return (m_val == lattice_val::bottom); }

  bool is_top() const { return (m_val == lattice_val::top); }

  void set_to_top() {
    clear();
    m_val = lattice_val::top;
  }

  void set_to_bottom() {
    clear();
    m_val = lattice_val::bottom;
  }

  void set(const variable_t &x, domain_t dom) {
    if (is_bottom() || is_top()) {
      return;
    }

    if (!contains(x)) {
      // Create a singleton equivalence class
      make(x, dom);
    } else {
      // Modify the abstract state of the whole equivalence class
      variable_t rep_x = find(x);
      equivalence_class_t &ec_x = m_classes.at(rep_x);
      ec_x.get_domain() = dom;
    }
  }

  // Pre-condition: contains(x)
  domain_t &get(const variable_t &x) {
    if (is_bottom()) {
      CRAB_ERROR("called union_find_domain::operator[] on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("called union_find_domain::operator[] on top");
    }
    if (!contains(x)) {
      CRAB_ERROR("union_find_domain::operator[] on a non-existing variable ", x,
                 " in ", *this);
    }

    variable_t rep_x = find(x);
    equivalence_class_t &ec_x = m_classes.at(rep_x);
    return ec_x.get_domain();
  }

  // Pre-condition: contains(x)
  const domain_t &get(const variable_t &x) const {
    if (is_bottom()) {
      CRAB_ERROR("called union_find_domain::operator[] on bottom");
    }
    if (is_top()) {
      CRAB_ERROR("called union_find_domain::operator[] on top");
    }
    if (!contains(x)) {
      CRAB_ERROR("union_find_domain::operator[] on a non-existing variable ", x,
                 " in ", *this);
    }

    variable_t rep_x = find(x);
    const equivalence_class_t &ec_x = m_classes.at(rep_x);
    return ec_x.get_domain();
  }  

  // Return true if *this is a refined partitioning of o
  //    forall x,y \in vars(o) :: o.find(x) != o.find(y) =>
  //                              this.find(x) != this.find(y)
  bool operator<=(const union_find_domain_t &o) const {
    if (is_bottom() || o.is_top()) {
      return true;
    } else if (is_top() || o.is_bottom()) {
      return false;
    } else {
      union_find_domain_t left(*this);
      union_find_domain_t right(o);
      equiv_class_vars_t right_equiv_classes = right.equiv_classes_vars();
      equiv_class_vars_t left_equiv_classes = left.equiv_classes_vars();
      for (auto &kv : right_equiv_classes) {
        variable_t right_repr = kv.first;
        variable_set_t &right_vars = kv.second;
        for (variable_t right_v :
             boost::make_iterator_range(right_vars.begin(), right_vars.end())) {
          if (!left.contains(right_v)) {
            return false;
          }
          variable_t left_repr = left.find(right_v);
          variable_set_t &left_vars = left_equiv_classes[left_repr];
          if (!(left_vars <= right_vars)) {
            return false;
          }
          const domain_t &left_dom =
              left.get_equiv_class(left_repr).get_domain();
          const domain_t &right_dom =
              right.get_equiv_class(right_repr).get_domain();
          if (!(left_dom <= right_dom)) {
            return false;
          }
        }
      }
      return true;
    }
  }

  union_find_domain_t operator|(const union_find_domain_t &o) const {
    return join_or_widening(o, true /*is_join*/);
  }

  union_find_domain_t operator&(const union_find_domain_t &o) const {
    return meet_or_narrowing(o, true /*is_meet*/);
  }

  union_find_domain_t operator||(const union_find_domain_t &o) const {
    return join_or_widening(o, false /*is_join*/);    
  }
  
  union_find_domain_t operator&&(const union_find_domain_t &o) const {
    return meet_or_narrowing(o, false /*is_meet*/);
  }

  // Add y into the equivalence class of x
  void add(const variable_t &x, const variable_t &y) {
    if (is_bottom() || is_top()) {
      return;
    }
    if (!contains(x)) {
      return;
    }
    if (contains(y)) {
      forget(y);
    }
    variable_t rep_x = find(x);
    m_parents.insert({y, rep_x});
  }

  void forget(const variable_t &v) {
    if (is_bottom() || is_top()) {
      return;
    }

    if (contains(v)) {
      if (m_parents.at(v) != v) {
        // v is not a representative
        variable_t rep_v = find(v);
        for (auto &kv : m_parents) {
          if (kv.second == v) {
            m_parents.at(kv.first) = rep_v;
          }
        }
      } else {
        // v is the representative of the equivalence class
        boost::optional<variable_t> new_rep;
        for (auto &kv : m_parents) {
          if (kv.first != v && kv.second == v) {
            if (!new_rep) {
              // choose the first one as representative
              new_rep = kv.first;
              m_classes.insert({*new_rep, m_classes.at(v)});
            }
            m_parents.at(kv.first) = *new_rep;
          }
        }
        m_classes.erase(v);
      }
      m_parents.erase(v);
    }
  }

  void project(const std::vector<variable_t> &variables) {
    if (is_bottom() || is_top()) {
      return;
    }

    // First, collect all variables to be forgotten
    std::vector<variable_t> sorted_variables(variables);
    std::sort(sorted_variables.begin(), sorted_variables.end());
    std::vector<variable_t> variables_to_forget;
    variables_to_forget.reserve(m_parents.size());
    for (auto &kv : m_parents) {
      auto lower = std::lower_bound(sorted_variables.begin(),
                                    sorted_variables.end(), kv.first);
      if (lower == sorted_variables.end() || kv.first < *lower) { // not found
        variables_to_forget.push_back(kv.first);
      }
    }
    // Forget variables
    for (auto &v : variables_to_forget) {
      forget(v);
    }
  }

  void rename(const std::vector<variable_t> &old_variables,
              const std::vector<variable_t> &new_variables) {
    if (is_top() || is_bottom()) {
      return;
    }
    if (old_variables.size() != new_variables.size()) {
      CRAB_ERROR(
          "union_find_domain::rename with input vectors of different sizes");
    }

    auto rename_var = [this](const variable_t &x, const variable_t &y) {
      if (contains(y)) {
        CRAB_ERROR("union_find_domain::rename assumes that ", y,
                   " does not exist");
      }
      if (contains(x)) {
        // tuples where x is the left operand
        std::vector<typename parents_map_t::value_type> left_x;
        for (auto &kv : m_parents) {
          if (kv.first == x) {
            left_x.push_back(kv);
          } else if (kv.second == x) {
            kv.second = y;
          }
        }
        for (auto &kv : left_x) {
          m_parents.erase(kv.first);
          m_parents.insert({y, (kv.second == x ? y : kv.second)});
        }
        auto it = m_classes.find(x);
        if (it != m_classes.end()) {
          equivalence_class_t equiv_class = it->second;
          m_classes.erase(it);
          m_classes.insert({y, std::move(equiv_class)});
        }
      }
    };

    for (unsigned i = 0, sz = old_variables.size(); i < sz; ++i) {
      rename_var(old_variables[i], new_variables[i]);
    }
  }

  void write(crab_os &o) const {
    if (is_top()) {
      o << "{}";
    } else if (is_bottom()) {
      o << "_|_";
    } else {

      // Sort everything to print always the same thing
      auto sorted_equiv_classes = [](const equiv_class_vars_t &unsorted_map) {
	     std::vector<std::pair<variable_t, std::vector<variable_t>>> sorted_eq_classes;
	     for (auto &kv: unsorted_map) {
	       std::vector<variable_t> sorted_eq_class(kv.second.begin(), kv.second.end());
	       std::sort(sorted_eq_class.begin(), sorted_eq_class.end());
	       sorted_eq_classes.emplace_back(std::make_pair(kv.first, sorted_eq_class));
	     }
	     std::sort(sorted_eq_classes.begin(), sorted_eq_classes.end(),
		       [](const std::pair<variable_t, std::vector<variable_t>> &p1,
			  const std::pair<variable_t, std::vector<variable_t>> &p2) {
			 return p1.first < p2.first;
		       });
	     return sorted_eq_classes;
      };

      auto print_vars_vector = [&o](const std::vector<variable_t>& vars) {
				 o << "{";
				 for (auto it = vars.begin(), et = vars.end(); it!=et;) {
				   o << *it;
				   ++it;
				   if (it != et) {
				     o << ",";
				   }
				 }
				 o << "}";
			       };
      
      union_find_domain_t tmp(*this);
      CRAB_LOG("union-find-print", tmp.print(o); o << "\n";);
      equiv_class_vars_t ec_vars = tmp.equiv_classes_vars();
      auto sorted_ec_vars = sorted_equiv_classes(ec_vars);
      o << "{";
      for (auto it = sorted_ec_vars.begin(), et = sorted_ec_vars.end(); it != et;) {
	auto &p = *it;
        variable_t &rep = tmp.find(p.first);
	print_vars_vector(p.second);
	o << "=>" << tmp.m_classes.at(rep).get_domain();
        ++it;
        if (it != et) {
          o << ",";
        }
      }
      o << "}";
    }
  }

  friend class crab::crab_os &operator<<(crab::crab_os &o,
                                         const union_find_domain_t &dom) {
    dom.write(o);
    return o;
  }
}; // end class union_find_domain

} // end namespace domains
} // end namespace crab
