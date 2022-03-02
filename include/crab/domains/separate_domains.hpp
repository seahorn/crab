/*******************************************************************************
 *
 * Generic implementation of non-relational domains.
 *
 * Author: Arnaud J. Venet (arnaud.j.venet@nasa.gov)
 *
 * Notices:
 *
 * Copyright (c) 2011 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 * Disclaimers:
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
 * ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
 * TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS,
 * ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 * OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE
 * ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
 * THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
 * ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
 * RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
 * RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
 * DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE,
 * IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST
 * THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL
 * AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS
 * IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH
 * USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM,
 * RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
 * AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
 * RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE,
 * UNILATERAL TERMINATION OF THIS AGREEMENT.
 *
 ******************************************************************************/

#pragma once

#include <crab/domains/patricia_trees.hpp>
#include  <crab/domains/discrete_domains.hpp>
#include <crab/support/debug.hpp>

#include <boost/optional.hpp>

#include <algorithm>
#include <set>
#include <vector>

namespace ikos {

/* Environment from Key to Value with all lattice operations */ 
template <typename Key, typename Value,
	  typename ValueEqual = std::equal_to<Value>>
class separate_domain {

private:
  using patricia_tree_t = patricia_tree<Key, Value, ValueEqual>;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

public:
  using unary_op_t = typename patricia_tree_t::unary_op_t;  
  using separate_domain_t = separate_domain<Key, Value, ValueEqual>;
  using iterator = typename patricia_tree_t::iterator;
  using key_type = Key;
  using value_type = Value;

private:
  bool _is_bottom;
  patricia_tree_t _tree;

  class join_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator|(y);
      if (z.is_top()) {
        return {false, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };

    virtual bool default_is_absorbing() override {
      return true;
    }
  }; // class join_op

  class widening_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator||(y);
      if (z.is_top()) {
        return {false, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };

    virtual bool default_is_absorbing() override {
      return true;
    }

  }; // class widening_op

  template <typename Thresholds>
  class widening_thresholds_op : public binary_op_t {
    const Thresholds &m_ts;
  public:
    widening_thresholds_op(const Thresholds &ts) : m_ts(ts) {}

    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.widening_thresholds(y, m_ts);
      if (z.is_top()) {
        return {false, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };

    virtual bool default_is_absorbing() override {
      return true;
    }

  }; // class widening_thresholds_op

  class meet_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator&(y);
      if (z.is_bottom()) {
        return {true, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };

    virtual bool default_is_absorbing() override {
      return false;
    }

  }; // class meet_op

  class narrowing_op : public binary_op_t {
    virtual std::pair<bool, boost::optional<Value>>
    apply(const Key &/*key*/, const Value &x, const Value &y) override {
      Value z = x.operator&&(y);
      if (z.is_bottom()) {
        return {true, boost::optional<Value>()};
      } else {
        return {false, boost::optional<Value>(z)};
      }
    };

    virtual bool default_is_absorbing() override {
      return false;
    }

  }; // class narrowing_op

  class domain_po : public partial_order_t {
    virtual bool leq(const Value &x, const Value &y) override {
      return x.operator<=(y);
    }

    virtual bool default_is_top() override {
      return true;
    }

  }; // class domain_po

  static patricia_tree_t apply_operation(binary_op_t &o, patricia_tree_t t1,
                                         const patricia_tree_t &t2,
                                         bool &is_bottom) {
    is_bottom = t1.merge_with(t2, o);
    return t1;
  }

  separate_domain(patricia_tree_t &&t) : _is_bottom(false), _tree(std::move(t)) {}

  separate_domain(bool b) : _is_bottom(!b) {}

public:
  separate_domain() : _is_bottom(false) {}
  separate_domain(const separate_domain_t &o) = default;
  separate_domain(separate_domain_t &&o) = default;    
  separate_domain_t &operator=(const separate_domain_t &o) = default;
  separate_domain_t &operator=(separate_domain_t &&o) = default;

  static separate_domain_t top() { return separate_domain_t(); }
  static separate_domain_t bottom() { return separate_domain_t(false); }
  
  iterator begin() const {
    if (is_bottom()) {
      CRAB_ERROR("Separate domain: trying to invoke iterator on bottom");
    } else {
      return _tree.begin();
    }
  }

  iterator end() const {
    if (is_bottom()) {
      CRAB_ERROR("Separate domain: trying to invoke iterator on bottom");
    } else {
      return _tree.end();
    }
  }

  bool is_bottom() const { return _is_bottom; }

  bool is_top() const {
    return (!is_bottom() && _tree.size() == 0);
  }

  bool operator<=(const separate_domain_t &e) const {
    if (is_bottom()) {
      return true;
    } else if (e.is_bottom()) {
      return false;
    } else {
      domain_po po;
      return _tree.leq(e._tree, po);
    }
  }

  bool operator==(const separate_domain_t &e) const {
    return (this->operator<=(e) && e.operator<=(*this));
  }

  // Join
  separate_domain_t operator|(const separate_domain_t &e) const {
    if (is_bottom()) {
      return e;
    } else if (e.is_bottom()) {
      return *this;
    } else {
      join_op o;
      bool is_bottom /*unused*/;
      patricia_tree_t res = apply_operation(o, _tree, e._tree, is_bottom);
      return separate_domain_t(std::move(res));
    }
  }

  // Meet
  separate_domain_t operator&(const separate_domain_t &e) const {
    if (is_bottom() || e.is_bottom()) {
      return bottom();
    } else {
      meet_op o;
      bool is_bottom;
      patricia_tree_t res = apply_operation(o, _tree, e._tree, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return separate_domain_t(std::move(res));
      }
    }
  }

  // Widening
  separate_domain_t operator||(const separate_domain_t &e) const {
    if (is_bottom()) {
      return e;
    } else if (e.is_bottom()) {
      return *this;
    } else {
      widening_op o;
      bool is_bottom /*unused*/;
      patricia_tree_t res = apply_operation(o, _tree, e._tree, is_bottom);
      return separate_domain_t(std::move(res));
    }
  }

  // Widening with thresholds
  template <typename Thresholds>
  separate_domain_t widening_thresholds(const separate_domain_t &e,
                                        const Thresholds &ts) const {
    if (is_bottom()) {
      return e;
    } else if (e.is_bottom()) {
      return *this;
    } else {
      widening_thresholds_op<Thresholds> o(ts);
      bool is_bottom /*unused*/;
      patricia_tree_t res = apply_operation(o, _tree, e._tree, is_bottom);
      return separate_domain_t(std::move(res));
    }
  }

  // Narrowing
  separate_domain_t operator&&(const separate_domain_t &e) const {
    if (is_bottom() || e.is_bottom()) {
      return separate_domain_t(false);
    } else {
      narrowing_op o;
      bool is_bottom;
      patricia_tree_t res = apply_operation(o, _tree, e._tree, is_bottom);
      if (is_bottom) {
        return bottom();
      } else {
        return separate_domain_t(std::move(res));
      }
    }
  }

  void set(const Key &k, const Value &v) {
    if (!is_bottom()) {
      if (v.is_bottom()) {
        set_to_bottom();
      } else if (v.is_top()) {
        _tree.remove(k);
      } else {
        _tree.insert(k, v);
      }
    }
  }

  void set_to_bottom() {
    _is_bottom = true;
    _tree = patricia_tree_t();
  }

  separate_domain_t &operator-=(const Key &k) {
    if (!is_bottom()) {
      _tree.remove(k);
    }
    return *this;
  }

  // return null if k is not found
  const Value* find(const Key &k) const {
    assert(!is_bottom());
    return _tree.find(k);
  }

  Value at(const Key &k) const {
    if (is_bottom()) {
      return Value::bottom();
    } else {
      boost::optional<Value> v = _tree.lookup(k);
      if (v) {
        return *v;
      } else {
        return Value::top();
      }
    }
  }

  std::size_t size() const {
    if (is_bottom()) {
      return 0;
    } else if (is_top()) {
      CRAB_ERROR("separate_domains::size() is undefined if top");
    } else {
      return _tree.size();
    }
  }

  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
    } else {
      o << "{";
      for (typename patricia_tree_t::iterator it = _tree.begin();
           it != _tree.end();) {
        Key k = it->first;
        k.write(o);
        o << " -> ";
        Value v = it->second;
        v.write(o);
        ++it;
        if (it != _tree.end()) {
          o << "; ";
        }
      }
      o << "}";
    }
  }

  void project(const std::vector<Key> &keys) {
    if (is_bottom() || is_top()) {
      return;
    }
    const int factor = 60;   // between 0 and 100
    const int small_env = 5; // size of a small environment
    int num_total_keys = (int)_tree.size();
    int num_keys = (int)keys.size();
    if (num_total_keys <= small_env ||
        num_keys < num_total_keys * factor / 100) {
      // project on less than factor% of keys: we copy
      separate_domain_t env;
      for (auto key : keys) {
        env.set(key, at(key));
      }
      std::swap(*this, env);
    } else {
      // project on more or equal then factor% of keys: that
      // might be too many copies so we remove instead.
      std::vector<Key> sorted_keys(keys);
      std::sort(sorted_keys.begin(), sorted_keys.end());
      std::vector<Key> project_out_keys;
      project_out_keys.reserve(num_total_keys);
      for (auto it = _tree.begin(); it != _tree.end(); ++it) {
        Key key = it->first;
        if (!std::binary_search(sorted_keys.begin(), sorted_keys.end(), key)) {
          project_out_keys.push_back(key);
        }
      }
      for (auto const &key : project_out_keys) {
        this->operator-=(key);
      }
    }
  }

  // Assume that from does not have duplicates.
  void rename(const std::vector<Key> &from, const std::vector<Key> &to) {
    if (is_top() || is_bottom()) {
      // nothing to rename
      return;
    }
    if (from.size() != to.size()) {
      CRAB_ERROR(
          "separate_domain::rename with input vectors of different sizes");
    }

    if (::crab::CrabSanityCheckFlag) {
      std::set<Key> s1, s2;
      s1.insert(from.begin(), from.end());
      if (s1.size() != from.size()) {
        CRAB_ERROR("separate_domain::rename expects no duplicates");
      }
      s2.insert(to.begin(), to.end());
      if (s2.size() != to.size()) {
        CRAB_ERROR("separate_domain::rename expects no duplicates");
      }
    }

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      Key k = from[i];
      Key new_k = to[i];
      if (k == new_k) { // nothing to rename
        continue;
      }
      if (::crab::CrabSanityCheckFlag) {
        if (_tree.lookup(new_k)) {
          CRAB_ERROR("separate_domain::rename assumes that  ", new_k,
                     " does not exist in ", *this);
        }
      }
      if (boost::optional<Value> k_val_opt = _tree.lookup(k)) {
        if (!(*k_val_opt).is_top()) {
          _tree.insert(new_k, *k_val_opt);
        }
        _tree.remove(k);
      }
    }
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const separate_domain<Key, Value, ValueEqual> &d) {
    d.write(o);
    return o;
  }
}; // class separate_domain

} // namespace ikos




namespace crab {
namespace domains {
//===================================================================//  
//        The NOSA license does not apply to this code  
//===================================================================//  

/* 
 * Environment from Key to discrete_domain<Element>
 *
 * We cannot use separate_domain because if the value of a key-value
 * pair is bottom then the whole environment becomes bottom. Here, if
 * the value is bottom then it means that the value is simply the
 * "empty set" so we want to allow entries whose values are
 * bottom. But, similar to separate_domain, the whole environment can
 * still be bottom meaning failure/undefined/unreachable.
 */ 
template <typename Key, typename Element> class separate_discrete_domain {
private:
  using discrete_domain_t = ikos::discrete_domain<Element>;
  using patricia_tree_t = ikos::patricia_tree<Key, discrete_domain_t>;
  using unary_op_t = typename patricia_tree_t::unary_op_t;
  using binary_op_t = typename patricia_tree_t::binary_op_t;
  using partial_order_t = typename patricia_tree_t::partial_order_t;

public:
  using key_type = Key;
  using value_type = discrete_domain_t;    
  using separate_discrete_domain_t = separate_discrete_domain<Key, Element>;
  using iterator = typename patricia_tree_t::iterator;

private:
  bool m_is_bottom;
  patricia_tree_t m_tree;

  static patricia_tree_t apply_operation(binary_op_t &o,
                                         const patricia_tree_t &t1,
                                         const patricia_tree_t &t2) {
    patricia_tree_t res(t1);
    res.merge_with(t2, o);
    return res;
  }
  /* begin patricia_tree API */
  class join_op : public binary_op_t {
    std::pair<bool, boost::optional<value_type>>
    apply(const key_type &/*key*/, const value_type &x, const value_type &y) override {
      value_type z = x.operator|(y);
      if (z.is_top()) {
	// special encoding for top: the patricia tree will not keep
	// top values.
        return {false, boost::optional<value_type>()};
      } else {
        return {false, boost::optional<value_type>(z)};
      }
    }
    bool default_is_absorbing() override { return true; }    
  }; // class join_op

  class meet_op : public binary_op_t {
    std::pair<bool, boost::optional<value_type>>
    apply(const key_type &/*key*/, const value_type &x, const value_type &y) override {
      value_type z = x.operator&(y);
      // Returning this pair means that if z is bottom do not treat it
      // special and just update the patricia tree with z.
      return {false, boost::optional<value_type>(z)};
    };
    bool default_is_absorbing() override { return false; }
  }; // class meet_op

  class domain_po : public partial_order_t {
    bool leq(const value_type &x, const value_type &y) { return x.operator<=(y); }
    bool default_is_top() { return true; }
  }; // class domain_po
  /* end patricia_tree API */

  
  separate_discrete_domain(patricia_tree_t &&t)
    : m_is_bottom(false), m_tree(std::move(t)) {}

public:
  
  static separate_discrete_domain_t top() {
    return separate_discrete_domain_t();
  }

  static separate_discrete_domain_t bottom() {
    return separate_discrete_domain_t(true);
  }

  // Default constructor returns a top environment
  separate_discrete_domain(bool is_bottom = false)
    : m_is_bottom(is_bottom), m_tree(patricia_tree_t()) {}
  
  separate_discrete_domain(const separate_discrete_domain_t &o) = default;
  separate_discrete_domain(separate_discrete_domain_t &&o) = default;
  separate_discrete_domain_t &operator=(const separate_discrete_domain_t &o)  = default;
  separate_discrete_domain_t &operator=(separate_discrete_domain_t &&o)  = default;
  
  iterator begin() const {
    if (is_bottom()) {
      CRAB_ERROR("Separate discrete domain: trying to invoke iterator on bottom");
    } else {
      return m_tree.begin();
    }
  }

  iterator end() const {
    if (is_bottom()) {
      CRAB_ERROR("Separate discrete domain: trying to invoke iterator on bottom");
    } else {
      return m_tree.end();
    }
  }

  bool is_bottom() const { return m_is_bottom; }

  bool is_top() const {
    return (!is_bottom()&& m_tree.size() == 0);
  }

  bool operator<=(const separate_discrete_domain_t &other) const {
    if (is_bottom() || other.is_top()) {
      return true;
    } else if (other.is_bottom() || is_top()) {
      return false;
    } else {
      domain_po po;
      return m_tree.leq(other.m_tree, po);
    }
  }
  
  separate_discrete_domain_t
  operator|(const separate_discrete_domain_t &o) const {
    CRAB_LOG("separate-domain",
	     crab::outs() << "Join " << *this << " and " << o << "=";);      
    
    if (is_bottom() || o.is_top()) {
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << o << "\n";);
      return o;
    } else if (o.is_bottom() || is_top()) {
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << *this << "\n";);      
      return *this;
    } else {
      join_op op;
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree);
      auto out = separate_discrete_domain_t(std::move(res));
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << out << "\n";);
      return out;
    }
  }

  separate_discrete_domain_t
  operator&(const separate_discrete_domain_t &o) const {
    CRAB_LOG("separate-domain",
	     crab::outs() << "Meet " << *this << " and " << o << "=";);      
    
    if (is_top() || o.is_bottom()) {
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << o << "\n";);
      
      return o;
    } else if (o.is_top() || is_bottom()) {
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << *this << "\n";);
      
      return *this;
    } else {
      meet_op op;
      patricia_tree_t res = apply_operation(op, m_tree, o.m_tree);
      auto out = separate_discrete_domain_t(std::move(res));
      CRAB_LOG("separate-domain",
	       crab::outs() << "Res=" << out << "\n";);
      return out;
    }
  }

  void set(const Key &k, value_type v) {
    // Note that we can store a key-value pair where the value is
    // bottom because bottom means empty set.
    if (!is_bottom()) {
      if (v.is_top()) {
	m_tree.remove(k);
      } else {
	m_tree.insert(k, v);
      }
    }
  }

  separate_discrete_domain_t &operator-=(const Key &k) {
    if (!is_bottom()) {
      m_tree.remove(k);
    }
    return *this;
  }

  value_type at(const Key &k) const {
    if (is_bottom()) {
      CRAB_ERROR("separate_discrete_domain::at is undefined on bottom");
    } else if (is_top()) {
      return value_type::top();
    } else {
      if (boost::optional<value_type> v = m_tree.lookup(k)) {
        return *v;
      } else {
        return value_type::top();
      }
    }
  }

  void project(const std::vector<Key> &keys) {
    if (is_bottom() || is_top()) {
      return;
    }
    const int factor = 60;   // between 0 and 100
    const int small_env = 5; // size of a small environment
    int num_total_keys = (int)m_tree.size();
    int num_keys = (int)keys.size();
    if (num_total_keys <= small_env ||
        num_keys < num_total_keys * factor / 100) {
      // project on less than factor% of keys: we copy
      separate_discrete_domain_t env;
      for (auto key : keys) {
        env.set(key, at(key));
      }
      std::swap(*this, env);
    } else {
      // project on more or equal then factor% of keys: that
      // might be too many copies so we remove instead.
      std::vector<Key> sorted_keys(keys);
      std::sort(sorted_keys.begin(), sorted_keys.end());
      std::vector<Key> project_out_keys;
      project_out_keys.reserve(num_total_keys);
      for (auto it = m_tree.begin(); it != m_tree.end(); ++it) {
        Key key = it->first;
        if (!std::binary_search(sorted_keys.begin(), sorted_keys.end(), key)) {
          project_out_keys.push_back(key);
        }
      }
      for (auto const &key : project_out_keys) {
        this->operator-=(key);
      }
    }
  }

  // Assume that from does not have duplicates.
  void rename(const std::vector<Key> &from, const std::vector<Key> &to) {
    if (is_top() || is_bottom()) {
      // nothing to rename
      return;
    }
    if (from.size() != to.size()) {
      CRAB_ERROR(
          "separate_discrete_domain::rename with input vectors of different sizes");
    }

    if (::crab::CrabSanityCheckFlag) {
      std::set<Key> s1, s2;
      s1.insert(from.begin(), from.end());
      if (s1.size() != from.size()) {
        CRAB_ERROR("separate_discrete_domain::rename expects no duplicates");
      }
      s2.insert(to.begin(), to.end());
      if (s2.size() != to.size()) {
        CRAB_ERROR("separate_discrete_domain::rename expects no duplicates");
      }
    }

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      Key k = from[i];
      Key new_k = to[i];
      if (k == new_k) { // nothing to rename
        continue;
      }
      if (::crab::CrabSanityCheckFlag) {
        if (m_tree.lookup(new_k)) {
          CRAB_ERROR("separate_discrete_domain::rename assumes that  ",
		     new_k, " does not exist in ", *this);
        }
      }
      if (boost::optional<value_type> val_opt = m_tree.lookup(k)) {
        if (!(*val_opt).is_top()) {
          m_tree.insert(new_k, *val_opt);
        }
        m_tree.remove(k);
      }
    }
  }

  
  void write(crab::crab_os &o) const {
    if (is_top()) {
      o << "{}";
    } else if (is_bottom()) {
      o << "_|_";
    } else {
      o << "{";
      for (typename patricia_tree_t::iterator it = m_tree.begin();
           it != m_tree.end();) {
        Key k = it->first;
        k.write(o);
        o << " -> ";
        value_type v = it->second;
        v.write(o);
        ++it;
        if (it != m_tree.end()) {
          o << "; ";
        }
      }
      o << "}";
    }
  }
}; // class separate_discrete_domain

template <typename Key, typename Element>
inline crab_os &operator<<(crab_os &o,
                           const separate_discrete_domain<Key, Element> &env) {
  env.write(o);
  return o;
}

} //end namespace domains  
} //end namespace crab
