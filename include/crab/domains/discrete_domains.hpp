/*******************************************************************************
 *
 * An implementation of discrete domains based on Patricia trees.
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
#include <crab/support/debug.hpp>

#include <boost/range.hpp>

namespace ikos {

template <typename Element> class discrete_domain {

private:
  using ptset_t = patricia_tree_set<Element>;

public:
  using discrete_domain_t = discrete_domain<Element>;
  using iterator = typename ptset_t::iterator;

private:
  bool m_is_top;
  ptset_t m_set;

private:
  discrete_domain(bool is_top) : m_is_top(is_top) {}
  discrete_domain(ptset_t set) : m_is_top(false), m_set(set) {}

public:
  // Return an empty set
  static discrete_domain_t bottom() { return discrete_domain_t(false); }
  // Return a set with all variables
  static discrete_domain_t top() { return discrete_domain_t(true); }

public:
  // Default constructor creates an empty set rather than top
  discrete_domain() : m_is_top(false) {}

  discrete_domain(const discrete_domain_t &other) = default;
  discrete_domain(discrete_domain_t &&other) = default;
  discrete_domain_t &operator=(const discrete_domain_t &other) = default;
  discrete_domain_t &operator=(discrete_domain_t &&other) = default;  

  discrete_domain(Element s) : m_is_top(false), m_set(s) {}

  template <typename Iterator>
  discrete_domain(Iterator eIt, Iterator eEt) : m_is_top(false) {
    for (auto e : boost::make_iterator_range(eIt, eEt)) {
      m_set += e;
    }
  }

  bool is_top() const { return m_is_top; }

  bool is_bottom() const { return (!m_is_top && m_set.empty()); }

  bool operator<=(const discrete_domain_t &other) const {
    return other.m_is_top || (!m_is_top && m_set <= other.m_set);
  }

  bool operator==(const discrete_domain_t &other) const {
    return (m_is_top && other.m_is_top) || (m_set == other.m_set);
  }

  void operator|=(const discrete_domain_t &other) { *this = *this | other; }

  discrete_domain_t operator|(const discrete_domain_t &other) const {
    if (m_is_top || other.m_is_top) {
      return discrete_domain_t(true);
    } else {
      return discrete_domain_t(m_set | other.m_set);
    }
  }

  discrete_domain_t operator&(const discrete_domain_t &other) const {
    if (is_bottom() || other.is_bottom()) {
      return discrete_domain_t(false);
    } else if (m_is_top) {
      return other;
    } else if (other.m_is_top) {
      return *this;
    } else {
      return discrete_domain_t(m_set & other.m_set);
    }
  }

  discrete_domain_t operator||(const discrete_domain_t &other) const {
    return operator|(other);
  }

  discrete_domain_t operator&&(const discrete_domain_t &other) const {
    return operator&(other);
  }

  discrete_domain_t &operator+=(Element s) {
    if (!m_is_top) {
      m_set += s;
    }
    return *this;
  }

  template <typename Range> discrete_domain_t &operator+=(Range es) {
    if (!m_is_top) {
      for (auto e : es) {
        m_set += e;
      }
    }
    return *this;
  }

  discrete_domain_t operator+(Element s) {
    discrete_domain_t r(*this);
    r.operator+=(s);
    return r;
  }

  template <typename Range> discrete_domain_t operator+(Range es) {
    discrete_domain_t r(*this);
    r.operator+=(es);
    return r;
  }

  discrete_domain_t &operator-=(Element s) {
    if (!m_is_top) {
      m_set -= s;
    }
    return *this;
  }

  template <typename Range> discrete_domain_t &operator-=(Range es) {
    if (!m_is_top) {
      for (auto e : es) {
        m_set -= e;
      }
    }
    return *this;
  }

  discrete_domain_t operator-(Element s) {
    discrete_domain_t r(*this);
    r.operator-=(s);
    return r;
  }

  template <typename Range> discrete_domain_t operator-(Range es) {
    discrete_domain_t r(*this);
    r.operator-=(es);
    return r;
  }

  std::size_t size() const {
    if (m_is_top) {
      assert(false);
      CRAB_ERROR("Size for discrete domain TOP is undefined");
    } else {
      return m_set.size();
    }
  }

  iterator begin() const {
    if (m_is_top) {
      assert(false);      
      CRAB_ERROR("Iterator for discrete domain TOP is undefined");
    } else {
      return m_set.begin();
    }
  }

  iterator end() const {
    if (m_is_top) {
      assert(false);      
      CRAB_ERROR("Iterator for discrete domain TOP is undefined");
    } else {
      return m_set.end();
    }
  }

  void write(crab::crab_os &o) const {
    if (m_is_top) {
      o << "{...}";
    } else if (m_set.empty()) {
      o << "{}";
    } else {
      o << m_set;
    }
  }

}; // class discrete_domain

template <typename Elem>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const discrete_domain<Elem> &d) {
  d.write(o);
  return o;
}

} // namespace ikos
