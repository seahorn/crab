/*******************************************************************************
 * Extend abstract domains with very specialized operations.
 *
 * FIXME: if we use the generic abtract domain (based on the
 * type-erasure idiom) we cannot extend abstract domains using traits
 * because we cannot access to the underlying domains.
 ******************************************************************************/

#pragma once
#include <crab/domains/generic_abstract_domain.hpp>

namespace crab {
namespace domains {

// Perform constraint simplifications depending on the abstract domain
template <typename Domain> class constraint_simp_domain_traits {
public:
  using number_t = typename Domain::number_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  typedef
      typename Domain::linear_constraint_system_t linear_constraint_system_t;

  // Convert an equality into two inequalities. This is not
  // possible for machine arithmetic domains.
  static void lower_equality(linear_constraint_t cst,
                             linear_constraint_system_t &csts) {
    if (cst.is_equality()) {
      csts += linear_constraint_t(cst.expression(),
                                  linear_constraint_t::INEQUALITY);
      csts += linear_constraint_t(cst.expression() * number_t(-1),
                                  linear_constraint_t::INEQUALITY);
    } else {
      csts += cst;
    }
  }
};

// Special operations needed by the checker
template <typename Domain> class checker_domain_traits {
public:
  using varname_t = typename Domain::varname_t;
  using number_t = typename Domain::number_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  typedef
      typename Domain::linear_constraint_system_t linear_constraint_system_t;
  using disjunctive_linear_constraint_system_t =
      typename Domain::disjunctive_linear_constraint_system_t;

private:
  struct entailment {
    Domain _dom;
    entailment(Domain dom) : _dom(dom) {}
    bool operator()(const linear_constraint_t &cst) {
      Domain dom(_dom); // copy is necessary
      linear_constraint_t neg_cst = cst.negate();
      dom += neg_cst;
      return dom.is_bottom();
    }
  };

  // Return true if (c1 or c2 or ... cn) entails (d1 and d2 and .. dn)
  static bool __entail(const disjunctive_linear_constraint_system_t &lhs,
                       const linear_constraint_system_t &rhs) {
    // -- trivial cases first
    if (rhs.is_false()) {
      return false;
    } else if (rhs.is_true()) {
      return true;
    } else if (lhs.is_false()) {
      return true;
    } else if (lhs.is_true()) {
      return false;
    }

    // -- return true if for all ci :: ci entails (d1 and d2 and .. dn)
    return std::all_of(
        lhs.begin(), lhs.end(), [&rhs](const linear_constraint_system_t &csts) {
          Domain lhs;
          lhs += csts;
          return std::all_of(
              rhs.begin(), rhs.end(),
              [&lhs](const linear_constraint_t &c) { return entail(lhs, c); });
        });
  }

public:
  /*
     Public API

     static bool entail(Domain&, const disjunctive_linear_constraint_system_t&);
     static bool entail(const disjunctive_linear_constraint_system_t&, Domain&);
     static bool entail(Domain&, const linear_constraint_t&);

     static bool intersect(Domain&, const linear_constraint_t&);
   */

  // Return true if lhs entails (c1 or c2 or ... cn)
  static bool entail(Domain &lhs,
                     const disjunctive_linear_constraint_system_t &rhs) {
    // -- trivial cases first
    if (rhs.is_false()) {
      return false;
    } else if (rhs.is_true()) {
      return true;
    } else if (lhs.is_bottom()) {
      return true;
    } else if (lhs.is_top()) {
      return false;
    }
    // -- return true if exists ci such that lhs entails ci
    for (linear_constraint_system_t csts : rhs) {
      if (std::all_of(csts.begin(), csts.end(),
                      [&lhs](const linear_constraint_t &c) {
                        return entail(lhs, c);
                      })) {
        return true;
      }
    }
    return false;
  }

  // Return true if (c1 or c2 or ... cn) entails rhs
  static bool entail(const disjunctive_linear_constraint_system_t &lhs,
                     Domain &rhs) {
    auto csts = rhs.to_linear_constraint_system();
    return __entail(lhs, csts);
  }

  // Return true if lhs entails rhs.
  static bool entail(Domain &lhs, const linear_constraint_t &rhs) {
    if (lhs.is_bottom())
      return true;
    if (rhs.is_tautology())
      return true;
    if (rhs.is_contradiction())
      return false;

    CRAB_LOG("checker-entailment", linear_constraint_t tmp(rhs);
             crab::outs() << "Checking whether\n"
                          << lhs << "\nentails " << tmp << "\n";);

    bool res;
    entailment op(lhs);
    if (rhs.is_equality()) {
      // try to convert the equality into inequalities so when it's
      // negated we do not have disequalities.
      linear_constraint_system_t inequalities;
      constraint_simp_domain_traits<Domain>::lower_equality(rhs, inequalities);
      res = std::all_of(inequalities.begin(), inequalities.end(), op);
    } else {
      res = op(rhs);
    }

    CRAB_LOG("checker-entailment",
             if (res) { crab::outs() << "\t**entailment holds.\n"; } else {
               crab::outs() << "\t**entailment does not hold.\n";
             });

    // Note: we cannot convert rhs into Domain and then use the <=
    //       operator. The problem is that we cannot know for sure
    //       whether Domain can represent precisely rhs. It is not
    //       enough to do something like
    //
    //       Dom dom = rhs;
    //       if (dom.is_top()) { ... }

    return res;
  }

  // Return true if inv intersects with cst.
  static bool intersect(Domain &inv, const linear_constraint_t &cst) {
    if (inv.is_bottom() || cst.is_contradiction())
      return false;
    if (inv.is_top() || cst.is_tautology())
      return true;

    Domain dom(inv);
    dom += cst;
    return !dom.is_bottom();
  }
};

// Special operations for applying reduction between domains.
template <typename Domain> class reduced_domain_traits {
public:
  using variable_t = typename Domain::variable_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  typedef
      typename Domain::linear_constraint_system_t linear_constraint_system_t;

  static_assert(
      !std::is_same<Domain, abstract_domain<typename Domain::variable_t>>::value,
      "reduced_domain_traits not supported for generic domain");
  static_assert(
      !std::is_same<Domain,
                   abstract_domain_ref<typename Domain::variable_t>>::value,
      "reduced_domain_traits not supported for generic domain");

  // extract linear constraints from dom involving x and store in ctsts
  static void extract(Domain &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    auto all_csts = dom.to_linear_constraint_system();
    for (auto const &cst : all_csts) {
      if (only_equalities && (!cst.is_equality())) {
        continue;
      }
      if (std::find(std::begin(cst.variables()), std::end(cst.variables()), x) !=
          std::end(cst.variables())) {
        csts += cst;
      }
    }
  }
};

// Some abstract domains have global state that should be reset by
// the client if multiple crab instances will be run. This is just
// a temporary hack. The proper solution is to avoid global state.
template <typename Domain> class special_domain_traits {
public:
  static void clear_global_state(void) {}
};

} // end namespace domains
} // end namespace crab
