/*******************************************************************************
 *
 * Standard domain of numerical congruences extended with bitwise
 * operations.
 *
 * Author: Alexandre C. D. Wimmers (alexandre.c.wimmers@nasa.gov)
 *
 * Contributors: Jorge A. Navas (jorge.a.navaslaserna@nasa.gov)
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

#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>

#include <boost/optional.hpp>

namespace ikos {

template <typename Number> class congruence {
  typedef interval<Number> interval_t;

public:
  typedef congruence<Number> congruence_t;

private:
  bool _is_bottom;

  /// A congruence is denoted by aZ + b, where b \in Z and a \in N.
  /// The abstract state aZ + b represents all numbers that are
  /// congruent to b modulo a.
  Number _a; // modulo
  Number _b; // remainder

  // Notes about the % operator
  //
  // The semantics of r = n % d is to set r to "n mod d". The sign of
  // the d is ignored and r is always non-negative.
  //
  // We assume that n % d (also n /d) raises a runtime error if d==0.

  void normalize(void) {
    // Set to standard form: 0 <= b < a for a != 0
    if (_a != 0) {
      _b = _b % _a;
    }
  }

  // if true then top (1Z + 0) else bottom
  congruence(bool b) : _is_bottom(!b), _a(1), _b(0) {}

  congruence(int n) : _is_bottom(false), _a(0), _b(n) {}

  congruence(Number a, Number b) : _is_bottom(false), _a(a), _b(b) {
    normalize();
  }

  Number abs(Number x) const { return x < 0 ? -x : x; }

  Number max(Number x, Number y) const { return x.operator<=(y) ? y : x; }

  Number min(Number x, Number y) const { return x.operator<(y) ? x : y; }

  Number gcd(Number x, Number y, Number z) const { return gcd(x, gcd(y, z)); }
  // Not to be called explicitly outside of gcd
  Number gcd_helper(Number x, Number y) const {
    return (y == 0) ? x : gcd_helper(y, x % y);
  }
  Number gcd(Number x, Number y) const { return gcd_helper(abs(x), abs(y)); }

  Number lcm(Number x, Number y) const {
    Number tmp = gcd(x, y);
    return abs(x * y) / tmp;
  }

  bool is_zero() const { return !is_bottom() && _a == 0 && _b == 0; }

  bool all_ones() const { return !is_bottom() && _a == 0 && _b == -1; }

  interval_t to_interval() const {
    assert(singleton());
    return interval_t(*(singleton()));
  }

public:
  static congruence_t top() { return congruence(true); }

  static congruence_t bottom() { return congruence(false); }

  congruence() : _is_bottom(false), _a(1), _b(0) {}

  congruence(Number n) : _is_bottom(false), _a(0), _b(n) {}

  congruence(const congruence_t &o)
      : _is_bottom(o._is_bottom), _a(o._a), _b(o._b) {}

  congruence_t operator=(congruence_t o) {
    _is_bottom = o._is_bottom;
    _a = o._a;
    _b = o._b;
    return *this;
  }

  bool is_bottom() const { return _is_bottom; }

  bool is_top() const { return _a == 1; }

  boost::optional<Number> singleton() const {
    if (!this->is_bottom() && _a == 0) {
      return boost::optional<Number>(_b);
    } else {
      return boost::optional<Number>();
    }
  }

  Number get_modulo() const { return _a; }

  Number get_remainder() const { return _b; }

  bool operator==(congruence_t o) const {
    return (is_bottom() == o.is_bottom() && _a == o._a && _b == o._b);
  }

  bool operator!=(congruence_t x) const { return !this->operator==(x); }

  /** Lattice Operations **/

  bool operator<=(congruence_t o) {
    if (is_bottom()) {
      return true;
    } else if (o.is_bottom()) {
      return false;
    } else if (_a == 0 && o._a == 0) {
      return (_b == o._b);
    } else if (_a == 0) {
      if ((_b % o._a) == (o._b % o._a)) {
        return true;
      }
    } else if (o._a == 0) {
      if (_b % _a == (o._b % _a)) {
        return false;
      }
    }
    return (_a % o._a == 0) && (_b % o._a == o._b % o._a);
  }

  congruence_t operator|(congruence_t o) {
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() || o.is_top()) {
      return top();
    } else {
      return congruence_t(gcd(_a, o._a, abs(_b - o._b)), min(_b, o._b));
    }
  }

  congruence_t operator&(congruence_t o) {
    if (is_bottom() || o.is_bottom()) {
      return bottom();
    }

    // lcm has meaning only if both a and o.a are not 0
    if (_a == 0 && o._a == 0) {
      if (_b == o._b) {
        return *this;
      } else {
        return bottom();
      }
    } else if (_a == 0) {
      // b & a'Z + b' iff \exists k such that a'*k + b' = b iff ((b - b') %a' ==
      // 0)
      if ((_b - o._b) % o._a == 0) {
        return *this;
      } else {
        return bottom();
      }
    } else if (o._a == 0) {
      // aZ+b & b' iff \exists k such that a*k+b  = b' iff ((b'-b %a) == 0)
      if ((o._b - _b) % _a == 0) {
        return o;
      } else {
        return bottom();
      }
    } else {
      // pre: a and o.a != 0
      Number x = gcd(_a, o._a);
      if (_b % x == (o._b % x)) {
        // the part max(b,o.b) needs to be verified. What we really
        // want is to find b'' such that
        // 1) b'' % lcm(a,a') == b  % lcm(a,a'), and
        // 2) b'' % lcm(a,a') == b' % lcm(a,a').
        // An algorithm for that is provided in Granger'89.
        return congruence_t(lcm(_a, o._a), max(_b, o._b));
      } else {
        return congruence_t::bottom();
      }
    }
  }

  congruence_t operator||(congruence_t o) {
    // Equivalent to join, domain is flat
    return *this | o;
  }

  congruence_t operator&&(congruence_t o) {
    // Simply refines top element
    return (is_top()) ? o : *this;
  }

  /** Arithmetic Operators **/

  congruence_t operator+(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else
      return congruence_t(gcd(_a, o._a), _b + o._b);
  }

  congruence_t operator-(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else
      return congruence_t(gcd(_a, o._a), _b - o._b);
  }

  congruence_t operator-() {
    if (this->is_bottom() || this->is_top())
      return *this;
    else
      return congruence_t(_a, -_b + _a);
  }

  congruence_t operator*(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if ((this->is_top() || o.is_top()) && _a != 0 && o._a != 0)
      return congruence_t::top();
    else
      return congruence_t(gcd(_a * o._a, _a * o._b, o._a * _b), _b * o._b);
  }

  // signed division
  congruence_t operator/(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (o == congruence(0))
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {

      /*
         aZ+b / 0Z+b':
            if b'|a then  (a/b')Z + b/b'
            else          top
      */
      if (o._a == 0) {
        if (_a % o._b == 0)
          return congruence_t(_a / o._b, _b / o._b);
        else
          return congruence_t::top();
      }

      /*
         0Z+b / a'Z+b':
            if N>0   (b div N)Z + 0
            else     0Z + 0

           where N = a'((b-b') div a') + b'
      */
      if (_a == 0) {
        Number n(o._a * (((_b - o._b) / o._a) + o._b));
        if (n > 0) {
          return congruence_t(_b / n, 0);
        } else {
          return congruence_t(0, 0);
        }
      }

      /*
        General case: no singleton
      */
      return congruence_t::top();
    }
  }

  // signed remainder operator
  congruence_t operator%(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (o == congruence(0))
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {
      /*
         aZ+b mod 0Z+b':
             if b'|a then  (a/b')Z + b/b'
             else          top
      */
      if (o._a == 0) {
        if (_a % o._b == 0) {
          return congruence_t(0, _b % o._b);
        } else {
          return congruence_t(gcd(_a, o._b), _b);
        }
      }

      /*
          0Z+b mod a'Z+b':
           if N<=0           then 0Z+b
           if (b div N) == 1 then gcd(b',a')Z + b
           if (b div N) >= 2 then N(b div N)Z  + b

         where N = a'((b-b') div a') + b'
      */
      if (_a == 0) {
        Number n(o._a * (((_b - o._b) / o._a) + o._b));
        if (n <= 0) {
          return congruence_t(_a, _b);
        } else if (_b == n) {
          return congruence_t(gcd(o._b, o._a), _b);
        } else if ((_b / n) >= 2) {
          return congruence_t(_b, _b);
        } else {
          CRAB_ERROR("unreachable");
        }
      }

      /*
          general case: no singleton
      */
      return congruence_t(gcd(_a, o._a, o._b), _b);
    }
  }

  /**
      Bitwise operators.
      They are very imprecise because we ignore bitwidth.

      Bitwise operation can be implemented more precisely based on
      Stefan Bygde's paper: Static WCET analysis based on abstract
      interpretation and counting of elements, Vasteras : School of
      Innovation, Design and Engineering, Malardalen University (2010).
   **/

  congruence_t And(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {
      if (is_zero() || o.is_zero()) {
        return congruence_t(0);
      } else if (all_ones()) {
        return o;
      } else if (o.all_ones()) {
        return *this;
      } else if (_a == 0 && o._a == 0) {
        return congruence_t(_b & o._b);
      } else {
        return top();
      }
    }
  }

  congruence_t Or(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {
      if (all_ones() || o.all_ones()) {
        return congruence_t(-1);
      } else if (is_zero()) {
        return o;
      } else if (o.is_zero()) {
        return *this;
      } else if (_a == 0 && o._a == 0) {
        return congruence_t(_b | o._b);
      } else {
        return top();
      }
    }
  }

  congruence_t Xor(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {
      if (is_zero()) {
        return o;
      } else if (o.is_zero()) {
        return *this;
      } else if (_a == 0 && o._a == 0) {
        return congruence_t(_b ^ o._b);
      } else {
        return top();
      }
    }
  }

  congruence_t Shl(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {

      if (o._a == 0) { // singleton

        if (o._b < 0) {
          return bottom();
        }

        // aZ + b << 0Z + b'  = (a*2^b')Z + b*2^b'
        Number x = Number(1) << o._b;
        return congruence_t(_a * x, _b * x);
      } else {

        Number x = Number(1) << o._b;
        Number y = Number(1) << o._a;
        // aZ + b << a'Z + b' = (gcd(a, b * (2^a' - 1)))*(2^b')Z + b*(2^b')
        return congruence_t(gcd(_a, _b * (y - 1)) * x, _b * x);
      }
    }
  }

  congruence_t AShr(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {

      if (o._a == 0) { // singleton
        // aZ + b >> 0Z + b'
        if (o._b < 0) {
          return congruence_t::bottom();
        }
      }

      if (singleton() && o.singleton()) {
        interval_t res = to_interval().AShr(o.to_interval());
        if (boost::optional<Number> n = res.singleton()) {
          return congruence(*n);
        }
      }

      return congruence_t::top();
    }
  }

  congruence_t LShr(congruence_t o) {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else {

      if (o._a == 0) {
        // aZ + b >> 0Z + b'
        if (o._b < 0) {
          return congruence_t::bottom();
        }
      }

      if (singleton() && o.singleton()) {
        interval_t res = to_interval().LShr(o.to_interval());
        if (boost::optional<Number> n = res.singleton()) {
          return congruence(*n);
        }
      }

      return congruence_t::top();
    }
  }

  // division and remainder operations

  congruence_t SDiv(congruence_t x) { return this->operator/(x); }

  congruence_t UDiv(congruence_t x) { return congruence_t::top(); }

  congruence_t SRem(congruence_t x) { return this->operator%(x); }

  congruence_t URem(congruence_t x) { return congruence_t::top(); }

  void write(crab::crab_os &o) const {
    if (is_bottom()) {
      o << "_|_";
      return;
    }

    if (_a == 0) {
      o << _b;
      return;
    }

    o << _a << "Z+" << _b;
  }
}; // end class congruence

template <typename Number>
inline crab::crab_os &operator<<(crab::crab_os &o,
                                 const congruence<Number> &c) {
  c.write(o);
  return o;
}

template <typename Number>
inline congruence<Number> operator+(Number c, congruence<Number> x) {
  return congruence<Number>(c) + x;
}

template <typename Number>
inline congruence<Number> operator+(congruence<Number> x, Number c) {
  return x + congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator*(Number c, congruence<Number> x) {
  return congruence<Number>(c) * x;
}

template <typename Number>
inline congruence<Number> operator*(congruence<Number> x, Number c) {
  return x * congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator/(Number c, congruence<Number> x) {
  return congruence<Number>(c) / x;
}

template <typename Number>
inline congruence<Number> operator/(congruence<Number> x, Number c) {
  return x / congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator-(Number c, congruence<Number> x) {
  return congruence<Number>(c) - x;
}

template <typename Number>
inline congruence<Number> operator-(congruence<Number> x, Number c) {
  return x - congruence<Number>(c);
}

template <typename Number, typename VariableName, typename CongruenceCollection>
class equality_congruence_solver {
  // TODO: check correctness of the solver. Granger provides a sound
  // and more precise solver for equality linear congruences (see
  // Theorem 4.4).
private:
  typedef congruence<Number> congruence_t;
  typedef variable<Number, VariableName> variable_t;
  typedef linear_expression<Number, VariableName> linear_expression_t;
  typedef linear_constraint<Number, VariableName> linear_constraint_t;
  typedef linear_constraint_system<Number, VariableName>
      linear_constraint_system_t;

  typedef std::vector<linear_constraint_t> cst_table_t;
  typedef typename linear_constraint_t::variable_set_t variable_set_t;

  std::size_t m_max_cycles;
  bool m_is_contradiction;
  cst_table_t m_cst_table;
  variable_set_t m_refined_variables;
  std::size_t m_op_count;

private:
  bool refine(variable_t v, congruence_t i, CongruenceCollection &env) {
    congruence_t old_i = env[v];
    congruence_t new_i = old_i & i;
    if (new_i.is_bottom()) {
      return true;
    }
    if (old_i != new_i) {
      env.set(v, new_i);
      m_refined_variables += v;
      ++(m_op_count);
    }
    return false;
  }

  congruence_t compute_residual(linear_constraint_t cst, variable_t pivot,
                                CongruenceCollection &env) {
    congruence_t residual(cst.constant());
    for (typename linear_constraint_t::iterator it = cst.begin();
         it != cst.end(); ++it) {
      variable_t v = it->second;
      if (!(v == pivot)) {
        residual = residual - (it->first * env[v]);
        ++(m_op_count);
      }
    }
    return residual;
  }

  bool propagate(linear_constraint_t cst, CongruenceCollection &env) {
    for (typename linear_constraint_t::iterator it = cst.begin();
         it != cst.end(); ++it) {
      Number c = it->first;
      variable_t pivot = it->second;
      congruence_t rhs = compute_residual(cst, pivot, env) / congruence_t(c);

      if (cst.is_equality()) {
        if (refine(pivot, rhs, env)) {
          return true;
        }
      } else if (cst.is_inequality() || cst.is_strict_inequality()) {
        // Inequations (>=, <=, >, and <) do not work well with
        // congruences because for any number n there is always x and y
        // \in gamma(aZ+b) such that n < x and n > y.
        //
        // The only cases we can catch is when all the expressions
        // are constants. We do not bother because any product
        // with intervals or constants should get those cases.
        continue;
      } else {
        // TODO: cst is a disequation
      }
    }
    return false;
  }

  bool solve_system(CongruenceCollection &env) {
    std::size_t cycle = 0;
    do {
      ++cycle;
      m_refined_variables.clear();
      for (typename cst_table_t::iterator it = m_cst_table.begin();
           it != m_cst_table.end(); ++it) {
        if (propagate(*it, env)) {
          return true;
        }
      }
    } while (m_refined_variables.size() > 0 && cycle <= m_max_cycles);
    return false;
  }

public:
  equality_congruence_solver(linear_constraint_system_t csts,
                             std::size_t max_cycles)
      : m_max_cycles(max_cycles), m_is_contradiction(false) {
    for (typename linear_constraint_system_t::iterator it = csts.begin();
         it != csts.end(); ++it) {
      linear_constraint_t cst = *it;
      if (cst.is_contradiction()) {
        m_is_contradiction = true;
        return;
      } else if (cst.is_tautology()) {
        continue;
      } else {
        m_cst_table.push_back(cst);
      }
    }
  }

  void run(CongruenceCollection &env) {
    if (m_is_contradiction) {
      env.set_to_bottom();
    } else {
      if (solve_system(env)) {
        env.set_to_bottom();
      }
    }
  }

}; // class equality_congruence_solver

template <typename Number, typename VariableName>
class congruence_domain final : public crab::domains::abstract_domain<
                                    congruence_domain<Number, VariableName>> {
public:
  typedef congruence<Number> congruence_t;

private:
  // note that this is assuming that all variables have the same bit
  // width which is unrealistic.
  typedef congruence_domain<Number, VariableName> congruence_domain_t;
  typedef crab::domains::abstract_domain<congruence_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  using typename abstract_domain_t::pointer_constraint_t;

private:
  typedef separate_domain<variable_t, congruence_t> separate_domain_t;
  typedef equality_congruence_solver<number_t, varname_t, separate_domain_t>
      solver_t;

public:
  typedef typename separate_domain_t::iterator iterator;

private:
  separate_domain_t _env;

private:
  congruence_domain(separate_domain_t env) : _env(env) {}

public:
  void set_to_top() {
    congruence_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    congruence_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  congruence_domain() : _env(separate_domain_t::top()) {}

  congruence_domain(const congruence_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  congruence_domain_t &operator=(const congruence_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }

  iterator begin() { return this->_env.begin(); }

  iterator end() { return this->_env.end(); }

  bool is_bottom() { return this->_env.is_bottom(); }

  bool is_top() { return this->_env.is_top(); }

  bool operator<=(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");
    return this->_env <= e._env;
  }

  void operator|=(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    this->_env = this->_env | e._env;
  }

  congruence_domain_t operator|(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");
    return this->_env | e._env;
  }

  congruence_domain_t operator&(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");
    return this->_env & e._env;
  }

  congruence_domain_t operator||(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");
    return this->_env || e._env;
  }

  congruence_domain_t
  widening_thresholds(congruence_domain_t other,
                      const crab::iterators::thresholds<number_t> &) {
    return (*this || other);
  }

  congruence_domain_t operator&&(congruence_domain_t e) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
    return this->_env && e._env;
  }

  void set(variable_t v, congruence_t i) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, i);
  }

  void set(variable_t v, number_t n) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");
    this->_env.set(v, congruence_t(n));
  }

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");
    this->_env -= v;
  }

  void operator-=(std::vector<variable_t> vs) {
    for (typename std::vector<variable_t>::iterator it = vs.begin(),
                                                    end = vs.end();
         it != end; ++it) {
      this->operator-= *it;
    }
  }

  congruence_t operator[](variable_t v) { return this->_env[v]; }

  congruence_t operator[](linear_expression_t expr) {
    congruence_t r(expr.constant());
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      congruence_t c(it->first);
      r = r + (c * this->_env[it->second]);
    }
    return r;
  }

  void add(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");
    const std::size_t threshold = 10;
    if (!this->is_bottom()) {
      solver_t solver(csts, threshold);
      solver.run(this->_env);
    }
  }

  void operator+=(linear_constraint_system_t csts) { this->add(csts); }

  congruence_domain_t operator+(linear_constraint_system_t csts) {
    congruence_domain_t e(this->_env);
    e += csts;
    return e;
  }

  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    congruence_t r = e.constant();
    for (typename linear_expression_t::iterator it = e.begin(); it != e.end();
         ++it) {
      r = r + (it->first * this->_env[it->second]);
    }
    this->_env.set(x, r);
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    congruence_t yi = this->_env[y];
    congruence_t zi = this->_env[z];
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    congruence_t yi = this->_env[y];
    congruence_t zi(k);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case OP_SDIV:
      xi = yi / zi;
      break;
    case OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case OP_SREM:
      xi = yi.SRem(zi);
      break;
    case OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  // backward operations
  void backward_assign(variable_t x, linear_expression_t e,
                       congruence_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign");

    crab::domains::BackwardAssignOps<congruence_domain_t>::assign(*this, x, e,
                                                                  inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      congruence_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<congruence_domain_t>::apply(*this, op, x,
                                                                 y, z, inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      congruence_domain_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<congruence_domain_t>::apply(*this, op, x,
                                                                 y, z, inv);
  }

  // cast operations

  void apply(crab::domains::int_conv_operation_t /*op*/, variable_t dst,
             variable_t src) {
    // ignore widths
    assign(dst, src);
  }

  // bitwise operations

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    congruence_t yi = this->_env[y];
    congruence_t zi = this->_env[z];
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default: {
      CRAB_ERROR("unreachable");
    }
    }
    this->_env.set(x, xi);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    congruence_t yi = this->_env[y];
    congruence_t zi(k);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default: {
      CRAB_ERROR("unreachable");
    }
    }
    this->_env.set(x, xi);
  }

  /*
     Begin unimplemented operations

     congruence_domain implements only standard abstract operations of
     a numerical domain.  The implementation of boolean, array, or
     pointer operations is empty because they should never be called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(crab::domains::bool_operation_t op, variable_t x,
                         variable_t y, variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                congruence_domain_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                congruence_domain_t invariant) {}
  void backward_apply_binary_bool(crab::domains::bool_operation_t op,
                                  variable_t x, variable_t y, variable_t z,
                                  congruence_domain_t invariant) {}
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           congruence_domain_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           congruence_domain_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            congruence_domain_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            congruence_domain_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  congruence_domain_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  congruence_domain_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             congruence_domain_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    separate_domain_t env;
    for (variable_t var : variables) {
      env.set(var, this->_env[var]);
    }
    std::swap(_env, env);
  }

  void expand(variable_t x, variable_t new_x) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set(new_x, this->_env[x]);
  }

  void normalize() {}

  void minimize() {}

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    _env.rename(from, to);    
  }
  
  /* begin intrinsics operations */  
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  congruence_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void write(crab::crab_os &o) {
    crab::CrabStats::count(getDomainName() + ".count.write");
    crab::ScopedCrabStats __st__(getDomainName() + ".write");

    this->_env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      variable_t v = it->first;
      congruence_t c = it->second;
      boost::optional<number_t> n = c.singleton();
      if (n) {
        csts += (v == *n);
      }
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  static std::string getDomainName() { return "Congruences"; }

}; // class congruence_domain

} // namespace ikos

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
struct abstract_domain_traits<ikos::congruence_domain<Number, VariableName>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

} // namespace domains
} // namespace crab
