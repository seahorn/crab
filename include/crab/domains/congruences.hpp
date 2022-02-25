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

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/variable.hpp>

#include <boost/optional.hpp>

namespace ikos {

template <typename Number> class congruence {
  using interval_t = interval<Number>;

public:
  using congruence_t = congruence<Number>;

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

  bool operator==(const congruence_t &o) const {
    return (is_bottom() == o.is_bottom() && _a == o._a && _b == o._b);
  }

  bool operator!=(const congruence_t &x) const { return !this->operator==(x); }

  /** Lattice Operations **/

  bool operator<=(const congruence_t &o) const {
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

  congruence_t operator|(const congruence_t &o) const {
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

  congruence_t operator&(const congruence_t &o) const {
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

  congruence_t operator||(const congruence_t &o) const {
    // Equivalent to join, domain is flat
    return *this | o;
  }

  congruence_t operator&&(const congruence_t &o) const {
    // Simply refines top element
    return (is_top()) ? o : *this;
  }

  /** Arithmetic Operators **/

  congruence_t operator+(const congruence_t &o) const {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else
      return congruence_t(gcd(_a, o._a), _b + o._b);
  }

  congruence_t operator-(const congruence_t &o) const {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if (this->is_top() || o.is_top())
      return congruence_t::top();
    else
      return congruence_t(gcd(_a, o._a), _b - o._b);
  }

  congruence_t operator-() const {
    if (this->is_bottom() || this->is_top())
      return *this;
    else
      return congruence_t(_a, -_b + _a);
  }

  congruence_t operator*(const congruence_t &o) const {
    if (this->is_bottom() || o.is_bottom())
      return congruence_t::bottom();
    else if ((this->is_top() || o.is_top()) && _a != 0 && o._a != 0)
      return congruence_t::top();
    else
      return congruence_t(gcd(_a * o._a, _a * o._b, o._a * _b), _b * o._b);
  }

  // signed division
  congruence_t operator/(const congruence_t &o) const {
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
  congruence_t operator%(const congruence_t &o) const {
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

  congruence_t And(const congruence_t &o) const {
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

  congruence_t Or(const congruence_t &o) const {
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

  congruence_t Xor(const congruence_t &o) const {
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

  congruence_t Shl(const congruence_t &o) const {
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

  congruence_t AShr(const congruence_t &o) const {
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

  congruence_t LShr(const congruence_t &o) const {
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

  congruence_t SDiv(const congruence_t &x) const { return this->operator/(x); }

  congruence_t UDiv(const congruence_t &x) const { return congruence_t::top(); }

  congruence_t SRem(const congruence_t &x) const { return this->operator%(x); }

  congruence_t URem(const congruence_t &x) const { return congruence_t::top(); }

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
inline congruence<Number> operator+(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) + x;
}

template <typename Number>
inline congruence<Number> operator+(const congruence<Number> &x, Number c) {
  return x + congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator*(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) * x;
}

template <typename Number>
inline congruence<Number> operator*(const congruence<Number> &x, Number c) {
  return x * congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator/(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) / x;
}

template <typename Number>
inline congruence<Number> operator/(const congruence<Number> &x, Number c) {
  return x / congruence<Number>(c);
}

template <typename Number>
inline congruence<Number> operator-(Number c, const congruence<Number> &x) {
  return congruence<Number>(c) - x;
}

template <typename Number>
inline congruence<Number> operator-(const congruence<Number> &x, Number c) {
  return x - congruence<Number>(c);
}

template <typename Number, typename VariableName, typename CongruenceCollection>
class equality_congruence_solver {
  // TODO: check correctness of the solver. Granger provides a sound
  // and more precise solver for equality linear congruences (see
  // Theorem 4.4).
private:
  using congruence_t = congruence<Number>;
  using variable_t = crab::variable<Number, VariableName>;
  using linear_expression_t = linear_expression<Number, VariableName>;
  using linear_constraint_t = linear_constraint<Number, VariableName>;
  using linear_constraint_system_t =
      linear_constraint_system<Number, VariableName>;

  using cst_table_t = std::vector<linear_constraint_t>;
  using variable_set_t = std::set<variable_t>;

  std::size_t m_max_cycles;
  bool m_is_contradiction;
  cst_table_t m_cst_table;
  variable_set_t m_refined_variables;
  std::size_t m_op_count;

private:
  bool refine(const variable_t &v, congruence_t i, CongruenceCollection &env) {
    congruence_t old_i = env.at(v);
    congruence_t new_i = old_i & i;
    if (new_i.is_bottom()) {
      return true;
    }
    if (old_i != new_i) {
      env.set(v, new_i);
      m_refined_variables.insert(v);
      ++(m_op_count);
    }
    return false;
  }

  congruence_t compute_residual(const linear_constraint_t &cst,
                                const variable_t &pivot,
                                CongruenceCollection &env) {
    congruence_t residual(cst.constant());
    for (auto kv : cst) {
      const variable_t &v = kv.second;
      if (!(v == pivot)) {
        residual = residual - (kv.first * env.at(v));
        ++(m_op_count);
      }
    }
    return residual;
  }

  bool propagate(const linear_constraint_t &cst, CongruenceCollection &env) {
    for (auto kv : cst) {
      Number c = kv.first;
      const variable_t &pivot = kv.second;
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
      for (const linear_constraint_t &cst : m_cst_table) {
        if (propagate(cst, env)) {
          return true;
        }
      }
    } while (m_refined_variables.size() > 0 && cycle <= m_max_cycles);
    return false;
  }

public:
  equality_congruence_solver(const linear_constraint_system_t &csts,
                             std::size_t max_cycles)
      : m_max_cycles(max_cycles), m_is_contradiction(false) {
    for (auto const &cst : csts) {
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
class congruence_domain final : public crab::domains::abstract_domain_api<
                                    congruence_domain<Number, VariableName>> {
public:
  using congruence_t = congruence<Number>;

private:
  // note that this is assuming that all variables have the same bit
  // width which is unrealistic.
  using congruence_domain_t = congruence_domain<Number, VariableName>;
  using abstract_domain_t =
      crab::domains::abstract_domain_api<congruence_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using number_t = Number;
  using varname_t = VariableName;
  using typename abstract_domain_t::reference_constraint_t;

private:
  using separate_domain_t = separate_domain<variable_t, congruence_t>;
  using solver_t =
      equality_congruence_solver<number_t, varname_t, separate_domain_t>;

public:
  using iterator = typename separate_domain_t::iterator;

private:
  separate_domain_t _env;

private:
  congruence_domain(separate_domain_t env) : _env(env) {}

public:
  congruence_domain_t make_top() const override {
    return congruence_domain_t(separate_domain_t::top());
  }

  congruence_domain_t make_bottom() const override {
    return congruence_domain_t(separate_domain_t::bottom());
  }

  void set_to_top() override {
    congruence_domain abs(separate_domain_t::top());
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    congruence_domain abs(separate_domain_t::bottom());
    std::swap(*this, abs);
  }

  congruence_domain() : _env(separate_domain_t::top()) {}

  congruence_domain(const congruence_domain_t &e) : _env(e._env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  congruence_domain_t &operator=(const congruence_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o)
      this->_env = o._env;
    return *this;
  }

  iterator begin() { return this->_env.begin(); }

  iterator end() { return this->_env.end(); }

  bool is_bottom() const override { return this->_env.is_bottom(); }

  bool is_top() const override { return this->_env.is_top(); }

  bool operator<=(const congruence_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");
    return this->_env <= e._env;
  }

  void operator|=(const congruence_domain_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    this->_env = this->_env | e._env;
  }

  congruence_domain_t operator|(const congruence_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");
    return this->_env | e._env;
  }

  congruence_domain_t operator&(const congruence_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");
    return this->_env & e._env;
  }

  congruence_domain_t operator||(const congruence_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    return this->_env || e._env;
  }

  congruence_domain_t widening_thresholds(
      const congruence_domain_t &other,
      const crab::iterators::thresholds<number_t> &) const override {
    return (*this || other);
  }

  congruence_domain_t operator&&(const congruence_domain_t &e) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");
    return this->_env && e._env;
  }

  void set(const variable_t &v, congruence_t i) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    this->_env.set(v, i);
  }

  void set(const variable_t &v, number_t n) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
    this->_env.set(v, congruence_t(n));
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");
    this->_env -= v;
  }

  virtual interval_t operator[](const variable_t &v) override {
    CRAB_WARN(domain_name(), "::operator[] not implemented");
    return interval_t::top();
  }

  virtual interval_t at(const variable_t &v) const override {
    CRAB_WARN(domain_name(), "::at not implemented");
    return interval_t::top();
  }  

  congruence_t to_congruence(const variable_t &v) { return this->_env.at(v); }

  congruence_t to_congruence(const linear_expression_t &expr) {
    congruence_t r(expr.constant());
    for (auto kv : expr) {
      congruence_t c(kv.first);
      r = r + (c * to_congruence(kv.second));
    }
    return r;
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");
    const std::size_t threshold = 10;
    if (!this->is_bottom()) {
      solver_t solver(csts, threshold);
      solver.run(this->_env);
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    congruence_t r = e.constant();
    for (auto kv : e) {
      r = r + (kv.first * this->_env.at(kv.second));
    }
    this->_env.set(x, r);
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    congruence_t yi = this->_env.at(y);
    congruence_t zi = this->_env.at(z);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case crab::domains::OP_ADDITION:
      xi = yi + zi;
      break;
    case crab::domains::OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case crab::domains::OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case crab::domains::OP_SDIV:
      xi = yi / zi;
      break;
    case crab::domains::OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case crab::domains::OP_SREM:
      xi = yi.SRem(zi);
      break;
    case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    congruence_t yi = this->_env.at(y);
    congruence_t zi(k);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case crab::domains::OP_ADDITION:
      xi = yi + zi;
      break;
    case crab::domains::OP_SUBTRACTION:
      xi = yi - zi;
      break;
    case crab::domains::OP_MULTIPLICATION:
      xi = yi * zi;
      break;
    case crab::domains::OP_SDIV:
      xi = yi / zi;
      break;
    case crab::domains::OP_UDIV:
      xi = yi.UDiv(zi);
      break;
    case crab::domains::OP_SREM:
      xi = yi.SRem(zi);
      break;
    case crab::domains::OP_UREM:
      xi = yi.URem(zi);
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    this->_env.set(x, xi);
  }

  // backward operations
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const congruence_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    crab::domains::BackwardAssignOps<congruence_domain_t>::assign(*this, x, e,
                                                                  inv);
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const congruence_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<congruence_domain_t>::apply(*this, op, x,
                                                                 y, z, inv);
  }

  void backward_apply(crab::domains::arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const congruence_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<congruence_domain_t>::apply(*this, op, x,
                                                                 y, z, inv);
  }

  // cast operations

  void apply(crab::domains::int_conv_operation_t /*op*/, const variable_t &dst,
             const variable_t &src) override {
    // ignore widths
    assign(dst, src);
  }

  // bitwise operations

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    congruence_t yi = this->_env.at(y);
    congruence_t zi = this->_env.at(z);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case crab::domains::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case crab::domains::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case crab::domains::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case crab::domains::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case crab::domains::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case crab::domains::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default: { CRAB_ERROR("unreachable"); }
    }
    this->_env.set(x, xi);
  }

  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    congruence_t yi = this->_env.at(y);
    congruence_t zi(k);
    congruence_t xi = congruence_t::bottom();

    switch (op) {
    case crab::domains::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case crab::domains::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case crab::domains::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case crab::domains::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case crab::domains::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case crab::domains::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default: { CRAB_ERROR("unreachable"); }
    }
    this->_env.set(x, xi);
  }

  DEFAULT_SELECT(congruence_domain_t)
  
  /// congruence_domain implements only standard abstract operations
  /// of a numerical domain so it is intended to be used as a leaf
  /// domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(congruence_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(congruence_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(congruence_domain_t)
    
  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t var : variables) {
      this->operator-=(var);
    }
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    _env.project(variables);
  }

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    set(new_x, this->_env.at(x));
  }

  void normalize() override {}

  void minimize() override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    _env.rename(from, to);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const congruence_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void write(crab::crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    this->_env.write(o);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    for (iterator it = this->_env.begin(); it != this->_env.end(); ++it) {
      const variable_t &v = it->first;
      congruence_t c = it->second;
      boost::optional<number_t> n = c.singleton();
      if (n) {
        csts += (v == *n);
      }
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  std::string domain_name() const override { return "Congruences"; }

}; // class congruence_domain

} // namespace ikos

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
struct abstract_domain_traits<ikos::congruence_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // namespace domains
} // namespace crab
