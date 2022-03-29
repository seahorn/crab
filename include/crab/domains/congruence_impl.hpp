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

#include <crab/domains/congruence.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>

namespace ikos {

/* Notes about the % operator
 *
 * The semantics of r = n % d is to set r to "n mod d". The sign of
 * the d is ignored and r is always non-negative.
 *
 * We assume that n % d (also n /d) raises a runtime error if d==0.
 */

template <typename Number> void congruence<Number>::normalize(void) {
  // Set to standard form: 0 <= b < a for a != 0
  if (m_a != 0) {
    m_b = m_b % m_a;
  }
}

// if true then top (1Z + 0) else bottom
template <typename Number>
congruence<Number>::congruence(bool b) : m_is_bottom(!b), m_a(1), m_b(0) {}

template <typename Number>
congruence<Number>::congruence(int n) : m_is_bottom(false), m_a(0), m_b(n) {}

template <typename Number>
congruence<Number>::congruence(Number a, Number b)
    : m_is_bottom(false), m_a(a), m_b(b) {
  normalize();
}

template <typename Number> Number congruence<Number>::abs(Number x) const {
  return x < 0 ? -x : x;
}

template <typename Number>
Number congruence<Number>::max(Number x, Number y) const {
  return (x <= y) ? y : x;
}

template <typename Number>
Number congruence<Number>::min(Number x, Number y) const {
  return (x < y) ? x : y;
}

template <typename Number>
Number congruence<Number>::gcd(Number x, Number y, Number z) const {
  return gcd(x, gcd(y, z));
}

// Not to be called explicitly outside of gcd
template <typename Number>
Number congruence<Number>::gcd_helper(Number x, Number y) const {
  return (y == 0) ? x : gcd_helper(y, x % y);
}

template <typename Number>
Number congruence<Number>::gcd(Number x, Number y) const {
  return gcd_helper(abs(x), abs(y));
}

template <typename Number>
Number congruence<Number>::lcm(Number x, Number y) const {
  Number tmp = gcd(x, y);
  return abs(x * y) / tmp;
}

template <typename Number> bool congruence<Number>::is_zero() const {
  return !is_bottom() && m_a == 0 && m_b == 0;
}

template <typename Number> bool congruence<Number>::all_ones() const {
  return !is_bottom() && m_a == 0 && m_b == -1;
}

template <typename Number> congruence<Number> congruence<Number>::top() {
  return congruence(true);
}

template <typename Number> congruence<Number> congruence<Number>::bottom() {
  return congruence(false);
}

template <typename Number>
congruence<Number>::congruence() : m_is_bottom(false), m_a(1), m_b(0) {}

template <typename Number>
congruence<Number>::congruence(Number n) : m_is_bottom(false), m_a(0), m_b(n) {}

template <typename Number> bool congruence<Number>::is_bottom() const {
  return m_is_bottom;
}

template <typename Number> bool congruence<Number>::is_top() const {
  return m_a == 1;
}

template <typename Number>
boost::optional<Number> congruence<Number>::singleton() const {
  if (!this->is_bottom() && m_a == 0) {
    return boost::optional<Number>(m_b);
  } else {
    return boost::optional<Number>();
  }
}

template <typename Number> Number congruence<Number>::get_modulo() const {
  return m_a;
}

template <typename Number> Number congruence<Number>::get_remainder() const {
  return m_b;
}

template <typename Number>
bool congruence<Number>::operator==(const congruence<Number> &o) const {
  return (is_bottom() == o.is_bottom() && m_a == o.m_a && m_b == o.m_b);
}

template <typename Number>
bool congruence<Number>::operator!=(const congruence<Number> &x) const {
  return !this->operator==(x);
}

template <typename Number>
bool congruence<Number>::operator<=(const congruence<Number> &o) const {
  if (is_bottom()) {
    return true;
  } else if (o.is_bottom()) {
    return false;
  } else if (m_a == 0 && o.m_a == 0) {
    return (m_b == o.m_b);
  } else if (m_a == 0) {
    if ((m_b % o.m_a) == (o.m_b % o.m_a)) {
      return true;
    }
  } else if (o.m_a == 0) {
    if (m_b % m_a == (o.m_b % m_a)) {
      return false;
    }
  }
  return (m_a % o.m_a == 0) && (m_b % o.m_a == o.m_b % o.m_a);
}

template <typename Number>
congruence<Number>
congruence<Number>::operator|(const congruence<Number> &o) const {
  if (is_bottom()) {
    return o;
  } else if (o.is_bottom()) {
    return *this;
  } else if (is_top() || o.is_top()) {
    return top();
  } else {
    return congruence<Number>(gcd(m_a, o.m_a, abs(m_b - o.m_b)),
                              min(m_b, o.m_b));
  }
}

template <typename Number>
congruence<Number>
congruence<Number>::operator&(const congruence<Number> &o) const {
  if (is_bottom() || o.is_bottom()) {
    return bottom();
  }

  // lcm has meaning only if both a and o.a are not 0
  if (m_a == 0 && o.m_a == 0) {
    if (m_b == o.m_b) {
      return *this;
    } else {
      return bottom();
    }
  } else if (m_a == 0) {
    // b & a'Z + b' iff \exists k such that a'*k + b' = b iff ((b - b') %a' ==
    // 0)
    if ((m_b - o.m_b) % o.m_a == 0) {
      return *this;
    } else {
      return bottom();
    }
  } else if (o.m_a == 0) {
    // aZ+b & b' iff \exists k such that a*k+b  = b' iff ((b'-b %a) == 0)
    if ((o.m_b - m_b) % m_a == 0) {
      return o;
    } else {
      return bottom();
    }
  } else {
    // pre: a and o.a != 0
    Number x = gcd(m_a, o.m_a);
    if (m_b % x == (o.m_b % x)) {
      // the part max(b,o.b) needs to be verified. What we really
      // want is to find b'' such that
      // 1) b'' % lcm(a,a') == b  % lcm(a,a'), and
      // 2) b'' % lcm(a,a') == b' % lcm(a,a').
      // An algorithm for that is provided in Granger'89.
      return congruence<Number>(lcm(m_a, o.m_a), max(m_b, o.m_b));
    } else {
      return congruence<Number>::bottom();
    }
  }
}

template <typename Number>
congruence<Number>
congruence<Number>::operator||(const congruence<Number> &o) const {
  // Equivalent to join, domain is flat
  return *this | o;
}

template <typename Number>
congruence<Number>
congruence<Number>::operator&&(const congruence<Number> &o) const {
  // Simply refines top element
  return (is_top()) ? o : *this;
}

/** Arithmetic Operators **/
template <typename Number>
congruence<Number>
congruence<Number>::operator+(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else
    return congruence<Number>(gcd(m_a, o.m_a), m_b + o.m_b);
}

template <typename Number>
congruence<Number>
congruence<Number>::operator-(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else
    return congruence<Number>(gcd(m_a, o.m_a), m_b - o.m_b);
}

template <typename Number>
congruence<Number> congruence<Number>::operator-() const {
  if (this->is_bottom() || this->is_top()) {
    return *this;
  } else {
    return congruence<Number>(m_a, -m_b + m_a);
  }
}

template <typename Number>
congruence<Number>
congruence<Number>::operator*(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if ((this->is_top() || o.is_top()) && m_a != 0 && o.m_a != 0)
    return congruence<Number>::top();
  else
    return congruence<Number>(gcd(m_a * o.m_a, m_a * o.m_b, o.m_a * m_b),
                              m_b * o.m_b);
}

// signed division
template <typename Number>
congruence<Number>
congruence<Number>::operator/(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (o == congruence(0))
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    /*
       aZ+b / 0Z+b':
          if b'|a then  (a/b')Z + b/b'
          else          top
    */
    if (o.m_a == 0) {
      if (m_a % o.m_b == 0)
        return congruence<Number>(m_a / o.m_b, m_b / o.m_b);
      else
        return congruence<Number>::top();
    }

    /*
         0Z+b / a'Z+b':
            if N>0   (b div N)Z + 0
            else     0Z + 0

           where N = a'((b-b') div a') + b'
    */
    if (m_a == 0) {
      Number n(o.m_a * (((m_b - o.m_b) / o.m_a) + o.m_b));
      if (n > 0) {
        return congruence<Number>(m_b / n, Number(0));
      } else {
        return congruence<Number>(Number(0), Number(0));
      }
    }

    /*
      General case: no singleton
    */
    return congruence<Number>::top();
  }
}

// signed remainder operator
template <typename Number>
congruence<Number>
congruence<Number>::operator%(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (o == congruence(0))
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    /*
         aZ+b mod 0Z+b':
             if b'|a then  (a/b')Z + b/b'
             else          top
    */
    if (o.m_a == 0) {
      if (m_a % o.m_b == 0) {
        return congruence<Number>(Number(0), m_b % o.m_b);
      } else {
        return congruence<Number>(gcd(m_a, o.m_b), m_b);
      }
    }
    /*
          0Z+b mod a'Z+b':
           if N<=0           then 0Z+b
           if (b div N) == 1 then gcd(b',a')Z + b
           if (b div N) >= 2 then N(b div N)Z  + b

         where N = a'((b-b') div a') + b'
    */
    if (m_a == 0) {
      Number n(o.m_a * (((m_b - o.m_b) / o.m_a) + o.m_b));
      if (n <= 0) {
        return congruence<Number>(m_a, m_b);
      } else if (m_b == n) {
        return congruence<Number>(gcd(o.m_b, o.m_a), m_b);
      } else if ((m_b / n) >= 2) {
        return congruence<Number>(m_b, m_b);
      } else {
        CRAB_ERROR("unreachable");
      }
    }

    /*
      general case: no singleton
    */
    return congruence<Number>(gcd(m_a, o.m_a, o.m_b), m_b);
  }
}

template <typename Number>
congruence<Number> congruence<Number>::And(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    if (is_zero() || o.is_zero()) {
      return congruence<Number>(0);
    } else if (all_ones()) {
      return o;
    } else if (o.all_ones()) {
      return *this;
    } else if (m_a == 0 && o.m_a == 0) {
      return congruence<Number>(m_b & o.m_b);
    } else {
      return top();
    }
  }
}

template <typename Number>
congruence<Number> congruence<Number>::Or(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    if (all_ones() || o.all_ones()) {
      return congruence<Number>(-1);
    } else if (is_zero()) {
      return o;
    } else if (o.is_zero()) {
      return *this;
    } else if (m_a == 0 && o.m_a == 0) {
      return congruence<Number>(m_b | o.m_b);
    } else {
      return top();
    }
  }
}

template <typename Number>
congruence<Number> congruence<Number>::Xor(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    if (is_zero()) {
      return o;
    } else if (o.is_zero()) {
      return *this;
    } else if (m_a == 0 && o.m_a == 0) {
      return congruence<Number>(m_b ^ o.m_b);
    } else {
      return top();
    }
  }
}

template <typename Number>
congruence<Number> congruence<Number>::Shl(const congruence<Number> &o) const {
  if (this->is_bottom() || o.is_bottom())
    return congruence<Number>::bottom();
  else if (this->is_top() || o.is_top())
    return congruence<Number>::top();
  else {
    if (o.m_a == 0) { // singleton
      if (o.m_b < 0) {
        return bottom();
      }
      // aZ + b << 0Z + b'  = (a*2^b')Z + b*2^b'
      Number x = Number(1) << o.m_b;
      return congruence<Number>(m_a * x, m_b * x);
    } else {
      Number x = Number(1) << o.m_b;
      Number y = Number(1) << o.m_a;
      // aZ + b << a'Z + b' = (gcd(a, b * (2^a' - 1)))*(2^b')Z + b*(2^b')
      return congruence<Number>(gcd(m_a, m_b * (y - 1)) * x, m_b * x);
    }
  }
}

template <typename Number>
congruence<Number> congruence<Number>::AShr(const congruence<Number> &o) const {
  auto to_interval = [](const congruence<Number> &c) -> interval<Number> {
    assert(c.singleton());
    return interval<Number>(*(c.singleton()));
  };

  if (this->is_bottom() || o.is_bottom()) {
    return congruence<Number>::bottom();
  } else if (this->is_top() || o.is_top()) {
    return congruence<Number>::top();
  } else {
    if (o.m_a == 0) { // singleton
      // aZ + b >> 0Z + b'
      if (o.m_b < 0) {
        return congruence<Number>::bottom();
      }
    }

    if (singleton() && o.singleton()) {
      auto res = to_interval(*this).AShr(to_interval(o));
      if (boost::optional<Number> n = res.singleton()) {
        return congruence(*n);
      }
    }
    return congruence<Number>::top();
  }
}

template <typename Number>
congruence<Number> congruence<Number>::LShr(const congruence<Number> &o) const {
  auto to_interval = [](const congruence<Number> &c) -> interval<Number> {
    assert(c.singleton());
    return interval<Number>(*(c.singleton()));
  };

  if (this->is_bottom() || o.is_bottom()) {
    return congruence<Number>::bottom();
  } else if (this->is_top() || o.is_top()) {
    return congruence<Number>::top();
  } else {
    if (o.m_a == 0) {
      // aZ + b >> 0Z + b'
      if (o.m_b < 0) {
        return congruence<Number>::bottom();
      }
    }
    if (singleton() && o.singleton()) {
      auto res = to_interval(*this).LShr(to_interval(o));
      if (boost::optional<Number> n = res.singleton()) {
        return congruence(*n);
      }
    }
    return congruence<Number>::top();
  }
}

// division and remainder operations
template <typename Number>
congruence<Number> congruence<Number>::SDiv(const congruence<Number> &x) const {
  return this->operator/(x);
}

template <typename Number>
congruence<Number> congruence<Number>::UDiv(const congruence<Number> &x) const {
  return congruence<Number>::top();
}

template <typename Number>
congruence<Number> congruence<Number>::SRem(const congruence<Number> &x) const {
  return this->operator%(x);
}

template <typename Number>
congruence<Number> congruence<Number>::URem(const congruence<Number> &x) const {
  return congruence<Number>::top();
}

template <typename Number>
void congruence<Number>::write(crab::crab_os &o) const {
  if (is_bottom()) {
    o << "_|_";
    return;
  }

  if (m_a == 0) {
    o << m_b;
    return;
  }

  o << m_a << "Z+" << m_b;
}

} // namespace ikos
