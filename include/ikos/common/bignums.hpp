/*******************************************************************************
 * Implementation of bignums based on GMP, the Gnu Multiple Precision
 * Arithmetic Library (http://gmplib.org).
 ******************************************************************************/

#ifndef IKOS_BIGNUMS_HPP
#define IKOS_BIGNUMS_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <gmpxx.h>
#include <ikos/common/types.hpp>

namespace ikos {
  
  class z_number: public writeable {
    
    friend class q_number;

  private:
    mpz_class _n;
    
  private:
    z_number();
    z_number(mpz_class n): _n(n) { }
    
  public:
    z_number(std::string s) { 
      try {
	this->_n = s;
      }
      catch (std::invalid_argument& e) {
	std::ostringstream buf;
	buf << "Integer bignum: invalid string in constructor '" << s << "'";
	throw error(buf.str());
      }
    }

    z_number(int n): _n(n) { }

    z_number operator+(z_number x) {
      mpz_class r = this->_n + x._n;
      return z_number(r);
    }

    z_number operator*(z_number x) {
      mpz_class r = this->_n * x._n;
      return z_number(r);
    }

    z_number operator-(z_number x) {
      mpz_class r = this->_n - x._n;
      return z_number(r);
    }

    z_number operator-() {
      mpz_class r = -this->_n;
      return z_number(r);
    }

    z_number operator/(z_number x) {
      mpz_class r = this->_n / x._n;
      return z_number(r);
    }

    z_number operator%(z_number x) {
      mpz_class r = this->_n % x._n;
      return z_number(r);
    }
    
    z_number& operator+=(z_number x) {
      this->_n += x._n;
      return *this;
    }
    
    z_number& operator*=(z_number x) {
      this->_n *= x._n;
      return *this;
    }
    
    z_number& operator-=(z_number x) {
      this->_n -= x._n;
      return *this;
    }
    
    z_number& operator/=(z_number x) {
      this->_n /= x._n;
      return *this;
    }

    z_number& operator%=(z_number x) {
      this->_n %= x._n;
      return *this;
    }

    z_number& operator--() {
      --(this->_n);
      return *this;
    }
    
    z_number& operator++() {
      ++(this->_n);
      return *this;
    }

    z_number operator++(int) {
      z_number r(*this);
      ++(*this);
      return r;
    }

    z_number operator--(int) {
      z_number r(*this);
      --(*this);
      return r;
    }
    
    bool operator==(z_number x) {
      return (this->_n == x._n);
    }

    bool operator!=(z_number x) {
      return (this->_n != x._n);
    }

    bool operator<(z_number x) {
      return (this->_n < x._n);
    }

    bool operator<=(z_number x) {
      return (this->_n <= x._n);
    }

    bool operator>(z_number x) {
      return (this->_n > x._n);
    }

    bool operator>=(z_number x) {
      return (this->_n >= x._n);
    }

    z_number operator&(z_number x) {    
      return z_number(this->_n & x._n);
    }

    z_number operator|(z_number x) {    
      return z_number(this->_n | x._n);
    }

    z_number operator^(z_number x) {    
      return z_number(this->_n ^ x._n);
    }

    z_number operator<<(z_number x) {    
      mpz_t tmp;
      mpz_init(tmp);
      mpz_mul_2exp(tmp, this->_n.get_mpz_t(), mpz_get_ui(x._n.get_mpz_t()));
      mpz_class result(tmp);
      return z_number(result);
    }

    z_number operator>>(z_number x) {    
      mpz_class tmp(this->_n);
      return z_number(tmp.operator>>=(mpz_get_ui(x._n.get_mpz_t())));
    }

    std::ostream& write(std::ostream& o) {
      o << this->_n;
      return o;
    }
    
  }; // class z_number

  class q_number: public writeable {
    
  private:
    mpq_class _n;
    
  private:
    q_number();
    q_number(mpq_class n): _n(n) { }
    
  public:
    q_number(std::string s) {
      try {
	this->_n = s;
	this->_n.canonicalize();
      }
      catch (std::invalid_argument& e) {
	std::ostringstream buf;
	buf << "Rational bignum: invalid string in constructor '" << s << "'";
	throw error(buf.str());
      }
    }
    
    q_number(int n): _n(n) { 
      this->_n.canonicalize();
    }
    
    q_number(z_number n): _n(n._n) { 
      this->_n.canonicalize();
    }
    
    q_number(z_number n, z_number d): _n(n._n, d._n) { 
      this->_n.canonicalize();
    }

    q_number operator+(q_number x) {
      mpq_class r = this->_n + x._n;
      return q_number(r);
    }

    q_number operator*(q_number x) {
      mpq_class r = this->_n * x._n;
      return q_number(r);
    }

    q_number operator-(q_number x) {
      mpq_class r = this->_n - x._n;
      return q_number(r);
    }

    q_number operator-() {
      mpq_class r = -this->_n;
      return q_number(r);
    }

    q_number operator/(q_number x) {
      mpq_class r = this->_n / x._n;
      return q_number(r);
    }

    q_number& operator+=(q_number x) {
      this->_n += x._n;
      return *this;
    }
    
    q_number& operator*=(q_number x) {
      this->_n *= x._n;
      return *this;
    }
    
    q_number& operator-=(q_number x) {
      this->_n -= x._n;
      return *this;
    }
    
    q_number& operator/=(q_number x) {
      this->_n /= x._n;
      return *this;
    }

    q_number& operator--() {
      --(this->_n);
      return *this;
    }
    
    q_number& operator++() {
      ++(this->_n);
      return *this;
    }

    q_number operator--(int) {
      q_number r(*this);
      --(*this);
      return r;
    }
    
    q_number operator++(int) {
      q_number r(*this);
      ++(*this);
      return r;      
    }

    bool operator==(q_number x) {
      return (this->_n == x._n);
    }

    bool operator!=(q_number x) {
      return (this->_n != x._n);
    }

    bool operator<(q_number x) {
      return (this->_n < x._n);
    }

    bool operator<=(q_number x) {
      return (this->_n <= x._n);
    }

    bool operator>(q_number x) {
      return (this->_n > x._n);
    }

    bool operator>=(q_number x) {
      return (this->_n >= x._n);
    }

    z_number numerator() {
      return z_number(this->_n.get_num());
    }
  
    z_number denominator() {
      return z_number(this->_n.get_den());
    }
    
    z_number round_to_upper() {
      z_number num = numerator();
      z_number den = denominator();
      z_number q = num / den;
      z_number r = num % den;
      if (r == 0 || *this < 0) {
	return q;
      } else {
	return (q + 1);
      }
    }
  
    z_number round_to_lower() {
      z_number num = numerator();
      z_number den = denominator();
      z_number q = num / den;
      z_number r = num % den;
      if (r == 0 || *this > 0) {
	return q;
      } else {
	return (q - 1);
      }
    }    
    
    std::ostream& write(std::ostream& o) {
      o << this->_n;
      return o;
    }
    
  }; // class q_number
  
}

#endif // IKOS_BIGNUMS_HPP

