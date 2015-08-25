/*******************************************************************************
 *
 * Standard domain of constants.
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

#ifndef IKOS_CONSTANTS_HPP
#define IKOS_CONSTANTS_HPP

#include <iostream>
#include <boost/optional.hpp>
#include <ikos/common/types.hpp>
#include <ikos/common/bignums.hpp>
#include <ikos/domains/intervals.hpp>

namespace ikos {

  template< typename Number >
  class constant: public writeable {
    
  public:
    typedef constant< Number > constant_t;

  private:
    typedef enum {
      BOTTOM,
      TOP,
      NUMBER
    } kind_t;

  private:
    kind_t _kind;
    Number _n;
    
  private:
    constant();

    constant(bool b): _kind(TOP), _n(0) {
      if(b) {
	this->_kind = BOTTOM;
      }
    }

    constant(kind_t kind, Number n): _kind(kind), _n(n) { }
    
  public:
    static constant_t bottom() {
      return constant_t(true);
    }
    
    static constant_t top() {
      return constant_t(false);
    }
    
  public:
    constant(Number n): _kind(NUMBER), _n(n) { }
    
    constant(int n): _kind(NUMBER), _n(n) { }

    constant(const constant_t& other): writeable(), _kind(other._kind), _n(other._n) { }

    constant_t& operator=(constant_t other) {
      this->_kind = other._kind;
      this->_n = other._n;
      return *this;
    }

    bool is_bottom() {
      return (this->_kind == BOTTOM);
    }

    bool is_top() {
      return (this->_kind == TOP);
    }

    bool is_number() {
      return (this->_kind = NUMBER);
    }

    boost::optional< Number > number() {
      if (this->_kind == NUMBER) {
	return boost::optional< Number >(this->_n);
      } else {
	return boost::optional< Number >();
      }
    }
    
    constant_t operator-() {
      return constant_t(this->_kind, -this->_n);
    }
    
    constant_t operator+(constant_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else if (this->is_top() || other.is_top()) {
	return top();
      } else {
	return constant_t(this->_n + other._n);
      }
    }

    constant_t operator-(constant_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else if (this->is_top() || other.is_top()) {
	return top();
      } else {
	return constant_t(this->_n - other._n);
      }
    }

    constant_t operator*(constant_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else if (this->is_top() || other.is_top()) {
	return top();
      } else {
	return constant_t(this->_n * other._n);
      }
    }

    constant_t operator/(constant_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else if (this->is_top() || other.is_top()) {
	return top();
      } else if (other._n == 0) {
	return bottom();
      } else {
	return constant_t(this->_n / other._n);
      }
    }

    constant_t& operator+=(constant_t other) {
      return this->operator=(this->operator+(other));
    }
    
    constant_t& operator-=(constant_t other) {
      return this->operator=(this->operator-(other));
    }

    constant_t& operator*=(constant_t other) {
      return this->operator=(this->operator*(other));
    }
    
    constant_t& operator/=(constant_t other) {
      return this->operator=(this->operator/(other));
    }

    bool operator<=(constant_t other) {
      if (this->_kind == BOTTOM) {
	return true;
      } else if (this->_kind == TOP) {
	return (other._kind == TOP);
      } else if (this->_kind == NUMBER) {
	return (other._kind == NUMBER) && (this->_n == other._n);
      } else {
	return false;
      }
    }

    bool operator==(constant_t other) {
      if (this->_kind == NUMBER) {
	return (other._kind == NUMBER) && (this->_n == other._n);
      } else {
	return (this->_kind == other._kind);
      }
    }
    
    constant_t operator|(constant_t other) {
      if (this->is_bottom()) {
	return other;
      } else if (other.is_bottom()) {
	return *this;
      } else if (this->_kind == NUMBER && other._kind == NUMBER && this->_n == other._n) {
	return *this;
      } else {
	return top();
      }
    }
    
    constant_t operator||(constant_t other) {
      return this->operator|(other);
    }
    
    constant_t operator&(constant_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else if (this->_kind == NUMBER && other._kind == NUMBER && this->_n == other._n) {
	return *this;
      } else {
	return bottom();
      }
    }

    constant_t operator&&(constant_t other) {
      return this->operator&(other);
    }

    std::ostream& write(std::ostream& o) {
      switch (this->_kind) {
      case BOTTOM: {
	o << "_|_";
	break;
      }
      case TOP: {
	o << "T";
	break;
      }
      case NUMBER: {
	o << this->_n;
	break;
      }
      }
      return o;
    }

  }; // class constant

  template< typename Number >
  inline constant< Number > operator+(Number c, constant< Number > x) {
    return constant< Number >(c) + x;
  }

  template< typename Number >
  inline constant< Number > operator+(constant< Number > x, Number c) {
    return x + constant< Number >(c);
  }

  template< typename Number >
  inline constant< Number > operator*(Number c, constant< Number > x) {
    return constant< Number >(c) * x;
  }

  template< typename Number >
  inline constant< Number > operator*(constant< Number > x, Number c) {
    return x * constant< Number >(c);
  }

  template< typename Number >
  inline constant< Number > operator/(Number c, constant< Number > x) {
    return constant< Number >(c) / x;
  }

  template< typename Number >
  inline constant< Number > operator/(constant< Number > x, Number c) {
    return x / constant< Number >(c);
  }

  template< typename Number >
  inline constant< Number > operator-(Number c, constant< Number > x) {
    return constant< Number >(c) - x;
  }

  template< typename Number >
  inline constant< Number > operator-(constant< Number > x, Number c) {
    return x - constant< Number >(c);
  }

  typedef constant< z_number > z_constant;
  typedef constant< q_number > q_constant;

  template< typename Number, typename VariableName, std::size_t max_reduction_cycles = 10 >
  class constant_domain: public writeable, numerical_domain< Number, VariableName > {

  public:
    typedef constant< Number > constant_t;
    typedef constant_domain< Number, VariableName > constant_domain_t;
    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    
  private:
    typedef separate_domain< VariableName, constant_t > separate_domain_t;
    typedef interval< Number > interval_t;
    typedef interval_domain< Number, VariableName > interval_domain_t;
    typedef typename linear_constraint_system_t::variable_set_t variable_set_t;

  public:
    typedef typename separate_domain_t::iterator iterator;
     
  private:
    separate_domain_t _env;

  private:
    constant_domain(separate_domain_t env): _env(env) { }
    
  public:
    static constant_domain_t top() {
      return constant_domain(separate_domain_t::top());
    }
    
    static constant_domain_t bottom() {
      return constant_domain(separate_domain_t::bottom());
    }

  public:
    constant_domain(): _env(separate_domain_t::top()) { }

    constant_domain(const constant_domain_t& e): writeable(), numerical_domain< Number, VariableName >(), _env(e._env) { }

    constant_domain_t& operator=(constant_domain_t e) {
      this->_env = e._env;
      return *this;
    }
    
    iterator begin() {
      return this->_env.begin();
    }

    iterator end() {
      return this->_env.end();
    }

    bool is_bottom() {
      return this->_env.is_bottom();
    }

    bool is_top() {
      return this->_env.is_top();
    }

    bool operator<=(constant_domain_t e) {
      return (this->_env <= e._env);
    }

    constant_domain_t operator|(constant_domain_t e) {
      return (this->_env | e._env);
    }

    constant_domain_t operator&(constant_domain_t e) {
      return (this->_env & e._env);
    }

    constant_domain_t operator||(constant_domain_t e) {
      return (this->_env || e._env);
    }

    constant_domain_t operator&&(constant_domain_t e) {
      return (this->_env && e._env);
    }

    void set(VariableName v, constant_t i) {
      this->_env.set(v, i);
    }

    void set(VariableName v, Number n) {
      this->_env.set(v, constant_t(n));
    }

    void operator-=(VariableName v) {
      this->_env -= v;
    }
    
    constant_t operator[](VariableName v) {
      return _env[v];
    }

    constant_t operator[](linear_expression_t expr) {
      constant_t r(expr.constant());
      for (typename linear_expression_t::iterator it = expr.begin(); it != expr.end(); ++it) {
	constant_t c(it->first);
	r += c * this->_env[it->second.name()];
      }
      return r;
    }
    
    void operator+=(linear_constraint_system_t csts) {
      this->add(csts);
    }
    
    void add(linear_constraint_system_t csts, std::size_t threshold = max_reduction_cycles) {
      if (!this->is_bottom()) {
	interval_domain_t e;
	variable_set_t vars = csts.variables();
	for (typename variable_set_t::iterator it = vars.begin(); it != vars.end(); ++it) {
	  VariableName n = it->name();
	  constant_t c = this->operator[](n);
	  boost::optional< Number > v = c.number();
	  if (v) {
	    e.set(n, *v);
	  }
	}
	e.add(csts, threshold);
	for (typename variable_set_t::iterator it = vars.begin(); it != vars.end(); ++it) {
	  VariableName n = it->name();
	  interval_t i = e[n];
	  if (i.is_bottom()) {
	    this->set(n, constant_t::bottom());
	    return;
	  } else {
	    boost::optional< Number > v = i.singleton();
	    if (v) {
	      this->set(n, *v);
	    }
	  }
	}
      }
    }
    
    constant_domain_t operator+(linear_constraint_system_t csts) {
      constant_domain_t e(this->_env);
      e += csts;
      return e;
    }

    void assign(VariableName x, linear_expression_t e) {
      constant_t r = e.constant();
      for (typename linear_expression_t::iterator it = e.begin(); it != e.end(); ++it) {
	r += it->first * this->_env[it->second.name()];
      }
      this->_env.set(x, r);
    }
    
    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
      constant_t yi = this->_env[y];
      constant_t zi = this->_env[z];
      constant_t xi = constant_t::bottom();
      
      switch (op) {
      case OP_ADDITION: {
	xi = yi + zi;
	break;
      }
      case OP_SUBTRACTION: {
	xi = yi - zi;
	break;
      }
      case OP_MULTIPLICATION: {
	xi = yi * zi;
	break;
      }
      case OP_DIVISION: {
	xi = yi / zi;
	break;
      }
      }
      this->_env.set(x, xi);
    }

    void apply(operation_t op, VariableName x, VariableName y, Number k) {
      constant_t yi = this->_env[y];
      constant_t zi(k);
      constant_t xi = constant_t::bottom();
      
      switch (op) {
      case OP_ADDITION: {
	xi = yi + zi;
	break;
      }
      case OP_SUBTRACTION: {
	xi = yi - zi;
	break;
      }
      case OP_MULTIPLICATION: {
	xi = yi * zi;
	break;
      }
      case OP_DIVISION: {
	xi = yi / zi;
	break;
      }
      }
      this->_env.set(x, xi);
    }
    
    std::ostream& write(std::ostream& o) {
      return this->_env.write(o);
    }

    const char* getDomainName () const {return "Constant Propagation";}

  }; // class constant_domain

} // namespace ikos

#endif // IKOS_CONSTANTS_HPP
