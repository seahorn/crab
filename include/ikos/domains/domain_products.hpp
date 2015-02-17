/*******************************************************************************
 * Products of abstract domains.
 ******************************************************************************/


#ifndef IKOS_DOMAIN_PRODUCTS_HPP
#define IKOS_DOMAIN_PRODUCTS_HPP

#include <iostream>
#include <ikos/common/types.hpp>
#include <ikos/domains/numerical_domains_api.hpp>

namespace ikos {

  template< typename Domain1, typename Domain2 >
  class domain_product2: public writeable {
    
  public:
    typedef domain_product2< Domain1, Domain2 > domain_product2_t;
    
  private:
    bool _is_bottom;
    Domain1 _first;
    Domain2 _second;

  public:
    static domain_product2_t top() {
      return domain_product2_t(Domain1::top(), Domain2::top());
    }
    
    static domain_product2_t bottom() {
      return domain_product2_t(Domain1::bottom(), Domain2::bottom());
    }
    
  private:
    void canonicalize() {
      if (!this->_is_bottom) {
	this->_is_bottom = this->_first.is_bottom() || this->_second.is_bottom();
	if (this->_is_bottom) {
	  this->_first = Domain1::bottom();
	  this->_second = Domain2::bottom();      
	}
      }
    }

  public:
    domain_product2(): _is_bottom(false), _first(Domain1::top()), _second(Domain2::top()) { }
    
    domain_product2(Domain1 first, Domain2 second): _is_bottom(false), _first(first), _second(second) {
      this->canonicalize();
    }
    
    domain_product2(const domain_product2_t& other): writeable(), _is_bottom(other._is_bottom), _first(other._first), _second(other._second) { }

    domain_product2_t& operator=(domain_product2_t other) {
      this->_is_bottom = other._is_bottom;
      this->_first = other._first;
      this->_second = other._second;
      return *this;
    }
    
    bool is_bottom() {
      this->canonicalize();
      return this->_is_bottom;
    }
    
    bool is_top() {
      return (this->_first.is_top() && this->_second.is_top());
    }
    
    Domain1& first() {
      this->canonicalize();
      return this->_first;
    }
    
    Domain2& second() {
      this->canonicalize();
      return this->_second;
    }

    bool operator<=(domain_product2_t other) {
      if (this->is_bottom()) {
	return true;
      } else if (other.is_bottom()) {
	return false;
      } else {
	return (this->_first <= other._first) && (this->_second <= other._second);
      }
    }
    
    bool operator==(domain_product2_t other) {
      return (this->operator<=(other) && other.operator<=(*this));
    }
    
    domain_product2_t operator|(domain_product2_t other) {
      if (this->is_bottom()) {
	return other;
      } else if (other.is_bottom()) {
	return *this;
      } else {
	return domain_product2_t(this->_first | other._first, this->_second | other._second);
      }
    }

    domain_product2_t operator||(domain_product2_t other) {
      if (this->is_bottom()) {
	return other;
      } else if (other.is_bottom()) {
	return *this;
      } else {
	return domain_product2_t(this->_first || other._first, this->_second || other._second);
      }
    }

    domain_product2_t operator&(domain_product2_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else {
	return domain_product2_t(this->_first & other._first, this->_second & other._second);
      }
    }

    domain_product2_t operator&&(domain_product2_t other) {
      if (this->is_bottom() || other.is_bottom()) {
	return bottom();
      } else {
	return domain_product2_t(this->_first && other._first, this->_second && other._second);
      }
    }
    
    std::ostream& write(std::ostream& o) {
      if (this->is_bottom()) {
	o << "_|_";
      } else {
	o << "(" << this->_first << ", " << this->_second << ")";
      }
      return o;
    }

  }; // class domain_product2

  template< typename Domain1, typename Domain2, typename Domain3 >
  class domain_product3: public writeable {

  public:
    typedef domain_product3< Domain1, Domain2, Domain3 > domain_product3_t;

  private:
    typedef domain_product2< Domain2, Domain3 > product23_t;
    typedef domain_product2< Domain1, product23_t > product123_t;
    
  private:
    product123_t _product;

  private:
    domain_product3(product123_t product): _product(product) { }

  public:
    static domain_product3_t top() {
      return domain_product3_t(product123_t::top());
    }
    
    static domain_product3_t bottom() {
      return domain_product3_t(product123_t::bottom());
    }
    
  public:
    domain_product3(): _product() { }

    domain_product3(Domain1 first, Domain2 second, Domain3 third): _product(first, product23_t(second, third)) { }

    domain_product3(const domain_product3_t& other): writeable(), _product(other._product) { }
    
    domain_product3_t& operator=(domain_product3_t other) {
      this->_product = other._product;
      return *this;
    }

    bool is_bottom() {
      return this->_product.is_bottom();
    }

    bool is_top() {
      return this->_product.is_top();
    }

    Domain1& first() {
      return this->_product.first();
    }

    Domain2& second() {
      return this->_product.second().first();
    }

    Domain3& third() {
      return this->_product.second().second();
    }

    bool operator<=(domain_product3_t other) {
      return (this->_product <= other._product);
    }

    bool operator==(domain_product3_t other) {
      return (this->_product == other._product);
    }

    domain_product3_t operator|(domain_product3_t other) {
      return domain_product3_t(this->_product | other._product);
    }

    domain_product3_t operator&(domain_product3_t other) {
      return domain_product3_t(this->_product & other._product);
    }

    domain_product3_t operator||(domain_product3_t other) {
      return domain_product3_t(this->_product || other._product);
    }

    domain_product3_t operator&&(domain_product3_t other) {
      return domain_product3_t(this->_product && other._product);
    }

    std::ostream& write(std::ostream& o) {
      if (this->is_bottom()) {
	o << "_|_";
      } else {
	o << "(" << this->first() << ", " << this->second() << ", " << this->third() << ")";
      }
      return o;
    }
    
  }; // class domain_product3

  template< typename Number, typename VariableName, typename Domain1, typename Domain2 >
  class numerical_domain_product2: public writeable, numerical_domain< Number, VariableName > {

  public:
    typedef numerical_domain_product2< Number, VariableName, Domain1, Domain2 > numerical_domain_product2_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;

  private:
    typedef domain_product2< Domain1, Domain2 > domain_product2_t;

  private:
    domain_product2_t _product;

  private:
    numerical_domain_product2(domain_product2_t product): _product(product) { }

    void reduce() {
      if (this->_product.first().is_bottom() || this->_product.second().is_bottom()) {
	_product = domain_product2_t::bottom();
      }
    }

  public:
    static numerical_domain_product2_t top() {
      return numerical_domain_product2_t(domain_product2_t::top());
    }
    
    static numerical_domain_product2_t bottom() {
      return numerical_domain_product2_t(domain_product2_t::bottom());
    }
    
  public:
    numerical_domain_product2(): _product() { }
    
    numerical_domain_product2(Domain1 first, Domain2 second): _product(domain_product2_t(first, second)) { }

    numerical_domain_product2(const numerical_domain_product2_t& other): writeable(), numerical_domain< Number, VariableName >(), _product(other._product) { }
    
    numerical_domain_product2_t& operator=(numerical_domain_product2_t other) {
      this->_product = other._product;
      return *this;
    }

    bool is_bottom() {
      return this->_product.is_bottom();
    }

    bool is_top() {
      return this->_product.is_top();
    }

    Domain1& first() {
      return this->_product.first();
    }

    Domain2& second() {
      return this->_product.second();
    }

    bool operator<=(numerical_domain_product2_t other) {
      return (this->_product <= other._product);
    }

    bool operator==(numerical_domain_product2_t other) {
      return (this->_product == other._product);
    }

    numerical_domain_product2_t operator|(numerical_domain_product2_t other) {
      return numerical_domain_product2_t(this->_product | other._product);
    }

    numerical_domain_product2_t operator&(numerical_domain_product2_t other) {
      return numerical_domain_product2_t(this->_product & other._product);
    }

    numerical_domain_product2_t operator||(numerical_domain_product2_t other) {
      return numerical_domain_product2_t(this->_product || other._product);
    }

    numerical_domain_product2_t operator&&(numerical_domain_product2_t other) {
      return numerical_domain_product2_t(this->_product && other._product);
    }

    void assign(VariableName x, linear_expression_t e) {
      this->_product.first().assign(x, e);
      this->_product.second().assign(x, e);
    }

    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
      this->_product.first().apply(op, x, y, z);
      this->_product.second().apply(op, x, y, z);
      this->reduce();
    }

    void apply(operation_t op, VariableName x, VariableName y, Number k) {
      this->_product.first().apply(op, x, y, k);
      this->_product.second().apply(op, x, y, k);
      this->reduce();
    }

    void operator+=(linear_constraint_system_t csts) {
      this->_product.first() += csts;
      this->_product.second() += csts;
      this->reduce();
    }

    void operator-=(VariableName v) {
      this->_product.first() -= v;
      this->_product.second() -= v;
    }
    
    std::ostream& write(std::ostream& o) {
      return this->_product.write(o);
    }

  }; // class numerical_domain_product2

  template< typename Number, typename VariableName, typename Domain1, typename Domain2, typename Domain3 >
  class numerical_domain_product3: public writeable, numerical_domain< Number, VariableName > {

  public:
    typedef numerical_domain_product3< Number, VariableName, Domain1, Domain2, Domain3 > numerical_domain_product3_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;

  private:
    typedef numerical_domain_product2< Number, VariableName, Domain2, Domain3 > product23_t;
    typedef numerical_domain_product2< Number, VariableName, Domain1, product23_t > product123_t;
    
  private:
    product123_t _product;
    
  private:
    numerical_domain_product3(product123_t product): _product(product) { }
    
  public:
    static numerical_domain_product3_t top() {
      return numerical_domain_product3_t(product123_t::top());
    }
    
    static numerical_domain_product3_t bottom() {
      return numerical_domain_product3_t(product123_t::bottom());
    }
    
  public:
    numerical_domain_product3(): _product() { }
    
    numerical_domain_product3(Domain1 first, Domain2 second, Domain3 third): _product(product123_t(first, product23_t(second, third))) { }

    numerical_domain_product3(const numerical_domain_product3_t& other): writeable(), numerical_domain< Number, VariableName >(), _product(other._product) { }
    
    numerical_domain_product3_t& operator=(numerical_domain_product3_t other) {
      this->_product = other._product;
      return *this;
    }
    
    bool is_bottom() {
      return this->_product.is_bottom();
    }
    
    bool is_top() {
      return this->_product.is_top();
    }
    
    Domain1& first() {
      return this->_product.first();
    }

    Domain2& second() {
      return this->_product.second().first();
    }

    Domain3& third() {
      return this->_product.second().second();
    }
    
    bool operator<=(numerical_domain_product3_t other) {
      return (this->_product <= other._product);
    }

    bool operator==(numerical_domain_product3_t other) {
      return (this->_product == other._product);
    }

    numerical_domain_product3_t operator|(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product | other._product);
    }

    numerical_domain_product3_t operator&(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product & other._product);
    }

    numerical_domain_product3_t operator||(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product || other._product);
    }

    numerical_domain_product3_t operator&&(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product && other._product);
    }

    void assign(VariableName x, linear_expression_t e) {
      this->_product.assign(x, e);
    }

    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
      this->_product.apply(op, x, y, z);
    }

    void apply(operation_t op, VariableName x, VariableName y, Number k) {
      this->_product.apply(op, x, y, k);
    }

    void operator+=(linear_constraint_system_t csts) {
      this->_product += csts;
    }

    void operator-=(VariableName v) {
      this->_product -= v;
    }
    
    std::ostream& write(std::ostream& o) {
      if (this->is_bottom()) {
	o << "_|_";
      } else {
	o << "(" << this->first() << ", " << this->second() << ", " << this->third() << ")";
      }
      return o;
    }

  }; // class numerical_domain_product3

} // namespace ikos

#endif // IKOS_DOMAIN_PRODUCTS_HPP
