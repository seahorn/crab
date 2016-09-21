/*******************************************************************************
 *
 * Products of abstract domains.
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

#ifndef IKOS_DOMAIN_PRODUCTS_HPP
#define IKOS_DOMAIN_PRODUCTS_HPP

#include <crab/common/types.hpp>
#include <crab/domains/operators_api.hpp>

namespace ikos {

  template< typename Domain1, typename Domain2 >
  class domain_product2: public writeable {
    
  public:
    typedef domain_product2< Domain1, Domain2 > domain_product2_t;
    typedef Domain1 first_type;
    typedef Domain2 second_type;

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
        this->_is_bottom = this->_first.is_bottom() || 
            this->_second.is_bottom();
        if (this->_is_bottom) {
          this->_first = Domain1::bottom();
          this->_second = Domain2::bottom();      
        }
      }
    }

  public:
    domain_product2(): 
        _is_bottom(false), 
        _first(Domain1::top()), _second(Domain2::top()) { }
    
    domain_product2(Domain1 first, Domain2 second): 
        _is_bottom(false), 
        _first(first), _second(second) {
      this->canonicalize();
    }
    
    domain_product2(const domain_product2_t& other): 
        writeable(), 
        _is_bottom(other._is_bottom), 
        _first(other._first), _second(other._second) { }
    
    domain_product2_t& operator=(const domain_product2_t& other) {
      if (this != &other) {
        this->_is_bottom = other._is_bottom;
        this->_first = other._first;
        this->_second = other._second;
      }
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
        return (this->_first <= other._first) && 
            (this->_second <= other._second);
      }
    }
    
    bool operator==(domain_product2_t other) {
      return (this->operator<=(other) && other.operator<=(*this));
    }

    void operator|=(domain_product2_t other) {
      if (this->is_bottom()) {
        *this = other;
      } else if (other.is_bottom()) {
        return;
      } else {
        this->_first |= other._first;
        this->_second |= other._second;
      }
    }
    
    domain_product2_t operator|(domain_product2_t other) {
      if (this->is_bottom()) {
        return other;
      } else if (other.is_bottom()) {
        return *this;
      } else {
        return domain_product2_t(this->_first | other._first, 
                                 this->_second | other._second);
      }
    }

    domain_product2_t operator||(domain_product2_t other) {
      if (this->is_bottom()) {
        return other;
      } else if (other.is_bottom()) {
        return *this;
      } else {
        return domain_product2_t(this->_first || other._first, 
                                 this->_second || other._second);
      }
    }

    domain_product2_t operator&(domain_product2_t other) {
      if (this->is_bottom() || other.is_bottom()) {
        return bottom();
      } else {
        return domain_product2_t(this->_first & other._first, 
                                 this->_second & other._second);
      }
    }

    domain_product2_t operator&&(domain_product2_t other) {
      if (this->is_bottom() || other.is_bottom()) {
        return bottom();
      } else {
        return domain_product2_t(this->_first && other._first, 
                                 this->_second && other._second);
      }
    }
    
    void write(crab::crab_os& o) {
      if (this->is_bottom()) {
        o << "_|_";
      } else {
        o << "(" << this->_first << ", " << this->_second << ")";
      }
    }

    static std::string getDomainName () { 
      std::string name = "Product(" +
          Domain1::getDomainName () + "," +
          Domain2::getDomainName () + ")";
      return name;
    }

  }; // class domain_product2

  template< typename Domain1, typename Domain2, typename Domain3 >
  class domain_product3: public writeable {

  public:
    typedef domain_product3< Domain1, Domain2, Domain3 > domain_product3_t;
    typedef Domain1 first_type;
    typedef Domain2 second_type;
    typedef Domain3 third_type;

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

    domain_product3(Domain1 first, Domain2 second, Domain3 third): 
        _product(first, product23_t(second, third)) { }

    domain_product3(const domain_product3_t& other): 
        writeable(), 
        _product(other._product) { }
    
    domain_product3_t& operator=(const domain_product3_t& other) {
      if (this != &other)
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

    void operator|=(domain_product3_t other) {
      *this = *this | other;
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

    void write(crab::crab_os& o) {
      if (this->is_bottom()) {
        o << "_|_";
      } else {
        o << "(" 
          << this->first() << ", " 
          << this->second() << ", " 
          << this->third() << ")";
      }
    }
    
  }; // class domain_product3

  // XXX: this domain product is not just numerical operations. It
  // also includes array and pointer operations.
  template< typename Number, typename VariableName, typename Domain1, typename Domain2 >
  class numerical_domain_product2: public writeable,
                                   public numerical_domain< Number, VariableName >,
                                   public bitwise_operators< Number, VariableName >, 
                                   public division_operators< Number, VariableName >,
                                   public crab::domains::array_operators< Number, VariableName >,
                                   public crab::domains::pointer_operators< Number, VariableName > {
  public:
    typedef numerical_domain_product2< Number, VariableName, Domain1, Domain2 > numerical_domain_product2_t;
    typedef Domain1 first_type;
    typedef Domain2 second_type;
    
    using typename numerical_domain< Number, VariableName >::linear_expression_t;
    using typename numerical_domain< Number, VariableName >::linear_constraint_t;
    using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
    using typename numerical_domain< Number, VariableName >::variable_t;
    using typename numerical_domain< Number, VariableName >::number_t;
    using typename numerical_domain< Number, VariableName >::varname_t;

    typedef crab::pointer_constraint<VariableName> ptr_cst_t;

  private:
    typedef domain_product2< Domain1, Domain2 > domain_product2_t;

  private:
    domain_product2_t _product;

  private:
    numerical_domain_product2(domain_product2_t product): _product(product) { }

    void reduce() {
      if (this->_product.first().is_bottom() || 
          this->_product.second().is_bottom()) {
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
    
    numerical_domain_product2(Domain1 first, Domain2 second): 
        _product(domain_product2_t(first, second)) { }

    numerical_domain_product2(const numerical_domain_product2_t& other): 
        writeable(), 
        numerical_domain< Number, VariableName >(), 
        _product(other._product) { }
    
    numerical_domain_product2_t& operator=(const numerical_domain_product2_t& other) {
      if (this != &other)
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

    void operator|=(numerical_domain_product2_t other) {
      this->_product |= other._product;
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

    template<typename Thresholds>
    numerical_domain_product2_t widening_thresholds (numerical_domain_product2_t other,
                                                     const Thresholds& ts) {
      return numerical_domain_product2_t (
          this->_product.first ().widening_thresholds (other._product.first (), ts),
          this->_product.second ().widening_thresholds (other._product.second (), ts));
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

    // bitwise_operators_api

    void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
      this->_product.first().apply(op, x, y, width);
      this->_product.second().apply(op, x, y, width);
      this->reduce();
    }
    
    void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
      this->_product.first().apply(op, x, k, width);
      this->_product.second().apply(op, x, k, width);
      this->reduce();
    }
    
    void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
      this->_product.first().apply(op, x, y, z);
      this->_product.second().apply(op, x, y, z);
      this->reduce();
    }
    
    void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
      this->_product.first().apply(op, x, y, k);
      this->_product.second().apply(op, x, y, k);
      this->reduce();
    }
    
    // division_operators_api
    
    void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
      this->_product.first().apply(op, x, y, z);
      this->_product.second().apply(op, x, y, z);
      this->reduce();
    }
    
    void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
      this->_product.first().apply(op, x, y, k);
      this->_product.second().apply(op, x, y, k);
      this->reduce();
    }

    // array_operators_api
    virtual void array_init (VariableName a, 
                             const vector<ikos::z_number>& values) override {
      this->_product.first().array_init (a, values);
      this->_product.second().array_init (a, values);
      this->reduce ();
    }

    virtual void array_assume (VariableName a,
                               boost::optional<Number> lb, boost::optional<Number> ub) override {
      this->_product.first().array_assume (a, lb, ub);
      this->_product.second().array_assume (a, lb, ub);
      this->reduce ();
    }

    virtual void array_load (VariableName lhs, VariableName a, 
                             VariableName i, ikos::z_number bytes) override {
      this->_product.first().array_load (lhs, a, i, bytes);
      this->_product.second().array_load (lhs, a, i, bytes);
      this->reduce ();
    }

    virtual void array_store (VariableName a, VariableName i,
                              linear_expression_t val, ikos::z_number bytes,
                              bool is_singleton) override {
      this->_product.first().array_store (a, i, val, bytes, is_singleton);
      this->_product.second().array_store (a, i, val, bytes, is_singleton);
      this->reduce ();
    }

    // pointer_operators_api
    virtual void pointer_load (VariableName lhs, VariableName rhs) override {
      this->_product.first().pointer_load (lhs, rhs);
      this->_product.second().pointer_load (lhs, rhs);
      this->reduce ();
    }

    virtual void pointer_store (VariableName lhs, VariableName rhs) override {
      this->_product.first().pointer_store (lhs, rhs);
      this->_product.second().pointer_store (lhs, rhs);
      this->reduce ();
    } 

    virtual void pointer_assign (VariableName lhs, VariableName rhs, linear_expression_t offset) override {
      this->_product.first().pointer_assign (lhs, rhs, offset);
      this->_product.second().pointer_assign (lhs, rhs, offset);
      this->reduce ();
    }

    virtual void pointer_mk_obj (VariableName lhs, ikos::index_t address) override {
      this->_product.first().pointer_mk_obj (lhs, address);
      this->_product.second().pointer_mk_obj (lhs, address);
      this->reduce ();
    }

    virtual void pointer_function (VariableName lhs, VariableName func) override {
      this->_product.first().pointer_function (lhs, func);
      this->_product.second().pointer_function (lhs, func);
      this->reduce ();
    }
    
    virtual void pointer_mk_null (VariableName lhs) override {
      this->_product.first().pointer_mk_null (lhs);
      this->_product.second().pointer_mk_null (lhs);
      this->reduce ();
    }
    
    virtual void pointer_assume (ptr_cst_t cst) override {
      this->_product.first().pointer_assume (cst);
      this->_product.second().pointer_assume (cst);
      this->reduce ();
    }    

    virtual void pointer_assert (ptr_cst_t cst) override {
      this->_product.first().pointer_assert (cst);
      this->_product.second().pointer_assert (cst);
      this->reduce ();
    }    

    void write(crab::crab_os& o) {
      this->_product.write(o);
    }

    static std::string getDomainName () { 
      return domain_product2_t::getDomainName ();
    }

  }; // class numerical_domain_product2

  template< typename Number, typename VariableName, 
            typename Domain1, typename Domain2, typename Domain3 >
  class numerical_domain_product3: 
      public writeable, 
      public numerical_domain< Number, VariableName > {

  public:
    typedef numerical_domain_product3< Number, VariableName, 
                                       Domain1, Domain2, Domain3 > 
    numerical_domain_product3_t;
    typedef Domain1 first_type;
    typedef Domain2 second_type;
    typedef Domain3 third_type;

    using typename numerical_domain< Number, VariableName >::linear_expression_t;
    using typename numerical_domain< Number, VariableName >::linear_constraint_t;
    using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
    using typename numerical_domain< Number, VariableName >::variable_t;
    using typename numerical_domain< Number, VariableName >::number_t;
    using typename numerical_domain< Number, VariableName >::varname_t;
  private:
    typedef numerical_domain_product2< Number, VariableName, Domain2, Domain3 > product23_t;
    typedef numerical_domain_product2< Number, VariableName, Domain1, product23_t > product123_t;
    
  private:
    product123_t _product;
    
  private:
    numerical_domain_product3(product123_t product): 
        _product(product) { }
    
  public:
    static numerical_domain_product3_t top() {
      return numerical_domain_product3_t(product123_t::top());
    }
    
    static numerical_domain_product3_t bottom() {
      return numerical_domain_product3_t(product123_t::bottom());
    }
    
  public:
    numerical_domain_product3(): _product() { }
    
    numerical_domain_product3(Domain1 first, Domain2 second, Domain3 third): 
        _product(product123_t(first, product23_t(second, third))) { }

    numerical_domain_product3(const numerical_domain_product3_t& other): 
        writeable(), 
        numerical_domain< Number, VariableName >(), 
        _product(other._product) { }
    
    numerical_domain_product3_t& 
    operator=(const numerical_domain_product3_t& other) {
      if (this != &other)
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

    void operator|=(numerical_domain_product3_t other) {
      *this = *this | other;
    }
    
    numerical_domain_product3_t 
    operator|(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product | other._product);
    }

    numerical_domain_product3_t 
    operator&(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product & other._product);
    }

    numerical_domain_product3_t 
    operator||(numerical_domain_product3_t other) {
      return numerical_domain_product3_t(this->_product || other._product);
    }

    template<typename Thresholds>
    numerical_domain_product3_t widening_thresholds (numerical_domain_product3_t other,
                                                     const Thresholds& ts) {
      return numerical_domain_product3_t (
          this->_product.widening_thresholds (other._product, ts));
    }

    numerical_domain_product3_t 
    operator&&(numerical_domain_product3_t other) {
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
    
    void write(crab::crab_os& o) {
      if (this->is_bottom()) {
        o << "_|_";
      } else {
        o << "(" << this->first() << ", " << this->second() << ", " << this->third() << ")";
      }
    }
    
  }; // class numerical_domain_product3

} // namespace ikos
#endif // IKOS_DOMAIN_PRODUCTS_HPP
