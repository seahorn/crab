/*******************************************************************************
 *
 * Reduced product of intervals and congruences.
 *
 * The reduce operator based on "Static Analysis of Arithmetical
 * Congruences" by P. Granger published in International Journal of
 * Computer Mathematics, 1989.
 ******************************************************************************/


#ifndef IKOS_INTERVALS_AND_CONGRUENCES_HPP
#define IKOS_INTERVALS_AND_CONGRUENCES_HPP

#include <iostream>
#include <ikos/common.hpp>
#include <ikos/numerical_domains_api.hpp>
#include <ikos/bitwise_operators_api.hpp>
#include <ikos/division_operators_api.hpp>
#include <ikos/domain_products.hpp>
#include <ikos/congruences.hpp>
#include <ikos/intervals.hpp>

namespace ikos {

template< typename Number, unsigned int typeSize = -1 >
class interval_congruence: public writeable {

 public:
  typedef interval_congruence< Number, typeSize > interval_congruence_t;

 private:
  typedef interval< Number >              interval_t;
  typedef bound< Number >                 bound_t;
  typedef congruence< Number, typeSize >  congruence_t;
  
 private:
  interval_t   _first;
  congruence_t _second;

 private:

  interval_congruence(bool is_bottom): 
     writeable(), 
     _first(is_bottom  ? interval_t::bottom()   : interval_t::top()), 
     _second(is_bottom ? congruence_t::bottom() : congruence_t::top()){ }

 public:

  static interval_congruence_t top(){
    return interval_congruence(false);
  }

  static interval_congruence_t bottom(){
    return interval_congruence(true);
  }

 private:
  
  inline Number abs(Number x){
    return x < 0 ? -x: x;
  }
  
  // R(c,a) is the least element of c greater or equal than a
  inline Number R(congruence_t c, Number a){
    Number m = c.get().first;
    Number p = c.get().second;
    return a + ((p-a) % abs(m));
  }
  
  // L(c,a) is the greatest element of c smaller or equal than a
  inline Number L(congruence_t c, Number a){
    Number m = c.get().first; // modulo
    Number p = c.get().second;
    return a - ((p-a) % abs(m));
  }
  
 public:

  interval_congruence(Number n): _first(interval_t(n)), _second(congruence_t(n)) { }

  interval_congruence(interval_t i, congruence_t c): _first(i), _second(c){ 
    this->reduce();
  }

  interval_congruence(interval_t i): _first(i), _second(congruence_t::top()){
    this->reduce();
  }

  interval_congruence(congruence_t c): _first(interval_t::top()), _second(c){
    this->reduce();
  }

  interval_congruence(const interval_congruence &other): 
     _first(other._first), _second(other._second){ }
  
  interval_congruence_t& operator=(interval_congruence_t other){
    this->_first  = other._first;
    this->_second = other._second;
    return *this;
  }

  interval_t&   first(){
    return this->_first;
  }

  congruence_t& second(){
    return this->_second;
  }

  /* 
     Let (i,c) be a pair of interval and congruence
     if (c.is_bottom() || i.is_bottom()) (bottom(), bottom());
     if (c = 0Z+a and a \notin i)        (bottom(), bottom());
     if (c = 0Z+a)                       ([a,a]   , c);
     if (i=[a,b] and R(c,a) > L(c,b))    (bottom(), bottom());
     if (i=[a,b])                        ([R(c,a), L(c,b)], c);
     if (i=[a,+oo])                      ([R(c,a), +oo], c);
     if (i=[-oo,b])                      ([-oo, L(c,b)], c);
     otherwise                           (i,c)
  */
  
  void reduce(){

    interval_t   i = first();
    congruence_t c = second();

    // congruence is top and interval is a singleton
    if (c.is_top()){
      boost::optional< Number> n = i.singleton();
      if (n)
        c = congruence_t(*n);
      return;
    }
    
    Number modulo = c.get().first;
    if (modulo == 0){
      // congruence is a singleton so we refine the interval
      interval_t a(c.get().second);
      if (!(a <= i)){
        i = interval_t::bottom();
        c = congruence_t::bottom();
      }
      else
        i = a;
    }
    else{
      // refine lower and upper bounds of the interval using
      // congruences
      bound_t lb = i.lb();
      bound_t ub = i.ub();

      if ( (lb.is_finite() && ub.is_finite())){
        Number x   = R(c,*(lb.number()));
        Number y   = L(c,*(ub.number()));
        if (x > y){
          i = interval_t::bottom();
          c = congruence_t::bottom();
        }
        else{
          i = interval_t(bound_t(x), bound_t(y));
        }
      }
      else if (lb.is_finite()){
        Number x   = R(c,*(lb.number()));
        i = interval_t(bound_t(x),bound_t("+oo"));
      }
      else if (ub.is_finite()){
        Number y   = L(c,*(ub.number()));
        i = interval_t(bound_t("-oo"), bound_t(y));
      }
      else{ // interval is top        
        }
    }
  }

  std::ostream& write(std::ostream& o) {
    o << "(" << this->_first << ", " << this->_second << ")";
    return o;
  }

 public:

    interval_congruence_t operator+(interval_congruence_t x) {
      interval_congruence_t y(this->_first.operator+(x.first()),
                              this->_second.operator+(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t operator-(interval_congruence_t x) {
      interval_congruence_t y(this->_first.operator-(x.first()),
                              this->_second.operator-(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t operator*(interval_congruence_t x) {
      interval_congruence_t y(this->_first.operator*(x.first()),
                              this->_second.operator*(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t operator/(interval_congruence_t x) {
      interval_congruence_t y(this->_first.operator/(x.first()),
                              this->_second.operator/(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t operator|(interval_congruence_t other){
      return interval_congruence_t(this->_first  | other._first,
                                   this->_second | other._second);
    }

    interval_congruence_t operator&(interval_congruence_t other){
      return interval_congruence_t(this->_first  & other._first,
                                   this->_second & other._second);
    }

 public:

    // division and remainder operations

    interval_congruence_t SDiv(interval_congruence_t x){
      interval_congruence_t y(this->_first.SDiv(x.first()),
                              this->_second.SDiv(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t UDiv(interval_congruence_t x){
      interval_congruence_t y(this->_first.UDiv(x.first()),
                              this->_second.UDiv(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t SRem(interval_congruence_t x){
      interval_congruence_t y(this->_first.SRem(x.first()),
                              this->_second.SRem(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t URem(interval_congruence_t x){
      interval_congruence_t y(this->_first.URem(x.first()),
                              this->_second.URem(x.second()));
      y.reduce();
      return y;
    }
   
    // bitwise operations

    interval_congruence_t Trunc(unsigned width){
      interval_congruence_t y(this->_first.Trunc(width),
                              this->_second.Trunc(width));
      y.reduce();
      return y;
    }

    interval_congruence_t ZExt(unsigned width){
      interval_congruence_t y(this->_first.ZExt(width),
                              this->_second.ZExt(width));
      y.reduce();
      return y;
    }
    
    interval_congruence_t SExt(unsigned width){
      interval_congruence_t y(this->_first.SExt(width),
                              this->_second.SExt(width));
      y.reduce();
      return y;
    }

    interval_congruence_t And(interval_congruence_t x){
      interval_congruence_t y(this->_first.And(x.first()),
                              this->_second.And(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t Or(interval_congruence_t x){
      interval_congruence_t y(this->_first.Or(x.first()),
                              this->_second.Or(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t Xor(interval_congruence_t x){
      interval_congruence_t y(this->_first.Xor(x.first()),
                              this->_second.Xor(x.second()));
      y.reduce();
      return y;
    }

    
    interval_congruence_t Shl(interval_congruence_t x){
      interval_congruence_t y(this->_first.Shl(x.first()),
                              this->_second.Shl(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t LShr(interval_congruence_t x){
      interval_congruence_t y(this->_first.LShr(x.first()),
                              this->_second.LShr(x.second()));
      y.reduce();
      return y;
    }

    interval_congruence_t AShr(interval_congruence_t x){
      interval_congruence_t y(this->_first.AShr(x.first()),
                              this->_second.AShr(x.second()));
      y.reduce();
      return y;
    }

}; 

  template< typename Number, typename VariableName, std::size_t max_reduction_cycles = 3, unsigned int typeSize = -1>
  class interval_congruence_domain: public writeable, 
                                    public numerical_domain< Number, VariableName >,
                                    public bitwise_operators< Number, VariableName >, 
                                    public division_operators< Number, VariableName > {

   public:
    // note that this is assuming that all variables have the same
    // bit width for the congruence domain which is unrealistic.
    typedef interval_congruence< Number, typeSize> interval_congruence_t;
    typedef interval_congruence_domain< Number, VariableName, max_reduction_cycles, typeSize > 
                                                             interval_congruence_domain_t;

   public:
    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName >        linear_expression_t;
    typedef linear_constraint< Number, VariableName >        linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    
   private:
    typedef interval_domain< Number, VariableName> interval_domain_t;
    typedef congruence_domain< Number, VariableName, typeSize> congruence_domain_t;
    typedef numerical_domain_product2< Number, VariableName, 
                                       interval_domain_t, congruence_domain_t > domain_product2_t;
                                       
   private:
    domain_product2_t _product;

   private:
    interval_congruence_domain(domain_product2_t product): _product(product) { }


    void reduce() {

      interval_domain_t   I = this->_product.first();
      congruence_domain_t C = this->_product.second();

      if (I.is_bottom() || C.is_bottom()) {
          _product = domain_product2_t::bottom();
          return;
      }

      // In general, the reduce operation may need a fixpoint.
      bool change=true;
      unsigned int num_iters = 0;
      while (change && (num_iters < max_reduction_cycles)){
        change=false;
        num_iters++;

        for(typename interval_domain_t::iterator it = I.begin(); it != I.end(); ++it){
          VariableName var_name = (*it).first;
          interval_congruence_t old_val((*it).second, C[var_name]);
          // An interval can only improve a congruence if the
          // congruence is top and the interval is a singleton.
          if (!old_val.first().is_top() && old_val.second().is_top()){
            interval_congruence_t new_val(old_val);
            new_val.reduce();
            change = (new_val.second() != old_val.second());
            if (change)
              C[var_name] = new_val.second();
          }
        } // end for

        for(typename congruence_domain_t::iterator it = C.begin(); it != C.end(); ++it){
          VariableName var_name = (*it).first;
          interval_congruence_t old_val(I[var_name], (*it).second);
          interval_congruence_t new_val(old_val);
          new_val.reduce();
          change = ( new_val.first()!=old_val.first() || new_val.second()!=old_val.second());
          if (change){
            I[var_name] = new_val.first();
            C[var_name] = new_val.second();
          }
        } // end for
      } // end while
    }

   public:
    static interval_congruence_domain_t top() {
      return interval_congruence_domain_t(domain_product2_t::top());
    }
    
    static interval_congruence_domain_t bottom() {
      return interval_congruence_domain_t(domain_product2_t::bottom());
    }
    
   public:

    interval_congruence_domain(): _product() { }
    
    interval_congruence_domain(const interval_congruence_domain_t& other): 
       writeable(), 
       numerical_domain< Number, VariableName >(), 
       bitwise_operators< Number, VariableName >(),
       division_operators< Number, VariableName > (),
       _product(other._product) { }
    
    interval_congruence_domain_t& operator=(interval_congruence_domain_t other) {
      this->_product = other._product;
      return *this;
    }

    bool is_bottom() {
      return this->_product.is_bottom();
    }

    bool is_top() {
      return this->_product.is_top();
    }

    interval_domain_t& first() {
      return this->_product.first();
    }

    congruence_domain_t& second() {
      return this->_product.second();
    }

    bool operator<=(interval_congruence_domain_t other) {
      return (this->_product <= other._product);
    }

    bool operator==(interval_congruence_domain_t other) {
      return (this->_product == other._product);
    }

    interval_congruence_domain_t operator|(interval_congruence_domain_t other) {
      return interval_congruence_domain_t(this->_product | other._product);
    }

    interval_congruence_domain_t operator&(interval_congruence_domain_t other) {
      return interval_congruence_domain_t(this->_product & other._product);
    }

    interval_congruence_domain_t operator||(interval_congruence_domain_t other) {
      return interval_congruence_domain_t(this->_product || other._product);
    }

    interval_congruence_domain_t operator&&(interval_congruence_domain_t other) {
      return interval_congruence_domain_t(this->_product && other._product);
    }

    // pre: x is already reduced
    void set(VariableName v, interval_congruence_t x) {
      this->_product.first().set(v, x.first());
      this->_product.second().set(v, x.second());
    }

    interval_congruence_t operator[](VariableName v) {
      return interval_congruence_t(this->_product.first()[v], this->_product.second()[v]);
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

    // bitwise_operators_api
    
    void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width){
      this->_product.first().apply(op, x, y, width);
      this->_product.second().apply(op, x, y, width);
      this->reduce();
    }

    void apply(conv_operation_t op, VariableName x, Number k, unsigned width){
      this->_product.first().apply(op, x, k, width);
      this->_product.second().apply(op, x, k, width);
      this->reduce();
    }

    void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z){
      this->_product.first().apply(op, x, y, z);
      this->_product.second().apply(op, x, y, z);
      this->reduce();
    }
    
    void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k){
      this->_product.first().apply(op, x, y, k);
      this->_product.second().apply(op, x, y, k);
      this->reduce();
    }
    
    // division_operators_api
    
    void apply(div_operation_t op, VariableName x, VariableName y, VariableName z){
      this->_product.first().apply(op, x, y, z);
      this->_product.second().apply(op, x, y, z);
      this->reduce();
    }

    void apply(div_operation_t op, VariableName x, VariableName y, Number k){
      this->_product.first().apply(op, x, y, k);
      this->_product.second().apply(op, x, y, k);
      this->reduce();

    }
    
    std::ostream& write(std::ostream& o) {
      return this->_product.write(o);
    }

  const char* getDomainName () const {return "Intervals+Congruences";}

  }; // class interval_congruence_domain

} // namespace ikos

#endif // IKOS_INTERVALS_AND_CONGRUENCES_HPP
