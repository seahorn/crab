/*******************************************************************************
 *
 * Standard domain of numerical congruences extended with bitwise
 * operations. 
 *
 * Bitwise extensions based on Stefan Bygde's paper: Static WCET
 * analysis based on abstract interpretation and counting of elements,
 * Vasteras : School of Innovation, Design and Engineering, Malardalen
 * University (2010)
 *
 * Comments:
 * - Assume that no overflow/underflow can occur.
 * - Bitwise operations are sound approximations for both signed and
 *   unsigned interpretations of bit strings.
 * - If typeSize > -1 then the domain assumes that all variables will
 *   have the same bit width which is unrealistic.
 * 
 ******************************************************************************/

#ifndef IKOS_CONGRUENCES_HPP
#define IKOS_CONGRUENCES_HPP

#include <iostream>
#include <ikos_domains/common.hpp>
#include <ikos_domains/separate_domains.hpp>
#include <ikos_domains/numerical_domains_api.hpp>
#include <ikos_domains/bitwise_operators_api.hpp>
#include <ikos_domains/division_operators_api.hpp>

using namespace std;

namespace ikos{
  
  template < typename Number , unsigned int typeSize = -1 >
  class congruence : public writeable {
  public:
    typedef congruence<Number, typeSize > congruence_t;
    static const unsigned int infty = -1;

  private:
    bool   _is_bottom;
    Number a, b;	// Of the form (aZ + b)
    
  private:

    congruence(bool b): writeable(), _is_bottom(!b), a(1), b(0) { }

    congruence(int x): writeable(), _is_bottom(false), a(0), b(x) { }

    congruence(Number a, Number b): writeable(), _is_bottom(false), 
        a(a), 
        // Set to standard form: 0 <= b < a for a != 0.
        b((a == 0) ? b : ((b % a) < 0 ? (b % a) + a : b % a)) { }	
    
    Number abs(Number x){
      return x < 0 ? -x: x;
    }

    static Number max(Number x, Number y) {
      return (x.operator<=(y) ? y : x);
    }

    static Number min(Number x, Number y) {
      return (x.operator<(y) ? x : y);
    }

    Number gcd(Number x, Number y, Number z){
      return gcd(x, gcd(y, z));
    }

    Number gcd(Number x, Number y){
      return gcd_helper(abs(x), abs(y));
    }

    // Not to be called explicitly outside of _gcd
    Number gcd_helper(Number x, Number y){	
      return (y == 0) ? x : gcd_helper(y, x % y);
    }

    Number lcm(Number x, Number y){
      Number tmp = gcd(x, y);
      if (tmp == 0) throw error("lcm causes division by zero");
      return abs(x* y)/tmp;
    }
    
    // Weakening operator
    congruence_t weakened(){	
      // Returns weakened version of _this, does not modify
      // _this. Returned _a is guaranteed to be a power of 2.
      if(a  == 0){
        return *this;
      }

      if (typeSize < 32){
        Number tmp = gcd(a, Number(1) << Number(typeSize)); 
        if (tmp == 0) throw error("weakened causes division by zero");
        return congruence_t(tmp, b % tmp);
      }
      else{
        throw error("if int is 32 bits shift can only range from 0 to 31");
      }
    }
    
    Number scanr1(Number x){	
      // Searches for the right-most bit containing a value of 1.
      Number count(1);
      for(Number idx(1); idx < typeSize+ 1; ){
        if((idx & x) > 0){
          return count;
        }
        ++count;
        idx = idx << 1;
      }
      return Number(-1); // to indicate not found, x==0
    }

    Number scanl1(Number x){	
      // Searches for the left-most bit containing a value of 1.
      if (typeSize < 32){
        Number count(sizeof(x));
        Number idx = Number(1) << Number((typeSize - 1));
        while (idx > 0){
          if( (idx & x) > 0){
            return count;
          }
          --count;
          idx = idx >> 1;
        }
        return Number(-1); // to indicate not found, x==0
      }
      else{
        throw error("if int is 32 bits shift can only range from 0 to 31");
      }
    }    
  public:
    
    static congruence_t top(){
      return congruence(true);
    }
    
    static congruence_t bottom(){
      return congruence(false);
    }
    
  public:

    congruence(): writeable(), _is_bottom(false), a(1), b(0) { }

    congruence(Number b): writeable(), _is_bottom(false), a(0), b(b) { }

    congruence(const congruence_t &o): writeable(), _is_bottom(o._is_bottom), a(o.a), b(o.b) { }
    
    congruence_t operator=(congruence_t o){
      _is_bottom= o._is_bottom;
      a= o.a;
      b= o.b;
      return *this;
    }
    
    bool is_bottom(){
      return _is_bottom;
    }
    bool is_top(){
      return (a ==1);
    }

    std::pair<Number,Number> get(){
      // aZ+b
      return std::make_pair(a,b);
    }

    boost::optional< Number > singleton() {
      if (!this->is_bottom() && a == 0) {
	return boost::optional< Number >(b);
      } else {
	return boost::optional< Number >();
      }
    }
    
    bool operator==(congruence_t o){
      return (a == o.a && b == o.b && is_bottom() == o.is_bottom());
    }
    
    bool operator!=(congruence_t x) {
      return !this->operator==(x);
    }

    //Lattice Operations

    bool operator<=(congruence_t o){	
      if(is_bottom()){
        return true;
      } 
      else if(o.is_bottom()){
        return false;
      } 
      else if(a == 0 && o.a == 0 && b == o.b){
        return true;
      } 
      else if(a == 0){ 
        if (o.a == 0) throw error("operator<= causes division by zero");
        if ((b % o.a < 0 ? (b % o.a)+ o.a : b % o.a) == o.b % o.a) 
          return true;
      } 
      else if(o.a == 0){ 
        if (a == 0) throw error("operator<= causes division by zero");
        if (b % a == (o.b % a < 0 ? (o.b % a)+ a : o.b % a))
          return false;
      }
      if (o.a == 0) throw error("operator<= causes division by zero");
      return (a % o.a == 0) && (b % o.a  == o.b % o.a);
    }

    congruence_t operator|(congruence_t o){	
      if(is_bottom()){
        return o;
      } else if(o.is_bottom()){
        return *this;
      }
      return congruence_t(gcd(a, o.a, abs(b - o.b)), min(b,o.b));
    }

    congruence_t operator&(congruence_t o){	

      if(is_bottom() || o.is_bottom()){
        return bottom();
      } 

      // note: lcm has meaning only if both a and o.a are not 0
      if(a == 0 && o.a == 0){
        if (b == o.b)
          return *this;
        else 
          return bottom();
      }
      else if (a == 0){
        // b & a'Z + b' iff \exists k such that a'*k + b' = b iff ((b - b') %a' == 0)
        if ((b - o.b) % o.a == 0)
          return *this;
        else 
          return bottom();
      }
      else if (o.a == 0){
        // aZ+b & b' iff \exists k such that a*k+b  = b' iff ((b'-b %a) == 0)
        if ((o.b - b) % a == 0)
          return o;
        else 
          return bottom();
      }
      else{
        // pre: a and and o.a != 0
        Number x = gcd(a, o.a);
        if (x == 0) throw error("operator& causes division by zero");
        if(b % x == (o.b % x)){	
          // the part max(b,o.b) needs to be verified. What we really
          // want is to find b'' such that 
          // 1) b'' % lcm(a,a') == b  % lcm(a,a'), and 
          // 2) b'' % lcm(a,a') == b' % lcm(a,a').
          // An algorithm for that is provided in Granger'89.
          return congruence_t(lcm(a, o.a), max(b, o.b));
        }
        else 
          return congruence_t::bottom();
      }

      ///
      // The below code is harder to verify. E.g. 2Z + 0 & 1 is
      // bottom, however, by the 4th rule it will return 1.
      ///
      // if(a == 0 && o.a == 0 && b == o.b){
      //   return *this;
      // } else if(a == 0 && o.a == 0){
      //   return bottom();
      // } else if(a == 0 && ((b % o.a < 0 ? (b % o.a)+ o.a : b % o.a) == o.b % o.a)){
      //   return *this;
      // } else if(o.a == 0 && ((o.b % a < 0 ? (o.b % a)+ a : o.b % a) == o.b % a)){
      //   return o;
      // } else if(b % gcd(a, o.a) == (o.b % gcd(a, o.a))){	
      //   // This case needs to be verified.
      //   return congruence_t(lcm(a, o.a), max(b, o.b));
      // }
      // return bottom();
    }
    congruence_t operator||(congruence_t o){ // Equivalent to join, domain is flat
      return *this | o;
    }
    congruence_t operator&&(congruence_t o){ // Simply refines top element
      return (is_top()) ? o : *this;
    }
    
    //Arithmetic Operators

    congruence_t operator+(congruence_t o){	
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else 
        return congruence_t(gcd(a,o.a), b+o.b);
    }
    
    congruence_t operator-(congruence_t o){
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else 
        return congruence_t(gcd(a,o.a), b-o.b);
    }

    congruence_t operator-(){
      if (this->is_bottom() || this->is_top())
        return *this;
      else 
        return congruence_t(a, -b+a);
    }
    
    congruence_t operator*(congruence_t o){
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else 
        return congruence_t(gcd(a*o.a, a*o.b, o.a*b), b*o.b);
    }

    // signed division
    congruence_t operator/(congruence_t o){	

      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (o == congruence(0))
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else{
        /* 
           (Special case 1) aZ+b / 0Z+b':         
           if b'|a then  (a/b')Z + b/b'
           else          top
           
           (Special case 2) 0Z+b / a'Z+b':
           if N>0   (b div N)Z + 0
           else     0Z + 0
           
           where N = a'((b-b') div a') + b'
        */

        if (o.a == 0){
          if (o.b == 0) throw error("operator/ causes division by zero");
          if (a % o.b == 0)
            return congruence_t(a/o.b, b/o.b);
          else 
            return congruence_t::top();
        }
        
        if (a == 0){
          Number n(o.a *(((b - o.b)/o.a) + o.b));
          if (n > 0)
          return congruence_t(b/n, 0);
          else
            return congruence_t(0, 0);          
        }
        
        /*
          General case: no singleton
        */
        return congruence_t::top();	
      }
    }      

    // signed remainder operator
    congruence_t operator%(congruence_t o){
      
      if (this->is_bottom() || o.is_bottom()) 
        return congruence_t::bottom();
      else if (o == congruence(0))
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else{
        /* 
         (Special case 1) aZ+b / 0Z+b':         
           if b'|a then  0Z + (b mod b')
           else          gcd(a,b')Z + b 

         (Special case 2) 0Z+b / a'Z+b':
           if N<=0           then 0Z+b
           if (b div N) == 1 then gcd(b',a')Z + b
           if (b div N) >= 2 then N(b div N)Z  + b

           where N = a'((b-b') div a') + b'
        */
        if (o.a == 0){
          if (o.b == 0) throw error("operator% causes division by zero");
          if (a % o.b == 0){
            return congruence_t(0, b % o.b);
          }
          else 
            return congruence_t(gcd(a, o.b), b);
        }
        if (a == 0){
          Number n(o.a *(((b - o.b)/o.a) + o.b));
          if (n <= 0){
            return congruence_t(a, b);
          }
          else if ( b == n){
            return congruence_t(gcd(o.b,o.a), b);
          }
          else if ( (b/n) >= 2){
            return congruence_t( b , b);
          }
          else 
            error("unreachable");
        }
        /* 
            general case: no singleton
        */
        return congruence_t(gcd(a, o.a, o.b), b);
      }
    }

    //Bitwise operators

   private:
    congruence_t Not(){	
      if (this->is_bottom())
        return congruence_t::bottom();
      else if (this->is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        if (typeSize < 32){
          Number two_to_n = Number(1) << Number(typeSize);
          return congruence_t(gcd(a, two_to_n), -b-1);
        }
        else {
          throw error("if int is 32 bits shift can only range from 0 to 31");
        }
      }
      // For simplicity, we don't implement the case when a==0:
      //   return (is_bottom()) ? bottom() :
      //   (a== 0) ? congruence_t(~b) : congruence_t(gcd(a, 1 << typeSize), -b- 1);
    }

   public:

    congruence_t And(congruence_t o){

      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        // first both numbers are approximated to power of 2
        congruence_t x(this->weakened());
        congruence_t y(o.weakened());
        
        if (x.a == 0 && y.a == 0)
          return congruence_t(0, x.b & y.b);
        
        if (x.a == 0) // (0Z + b) & (a'Z+b')
          return congruence_t((x.b < y.a) ? 0 : y.a, x.b & y.b);
        
        if (y.a == 0) // (aZ + b) & (0Z+b')
          return congruence_t((y.b < x.a) ? 0 : x.a, x.b & y.b);
        
        if (scanr1(y.b) > scanr1(y.a)) // choose the smallest power of 2
          return congruence_t(scanr1(y.a), x.b & y.b);
        else 
          return congruence_t(scanr1(x.a), x.b & y.b);        
      }
    }
    
    congruence_t Or(congruence_t o){	
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        return this->Not().And(o.Not()).Not();
      }
    }
    
    congruence_t Xor(congruence_t o){	
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        return this->And(o.Not()).Or(this->Not().And(o));
      }
    }

    congruence_t Shl(congruence_t o){
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        if (o.a == 0){ // singleton
          // aZ + b << 0Z + b'  = (a*2^b')Z + b*2^b'
          Number x(Number(2) << o.b);
          return congruence_t(a * x, b * x);
        }
        else{
          Number x(Number(2) << o.b);
          Number y(Number(2) << o.a);
          // aZ + b << a'Z + b' = (gcd(a, b * (2^a' - 1)))*(2^b')Z + b*(2^b')
          return congruence_t(gcd(a, b*(y - 1)) * x , b * x);
        }
      }
    }


    congruence_t AShr(congruence_t o){	
      if (this->is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else if (this->is_top() || o.is_top())
        return congruence_t::top();
      else if(typeSize == infty)
        return congruence_t::top();
      else{
        Number t= scanr1(a);
        if (o.a == 0 && scanr1(a) == scanl1(a) && o.b < scanr1(a))
          return congruence_t(Number(2) << (t - o.b), b >> o.b);
        else
          return congruence_t::top();
      }
    }

    congruence_t LShr(congruence_t  o){
      if (is_bottom() || o.is_bottom())
        return congruence_t::bottom();
      else
        return congruence_t::top();
    }


    // division and remainder operations

    congruence_t SDiv(congruence_t x){
      return this->operator/(x);
    }

    congruence_t UDiv(congruence_t x){
      return congruence_t::top();
    }

    congruence_t SRem(congruence_t x){
      return this->operator%(x);
    }

    congruence_t URem(congruence_t x){
      return congruence_t::top();
    }
   
    // cast operations

    congruence_t Trunc(unsigned width){
      if (this->is_bottom())
        return congruence_t::bottom();
      else if (this->is_top())
        return congruence_t::top();
      else
        return congruence_t( gcd(1 << width, a), b);
    }
    
    congruence_t ZExt(unsigned width){
      return *this;
    }
    
    congruence_t SExt(unsigned width){
      return *this;
    }

    
    std::ostream & write(std::ostream &o){

      if(is_bottom()){
        o << "_|_";
        return o;
      }

      if(a== 0){
        o << b ;
        return o;
      }
      
      bool positive = true;
      if(b < 0){
        positive = false;
        b*= -1;
      }
      o << a << "Z" << ((positive) ? "+ " : "- ") << b ;
      return o;
    }
  }; // end class congruence

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator+(Number c, congruence< Number, typeSize > x) {
    return congruence< Number, typeSize >(c) + x;
  }
  
  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator+(congruence< Number, typeSize > x, Number c) {
    return x + congruence< Number, typeSize >(c);
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator*(Number c, congruence< Number, typeSize > x) {
    return congruence< Number, typeSize >(c) * x;
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator*(congruence< Number, typeSize > x, Number c) {
    return x * congruence< Number, typeSize >(c);
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator/(Number c, congruence< Number, typeSize > x) {
    return congruence< Number, typeSize >(c) / x;
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator/(congruence< Number, typeSize > x, Number c) {
    return x / congruence< Number, typeSize >(c);
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator-(Number c, congruence< Number, typeSize > x) {
    return congruence< Number, typeSize >(c) - x;
  }

  template< typename Number, unsigned int typeSize>
  inline congruence< Number, typeSize > operator-(congruence< Number, typeSize > x, Number c) {
    return x - congruence< Number, typeSize >(c);
  }


  template< typename Number, typename VariableName, typename CongruenceCollection, unsigned int typeSize = -1 >
  class equality_congruence_solver {
    // TODO: check correctness of the solver. Granger provides a sound
    // and more precise solver for equality linear congruences (see
    // Theorem 4.4).
  private:
    typedef congruence< Number, typeSize > congruence_t;
    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    
  private:
    typedef std::vector< linear_constraint_t > cst_table_t;
    typedef typename linear_constraint_t::variable_set_t variable_set_t;

  private:
    class bottom_found { };

  private:
    std::size_t _max_cycles;
    bool _is_contradiction;
    cst_table_t _cst_table;
    variable_set_t _refined_variables;
    std::size_t _op_count;
    
  private:
    void refine(variable_t v, congruence_t i, CongruenceCollection& env) {
      VariableName n = v.name();
      congruence_t old_i = env[n];
      congruence_t new_i = old_i & i;
      if (new_i.is_bottom()) {
        throw bottom_found();
      }
      if (old_i != new_i) {
	env.set(n, new_i);
	this->_refined_variables += v;
	++(this->_op_count);
      }
    }

    congruence_t compute_residual(linear_constraint_t cst, variable_t pivot, CongruenceCollection& env) {
      congruence_t residual(cst.constant());
      for (typename linear_constraint_t::iterator it = cst.begin(); it != cst.end(); ++it) {
	variable_t v = it->second;
	if (v.index() != pivot.index()) {
	  residual = residual - (it->first * env[v.name()]);
	  ++(this->_op_count);
	}
      }
      return residual;
    }
    
    void propagate(linear_constraint_t cst, CongruenceCollection& env) {
      for (typename linear_constraint_t::iterator it = cst.begin(); it != cst.end(); ++it) {
	Number c         = it->first;
	variable_t pivot = it->second;
	congruence_t rhs = this->compute_residual(cst, pivot, env) / congruence_t(c);
	if (cst.is_equality()) {
	  this->refine(pivot, rhs, env);
	} 
          else if (cst.is_inequality()) {
            // Inequations (>=, <=, >, and <) do not work well with
            // congruences because for any number n there is always x and y
            // \in gamma(aZ+b) such that n < x and n > y.
            // 
            // The only cases we can catch is when all the expressions
            // are constants. We do not bother because any product
            // with intervals or constants should get those cases.
            continue; 
	} 
          else {
	  // TODO: cst is a disequation
	}
      }
    }
    
    void solve_system(CongruenceCollection& env) {
      std::size_t cycle = 0;
      do {
        ++cycle;
        this->_refined_variables.clear();
        for (typename cst_table_t::iterator it = this->_cst_table.begin(); it != this->_cst_table.end(); ++it) {
          this->propagate(*it, env);
        }
      }
      while (!this->_refined_variables.is_empty() && cycle <= this->_max_cycles);
    }
    
    
  public:

    equality_congruence_solver(linear_constraint_system_t csts, std::size_t max_cycles): 
        _max_cycles(max_cycles), _is_contradiction(false)  {

      for (typename linear_constraint_system_t::iterator it = csts.begin(); it != csts.end(); ++it) {
        linear_constraint_t cst = *it;
        if (cst.is_contradiction()) {
          this->_is_contradiction = true;
          return;
        } else if (cst.is_tautology()) {
          continue;
        } else {
          this->_cst_table.push_back(cst);
        }
      }
    }
    
    void run(CongruenceCollection& env) {
      if (this->_is_contradiction) {
        env.set_to_bottom();
      } 
      else {
        try {
          this->solve_system(env);
        }
        catch (bottom_found& e) {
          env.set_to_bottom();  
        }
      }
    }
    
  }; // class equality_congruence_solver



  template< typename Number, typename VariableName, unsigned int typeSize = -1 >
  class congruence_domain: public writeable, 
                           public numerical_domain< Number, VariableName >,
                           public bitwise_operators< Number, VariableName >, 
                           public division_operators< Number, VariableName > {

   public:
    typedef congruence< Number, typeSize >  congruence_t;
    typedef congruence_domain< Number, VariableName, typeSize > congruence_domain_t;

   public:
    // note that this is assuming that all variables have the same bit
    // width which is unrealistic.
    typedef variable< Number, VariableName > variable_t;
    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    
   private:
    typedef separate_domain< VariableName, congruence_t > separate_domain_t;
    typedef equality_congruence_solver< Number, VariableName, separate_domain_t, typeSize > solver_t;

   public:
    typedef typename separate_domain_t::iterator iterator;
     
   private:
    separate_domain_t _env;

   private:
    congruence_domain(separate_domain_t env): _env(env) { }
    
   public:
    static congruence_domain_t top() {
      return congruence_domain(separate_domain_t::top());
    }
    
    static congruence_domain_t bottom() {
      return congruence_domain(separate_domain_t::bottom());
    }
    
   public:
    congruence_domain(): _env(separate_domain_t::top()) { }

    congruence_domain(const congruence_domain_t& e): 
         writeable(), 
         numerical_domain< Number, VariableName >(), 
         bitwise_operators< Number, VariableName >(),
         division_operators< Number, VariableName >(),
         _env(e._env) { }

    congruence_domain_t& operator=(congruence_domain_t e) {
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

    bool operator<=(congruence_domain_t e) {
      return (this->_env <= e._env);
    }

    congruence_domain_t operator|(congruence_domain_t e) {
      return (this->_env | e._env);
    }

    congruence_domain_t operator&(congruence_domain_t e) {
      return (this->_env & e._env);
    }

    congruence_domain_t operator||(congruence_domain_t e) {
      return (this->_env || e._env);
    }

    congruence_domain_t operator&&(congruence_domain_t e) {
      return (this->_env && e._env);
    }

    void set(VariableName v, congruence_t i) {
      this->_env.set(v, i);
    }

    void set(VariableName v, Number n) {
      this->_env.set(v, congruence_t(n));
    }

    void operator-=(VariableName v) {
      this->_env -= v;
    }
    
    congruence_t operator[](VariableName v) {
      return this->_env[v];
    }

    congruence_t operator[](linear_expression_t expr) {
      congruence_t r(expr.constant());
      for (typename linear_expression_t::iterator it = expr.begin(); it != expr.end(); ++it) {
          congruence_t c(it->first);
          r = r + (c * this->_env[it->second.name()]);
      }
      return r;
    }
    
    void add(linear_constraint_system_t csts) {
      const std::size_t threshold = 10;
      if (!this->is_bottom()) {
	solver_t solver(csts, threshold);
	solver.run(this->_env);
      }
    }

    void operator+=(linear_constraint_system_t csts) {
      this->add(csts);
    }

    congruence_domain_t operator+(linear_constraint_system_t csts) {
      congruence_domain_t e(this->_env);
      e += csts;
      return e;
    }

    
    void assign(VariableName x, linear_expression_t e) {
      congruence_t r = e.constant();
      for (typename linear_expression_t::iterator it = e.begin(); it != e.end(); ++it) {
        r = r + (it->first * this->_env[it->second.name()]);
      }
      this->_env.set(x, r);
    }

    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
      congruence_t yi = this->_env[y];
      congruence_t zi = this->_env[z];
      congruence_t xi = congruence_t::bottom();
      
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
        case OP_DIVISION: {// signed division
          xi = yi / zi;
          break;
        }
      }
      this->_env.set(x, xi);
    }

    void apply(operation_t op, VariableName x, VariableName y, Number k) {

      congruence_t yi = this->_env[y];
      congruence_t zi(k);
      congruence_t xi = congruence_t::bottom();
      
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
        case OP_DIVISION: { // signed division
          xi = yi / zi;
          break;
        }
      }
      this->_env.set(x, xi);
    }

    // bitwise_operators_api
    
    void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width){
      congruence_t yi = this->_env[y];
      congruence_t xi = congruence_t::bottom();
      switch(op){
        case OP_TRUNC:
          xi = yi.Trunc(width);
          break;  
        case OP_ZEXT:
          xi = yi.ZExt(width);
          break;        
        case OP_SEXT:
          xi = yi.SExt(width);
          break;
        default: 
          throw error("unreachable");
      }      
      this->_env.set(x, xi);
    }

    void apply(conv_operation_t op, VariableName x, Number k, unsigned width){
      congruence_t yi(k);
      congruence_t xi = congruence_t::bottom();
      switch(op){
        case OP_TRUNC:
          xi = yi.Trunc(width);
          break;  
        case OP_ZEXT:
          xi = yi.ZExt(width);
          break;        
        case OP_SEXT:
          xi = yi.SExt(width);
          break;
        default: 
          throw error("unreachable");
      }      
      this->_env.set(x, xi);
    }
    void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z){
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
        default: 
          throw error("unreachable");
      }
      this->_env.set(x, xi);
    }
    
    void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k){
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
        default: 
          throw error("unreachable");
      }
      this->_env.set(x, xi);
    }
    
    // division_operators_api
    
    void apply(div_operation_t op, VariableName x, VariableName y, VariableName z){
      congruence_t yi = this->_env[y];
      congruence_t zi = this->_env[z];
      congruence_t xi = congruence_t::bottom();
      
      switch (op) {
        case OP_SDIV: {
	xi = yi / zi;
	break;
        }
        case OP_UDIV: {
	xi = yi.UDiv(zi);
	break;
        }
        case OP_SREM: {
	xi = yi.SRem(zi);
	break;
        }
        case OP_UREM: {
	xi = yi.URem(zi);
	break;
        }
        default: 
          throw error("unreachable");
      }
      this->_env.set(x, xi);

    }
    void apply(div_operation_t op, VariableName x, VariableName y, Number k){
      congruence_t yi = this->_env[y];
      congruence_t zi(k);
      congruence_t xi = congruence_t::bottom();
      switch (op) {
        case OP_SDIV: {
	xi = yi / zi;
	break;
        }
        case OP_UDIV: {
	xi = yi.UDiv(zi);
	break;
        }
        case OP_SREM: {
	xi = yi.SRem(zi);
	break;
        }
        case OP_UREM: {
	xi = yi.URem(zi);
	break;
        }
        default: 
          throw error("unreachable");
      }
      this->_env.set(x, xi);
    }


    std::ostream& write(std::ostream& o) {
      return this->_env.write(o);
    }

    const char* getDomainName () const {return "Congruences";}


  }; // class congruence_domain


} // namespace ikos

#endif // IKOS_CONGRUENCES_HPP
