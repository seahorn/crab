/*******************************************************************************
 * Symbolic expressions 
 * 
 * Based on Section 5.1 from paper "Symbolic Methods to Enhance the
 * Precision of Numerical Abstract Domains" by Mine.
 * 
 ******************************************************************************/

#ifndef CRAB_SYMBOLIC_EXPRESSIONS_HPP
#define CRAB_SYMBOLIC_EXPRESSIONS_HPP

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <crab/domains/patricia_trees.hpp>

namespace crab { 

using namespace std;
using namespace ikos;

template< typename VariableName, typename Number >
class expression_visitor;

template <typename VariableName, typename Number>    
struct ExprMarshal;

template<typename VariableName, typename Number>
class _expression: public writeable{

 public:
  typedef _expression<VariableName, Number>  _expression_t;  
  typedef boost::shared_ptr< _expression_t > _expression_ptr;  
  typedef patricia_tree_set<VariableName> variable_set_t;

 public:
  typedef expression_visitor< VariableName, Number > expression_visitor_t;
  typedef boost::shared_ptr<expression_visitor_t> expression_visitor_ptr;

 protected:
  _expression(): writeable(){ }

 public:
  virtual size_t hash() const  = 0;
  virtual variable_set_t variables() const = 0;
  virtual void write(ostream& o) = 0;
  virtual void accept(expression_visitor_ptr visitor) = 0;
  virtual ~_expression() { }

}; // end class _expression


template<typename VariableName, typename Number>
class expression_number: public _expression<VariableName,Number> {
 public:
  typedef typename _expression<VariableName,Number>::expression_visitor_ptr expression_visitor_ptr;

 private:
  typedef typename _expression<VariableName,Number>::variable_set_t variable_set_t;

  Number _n;

 public:

  expression_number(Number n): 
      _expression<VariableName,Number>(), _n(n) { }

  size_t hash() const{
    boost::hash<string> hasher; 
    return hasher(boost::lexical_cast<string>(this->_n));
  }

  variable_set_t variables() const{
    return variable_set_t();
  }

  Number get_num () const {
    return _n;
  }
  
  void write(ostream& o){
    o << this->_n;
  }

  void accept(expression_visitor_ptr visitor) {  
    visitor->visit(*this);
  }

}; // end class expression_number

template<typename VariableName, typename Number>
class expression_variable: public _expression<VariableName,Number> {
 public:
  typedef typename _expression<VariableName,Number>::expression_visitor_ptr expression_visitor_ptr;

 private:
  typedef typename _expression<VariableName,Number>::variable_set_t variable_set_t;

  VariableName _v;

 public:
  expression_variable(VariableName v): 
      _expression<VariableName,Number>(), _v(v) { }

  size_t hash() const{
    boost::hash<string> hasher; 
    return hasher(boost::lexical_cast<string>(this->_v));
  }

  VariableName name() const{
    return this->_v;
  }

  variable_set_t variables() const{
    variable_set_t vars;
    vars += this->_v;
    return vars;
  }

  void write(ostream& o){
    o << this->_v;
  }

  void accept(expression_visitor_ptr visitor) {  
    visitor->visit(*this);
  }

}; // end class expression_variable

enum expression_op { add, sub, mul, udiv, sdiv, urem, srem, 
                     fadd, fsub, fmul, fdiv, frem, 
                     _shl, _lshr, _ashr, _and, _or, _xor };

inline std::ostream& operator<<(std::ostream& o, expression_op op){
  switch (op){
    case add: o << "+"; break;
    case sub: o << "-"; break;
    case mul: o << "*"; break;
    case sdiv: o << "/"; break;
    case udiv: o << "/_u"; break;
    case srem: o << "%"; break;
    case urem: o << "%_u"; break;
    case _shl: o << "<<"; break;
    case _lshr: o << ">>_l"; break;
    case _ashr: o << ">>"; break;
    case _and: o << "&"; break;
    case _or: o << "|"; break;
    case _xor: o << "^"; break;
    default: CRAB_ERROR ("symbolic expression operation not supported");
  }
  return o;
}

template<typename VariableName, typename Number>
class expression_binary_op: public _expression<VariableName, Number>{

 public:
  typedef typename _expression<VariableName,Number>::expression_visitor_ptr expression_visitor_ptr;

 private:
  typedef _expression<VariableName, Number> _expression_t;
  typedef typename _expression_t::_expression_ptr _expression_ptr;
  typedef typename _expression_t::variable_set_t  variable_set_t;

 private:
  _expression_ptr _left;
  _expression_ptr _right;
  expression_op   _op;

 public:
  expression_binary_op(expression_op op, 
                       _expression_ptr left, 
                       _expression_ptr right): 
      _expression<VariableName, Number>(), 
      _left(left), _right(right), _op(op){ }

  size_t hash() const{
    size_t res;
    boost::hash<int> hasher; 
    res = hasher(_op); 
    if (_left)
      boost::hash_combine (res, _left);
    if (_right)
      boost::hash_combine (res, _right);
    return res;
  }

  variable_set_t variables() const{
    return variable_set_t(_left->variables() | _right->variables());
  }

  expression_op get_op() const { return _op; }
  
  _expression_ptr get_left() const { return _left; }
  
  _expression_ptr get_right() const {  return _right; }

  virtual void accept(expression_visitor_ptr visitor) {  
    _left->accept(visitor);
    visitor->visit(*this);
    _right->accept(visitor);
  }

  void write(ostream& o){
    o << _op << "(" << *(_left) << "," << *(_right) << ")";
  }

}; // end expression_binary_op

template<typename VariableName, typename Number>
class expression;

template <typename VariableName, typename Number>    
struct ExprSubst {

  typedef expression_variable<VariableName, Number>  expression_variable_t;
  typedef expression_number<VariableName, Number>    expression_number_t;
  typedef expression_binary_op<VariableName, Number> expression_binary_op_t;
  typedef typename _expression<VariableName, Number>::_expression_ptr Expr;

  typedef boost::unordered_map<Expr, Expr> expr_expr_map;

  // replace in e all occurrences of x with y
  static Expr substitute (Expr e, VariableName x, Expr y) {
    expr_expr_map seen;
    return substitute (e, x, y, seen);
  }

  static Expr substitute (Expr e, VariableName x, Expr y, 
                          expr_expr_map &seen){

    auto  it = seen.find (e);
    if (it != seen.end ()) return it->second;

    if (boost::dynamic_pointer_cast<expression_number_t>(e)){
      return e;
    }
    
    if (boost::shared_ptr<expression_variable_t> v = 
        boost::dynamic_pointer_cast<expression_variable_t>(e)){
      if (v->name() == x)
        return y;
      else
        return e;
    }
    
    // otherwise must be a binary expression
    boost::shared_ptr< expression_binary_op_t > bin_e = 
        boost::static_pointer_cast<expression_binary_op_t>(e);
    Expr e1 = substitute (bin_e->get_left() , x, y, seen);
    Expr e2 = substitute (bin_e->get_right(), x, y, seen);

    Expr res(new expression_binary_op_t (bin_e->get_op(), e1, e2));
    seen.insert (typename expr_expr_map::value_type (e, res ));
    return res;
  }
};

template<typename VariableName, typename Number>
class expression : public writeable {

  friend class ExprMarshal<VariableName, Number>;

 public:
  typedef typename _expression<VariableName, Number>::variable_set_t  variable_set_t;

 private:
  typedef expression<VariableName, Number> expression_t;
  typedef expression_variable<VariableName, Number>  expression_var_t;
  typedef expression_number<VariableName, Number> expression_number_t;
  typedef expression_binary_op<VariableName, Number> expression_binary_op_t;

  typedef typename _expression<VariableName,Number>::_expression_ptr  _expression_ptr;

 private:
  _expression_ptr _ptr;
  bool _is_bottom;

  expression(bool is_bottom): 
      writeable(), _is_bottom(is_bottom)  { }
  
  expression(_expression_ptr ptr): 
      writeable(), _ptr(ptr), _is_bottom(false) { }

  _expression_ptr get () { return _ptr; }

 public:
  
  static expression_t top() {
    return expression_t(false);
  }

  static expression_t bottom() {
    return expression_t(true);
  }

 public:

  expression(VariableName v): 
      writeable(), 
      _ptr (_expression_ptr (new expression_var_t(v))), 
      _is_bottom(false) { }

  expression(Number n): 
      writeable(), 
      _ptr (_expression_ptr(new expression_number_t(n))), 
      _is_bottom(false) { }

  bool is_bottom() const { 
    return _is_bottom; 
  }

  bool is_top() const { 
    return (!is_bottom() && (!this->_ptr)); 
  }

  friend size_t hash_value(const expression_t &e){
    boost::hash<string> hasher; 
    if (e.is_top())
      return hasher("top");
    else if (e.is_bottom())
      return hasher("bottom");
    else
      return e._ptr->hash();
  }

  bool operator==(expression_t other) const {
    if (is_bottom() && other.is_bottom())
      return true;
    else if (is_top() && other.is_top())
      return true;
    else if (is_bottom() || is_top())
      return false;
    else if (other.is_bottom() || other.is_top())
      return false;
    else
      return (_ptr->hash() == other._ptr->hash());
  }
  
  bool operator<=(expression_t other) {
    if (other.is_top())
      return true;
    else if (is_bottom())
      return true;
    else
      return *this == other; 
  }

  expression_t operator| (expression_t other) {
    if (*this == other)
      return expression_t(this->_ptr);
    else 
      return top();
  }

  expression_t operator& (expression_t other) {
    return *this;
  }

  // Since we satisfy ACC we can replace widening with join
  expression_t operator|| (expression_t other){
    return *this | other;
  }

  // Since we satisfy ACC we can replace narrowing with meet
  expression_t operator&&(expression_t other){
    return *this & other;
  }

  variable_set_t occ() const {
    if (this->is_bottom()) {
      CRAB_ERROR ("symbolic expression's occurences of a bottom element");
    }
    else if (this->is_top()) {
      return variable_set_t(); //empty set
    }
    else {
      return this->_ptr->variables();
    }
  }

  expression_t substitute (VariableName v, expression_t e) const {
    if (is_bottom() || e.is_bottom()) {
      CRAB_ERROR ("symbolic expression substitution of a bottom element");
    }
    else if (is_top()){
      return top();
    }
    else if (e.is_top()){
      if (!e.occ()[v]) {
        return expression_t(this->_ptr);
      }
      else{
        return top();
      }
    }
    else{
      return expression_t(ExprSubst<VariableName, Number>::
                          substitute(this->_ptr, v, e._ptr ));
    }
  }

  void write (ostream& o){
    if (this->is_top()){
      o << "top";
    }
    else if (this->is_bottom()){
      o << "_|_";
    }
    else
      this->_ptr->write(o);
  }

  // Helper to build expressions
  expression_t combine(expression_op op, const expression_t &other) const {
                       
    if (is_top() || other.is_top())
      return top();
    else if (is_bottom() || other.is_bottom())
      return bottom();
    else
      return expression_t(_expression_ptr(new expression_binary_op_t(op, 
                                                                     _ptr, 
                                                                     other._ptr)));
  }


  void accept(boost::shared_ptr<expression_visitor< VariableName, Number > > visitor) {
    if (this->is_bottom())
      CRAB_ERROR ("symbolic expression accept method of a bottom element");
    
    if (!this->is_top())
      this->_ptr->accept(visitor);
  }
  
}; // end class expression

// To build common expressions
template<typename VariableName, typename Number>
expression<VariableName, Number> operator+(expression<VariableName,Number> x, 
                                           expression<VariableName,Number> y){
  return x.combine(add, y);
}

template<typename VariableName, typename Number>
expression<VariableName, Number> operator-(expression<VariableName,Number> x, 
                                           expression<VariableName,Number> y){
  return x.combine(sub, y);
}

template<typename VariableName, typename Number>
expression<VariableName, Number> operator*(expression<VariableName,Number> x, 
                                           expression<VariableName,Number> y){
  return x.combine(mul, y);
}

template<typename VariableName, typename Number>
expression<VariableName, Number> operator/(expression<VariableName,Number> x, 
                                           expression<VariableName,Number> y){
  return x.combine(sdiv, y);
}

template<typename VariableName, typename Number>
class expression_visitor {

  template< typename Any1, typename Any2 >
  friend class expression;
  typedef _expression<VariableName, Number> _expression_t;

 public:
  typedef expression<VariableName, Number> expression_t;
  typedef expression_binary_op<VariableName, Number> expression_binary_op_t;
  typedef expression_variable<VariableName, Number> expression_variable_t;
  typedef expression_number<VariableName, Number> expression_number_t;

  virtual ~expression_visitor() { }
  virtual void visit(expression_t &) = 0;
  virtual void visit(expression_binary_op_t &) = 0;
  virtual void visit(expression_variable_t &) = 0;
  virtual void visit(expression_number_t &) = 0;

 private:
  void visit(_expression_t & exp){ }

}; // end class expression_visitor


// Marshalling/UnMarshalling between linear and nonlinear expressions
template <typename VariableName, typename Number>    
struct ExprMarshal {
  
  typedef expression<VariableName, Number>  expression_t;
  typedef expression_number<VariableName, Number>  expression_number_t;
  typedef expression_variable<VariableName, Number>  expression_variable_t;
  typedef expression_binary_op<VariableName, Number>  expression_binop_t;
  typedef linear_expression< Number, VariableName > linear_expression_t;
  typedef typename linear_expression_t::variable_t variable_t;

  static expression_t marshal (const linear_expression_t &e) {
    if (e.is_constant())
      return expression_t(e.constant());
    else {
      expression_t y (e.constant());
      for (auto l : e) {      
        Number n     = l.first;
        variable_t v = l.second;
           
        if (n == 1){
          y = y + expression_t (v.name()); 
        }
        else if (n > 1){
          y = y + expression_t (expression_t(n) * expression_t(v.name())); 
        }
        else if (n == -1){
          y = y - expression_t (v.name()); 
        }
        else if (n < -1){
          Number abs_n = n * Number(-1);
          y = y - expression_t (expression_t(abs_n) * expression_t(v.name())); 
        }
      }
      return y;
    }

  }

  typedef boost::unordered_map <expression_t, boost::optional<linear_expression_t> > cache_t;
  
  static boost::optional<linear_expression_t> translate (expression_t e, cache_t &seen) {

    auto it = seen.find (e);
    if (it != seen.end ()) {
      return it->second;
    }
    
    
    boost::optional<linear_expression_t> res;
    
    auto x = e.get ();
    if (expression_number_t * n = dynamic_cast<expression_number_t*>  (&(*x))) {
      res = linear_expression_t(n->get_num ());
    }
    else if (expression_variable_t * v = dynamic_cast<expression_variable_t*>  (&(*x))) {
      res = linear_expression_t(variable_t (v->name ()));
    }
    else {
      expression_binop_t * f = static_cast<expression_binop_t*>  (&(*x));
      if (f->get_op () == add || f->get_op () == sub || f->get_op () == mul){
        expression_t left (f->get_left ());
        expression_t right (f->get_right ());
        auto e1 = translate (left, seen);
        auto e2 = translate (right, seen);
        if (e1 && e2) {
          if (f->get_op () == add) 
            res = *e1 + *e2;
          else if (f->get_op () == sub) 
            res = *e1 - *e2;
          else {
            if((*e1).is_constant ())
              res = (*e1).constant () * (*e2);
            else if ((*e2).is_constant ())
              res = (*e1) * (*e2).constant (); 
          }
        }
      }
    }
    seen.insert (make_pair (e,res));
    return res;
  }

  static boost::optional<linear_expression_t> unmarshal (expression_t e) {
    cache_t seen;
    return translate (e, seen);
  }

};

} 

#endif 
