/*******************************************************************************
 *
 * Symbolic Constant Domain
 *
 * Based on the paper "Symbolic Methods to Enhance the Precision of
 * Numerical Abstract Domains" by A. Mine, VMCAI'06.
 *
 ******************************************************************************/

#ifndef SYMBOLIC_CONSTANT_DOMAIN_HPP
#define SYMBOLIC_CONSTANT_DOMAIN_HPP

// Uncomment for enabling debug information
// #include <crab/common/dbg.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

#include <crab/common/types.hpp>
#include <crab/domains/symexpr/symbolic_expressions.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/domain_products.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>
#include <crab/domains/division_operators_api.hpp>


namespace crab { 

 namespace domains {

   using namespace ikos;
   using namespace std;

   template<typename Number, typename VariableName>
   class symbolic_constant_domain: public writeable,
                                   public numerical_domain< Number, VariableName >,
                                   public bitwise_operators< Number, VariableName >, 
                                   public division_operators< Number, VariableName > {

    public:
     using typename numerical_domain< Number, VariableName >::linear_expression_t;
     using typename numerical_domain< Number, VariableName >::linear_constraint_t;
     using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
     using typename numerical_domain< Number, VariableName >::variable_t;
     using typename numerical_domain< Number, VariableName >::number_t;
     using typename numerical_domain< Number, VariableName >::varname_t;

    public:
     typedef expression<VariableName, Number> expression_t;
     typedef symbolic_constant_domain< Number, VariableName > sym_constant_domain_t;
     typedef typename expression_t::variable_set_t expr_variable_set_t;

    private:
     typedef separate_domain< VariableName, expression_t > separate_domain_t;
     
    public:
     typedef typename separate_domain_t::iterator iterator;
     
    private:

     // --- We ensure that there is no cyclic dependencies in _env.
     //     We have a cycle if given variable v1,...,vn, for all i. vi
     //     \in occ (_env [v_i+1]),..., vn \in occ (_env [v1])
     separate_domain_t _env;

     symbolic_constant_domain(separate_domain_t env): _env(env) { }
     
    public:
     
     static sym_constant_domain_t top() {
       return symbolic_constant_domain(separate_domain_t::top());
     }
     
     static sym_constant_domain_t bottom() {
       return symbolic_constant_domain(separate_domain_t::bottom());
     }
     
    public:

     template<typename T>
     static expression_t 
     mkExp (bitwise_operation_t op, VariableName x, T y) {
       expression_t e1 (x);
       expression_t e2 (y);
       switch (op) {
         case OP_AND:  return e1.combine(_and, e2); break;
         case OP_OR:   return e1.combine(_or, e2); break;
         case OP_XOR:  return e1.combine(_xor, e2); break;
         case OP_SHL:  return e1.combine(_shl, e2); break;
         case OP_LSHR: return e1.combine(_lshr, e2); break;
         default:      return e1.combine(_ashr, e2); break;
       }
     }

     template<typename T>
     static expression_t 
     mkExp (operation_t op, VariableName x, T y) { 
       expression_t e1 (x);
       expression_t e2 (y);
       switch (op) {
         case OP_ADDITION: return e1.combine(add, e2); break; 
         case OP_SUBTRACTION: return e1.combine(sub, e2); break;
         case OP_MULTIPLICATION: return e1.combine(mul, e2); break;
         default: return e1.combine(sdiv, e2); break;
       }
     }

     template<typename T>
     static expression_t
     mkExp (div_operation_t op, VariableName x, T y) {
       expression_t e1 (x);
       expression_t e2 (y);
       switch (op) {
         case OP_SDIV: return e1.combine (sdiv, e2); break;
         case OP_UDIV: return e1.combine (udiv, e2); break;
         case OP_SREM: return e1.combine (srem, e2); break;
         default:      return e1.combine (urem, e2); break;
       }
     }

    public:
  
     symbolic_constant_domain(): 
         _env(separate_domain_t::top()) { }
     
     symbolic_constant_domain(const sym_constant_domain_t& e): 
         writeable(), _env(e._env) { }
     
     sym_constant_domain_t& operator=(sym_constant_domain_t e) {
       _env = e._env;
       return *this;
     }
     
     iterator begin() { return _env.begin(); }
     
     iterator end() { return _env.end(); }
     
     bool is_bottom() { return _env.is_bottom(); }
     
     bool is_top() { return _env.is_top(); }

     bool operator<=(sym_constant_domain_t e) {
       return (_env <= e._env);
     }
  
     sym_constant_domain_t operator|(sym_constant_domain_t e) {
       return (_env | e._env);
     }
     
     sym_constant_domain_t operator&(sym_constant_domain_t e) {
       return (_env & e._env);
     }
  
     sym_constant_domain_t operator||(sym_constant_domain_t e) {
       return (_env || e._env);
     }
     
     sym_constant_domain_t operator&&(sym_constant_domain_t e) {
       return (_env && e._env);
     }
     
     void set(VariableName v, expression_t e) {
       _env.set(v, e);
     }
     
     void operator-=(VariableName v) {
       _env -= v;
     }
  
     expression_t operator[](VariableName v) {
       return _env[v];
     }
     
     void apply (VariableName x, expression_t e) {
       if (is_bottom ())
         return;

       // update x
       set (x, e.substitute( x, this->operator[](x)));

       // update the rest of variables with the new value of x
       expression_t y = this->operator[](x);
       for (auto p: *this) {
         if (!(p.first == x)) {
           set (p.first, p.second.substitute (x, y));
         }
       }
     }
     
     // --- numerical_domains_api
     
     void operator+=(linear_constraint_t csts) { }

     void operator+=(linear_constraint_system_t csts) { }

     void assign(VariableName x, linear_expression_t e) { 
       if (is_bottom()) 
         return;
       
       auto vars = e.variables ();
       typedef typename linear_expression_t::variable_set_t variable_set_t;
       if (variable_set_t (variable_t (x)) <= vars)
         return;

       apply  (x, ExprMarshal <VariableName, Number>::marshal (e));
     }

     void apply(operation_t op, VariableName x, VariableName y, VariableName z) { 
       if (is_bottom ())
         return;

       if (x == y || x == z)
         return;

       apply (x, mkExp (op, y, z));
     }
     
     void apply(operation_t op, VariableName x, VariableName y, Number k) { 
       if (is_bottom ())
         return;

       if (x == y)
         return;

       apply (x, mkExp (op, y, k));
     }

     // --- bitwise_operators_api

     void apply(conv_operation_t op, 
                VariableName x, VariableName y, unsigned width) { 
       // we ignore width since we assume unlimited precision
       assign (x, linear_expression_t (y));
     }
     
     void apply(conv_operation_t op, 
                VariableName x, Number k, unsigned width){ 
       // we ignore width since we assume unlimited precision
       assign (x, linear_expression_t (k));
     }
     
     void apply(bitwise_operation_t op, 
                VariableName x, VariableName y, VariableName z) { 
       if (is_bottom ())
         return;

       if (x == y || x == z)
         return;

       apply (x, mkExp (op, y, z));
     }
    
     void apply(bitwise_operation_t op, 
                VariableName x, VariableName y, Number k) { 
       if (is_bottom () || x==y) 
         return;
         
       apply (x, mkExp (op, y, k));
     }

     // --- division_operators_api

     void apply(div_operation_t op, 
                VariableName x, VariableName y, VariableName z) { 
       if (is_bottom ())
         return;

       if (x == y || x == z)
         return;

       apply (x, mkExp (op, y, z));
     }

     void apply(div_operation_t op, 
                VariableName x, VariableName y, Number k) { 
       if (is_bottom ())
         return;

       if (x == y)
         return;

       apply (x, mkExp (op, y, k));
     }

     void write(ostream& o) { return _env.write(o); }
     
     const char* getDomainName () const { 
       return "Symbolic constants";
     }
     
   }; // end class symbolic_constant_domain


   // Reduced product of symbolic constant domain and a numerical
   // domain.
   //
   // FIXME: although NumDomain can be any numerical domain the
   // symbolic substitutions quickly yield assignments or arithmetic
   // operations with multiple variables (more than two) which domains
   // like DBM cannot handle unless we do some preprocessing (i.e.,
   // split into several operations).
   template< typename NumDomain, typename Number, typename VariableName>
   class num_sym_constant_domain: public writeable, 
                                  public numerical_domain< Number, VariableName >,
                                  public bitwise_operators< Number, VariableName >, 
                                  public division_operators< Number, VariableName > {
     
    public:
     using typename numerical_domain< Number, VariableName >::linear_expression_t;
     using typename numerical_domain< Number, VariableName >::linear_constraint_t;
     using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
     using typename numerical_domain< Number, VariableName >::variable_t;
     using typename numerical_domain< Number, VariableName >::number_t;
     using typename numerical_domain< Number, VariableName >::varname_t;
     typedef num_sym_constant_domain < NumDomain, Number, VariableName > num_sym_constant_domain_t;
     
    private:
     typedef expression<VariableName, Number> expression_t;
     typedef typename expression_t::variable_set_t expr_variable_set_t;
     typedef symbolic_constant_domain< Number, VariableName> sym_constant_domain_t;
     typedef numerical_domain_product2< Number, VariableName, 
                                        NumDomain, sym_constant_domain_t > domain_product2_t;
     typedef interval< Number > interval_t;
     typedef ExprMarshal <VariableName, Number> expr_marshal_t;

     domain_product2_t _product;
  
    private:
     num_sym_constant_domain (domain_product2_t product): 
         _product(product) { }
     
     // Substitution strategy (Section 5.3 from VMCAI'06 paper)
     expression_t  strat (expression_t e) {
       expression_t r (e);
       auto variables = e.occ();
       for(auto v: variables) { 
         expression_t y = _product.second () [v];
         if (!y.is_top ()) {
           r = r.substitute(v, y);
         }
       }
       return r;
     }

    public:
     
     static num_sym_constant_domain_t top() {
       return num_sym_constant_domain_t(domain_product2_t::top());
     }
     
     static num_sym_constant_domain_t bottom() {
       return num_sym_constant_domain_t(domain_product2_t::bottom());
     }
     
    public:
     
     num_sym_constant_domain(): _product() { }
     
     num_sym_constant_domain (const num_sym_constant_domain_t& other): 
         writeable(), 
         numerical_domain< Number, VariableName >(), 
         bitwise_operators< Number, VariableName >(),
         division_operators< Number, VariableName > (),
         _product(other._product) { }
     
     num_sym_constant_domain_t& 
     operator=(num_sym_constant_domain_t other) {
       this->_product = other._product;
       return *this;
     }
     
     bool is_bottom() {
       return _product.is_bottom();
     }
     
     bool is_top() {
       return _product.is_top();
     }
     
     NumDomain& first() {
       return _product.first();
     }
     
     sym_constant_domain_t& second() {
       return _product.second();
     }
  
     bool operator<=(num_sym_constant_domain_t other) {
       return (_product <= other._product);
     }

     bool operator==(num_sym_constant_domain_t other) {
       return (_product == other._product);
     }
  
     num_sym_constant_domain_t operator|(num_sym_constant_domain_t other) {
       return num_sym_constant_domain_t(_product | other._product);
     }
     
     num_sym_constant_domain_t operator&(num_sym_constant_domain_t other) {
       return num_sym_constant_domain_t(_product & other._product);
     }
     
     num_sym_constant_domain_t operator||(num_sym_constant_domain_t other) {
       return num_sym_constant_domain_t(_product || other._product);
     }
     
     num_sym_constant_domain_t operator&&(num_sym_constant_domain_t other) {
       return num_sym_constant_domain_t (_product && other._product);
     }
  
     void set(VariableName v, interval_t x) {
       _product.first().set(v, x);

       _product.second() -= v;
       boost::optional <Number> n = x.singleton ();
       if (n) {
         _product.second ().set (v, expression_t (*n));
       }
     }

     interval_t operator[](VariableName v) {
       return _product.first()[v].first();
     }
     
     void operator+=(linear_constraint_t cst) {
       if (is_bottom ())
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       linear_expression_t x = cst.expression ();
       expression_t e1 = expr_marshal_t::marshal (x);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> y = expr_marshal_t::unmarshal (e2)) {
         linear_constraint_t subst_cst (*y, cst.kind ());
         _product += subst_cst;
         CRAB_DEBUG("Added cst ", subst_cst, ": ", *this);
       }
       //else {
       if (!is_bottom ()) {
         _product += cst;
         CRAB_DEBUG("Added cst ", cst, ": ", *this);
       }
       // }
     }
  
     void operator+=(linear_constraint_system_t csts) {
       if (is_bottom ()) 
         return;
       
       for (auto cst : csts) {
         this->operator+= (cst);
       }
     }

     void operator-=(VariableName v) {
       _product -= v;
     }
     
     void assign(VariableName x, linear_expression_t e) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = ExprMarshal <VariableName, Number>::marshal (e);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> y = expr_marshal_t::unmarshal (e2)) {
         _product.assign (x, *y);
         CRAB_DEBUG("After ",x, ":=", *y, ": ", *this);
       }
       else {
         _product.assign (x, e);
         CRAB_DEBUG("After ",x, ":=", e, ": ", *this);
       }
     }
     
     void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, z);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2)) {
         _product.assign (x, *w);
         CRAB_DEBUG("After ",x, ":=", *w, ": ", *this);
       }
       else {
         _product.apply(op, x, y, z);
         CRAB_DEBUG("After ",x, ":=", y, op, z ,": ", *this);
       }
     }
     
     void apply(operation_t op, VariableName x, VariableName y, Number k) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, k);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2)) {
         _product.assign (x, *w);
         CRAB_DEBUG("After ",x, ":=", *w, ": ", *this);
       }
       else {
         _product.apply(op, x, y, k);
         CRAB_DEBUG("After ",x, ":=", y, op, k ,": ", *this);
       }
     }
     
     // bitwise_operators_api
     
     void apply(conv_operation_t op, 
                VariableName x, VariableName y, unsigned width) {
       assign (x, variable_t (y));
     }
     
     void apply(conv_operation_t op, 
                VariableName x, Number k, unsigned width) {
       assign (x, k);
     }
     
     void apply(bitwise_operation_t op, 
                VariableName x, VariableName y, VariableName z) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, z);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2))
         _product.assign (x, *w);
       else
         _product.apply(op, x, y, z);
     }
     
     void apply(bitwise_operation_t op, 
                VariableName x, VariableName y, Number k) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, k);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2))
         _product.assign (x, *w);
       else
         _product.apply(op, x, y, k);
     }
     
     // division_operators_api
     
     void apply(div_operation_t op, 
                VariableName x, VariableName y, VariableName z) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, z);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2))
         _product.assign (x, *w);
       else
         _product.apply(op, x, y, z);
     }
     
     void apply(div_operation_t op, 
                VariableName x, VariableName y, Number k) {
       if (is_bottom ()) 
         return;

       // --- reduction between the symbolic constant domain and the
       //     numerical domain
       expression_t e1 = sym_constant_domain_t::mkExp (op, y, k);
       expression_t e2 = strat (e1);
       if (boost::optional<linear_expression_t> w = expr_marshal_t::unmarshal (e2))
         _product.assign (x, *w);
       else
         _product.apply(op, x, y, k);
     }

     linear_constraint_system_t to_linear_constraint_system () {
       linear_constraint_system_t csts = 
           _product.first ().to_linear_constraint_system ();
       
       // Add equivalences from the symbolic domain
       for (auto p:_product.second ()) {
         if (boost::optional<linear_expression_t> e = 
             expr_marshal_t::unmarshal (p.second)) {
           csts += (p.first == *e);
         }
       } 
       return csts;
     }
     
     void write(ostream& o) {
       _product.write(o);
     }

     const char* getDomainName () const { 
       return _product.getDomainName ();
     }     

   }; 

 } // end namespace
} // end namespace

#endif 
