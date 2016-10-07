#ifndef CRAB_COMBINED_DOMAINS_HPP
#define CRAB_COMBINED_DOMAINS_HPP

/************************************************************************** 
 * Combination of domains:
 * 
 * (1) reduced product of two arbitrary domains with only lattice
 *     operations.
 *
 * (2) reduced product of two arbitrary domains with all operations.
 * 
 * The reduction in (1) and (2) is simply done by making bottom the
 * abstract state if one of them is bottom.
 * 
 * (3) reduced product of two numerical domains. The reduction is done
 *     via a special push operation that must be defined by the
 *     domains. This is often more precise than (2).
 *
 * (4) reduced product of an arbitrary numerical domain and congruences.
 **************************************************************************/

#include <crab/common/types.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/congruences.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/operators_api.hpp>

namespace crab {
  namespace domains {

    using namespace ikos;


    // Provided by Ikos
    // Reduced product of two arbitrary domains with only lattice
    // operations.
    template< typename Domain1, typename Domain2 >
    class basic_domain_product2: public writeable {
    
     public:
      typedef basic_domain_product2< Domain1, Domain2 > basic_domain_product2_t;
      typedef Domain1 first_type;
      typedef Domain2 second_type;
      
     private:
      bool _is_bottom;
      Domain1 _first;
      Domain2 _second;
      
     public:
      
      static basic_domain_product2_t top() {
        return basic_domain_product2_t(Domain1::top(), Domain2::top());
      }
      
      static basic_domain_product2_t bottom() {
        return basic_domain_product2_t(Domain1::bottom(), Domain2::bottom());
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
      basic_domain_product2(): 
          _is_bottom(false), 
          _first(Domain1::top()), _second(Domain2::top()) { }
      
      basic_domain_product2(Domain1 first, Domain2 second): 
          _is_bottom(false), 
          _first(first), _second(second) {
        this->canonicalize();
      }
      
      basic_domain_product2(const basic_domain_product2_t& other): 
          writeable(), 
          _is_bottom(other._is_bottom), 
          _first(other._first), _second(other._second) { }
      
      basic_domain_product2_t& operator=(const basic_domain_product2_t& other) {
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
      
      bool operator<=(basic_domain_product2_t other) {
        if (this->is_bottom()) {
          return true;
        } else if (other.is_bottom()) {
          return false;
        } else {
          return (this->_first <= other._first) && 
              (this->_second <= other._second);
        }
      }
      
      bool operator==(basic_domain_product2_t other) {
        return (this->operator<=(other) && other.operator<=(*this));
      }
      
      void operator|=(basic_domain_product2_t other) {
        if (this->is_bottom()) {
          *this = other;
        } else if (other.is_bottom()) {
          return;
        } else {
          this->_first |= other._first;
          this->_second |= other._second;
        }
      }
      
      basic_domain_product2_t operator|(basic_domain_product2_t other) {
        if (this->is_bottom()) {
          return other;
        } else if (other.is_bottom()) {
          return *this;
        } else {
          return basic_domain_product2_t(this->_first | other._first, 
                                         this->_second | other._second);
        }
      }
      
      basic_domain_product2_t operator||(basic_domain_product2_t other) {
        if (this->is_bottom()) {
          return other;
        } else if (other.is_bottom()) {
          return *this;
        } else {
          return basic_domain_product2_t(this->_first || other._first, 
                                         this->_second || other._second);
        }
      }
      
      basic_domain_product2_t operator&(basic_domain_product2_t other) {
        if (this->is_bottom() || other.is_bottom()) {
          return bottom();
        } else {
          return basic_domain_product2_t(this->_first & other._first, 
                                         this->_second & other._second);
        }
      }
      
      basic_domain_product2_t operator&&(basic_domain_product2_t other) {
        if (this->is_bottom() || other.is_bottom()) {
          return bottom();
        } else {
          return basic_domain_product2_t(this->_first && other._first, 
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
    }; // class basic_domain_product2


    // Provided by Ikos
    // Reduced product of two arbitrary domains with all operations.
    template< typename Number, typename VariableName, typename Domain1, typename Domain2 >
    class domain_product2: public writeable,
                           public numerical_domain< Number, VariableName >,
                           public bitwise_operators< Number, VariableName >, 
                           public division_operators< Number, VariableName >,
                           public crab::domains::array_operators< Number, VariableName >,
                           public crab::domains::pointer_operators< Number, VariableName > {
     public:
      typedef domain_product2< Number, VariableName, Domain1, Domain2 > domain_product2_t;
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
      typedef basic_domain_product2< Domain1, Domain2 > basic_domain_product2_t;
      
     private:
      basic_domain_product2_t _product;
      
     private:
      domain_product2(basic_domain_product2_t product): _product(product) { }
      
      void reduce() {
        if (this->_product.first().is_bottom() || 
            this->_product.second().is_bottom()) {
          _product = basic_domain_product2_t::bottom();
        }
      }
      
     public:
      static domain_product2_t top() {
        return domain_product2_t(basic_domain_product2_t::top());
      }
      
      static domain_product2_t bottom() {
        return domain_product2_t(basic_domain_product2_t::bottom());
    }
      
     public:
      domain_product2(): _product() { }
      
      domain_product2(Domain1 first, Domain2 second): 
          _product(basic_domain_product2_t(first, second)) { }
      
      domain_product2(const domain_product2_t& other): 
          writeable(), 
          numerical_domain< Number, VariableName >(), 
          _product(other._product) { }
      
      domain_product2_t& operator=(const domain_product2_t& other) {
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
      
      bool operator<=(domain_product2_t other) {
        return (this->_product <= other._product);
      }
      
      bool operator==(domain_product2_t other) {
        return (this->_product == other._product);
      }
      
      void operator|=(domain_product2_t other) {
        this->_product |= other._product;
      }
      
      domain_product2_t operator|(domain_product2_t other) {
        return domain_product2_t(this->_product | other._product);
      }
      
      domain_product2_t operator&(domain_product2_t other) {
        return domain_product2_t(this->_product & other._product);
      }
      
      domain_product2_t operator||(domain_product2_t other) {
        return domain_product2_t(this->_product || other._product);
      }
      
      template<typename Thresholds>
      domain_product2_t widening_thresholds (domain_product2_t other,
                                             const Thresholds& ts) {
        return domain_product2_t (
            this->_product.first ().widening_thresholds (other._product.first (), ts),
            this->_product.second ().widening_thresholds (other._product.second (), ts));
      }
      
      domain_product2_t operator&&(domain_product2_t other) {
        return domain_product2_t(this->_product && other._product);
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
      
      virtual void array_assume (VariableName a, crab::variable_type a_ty, 
                                 VariableName var) override {
        this->_product.first().array_assume (a, a_ty, var);
        this->_product.second().array_assume (a, a_ty, var);
        this->reduce ();
      }
      
      virtual void array_load (VariableName lhs, VariableName a, crab::variable_type a_ty, 
                               linear_expression_t i, ikos::z_number bytes) override {
        
        this->_product.first().array_load (lhs, a, a_ty, i, bytes);
        this->_product.second().array_load (lhs, a, a_ty, i, bytes);
        this->reduce ();
      }
      
      virtual void array_store (VariableName a, crab::variable_type a_ty, 
                                linear_expression_t i, VariableName val, 
                                ikos::z_number bytes, bool is_singleton) override {
        this->_product.first().array_store (a, a_ty, i, val, bytes, is_singleton);
        this->_product.second().array_store (a, a_ty, i, val, bytes, is_singleton);
        this->reduce ();
      }
      
      virtual void array_assign (VariableName lhs, VariableName rhs, crab::variable_type ty) override {
        this->_product.first().array_assign (lhs, rhs, ty);
        this->_product.second().array_assign (lhs, rhs, ty);
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
        return basic_domain_product2_t::getDomainName ();
      }
      
    }; // class domain_product2
  

    // This domain is similar to domain_product2 but it uses a more
    // precise reduction operation.
    template<typename Domain1, typename Domain2 >
    class reduced_numerical_domain_product2: 
       public writeable,
       public numerical_domain<typename Domain1::number_t, typename Domain1::varname_t>,
       public bitwise_operators<typename Domain1::number_t, typename Domain1::varname_t>, 
       public division_operators<typename Domain1::number_t, typename Domain1::varname_t>,
       public array_operators<typename Domain1::number_t, typename Domain1::varname_t>,
       public pointer_operators<typename Domain1::number_t, typename Domain1::varname_t> {
      
     public:
      // Assume that Domain1 and Domain2 have the same types for
      // number_t and varname_t
      typedef typename Domain1::number_t Number;
      typedef typename Domain1::varname_t VariableName;

      using typename numerical_domain< Number, VariableName >::linear_expression_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
      using typename numerical_domain< Number, VariableName >::variable_t;
      using typename numerical_domain< Number, VariableName >::number_t;
      using typename numerical_domain< Number, VariableName >::varname_t;
      
      typedef interval<number_t> interval_t;
      typedef reduced_numerical_domain_product2<Domain1,Domain2> reduced_numerical_domain_product2_t;
      
     private:
      
      typedef patricia_tree_set<VariableName> variable_set_t;
      typedef domain_product2<Number, VariableName, Domain1, Domain2> domain_product2_t; 
      
      domain_product2_t _product;
      
      reduced_numerical_domain_product2(const domain_product2_t& product):
          _product(product) {}

      linear_constraint_system_t 
      to_linear_constraints (variable_t v, interval_t i) const {
        linear_constraint_system_t csts;
        if (i.lb ().is_finite () && i.ub ().is_finite ()) {
          auto lb = *(i.lb ().number());
          auto ub = *(i.ub ().number());
          if (lb == ub) {
            csts += (v == lb);
          } else {
            csts += (v >= lb);
            csts += (v <= ub);
          }
        } else if (i.lb ().is_finite ()) {
          auto lb = *(i.lb ().number());
          csts += (v >= lb);
        } else if (i.ub ().is_finite ()) {
          auto ub = *(i.ub ().number());
          csts += (v <= ub);
        }
        return csts;
      }

      void reduce_variable(const VariableName& v) {
        crab::CrabStats::count (getDomainName() + ".count.reduce");
        crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

        if (!is_bottom()) {
          
          // We just propagate from one domain to another.  We could
          // propagate in the other direction ... and repeat it for a
          // while ...

          Domain1& inv1 = _product.first ();
          Domain2& inv2 = _product.second ();

          // TODO?: having knowledge about the underlying domains we
          // might want to push only certain constraints (e.g.,
          // only equalities, only inequalities, etc). In that way, we
          // can avoid propagating redundant constraints.

          //////
          // propagate interval constraints between domains
          //////
          interval_t i1 = inv1[v];
          if (!i1.is_top ())
            inv2 += to_linear_constraints (variable_t(v), i1);
          interval_t i2 = inv2[v];
          if (!i2.is_top ())
            inv1 += to_linear_constraints (variable_t(v), i2);
          
          //////
          // propagate other constraints expressed by the domains
          //////

          crab::domains::product_domain_traits<Domain1,Domain2>::push (v, inv1, inv2);
          crab::domains::product_domain_traits<Domain2,Domain1>::push (v, inv2, inv1);
        }
      }
      
     public:
    
      static reduced_numerical_domain_product2_t top() {
        return reduced_numerical_domain_product2_t (domain_product2_t::top());
      }
      
      static reduced_numerical_domain_product2_t bottom() {
        return reduced_numerical_domain_product2_t(domain_product2_t::bottom());
      }
    
     public:
      
      reduced_numerical_domain_product2() : 
          writeable(),
          numerical_domain<number_t,varname_t>(),
          bitwise_operators<number_t,varname_t>(),
          division_operators<number_t,varname_t>(),
          array_operators<number_t,varname_t>(),
          pointer_operators<number_t,varname_t>(),
          _product() {}
      
      reduced_numerical_domain_product2(const reduced_numerical_domain_product2_t& other) :
          writeable(),
          numerical_domain<number_t,varname_t>(),
          bitwise_operators<number_t,varname_t>(),
          division_operators<number_t,varname_t>(),
          array_operators<number_t,varname_t>(),
          pointer_operators<number_t,varname_t>(),
          _product(other._product) { }
      
      reduced_numerical_domain_product2_t& 
      operator=(const reduced_numerical_domain_product2_t& other) {
        if (this != &other)
          this->_product = other._product;
        
         return *this;
      }
      
      bool is_bottom() { return this->_product.is_bottom(); }
      
      bool is_top() { return this->_product.is_top(); }
      
      Domain1& first() { return this->_product.first(); }
      
      Domain2& second() { return this->_product.second(); }
      
      bool operator<=(reduced_numerical_domain_product2_t other) {
        return this->_product <= other._product;
      }
      
      void operator|=(reduced_numerical_domain_product2_t other) {
        this->_product |= other._product;
      }
      
      reduced_numerical_domain_product2_t 
      operator|(reduced_numerical_domain_product2_t other) {
        return reduced_numerical_domain_product2_t(this->_product | other._product);
      }
      
      reduced_numerical_domain_product2_t 
      operator&(reduced_numerical_domain_product2_t other) {
        return reduced_numerical_domain_product2_t(this->_product & other._product);
      }
      
      reduced_numerical_domain_product2_t 
      operator||(reduced_numerical_domain_product2_t other) {
        return reduced_numerical_domain_product2_t(this->_product || other._product);
      }
      
      template<typename Thresholds>
      reduced_numerical_domain_product2_t widening_thresholds 
      (reduced_numerical_domain_product2_t other, const Thresholds& ts) {
        return reduced_numerical_domain_product2_t
            (this->_product.widening_thresholds (other._product, ts));
      }
      
      reduced_numerical_domain_product2_t 
      operator&&(reduced_numerical_domain_product2_t other) {
        return reduced_numerical_domain_product2_t(this->_product && other._product);
      }
      
      void set (varname_t v, interval_t x) {
        this->_product.first().set(v, x);
        this->_product.second().set(v, x);
      }
      
      interval_t operator[](varname_t v) {
        // We can choose either first or second domain
        return this->_product.second()[v];
      }
      
      void operator+=(linear_constraint_system_t csts) {
        this->_product += csts;
        for (auto v: csts.variables()) {
          reduce_variable(v.name());
        }
      }
      
      void operator-=(varname_t v) { this->_product -= v; }
            
      void assign(varname_t x, linear_expression_t e) {
        this->_product.assign(x, e);
        this->reduce_variable(x);
      }
      
      void apply(operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
      }
      
      void apply(operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      // bitwise_operators_api
      
      void apply(conv_operation_t op, varname_t x, varname_t y, unsigned width) {
        this->_product.apply(op, x, y, width);
        this->reduce_variable(x);
      }
      
      void apply(conv_operation_t op, varname_t x, number_t k, unsigned width) {
        this->_product.apply(op, x, k, width);
        this->reduce_variable(x);
      }
      
      void apply(bitwise_operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
      }
      
      void apply(bitwise_operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      // division_operators_api
      
      void apply(div_operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
      }
      
      void apply(div_operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      
      // domain_traits_api
      
      void expand(VariableName x, VariableName new_x) {
        crab::domains::domain_traits<Domain1>::expand (this->_product.first(), 
                                                       x, new_x);
        crab::domains::domain_traits<Domain2>::expand (this->_product.second(), 
                                                       x, new_x);
      }
      
      void normalize() {
        crab::domains::domain_traits<Domain1>::normalize(this->_product.first());
        crab::domains::domain_traits<Domain2>::normalize(this->_product.second());
      }
      
      template <typename Range>
      void forget(Range vars){
        crab::domains::domain_traits<Domain1>::forget(this->_product.first(), 
                                                      vars.begin (), vars.end());
        crab::domains::domain_traits<Domain2>::forget(this->_product.second(), 
                                                      vars.begin (), vars.end());
      }
      
      template <typename Range>
      void project(Range vars) {
        crab::domains::domain_traits<Domain1>::project(this->_product.first(), 
                                                       vars.begin(), vars.end());
        crab::domains::domain_traits<Domain2>::project(this->_product.second(), 
                                                       vars.begin(), vars.end());
      }
      
      void write(crab_os& o) { 
        this->_product.write(o); 
      }
      
      linear_constraint_system_t to_linear_constraint_system() {
        // We can choose either first or second domain
        return this->_product.second().to_linear_constraint_system();
      }
      
      static std::string getDomainName() { 
        std::string name = "ReducedProduct(" +
                           Domain1::getDomainName () +
                           "," + 
                           Domain2::getDomainName () + 
                           ")";
        return name;
      }
      
    }; // class reduced_numerical_domain_product2


   /*
    *  The reduce operator based on "Static Analysis of Arithmetical
    *  Congruences" by P. Granger published in International Journal of
    *  Computer Mathematics, 1989.
    */
    template < typename Number, int typeSize = -1 >
    class interval_congruence : public writeable {
     public:
      typedef interval_congruence< Number, typeSize > interval_congruence_t;
      
     private:
      typedef interval< Number > interval_t;
      typedef bound< Number > bound_t;
      typedef congruence< Number, typeSize > congruence_t;
      
     private:
      interval_t _first;
      congruence_t _second;
      
     private:
      interval_congruence(bool is_bottom)
          : writeable(),
            _first(is_bottom ? interval_t::bottom() : interval_t::top()),
            _second(is_bottom ? congruence_t::bottom() : congruence_t::top()) {}
      
     public:
      static interval_congruence_t top() { return interval_congruence(false); }
      
      static interval_congruence_t bottom() { return interval_congruence(true); }
      
     private:
      inline Number abs(Number x) { return x < 0 ? -x : x; }
      
      // operator % can return a negative number
      // mod(a, b) always returns a positive number
      inline Number mod(Number a, Number b) {
        Number m = a % b;
        if (m < 0)
          return m + b;
        else
          return m;
      }
      
      // R(c,a) is the least element of c greater or equal than a
      inline Number R(congruence_t c, Number a) {
        Number m = c.get().first;
        Number p = c.get().second;
        return a + mod(p - a, abs(m));
      }
      
      // L(c,a) is the greatest element of c smaller or equal than a
      inline Number L(congruence_t c, Number a) {
        Number m = c.get().first;
        Number p = c.get().second;
        return a - mod(a - p, abs(m));
      }
      
     public:
       interval_congruence(Number n)
           : _first(interval_t(n)), _second(congruence_t(n)) {}
      
      interval_congruence(interval_t i, congruence_t c) : 
          _first(i), _second(c) {
        this->reduce();
      }
      
      interval_congruence(interval_t i) : 
          _first(i), _second(congruence_t::top()) {
        this->reduce();
      }
      
      interval_congruence(congruence_t c) : 
          _first(interval_t::top()), _second(c) {
        this->reduce();
      }
      
      interval_congruence(const interval_congruence& other)
          : _first(other._first), _second(other._second) {}
      
      interval_congruence_t& operator=(interval_congruence_t other) {
        this->_first = other._first;
        this->_second = other._second;
        return *this;
      }
      
      bool is_bottom() {
        return this->_first.is_bottom() || this->_second.is_bottom();
      }
      
      bool is_top() { return this->_first.is_top() && this->_second.is_top(); }
      
      interval_t& first() { return this->_first; }
      
      congruence_t& second() { return this->_second; }
      
      /*
         Let (i,c) be a pair of interval and congruence these are the
         main rules described by Granger:

         if (c.is_bottom() || i.is_bottom()) (bottom(), bottom());
         if (c = 0Z+a and a \notin i)        (bottom(), bottom());
         if (c = 0Z+a)                       ([a,a]   , c);
         if (i=[a,b] and R(c,a) > L(c,b))    (bottom(), bottom());
         if (i=[a,b])                        ([R(c,a), L(c,b)], c);
         if (i=[a,+oo])                      ([R(c,a), +oo], c);
         if (i=[-oo,b])                      ([-oo, L(c,b)], c);
         otherwise                           (i,c)
       */
       
      void reduce() {
        interval_t& i = first();
        congruence_t& c = second();
        
        if (i.is_bottom() || c.is_bottom()) {
          i = interval_t::bottom();
          c = congruence_t::bottom();
        }
        
        // congruence is top and interval is a singleton
        if (c.is_top()) {
          boost::optional< Number > n = i.singleton();
          if (n) {
            c = congruence_t(*n);
          }
          return;
         }
        
        Number modulo = c.get().first;
        if (modulo == 0) {
          // congruence is a singleton so we refine the interval
          interval_t a(c.get().second);
          if (!(a <= i)) {
            i = interval_t::bottom();
             c = congruence_t::bottom();
          } else {
            i = a;
          }
        } else {
          // refine lower and upper bounds of the interval using
          // congruences
          bound_t lb = i.lb();
          bound_t ub = i.ub();
          
          if (lb.is_finite() && ub.is_finite()) {
            Number x = R(c, *(lb.number()));
            Number y = L(c, *(ub.number()));
            if (x > y) {
              i = interval_t::bottom();
              c = congruence_t::bottom();
            } else if (x == y) {
              i = interval_t(x);
              c = congruence_t(x);
            } else {
              i = interval_t(bound_t(x), bound_t(y));
            }
          } else if (lb.is_finite()) {
            Number x = R(c, *(lb.number()));
             i = interval_t(bound_t(x), bound_t::plus_infinity());
          } else if (ub.is_finite()) {
            Number y = L(c, *(ub.number()));
            i = interval_t(bound_t::minus_infinity(), bound_t(y));
          } else {
            // interval is top
          }
        }
      }
      
      void write(crab_os& o) {
        o << "(" << this->_first << ", " << this->_second << ")";
      }
      
     public:
      interval_congruence_t operator+(interval_congruence_t x) {
        return interval_congruence_t(this->_first.operator+(x.first()),
                                     this->_second.operator+(x.second()));
       }
      
      interval_congruence_t operator-(interval_congruence_t x) {
        return interval_congruence_t(this->_first.operator-(x.first()),
                                     this->_second.operator-(x.second()));
      }
      
      interval_congruence_t operator*(interval_congruence_t x) {
        return interval_congruence_t(this->_first.operator*(x.first()),
                                     this->_second.operator*(x.second()));
       }
      
      interval_congruence_t operator/(interval_congruence_t x) {
        return interval_congruence_t(this->_first.operator/(x.first()),
                                     this->_second.operator/(x.second()));
      }
      
      interval_congruence_t operator|(interval_congruence_t other) {
        return interval_congruence_t(this->_first | other._first,
                                     this->_second | other._second);
      }
      
      interval_congruence_t operator&(interval_congruence_t other) {
        return interval_congruence_t(this->_first & other._first,
                                     this->_second & other._second);
      }
      
     public:
      // division and remainder operations
      
      interval_congruence_t SDiv(interval_congruence_t x) {
        return interval_congruence_t(this->_first.SDiv(x.first()),
                                     this->_second.SDiv(x.second()));
      }
      
      interval_congruence_t UDiv(interval_congruence_t x) {
        return interval_congruence_t(this->_first.UDiv(x.first()),
                                     this->_second.UDiv(x.second()));
      }
      
      interval_congruence_t SRem(interval_congruence_t x) {
        return interval_congruence_t(this->_first.SRem(x.first()),
                                     this->_second.SRem(x.second()));
      }
      
      interval_congruence_t URem(interval_congruence_t x) {
        return interval_congruence_t(this->_first.URem(x.first()),
                                     this->_second.URem(x.second()));
       }
      
      // bitwise operations
      
      interval_congruence_t Trunc(unsigned width) {
        return interval_congruence_t(this->_first.Trunc(width),
                                     this->_second.Trunc(width));
      }
      
      interval_congruence_t ZExt(unsigned width) {
        return interval_congruence_t(this->_first.ZExt(width),
                                     this->_second.ZExt(width));
      }
      
      interval_congruence_t SExt(unsigned width) {
        return interval_congruence_t(this->_first.SExt(width),
                                     this->_second.SExt(width));
      }
      
      interval_congruence_t And(interval_congruence_t x) {
        return interval_congruence_t(this->_first.And(x.first()),
                                     this->_second.And(x.second()));
      }
      
      interval_congruence_t Or(interval_congruence_t x) {
        return interval_congruence_t(this->_first.Or(x.first()),
                                     this->_second.Or(x.second()));
      }
      
      interval_congruence_t Xor(interval_congruence_t x) {
        return interval_congruence_t(this->_first.Xor(x.first()),
                                     this->_second.Xor(x.second()));
      }
      
      interval_congruence_t Shl(interval_congruence_t x) {
        return interval_congruence_t(this->_first.Shl(x.first()),
                                     this->_second.Shl(x.second()));
      }
      
      interval_congruence_t LShr(interval_congruence_t x) {
        return interval_congruence_t(this->_first.LShr(x.first()),
                                     this->_second.LShr(x.second()));
       }
      
      interval_congruence_t AShr(interval_congruence_t x) {
        return interval_congruence_t(this->_first.AShr(x.first()),
                                     this->_second.AShr(x.second()));
      }
    };

    // Reduced product of a numerical domain with congruences.
    // It assumes that all variables have the same bitwdith which is
    // not realistic.
    template < typename NumAbsDom, int typeSize=-1 >
    class numerical_congruence_domain:
        public writeable,
        public numerical_domain<typename NumAbsDom::number_t, 
                                typename NumAbsDom::varname_t>,
        public bitwise_operators<typename NumAbsDom::number_t, 
                                 typename NumAbsDom::varname_t>,
        public division_operators<typename NumAbsDom::number_t,
                                  typename NumAbsDom::varname_t>,
        public array_operators<typename NumAbsDom::number_t,
                               typename NumAbsDom::varname_t>,
        public pointer_operators<typename NumAbsDom::number_t,
                                 typename NumAbsDom::varname_t> {
     private:
      typedef typename NumAbsDom::number_t N;
      typedef typename NumAbsDom::varname_t V;
      
     public:
      typedef numerical_congruence_domain< NumAbsDom, typeSize > rnc_domain_t;
      using typename numerical_domain< N, V>::linear_expression_t;
      using typename numerical_domain< N, V>::linear_constraint_t;
      using typename numerical_domain< N, V>::linear_constraint_system_t;
      using typename numerical_domain< N, V>::variable_t;
      using typename numerical_domain< N, V>::number_t;
      using typename numerical_domain< N, V>::varname_t;

      typedef congruence_domain<number_t, varname_t, typeSize> congruence_domain_t;
      typedef interval_congruence<number_t, typeSize> interval_congruence_t;
      typedef interval<number_t> interval_t;
      
     private:
      typedef patricia_tree_set<variable_t > variable_set_t;
      typedef domain_product2<number_t, varname_t, 
                                        NumAbsDom, congruence_domain_t > domain_product2_t; 
      
      domain_product2_t _product;
      
      numerical_congruence_domain(const domain_product2_t& product):
          _product(product) {}
      
      void reduce_variable(const varname_t& v) {
        crab::CrabStats::count (getDomainName() + ".count.reduce");
        crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

        if (is_bottom())
          return;
        
        auto i = this->_product.first()[v]; // project on intervals
        auto c = this->_product.second()[v];
        interval_congruence_t val(i, c);
        
        if (val.is_bottom()) {
          *this = bottom();
        } else {
          if (val.first() != i) { 
            // FIXME: this is imprecise for relational domains
            this->_product.first().set(v, val.first());
          }
          
          if (val.second() != c) 
            this->_product.second().set(v, val.second());
        }
      }
      
      void reduce_variables(variable_set_t variables) {
        for (typename variable_set_t::iterator it = variables.begin();
             !is_bottom() && it != variables.end();
              ++it) {
          this->reduce_variable((*it).name());
        }
      }
      
     public:
      
      static rnc_domain_t top() {
        return rnc_domain_t (domain_product2_t::top());
      }
      
      static rnc_domain_t bottom() {
        return rnc_domain_t(domain_product2_t::bottom());
      }
      
     public:
       
      numerical_congruence_domain() : 
          writeable(),
          numerical_domain<number_t,varname_t>(),
          bitwise_operators<number_t,varname_t>(),
          division_operators<number_t,varname_t>(),
          array_operators<number_t,varname_t>(),
          pointer_operators<number_t,varname_t>(),
          _product() {}
      
      numerical_congruence_domain(const rnc_domain_t& other) :
          writeable(),
          numerical_domain<number_t,varname_t>(),
          bitwise_operators<number_t,varname_t>(),
          division_operators<number_t,varname_t>(),
          array_operators<number_t,varname_t>(),
          pointer_operators<number_t,varname_t>(),
          _product(other._product) { }
      
      rnc_domain_t& operator=(const rnc_domain_t& other) {
        if (this != &other)
          this->_product = other._product;
         
        return *this;
      }
      
      bool is_bottom() { return this->_product.is_bottom(); }
      
      bool is_top() { return this->_product.is_top(); }
      
      NumAbsDom& first() { return this->_product.first(); }
      
      congruence_domain_t& second() { return this->_product.second(); }
      
      bool operator<=(rnc_domain_t other) {
        return this->_product <= other._product;
      }
      
      bool operator==(rnc_domain_t other) {
        return this->_product == other._product;
      }
      
      void operator|=(rnc_domain_t other) {
        this->_product |= other._product;
      }

      rnc_domain_t operator|(rnc_domain_t other) {
         return rnc_domain_t(this->_product | other._product);
      }
      
      rnc_domain_t operator&(rnc_domain_t other) {
        return rnc_domain_t(this->_product & other._product);
      }
      
      rnc_domain_t operator||(rnc_domain_t other) {
        return rnc_domain_t(this->_product || other._product);
      }
       
      template<typename Thresholds>
      rnc_domain_t widening_thresholds (rnc_domain_t other,
                                        const Thresholds& ts) {
        return rnc_domain_t(this->_product.widening_thresholds (other._product, ts));
      }
      
      rnc_domain_t operator&&(rnc_domain_t other) {
        return rnc_domain_t(this->_product && other._product);
      }
      
      // pre: x is already reduced
      void set (varname_t v, interval_congruence_t x) {
        this->_product.first().set(v, x.first());
        this->_product.second().set(v, x.second());
       }
      
      interval_congruence_t get(varname_t v) {
        return interval_congruence_t(this->_product.first()[v],
                                     this->_product.second()[v]);
      }
      
      interval_t operator[](varname_t v) {
        interval_congruence_t x = get (v);
        return x.first ();
      }
      
      void operator+=(linear_constraint_system_t csts) {
        this->_product += csts;
        this->reduce_variables(csts.variables());
      }
      
      void operator-=(varname_t v) { this->_product -= v; }
      
      void assign(varname_t x, linear_expression_t e) {
        this->_product.assign(x, e);
        this->reduce_variable(x);
      }
      
      void apply(operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
       }
      
      void apply(operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      // bitwise_operators_api
      
      void apply(conv_operation_t op, varname_t x, varname_t y, unsigned width) {
         this->_product.apply(op, x, y, width);
         this->reduce_variable(x);
      }
      
      void apply(conv_operation_t op, varname_t x, number_t k, unsigned width) {
        this->_product.apply(op, x, k, width);
        this->reduce_variable(x);
      }
      
      void apply(bitwise_operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
      }
      
      void apply(bitwise_operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      // division_operators_api
      
      void apply(div_operation_t op, varname_t x, varname_t y, varname_t z) {
        this->_product.apply(op, x, y, z);
        this->reduce_variable(x);
      }
      
      void apply(div_operation_t op, varname_t x, varname_t y, number_t k) {
        this->_product.apply(op, x, y, k);
        this->reduce_variable(x);
      }
      
      void write(crab_os& o) { 
        this->_product.write(o); 
       }
      
      linear_constraint_system_t to_linear_constraint_system() {
        return this->_product.first().to_linear_constraint_system();
      }
      
      static std::string getDomainName() { 
        return domain_product2_t::getDomainName (); 
      }
      
    }; // class numerical_congruence_domain

    template<typename Domain1, typename Domain2>
    class domain_traits <reduced_numerical_domain_product2<Domain1,Domain2> > {
     public:

      typedef reduced_numerical_domain_product2<Domain1,Domain2> product_t;
      typedef typename product_t::varname_t VariableName;

      static void normalize (product_t& inv) {
        inv.normalize();
      }

      static void expand (product_t& inv, VariableName x, VariableName new_x) {
        inv.expand (x, new_x);
      }
      
      template <typename Iter>
      static void forget (product_t& inv, Iter it, Iter end){
        inv.forget (boost::make_iterator_range (it, end));
      }
      
      template <typename Iter>
      static void project (product_t& inv, Iter it, Iter end) {
        inv.project (boost::make_iterator_range (it, end));
      }
    };

  } // end namespace domains
} // namespace crab

#endif 
