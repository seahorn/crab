#include "crab/config.h"

#ifndef HAVE_MDD
int main(int argc, char**argv) {
  return 0;
}
#else
#include "include/MDD.hh"
#include "include/MDD.hpp"
#include "include/mdd_builtin_cache.hh"
#include "include/mdd_ext_cache.hh"
#include "include/MDD_ops.hh"
#include "include/MDD_vis.hh"
#include "include/MDD_visit.hh"
#include "include/MDD_arith.hh"
#include "include/interval.hh"

#include <boost/optional.hpp>
#include <vector>

using namespace mdd_boxes;

typedef mdd_mgr<int64_t> mgr_t;
typedef mdd_node<int64_t> node_t;
typedef mdd_ref<int64_t> ref_t;

template<class T>
struct ShowNop {
  T operator()(T k) { return k; }
};


template<class V, class ShowV>
class MDD_print {
public:
  
  typedef mdd_mgr<V>  mgr_t;
  typedef mdd_node<V> node_t;

private:
  
  class interval {
    boost::optional<V> lb;
    boost::optional<V> ub;

    // create an interval [lb, ub)
    interval(boost::optional<V> _lb, boost::optional<V> _ub)
      : lb(_lb), ub(_ub) {}
    
    public:
    // create an interval [lb, ub)
    interval(V _lb, V _ub)
      : lb(_lb), ub(_ub) {}

    // return (-oo, ub]
    static interval lower_half_line(V ub) {
      return interval(boost::optional<V>(), ub);
    }

    // return [lb, +oo)
    static interval upper_half_line(V lb) {
      return interval(lb, boost::optional<V>());
    }

    void write(std::ostream& o) const {
      if (lb && ub) {
	o << "[" << *lb << "," << *ub << ")";
      } else if (!lb && ub) {
	o << "(" << "-oo" << "," << *ub << "]";
      } else  {
	assert (lb && !ub);
	o << "[" << *lb << "," << "+oo" << ")";
      }
    }
  };

  typedef std::vector<std::pair<node_t*, interval>> stack_t;

  void rec_mdd_print(node_t* n, stack_t& stack) {
    
    if (n == MDD_ops<V>::MDD_TRUE) {
      // end of the path: we print the stack
      for (auto it = stack.begin(), et = stack.end(); it!=et; ) {
	auto p = *it;
	o << show_v(p.first->var);
	o << "=";
	p.second.write(o);
	++it;
	if (it != et) {
	  o << " and ";
	}
      }
      o << "\n";
      return;
    }
    
    auto it = n->edges.begin();
    auto en = n->edges.end();
    if(n->dest0 != MDD_ops<V>::MDD_FALSE) {
      // there is an edge  from n to dest0 with interval (-oo, (*it).lb]
      stack.push_back(std::make_pair(n, interval::lower_half_line((*it).lb)));
      rec_mdd_print(n->dest0, stack);
      stack.pop_back();      
    }
    
    V lb((*it).lb);
    node_t* dest((*it).dest);    
    for(++it; it != en; ++it) {
      if(dest != MDD_ops<V>::MDD_FALSE) {
	// there is an edge from n to dest with interval [lb, (*it).lb)
	stack.push_back(std::make_pair(n, interval(lb, (*it).lb)));	
	rec_mdd_print(dest, stack);
	stack.pop_back();
      }
      lb = (*it).lb;
      dest = (*it).dest;
    }

    if (dest != MDD_ops<V>::MDD_FALSE) {
      // there is an edge from n to dest with interval [(*it).lb, +oo)
      stack.push_back(std::make_pair(n, interval::upper_half_line(lb)));      
      rec_mdd_print(dest, stack);
      stack.pop_back();
    }
  }
  
public:
  MDD_print(std::ostream& _o, ShowV& _show_v)
    : o(_o), show_v(_show_v) { }
  
  void operator()(node_t* n) {
    if(n == MDD_ops<V>::MDD_FALSE) {
      o << "_|_";
    } else if(n == MDD_ops<V>::MDD_TRUE) {
      o << "top";
    } else {
      stack_t stack;
      rec_mdd_print(n, stack);
    }
  }

protected:
  std::ostream& o;
  ShowV& show_v;
};

ref_t make_range(mgr_t& m, size_t v, int64_t l, int64_t u) {
  auto b(m.begin_node(v, m.mdd_false().p));
  b.push(l, m.mdd_true().p);
  b.push(u, m.mdd_false().p);

  node_t* r(b.finalize());
  return ref_t(&m, r);
}

template<typename V,typename R>
std::ostream& operator<<(std::ostream& o, num_interval<V,R> i) {
  if (i.is_empty()) {
    o << "_|_";
  } else if (i.is_top()) {
    o << "[-oo,+oo]";
  } else if (i.lb_is_finite() && i.ub_is_finite()) {
    o << "[" << i.lb() << "," << i.ub() << ")";
  } else if (!i.ub_is_finite()) {
    o << "[" << i.lb() << ", +oo]";    
  } else {
    o << "[-oo," << i.ub() << "]";        
  }
  return o;
}

template<typename V,typename R>
std::ostream& operator<<(std::ostream& o, num_interval_NE<V,R> i) {
  if (i.is_top()) {
    o << "[-oo,+oo]";
  } else if (i.lb_is_finite() && i.ub_is_finite()) {
    o << "[" << i.lb() << "," << i.ub() << ")";
  } else if (!i.ub_is_finite()) {
    o << "[" << i.lb() << ", +oo]";    
  } else {
    o << "[-oo," << i.ub() << "]";        
  }
  return o;
}

typedef MDD_dot<int64_t, ShowNop<int64_t>> MDD_dot_vis;
typedef MDD_print<int64_t, ShowNop<int64_t>> MDD_printer_t;

template<typename MDD_printer>
void test_meet(MDD_printer& v, mgr_t& m) {
  std::cout << "== Begin test_meet\n";  
  ref_t p(m.var_lt(10, 100));
  ref_t q(m.var_ge(10, 50));
  ref_t r(m.var_ge(15, 80));

  ref_t x(&m, MDD_ops<int64_t>::meet(&m, p.get(), r.get()));
  std::cout << "\"meet: v_10 < 100 and v_15 >= 80\"\n";
  v(x.get());
  std::cout << "\"v_10 >= 50\"\n";
  v(q.get());
  ref_t y(&m, MDD_ops<int64_t>::meet(&m, x.get(), q.get()));

  //fprintf(stdout, "Conj: %p\n", y.get());
  std::cout << "\"meet v_10 < 100 and v_15 >= 80 and v_10 >= 50\"\n";
  v(y.get());
  std::cout << "== End test_meet\n";    
}

template<typename MDD_printer>
void test_widen(MDD_printer& v, mgr_t& m) {
  std::cout << "== Begin test_widen\n";
  ref_t x1(make_range(m, 0, 0, 2));
  ref_t x2(make_range(m, 0, 0, 3));
  ref_t x3(make_range(m, 0, 4, 6));

  ref_t y1(make_range(m, 1, 0, 2));
  ref_t y2(make_range(m, 1, 0, 3));
  ref_t y3(make_range(m, 1, 2, 5));

  ref_t p1(&m, MDD_ops<int64_t>::meet(&m, x1.get(), y1.get()));
  ref_t p2(&m, MDD_ops<int64_t>::meet(&m, x3.get(), y3.get()));
  ref_t q1(&m, MDD_ops<int64_t>::meet(&m, x2.get(), y2.get()));

  ref_t p(&m, MDD_ops<int64_t>::join(&m, p1.get(), p2.get()));
  ref_t q(&m, MDD_ops<int64_t>::join(&m, q1.get(), p2.get()));
  std::cout << "\"(0<= v_0 < 2 and 0<= v_1 < 2) or (4 <= v_0 < 6 and 2 <= v_1 < 5)\"\n";
  v(p.get());
  std::cout << "\"(0<= v_0 < 3 and 0<= v_1 < 3) or (4 <= v_0 < 6 and 2 <= v_1 < 5)\"\n";  
  v(q.get());
  std::cout << "\"widening of the last two:\"\n";
  ref_t w(&m, MDD_ops<int64_t>::widen(&m, p.get(), q.get()));
  v(w.get());
  std::cout << "== End test_widen\n";  
}

template<typename MDD_printer>
void test_project(MDD_printer& v, mgr_t& m) {
  std::cout << "== Begin test_project\n";    
  
  ref_t x1(make_range(m, 3, 0, 2));
  ref_t x3(make_range(m, 3, 4, 6));

  ref_t y1(make_range(m, 5, 0, 2));
  ref_t y3(make_range(m, 5, 3, 5));

  ref_t p1(&m, MDD_ops<int64_t>::meet(&m, x1.get(), y1.get()));
  ref_t p2(&m, MDD_ops<int64_t>::meet(&m, x3.get(), y3.get()));

  ref_t p(&m, MDD_ops<int64_t>::join(&m, p1.get(), p2.get()));

  std::cout << "\"(0<= v_3  <2 and 0 <= v_5 < 2) or (4 <= v_3 < 6 and 3 <= v_5 < 5)\"\n";
  vec<size_t> pvars { 2, 3, 6 };
  ref_t w(&m, MDD_ops<int64_t>::forget(&m, p.get(), pvars.begin(), pvars.end()));
  std::cout << "After forgetting v_2, v_3, and v_6\n";
  v(w.get());
  std::cout << "== End test_project\n";      
}

struct CheckInt64 {
  bool operator()(int64_t k) const { return true; }
};

typedef mdd_ext_cache<int64_t, CheckInt64, mdd_node<int64_t>* > unary_ecache_t;

/*
struct int64_interval {
  static int64_t pred(int64_t k) { return k-1; }
  static int64_t succ(int64_t k) { return k+1; }
  bool lb_is_finite(void) const { return true; }
  bool ub_is_finite(void) const { return true; }
  int64_t lb(void) const { return _lb; }
  int64_t ub(void) const { return _ub; }
  size_t hash(size_t h) const {
    h = (h<<5) + h + _lb;
    h = (h<<5) + h + _ub;
    return h;
  }
  bool operator==(const int64_interval& o) const { return _lb == o._lb && _ub == o._ub; }

  int64_t _lb;
  int64_t _ub; 
};
*/

static const unary_ecache_t::op_tag ec_tag = unary_ecache_t::bind_op();
void report_ecache(const char* s, unary_ecache_t& c, mdd_node<int64_t>* m) {
  int64_t* r(c.lookup(ec_tag, m));
  fprintf(stdout, "Checking %s: ", s);
  if(r)
    fprintf(stdout, "hit (%ld)", *r);
  else
    fprintf(stdout, "miss");
  fprintf(stdout, "\n");
}

void test_ecache(mgr_t& m) {
  unary_ecache_t c;

  ref_t x1(make_range(m, 3, 0, 2)); 
  ref_t y1(make_range(m, 5, 0, 2));
  ref_t y2(make_range(m, 5, 3, 5));

  fprintf(stdout, "Inserting x1 -> 8\n");
  c.insert(ec_tag, x1.get(), 8);

  report_ecache("x1", c, x1.get());
  report_ecache("y1", c, y1.get());
  report_ecache("y2", c, y2.get());

  fprintf(stdout, "Inserting y1 -> 3\n");
  c.insert(ec_tag, y1.get(), 3);

  report_ecache("x1", c, x1.get());
  report_ecache("y1", c, y1.get());
  report_ecache("y2", c, y2.get());
}


struct int64_repr {
  int64_repr(void) : x(INT64_MIN) { }
  int64_repr(int64_t _x) : x(_x) { }
  bool operator==(const int64_repr& o) const { return x == o.x; }

  bool is_finite(void) const { return x != INT64_MIN; }
  int64_repr operator+(int64_t r) const { return int64_repr(x+r); }
  int64_repr operator*(int64_t c) const { return int64_repr(x*c); }
  bool operator<(const int64_repr& o) const { return x < o.x; }

  static int64_repr minus_infty(void) { return INT64_MIN; }
  static int64_t succ(int64_t k) { return k+1; }
  static int64_t pred(int64_t k) { return k-1; }
  int64_t value(void) const { return x; }
  
  int64_t x;
};

namespace std {
  template<>
  struct hash<int64_repr> {
    size_t operator()(const int64_repr& r) { return r.x; }
  };
};

typedef num_interval<int64_t, int64_repr> int64_interval;
typedef mdd_transformer<int64_t, mdd_assign_interval<int64_t, int64_interval> > assign_t;

typedef mdd_transformer<int64_t, mdd_linexpr_lb<int64_t, int64_repr> > linex_lb;
typedef mdd_transformer<int64_t, mdd_eval_linexpr<int64_t, int64_repr> > linex_eval;

template<typename MDD_printer>
void test_assign(MDD_printer& vis, mgr_t* mgr) {
  std::cout << "== Begin test_assign\n";        
  ref_t ttt(mgr->mdd_true()); 

  int64_interval a(int64_interval::range(10, 1000));
  int64_interval b(int64_interval::range(80, 90));
  int64_interval c(int64_interval::range(200, 300));

  /*
  ref_t p(mgr, assign_t::apply(mgr, ttt.get(), 3, a));
  ref_t q(mgr, assign_t::apply(mgr, p.get(), 5, b));
  vis(q.get());
  */
  ref_t p1(mgr, assign_t::apply(mgr, ttt.get(), 5, b));
  ref_t p2(mgr, assign_t::apply(mgr, ttt.get(), 5, c));
  ref_t p(mgr, MDD_ops<int64_t>::join(mgr, p1.get(), p2.get()));

  std::cout << "\" 80 <= v_5 <= 90 or 200 <= v_5 <= 300 \"\n";
  vis(p.get());
  std::cout << "\"(assign interval) v_3 = [10,1000]\"\n";
  ref_t q(mgr, assign_t::apply(mgr, p.get(), 3, a));
  vis(q.get());
  
  vec< linterm<int64_t> > xs;
  xs.push(linterm<int64_t> { 3, 3 });
  xs.push(linterm<int64_t> { 8, 5 });
  
  // auto pt = hc::lookup(xs);
  // {
  //   int64_t r = linex_lb::apply(mgr, q.get(), hc::lookup(xs)).x;
  //   fprintf(stdout, "%% lb(3 [v_3] + 8 [v_5])[%p] = %ld\n", pt, r);
  // }

  // xs.push(linterm<int64_t> { 1, 8 }); 
  // {
  //   int64_t s = linex_lb::apply(mgr, q.get(), hc::lookup(xs)).x;
  //   fprintf(stdout, "%% lb(3 [v_3] + 8 [v_5] + 1 [v_8]) = %ld\n", s);
  // }

  // xs.pop();
  auto qt = hc::lookup(xs);
  // {
  //   int64_t t = linex_lb::apply(mgr, q.get(), qt).x;
  //   fprintf(stdout, "%% lb(3 [v_3] + 8 [v_5])[%p] = %ld\n", qt, t);
  // }

  {
    std::cout << "== Begin evaluation linear expression \n";
    vis(q.get());
    auto i = linex_eval::apply(mgr, q.get(), qt);
    std::cout << "\"Eval lin exp 3*v_3 + 8*v_5=" << i << "\n";
    std::cout << "== End evaluation linear expression \n";        
  }

  {
    std::cout << "== Begin apply linear expression \n";
    vis(q.get());    
    ref_t r(mgr, mdd_transformer< int64_t, mdd_assign_linexpr<int64_t, int64_repr> >::apply(mgr, q.get(), 4, qt, 10));
    std::cout << "\"(apply lin expr) v_4 = 3*v_3 + 8*v_5 + 10\"\n";
    vis(r.get());
    std::cout << "== End apply linear expression \n";    
  }

  {
    std::cout << "== Begin apply linear expression \n";
    vis(q.get());    
    std::cout << "\"(apply lin expr) v_4 = 3*v_3 + 8*v_5 + 10\"\n";    
    ref_t r(mgr, mdd_transformer< int64_t, mdd_assign_linexpr<int64_t, int64_repr, num_interval_NE<int64_t, int64_repr>, true> >::apply(mgr, q.get(), 4, qt, 10));
    vis(r.get());
    std::cout << "== End apply linear expression \n";        
  }
      
  {
    std::cout << "== Begin add linear inequality \n";
    vis(q.get());    
    ref_t r(mgr, mdd_transformer< int64_t, mdd_lin_leq<int64_t, int64_repr> >::apply(mgr, q.get(), qt, 2000));
    std::cout << "\" 3*v_3 + 8*v_5 <= 2000 \"\n";
    vis(r.get());
    std::cout << "== End add linear inequality \n";        
  }

  {    
    std::cout << "== Begin multiplication \n";
    vis(q.get());
    vec<unsigned int> vs { 3, 5 };
    hc::hc_list<unsigned int>* xs(hc::lookup(vs));    
    typedef mdd_eval_prod<// type for number
                          int64_t,
			  // number augmented with -oo, +oo
			  int64_repr,
			  // return type
                          num_interval<int64_t, int64_repr>> eval_prod_t;
    auto i = mdd_transformer<int64_t, eval_prod_t>::apply(mgr, q.get(), xs);
    std::cout << "Evaluation of (v_3 * v_5)=" << i << "\n";
    std::cout << "\"v_12 =  (3*(v_3 * v_5)) + 1\"\n";
    ref_t r(mgr, mdd_transformer< int64_t, mdd_assign_prod<int64_t, int64_repr, num_interval<int64_t, int64_repr> > >::apply(mgr, q.get(), 12, xs, 3));
    vis(r.get());
    std::cout << "== End multiplication \n";                
  }

  {
    std::cout << "== Begin division \n";
    vis(q.get());        
    std::cout << "\" v_4 =  (v_5 / v_3)\"\n";
    ref_t r(mgr, mdd_transformer< int64_t,
            mdd_assign_div<int64_t, int64_repr, num_interval<int64_t, int64_repr> > >::
	    apply(mgr, q.get(), 4, 5, 3));
    vis(r.get());
    std::cout << "== End division \n";    
  }

  {
    // Extract bounds from the MDD.
    vec< var_lb<int64_t> > lbs;
    vec< var_ub<int64_t> > ubs;
    mdd_lower_bounds<int64_t>::eval(q.get(), lbs);
    mdd_upper_bounds<int64_t>::eval(q.get(), ubs);

    fprintf(stdout, "Lower bounds:");
    for(var_lb<int64_t> b : lbs)
      fprintf(stdout, " [v_%d >= %ld]", b.var, b.lb);
    fprintf(stdout, "\n");

    fprintf(stdout, "Upper bounds:");
    for(var_ub<int64_t> b : ubs)
      fprintf(stdout, " [v_%d < %ld]", b.var, b.ub);
    fprintf(stdout, "\n");
  }
  std::cout << "== End test_assign\n";          
}

template<typename MDD_printer>
void test_convexify(MDD_printer& vis, mgr_t* mgr) {
  std::cout << "== Begin test_convexify\n";            
  ref_t ttt(mgr->mdd_true()); 

  int64_interval a(int64_interval::range(10, 20));
  int64_interval b(int64_interval::range(80, 100));
  int64_interval c(int64_interval::range(95, 300));

  ref_t p(mgr, assign_t::apply(mgr, assign_t::apply(mgr, ttt.get(), 1, a), 2, b));
  ref_t q(mgr, assign_t::apply(mgr, assign_t::apply(mgr, ttt.get(), 1, b), 2, c));

  ref_t r(mgr, MDD_ops<int64_t>::join(mgr, p.get(), q.get()));

  std::cout << "\"(v_1 = [10,20) and v_2 = [80,100)) or (v_1 = [80,100) and v_2 = [95,300)) \"\n";
  vis(r.get());
  std::cout << "After convexify:\n";
  ref_t s(mgr, convexify<int64_t>(mgr, r.get()));
  vis(s.get());
  std::cout << "== End test_convexify\n";              
}

typedef mdd_transformer<int64_t, mdd_rename<int64_t, int64_repr> > rename_t;

template<typename MDD_printer>
void test_rename(MDD_printer& vis, mgr_t* mgr) {
  std::cout << "== Begin test_rename\n";                
  ref_t ttt(mgr->mdd_true());

  int64_interval a(int64_interval::range(10, 20));
  int64_interval b(int64_interval::range(80, 100));
  int64_interval c(int64_interval::range(95, 300));

  ref_t p(mgr, assign_t::apply(mgr, assign_t::apply(mgr, ttt.get(), 1, a), 2, b));
  ref_t q(mgr, assign_t::apply(mgr, assign_t::apply(mgr, ttt.get(), 1, b), 2, c));

  ref_t r(mgr, MDD_ops<int64_t>::join(mgr, p.get(), q.get()));
  std::cout << "\"(v_1 = [10,20) and v_2 = [80,100)) or (v_1 = [80,100) and v_2 = [95,300)) \"\n";
  vis(r.get());

  vec<var_pair> pi;
  pi.push(var_pair { 1, 3 });
  pi.push(var_pair { 2, 1 });
  pi.push(var_pair { 3, 2 });

  ref_t s(mgr, rename_t::apply(mgr, r.get(), hc::lookup(pi)));
  std::cout << "\"After renaming v_1 <-> v_3, v_2 <-> v_1, and v_3 <-> v_2\"\n";
  vis(s.get());
  std::cout << "== End test_rename\n";                  

}

template<typename MDD_printer>
void test_copy(MDD_printer& v, mgr_t* m) {
  std::cout << "== Begin test_copy\n";                    
  ref_t p0(m->var_lt(10, 100));
  ref_t q0(p0); // make a copy

  std::cout << "mdd1: \" v_10 < 100 \"\n";
  v(p0.get());
  std::cout << "mdd2: \" v_10 < 100 \"\n";  
  v(q0.get());

  int64_interval a(int64_interval::range(80, 90));
  ref_t q1(m, assign_t::apply(m, q0.get(), 5, a));

  std::cout << "mdd1: \" v_10 < 100 \"\n";
  v(p0.get());
  std::cout << "mdd2: \" v_10 < 100 and  v_5 = [80,91)\"\n";  
  v(q1.get());
  std::cout << "== End test_copy\n";                      
}

template<typename MDD_printer>
void test_copy_assign(MDD_printer& v, mgr_t* m) {
  std::cout << "== Begin test_copy_assign\n";                    
  ref_t p0(m->var_lt(10, 100));
  ref_t q0 = p0; 
  std::cout << "mdd1: \" v_10 < 100 \"\n";
  v(p0.get());
  std::cout << "mdd2: \" v_10 < 100 \"\n";  
  v(q0.get());

  int64_interval a(int64_interval::range(80, 90));
  ref_t q1(m, assign_t::apply(m, q0.get(), 5, a));

  std::cout << "mdd1: \" v_10 < 100 \"\n";
  v(p0.get());
  std::cout << "mdd2: \" v_10 < 100 and  v_5 = [80,91)\"\n";  
  v(q1.get());
  std::cout << "== End test_copy\n";                      
}

int main(int argc, char** argv) {
  mgr_t m;
  
  ShowNop<int64_t> nop;
  MDD_printer_t vis(std::cout, nop);
  test_meet(vis, m);
  test_widen(vis, m);
  test_project(vis, m);
  //test_ecache(m);
  test_assign(vis, &m);
  test_convexify(vis, &m);
  test_rename(vis, &m);
  test_copy(vis, &m);
  test_copy_assign(vis, &m);  
  return 0;
}
#endif
