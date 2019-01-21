#pragma once

/*******************************************************************************
 * Implementation of Octagons based on split normal form.
 * Author: Zequn Ma.
 *******************************************************************************/

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/graphs/adapt_sgraph.hpp>
#include <crab/domains/graphs/sparse_graph.hpp>
#include <crab/domains/graphs/ht_graph.hpp>
#include <crab/domains/graphs/pt_graph.hpp>
#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>

#include <type_traits>

#include <boost/optional.hpp>
#include <boost/unordered_set.hpp>
#include <boost/container/flat_map.hpp>

#define CLOSE_BOUNDS_INLINE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {
  /* 
     It was originally in graphs/graph_ops.hpp
     This is similar to SubGraph but for DBMs adapted for Octagos.
   */
  
  // View of a graph, omitting a given vertex
  template<class G>
  class SplitGraph {
  public:
    typedef typename G::vert_id vert_id;
    typedef typename G::Wt Wt;

    typedef typename G::pred_range g_pred_range;
    typedef typename G::succ_range g_succ_range;

    typedef typename G::e_pred_range g_e_pred_range;
    typedef typename G::e_succ_range g_e_succ_range;

    typedef typename G::mut_val_ref_t mut_val_ref_t;

    SplitGraph(G& _g)
      : g(_g)
    { }

     
    bool elem(vert_id x, vert_id y) const {
      return ((x / 2 != y / 2) && g.elem(x, y));
    }

    bool lookup(vert_id x, vert_id y, mut_val_ref_t* w) const {
      return ((x / 2 != y / 2) && g.lookup(x, y, w));
    }

    Wt edge_val(vert_id x, vert_id y) const {
      return g.edge_val(x, y);  
    }

    // Precondition: elem(x, y) is true.
    Wt operator()(vert_id x, vert_id y) const {
      return g(x, y);
    }

    void clear_edges(void) { g.clear_edges(); }

    void clear(void)
    {
      assert(0 && "SplitGraph::clear not implemented.");       
    }

    // Number of allocated vertices
    int size(void) const {
      return g.size();
    }

    // Assumption: (x, y) not in mtx
    void add_edge(vert_id x, Wt wt, vert_id y)
    {
//      assert(x != v_ex && y != v_ex);
      g.add_edge(x, wt, y);
    }

    void set_edge(vert_id s, Wt w, vert_id d)
    {
//      assert(s != v_ex && d != v_ex);
      g.set_edge(s, w, d);
    }

    template<class Op>
    void update_edge(vert_id s, Wt w, vert_id d, Op& op)
    {
//      assert(s != v_ex && d != v_ex);
      g.update_edge(s, w, d, op);
    }

    class vert_iterator {
    public:
      vert_iterator(const typename G::vert_iterator& _iG)
        : iG(_iG)
      { }

      // Skipping of v_ex is done entirely by !=.
      // So we _MUST_ test it != verts.end() before dereferencing.
      vert_id operator*(void) { return *iG; }
      vert_iterator operator++(void) { ++iG; return *this; } 
      bool operator!=(const vert_iterator& o) {
        return iG != o.iG;
      }
      
      typename G::vert_iterator iG; 
    };
    class vert_range {
    public:
      vert_range(const typename G::vert_range& _rG)
        : rG(_rG)
      { }

      vert_iterator begin(void) const { return vert_iterator(rG.begin()); }
      vert_iterator end(void) const { return vert_iterator(rG.end()); }

      typename G::vert_range rG;
    };
    vert_range verts(void) const { return vert_range(g.verts()); }

    template<class It>
    class adj_iterator {
    public:
      adj_iterator(const It& _iG, vert_id _v_ex)
        : iG(_iG), v_ex(_v_ex)
      { }
      vert_id operator*(void) const { return *iG; }
      adj_iterator& operator++(void) { ++iG; return *this; }
      bool operator!=(const adj_iterator& o)
      {
        if(iG != o.iG && (*iG) == v_ex)
          ++iG;
        return iG != o.iG;
      }

      It iG;
      vert_id v_ex;
    };

    template<class It>
    class e_adj_iterator {
    public:
      typedef typename It::edge_ref edge_ref;

      e_adj_iterator(const It& _iG, vert_id _v_ex)
        : iG(_iG), v_ex(_v_ex)
      { }
      edge_ref operator*(void) const { return *iG; }
      e_adj_iterator& operator++(void) { ++iG; return *this; }
      bool operator!=(const e_adj_iterator& o)
      {
        if(iG != o.iG && (*iG).vert == v_ex)
          ++iG;
        return iG != o.iG;
      }

      It iG;
      vert_id v_ex;
    };

    template<class R, class It>
    class adj_list {
    public: 
      typedef typename R::iterator g_iter;
      typedef It iterator;

      adj_list(const R& _rG, vert_id _v_ex)
        : rG(_rG), v_ex(_v_ex)
      { }
      iterator begin() const { return iterator(rG.begin(), v_ex); }
      iterator end() const { return iterator(rG.end(), v_ex); }
      
    protected:
      R rG;
      vert_id v_ex;
    };
    typedef adj_list<g_pred_range,
              adj_iterator<typename g_pred_range::iterator> > pred_range;
    typedef adj_list<g_succ_range,
              adj_iterator<typename g_succ_range::iterator> > succ_range;

    typedef adj_list<g_e_pred_range,
              e_adj_iterator<typename g_e_pred_range::iterator> > e_pred_range;
    typedef adj_list<g_e_succ_range,
              e_adj_iterator<typename g_e_succ_range::iterator> > e_succ_range;

    succ_range succs(vert_id v) {
      // assert(v != v_ex);
      vert_id v_opp;
      if (v % 2 == 0) v_opp = v + 1;
      else v_opp = v - 1;
      return succ_range(g.succs(v), v_opp);
    }
    pred_range preds(vert_id v) {
      // assert(v != v_ex);
      vert_id v_opp;
      if (v % 2 == 0) v_opp = v + 1;
      else v_opp = v - 1;
      return pred_range(g.preds(v), v_opp);
    }
    e_succ_range e_succs(vert_id v) {
      vert_id v_opp;
      if (v % 2 == 0) v_opp = v + 1;
      else v_opp = v - 1;
      return e_succ_range(g.e_succs(v), v_opp);
    }
    e_pred_range e_preds(vert_id v) {
      vert_id v_opp;
      if (v % 2 == 0) v_opp = v + 1;
      else v_opp = v - 1;
      return e_pred_range(g.e_preds(v), v_opp);
    }

    G& g;
  };
    
  namespace domains {

    namespace split_octagons_impl {

      template<typename Number, typename Wt>
      class NtoV {
      public:
	static Wt ntov(const Number& n) {
	  return (Wt) n;
	}
      };

      enum GraphRep {
	ss = 1,
	adapt_ss = 2,
	pt = 3,
	ht = 4
      };

      // enum Sign {
      // 	pos = 0,
      // 	neg = 1,
      // };

      template<typename Number, GraphRep Graph = GraphRep::adapt_ss>
      class DefaultParams {
      public:
	enum { chrome_dijkstra = 1 };
	enum { widen_restabilize = 1 };
	enum { special_assign = 1 };

	//typedef Number Wt;
	typedef long Wt;

	typedef typename std::conditional< 
	  (Graph == ss), SparseWtGraph<Wt>,
	  typename std::conditional< 
	    (Graph == adapt_ss), AdaptGraph<Wt>,
	    typename std::conditional< 
	      (Graph == pt), PtGraph<Wt>, HtGraph<Wt> 
	      >::type 
	    >::type 
	  >::type graph_t;
      };

      template<typename Number, GraphRep Graph = GraphRep::adapt_ss>
      class SimpleParams {
      public:
	enum { chrome_dijkstra = 0 };
	enum { widen_restabilize = 0 };
	enum { special_assign = 0 };

	typedef long Wt;

	typedef typename std::conditional< 
	  (Graph == ss), SparseWtGraph<Wt>,
	  typename std::conditional< 
	    (Graph == adapt_ss), AdaptGraph<Wt>,
	    typename std::conditional< 
	      (Graph == pt), PtGraph<Wt>, HtGraph<Wt> 
	      >::type 
	    >::type 
	  >::type graph_t;
      };
    }; // end namespace split_octagons_impl

    template<class Number, class VariableName,
	     class Params = split_octagons_impl::DefaultParams<Number>>
    class split_oct_domain :
      public abstract_domain<split_oct_domain<Number, VariableName, Params>> {
      typedef split_oct_domain<Number, VariableName, Params> split_oct_domain_t;
      typedef abstract_domain<split_oct_domain_t> abstract_domain_t;
			
    public:
      using typename abstract_domain_t::linear_expression_t;
      using typename abstract_domain_t::linear_constraint_t;
      using typename abstract_domain_t::linear_constraint_system_t;
      using typename abstract_domain_t::disjunctive_linear_constraint_system_t;      
      using typename abstract_domain_t::variable_t;
      using typename abstract_domain_t::variable_vector_t;
      using typename abstract_domain_t::pointer_constraint_t;
      typedef Number number_t;
      typedef VariableName varname_t;
      typedef typename linear_constraint_t::kind_t constraint_kind_t;
      typedef interval<number_t>  interval_t;

    private:
      typedef bound<number_t> bound_t;
      typedef typename Params::Wt Wt;
      typedef typename Params::graph_t graph_t;
      typedef split_octagons_impl::NtoV<number_t, Wt> ntov;
      typedef typename graph_t::vert_id vert_id;
      // (variable: (vert_id of pos, vert_id of neg))
      typedef boost::container::flat_map<variable_t, std::pair<vert_id, vert_id>> vert_map_t;
      typedef typename vert_map_t::value_type vmap_elt_t;
      typedef std::vector<boost::optional<variable_t>> rev_map_t;
      typedef GraphOps<graph_t> GrOps;
      typedef GraphPerm<graph_t> GrPerm;
      typedef typename GrOps::edge_vector edge_vector;
      typedef std::pair<std::pair<variable_t, variable_t>, Wt> diffcst_t;
      typedef boost::unordered_set<vert_id> vert_set_t;

      vert_map_t vert_map;
      rev_map_t rev_map;
      graph_t graph;
      std::vector<Wt> potential;
      vert_set_t unstable;
      bool _is_bottom;

      class Wt_max {
      public:
	Wt_max() { }
	Wt apply(const Wt& x, const Wt& y) { return std::max(x, y); }
	bool default_is_absorbing() { return true; }
      };
      class Wt_min {
      public:
	Wt_min() { }
	Wt apply(const Wt& x, const Wt& y) { return std::min(x, y); }
	bool default_is_absorbing() { return false; }
      };
      
      // construct difference constraint (i, (j, k)) from a linear constraint
      boost::optional<std::pair<vert_id, std::pair<vert_id, Wt>>>
      diffcst_of_leq(linear_constraint_t cst) {
	assert (cst.size() > 0);
	assert(cst.is_inequality());

	std::vector<std::pair<vert_id,std::pair<vert_id, Wt> > > diffcsts;
	Wt weight = Wt(cst.constant());

	typename linear_expression_t::iterator it1 = cst.begin();
	typename linear_expression_t::iterator it2 = ++cst.begin();
	vert_id i, j;

	if (cst.size() == 1 && it1->first == 1) {
	  i = get_vert(it1->second.name());
	  j = i + 1;
	  weight = 2 * weight;
	} else if (cst.size() == 1 && it1->first == -1) {
	  i = get_vert(it1->second.name()) + 1;
	  j = i - 1;
	  weight = 2 * weight;
	} else if (cst.size() == 2 && it1->first == 1 && it2->first == -1) {
	  i = get_vert(it1->second.name());
	  j = get_vert(it2->second.name());
	} else if (cst.size() == 2 && it1->first == -1 && it2->first == 1) {
	  i = get_vert(it2->second.name());
	  j = get_vert(it1->second.name());
	} else if (cst.size() == 2 && it1->first == 1 && it2->first == 1) {
	  i = get_vert(it1->second.name());
	  j = get_vert(it2->second.name()) + 1;
	} else if (cst.size() == 2 && it1->first == -1 && it2->first == -1) {
	  i = get_vert(it2->second.name()) + 1;
	  j = get_vert(it1->second.name());
	} else {
	  return boost::none;
	}

	if (abs(i) > abs(j)) {
	  auto tmp = i;
	  i = j;
	  j = tmp;
	}

	return std::make_pair(j, std::make_pair(i, Wt(cst.constant())));
      }

      // get vert_id of both v+ and v- after creating verices in graph
      // return vert_id of positive node, negative node is positive + 1
      vert_id get_vert(variable_t v) {
	auto it = vert_map.find(v);
	if (it != vert_map.end()) return(*it).second.first;

	vert_id vert_pos(graph.new_vertex());
	vert_id vert_neg(graph.new_vertex());
	if (vert_pos > vert_neg) {
	  vert_id tmp = vert_pos;
	  vert_pos = vert_neg;
	  vert_neg = tmp;
	}
	vert_map.insert(vmap_elt_t(v, std::make_pair(vert_pos, vert_neg)));
	assert(vert_pos <= rev_map.size());
	assert(vert_neg <= rev_map.size()+1);

	if (vert_pos < rev_map.size()) {
	  potential[vert_pos] = Wt(0);
	  rev_map[vert_pos] = v;
	} else {
	  potential.push_back(Wt(0));
	  rev_map.push_back(v);
	}
		    	
	if (vert_neg < rev_map.size()) {
	  potential[vert_neg] = Wt(0);
	  rev_map[vert_neg] = v;
	} else {
	  potential.push_back(Wt(0));
	  rev_map.push_back(v);
	}
                        
	//vert_map.insert(vmap_elt_t(v, std::make_pair(vert_pos, vert_neg)));
                            
	return vert_pos;
      }

      vert_id get_vert(graph_t& g, vert_map_t& vmap, rev_map_t& rmap,
		       std::vector<Wt>& pot, variable_t v) {
	auto it = vmap.find(v);
	if (it != vmap.end()) return (*it).second.first;

	vert_id vert_pos(g.new_vertex());
	vert_id vert_neg(g.new_vertex());
	if (vert_pos > vert_neg) {
	  //CRAB_LOG ("octagon-split", crab::outs() << "\tfixing pos and neg" << "\n");
	  vert_id tmp = vert_pos;
	  vert_pos = vert_neg;
	  vert_neg = tmp;
	}
	vmap.insert(vmap_elt_t(v, std::make_pair(vert_pos, vert_neg)));
	assert(vert_pos <= rev_map.size());
	assert(vert_neg <= rev_map.size()+1);

	if (vert_pos < rmap.size()) {
	  pot[vert_pos] = Wt(0);
	  rmap[vert_pos] = v;
	} else {
	  pot.push_back(Wt(0));
	  rmap.push_back(v);
	}

	if (vert_neg < rmap.size()) {
	  pot[vert_neg] = Wt(0);
	  rmap[vert_neg] = v;
	} else {
	  pot.push_back(Wt(0));
	  rmap.push_back(v);
	}

	//vmap.insert(vmap_elt_t(v, std::make_pair(vert_pos, vert_neg)));

	return vert_pos;
      }

      class vert_set_wrap_t {
      public:
	const vert_set_t& vs;
	vert_set_wrap_t(const vert_set_t& _vs): vs(_vs) { }
	bool operator[](vert_id v) const { return vs.find(v) != vs.end(); }
      };

      template<class G, class P>
      inline bool check_potential(G& g, P& p) {
        #ifdef CHECK_POTENTIAL
	for(vert_id v : g.verts()) {
	  for(vert_id d : g.succs(v)) {
	    if(p[v] + g.edge_val(v, d) - p[d] < Wt(0)) {
	      assert(0 && "Invalid potential.");
	      return false;
	    }
	  }
	}
        #endif
	return true;
      }


      // get potential of a variable
      Wt pot_value(variable_t v) {
	auto it = vert_map.find(v);
	if (it != vert_map.end()) return potential[(*it).second.first];
	return ((Wt)0);
      }
      Wt pot_value(variable_t v, std::vector<Wt>& potential) {
	auto it = vert_map.find(v);
	if (it != vert_map.end()) return potential[(*it).second.first];
	return ((Wt)0);
      }

      Wt eval_expression(linear_expression_t e) {
	Wt v(ntov::ntov(e.constant())); 
	for(auto p : e) {
	  v += (pot_value(p.second))*(ntov::ntov(p.first));
	}
	return v;
      }

      interval_t eval_interval(linear_expression_t e) {
	interval_t r = e.constant();
	for (auto p : e) r += p.first * operator[](p.second);
	return r;
      }

      // Turn an assignment into a set of difference constraints.
      void diffcsts_of_assign(variable_t x, linear_expression_t exp,
			      std::vector<std::pair<variable_t, Wt>>& lb,
			      std::vector<std::pair<variable_t, Wt>>& ub) {
	{
	  boost::optional<variable_t> unbounded_ubvar;
	  Wt exp_ub(ntov::ntov(exp.constant()));
	  std::vector<std::pair<variable_t, Wt>> ub_terms;
	  for (auto p : exp) {
	    Wt coeff(ntov::ntov(p.first));
	    if (p.first < Wt(0)) {
	      bound_t y_lb = operator[](p.second).lb();
	      if (y_lb.is_infinite()) goto assign_ub_finish;
	      exp_ub += ntov::ntov(*(y_lb.number()))*coeff;
	    } else {
	      variable_t y(p.second);
	      bound_t y_ub = operator[](y).ub();
	      if (y_ub.is_infinite()) {
		if (unbounded_ubvar || coeff != Wt(1))
		  goto assign_ub_finish;
		unbounded_ubvar = y;
	      } else {
		Wt ymax(ntov::ntov(*(y_ub.number())));
		exp_ub += ymax*coeff;
		ub_terms.push_back(std::make_pair(y, ymax));
	      }
	    }
	  }
	  if (unbounded_ubvar) {
	    ub.push_back(std::make_pair(*unbounded_ubvar, exp_ub));
	  } else {
	    for (auto p : ub_terms)
	      ub.push_back(std::make_pair(p.first, exp_ub - p.second));
	  }
	}

      assign_ub_finish:
	{
	  boost::optional<variable_t> unbounded_lbvar;
	  Wt exp_lb(ntov::ntov(exp.constant()));
	  std::vector< std::pair<variable_t, Wt> > lb_terms;
	  for(auto p : exp) {
	    Wt coeff(ntov::ntov(p.first));
	    if(p.first < Wt(0)) {
	      // Again, can't do anything with negative coefficients.
	      bound_t y_ub = operator[](p.second).ub();
	      if(y_ub.is_infinite())
		goto assign_lb_finish;
	      exp_lb += (ntov::ntov(*(y_ub.number())))*coeff;
	    } else {
	      variable_t y(p.second);
	      bound_t y_lb = operator[](y).lb(); 
	      if(y_lb.is_infinite()) {
		if(unbounded_lbvar || coeff != Wt(1))
		  goto assign_lb_finish;
		unbounded_lbvar = y;
	      } else {
		Wt ymin(ntov::ntov(*(y_lb.number())));
		exp_lb += ymin*coeff;
		lb_terms.push_back(std::make_pair(y, ymin));
	      }
	    }
	  }

	  if(unbounded_lbvar) {
	    lb.push_back(std::make_pair(*unbounded_lbvar, exp_lb));
	  } else {
	    for(auto p : lb_terms) {
	      lb.push_back(std::make_pair(p.first, exp_lb - p.second));
	    }
	  }
	}

      assign_lb_finish:
	return;
      }

      void oct_diffcsts_of_lin_leq(const linear_expression_t& exp,
				   std::vector<diffcst_t>& csts, bool add) {
	Wt unbounded_lbcoeff;
	Wt unbounded_ubcoeff;
	boost::optional<variable_t> unbounded_lbvar;
	boost::optional<variable_t> unbounded_ubvar;
	Wt exp_ub = - (ntov::ntov(exp.constant()));
	std::vector< std::pair< std::pair<Wt, variable_t>, Wt> > pos_terms;
	std::vector< std::pair< std::pair<Wt, variable_t>, Wt> > neg_terms;
                            
	for (auto p : exp) {
	  Wt coeff(ntov::ntov(p.first));
	  if(coeff > Wt(0)) {
	    variable_t y(p.second);
	    bound_t y_lb = operator[](y).lb();
	    if(y_lb.is_infinite()) {
	      if(unbounded_lbvar) return;
	      unbounded_lbvar = y;
	      unbounded_lbcoeff = coeff;
	    } else {
	      Wt ymin(ntov::ntov(*(y_lb.number())));
	      // Coeff is negative, so it's still add
	      exp_ub -= ymin*coeff;
	      pos_terms.push_back(std::make_pair(std::make_pair(coeff, y), ymin));
	    }
	  } else {
	    variable_t y(p.second);
	    bound_t y_ub = operator[](y).ub(); 
	    if(y_ub.is_infinite()) {
	      if(unbounded_ubvar) return;
	      unbounded_ubvar = y;
	      unbounded_ubcoeff = -(ntov::ntov(coeff));
	    } else {
	      Wt ymax(ntov::ntov(*(y_ub.number())));
	      exp_ub -= ymax*coeff;
	      neg_terms.push_back(std::make_pair(std::make_pair(-coeff, y), ymax));
	    }
	  }
	}
                            
	if (neg_terms.size() == 0 && add){
	  for(auto pu1 : pos_terms)
	    for(auto pu2 : pos_terms){
	      if (pu1.first.second == pu2.first.second) continue;
	      csts.push_back(std::make_pair(std::make_pair(pu1.first.second, pu2.first.second),
					    exp_ub + pu1.second + pu2.second));
	    }
	}
                            
	if (pos_terms.size() == 0 && !add){
	  for(auto pu1 : neg_terms)
	    for(auto pu2 : neg_terms){
	      if (pu1.first.second == pu2.first.second) continue;
	      csts.push_back(std::make_pair(std::make_pair(pu1.first.second, pu2.first.second),
					    exp_ub - pu1.second - pu2.second));
	    }
	}
      }
                        
			
      void diffcsts_of_lin_leq(const linear_expression_t& exp, std::vector<diffcst_t>& csts,
			       std::vector<std::pair<variable_t, Wt> >& lbs,
			       std::vector<std::pair<variable_t, Wt> >& ubs) {

	Wt unbounded_lbcoeff;
	Wt unbounded_ubcoeff;
	boost::optional<variable_t> unbounded_lbvar;
	boost::optional<variable_t> unbounded_ubvar;
	Wt exp_ub = - (ntov::ntov(exp.constant()));
	std::vector< std::pair< std::pair<Wt, variable_t>, Wt> > pos_terms;
	std::vector< std::pair< std::pair<Wt, variable_t>, Wt> > neg_terms;

	for(auto p : exp) {
	  Wt coeff(ntov::ntov(p.first));
	  if(coeff > Wt(0)) {
	    variable_t y(p.second);
	    bound_t y_lb = operator[](y).lb();
	    if(y_lb.is_infinite()) {
	      if(unbounded_lbvar) goto diffcst_finish;
	      unbounded_lbvar = y;
	      unbounded_lbcoeff = coeff;
	    } else {
	      Wt ymin(ntov::ntov(*(y_lb.number())));
	      // Coeff is negative, so it's still add
	      exp_ub -= ymin*coeff;
	      pos_terms.push_back(std::make_pair(std::make_pair(coeff, y), ymin));
	    }
	  } else {
	    variable_t y(p.second);
	    bound_t y_ub = operator[](y).ub(); 
	    if(y_ub.is_infinite()) {
	      if(unbounded_ubvar) goto diffcst_finish;
	      unbounded_ubvar = y;
	      unbounded_ubcoeff = -(ntov::ntov(coeff));
	    } else {
	      Wt ymax(ntov::ntov(*(y_ub.number())));
	      exp_ub -= ymax*coeff;
	      neg_terms.push_back(std::make_pair(std::make_pair(-coeff, y), ymax));
	    }
	  }
	}

	if(unbounded_lbvar) {
	  variable_t x(*unbounded_lbvar);
	  if(unbounded_ubvar) {
	    if(unbounded_lbcoeff != Wt(1) || unbounded_ubcoeff != Wt(1))
	      goto diffcst_finish;
	    variable_t y(*unbounded_ubvar);
	    csts.push_back(std::make_pair(std::make_pair(x, y), exp_ub));
	  } else {
	    if(unbounded_lbcoeff == Wt(1)) {
	      for(auto p : neg_terms)
		csts.push_back(std::make_pair(std::make_pair(x, p.first.second),
					      exp_ub - p.second));
	    }
	    // Add bounds for x
	    ubs.push_back(std::make_pair(x, exp_ub/unbounded_lbcoeff));
	  }
	} else {
	  if(unbounded_ubvar) {
	    variable_t y(*unbounded_ubvar);
	    if(unbounded_ubcoeff == Wt(1)) {
	      for(auto p : pos_terms)
		csts.push_back(std::make_pair(std::make_pair(p.first.second, y),
					      exp_ub + p.second));
	    }
	    // Bounds for y
	    lbs.push_back(std::make_pair(y, -exp_ub/unbounded_ubcoeff));
	  } else {
	    for(auto pl : neg_terms)
	      for(auto pu : pos_terms)
		csts.push_back(std::make_pair(std::make_pair(pu.first.second, pl.first.second),
					      exp_ub - pl.second + pu.second));
	    for(auto pl : neg_terms)
	      lbs.push_back(std::make_pair(pl.first.second,
					   -exp_ub/pl.first.first + pl.second));
	    for(auto pu : pos_terms)
	      ubs.push_back(std::make_pair(pu.first.second,
					   exp_ub/pu.first.first + pu.second));
	  }
	}
      diffcst_finish:
	return;
      }

      bool add_linear_leq(const linear_expression_t& exp) {
	normalize();
	CRAB_LOG("octagon-split", linear_expression_t exp_tmp(exp);
		 crab::outs() << "Adding: "<< exp_tmp << "<= 0" <<"\n");
	CRAB_LOG("octagon-add", linear_expression_t exp_tmp(exp);
		 crab::outs() << "Adding: "<< exp_tmp << "<= 0" <<" to:\n"<<*this<<"\n");

	std::vector< std::pair<variable_t, Wt> > lbs;
	std::vector< std::pair<variable_t, Wt> > ubs;
	std::vector<diffcst_t> csts;
	std::vector<diffcst_t> csts_add;
	std::vector<diffcst_t> csts_min;
	diffcsts_of_lin_leq(exp, csts, lbs, ubs);
	oct_diffcsts_of_lin_leq(exp, csts_add, true);
	oct_diffcsts_of_lin_leq(exp, csts_min, false);

	// assert(check_potential(graph, potential));

	Wt_min min_op;
	typename graph_t::mut_val_ref_t w;
	for(auto p : lbs) {
	  CRAB_LOG("octagon-split", crab::outs() << p.first<< ">="<< p.second <<"\n");
	  variable_t x(p.first);
	  // CRAB_LOG("octagon-split", crab::outs() << "declare variable " << p.first << "\n");
	  vert_id v = get_vert(p.first);
	  // CRAB_LOG("octagon-split", crab::outs() << "got vertex for " << p.first << " at " << v << "\n");
	  if(graph.lookup(v, v+1, &w) && w <= -2*p.second) continue;
	  graph.set_edge(v, -2*p.second, v+1);
	  //CRAB_LOG("octagon-split", crab::outs() << "added edge for lower bound of " << p.first << "\n" << *this << "\n");
	  if(!repair_potential(v, v+1)) {
	    set_to_bottom();
	    return false;
	  }
	  //CRAB_LOG("octagon-split", crab::outs() << "added edge for lower bound of " << p.first << "\n" << *this << "\n");

	  assert(check_potential(graph, potential));
	  // Compute other updated bounds
#ifdef CLOSE_BOUNDS_INLINE
	  //for(auto e : graph.e_preds(v)) {
	  //	if(e.vert == v+1) continue;
	  //	graph.update_edge(e.vert, e.val - 2 * p.second, v+1, min_op);
	  //	if(!repair_potential(e.vert, v+1)) {
	  //		set_to_bottom();
	  //		return false;
	  //	}
	  //	assert(check_potential(graph, potential));
	  //}
	  for(auto e : graph.e_preds(v)) {
	    if (e.vert % 2 == 0) continue;
	    for (auto f : graph.e_succs(v+1)) {
	      if (f.vert + 1 != e.vert) continue;
	      graph.update_edge(f.vert, f.val + 2*p.second + e.val, e.vert, min_op);
	      if(!repair_potential(f.vert, e.vert)) {
		set_to_bottom();
		return false;
	      }
	      assert(check_potential(graph, potential));
	    }
	  }

#endif
	}

	for(auto p : ubs) {
	  CRAB_LOG("octagon-split", crab::outs() << p.first<< "<="<< p.second <<"\n");
	  CRAB_LOG("octagon-add", crab::outs() << p.first<< "<="<< p.second <<"\n");
	  variable_t x(p.first);
	  vert_id v = get_vert(p.first);
	  if(graph.lookup(v+1, v, &w) && w <= 2 * p.second) {
	    //CRAB_LOG("octagon-split", crab::outs() << "edge exists for upper bound of " << p.first << "with wt of " << w << "\n");
	    continue;
	  }

	  graph.set_edge(v+1, 2 * p.second, v);
	  graph.lookup(v+1, v, &w);
	  //CRAB_LOG("octagon-split", crab::outs() << "added edge for upper bound of " << p.first << "with wt of " << w << "\n");
	  if(!repair_potential(v+1, v)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));
#ifdef CLOSE_BOUNDS_INLINE
	  //for(auto e : graph.e_succs(v)) {
	  //	if(e.vert == v+1) continue;
	  //	graph.update_edge(v+1, e.val + 2 * p.second, e.vert, min_op);
	  //	if(!repair_potential(v, e.vert)) {
	  //		set_to_bottom();
	  //		return false;
	  //	}
	  //	assert(check_potential(graph, potential));
	  //}
	  for(auto e : graph.e_succs(v)) {
	    if (e.vert % 2 != 0) continue;
	    for (auto f : graph.e_preds(v+1)) {
	      if (f.vert - 1 != e.vert) continue;
	      graph.update_edge(f.vert, f.val + 2*p.second + e.val, e.vert, min_op);
	      if(!repair_potential(f.vert, e.vert)) {
		set_to_bottom();
		return false;
	      }
	      assert(check_potential(graph, potential));
	    }
	  }
#endif
	}

	for(auto diff : csts) {
	  CRAB_LOG("octagon-split", crab::outs() << diff.first.first<< "-"<< diff.first.second<< "<="
		   << diff.second <<"\n");
	  CRAB_LOG("octagon-add", crab::outs() << diff.first.first<< "-"<< diff.first.second<< "<="
		   << diff.second <<"\n");

	  vert_id src = get_vert(diff.first.second);
	  vert_id dest = get_vert(diff.first.first);
	  graph.update_edge(src, diff.second, dest, min_op);
	  graph.update_edge(dest+1, diff.second, src+1, min_op);
	  if(!repair_potential(src, dest)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));

	  if(!repair_potential(dest+1, src+1)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));
	  close_over_edge(src, dest);
	  close_over_edge(dest+1, src+1);
	  //assert(check_potential(graph, potential));
	}
				
	for(auto diff : csts_add) {
	  CRAB_LOG("octagon-split", crab::outs() << diff.first.first<< "+"<< diff.first.second<< "<="
		   << diff.second <<"\n");
	  CRAB_LOG("octagon-add", crab::outs() << diff.first.first<< "+"<< diff.first.second<< "<="
		   << diff.second <<"\n");

	  vert_id src = get_vert(diff.first.second);
	  vert_id dest = get_vert(diff.first.first);
	  graph.update_edge(src, diff.second, dest+1, min_op);
	  graph.update_edge(dest, diff.second, src+1, min_op);
	  if(!repair_potential(src, dest+1)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));

	  if(!repair_potential(dest, src+1)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));
	  close_over_edge(src, dest+1);
	  close_over_edge(dest, src+1);
	  //assert(check_potential(graph, potential));
	}
				
	for(auto diff : csts_min) {
	  CRAB_LOG("octagon-split", crab::outs() << "-"<<diff.first.first<< "-"<< diff.first.second<< "<="
		   << diff.second <<"\n");
	  CRAB_LOG("octagon-add", crab::outs() <<"-"<< diff.first.first<< "-"<< diff.first.second<< "<="
		   << diff.second <<"\n");

	  vert_id src = get_vert(diff.first.second);
	  vert_id dest = get_vert(diff.first.first);
	  graph.update_edge(src+1, diff.second, dest, min_op);
	  graph.update_edge(dest+1, diff.second, src, min_op);
	  if(!repair_potential(src+1, dest)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));

	  if(!repair_potential(dest+1, src)) {
	    set_to_bottom();
	    return false;
	  }
	  assert(check_potential(graph, potential));
	  close_over_edge(src+1, dest);
	  close_over_edge(dest+1, src);
	  //assert(check_potential(graph, potential));
	}
	// Collect bounds
	// GKG: Now done in close_over_edge

	CRAB_LOG("octagon-add", crab::outs() << "after adding: " << *this << "\n" << "With graph" << graph << "\n");
#ifndef CLOSE_BOUNDS_INLINE
	edge_vector delta;
	for (vert_id v : gragh.verts()) {
	  if (v % 2 != 0) continue;
	  GrOps::close_after_assign(graph, potential, v, delta);
	  GrOps::update_delta(graph, delta);
	}
#endif
        /* */

	assert(check_potential(graph, potential));
	// CRAB_WARN("SplitDBM::add_linear_leq not yet implemented.");
	CRAB_LOG("octagon-split", crab::outs() << "after adding: " << *this << "\n" << "With graph" << graph << "\n");
	CRAB_LOG("octagon-add", crab::outs() << "after adding: " << *this << "\n" << "With graph" << graph << "\n");
	return true;  
      }

      void add_univar_disequation(variable_t x, number_t n) {
	interval_t i = get_interval(x);
	interval_t new_i = linear_interval_solver_impl::trim_interval<interval_t>(i, interval_t(n));
	CRAB_LOG("octagon-split", crab::outs() << "Adding disequation: " << x << "!=" << n << "\n" << new_i << "\n");
	if (new_i.is_bottom()) set_to_bottom();
	else if (!new_i.is_top() && (new_i <= i)) {
	  vert_id v = get_vert(x);
	  typename graph_t::mut_val_ref_t w;
	  if (new_i.lb().is_finite()) {
	    Wt lb_val = ntov::ntov(-(*(new_i.lb().number())));
	    if (graph.lookup(v, v+1, &w) && 2*lb_val < w) {
	      graph.set_edge(v, 2*lb_val, v+1);
	      if (!repair_potential(v, v+1)) {
		set_to_bottom();
		return;
	      }
	      assert(check_potential(graph, potential));
	    }
	  }

	  if(new_i.ub().is_finite()) {	    
	    // strengthen ub
	    Wt ub_val = ntov::ntov(*(new_i.ub().number()));
	    if(graph.lookup(v+1, v, &w) && (2*ub_val < w)) {
	      graph.set_edge(v+1, 2*ub_val, v);
	      if(!repair_potential(v+1, v)) {
		set_to_bottom();
		return;
	      }
	      assert(check_potential(graph, potential));
	    }
	  }

	}
      }

      interval_t compute_residual(linear_expression_t e, variable_t pivot) {
	interval_t residual(-e.constant());
	for (typename linear_expression_t::iterator it = e.begin(); it != e.end(); ++it) {
	  variable_t v = it->second;
	  if (v.index() != pivot.index()) {
	    residual = residual - (interval_t(it->first) * this->operator[](v));
	  }
	}
	return residual;
      }
			
      void add_disequation(linear_expression_t e) {
	for (typename linear_expression_t::iterator it = e.begin(); it != e.end(); ++it) {
	  variable_t pivot = it->second;
	  interval_t i = compute_residual(e, pivot) / interval_t(it->first);
	  if (auto k = i.singleton()){
	    add_univar_disequation(pivot, *k);
	  }
	}
	return;
      }

      interval_t get_interval(variable_t x) { return get_interval(vert_map, graph, x); }
      interval_t get_interval(vert_map_t& m, graph_t& r, variable_t x) {
	auto it = m.find(x);
	if (it == m.end()) return interval_t::top();
	vert_id v = (*it).second.first;
	interval_t x_out = interval_t(
				      r.elem(v, v+1) ? -number_t(r.edge_val(v, v+1))/2 : bound_t::minus_infinity(),
				      r.elem(v+1, v) ? number_t(r.edge_val(v+1, v))/2 : bound_t::plus_infinity());
	return x_out;
      }

      bool repair_potential(vert_id src, vert_id dest) {
	return GrOps::repair_potential(graph, potential, src, dest);
      }

      void close_over_edge(vert_id ii, vert_id jj) {
	Wt_min min_op;

	assert(ii / 2 != jj / 2);
	SplitGraph<graph_t> g_excl(graph);

	Wt c = g_excl.edge_val(ii,jj);

	typename graph_t::mut_val_ref_t w;

        #ifdef CLOSE_BOUNDS_INLINE
	if (ii % 2 == 0 && jj % 2 == 0) {
	  if(graph.lookup(jj+1, ii, &w))
	    graph.update_edge(jj+1, w + c, jj, min_op); //ub of jj
	  if(graph.lookup(jj, ii+1, &w))
	    graph.update_edge(ii, w + c, ii+1, min_op); //lb of ii
	} else if (ii % 2 == 0 && jj % 2 != 0) {
	  if(graph.lookup(jj-1, ii, &w))
	    graph.update_edge(jj-1, w + c, jj, min_op); //lb of jj
	  if(graph.lookup(jj, ii+1, &w))
	    graph.update_edge(ii, w + c, ii+1, min_op); //lb of ii
	} else if (ii % 2 != 0 && jj % 2 == 0) {
	  if(graph.lookup(jj+1, ii, &w))
	    graph.update_edge(jj+1, w + c, jj, min_op); //ub of jj
	  if(graph.lookup(jj, ii-1, &w))
	    graph.update_edge(ii, w + c, ii-1, min_op); //ub of ii
	} else if (ii % 2 != 0 && jj % 2 != 0) {
	  if(graph.lookup(jj-1, ii, &w))
	    graph.update_edge(jj-1, w + c, jj, min_op); //lb of jj
	  if(graph.lookup(jj, ii-1, &w))
	    graph.update_edge(ii, w + c, ii-1, min_op); //ub of ii
	}
        #endif

	std::vector<std::pair<vert_id, Wt> > src_dec;
	for(auto edge : g_excl.e_preds(ii)) {
	  vert_id se = edge.vert;
	  Wt wt_sij = edge.val + c;

	  assert(g_excl.succs(se).begin() != g_excl.succs(se).end());
	  if(se != jj) {
	    if(g_excl.lookup(se, jj, &w)) {
	      if(w <= wt_sij) continue;
	      w = wt_sij;
	    } else {
	      g_excl.add_edge(se, wt_sij, jj);
	    }
	    src_dec.push_back(std::make_pair(se, edge.val));  
            #ifdef CLOSE_BOUNDS_INLINE
	    if (se % 2 == 0 && jj % 2 == 0) {
	      if(graph.lookup(jj+1, se, &w))
		graph.update_edge(jj+1, w + wt_sij, jj, min_op); //ub of jj
	      if(graph.lookup(jj, se+1, &w))
		graph.update_edge(se, w + wt_sij, se+1, min_op); //lb of se
	    } else if (se % 2 != 0 && jj % 2 == 0) {
	      if(graph.lookup(jj+1, se, &w))
		graph.update_edge(jj+1, w + wt_sij, jj, min_op); //ub of jj
	      if(graph.lookup(jj, se-1, &w))
		graph.update_edge(se, w + wt_sij, se-1, min_op); //ub of se
	    } else if (se % 2 == 0 && jj % 2 != 0) {
	      if(graph.lookup(jj-1, se, &w))
		graph.update_edge(jj-1, w + wt_sij, jj, min_op); //lb of jj
	      if(graph.lookup(jj, se+1, &w))
		graph.update_edge(se, w + wt_sij, se+1, min_op); //lb of se
	    } else if (se % 2 != 0 && jj % 2 != 0) {
	      if(graph.lookup(jj-1, se, &w))
		graph.update_edge(jj-1, w + wt_sij, jj, min_op); //lb of jj
	      if(graph.lookup(jj, se-1, &w))
		graph.update_edge(se, w + wt_sij, se-1, min_op); //ub of jj
	    }
            #endif
	  }
	}

	std::vector<std::pair<vert_id, Wt> > dest_dec;
	for(auto edge : g_excl.e_succs(jj)) {
	  vert_id de = edge.vert;
	  Wt wt_ijd = edge.val + c;
	  if(de != ii) {
	    if(g_excl.lookup(ii, de, &w)) {
	      if(w <= wt_ijd) continue;
	      w = wt_ijd;
	    } else {
	      g_excl.add_edge(ii, wt_ijd, de);
	    }
	    dest_dec.push_back(std::make_pair(de, edge.val));
            #ifdef CLOSE_BOUNDS_INLINE
	    if (de % 2 == 0 && ii % 2 == 0) {
	      if(graph.lookup(de+1, ii, &w))
		graph.update_edge(de+1, w + wt_ijd, de, min_op); //ub of de
	      if(graph.lookup(de, ii+1, &w))
		graph.update_edge(ii, w + wt_ijd, ii+1, min_op); //lb of ii
	    } else if (de % 2 != 0 && ii % 2 == 0) {
	      if(graph.lookup(de-1, ii, &w))
		graph.update_edge(de-1, w + wt_ijd, de, min_op); //lb of de
	      if(graph.lookup(de, ii+1, &w))
		graph.update_edge(ii, w + wt_ijd, ii+1, min_op); //lb of ii
	    } else if (de % 2 == 0 && ii % 2 != 0) {
	      if(graph.lookup(de+1, ii, &w))
		graph.update_edge(de+1, w + wt_ijd, de, min_op); //ub of de
	      if(graph.lookup(de, ii-1, &w))
		graph.update_edge(ii, w + wt_ijd, ii-1, min_op); //ub of ii
	    } else if (de % 2 != 0 && ii % 2 != 0) {
	      if(graph.lookup(de-1, ii, &w))
		graph.update_edge(de-1, w + wt_ijd, de, min_op); //lb of de
	      if(graph.lookup(de, ii-1, &w))
		graph.update_edge(ii, w + wt_ijd, ii-1, min_op); //ub of ii
	    }
            #endif
	  }
	}

	for(auto s_p : src_dec) {
	  vert_id se = s_p.first;
	  Wt wt_sij = c + s_p.second;
	  for(auto d_p : dest_dec) {
	    vert_id de = d_p.first;
	    Wt wt_sijd = wt_sij + d_p.second; 
	    if(graph.lookup(se, de, &w)) {
	      if(w <= wt_sijd) continue;
	      w = wt_sijd;
	    } else {
	      graph.add_edge(se, wt_sijd, de);
	    }
            #ifdef CLOSE_BOUNDS_INLINE
	    if (se % 2 == 0 && de % 2 == 0) {
	      if(graph.lookup(de+1, se, &w))
		graph.update_edge(de+1, w + wt_sijd, de, min_op); //ub of de
	      if(graph.lookup(de, se+1, &w))
		graph.update_edge(se, w + wt_sijd, se+1, min_op); //lb of se
	    } else if (se % 2 != 0 && de % 2 == 0) {
	      if(graph.lookup(de+1, se, &w))
		graph.update_edge(de+1, w + wt_sijd, de, min_op); //ub of de
	      if(graph.lookup(de, se-1, &w))
		graph.update_edge(se, w + wt_sijd, se-1, min_op); //ub of se
	    } else if (se % 2 == 0 && de % 2 != 0) {
	      if(graph.lookup(de-1, se, &w))
		graph.update_edge(de-1, w + wt_sijd, de, min_op); //lb of de
	      if(graph.lookup(de, se+1, &w))
		graph.update_edge(se, w + wt_sijd, se+1, min_op); //lb of se
	    } else if (se % 2 != 0 && de % 2 != 0) {
	      if(graph.lookup(de-1, se, &w))
		graph.update_edge(de-1, w + wt_sijd, de, min_op); //lb of de
	      if(graph.lookup(de, se-1, &w))
		graph.update_edge(se, w + wt_sijd, se-1, min_op); //ub of se
	    }
            #endif
	  }
	}
      }

      
    public:
      
      split_oct_domain(bool is_bottom = false)
	: _is_bottom(is_bottom) {}
      
      split_oct_domain(const split_oct_domain_t& o)
	: vert_map(o.vert_map), rev_map(o.rev_map), graph(o.graph),
	  potential(o.potential), unstable(o.unstable), _is_bottom(false) {
	crab::CrabStats::count(getDomainName() + ".count.copy");
	crab::ScopedCrabStats __st__(getDomainName() + ".copy");
	if (o._is_bottom) set_to_bottom();
	if (!_is_bottom) assert(graph.size() >= 0);
      }
      split_oct_domain(split_oct_domain_t&& o)
	: vert_map(std::move(o.vert_map)), rev_map(std::move(o.rev_map)),
	  graph(std::move(o.graph)), potential(std::move(o.potential)),
	  unstable(std::move(o.unstable)), _is_bottom(o._is_bottom)
      { }
      split_oct_domain(vert_map_t& _vert_map, rev_map_t& _rev_map, graph_t& _graph,
		       std::vector<Wt>& _potential, vert_set_t& _unstable)
	: vert_map(_vert_map), rev_map(_rev_map), graph(_graph),
	  potential(_potential), unstable(_unstable), _is_bottom(false) {
	CRAB_WARN("Non-moving constructor.");
	assert(graph.size() >= 0);
      }
      split_oct_domain(vert_map_t&& _vert_map, rev_map_t&& _rev_map, graph_t&& _graph,
		       std::vector<Wt>&& _potential, vert_set_t&& _unstable)
	: vert_map(std::move(_vert_map)), rev_map(std::move(_rev_map)), graph(std::move(_graph)),
	  potential(std::move(_potential)), unstable(std::move(_unstable)), _is_bottom(false)
      { assert(graph.size() >= 0); }

      split_oct_domain& operator=(const split_oct_domain& o) {
	crab::CrabStats::count(getDomainName() + ".count.copy");
	crab::ScopedCrabStats __st__(getDomainName() + ".copy");

	if (this != &o) {
	  if (o._is_bottom) set_to_bottom();
	  else {
	    _is_bottom = false;
	    vert_map = o.vert_map;
	    rev_map = o.rev_map;
	    graph = o.graph;
	    potential = o.potential;
	    unstable = o.unstable;
	    assert(graph.size() >= 0);
	  }
	}
	return *this;
      }

      split_oct_domain& operator=(split_oct_domain&& o) {
	if (o._is_bottom) {
	  set_to_bottom();
	} else {
	  _is_bottom = false;
	  vert_map = std::move(o.vert_map);
	  rev_map = std::move(o.rev_map);
	  graph = std::move(o.graph);
	  potential = std::move(o.potential);
	  unstable = std::move(o.unstable);
	}
	return *this;
      }
      
      void set_to_top() {
	split_oct_domain_t tmp(false);
	std::swap(*this, tmp);	
      }
      
      void set_to_bottom() {
	vert_map.clear();
	rev_map.clear();
	graph.clear();
	potential.clear();
	unstable.clear();
	_is_bottom = true;
      }

      bool is_bottom() {
	return _is_bottom;
      }
      
      bool is_top() {
	if (_is_bottom) return false;
	return graph.is_empty();
      }
      
      bool operator<=(split_oct_domain_t o) {
	crab::CrabStats::count(getDomainName() + ".count.leq");
	crab::ScopedCrabStats __st__(getDomainName() + ".leq");
	
	if (is_bottom()) return true;
	else if (o.is_bottom()) return false;
	else if (o.is_top()) return true;
	else if (is_top()) return false;
	else {
	  normalize();

	  if (vert_map.size() < o.vert_map.size()) return false;
	  typename graph_t::mut_val_ref_t wx;
	  typename graph_t::mut_val_ref_t wy;
	  typename graph_t::mut_val_ref_t wz;
	  bool broken = false;

	  std::vector<unsigned int> vert_renaming(o.graph.size(), -1);
	  for (auto p : o.vert_map) {
	    if (o.graph.succs(p.second.first).size() == 0
		&& o.graph.succs(p.second.second).size() == 0
		&& o.graph.preds(p.second.first).size() == 0
		&& o.graph.preds(p.second.second).size() == 0)
	      continue;

	    auto it = vert_map.find(p.first);
	    if (it == vert_map.end()) return false;
	    vert_renaming[p.second.first] = (*it).second.first;
	    vert_renaming[p.second.second] = (*it).second.second;
	  }

	  assert(graph.size() >= 0);

	  for (vert_id ox : o.graph.verts()) {
	    if (o.graph.succs(ox).size() == 0) continue;

	    assert(vert_renaming[ox] != -1);
	    vert_id x = vert_renaming[ox];
	    for (auto edge : o.graph.e_succs(ox)) {
	      broken = false;
	      vert_id oy = edge.vert;
	      assert(vert_renaming[oy] != -1);
	      if (ox == oy) continue;
	      vert_id y = vert_renaming[oy];
	      Wt ow = edge.val;

	      if(graph.lookup(x, y, &wx) && (wx <= ow)) continue;

	      for (vert_id gx : graph.verts()) {
		if (gx % 2 != 0) continue;
		if (graph.elem(gx, gx+1)) {
		  if (graph.lookup(x, gx, &wx) && graph.lookup(gx, gx+1, &wy) &&
		      graph.lookup(gx+1, y, &wz)) {
		    if (wx + wy + wz <= ow) {
		      broken = true;
		      break;
		    }
		  }
		}
		if (graph.elem(gx+1, gx)) {
		  if (graph.lookup(x, gx+1, &wx) && graph.lookup(gx+1, gx, &wy) &&
		      graph.lookup(gx, y, &wz))
		    if (wx + wy + wz <= ow) {
		      broken = true;
		      break;
		    }
		}
	      }

	      if (!broken) return false;
	    }
	  }
	  return true;
	}
      }

      void operator|=(split_oct_domain_t o) {
	*this = *this | o;
      }

      split_oct_domain_t operator|(split_oct_domain_t o) {
	crab::CrabStats::count(getDomainName() + ".count.join");
	crab::ScopedCrabStats __st__(getDomainName() + ".join");
		    	
	if (is_bottom() || o.is_top ()) return o;
	else if (is_top() || o.is_bottom()) return *this;
	else {
	  CRAB_LOG("octagon-split", 
		    crab::outs() << "Before join:\n"<<"DBM 1\n"<<*this<<"\n"<<graph<<"\n"<<"DBM 2\n"
		    << o <<"\n"<<o.graph<<"\n");
	  CRAB_LOG("octagon-join", 
		    crab::outs() << "Before join:\n"<<"DBM 1\n"<<*this<<"\n"<<graph<<"\n"<<"DBM 2\n"
		    << o <<"\n"<<o.graph<<"\n");

	  normalize();
	  o.normalize();

	  assert(check_potential(graph, potential));
	  assert(check_potential(o.graph, o.potential));

	  std::vector<vert_id> perm_x;
	  std::vector<vert_id> perm_y;
	  std::vector<variable_t> perm_inv;

	  std::vector<Wt> pot_rx;
	  std::vector<Wt> pot_ry;
	  vert_map_t out_vmap;
	  rev_map_t out_revmap;

	  for (auto p : vert_map) {
	    auto it = o.vert_map.find(p.first);
	    if(it != o.vert_map.end()) {
	      out_vmap.insert(vmap_elt_t(p.first, std::make_pair(perm_x.size(), perm_x.size()+1)));
	      out_revmap.push_back(p.first);
	      out_revmap.push_back(p.first);

	      pot_rx.push_back(potential[p.second.first]);
	      pot_rx.push_back(potential[p.second.second]);
	      pot_ry.push_back(o.potential[(*it).second.first]);
	      pot_ry.push_back(o.potential[(*it).second.second]);
	      perm_inv.push_back(p.first);
	      perm_x.push_back(p.second.first);
	      perm_x.push_back(p.second.second);
	      perm_y.push_back((*it).second.first);
	      perm_y.push_back((*it).second.second);
	    }
	  }
	  
	  unsigned int sz = perm_x.size();
					
	  // Build the permuted view of x and y.
	  assert(graph.size() > 0);
	  GrPerm gx(perm_x, graph);
	  assert(o.graph.size() > 0);
	  GrPerm gy(perm_y, o.graph);

	  // Compute the deferred relations
	  graph_t g_ix_ry;
	  g_ix_ry.growTo(sz);
	  SplitGraph<GrPerm> gy_excl(gy);
	  for (vert_id s : gy_excl.verts()) {
	    for (vert_id d : gy_excl.succs(s)) {
	      typename graph_t::mut_val_ref_t ws;
	      typename graph_t::mut_val_ref_t wd;
	      typename graph_t::mut_val_ref_t wx;
	      CRAB_LOG("octagon-join", crab::outs()<<"s,d: "<<s<<", "<<d<<"\n");
	      if (s%2==0 && d%2==0){ // both pos
		if (gx.lookup(s, s+1, &ws) && gx.lookup(d+1, d, &wd)){
		  g_ix_ry.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2==0 && d%2!=0){ // s pos
		if (gx.lookup(s, s+1, &ws) && gx.lookup(d-1, d, &wd)){
		  g_ix_ry.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2!=0 && d%2==0){ // d pos
		if (gx.lookup(s, s-1, &ws) && gx.lookup(d+1, d, &wd)){
		  g_ix_ry.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2!=0 && d%2!=0){ // both neg
		if (gx.lookup(s, s-1, &ws) && gx.lookup(d-1, d, &wd)){
		  g_ix_ry.add_edge(s, (ws+wd)/2 ,d);
		}
	      }

	    }
	  }

	  // Apply the deferred relations, and re-close.
	  edge_vector delta;
	  bool is_closed;
	  graph_t g_rx(GrOps::meet(gx, g_ix_ry, is_closed));
	  assert(check_potential(g_rx, pot_rx));
	  if (!is_closed) {
	    SplitGraph<graph_t> g_rx_excl(g_rx);
	    GrOps::close_after_meet(g_rx_excl, pot_rx, gx, g_ix_ry, delta);
	    GrOps::update_delta(g_rx, delta);
	  }

	  // Compute the deferred relations
	  graph_t g_rx_iy;
	  g_rx_iy.growTo(sz);
	  SplitGraph<GrPerm> gx_excl(gx);
	  for (vert_id s : gx_excl.verts()) {
	    for (vert_id d : gx_excl.succs(s)) {
	      typename graph_t::mut_val_ref_t ws;
	      typename graph_t::mut_val_ref_t wd;
	      typename graph_t::mut_val_ref_t wy;
	      if (s%2==0 && d%2==0){ // both pos
		if (gy.lookup(s, s+1, &ws) && gy.lookup(d+1, d, &wd)){
		  CRAB_LOG("octagon-join", crab::outs()<<"s -min-> d: "<<s<<" -"<<(ws+wd)/2<<"-> "<<d<<"\n");
		  g_rx_iy.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2==0 && d%2!=0){ // s pos
		if (gy.lookup(s, s+1, &ws) && gy.lookup(d-1, d, &wd)){
		  CRAB_LOG("octagon-join", crab::outs()<<"s -min-> d: "<<s<<" -"<<(ws+wd)/2<<"-> "<<d<<"\n");
		  g_rx_iy.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2!=0 && d%2==0){ // d pos
		if (gy.lookup(s, s-1, &ws) && gy.lookup(d+1, d, &wd)){
		  CRAB_LOG("octagon-join", crab::outs()<<"s -min-> d: "<<s<<" -"<<(ws+wd)/2<<"-> "<<d<<"\n");
		  g_rx_iy.add_edge(s, (ws+wd)/2 ,d);
		}
	      } else if (s%2!=0 && d%2!=0){ // both neg
		if (gy.lookup(s, s-1, &ws) && gy.lookup(d-1, d, &wd)){
		  CRAB_LOG("octagon-join", crab::outs()<<"s -min-> d: "<<s<<" -"<<(ws+wd)/2<<"-> "<<d<<"\n");
		  g_rx_iy.add_edge(s, (ws+wd)/2 ,d);
		}
	      }
	    }
	  }

	  // Apply the deferred relations, and re-close.
	  delta.clear();
	  graph_t g_ry(GrOps::meet(gy, g_rx_iy, is_closed));
	  CRAB_LOG("octagon-join", crab::outs()<<"\tMET g_ry:\n\t"<<g_ry<<"\n");
	  if (!is_closed) {
	    SplitGraph<graph_t> g_ry_excl(g_ry);
	    GrOps::close_after_meet(g_ry_excl, pot_ry, gy, g_rx_iy, delta);
	    GrOps::update_delta(g_ry, delta);
	  }

	  // We now have the relevant set of relations. Because g_rx and g_ry are closed,
	  // the result is also closed.
	  Wt_min min_op;
	  CRAB_LOG("octagon-join",crab::outs()<<"\tBefore joined:\n\tg_rx:"<<g_rx<<"\n\tg_ry:"<<g_ry<<"\n");
	  graph_t join_g(GrOps::join(g_rx, g_ry));
	  CRAB_LOG("octagon-join", crab::outs()<<"Joined graph:\n"<<join_g<<"\n");

	  std::vector<vert_id> lb_up;
	  std::vector<vert_id> lb_down;
	  std::vector<vert_id> ub_up;
	  std::vector<vert_id> ub_down;

	  typename graph_t::mut_val_ref_t wx;
	  typename graph_t::mut_val_ref_t wy;
	  for (vert_id v : gx.verts()) {
	    if (v % 2 != 0) continue;
	    if (gx.lookup(v+1, v, &wx) && gy.lookup(v+1, v, &wy)) {
	      if (wx < wy) ub_up.push_back(v);
	      if (wy < wx) ub_down.push_back(v);
	    }
	    if (gx.lookup(v, v+1, &wx) && gy.lookup(v, v+1, &wy)) {
	      if (wx < wy) lb_down.push_back(v);
	      if (wy < wx) lb_up.push_back(v);
	    }
	  }

	  for(vert_id s : lb_up) {
	    Wt dx_s = gx.edge_val(s, s+1)/2;
	    Wt dy_s = gy.edge_val(s, s+1)/2;
	    for(vert_id d : ub_up) {
	      if(s == d) continue;
	      join_g.update_edge(s, std::max(dx_s + gx.edge_val(d+1, d)/2,
					     dy_s + gy.edge_val(d+1, d)/2),
				 d, min_op);
	    }
	    for(vert_id d : lb_down) {
	      if(s == d) continue;
	      join_g.update_edge(d, std::max(dx_s + gx.edge_val(d, d+1)/2,
					     dy_s + gy.edge_val(d, d+1)/2),
				 s+1, min_op);
	    }
	  }

	  for(vert_id s : lb_down) {
	    Wt dx_s = gx.edge_val(s, s+1)/2;
	    Wt dy_s = gy.edge_val(s, s+1)/2;
	    for(vert_id d : ub_down) {
	      if(s == d) continue;
	      join_g.update_edge(s, std::max(dx_s + gx.edge_val(d+1, d)/2,
					     dy_s + gy.edge_val(d+1, d)/2),
				 d, min_op);
	    }
	  }
                                        
	  for(vert_id s : ub_up) {
	    Wt dx_s = gx.edge_val(s+1, s)/2;
	    Wt dy_s = gy.edge_val(s+1, s)/2;
	    for(vert_id d : ub_down) {
	      if(s == d) continue;
	      join_g.update_edge(s+1, std::max(dx_s + gx.edge_val(d+1, d)/2,
					       dy_s + gy.edge_val(d+1, d)/2),
				 d, min_op);
	    }
	  }

	  // Conjecture: join_g remains closed.
	  // Now garbage collect any unused vertices
	  CRAB_LOG("octagon-join", crab::outs()<<"Joined graph:\n"<<join_g<<"\n");
	  //CRAB_LOG("octagon-join", crab::outs()<<"Start garbage collection:\n");
	  for(vert_id v : join_g.verts()) {
	    //CRAB_LOG("octagon-join", crab::outs()<<"\tCurrent vertex: "<<v<<"\n");
	    if (v % 2 != 0) continue;
	    if(join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0 && join_g.succs(v+1).size() == 0 && join_g.preds(v+1).size() == 0) {
	      join_g.forget(v);
	      join_g.forget(v+1);
	      if(out_revmap[v]) {
		out_vmap.erase(*(out_revmap[v]));
		out_revmap[v] = boost::none;
		out_revmap[v+1] = boost::none;
	      }
	    }
	  }
          
	  // split_oct_domain_t res(join_range, out_vmap, out_revmap, join_g, join_pot);
	  split_oct_domain_t res(std::move(out_vmap), std::move(out_revmap), std::move(join_g), 
				 std::move(pot_rx), vert_set_t());
	  //join_g.check_adjs();
	  CRAB_LOG("octagon-split", crab::outs() << "Result join:\n"<<res <<"\n");
	  CRAB_LOG("octagon-join", crab::outs() << "Result join:\n"<<res <<"\n");

	  return res;
	}
      }

      // widening
      split_oct_domain_t operator||(split_oct_domain_t o) {
	crab::CrabStats::count(getDomainName() + ".count.widening");
	crab::ScopedCrabStats __st__(getDomainName() + ".widening");
		    	
	if (is_bottom()) return o;
	else if (o.is_bottom()) return *this;
	else {
	  CRAB_LOG ("octagon-split",
		    crab::outs() << "Before widening:\n"<<"DBM 1\n"<<*this<<"\n"<<"DBM 2\n"
		    <<o <<"\n");
	  o.normalize();

	  // Figure out the common renaming
	  std::vector<vert_id> perm_x;
	  std::vector<vert_id> perm_y;
	  vert_map_t out_vmap;
	  rev_map_t out_revmap;
	  std::vector<Wt> widen_pot;
	  vert_set_t widen_unstable(unstable);

	  assert(potential.size() > 0);
	  for(auto p : vert_map) {
	    auto it = o.vert_map.find(p.first); 
	    // Variable exists in both
	    if(it != o.vert_map.end()) {
	      out_vmap.insert(vmap_elt_t(p.first,
					 std::make_pair(perm_x.size(), perm_x.size() + 1)));
	      out_revmap.push_back(p.first);
	      out_revmap.push_back(p.first);

	      widen_pot.push_back(potential[p.second.first]);
	      widen_pot.push_back(potential[p.second.second]);
	      perm_x.push_back(p.second.first);
	      perm_x.push_back(p.second.second);
	      perm_y.push_back((*it).second.first);
	      perm_y.push_back((*it).second.second);
	    }
	  }

	  // Build the permuted view of x and y.
	  assert(graph.size() > 0);
	  GrPerm gx(perm_x, graph);            
	  assert(o.graph.size() > 0);
	  GrPerm gy(perm_y, o.graph);

	  // Now perform the widening 
	  std::vector<vert_id> destabilized;
	  graph_t widen_g(GrOps::widen(gx, gy, destabilized));
	  for(vert_id v : destabilized) widen_unstable.insert(v);

	  split_oct_domain_t res(std::move(out_vmap), std::move(out_revmap), std::move(widen_g), 
		    std::move(widen_pot), std::move(widen_unstable));

	  CRAB_LOG("octagon-split", crab::outs() << "Result widening:\n"<<res <<"\n");
	  return res;
	}
      }

      // MEET
      split_oct_domain_t operator&(split_oct_domain_t o) {
	crab::CrabStats::count(getDomainName() + ".count.meet");
	crab::ScopedCrabStats __st__(getDomainName() + ".meet");
		    	
	if (is_bottom() || o.is_bottom()) return split_oct_domain_t::bottom();
	else if (is_top()) return o;
	else if (o.is_top()) return *this;
	else{
	  CRAB_LOG("octagon-split",
		    crab::outs() << "Before meet:\n"<<"DBM 1\n"<<*this<<"\n"<<"DBM 2\n"<<o
		    <<"\n");
	  normalize();
	  o.normalize();

	  vert_map_t meet_verts;
	  rev_map_t meet_rev;

	  std::vector<vert_id> perm_x;
	  std::vector<vert_id> perm_y;
	  std::vector<Wt> meet_pi;

	  for (auto p : vert_map) {
	    vert_id vv = perm_x.size();
	    meet_verts.insert(vmap_elt_t(p.first, std::make_pair(vv, vv+1)));
	    meet_rev.push_back(p.first);
	    meet_rev.push_back(p.first);

	    perm_x.push_back(p.second.first);
	    perm_x.push_back(p.second.second);
	    perm_y.push_back(-1);
	    perm_y.push_back(-1);
	    meet_pi.push_back(potential[p.second.first]);
	    meet_pi.push_back(potential[p.second.second]);
	  }
	  
	  // Add missing mappings from the right operand.
	  for(auto p : o.vert_map) {
	    auto it = meet_verts.find(p.first);

	    if(it == meet_verts.end()) {
	      vert_id vv = perm_y.size();
	      meet_rev.push_back(p.first);
	      meet_rev.push_back(p.first);

	      perm_y.push_back(p.second.first);
	      perm_y.push_back(p.second.second);
	      perm_x.push_back(-1);
	      perm_x.push_back(-1);
	      meet_pi.push_back(o.potential[p.second.first]);
	      meet_pi.push_back(o.potential[p.second.second]);
	      meet_verts.insert(vmap_elt_t(p.first, std::make_pair(vv, vv+1)));
	    } else {
	      perm_y[(*it).second.first] = p.second.first;
	      perm_y[(*it).second.second] = p.second.second;
	    }
	  }

	  // Build the permuted view of x and y.
	  assert(graph.size() > 0);
	  GrPerm gx(perm_x, graph);
	  assert(o.graph.size() > 0);
	  GrPerm gy(perm_y, o.graph);

	  // Compute the syntactic meet of the permuted graphs.
	  bool is_closed;
	  graph_t meet_g(GrOps::meet(gx, gy, is_closed));
           
	  // Compute updated potentials on the zero-enriched graph
	  //vector<Wt> meet_pi(meet_g.size());
	  // We've warm-started pi with the operand potentials
	  if(!GrOps::select_potentials(meet_g, meet_pi)) {
	    // Potentials cannot be selected -- state is infeasible.
	    return split_oct_domain_t::bottom();
	  }

	  if(!is_closed) {
	    edge_vector delta;
	    SplitGraph<graph_t> meet_g_excl(meet_g);
	    // GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);

	    if(Params::chrome_dijkstra) {
	      GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
	    } else {
	      GrOps::close_johnson(meet_g_excl, meet_pi, delta);
	    }
	    
	    GrOps::update_delta(meet_g, delta);
	    
	    // Recover updated LBs and UBs.
            #ifdef CLOSE_BOUNDS_INLINE
	    Wt_min min_op;
	    for(auto e : delta) {
	      if (e.first.first % 2 == 2) {
		if(meet_g.elem(e.first.first+1, e.first.first))
		  meet_g.update_edge(e.first.first+1,
				     meet_g.edge_val(e.first.first+1, e.first.first) + 2*e.second,
				     e.first.second, min_op);
		if(meet_g.elem(e.first.second, e.first.first+1))
		  meet_g.update_edge(e.first.first,
				     meet_g.edge_val(e.first.second, e.first.first+1) + 2*e.second,
				     e.first.first+1, min_op);
	      } else if (e.first.first % 2 == 2) {
		if(meet_g.elem(e.first.first, e.first.first+1))
		  meet_g.update_edge(e.first.first,
				     meet_g.edge_val(e.first.first, e.first.first+1) + 2*e.second,
				     e.first.second, min_op);
		if(meet_g.elem(e.first.second, e.first.first+1))
		  meet_g.update_edge(e.first.first+1,
				     meet_g.edge_val(e.first.second+1, e.first.first) + 2*e.second,
				     e.first.first, min_op);
	      }
	    }
            #else
	    delta.clear();
	    GrOps::close_after_assign(meet_g, meet_pi, 0, delta);
	    GrOps::update_delta(meet_g, delta);
            #endif
	  }
					
	  assert(check_potential(meet_g, meet_pi)); 
	  split_oct_domain_t res(std::move(meet_verts), std::move(meet_rev), std::move(meet_g), 
				 std::move(meet_pi), vert_set_t());
	  CRAB_LOG("octagon-split",
		    crab::outs() << "Result meet:\n"<<res <<"\n");
	  return res;

	}
      }

      template<typename Thresholds>
      split_oct_domain_t widening_thresholds(split_oct_domain_t o, const Thresholds& ts) {
	return (*this || o);
      }
		    
      
      // NARROWING
      split_oct_domain_t operator&&(split_oct_domain_t o) {
	crab::CrabStats::count(getDomainName() + ".count.narrowing");
	crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");
		    	
	if (is_bottom() || o.is_bottom()) return split_oct_domain_t::bottom();
	else if (is_top()) return o;
	else{
	  CRAB_LOG ("octagon-split",
		    crab::outs() << "Before narrowing:\n"<<"DBM 1\n"<<*this<<"\n"<<"DBM 2\n"
		    << o <<"\n");

	  // FIXME: Implement properly
	  // Narrowing as a no-op should be sound.
	  normalize();
	  split_oct_domain_t res(*this);

	  CRAB_LOG("octagon-split",
		    crab::outs() << "Result narrowing:\n"<<res <<"\n");
	  return res;
	}
      }

      void operator-=(variable_t v) {
	if (is_bottom()) return;
	normalize();
	auto it = vert_map.find(v);
	if (it != vert_map.end()) {
	  graph.forget((*it).second.first);
	  graph.forget((*it).second.second);
	  rev_map[(*it).second.first] = boost::none;
	  rev_map[(*it).second.second] = boost::none;
	  vert_map.erase(v);
	}
      }

      void operator+=(linear_constraint_t cst) {
	crab::CrabStats::count(getDomainName() + ".count.add_constraints");
	crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");
	if (is_bottom()) return;
	normalize();
	if (cst.is_tautology()) return;
	if (cst.is_contradiction()) {
	  set_to_bottom();
	  return;
	}

	if (cst.is_inequality()) {
	  if (!add_linear_leq(cst.expression())) set_to_bottom();
	  CRAB_LOG("octagon-split",
		   crab::outs() << "--- "<< cst<< "\n"<< *this <<"\n");
	  return;
	}

	if (cst.is_equality()) {
	  linear_expression_t exp = cst.expression();
	  if (!add_linear_leq(exp) || !add_linear_leq(-exp)) {
	    CRAB_LOG("octagon-split", crab::outs() << " ~~> _|_" <<"\n");
	    set_to_bottom();
	  }
	  CRAB_LOG("octagon-split",
		   crab::outs() << "--- "<< cst<< "\n"<< *this <<"\n");
	  return;
	}

	if (cst.is_disequation()) {
	  add_disequation(cst.expression());
	  CRAB_LOG("octagon-split",
		   crab::outs() << "--- "<< cst<< "\n"<< *this <<"\n");
	  return;
	}
        		
	CRAB_WARN("Unhandled constraint in SplitOct");

	CRAB_LOG("octagon-split", crab::outs() << "---"<< cst<< "\n"<< *this <<"\n");
	return;
      }

      void operator+=(linear_constraint_system_t csts) {
	if (is_bottom()) return;
	for (auto cst : csts) {
	  operator+=(cst);
	}
      }

      interval_t operator[](variable_t x) {
	crab::CrabStats::count(getDomainName() + ".count.to_intervals");
	crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");
	if (is_bottom()) return interval_t::bottom();
	if (this->is_bottom()) {
	  return interval_t::bottom();
	} else {
	  //return get_interval(ranges, x);
	  return get_interval(vert_map, graph, x);
	}
      }

      void normalize() {
	Wt_min min_op;
	for (vert_id v : graph.verts()){
	  for (vert_id w : graph.succs(v)){
	    if (v / 2 == w / 2) continue;
	    Wt current = graph.edge_val(v, w);
	    typename graph_t::mut_val_ref_t mirror;	
	    if (v % 2 == 0 and w % 2 == 0) {
	      if (graph.lookup(w+1, v+1, &mirror)){
		graph.update_edge(w+1, (mirror < current ? mirror : current), v+1, min_op);
		graph.update_edge(v, (mirror < current ? mirror : current), w, min_op);
	      } else {
		graph.add_edge(w+1, current, v+1);
		if(!repair_potential(w+1, v+1)) {set_to_bottom();}
	      }
	    } else if (v % 2 != 0 and w % 2 == 0) {
	      if (graph.lookup(w+1, v-1, &mirror)){
		graph.update_edge(w+1, (mirror < current ? mirror : current), v-1, min_op);
		graph.update_edge(v, (mirror < current ? mirror : current), w, min_op);
	      } else {
		graph.add_edge(w+1, current, v-1);
		if(!repair_potential(w+1, v-1)) {set_to_bottom();}
	      }
	    } else if (v % 2 == 0 and w % 2 != 0) {
	      if (graph.lookup(w-1, v+1, &mirror)){
		graph.update_edge(w-1, (mirror < current ? mirror : current), v+1, min_op);
		graph.update_edge(v, (mirror < current ? mirror : current), w, min_op);
	      } else {
		graph.add_edge(w-1, current, v+1);
		if(!repair_potential(w-1, v+1)) {set_to_bottom();}
	      }
	    } else if (v % 2 != 0 and w % 2 != 0) {
	      if (graph.lookup(w-1, v-1, &mirror)){
		graph.update_edge(w-1, (mirror < current ? mirror : current), v-1, min_op);
		graph.update_edge(v, (mirror < current ? mirror : current), w, min_op);
	      } else {
		graph.add_edge(w-1, current, v-1);
		if(!repair_potential(w-1, v-1)) {set_to_bottom();}
	      }
	    }
	  }
	}
	//CRAB_LOG("octagon-split", crab::outs() << "--- graph after normalize " << graph <<"\n");
        #ifdef SOCT_NO_NORMALIZE
	return;
        #endif
	
	CRAB_LOG("octagon-unstable",
		 crab::outs() << "Size of unstable list: " << unstable.size()
		              << "\n\tWith domaim:\n\t" << *this << "\n\tand graph:\n\t" << graph <<"\n");
	if (unstable.size() == 0) return;
		    	
	edge_vector delta;
	SplitGraph<graph_t> g_excl(graph);
	if (Params::widen_restabilize)
	  GrOps::close_after_widen(g_excl,potential, vert_set_wrap_t(unstable), delta);
	else GrOps::close_johnson(g_excl, potential, delta);
	GrOps::update_delta(graph, delta);

	unstable.clear();
      }

      // set a variable to an interval
      void set(variable_t x, interval_t intv) {
	crab::CrabStats::count(getDomainName() + ".count.assign");
	crab::ScopedCrabStats __st__(getDomainName() + ".assign");
	if (is_bottom()) return;
	this->operator-=(x);
	if (intv.is_top()) return;
	vert_id v = get_vert(x);
	if (intv.ub().is_finite()) {
	  Wt ub = ntov::ntov(*(intv.ub().number()));
	  potential[v] = ub;
	  potential[v+1] = -ub;
	  graph.set_edge(v+1, 2*ub, v);
	}
	if (intv.lb().is_finite()) {
	  Wt lb = ntov::ntov(*(intv.lb().number()));
	  potential[v] = lb;
	  potential[v+1] = -lb;
	  graph.set_edge(v, -2*lb, v+1);
	}
      }

      // assign am exact expression to a variable
      void assign(variable_t x, linear_expression_t e) {
	crab::CrabStats::count(getDomainName() + ".count.assign");
	crab::ScopedCrabStats __st__(getDomainName() + ".assign");

	if (is_bottom()) return;
	CRAB_LOG("octagon-split", crab::outs() << "Before assign: "<< *this <<"\n");
	CRAB_LOG("octagon-split", crab::outs() << x<< ":="<< e <<"\n");
	CRAB_LOG("octagon-assign", crab::outs() << "Before assign: "<< *this <<"\n");
	CRAB_LOG("octagon-assign", crab::outs() << x<< ":="<< e <<"\n");
	normalize();
	//assert(check_potential(graph, potential));

	if (e.is_constant()) set(x, e.constant());
	else {
	  interval_t x_int = eval_interval(e);
	  CRAB_LOG("octagon-assign", crab::outs() << "Interval of e is " << x_int <<"\n");
	  std::vector<std::pair<variable_t, Wt>> diffs_lb;
	  std::vector<std::pair<variable_t, Wt>> diffs_ub;
	  diffcsts_of_assign(x, e, diffs_lb, diffs_ub);
	  if (diffs_lb.size() > 0 || diffs_ub.size() > 0) {
	    if (Params::special_assign) {
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning " << x << " to " << x_int <<"\n");
	      vert_id v = graph.new_vertex();
	      vert_id w = graph.new_vertex();
                                                        
	      if (w < v){
		vert_id tmp = v;
		v = w;
		w = tmp;
	      }
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning pos " << v << ", neg " << w <<"\n");
	      // assert(w <= rev_map.size());
	      if (v == rev_map.size()) {
		rev_map.push_back(x);
		potential.push_back(eval_expression(e));
	      } else {
		potential[v] = eval_expression(e);
		rev_map[v] = x;
	      }
	      if (w == rev_map.size()) {
		rev_map.push_back(x);
		potential.push_back(-eval_expression(e));
	      } else {
		potential[w] = -eval_expression(e);
		rev_map[w] = x;
	      }

	      edge_vector delta;
	      for (auto diff : diffs_lb) {
		delta.push_back(std::make_pair(std::make_pair(v, get_vert(diff.first)),
					       -diff.second));
	      }

	      for (auto diff : diffs_ub) {
		delta.push_back(std::make_pair(std::make_pair(get_vert(diff.first), v),
					       diff.second));
	      }
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning before meet " << *this <<"\n");
	      GrOps::apply_delta(graph, delta);
	      delta.clear();
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning appied delta " << *this <<"\n");
	      SplitGraph<graph_t> g_excl(graph);
	      GrOps::close_after_assign(g_excl, potential, v, delta);
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning closed " << *this <<"\n");
	      GrOps::update_delta(graph, delta);
	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning appied closed delta " << *this << "\nwith graph: "<< graph <<"\n");

	      Wt_min min_op;
							
                                                        
	      if (x_int.lb().is_finite()) {
		CRAB_LOG("octagon-assign", crab::outs() << "Assigning " << x << " lb is " << ntov::ntov(-(*(x_int.lb().number()))) <<"\n");
		CRAB_LOG("octagon-assign", crab::outs() << "Assigning pos " << v << ", neg " << w <<"\n");
		graph.update_edge(v, 2*ntov::ntov(-(*(x_int.lb().number()))), w, min_op);
	      }
	      if (x_int.ub().is_finite()) {
		CRAB_LOG("octagon-assign", crab::outs() << "Assigning " << x << " ub is " << ntov::ntov(*(x_int.ub().number())) <<"\n");
		graph.update_edge(w, 2*ntov::ntov(*(x_int.ub().number())), v, min_op);
	      }

	      CRAB_LOG("octagon-assign", crab::outs() << "Assigning " << x << " to " << x_int << ", edge updated\n");
	      operator-=(x);
	      CRAB_LOG("octagon-assign", crab::outs() << "Add var to vert_map"<<"\n");
	      vert_map.insert(vmap_elt_t(x, std::make_pair(v, w)));
	      CRAB_LOG("octagon-assign", crab::outs() << "Added var to vert_map"<<"\n");
	    } else {
	      vert_id v = graph.new_vertex();
	      vert_id w = graph.new_vertex();
	      assert(w <= rev_map.size());
	      if (v == rev_map.size()) {
		rev_map.push_back(x);
		potential.push_back(Wt(0));
	      } else {
		potential[v] = Wt(0);
		rev_map[v] = x;
	      }
	      if (w == rev_map.size()) {
		rev_map.push_back(x);
		potential.push_back(Wt(0));
	      } else {
		potential[w] = Wt(0);
		rev_map[w] = x;
	      }
	      Wt_min min_op;
	      edge_vector cst_edges;
	      for (auto diff : diffs_lb) {
		cst_edges.push_back(std::make_pair(std::make_pair(v, get_vert(diff.first)),
						   -diff.second));
	      }

	      for (auto diff : diffs_ub) {
		cst_edges.push_back(std::make_pair(std::make_pair(get_vert(diff.first), v),
						   diff.second));
	      }

	      for(auto diff : cst_edges) {
		vert_id src = diff.first.first;
		vert_id dest = diff.first.second;
		graph.update_edge(src, diff.second, dest, min_op);
		if(!repair_potential(src, dest)) {
		  // assert(0 && "Unreachable");
		  set_to_bottom();
		}
		assert(check_potential(graph, potential));

		close_over_edge(src, dest);
		//assert(check_potential(graph, potential));
	      }

	      if(x_int.lb().is_finite())
		graph.update_edge(v, ntov::ntov(-(*(x_int.lb().number()))), w, min_op);
	      if(x_int.ub().is_finite())
		graph.update_edge(w, ntov::ntov(*(x_int.ub().number())), v, min_op);

              // Clear the old x vertex
	      operator-=(x);
	      vert_map.insert(vmap_elt_t(x, std::make_pair(v, w)));
	    }
	  } else {
	    set(x, x_int);
	  }
	}
	assert(check_potential(graph, potential));
	CRAB_LOG("octagon-split", crab::outs() << "---"<< x<< ":="<< e<<"\n"<<*this <<"\n");
	CRAB_LOG("octagon-assign", crab::outs() << "---"<< x<< ":="<< e<<"\n"<<*this <<"\n");
      }

      void forget(const variable_vector_t& vars) {
	if (is_bottom()) return;
	for (auto v : vars) {
	  auto it = vert_map.find(v);
	  if (it != vert_map.end()) {
	    operator-=(v);
	  }
	}
      }

      void project (const variable_vector_t& vars) {
	crab::CrabStats::count(getDomainName() + ".count.project");
	crab::ScopedCrabStats __st__(getDomainName() + ".project");
	if (is_bottom()) return;
	normalize();
	std::vector<bool> save(rev_map.size(), false);
	for (auto x : vars) {
	  auto it = vert_map.find(x);
	  if (it != vert_map.end()) {
	    save[(*it).second.first] = true;
	    save[(*it).second.second] = true;
	  }
	}
	for(vert_id v = 0; v < rev_map.size(); v++) {
	  if(!save[v] && rev_map[v])
	    operator-=((*rev_map[v]));
	}
      }

      void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
	crab::CrabStats::count(getDomainName() + ".count.apply");
	crab::ScopedCrabStats __st__(getDomainName() + ".apply");
	if (is_bottom()) return;

	normalize();

        switch(op) {
          case OP_ADDITION:
            assign(x, y + z);
            break;
          case OP_SUBTRACTION:
            assign(x, y - z);
            break;
          // For the rest of operations, we fall back on intervals.
          case OP_MULTIPLICATION:
            set(x, get_interval(y)*get_interval(z));
            break;
          case OP_SDIV:
            set(x, get_interval(y)/get_interval(z));
            break;
          case OP_UDIV:
            set(x, get_interval(y).UDiv(get_interval(z)));
            break;
          case OP_SREM:
            set(x, get_interval(y).SRem(get_interval(z)));
            break;
          case OP_UREM:
            set(x, get_interval(y).URem(get_interval(z)));
            break;
	  default:
	    CRAB_ERROR("Operation ", op, " not supported");
        }
	
	CRAB_LOG("octagon-split",
		 crab::outs() << "---"<< x<< ":="<< y<< op<< z<<"\n"<< *this <<"\n");

      }

      void apply(operation_t op, variable_t x, variable_t y, number_t k) {
	crab::CrabStats::count(getDomainName() + ".count.apply");
	crab::ScopedCrabStats __st__(getDomainName() + ".apply");
	if (is_bottom()) return;
	normalize();

        switch(op) {
          case OP_ADDITION:
            assign(x, y + k);
            break;
          case OP_SUBTRACTION:
            assign(x, y - k);
            break;
          // For the rest of operations, we fall back on intervals.
          case OP_MULTIPLICATION:
            set(x, get_interval(y)*interval_t(k));
            break;
          case OP_SDIV:
            set(x, get_interval(y)/interval_t(k));
            break;
          case OP_UDIV:
            set(x, get_interval(y).UDiv(interval_t(k)));
            break;
          case OP_SREM:
            set(x, get_interval(y).SRem(interval_t(k)));
            break;
          case OP_UREM:
            set(x, get_interval(y).URem(interval_t(k)));
            break;
	  default:
	    CRAB_ERROR("Operation ", op, " not supported");
        }
	CRAB_LOG("octagon-split",
		 crab::outs() << "---"<< x<< ":="<< y<< op<< k<<"\n"<< *this <<"\n");
      }

      void apply(int_conv_operation_t /*op*/, variable_t dst, variable_t src) {
	// since reasoning about infinite precision we simply assign and
	// ignore the widths.
	assign(dst, src);
      }

      void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
	crab::CrabStats::count(getDomainName() + ".count.apply");
	crab::ScopedCrabStats __st__(getDomainName() + ".apply");
	normalize();
	this->operator-=(x);
        		
	interval_t yi = operator[](y);
	interval_t zi = operator[](z);
	interval_t xi = interval_t::bottom();

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
	default: CRAB_ERROR("Oct: unreachable");
	}
	set(x, xi);
      }

      void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
	crab::CrabStats::count(getDomainName() + ".count.apply");
	crab::ScopedCrabStats __st__(getDomainName() + ".apply");
	// Convert to intervals and perform the operation
	normalize();
	interval_t yi = operator[](y);
	interval_t zi(k);
	interval_t xi = interval_t::bottom();

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
	default: CRAB_ERROR("Oct: unreachable");
	}
	set(x, xi);
      }

      void backward_assign(variable_t x, linear_expression_t e, split_oct_domain_t inv) {
	crab::domains::BackwardAssignOps<split_oct_domain_t>::assign(*this, x, e, inv);
      }
      
      void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
			  split_oct_domain_t inv) {
	crab::domains::BackwardAssignOps<split_oct_domain_t>::apply(*this, op, x, y, z, inv);
      }
      
      void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
			  split_oct_domain_t inv) {
	crab::domains::BackwardAssignOps<split_oct_domain_t>::apply(*this, op, x, y, z, inv);
      }
      

      void expand(variable_t x, variable_t y) {
	crab::CrabStats::count(getDomainName() + ".count.expand");
	crab::ScopedCrabStats __st__(getDomainName() + ".expand");
		        
	if (is_bottom()) return;
	CRAB_LOG("octagon-split",
		  crab::outs() << "Before expand " << x << " into " << y << ":\n"
		  << *this <<"\n");

	auto it = vert_map.find(variable_t(y));
	if(it != vert_map.end()) {
	  CRAB_ERROR("split_dbm expand operation failed because y already exists");
	}

	vert_id ii = get_vert(x);
	vert_id jj = get_vert(y);

	for (auto edge : graph.e_preds(ii))
	  graph.add_edge(edge.vert, edge.val, jj);

	for (auto edge : graph.e_succs(ii))
	  graph.add_edge(jj, edge.val, edge.vert);

	for (auto edge : graph.e_preds(ii+1))
	  graph.add_edge(edge.vert, edge.val, jj+1);

	for (auto edge : graph.e_succs(ii+1))
	  graph.add_edge(jj+1, edge.val, edge.vert);

	potential[jj] = potential[ii];
	potential[jj+1] = potential[ii+1];

	CRAB_LOG("octagon-split",
		  crab::outs() << "After expand " << x << " into " << y << ":\n"
		  << *this <<"\n");
      }

      void rename(const variable_vector_t &from, const variable_vector_t &to) {
	if (is_top() || is_bottom()) return;
		    	
	CRAB_LOG("octagon-split",
		 crab::outs() << "Replacing {";
		 for (auto v: from) crab::outs() << v << ";";
		 crab::outs() << "} with ";
		 for (auto v: to) crab::outs() << v << ";";
		 crab::outs() << "}:\n";
		 crab::outs() << *this << "\n";);

	vert_map_t new_vert_map;
	for (auto kv : vert_map) {
	  ptrdiff_t pos = std::distance(from.begin(),
		    			std::find(from.begin(), from.end(), kv.first));
	  if (pos < from.size()) {
	    variable_t new_v(to[pos]);
	    new_vert_map.insert(vmap_elt_t(new_v, kv.second));
	    rev_map[kv.second.first] = new_v;
	    rev_map[kv.second.second] = new_v;
	  } else {
	    new_vert_map.insert(kv);
	  }
	}
	std::swap(vert_map, new_vert_map);
	CRAB_LOG("octagon-split", crab::outs() << "RESULT=" << *this << "\n");
      }

      template<typename OtherNumDom>
      void push(const VariableName& x, OtherNumDom&inv) {
	crab::CrabStats::count(getDomainName() + ".count.push");
	crab::ScopedCrabStats __st__(getDomainName() + ".push");
	normalize();
	if (is_bottom() || inv.is_bottom()) return;
	linear_constraint_system_t csts;

	auto it = vert_map.find(x);
	if (it != vert_map.end()) {
	  vert_id s = (*it).second;
	  if (rev_map[s]) {
	    variable_t vs = *rev_map[s];
	    SplitGraph<graph_t> g_excl(graph);
	    for (vert_id d : g_excl.vert()) {
	      if (rev_map[d]) {
		variable_t vd = *rev_map[d];
		if (g_excl.elem(s, d) && g_excl.elem(d, s) &&
		    g_excl.edge_val(s, d) == 0 &&
		    g_excl.edge_val(d, s) == 0) {
		  linear_constraint_t cst(vs == vd);
		  csts += cst;
		  continue;
		}
		
		if (g_excl.elem(s, d)) {
		  linear_constraint_t cst(vd - vs <= g_excl.edge_val(s, d));
		  csts += cst;
		}
		
		if (g_excl.elem(d, s)) {
		  linear_constraint_t cst(vs - vd <= g_excl.edge_val(d, s));
		  csts += cst;
		}
	      }
	    }
	  }
	}
	inv += csts;
      }

      /* Begin unimplemented operations */
      // boolean operations
      void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
      void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
      void apply_binary_bool(bool_operation_t op, variable_t x,variable_t y,variable_t z) {}
      void assume_bool(variable_t v, bool is_negated) {}
      // backward boolean operations
      void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
				    split_oct_domain_t invariant){}
      void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
				      split_oct_domain_t invariant) {}
      void backward_apply_binary_bool(bool_operation_t op,
				      variable_t x,variable_t y,variable_t z,
				      split_oct_domain_t invariant) {}
      // array operations
      void array_init(variable_t a, linear_expression_t elem_size,
		      linear_expression_t lb_idx, linear_expression_t ub_idx, 
		      linear_expression_t val) {}      
      void array_load(variable_t lhs,
		      variable_t a, linear_expression_t elem_size,
		      linear_expression_t i) {}
      void array_store(variable_t a, linear_expression_t elem_size,
		       linear_expression_t i, linear_expression_t v, 
		       bool is_singleton) {}
      void array_store_range(variable_t a, linear_expression_t elem_size,
			     linear_expression_t i, linear_expression_t j,
			     linear_expression_t v) {}                  
      void array_assign(variable_t lhs, variable_t rhs) {}
      // pointer operations
      void pointer_load(variable_t lhs, variable_t rhs)  {}
      void pointer_store(variable_t lhs, variable_t rhs) {} 
      void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
      void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
      void pointer_function(variable_t lhs, varname_t func) {}
      void pointer_mk_null(variable_t lhs) {}
      void pointer_assume(pointer_constraint_t cst) {}
      void pointer_assert(pointer_constraint_t cst) {}
      /* End unimplemented operations */
      
      void write(crab_os& o) {
	typename graph_t::mut_val_ref_t w;
	if (is_bottom()) {
	  o << "_|_";
	  return;
	} else if (is_top()) {
	  o << "{}";
	  return;
	} else {
	  bool first = true;
	  o << "{";
	  SplitGraph<graph_t> g_excl(graph);
	  for (vert_id v : graph.verts()) {
	    // CRAB_LOG("octagon-split", crab::outs() << "writing vertex " << v << "\n");
	    if (v%2 != 0) continue;
	    if (!rev_map[v]) continue;
	    if (!graph.elem(v, v+1) && !graph.elem(v+1, v)) continue;
	    interval_t v_out = interval_t(graph.elem(v, v+1) ?
					  -number_t(graph.edge_val(v, v+1))/2 :
					  bound_t::minus_infinity(),
					  graph.elem(v+1, v) ?
					  number_t(graph.edge_val(v+1, v))/2 :
					  bound_t::plus_infinity());
		    			
	    if (first) first = false;
	    else o << ", ";
	    o << *(rev_map[v]) << " -> " << v_out;
	  }
	  for (vert_id s : graph.verts()) {
	    if (!rev_map[s]) continue;
	    variable_t vs = *(rev_map[s]);
	    for (vert_id d : g_excl.succs(s)) {
	      if (!rev_map[d]) continue;
	      variable_t vd = *(rev_map[d]);
	      if (first) first = false;
	      else o << ", ";

	      if (s % 2 == 0 && d % 2 == 0) {
		o << vd << "-" << vs << "<=" << g_excl.edge_val(s, d);
	      } else if (s % 2 != 0 && d % 2 == 0) {
		o << vd << "+" << vs << "<=" << g_excl.edge_val(s, d);
	      } else if (s % 2 == 0 && d % 2 != 0) {
		o << "-" << vd << "-" << vs << "<=" << g_excl.edge_val(s, d);
	      } else if (s % 2 != 0 && d % 2 != 0) {
		o << "-" << vd << "+" << vs << "<=" << g_excl.edge_val(s, d);
	      }
		    				
	    }
	  }
	  o << '}';
	}
      }

      linear_constraint_system_t to_linear_constraint_system() {
	normalize();
	linear_constraint_system_t csts;
	if (is_bottom()) {
          csts += linear_constraint_t::get_false();	  
	  return csts;
	}
	SplitGraph<graph_t> g_excl(graph);
	for (vert_id v : g_excl.verts()) {
	  if (v%2 != 0) continue;
	  if (!rev_map[v]) continue;
	  if (graph.elem(v, v+1))
	    csts += linear_constraint_t(linear_expression_t(*rev_map[v]) >= -(graph.edge_val(v, v+1)/2));
	  if (graph.elem(v+1, v))
	    csts += linear_constraint_t(linear_expression_t(*rev_map[v]) <= (graph.edge_val(v+1, v)/2));
	}
	for (vert_id s : g_excl.verts()) {
	  if (!rev_map[s]) continue;
	  variable_t vs = *rev_map[s];
	  auto s_exp = linear_expression_t(vs);
	  for (vert_id d : g_excl.succs(s)){
	    if (!rev_map[d]) continue;
	    variable_t vd = *rev_map[d];
	    auto d_exp = linear_expression_t(vd);
	    auto w = g_excl.edge_val(s, d);

	    if (s % 2 == 0 && d % 2 == 0) {
	      csts += linear_constraint_t(d_exp - s_exp <= w);
	    } else if (s % 2 != 0 && d % 2 == 0) {
	      csts += linear_constraint_t(d_exp + s_exp <= w);
	    } else if (s % 2 == 0 && d % 2 != 0) {
	      csts += linear_constraint_t(- d_exp - s_exp <= w);
	    } else if (s % 2 != 0 && d % 2 != 0) {
	      csts += linear_constraint_t(- d_exp - d_exp <= w);
	    }
	  }
	}
	return csts;
      }

      disjunctive_linear_constraint_system_t to_disjunctive_linear_constraint_system() {
	auto lin_csts = to_linear_constraint_system();
	if (lin_csts.is_false()) {
	  return disjunctive_linear_constraint_system_t(true /*is_false*/); 
	} else if (lin_csts.is_true()) {
	  return disjunctive_linear_constraint_system_t(false /*is_false*/);
	} else {
	  return disjunctive_linear_constraint_system_t(lin_csts);
	}
      }

      /* Begin array_sgraph_domain_traits */ 
      // return true iff inequality cst is unsatisfiable.
      bool is_unsat(linear_constraint_t cst) {
	if (is_bottom() || cst.is_contradiction()) return true;
	if (is_top() || cst.is_tautology()) return false;

	if (!cst.is_inequality()) {
	  return false;
	}

	auto diffcst = diffcst_of_leq(cst);
	if (!diffcst) return false;

	// x - y <= k ?
	auto x = (*diffcst).first;
	auto y = (*diffcst).second.first;
	auto k = (*diffcst).second.second;

	typename graph_t::mut_val_ref_t w;
	if (graph.lookup(y, x, &w)){
	  return ((w + k) < 0);
	} else {
	  interval_t intv_x = interval_t::top();
	  interval_t intv_y = interval_t::top();
	  if (x % 2 == 0) {
	    if (graph.elem(x, x+1) || graph.elem(x+1, x)){
	      intv_x = interval_t(graph.elem(x+1, x) ?
				  -number_t(graph.edge_val(x+1, x))/2 :
				  bound_t::minus_infinity(),
				  graph.elem(x, x+1) ?
				  number_t(graph.edge_val(x, x+1))/2 :
				  bound_t::plus_infinity());
	    }
	  } else {
	    if (graph.elem(x, x-1) || graph.elem(x-1, x)){
	      intv_x = interval_t(graph.elem(x, x-1) ?
				  -number_t(graph.edge_val(x, x-1))/2 :
				  bound_t::minus_infinity(),
				  graph.elem(x-1, x) ?
				  number_t(graph.edge_val(x-1, x))/2 :
				  bound_t::plus_infinity());
	    }
	  }

	  if (y % 2 == 0) {
	    if (graph.elem(y, y+1) || graph.elem(y+1, y)){
	      intv_y = interval_t(graph.elem(y+1, y) ?
				  -number_t(graph.edge_val(y+1, y))/2 :
				  bound_t::minus_infinity(),
				  graph.elem(y, y+1) ?
				  number_t(graph.edge_val(y, y+1))/2 :
				  bound_t::plus_infinity());
	    }
	  } else {
	    if (graph.elem(y, y-1) || graph.elem(y-1, y)){
	      intv_y = interval_t(graph.elem(y, y-1) ?
				  -number_t(graph.edge_val(y, y-1))/2 :
				  bound_t::minus_infinity(),
				  graph.elem(y-1, y) ?
				  number_t(graph.edge_val(y-1, y))/2 :
				  bound_t::plus_infinity());
	    }
	  }

	  if (intv_x.is_top() || intv_y.is_top()) {
	    return false;
	  } else {
	    return (!((intv_y - intv_x).lb() <= k));
	  }
	}
      }
		    
      void active_variables(std::vector<variable_t>& out) const {
	out.reserve(graph.size());
	for (auto v: graph.verts()) {
	  if (rev_map[v]) 
	    out.push_back((*(rev_map[v])));
	}
      }
      /* End array_sgraph_domain_traits */
      
      static std::string getDomainName() {
	return "SplitOctagons";
      }

    }; // end class split_oct_domain

    template<typename Number, typename VariableName, typename Params>    
    struct abstract_domain_traits<split_oct_domain<Number, VariableName, Params>> {
      typedef Number number_t;
      typedef VariableName varname_t;       
    };        
    
    template<typename Number, typename VariableName, typename Params>
    struct array_sgraph_domain_helper_traits<split_oct_domain<Number, VariableName, Params>> {
      typedef split_oct_domain<Number, VariableName, Params> soct_domain_t;
      typedef typename soct_domain_t::linear_constraint_t linear_constraint_t;
      typedef ikos::variable<Number, VariableName> variable_t;

      static bool is_unsat(soct_domain_t& inv, linear_constraint_t cst) {
	return inv.is_unsat(cst);
      }

      static void active_variables(soct_domain_t &inv, std::vector<variable_t>& out) {
	inv.active_variables(out);
      }
    };

  } // end namespace domains
} // end namespace crab

#pragma GCC diagnostic pop
