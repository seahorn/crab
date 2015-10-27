/*******************************************************************************
 *
 * A re-engineered implementation of the Difference Bound Matrix domain,
 * which maintains bounds and relations separately.
 *
 * Closure operations based on the paper "Fast and Flexible Difference Constraint
 * Propagation for DPLL(T)" by Cotton and Maler.
 *
 * Graeme Gange (gkgange@unimelb.edu.au)
 ******************************************************************************/

#ifndef SPLIT_DBM_HPP
#define SPLIT_DBM_HPP

// Uncomment for enabling debug information
//#include <crab/common/dbg.hpp>

#include <crab/common/types.hpp>
#include <crab/common/sparse_graph.hpp>
#include <crab/common/graph_ops.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>
#include <crab/domains/division_operators_api.hpp>

#include <boost/optional.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/container/flat_map.hpp>

using namespace boost;
using namespace std;

namespace crab {

  namespace domains {

     using namespace ikos;

    template<class Number, class VariableName>
    class SplitDBM: public writeable,
               public numerical_domain<Number, VariableName >,
               public bitwise_operators<Number,VariableName >,
               public division_operators<Number, VariableName >{
     public:
      using typename numerical_domain< Number, VariableName >::linear_expression_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_t;
      using typename numerical_domain< Number, VariableName >::linear_constraint_system_t;
      using typename numerical_domain< Number, VariableName >::variable_t;
      using typename numerical_domain< Number, VariableName >::number_t;
      using typename numerical_domain< Number, VariableName >::varname_t;
      
      typedef typename linear_constraint_t::kind_t constraint_kind_t;
      typedef interval<Number>  interval_t;

     private:
      typedef bound<Number>  bound_t;
      // Can't use separate_domain directly, as we need to
      // retrofit some operations onto the join.
      typedef patricia_tree< VariableName, interval_t > ranges_t;
      typedef typename ranges_t::key_binary_op_t key_binary_op_t;
      
      // Eventually break this out into a template param
      typedef Number Wt;
      typedef SparseWtGraph<Wt> graph_t;
      typedef typename graph_t::vert_id vert_id;
      typedef boost::container::flat_map<variable_t, vert_id> vert_map_t;
      typedef typename vert_map_t::value_type vmap_elt_t;
      typedef vector< boost::optional<variable_t> > rev_map_t;

      typedef SplitDBM<Number, VariableName> DBM_t;

      typedef GraphOps<Wt> GrOps;
      typedef GraphPerm<graph_t> GrPerm;
      typedef typename GrOps::edge_vector edge_vector;
      // < <x, y>, k> == x - y <= k.
      typedef pair< pair<VariableName, VariableName>, Wt > diffcst_t;

      protected:
        
      //================
      // Domain data
      //================
      ranges_t ranges; // Intervals for each variable
      vert_map_t vert_map; // Mapping from variables to vertices
      rev_map_t rev_map;
      graph_t g; // The underlying relation graph
      vector<Wt> potential; // Stored potential for the vertex

      bool _is_bottom;

   public:
      SplitDBM(bool is_bottom = false):
        writeable(), _is_bottom(is_bottom)
      {
        g.growTo(1);  // Allocate the zero vector
        potential.push_back(Wt(0));
        rev_map.push_back(none);
      }

      // FIXME: Rewrite to avoid copying if o is _|_
      SplitDBM(const DBM_t& o)
        : writeable(),
          numerical_domain<Number, VariableName >(),
          bitwise_operators< Number, VariableName >(),
          division_operators< Number, VariableName >(),
          ranges(o.ranges), 
          vert_map(o.vert_map),
          rev_map(o.rev_map),
          g(o.g),
          potential(o.potential),
          _is_bottom(false)
      {
        if(o._is_bottom)
          set_to_bottom();

        /*
        for(vert_id v : g.verts())
          if(g.succs(v).begin() != g.succs(v).end() || g.preds(v).begin() != g.preds(v).end())
            if(v != 0)
              assert(rev_map[v]);
              */
      }

      // We should probably use the magical rvalue ownership semantics stuff.
      SplitDBM(ranges_t& _ranges, vert_map_t& _vert_map, rev_map_t& _rev_map, graph_t& _g, vector<Wt>& _potential)
        : writeable(),
          numerical_domain<Number, VariableName >(),
          bitwise_operators< Number, VariableName >(),
          division_operators< Number, VariableName >(),
          ranges(_ranges), vert_map(_vert_map), rev_map(_rev_map), g(_g), potential(_potential),
          _is_bottom(false)
      { }

      // FIXME: Add a move constructor
      SplitDBM& operator=(const SplitDBM& o)
      {
        if(this != &o)
        {
          if(o._is_bottom)
            set_to_bottom();
          else {
            _is_bottom = false;
            ranges = o.ranges;
            vert_map = o.vert_map;
            rev_map = o.rev_map;
            g = o.g;
            potential = o.potential;
          }
        }
        return *this;
      }
       
     private:

      /*
      void forget(vector<int> idxs) {
        dbm ret = NULL;
        ret = dbm_forget_array(&idxs[0], idxs.size(), _dbm);
        dbm_dealloc(_dbm);
        swap(_dbm, ret);
      }
      */

      void set_to_bottom() {
        ranges.clear();
        vert_map.clear();
        rev_map.clear();
        g.clear();
        potential.clear();
        _is_bottom = true;
      }

      // check satisfiability of cst using intervals
      // Only to be used if cst is too hard for dbm
      bool intervals_check_sat(linear_constraint_t cst)  {
        if (is_top())    return true;
        if (is_bottom()) return false;

        auto vars = cst.variables();
        ranges_t inv;
        for(auto v: vars)
          inv.set (v.name(), operator[](v.name));
        inv += cst;
        return !inv.is_bottom();
      }

     public:

      static DBM_t top() { return SplitDBM(false); }
    
      static DBM_t bottom() { return SplitDBM(true); }
    
     public:

      bool is_bottom() {
//        if(!_is_bottom && g.has_negative_cycle())
//          _is_bottom = true;
        return _is_bottom;
      }
    
      bool is_top() {
        if(_is_bottom)
          return false;
        if(ranges.size() != 0)
         return false;
        // GKG: Come up with a cheaper approach for this
        for(auto p : vert_map)
          if(g.succs(p.second).size() != 0)
            return false;

        return true;
      }
    
      bool operator<=(DBM_t& o)  {
        // cover all trivial cases to avoid allocating a dbm matrix
        if (is_bottom()) 
          return true;
        else if(o.is_bottom())
          return false;
        else if (o.is_top ())
          return true;
        else if (is_top ())
          return false;
        else {
          interval_po po;
          if(!ranges.leq(o.ranges, po))
            return false;

          // Set up a mapping from o to this.
          vector<unsigned int> vert_renaming(o.g.size());
          for(auto p : o.vert_map)
          {
            auto it = vert_map.find(p.first);
            // We can't have this <= o if we're missing some
            // vertex.
            if(it == vert_map.end())
              return false;
            vert_renaming[p.second] = (*it).second;
          }

          assert(g.size() > 0);
          GrPerm g_perm(vert_renaming, g);
          // FIXME: Check incorporation of bounds
          //        (and come up with a more efficient method)
          for(vert_id ox : o.g.verts())
          {
            vert_id x = vert_renaming[ox];
            for(vert_id oy : o.g.succs(ox))
            {
              vert_id y = vert_renaming[oy];
              Wt ow = o.g.edge_val(ox, oy);

              // Is the edge implied by ranges?
              if(rev_map[x] && rev_map[y])
              {
                if(get_interval(ranges, *rev_map[x]).ub()
                    - get_interval(ranges, *rev_map[y]).lb() <= bound_t(ow))
                  continue;
              }
              // Edge not present
              if(!g_perm.elem(x, y))
                return false;
              if(!(g_perm.edge_val(x, y) <= ow))
                return false;
            }
          }
          return true;
        }
      }

      class Wt_max {
      public:
       Wt_max() { } 
       Wt apply(const Wt& x, const Wt& y) { return max(x, y); }
       bool default_is_absorbing() { return true; }
      };

      class Wt_min {
      public:
        Wt_min() { }
        Wt apply(const Wt& x, const Wt& y) { return min(x, y); }
        bool default_is_absorbing() { return false; }
      };

      // Perform a join on intervals, but record which
      // lower/upper bounds move in each direction.
      class interval_join_t : public key_binary_op_t {
      public:
        interval_join_t(
            vector<VariableName>& _lb_up,
            vector<VariableName>& _lb_down,
            vector<VariableName>& _ub_up,
            vector<VariableName>& _ub_down)
          : lb_up(_lb_up), lb_down(_lb_down),
            ub_up(_ub_up), ub_down(_ub_down)
        { }

        boost::optional<interval_t> apply(
            VariableName v, interval_t x, interval_t y)
        {
          interval_t z = x.operator|(y);
          if (z.is_top()) {
            return boost::optional< interval_t >();
          } else {
            if(z.lb().is_finite())
            {
              if(x.lb() < y.lb())
                lb_up.push_back(v);
              if(y.lb() < x.lb())
                lb_down.push_back(v);
            }
            if(z.ub().is_finite())
            {
              if(x.ub() < y.ub())
                ub_up.push_back(v);
              if(y.ub() < x.ub())
                ub_down.push_back(v);
            }
            return boost::optional< interval_t >(z);
          }
        }

        bool default_is_absorbing() { return true; }
        vector<VariableName>& lb_up;
        vector<VariableName>& lb_down;
        vector<VariableName>& ub_up;
        vector<VariableName>& ub_down;
      }; // class binary_op

      typedef typename separate_domain<VariableName, interval_t>::meet_op interval_meet_t;
      typedef typename separate_domain<VariableName, interval_t>::widening_op interval_widen_t;
      typedef typename separate_domain<VariableName, interval_t>::domain_po interval_po;
      typedef typename separate_domain<VariableName, interval_t>::bottom_found bottom_found;

      vert_id get_vert(VariableName v)
      {
        auto it = vert_map.find(variable_t(v));
        if(it != vert_map.end())
          return (*it).second;

        Wt v_pot(pot_value(v));
        vert_id vert(g.new_vertex());
        vert_map.insert(vmap_elt_t(variable_t(v), vert)); 
        // Initialize 
        assert(vert <= rev_map.size());
        if(vert < rev_map.size())
        {
          potential[vert] = v_pot;
          rev_map[vert] = v;
        } else {
          potential.push_back(v_pot);
          rev_map.push_back(variable_t(v));
        }
        vert_map.insert(vmap_elt_t(v, vert));
        return vert;
      }

      vert_id get_vert(graph_t& g, vert_map_t& vmap, rev_map_t& rmap,
          ranges_t& ranges, vector<Wt>& pot, VariableName v)
      {
        auto it = vmap.find(variable_t(v));
        if(it != vmap.end())
          return (*it).second;

        Wt v_pot(pot_value(v, ranges, pot));
        vert_id vert(g.new_vertex());
        vmap.insert(vmap_elt_t(variable_t(v), vert)); 
        // Initialize 
        assert(vert <= rmap.size());
        if(vert < rmap.size())
        {
          pot[vert] = v_pot;
          rmap[vert] = v;
        } else {
          pot.push_back(v_pot);
          rmap.push_back(variable_t(v));
        }
        vmap.insert(vmap_elt_t(v, vert));
        return vert;
      }

      DBM_t operator|(DBM_t o) {
        if (is_bottom() || o.is_top ())
          return o;
        else if (is_top () || o.is_bottom())
          return *this;
        else {
          CRAB_DEBUG ("Before join:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

          // First, join the intervals, collecting change directions
          vector<VariableName> lb_up;
          vector<VariableName> lb_down;
          vector<VariableName> ub_up;
          vector<VariableName> ub_down;
          interval_join_t join_op(lb_up, lb_down, ub_up, ub_down);
          
          ranges_t join_range(ranges);
          join_range.merge_with(o.ranges, join_op);

          // Figure out the common renaming, initializing the
          // resulting potentials as we go.
          vector<vert_id> perm_x;
          vector<vert_id> perm_y;
          vector<variable_t> perm_inv;

          vector<Wt> join_pot;
          vert_map_t out_vmap;
          rev_map_t out_revmap;
          // Add the zero vertex
          assert(potential.size() > 0);
          join_pot.push_back(potential[0]);
          perm_x.push_back(0);
          perm_y.push_back(0);
          out_revmap.push_back(none);
          for(auto p : vert_map)
          {
            auto it = o.vert_map.find(p.first); 
            // Variable exists in both
            if(it != o.vert_map.end())
            {
              out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
              out_revmap.push_back(p.first);

              join_pot.push_back(potential[p.second]);
              perm_inv.push_back(p.first);
              perm_x.push_back(p.second);
              perm_y.push_back((*it).second);
            }
          }
          unsigned int sz = perm_x.size();

          // Build the permuted view of x and y.
          assert(g.size() > 0);
          GrPerm gx(perm_x, g);
          assert(o.g.size() > 0);
          GrPerm gy(perm_y, o.g);

          // Compute the deferred relations
          graph_t g_ix_ry;
          g_ix_ry.growTo(sz);
          SubGraph<GrPerm> gy_excl(gy, 0);
          for(vert_id s : gy_excl.verts())
          {
            for(vert_id d : gy_excl.succs(s))
            {
              // Assumption: gx.elem(s, d) -> gx.edge_val(s, d) <= ranges[var(s)].ub() - ranges[var(d)].lb()
              // That is, if the relation exists, it's at least as strong as the bounds.
              if(!gx.elem(s, d))
              {
                // Check the bounds implied by o.ranges
                bound_t b = get_interval(ranges,perm_inv[d]).ub() - get_interval(ranges,perm_inv[s]).lb();
                if(b.is_finite())
                {
                  g_ix_ry.add_edge(s, *(b.number()), d);
                }
              }
            }
          }
          // Apply the deferred relations, and re-close.
          edge_vector delta;
          graph_t g_rx(GrOps::meet(gx, g_ix_ry));
          SubGraph<graph_t> g_rx_excl(g_rx, 0);
          GrOps::close_after_meet(g_rx_excl, potential, gx, g_ix_ry, delta);
          GrOps::apply_delta(g_rx, delta);

          graph_t g_rx_iy;
          g_rx_iy.growTo(sz);

          SubGraph<GrPerm> gx_excl(gx, 0);
          for(vert_id s : gx_excl.verts())
          {
            for(vert_id d : gx_excl.succs(s))
            {
              // Assumption: gx.mem(s, d) -> gx.edge_val(s, d) <= ranges[var(s)].ub() - ranges[var(d)].lb()
              // That is, if the relation exists, it's at least as strong as the bounds.
              if(!gy.elem(s, d))
              {
                // Check the bounds implied by o.ranges
                bound_t b = get_interval(o.ranges,perm_inv[d]).ub() - get_interval(o.ranges,perm_inv[s]).lb();
                if(b.is_finite())
                {
                  g_rx_iy.add_edge(s, *(b.number()), d);
                }
              }
            }
          }
          delta.clear();
          // Similarly, should use a SubGraph view.
          graph_t g_ry(GrOps::meet(gy, g_rx_iy));
          SubGraph<graph_t> g_ry_excl(g_ry, 0);
          GrOps::close_after_meet(g_ry_excl, o.potential, gy, g_rx_iy, delta);
          GrOps::apply_delta(g_ry, delta);
           
          // We now have the relevant set of relations. Because g_rx and g_ry are closed,
          // the result is also closed.
          Wt_min min_op;
          graph_t join_g(GrOps::join(g_rx, g_ry));
          // Now reapply the missing independent relations.
          // Need to derive vert_ids from lb_up/lb_down, and make sure the vertices exist
          for(VariableName vs : lb_up)
          {
            bound_t ls_0 = get_interval(ranges, vs).lb();  
            bound_t ls_1 = get_interval(o.ranges, vs).lb();
            for(VariableName vd : ub_up)
            {
              if(vs == vd)
                continue;

              bound_t ud_0 = get_interval(ranges, vd).ub();
              bound_t ud_1 = get_interval(o.ranges, vd).ub();

              bound_t b = max(ud_0 - ls_0, ud_1 - ls_1);

              vert_id s(get_vert(join_g, out_vmap, out_revmap, join_range, join_pot, vs));
              vert_id d(get_vert(join_g, out_vmap, out_revmap, join_range, join_pot, vd));
              join_g.update_edge(s, *(b.number()), d, min_op);
            }
          }

          for(VariableName vs : lb_down)
          {
            bound_t ls_0 = get_interval(ranges, vs).lb();  
            bound_t ls_1 = get_interval(o.ranges, vs).lb();
            for(VariableName vd : ub_down)
            {
              if(vs == vd)
                continue;

              bound_t ud_0 = get_interval(ranges, vd).ub();
              bound_t ud_1 = get_interval(o.ranges, vd).ub();

              bound_t b = max(ud_0 - ls_0, ud_1 - ls_1);

              vert_id s(get_vert(join_g, out_vmap, out_revmap, join_range, join_pot, vs));
              vert_id d(get_vert(join_g, out_vmap, out_revmap, join_range, join_pot, vd));
              join_g.update_edge(s, *(b.number()), d, min_op);
            }
          }

          // Conjecture: join_g remains closed.
          
          DBM_t res(join_range, out_vmap, out_revmap, join_g, join_pot);
          CRAB_DEBUG ("Result join:\n",res);
           
          return res;
        }
      }

      DBM_t operator||(DBM_t o) {	
        if (is_bottom())
          return o;
        else if (o.is_bottom())
          return *this;
        else {
          CRAB_DEBUG ("Before widening:\n","DBM 1\n",*this,"\n","DBM 2\n",o);
          interval_widen_t widen_op;
          ranges_t widen_range(ranges);
          widen_range.merge_with(o.ranges, widen_op);
          
          // Figure out the common renaming
          vector<vert_id> perm_x;
          vector<vert_id> perm_y;
          vert_map_t out_vmap;
          rev_map_t out_revmap;
          vector<Wt> widen_pot;
          assert(potential.size() > 0);
          widen_pot.push_back(potential[0]);
          perm_x.push_back(0);
          perm_y.push_back(0);
          out_revmap.push_back(none);
          for(auto p : vert_map)
          {
            auto it = o.vert_map.find(p.first); 
            // Variable exists in both
            if(it != o.vert_map.end())
            {
              out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
              out_revmap.push_back(p.first);

              widen_pot.push_back(potential[p.second]);
              perm_x.push_back(p.second);
              perm_y.push_back((*it).second);
            }
          }

          // Build the permuted view of x and y.
          assert(g.size() > 0);
          GrPerm gx(perm_x, g);            
          assert(o.g.size() > 0);
          GrPerm gy(perm_y, o.g);
         
          // Now perform the widening 
          graph_t widen_g(GrOps::widen(gx, gy));
          DBM_t res(widen_range, out_vmap, out_revmap, widen_g, widen_pot);

          // GKG: need to mark changes so we can restore closure
           
          CRAB_DEBUG ("Result widening:\n",res);
          return res;
        }
      }

      DBM_t operator&(DBM_t o) {
        if (is_bottom() || o.is_bottom())
          return bottom();
        else if (is_top())
          return o;
        else if (o.is_top())
          return *this;
        else{
          CRAB_DEBUG ("Before meet:\n","DBM 1\n",*this,"\n","DBM 2\n",o);
          interval_meet_t meet_op;

          ranges_t meet_range(ranges);
          try {
            meet_range.merge_with(o.ranges, meet_op);
          } 
          catch (bottom_found& exc) {
            return bottom();
          }

          // We map vertices in the left operand onto a contiguous range.
          // This will often be the identity map, but there might be gaps.
          vert_map_t meet_verts;
          rev_map_t meet_rev;

          vector<vert_id> perm_x;
          vector<vert_id> perm_y;
          perm_x.push_back(0);
          perm_y.push_back(0);
          meet_rev.push_back(none);
          for(auto p : vert_map)
          {
            vert_id vv = perm_x.size();
            meet_verts.insert(vmap_elt_t(p.first, vv));
            meet_rev.push_back(p.first);

            perm_x.push_back(p.second);
            perm_y.push_back(-1);
          }

          // Add missing mappings from the right operand.
          for(auto p : o.vert_map)
          {
            auto it = meet_verts.find(p.first);

            if(it == meet_verts.end())
            {
              vert_id vv = perm_y.size();
              meet_rev.push_back(p.first);

              perm_y.push_back(p.second);
              perm_x.push_back(-1);
              meet_verts.insert(vmap_elt_t(p.first, vv));
            } else {
              perm_y[(*it).second] = p.second;
            }
          }

          // Build the permuted view of x and y.
          assert(g.size() > 0);
          GrPerm gx(perm_x, g);
          assert(o.g.size() > 0);
          GrPerm gy(perm_y, o.g);

          // Compute the syntactic meet of the permuted graphs.
          graph_t meet_g(GrOps::meet(gx, gy));
           
          // Compute updated potentials on the zero-enriched graph
          vector<Wt> meet_pi(meet_g.size());
          if(!GrOps::select_potentials(meet_g, meet_pi))
          {
            // Potentials cannot be selected -- state is infeasible.
            return bottom();
          }

          edge_vector delta;
          SubGraph<graph_t> meet_g_excl(meet_g, 0);
          GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
          GrOps::apply_delta(meet_g, delta);

          for(auto e : delta)
          {
            vert_id s = e.first.first;
            vert_id d = e.first.second; 
            Wt w = e.second;

            // d - s <= w --> d <= w + s, s >= d - w
            // Potentially update ub(d) and lb(s)
            /*
             * // FIXME: Again, map vertices back to variables
            interval_t b_s0 = get_interval(meet_range,s);
            interval_t b_d0 = get_interval(meet_range,d);
            if(b_s0.ub().is_finite())
            {
              bound_t ub_d = bound_t(w) + b_s0.ub();
              if(ub_d <= b_d0.ub())
                meet_range.insert(d, interval_t(b_d0.lb(), ub_d));
            }

            if(b_d0.lb().is_finite())
            {
              bound_t lb_s = bound_t(w) - b_d0.lb();
              if(b_s0.lb() <= lb_s)
                meet_range.insert(d, interval_t(lb_s, b_s0.ub()));
            }
            */
          }
          DBM_t res(meet_range, meet_verts, meet_rev, meet_g, meet_pi);
          CRAB_DEBUG ("Result meet:\n",res);
          return res;
        }
      }
    
      DBM_t operator&&(DBM_t o) {
        if (is_bottom() || o.is_bottom())
          return bottom();
        else if (is_top ())
          return o;
        else{
          CRAB_DEBUG ("Before narrowing:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

          // FIXME: Implement properly
          // Narrowing as a no-op should be sound.
          DBM_t res(*this);

          /*
          id_t id = 0;
          // Build new map
          var_map_t var_map = merge_var_map (_var_map, _dbm, 
                                             o._var_map, o._dbm,
                                             id, subs_x, subs_y);
          // Build the reverse map
          rev_map_t rev_map;
          for(auto p: var_map)
            rev_map.insert (make_pair (p.second, p.first));

          */
           

          CRAB_DEBUG ("Result narrowing:\n",res);
          return res;
        }
      }	

      void normalize() {
        // dbm_canonical(_dbm);
        // Always maintained in normal form, except for widening
        // FIXME: Handle closing after widening
      }

      void operator-=(VariableName v) {
        if (is_bottom ())
          return;

        ranges.remove(v);
        auto it = vert_map.find (v);
        if (it != vert_map.end ()) {
          CRAB_DEBUG("Before forget ", it->second, ": ", g);
          g.forget(it->second);
          CRAB_DEBUG("After: ", g);
          rev_map[it->second] = boost::none;
          vert_map.erase(v);
        }
      }

      template<typename Iterator>
      void forget (Iterator vIt, Iterator vEt) {
        if (is_bottom ())
          return;
        CRAB_WARN("forget not implemented.");
        /*
        vector<int> idxs;
        for (auto v: boost::make_iterator_range (vIt,vEt)) {
          auto it = vert_map.find (v);
          if (it != vert_map.end ()) {
            idxs.push_back  ((int) it->second);
            _var_map.erase (v);
            _rev_map.erase (it->second);
          }
        }
        if (!idxs.empty()) {
          // --- forget all at once is more efficient than calling
          //     operator-= multiple times.
          forget (idxs);
          if (!dbm_compact ())
            resize (Delta);
        }
        */
      }

      // Evaluate the potential value of a variable.
      Wt pot_value(variable_t v)
      {
        auto it = vert_map.find(v); 
        if(it != vert_map.end())
          return potential[(*it).second];

        interval_t r(get_interval(ranges,v));
        if(r.lb().is_finite())
          return (Wt) (*(r.lb().number()));
        if(r.ub().is_finite())
          return (Wt) (*(r.ub().number()));
        return ((Wt) 0);
      }

      Wt pot_value(variable_t v, ranges_t& ranges, vector<Wt>& potential)
      {
        auto it = vert_map.find(v); 
        if(it != vert_map.end())
          return potential[(*it).second];

        interval_t r(get_interval(ranges,v));
        if(r.lb().is_finite())
          return (Wt) (*(r.lb().number()));
        if(r.ub().is_finite())
          return (Wt) (*(r.ub().number()));
        return ((Wt) 0);
      }

      // Evaluate an expression under the chosen potentials
      Wt eval_expression(linear_expression_t e)
      {
        Wt v((Wt) e.constant()); 
        for(auto p : e)
        {
          v += pot_value(p.second)*((Wt) p.first);
        }
        return v;
      }
      
      interval_t eval_interval(linear_expression_t e)
      {
        interval_t r = e.constant();
        for (auto p : e)
          r += p.first * operator[](p.second.name());

        return r;
      }

      // Turn an assignment into a set of difference constraints.
      void diffcsts_of_assign(VariableName x, linear_expression_t exp,
          vector<pair<VariableName, Wt> >& lb, vector<pair<VariableName,Wt> >& ub)
      {
        {
          // Process upper bounds.
          optional<VariableName> unbounded_ubvar;
          Wt exp_ub = exp.constant();
          vector< pair<VariableName, Wt> > ub_terms;
          for(auto p : exp)
          {
            Wt coeff = p.first;
            if(p.first < Wt(0))
            {
              // Can't do anything with negative coefficients.
              bound_t y_lb = operator[](p.second.name()).lb();
              if(y_lb.is_infinite())
                goto assign_ub_finish;
              exp_ub += (*(y_lb.number()))*coeff;
            } else {
              VariableName y(p.second.name());
              bound_t y_ub = operator[](y).ub(); 
              if(y_ub.is_infinite())
              {
                if(unbounded_ubvar || coeff != Wt(1))
                  goto assign_ub_finish;
                unbounded_ubvar = y;
              } else {
                Wt ymax(*(y_ub.number()));
                exp_ub += ymax*coeff;
                ub_terms.push_back(make_pair(y, ymax));
              }
            }
          }

          if(unbounded_ubvar)
          {
            // There is exactly one unbounded variable. 
            ub.push_back(make_pair(*unbounded_ubvar, exp_ub));
          } else {
            for(auto p : ub_terms)
            {
              ub.push_back(make_pair(p.first, exp_ub - p.second));
            }
          }
        }
      assign_ub_finish:

        {
          optional<VariableName> unbounded_lbvar;
          Wt exp_lb = exp.constant();
          vector< pair<VariableName, Wt> > lb_terms;
          for(auto p : exp)
          {
            Wt coeff = p.first;
            if(p.first < Wt(0))
            {
              // Again, can't do anything with negative coefficients.
              bound_t y_ub = operator[](p.second.name()).ub();
              if(y_ub.is_infinite())
                goto assign_lb_finish;
              exp_lb += (*(y_ub.number()))*coeff;
            } else {
              VariableName y(p.second.name());
              bound_t y_lb = operator[](y).lb(); 
              if(y_lb.is_infinite())
              {
                if(unbounded_lbvar || coeff != Wt(1))
                  goto assign_lb_finish;
                unbounded_lbvar = y;
              } else {
                Wt ymin(*(y_lb.number()));
                exp_lb += ymin*coeff;
                lb_terms.push_back(make_pair(y, ymin));
              }
            }
          }

          if(unbounded_lbvar)
          {
            lb.push_back(make_pair(*unbounded_lbvar, -exp_lb));
          } else {
            for(auto p : lb_terms)
            {
              lb.push_back(make_pair(p.first, p.second - exp_lb));
            }
          }
        }
    assign_lb_finish:
        return;
      }
   
      // GKG: I suspect there're some sign/bound direction errors in the 
      // following.
      void diffcsts_of_lin_leq(const linear_expression_t& exp, vector<diffcst_t>& csts,
          vector<pair<VariableName, Wt> >& lbs, vector<pair<VariableName, Wt> >& ubs)
      {
        // Process upper bounds.
        Wt unbounded_lbcoeff;
        Wt unbounded_ubcoeff;
        optional<VariableName> unbounded_lbvar;
        optional<VariableName> unbounded_ubvar;
        Wt exp_ub = -exp.constant();
        vector< pair< pair<Wt, VariableName>, Wt> > pos_terms;
        vector< pair< pair<Wt, VariableName>, Wt> > neg_terms;
        for(auto p : exp)
        {
          Wt coeff = p.first;
          if(coeff > Wt(0))
          {
            VariableName y(p.second.name());
            bound_t y_lb = operator[](y).lb();
            if(y_lb.is_infinite())
            {
              if(unbounded_lbvar);
                goto diffcst_finish;
              unbounded_lbvar = y;
              unbounded_lbcoeff = coeff;
            } else {
              Wt ymin = (*(y_lb.number()));
              // Coeff is negative, so it's still add
              exp_ub += ymin*coeff;
              pos_terms.push_back(make_pair(make_pair(coeff, y), ymin));
            }
          } else {
            VariableName y(p.second.name());
            bound_t y_ub = operator[](y).ub(); 
            if(y_ub.is_infinite())
            {
              if(unbounded_ubvar)
                goto diffcst_finish;
              unbounded_ubvar = y;
              unbounded_ubcoeff = -coeff;
            } else {
              Wt ymax(*(y_ub.number()));
              exp_ub += ymax*coeff;
              neg_terms.push_back(make_pair(make_pair(-coeff, y), ymax));
            }
          }
        }

        if(unbounded_lbvar)
        {
          VariableName x(*unbounded_lbvar);
          if(unbounded_ubvar)
          {
            if(unbounded_lbcoeff != Wt(1) || unbounded_ubcoeff != Wt(1))
              goto diffcst_finish;
            VariableName y(*unbounded_ubvar);
            csts.push_back(make_pair(make_pair(x, y), exp_ub));
          } else {
            if(unbounded_lbcoeff == Wt(1))
            {
              for(auto p : pos_terms)
                csts.push_back(make_pair(make_pair(x, p.first.second), exp_ub - p.second));
            }
            // Add bounds for x
            ubs.push_back(make_pair(x, exp_ub/unbounded_lbcoeff));
          }
        } else {
          if(unbounded_ubvar)
          {
            VariableName y(*unbounded_ubvar);
            if(unbounded_ubcoeff == Wt(1))
            {
              for(auto p : neg_terms)
                csts.push_back(make_pair(make_pair(p.first.second, y), exp_ub + p.second));
            }
            // Bounds for y
            lbs.push_back(make_pair(y, -exp_ub/unbounded_ubcoeff));
          } else {
            for(auto pl : neg_terms)
              for(auto pu : pos_terms)
                csts.push_back(make_pair(make_pair(pl.first.second, pu.first.second), exp_ub + pl.second - pu.second));
            for(auto pl : neg_terms)
              lbs.push_back(make_pair(pl.first.second, -pl.second - exp_ub/pl.first.first));
            for(auto pu : pos_terms)
              ubs.push_back(make_pair(pu.first.second, exp_ub/pu.first.first - pu.second));
          }
        }
    diffcst_finish:
        return;
      }

      // Assumption: state is currently feasible.
      void assign(VariableName x, linear_expression_t e) {

        if(is_bottom())
          return;
        CRAB_DEBUG(x, ":=", e);

        // If it's a constant, just assign the interval.
        if (e.is_constant()){
          set(x, e.constant());
        } else {
          interval_t x_int = eval_interval(e);
          vector<pair<VariableName, Wt> > diffs_lb;
          vector<pair<VariableName, Wt> > diffs_ub;
          // Construct difference constraints from the assignment
          diffcsts_of_assign(x, e, diffs_lb, diffs_ub);
          if(diffs_lb.size() > 0 || diffs_ub.size() > 0)
          {
            // Allocate a new vertex for x
            vert_id v = g.new_vertex();
            assert(v <= rev_map.size());
            if(v == rev_map.size())
            {
              rev_map.push_back(variable_t(x));
              potential.push_back(eval_expression(e));
            } else {
              potential[v] = eval_expression(e);
              rev_map[v] = x;
            }
            
            edge_vector delta;
            for(auto diff : diffs_lb)
            {
              delta.push_back(make_pair(make_pair(get_vert(diff.first), v), diff.second));
            }

            for(auto diff : diffs_ub)
            {
              delta.push_back(make_pair(make_pair(v, get_vert(diff.first)), diff.second));
            }
               
            for(auto diff : delta)
            {
              vert_id s = diff.first.first;
              vert_id d = diff.first.second;

              CRAB_DEBUG("|- ", (*rev_map[s]).name(), "-", (*rev_map[d]).name(), "<=", diff.second);
            }
            GrOps::apply_delta(g, delta);
            delta.clear();
            SubGraph<graph_t> g_excl(g, 0);
            GrOps::close_after_assign(g_excl, potential, v, delta);
            GrOps::apply_delta(g, delta);

            // Clear the old x vertex
            operator-=(x);
            ranges.insert(x, x_int);
            vert_map.insert(vmap_elt_t(variable_t(x), v));
          } else {
            set(x, x_int);
          }
          // CRAB_WARN("DBM only supports a cst or var on the rhs of assignment");
          // this->operator-=(x);
        }

        // Only need to update potentials and paths for x
            

        CRAB_DEBUG("---", x, ":=", e,"\n",*this);
      }

      void apply(operation_t op, VariableName x, VariableName y, VariableName z){	
        if(is_bottom())
          return;

        switch(op)
        {
          case OP_ADDITION:
          {
            linear_expression_t e(linear_expression_t(y) + linear_expression_t(z));
            assign(x, e);
            break;
          }
          case OP_SUBTRACTION:
          {
            linear_expression_t e(linear_expression_t(y) - linear_expression_t(z));
            assign(x, e);
            break;
          }
          // For mul and div, we fall back on intervals.
          case OP_MULTIPLICATION:
          {
            set(x, get_interval(ranges, y)*get_interval(ranges, z));
            break;
          }
          case OP_DIVISION:
          {
            interval_t xi(get_interval(ranges, y)/get_interval(ranges, z));
            if(xi.is_bottom())
              set_to_bottom();
            else
              set(x, xi);
            break;
          }
        }
        /*
        if (x == y){
          // --- ensure lhs does not appear on the rhs
          assign_tmp(y); 
          apply(op, x, get_tmp(), get_dbm_index(z));
          forget(get_tmp());
        }
        else if (x == z){
          // --- ensure lhs does not appear on the rhs
          assign_tmp(z); 
          apply(op, x, get_dbm_index (y), get_tmp());
          forget(get_tmp());
        }
        else{
          if (x == y && y == z)
            CRAB_ERROR("DBM: does not support x := x + x ");
          else
            apply(op, x, get_dbm_index(y), get_dbm_index(z));
        }
      */
        CRAB_DEBUG("---", x, ":=", y, op, z,"\n", *this);
      }

    
      void apply(operation_t op, VariableName x, VariableName y, Number k) {	
        if(is_bottom())
          return;

        // FIXME: Implement
        switch(op)
        {
          case OP_ADDITION:
          {
            linear_expression_t e(linear_expression_t(y) + linear_expression_t(k));
            assign(x, e);
            break;
          }
          case OP_SUBTRACTION:
          {
            linear_expression_t e(linear_expression_t(y) - linear_expression_t(k));
            assign(x, e);
            break;
          }
          // For mul and div, we fall back on intervals.
          case OP_MULTIPLICATION:
          {
            set(x, get_interval(ranges, y)*k);

            break;
          }
          case OP_DIVISION:
          {
            // FIXME: Implement
            if(k == Wt(0))
              set_to_bottom();
            else
              set(x, get_interval(ranges, y)/k);

            break;
          }
        }

        CRAB_DEBUG("---", x, ":=", y, op, k,"\n", *this);
      }
      
      bool add_linear_leq(const linear_expression_t& exp)
      {
        CRAB_DEBUG("Adding: ", exp);
        vector< pair<VariableName, Wt> > lbs;
        vector< pair<VariableName, Wt> > ubs;
        vector<diffcst_t> csts;
        diffcsts_of_lin_leq(exp, csts, lbs, ubs);

        for(auto p : lbs)
        {
          CRAB_DEBUG(p.first, ">=", p.second);
          VariableName x(p.first);
          auto it = vert_map.find(x);
          // Vertex already exists; add the edge and repair
          if(it != vert_map.end())
          {
            vert_id v = (*it).second;
            g.add_edge(v, -p.second, 0); 
            if(!repair_potential(v, 0))
            {
              set_to_bottom();
              return false;
            }
          }
          // Apply the interval.
          // GKG: Check for bottom.
          interval_t x_out = get_interval(ranges, x)&
            interval_t(p.second, bound_t::plus_infinity());
          if(x_out.is_bottom())
          {
            set_to_bottom(); return false;
          }
          ranges.insert(x, x_out);
        }
        for(auto p : ubs)
        {
          CRAB_DEBUG(p.first, "<=", p.second);
          VariableName x(p.first);
          auto it = vert_map.find(x);
          if(it != vert_map.end())
          {
            vert_id v = (*it).second;
            g.add_edge(0, p.second, v);
            if(!repair_potential(0, v))
            {
              set_to_bottom();
              return false;
            }
          }
          // Apply the interval.
          interval_t x_out = get_interval(ranges, x)&
              interval_t(bound_t::minus_infinity(), p.second);
          if(x_out.is_bottom())
          {
            set_to_bottom(); return false;
          }
          ranges.insert(x, x_out);
        }

        for(auto diff : csts)
        {
          CRAB_DEBUG(diff.first.first, "-", diff.first.second, "<=", diff.second);

          vert_id src = get_vert(diff.first.second);
          vert_id dest = get_vert(diff.first.first);
          g.add_edge(src, diff.second, dest);
          if(!repair_potential(src, dest))
          {
            set_to_bottom();
            return false;
          }
          
          close_over_edge(src, dest);
        }

        // CRAB_WARN("SplitDBM::add_linear_leq not yet implemented.");
        return true;  
      }
   
      void add_disequation(linear_expression_t exp)
      {
        return;
        /*
        // Can only exploit \sum_i c_i x_i \neq k if:
        // (1) exactly one x_i is unfixed
        // (2) lb(x_i) or ub(x_i) = k - \sum_i' c_i' x_i'
        Wt k = exp.constant();
        auto it = exp.begin();
        for(; it != exp.end(); ++it)
        {
          if(!var_is_fixed((*it).second)) 
            break;
          k -= (*it).first*get_value((*it).second);
        }

        // All variables are fixed
        if(it == exp.end())
        {
          if(k == Wt(0))
            set_to_bottom();
          return;
        }

        // Found one unfixed variable; collect the rest.
        Wt ucoeff = (*it).first;
        VariableName uvar((*it).second;
        interval_t u_int = get_interval(ranges, uvar);
        // We need at least one side of u to be finite.
        if(u_int.lb().is_infinite() && u_int.ub().is_infinite())
          return;

        for(++it; it != exp.end(); ++it)
        {
          // Two unfixed variables; nothing we can do.
          if(!var_is_fixed((*it).second))
            return;
          k -= (*it).first*get_value((*it).second);
        }
        */
      }

      void operator+=(linear_constraint_t cst) {
        if(is_bottom())
          return;

        if (cst.is_tautology())
          return;
      
        if (cst.is_contradiction()){
          set_to_bottom();
          return ;
        }

        if (cst.is_inequality())
        {
          if(!add_linear_leq(cst.expression()))
            set_to_bottom();
          return;
        }

        if (cst.is_equality())
        {
          linear_expression_t exp = cst.expression();
          if(!add_linear_leq(exp) || !add_linear_leq(-exp))
            set_to_bottom();
          return;
        }

        if (cst.is_disequation())
        {
          add_disequation(cst.expression());
        }

        CRAB_WARN("Unhandled constraint in SplitDBM");

        CRAB_DEBUG("---", cst, "\n", *this);
        return;
      }
    
      void operator+=(linear_constraint_system_t csts) {  
        if(is_bottom()) return;

        for(auto cst: csts) {
          operator+=(cst);
        }
      }

      interval_t get_interval(ranges_t& r, variable_t x) {
        return get_interval(r, x.name());
      }

      interval_t get_interval(ranges_t& r, VariableName x) {
        boost::optional< interval_t > v = r.lookup(x);
        if(v)
          return *v;
        else
          return interval_t::top();
      }

      interval_t operator[](VariableName x) { 
        if (is_top())    return interval_t::top();
        if (is_bottom()) return interval_t::bottom();

        if (this->is_bottom()) {
            return interval_t::bottom();
        } else {
          return get_interval(ranges, x);
        }
      }

      void set(VariableName x, interval_t intv) {
        if(is_bottom())
          return;

        this->operator-=(x);
        ranges.insert(x, intv);
      }

      // bitwise_operators_api
      void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
        // since reasoning about infinite precision we simply assign and
        // ignore the width.
        assign(x, linear_expression_t(y));
      }

      void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
        // since reasoning about infinite precision we simply assign
        // and ignore the width.
        assign(x, k);
      }

      void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
        // Convert to intervals and perform the operation
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
          default: 
            CRAB_ERROR("DBM: unreachable");
        }
        set(x, xi);
      }
    
      void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
        // Convert to intervals and perform the operation
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
          default: 
            CRAB_ERROR("DBM: unreachable");
        }
        set(x, xi);
      }
    
      // division_operators_api
    
      void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
        if (op == OP_SDIV){
          apply(OP_DIVISION, x, y, z);
        }
        else{
          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi = operator[](z);
          interval_t xi = interval_t::bottom();
      
          switch (op) {
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
              CRAB_ERROR("spDBM: unreachable");
          }
          set(x, xi);
        }
      }

      void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
        if (op == OP_SDIV){
          apply(OP_DIVISION, x, y, k);
        }
        else{
          // Convert to intervals and perform the operation
          interval_t yi = operator[](y);
          interval_t zi(k);
          interval_t xi = interval_t::bottom();
      
          switch (op) {
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
              CRAB_ERROR("DBM: unreachable");
          }
          set(x, xi);
        }
      }

      // Resore potential after an edge addition
      bool repair_potential(vert_id src, vert_id dest)
      {
        return GrOps::repair_potential(g, potential, src, dest);
      }

      // Restore closure after a single edge addition
      void close_over_edge(vert_id ii, vert_id jj)
      {
        assert(ii != 0 && jj != 0);
        SubGraph<graph_t> g_excl(g, 0);

        Wt c = g_excl.edge_val(ii,jj);

        // There may be a cheaper way to do this.
        for(vert_id se : g_excl.preds(ii))
        {
          Wt wt_sij = g_excl.edge_val(se,ii) + c;

          // assert(g_excl.preds(se).size() > 0);
          assert(g_excl.preds(se).begin() != g_excl.preds(se).end());
          if(se != jj)
          {
            if(g_excl.elem(se, jj))
            {
              if(g_excl.edge_val(se,jj) <= wt_sij)
                continue;

              g_excl.edge_val(se,jj) = wt_sij;
            } else {
              g_excl.add_edge(se, wt_sij, jj);
            }
            
            for(vert_id de : g_excl.succs(jj))
            {
              if(se != de)
              {
                Wt wt_sijd = wt_sij + g_excl.edge_val(jj, de);
                if(g_excl.elem(ii, jj))
                {
                  g_excl.edge_val(se, de) = min(g_excl.edge_val(se, de), wt_sijd);
                } else {
                  g_excl.add_edge(se, wt_sijd, de);
                }
              }
            }
          }
        }

        for(vert_id de : g_excl.succs(jj))
        {
          Wt wt_ijd = g_excl.edge_val(jj, de) + c;
          if(de != ii)
          {
            if(g_excl.elem(ii, de))
            {
              g_excl.edge_val(ii, de) = min(g_excl.edge_val(ii, de), wt_ijd);
            } else {
              g_excl.add_edge(ii, wt_ijd, de);
            }
          }
        }
        // Closure is now updated.
      }
    
      // Restore closure after a variable assignment
      // Assumption: x = f(y_1, ..., y_n) cannot induce non-trivial
      // relations between (y_i, y_j)
      bool close_after_assign(vert_id v)
      {
        // Run Dijkstra's forward to collect successors of v,
        // and backward to collect predecessors
        edge_vector delta; 
        if(!GrOps::close_after_assign(g, potential, v, delta))
          return false;
        GrOps::apply_delta(g, delta);
        return true; 
      }

      bool closure(void)
      {
        // Full Johnson-style all-pairs shortest path
        CRAB_ERROR("SparseWtGraph::closure not yet implemented."); 
      }

      //! copy of x into a new fresh variable y
      void expand (VariableName x, VariableName y) {
        if(is_bottom())
          return;

        assign(x, linear_expression_t(y));
//        dbm ret = NULL;      
//        ret = dbm_expand(get_dbm_index(x), get_dbm_index(y), _dbm);
//        dbm_dealloc(_dbm);
//        swap(_dbm, ret);
      }

      // Output function
      void write(ostream& o) {

        normalize ();

#if 0
        o << "edges={";
        for(vert_id v : g.verts())
        {
          for(vert_id d : g.succs(v))
          {
            if(!rev_map[v] || !rev_map[d])
            {
              CRAB_WARN("Edge incident to un-mapped vertex.");
              continue;
            }
            o << "(" << (*(rev_map[v])) << "," << (*(rev_map[d])) << ":" << g.edge_val(v,d) << ")";
          }
        }
        o << "}";
#endif

#if 0
        cout << "var_map={";
        for (auto &p: _var_map) 
          cout << p.first << "(" << p.first.index () << ") " << "->" << p.second << ";";
        cout << "}\n";
        cout << "rev_map={";
        for (auto &p: _rev_map) 
          cout << p.first << "->" << p.second << ";";
        cout << "}\n";
        cout << "matrix:\n";
        dbm_print_to(cout, _dbm);
#endif 

        if(is_bottom()){
          o << "_|_";
          return;
        }
        else if (is_top()){
          o << "{}";
          return;
        }
        else
        {
          // Intervals
          bool first = true;
          o << "{";
          for(auto p : ranges)
          {
            if(first)
              first = false;
            else
              o << ", ";
            o << p.first << " -> " << p.second;
          }
         
          // Extract all the edges
          SubGraph<graph_t> g_excl(g, 0);
          for(vert_id s : g_excl.verts())
          {
            if(!rev_map[s])
              continue;
            variable_t vs = *rev_map[s];
            for(vert_id d : g_excl.succs(s))
            {
              if(!rev_map[d])
                continue;
              variable_t vd = *rev_map[d];
              if(first)
                first = false;
              else
                o << ", ";
              o << vd << "-" << vs << "<=" << g_excl.edge_val(s, d);
            }
          }
          o << "}";

//          linear_constraint_system_t inv = to_linear_constraint_system ();
//          o << inv;
        }
      }

      linear_constraint_system_t to_linear_constraint_system () {

        normalize ();

        linear_constraint_system_t csts;
    
        if(is_bottom ()) {
          csts += linear_constraint_t (linear_expression_t (Number(1)) == 
                                       linear_expression_t (Number(0)));
          return csts;
        }

        // Extract the intervals
        for(auto p : ranges)
        {
          variable_t x = p.first;
          interval_t b = p.second;

          if(b.lb().is_finite())
            csts += linear_constraint_t(linear_expression_t(x) >= linear_expression_t(*(b.lb().number())));
          if(b.ub().is_finite())
            csts += linear_constraint_t(linear_expression_t(x) <= linear_expression_t(*(b.ub().number())));
        }

        // Extract all the edges
        SubGraph<graph_t> g_excl(g, 0);
        for(vert_id s : g_excl.verts())
        {
          if(!rev_map[s])
            continue;
          variable_t vs = *rev_map[s];
          for(vert_id d : g_excl.succs(s))
          {
            if(!rev_map[d])
              continue;
            variable_t vd = *rev_map[d];
            csts += linear_constraint_t(linear_expression_t(vd) - linear_expression_t(vs) <= linear_expression_t(g_excl.edge_val(s, d)));
          }
        }

        return csts;
      }

      const char* getDomainName () const {return "spDBM";}

    }; // class SplitDBM

//    template<class Var, class Num>
//    int SplitDBM<Var, Num>::count = 0;
  } // namespace domains


  namespace domain_traits {

       using namespace domains;

       template <typename Number, typename VariableName>
       void expand (SplitDBM<Number,VariableName>& inv, 
                    VariableName x, VariableName new_x) {
         inv.expand (x, new_x);
       }
    
       template <typename Number, typename VariableName>
       void normalize (SplitDBM<Number,VariableName>& inv) {
         inv.normalize();
       }
    
       template <typename Number, typename VariableName, typename Iterator >
       void forget (SplitDBM<Number,VariableName>& inv, Iterator it, Iterator end){
         inv.forget (it, end);
       }
    } // namespace domain_traits

} // namespace crab


#endif // SPLIT_DBM_HPP
