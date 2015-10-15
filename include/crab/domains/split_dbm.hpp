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
      typedef patricia_tree< VariableName, interval_t > intervals_t;
      typedef intervals_t::key_binary_op_t key_binary_op_t;
      
      // Eventually break this out into a template param
      typedef Number Wt;
      typedef SparseWtGraph<Wt> graph_t;
      typedef boost::container::flat_map<variable_t, graph_t::vert_id> vert_map_t;

      typedef SplitDBM<Number, VariableName> DBM_t;

      typedef GraphOps<Wt>::edge_vector edge_vector;

      //================
      // Domain data
      //================
      intervals_t ranges; // Intervals for each variable
      vert_map_t vert_map; // Mapping from variables to vertices
      graph_t g; // The underlying relation graph
      vector<Wt> potential; // Stored potential for the vertex

      bool is_bottom;

     private:
      SplitDBM(bool is_bottom):
        writeable(), is_bottom(_is_bottom)
      {

      }

      // FIXME: Rewrite to avoid copying if o is _|_
      SplitDBM(DBM_t& o)
        : writeable()
          numerical_domain<Number, VariableName >(),
          bitwise_operators< Number, VariableName >(),
          division_operators< Number, VariableName >(),
          ranges(o.ranges), 
          vert_map(o.vert_map),
          g(o.g)
      {
        if(o._is_bottom)
          set_to_bottom();
      }

      // We should probably use the magical rvalue ownership semantics stuff.
      SplitDBM(intervals_t& _ranges, vert_map_t& _vert_map, graph_t& _g, vector<Wt>& _potential)
        : writeable()
          numerical_domain<Number, VariableName >(),
          bitwise_operators< Number, VariableName >(),
          division_operators< Number, VariableName >(),
          ranges(_ranges), vert_map(_vert_map), g(_g), potential(_potential)
      { }

      SplitDBM& operator=(SplitDBM& o)
      {
        if(this != &o)
        {
          ranges = o.ranges;
          vert_map = o.vert_map;
          g = o.g;
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
        g.clear();
        _is_bottom = true;
      }

      void apply(operation_t op, VariableName x, id_t y, id_t z) {
        switch(op){
          case OP_ADDITION:
            {
              exp_t exp = exp_add(exp_var(y), exp_var(z));
              assign(x, exp);
            }
            break;
          case OP_SUBTRACTION:
            {
              exp_t exp = exp_sub(exp_var(y), exp_var(z));
              assign(x, exp);
            }
            break;
          case OP_MULTIPLICATION:
            {
              exp_t exp = exp_mul(exp_var(y), exp_var(z));
              assign(x, exp);
            }
            break;
          case OP_DIVISION:
            {
              exp_t exp = exp_div(exp_var(y), exp_var(z));
              assign(x, exp);
            }
            break;
          default:;;
        }
      }	

      void apply(operation_t op, VariableName x, id_t y, Number k) {	
        int z = DBM_impl::ntoi<Number>(k);
        switch(op){
          case OP_ADDITION:
            {
              exp_t exp = exp_add(exp_var(y), exp_const(z));
              assign(x, exp);
            }
            break;
          case OP_SUBTRACTION:
            {
              exp_t exp = exp_sub(exp_var(y), exp_const(z));
              assign(x, exp);
            }
            break;
          case OP_MULTIPLICATION:
            {
              exp_t exp = exp_mul(exp_var(y), exp_const(z));
              assign(x, exp);
            }
            break;
          case OP_DIVISION:
            {
              exp_t exp = exp_div(exp_var(y), exp_const(z));
              assign(x, exp);          
            }
            break;
          default:;;
        }
      }

      // check satisfiability of cst using intervals
      // Only to be used if cst is too hard for dbm
      bool intervals_check_sat(linear_constraint_t cst)  {
        if (is_top())    return true;
        if (is_bottom()) return false;

        auto vars = cst.variables();
        intervals_t inv;
        for(auto v: vars)
          inv.set (v.name(), operator[](v.name()));
        inv += cst;
        return !inv.is_bottom();
      }

     public:

      static DBM_t top() { return SplitDBM(false); }
    
      static DBM_t bottom() { return SplitDBM(true); }
    
     public:

      bool is_bottom() {
        if(!_is_bottom && g.has_negative_cycle())
          _is_bottom = true;
        return _is_bottom;
      }
    
      bool is_top() {
        return ranges.size() == 0 && g.is_empty();
      }
    
      bool operator<=(DBM_t o)  {
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
          if(!(ranges <= o.ranges))
            return false;

          // Set up a mapping from o to this.
          vector<unsigned int> vert_renaming(o.g.size());
          for(auto p : o.vert_map)
          {
            vert_map_t::iterator it = vert_map.find(p.first);
            // We can't have this <= o if we're missing some
            // vertex.
            if(it == vert_map.end())
              return false;
            vert_renaming[p.second] = (*it).second;
          }

          PermView<graph_t> g_perm(vert_renaming, g);
          // FIXME: Check incorporation of bounds
          for(auto e : o.edges())
          {
            vert_id x = vert_renaming[e.src];
            vert_id y = vert_renaming[e.dest];
            // Is the edge implied by ranges?
            if(ranges[y].ub() - ranges[x].lb() <= bound_t(e.weight))
              continue;
            // Edge not present
            if(!g_perm.mem(e.src, e.dest))
              return false;
            if(!(g_perm.get(e.src, e.dest) <= e.weight))
              return false;
          }
          return true;
        }
      }

      class Wt_max {
      public:
       Wt_max() { } 
       Wt apply(Wt& x, Wt& y) { return max(x, y); }
       bool default_is_absorbing() { return true; }
      };

      class Wt_min {
      public:
        Wt_min() { }
        Wt apply(Wt& x, Wt& y) { return min(x, y); }
        bool default_is_absorbing() { return false; }
      }

      // Perform a join on intervals, but record which
      // lower/upper bounds move in each direction.
      class interval_join_t : public key_binary_op_t {
      public:
        interval_join(
            vector<vert_id>& _lb_up,
            vector<vert_id>& _lb_down,
            vector<vert_id>& _ub_up,
            vector<vert_id>& _ub_down)
          : lb_up(_lb_up), lb_down(_lb_down),
            ub_up(_ub_up), ub_down(_ub_down)
        { }

        boost::optional<interval_t> apply(
            VariableName v, interval_t x, interval_t y)
        {
          interval_t z = x.operator|(y);
          if (z.is_top()) {
            return boost::optional< Value >();
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
            return boost::optional< Value >(z);
          }
        }

        bool default_is_absorbing() { return true; }
        vector<vert_id>& lb_up;
        vector<vert_id>& lb_down;
        vector<vert_id>& ub_up;
        vector<vert_id>& ub_down;
      }; // class binary_op
      typedef typename separate_domain<Number, VariableName>::meet_op interval_meet_t;
      typedef typename separate_domain<Number, VariableName>::widening_op interval_widen_t;

      DBM_t operator|(DBM_t o) {
        if (is_bottom() || o.is_top ())
          return o;
        else if (is_top () || o.is_bottom())
          return *this;
        else {
          CRAB_DEBUG ("Before join:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

          // First, join the intervals, collecting change directions
          vector<vert_id> lb_up;
          vector<vert_id> lb_down;
          vector<vert_id> ub_up;
          vector<vert_id> ub_down;
          interval_join_t join_op(lb_up, lb_down, ub_up, ub_down);
          
          intervals_t join_range = ranges.key_merge(o.ranges, join_op);

          // Figure out the common renaming, initializing the
          // resulting potentials as we go.
          vector<vert_id> perm_x;
          vector<vert_id> perm_y;
          vector<VariableName> perm_inv;

          vector<Wt> join_pot;
          vert_map_t out_vmap;
          for(auto p : vert_map)
          {
            auto it = o.vert_map.find(p.first); 
            // Variable exists in both
            if(it != o.vert_map.end())
            {
              out_vmap.insert(p.first, perm_x.size());
              join_pot.push_back(potential[p.second]);
              perm_inv.push_back(p.first);
              perm_x.push_back(p.second);
              perm_y.push_back((*it).second);
            }
          }
          unsigned int sz = perm_x.size();

          // Build the permuted view of x and y.
          PermGraph gx(perm_x, g);
          PermGraph gy(perm_y, o.g);

          // Compute the deferred relations
          graph_t g_ix_ry(sz);
          for(vert_id s : gy.verts())
          {
            for(vert_id d : gy.succs(s))
            {
              // Assumption: gx.mem(s, d) -> gx.edge_val(s, d) <= ranges[var(s)].ub() - ranges[var(d)].lb()
              // That is, if the relation exists, it's at least as strong as the bounds.
              if(!gx.mem(s, d))
              {
                // Check the bounds implied by o.ranges
                bound_t b = ranges[perm_inv[d]].ub() - ranges[perm_inv[s]].lb();
                if(b.is_finite())
                {
                  g_ix_ry.add_edge(s, *(b.number()), d);
                }
              }
            }
          }
          // Apply the deferred relations, and re-close.
          // Assumption: potential is valid for relations _and bounds_.
          // GKG: Make sure this is preserved, or fix.
          edge_vector delta;
          graph_t g_rx(GraphOps<Wt>::meet(gx, g_ix_ry));
          GraphOps<Wt>::close_after_meet(g_rx, potential, gx, g_ix_ry, delta);
          GraphOps<Wt>::apply_delta(g_rx, delta);          

          graph_t g_rx_iy(sz);
          for(vert_id s : gx.verts())
          {
            for(vert_id d : gx.succs(s))
            {
              // Assumption: gx.mem(s, d) -> gx.edge_val(s, d) <= ranges[var(s)].ub() - ranges[var(d)].lb()
              // That is, if the relation exists, it's at least as strong as the bounds.
              if(!gy.mem(s, d))
              {
                // Check the bounds implied by o.ranges
                bound_t b = o.ranges[perm_inv[d]].ub() - o.ranges[perm_inv[s]].lb();
                if(b.is_finite())
                {
                  g_iy_rx.add_edge(s, *(b.number()), d);
                }
              }
            }
          }
          delta.clear();
          graph_t g_ry(GraphOps<Wt>::meet(gy, g_iy_rx));
          GraphOps<Wt>::close_after_meet(g_ry, o.potential, gy, g_iy_rx, delta);
          GraphOps<Wt>::apply_delta(g_ry, delta);
           
          // We now have the relevant set of relations. Because g_rx and g_ry are closed,
          // the result is also closed.
          Wt_min min_op;
          graph_t join_g(g_rx, g_ry);
          // Now reapply the missing independent relations.
          for(vert_id s : ub_up)
          {
            for(vert_id d : lb_up)
            {
              if(s == d)
                continue;
              bound_t b = max(ranges[perm_inv[d]].ub() - ranges[perm_inv[s]].lb(),
                        o.ranges[perm_inv[d]].ub() - o.ranges[perm_inv[s]].lb());
              g_join.update_edge(s, *(b.number()), d, min_op);
            }
          }
          // Conjecture: g_join remains closed.
          
          DBM_t res(join_range, join_map, join_g, join_pot);
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
          intervals_t widen_range = ranges.merge(o.ranges, widen_op);    
          
          // Figure out the common renaming
          vector<vert_id> perm_x;
          vector<vert_id> perm_y;
          vert_map_t out_vmap;
          vector<Wt> widen_pot;
          for(auto p : vert_map)
          {
            auto it = o.vert_map.find(p.first); 
            // Variable exists in both
            if(it != o.vert_map.end())
            {
              out_vmap.insert(p.first, perm_x.size());
              widen_pot.push(potential[p.second]);
              perm_x.push(p.second);
              perm_y.push((*it).second);
            }
          }

          // Build the permuted view of x and y.
          PermGraph gx(perm_x, g);            
          PermGraph gy(perm_y, o.g);
         
          // Now perform the widening 
          graph_t widen_g(GraphOps::widen(gx, gy));
          DBM_t res(widen_range, widen_map, widen_g, widen_pot);
           
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

          interval_t meet_range = ranges.merge(o.ranges, meet_op);

          // Vertices in the left operand map to themselves
          vert_map_t meet_verts(vert_map);
          vector<variable_t> meet_rev(var_rev);

          // Set up mapping.
          vector<vert_id> r_map;
          for(vert_id vr = 0 : o.g.verts())
          {
            variable_t v = o.var_rev[vr]; 
            auto it = meet_verts.find(o.var_rev[vr]);

            vert_id vv;
            if(it == meet_verts.end())
            {
              vv = meet_rev.size();
              meet_verts.insert(v, vv);
              meet_rev.push_back(v);
            } else {
              vv = (*it).second();
            }
            r_map.push_back(vv);
          }

          // Build the permuted view of x and y.
          PermGraph gx(perm_x, g);            
          PermGraph gy(perm_y, o.g);

          // Compute the syntactic meet of the permuted graphs.
          graph_t meet_g(GraphOps::meet(gx, gy));
          
          // Compute updated potentials
          vector<Wt> meet_pi(meet_g.size());
          if(!GraphOps::select_potentials(meet_g, meet_pi))
          {
            // Potentials cannot be selected -- state is infeasible.
            return bottom();
          }

          edge_vector delta;
          GraphOps::close_after_meet(meet_g, meet_pi, gx, gy, delta);
          GraphOps::apply_delta(meet_g, delta);
          for(vert_id v = 0; v < delta.size(); v++)
          {
            for(auto p : delta[v])
            {
              // Update intervals
                        
            }
          }
          
          CRAB_DEBUG ("Result meet:\n",res);
          return DBM_t(meet_range, meet_verts, meet_g, meet_pi);
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
          return DBM_t(*this);

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

          dbm dbm_x = dbm_copy (this->_dbm);    
          dbm dbm_y = dbm_copy (o._dbm);    
          
          // Resize before renaming          
          int max_sz = max (dbm_x->sz, dbm_y->sz);
          if (id >= max_sz - 2)
            max_sz = id + Delta;

          if (dbm_x->sz != max_sz) {
            dbm tmp = resize (dbm_x, max_sz);
            dbm_dealloc(dbm_x);    
            swap(dbm_x, tmp);
          }

          if (dbm_y->sz != max_sz) {
            dbm tmp = resize (dbm_y, max_sz);
            dbm_dealloc(dbm_y);    
            swap(dbm_y, tmp);
          }

          assert (dbm_x->sz == dbm_y->sz && dbm_y->sz == max_sz);

          dbm dbm_xx = dbm_rename (&subs_x[0], subs_x.size(), dbm_x);
          dbm dbm_yy = dbm_rename (&subs_y[0], subs_y.size(), dbm_y);

          dbm_dealloc(dbm_x);
          dbm_dealloc(dbm_y);

          DBM_t res (dbm_narrowing(dbm_xx, dbm_yy), 
                     id, var_map, rev_map);

          dbm_dealloc(dbm_xx);
          dbm_dealloc(dbm_yy);
          */
           

          CRAB_DEBUG ("Result narrowing:\n",res);
          return res;
        }
      }	

      void normalize() {
        // dbm_canonical(_dbm);
        // Always maintained in normal form, except for widening
      }

      void operator-=(VariableName v) {
        if (is_bottom ())
          return;

        auto it = _var_map.find (v);
        if (it != _var_map.end ()) {
          forget (it->second);
          _var_map.erase (v);
          _rev_map.erase (it->second);
        }
      }

      template<typename Iterator>
      void forget (Iterator vIt, Iterator vEt) {
        if (is_bottom ())
          return;
        vector<int> idxs;
        for (auto v: boost::make_iterator_range (vIt,vEt)) {
          auto it = _var_map.find (v);
          if (it != _var_map.end ()) {
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
      }

      // Evaluate the potential value of a variable.
      Wt pot_value(variable_t v)
      {
        vert_map_t::iterator it = vert_map.find(v); 
        if(it != vert_map.end())
          return potential[(*it).second];

        interval_t r(ranges[v]);
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

      void assign(VariableName x, linear_expression_t e) {

        if(is_bottom())
          return;

        // FIXME: Implement
        this->operator-=(x);
        /*
        if (e.is_constant()){
          exp_t exp = exp_const(DBM_impl::ntoi<Number>(e.constant()));
          assign (x,exp);
        }
        else if  (optional<variable_t> v = e.get_variable()){
          VariableName y = (*v).name();
          if (!(x==y)){
            exp_t exp = exp_var (get_dbm_index(y));        
            assign (x,exp);
          }
        }
        else {
          CRAB_WARN("DBM only supports a cst or var on the rhs of assignment");
          this->operator-=(x);
        }
        */

        // Only need to update potentials and paths for x
          

        CRAB_DEBUG("---", x, ":=", e,"\n",*this);
      }

      void apply(operation_t op, VariableName x, VariableName y, VariableName z){	
        if(is_bottom())
          return;

        // FIXME: implement
        this->operator-=(x);
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

        CRAB_DEBUG("---", x, ":=", y, op, z,"\n", *this);
      }
      */

    
      void apply(operation_t op, VariableName x, VariableName y, Number k) {	

        if(is_bottom())
          return;

        // FIXME: Implement
        this->operator-=(x);
        /*
        if (x == y){
          // to make sure that lhs does not appear on the rhs
          graph_t::vert_id v = g.new_vertex();
          
          apply(op, x, get_tmp(), k);
          forget(get_tmp());
        }
        else{
          apply(op, x, get_dbm_index(y), k);
        }
        */

        CRAB_DEBUG("---", x, ":=", y, op, k,"\n", *this);
      }
   
      // Assumption: state is currently feasible.
      bool add_linear_leq(linear_expression_t& exp)
      {
        return true;  
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
          return add_linear_leq(cst.expression());

        if (cst.is_equality())
        {
          linear_expression_t = cst.expression();
          return add_linear_leq(exp) && add_linear_leq(-exp);
        }

        if (cst.is_disequality())
        {
          // FIXME: add handling  
        }

        CRAB_WARN("Unhandled constraint in SplitDBM");

        CRAB_DEBUG("---", cst, "\n", *this);
        return true;
      }
    
      void operator+=(linear_constraint_system_t csts) {  
        if(is_bottom()) return;

        for(auto cst: csts) {
          operator+=(cst);
        }
      }

      interval_t operator[](VariableName x) { 
        if (is_top())    return interval_t::top();
        if (is_bottom()) return interval_t::bottom();

        return ranges[x];
      }

      void set(VariableName x, interval_t intv) {
        if(is_bottom())
          return;

        this->operator-=(x);
        ranges.set(x, intv);
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
        /*
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
        */
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
              CRAB_ERROR("DBM: unreachable");
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

      // Restore closure after a single edge addition
      bool close_over_edge(vert_id ii, vert_id jj)
      {
        Wt c = edge_val(ii,jj);

        // There may be a cheaper way to do this.
        for(vert_id se : preds(ii))
        {
          Wt wt_sij = edge_val(se,ii) + c;

          assert(preds(se).size() > 0);
          if(se != jj)
          {
            if(elem(se, jj))
            {
              if(edge_val(se,jj) <= wt_sij)
                continue;

              edge_val(se,jj) = wt_sij;
            } else {
              add_edge(se, jj, wt_sij);
            }
            
            for(vert_id de : succs(jj))
            {
              if(se != de)
              {
                int wt_sijd = wt_sij + edge_val(jj, de);
                if(elem(ii, jj))
                {
                  edge_val(se, de) = min(edge_val(se, de), wt_sijd);
                } else {
                  add_edge(se, de, wt_sijd);
                }
              }
            }
          }
        }

        for(vert_id de : succs(jj))
        {
          int wt_ijd = edge_val(jj, de) + c;
          if(de != ii)
          {
            if(elem(ii, de))
            {
              edge_val(ii, de) = min(edge_val(ii, de), wt_ijd);
            } else {
              add_edge(ii, de, wt_ijd);
            }
          }
        }
        // Closure is now updated.
      }
    
      // Restore closure after a variable assignment
      // Assumption: x = f(y_1, ..., y_n) cannot induce non-trivial
      // relations between (y_i, y_j)
      bool close_assign(vert_id v)
      {
        // Run Dijkstra's forward to collect successors of v,
        // and backward to collect predecessors
        CRAB_ERROR("SparseWtGraph::close_assign not yet implemented."); 
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

        dbm ret = NULL;      
        ret = dbm_expand(get_dbm_index(x), get_dbm_index(y), _dbm);
        dbm_dealloc(_dbm);
        swap(_dbm, ret);
      }

      // Output function
      void write(ostream& o) {

        normalize ();
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
          linear_constraint_system_t inv = to_linear_constraint_system ();
          o << inv;
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

        // GKG: Implement me.
        return csts;
      }

      const char* getDomainName () const {return "spDBM";}

    }; // class SplitDBM

  } // namespace domains


  namespace domain_traits {

       using namespace domains;

       template <typename Number, typename VariableName>
       void expand (DBM<Number,VariableName>& inv, 
                    VariableName x, VariableName new_x) {
         inv.expand (x, new_x);
       }
    
       template <typename Number, typename VariableName>
       void normalize (DBM<Number,VariableName>& inv) {
         inv.normalize();
       }
    
       template <typename Number, typename VariableName, typename Iterator >
       void forget (DBM<Number,VariableName>& inv, Iterator it, Iterator end){
         inv.forget (it, end);
       }
    } // namespace domain_traits

} // namespace crab


#endif SPLIT_DBM_HPP
