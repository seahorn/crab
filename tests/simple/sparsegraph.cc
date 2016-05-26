#include "../common.hpp"
#include <crab/domains/graphs/sparse_graph.hpp>
#include <crab/domains/graphs/pt_graph.hpp>
#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/domains/graphs/adapt_sgraph.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab;

typedef linear_constraint<z_number, varname_t> linear_constraint_t;
typedef linear_expression<z_number, varname_t> linear_expression_t;

//typedef SparseWtGraph<z_number> graph_t;
typedef PtGraph<int64_t> graph_t;
typedef graph_t::vert_id vert_id;

typedef AdaptGraph<int64_t> agraph_t;

typedef GraphOps<graph_t> GrOps;

template<class G>
void check_graph(G& g)
{
  for(vert_id v : g.verts())
  {
    for(vert_id d : g.succs(v))
      assert(g.elem(v, d));
    for(vert_id s : g.preds(v))
      assert(g.elem(s, v));
  }
}

template<class V>
void write_vec(ostream& o, V& v)
{
  o << "[";
  bool first = true;
  for(auto e : v)
  {
    if(first)
      first = false;
    else
      o << ", ";
    o << e;
  }
  o << "]";
}

void test_sdbm(void)
{
  VariableFactory vfac;

  varname_t x = vfac["x"];
  varname_t A = vfac["A"];
  varname_t x_prime = vfac["x\'"];

  sdbm_domain_t dbm = sdbm_domain_t::top ();
  // for all i. A[i] >= 0
  dbm += linear_constraint_t ( linear_expression_t (A) >= z_number (0));
  // x = A[..];
  dbm.expand (A, x_prime);
  dbm.assign (x, z_var (x_prime));
  dbm -= x_prime;
  // if (x <= 0)
  dbm += linear_constraint_t ( linear_expression_t (x) <= z_number (0));
  cout << dbm << "\n";
}

void test_sgraph(void)
{
  agraph_t x;
  for(vert_id v = 0; v < 10; v++)
    x.new_vertex();
  //x.growTo(10);

  x.add_edge(0, -1, 2);
  x.add_edge(1, 5, 2);
  x.add_edge(0, 2, 1);
  x.add_edge(2, 1, 1);

  graph_t y;
  for(vert_id v = 0; v < 10; v++)
    y.new_vertex();
//  y.growTo(10);

  y.add_edge(0, 3, 2);
  y.add_edge(0, 3, 1);
  cout << x << y << "\n";

  check_graph(x);
  check_graph(y);

  vector<int64_t> x_pot(x.size());
  if(!GrOps::select_potentials(x, x_pot))
    assert(0 && "Should be feasible.");
  cout << "x model: "; write_vec(cout, x_pot); cout << "\n";

  graph_t g = GrOps::join(x, y);
  cout << g << "\n";

  x.clear_edges();
  x.add_edge(0, -1, 1);

  y.clear_edges();
  y.add_edge(1, -1, 2);

  // Compute the syntactic meet
  bool is_closed = false;
//  graph_t gm = GrOps::meet(x, y, is_closed);
  graph_t gm = GrOps::meet(x, y);
  vector<int64_t> gm_pot(gm.size());

  // Close
  GrOps::select_potentials(gm, gm_pot);
  if(!is_closed)
  {
    GrOps::edge_vector delta;
    GrOps::close_after_meet(gm, gm_pot, x, y, delta);
    GrOps::apply_delta(gm, delta);
  }

  cout << gm << "\n";

}
int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  test_sgraph();
  test_sdbm();
  return 0;
}
