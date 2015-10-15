#include "../common.hpp"
#include <crab/common/sparse_graph.hpp>
#include <crab/common/graph_ops.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab;

typedef linear_constraint<z_number, varname_t> linear_constraint_t;
typedef linear_expression<z_number, varname_t> linear_expression_t;

typedef SparseWtGraph<z_number> graph_t;

void check_graph(graph_t& g)
{
  for(int v = 0; v < g.size(); v++)
  {
    for(int d : g.succs(v))
      assert(g.elem(v, d));
    for(int s : g.preds(v))
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

int main (int argc, char** argv )
{
  graph_t x(10);
  x.add_edge(0, -1, 2);
  x.add_edge(1, 5, 2);
  x.add_edge(0, 2, 1);
  x.add_edge(2, 1, 1);

  graph_t y(10);
  y.add_edge(0, 3, 2);
  y.add_edge(0, 3, 1);
  cout << x << y << endl;

  check_graph(x);
  check_graph(y);

  vector<z_number> x_pot(x.size());
  if(!GraphOps<z_number>::select_potentials(x, x_pot))
    assert(0 && "Should be feasible.");
  cout << "x model: "; write_vec(cout, x_pot); cout << endl;

  graph_t g = GraphOps<z_number>::join(x, y);
  cout << g << endl;

  return 0;
}
