#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/dataflow/assertion_crawler.hpp>

#include <boost/range/iterator_range.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

z_cfg_t *prog(variable_factory_t &vfac) {
  /*
     i := 0;
     x := 1;
     y := 0;
     z := 3;
     w := 3;
     while (i < 100) {
       x  := x + y;
       y  := y + 1;
       nd := *;
       z  := z xor nd;
       w  := w xor nd;
       i  := i + 1;
     }
     assert(i >= 100);
   */
  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var nd1(vfac["nd1"], crab::INT_TYPE, 32);
  z_var nd2(vfac["nd2"], crab::INT_TYPE, 32);
  // entry and exit block
  function_decl<z_number, varname_t> decl("main", {}, {});  
  z_cfg_t *cfg = new z_cfg_t("entry", "ret", decl);
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &exit = cfg->insert("exit");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> exit;
  exit >> ret;
  // adding statements
  entry.assign(i, 0);
  entry.assign(x, 1);
  entry.assign(y, 0);
  entry.assign(z, 3);
  entry.assign(w, 3);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.callsite("bar",{x,y}, {x,y,z,w});
  bb2.add(i, i, 1);
  exit.assume(x <= y);
  exit.assertion(w >= 0);
  exit.assertion(x >= y);
  exit.assertion(i >= 0);    
  return cfg;
}

z_cfg_t *foo_cfg(variable_factory_t &vfac) {
  // Definining program variables
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var i3(vfac["i3"], crab::INT_TYPE, 32);
  z_var i4(vfac["i4"], crab::INT_TYPE, 32);
  z_var o1(vfac["o1"], crab::INT_TYPE, 32);
  z_var o2(vfac["o2"], crab::INT_TYPE, 32);  

  z_var tmp1(vfac["tmp1"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);  
  
  // entry and exit block
   z_cfg_t *cfg = new z_cfg_t("entry", "exit",
			      function_decl<z_number, varname_t>("foo", {i1,i2,i3,i4}, {o1,o2}));
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");  

  entry >> exit;
  entry.add(tmp1, i1,i2);
  entry.add(tmp2, i3,i4);
  
  exit.assign(o1, tmp1);
  exit.assign(o2, tmp2);
  
  return cfg;
}

z_cfg_t *bar_cfg(variable_factory_t &vfac) {
#if 0   
  // Definining program variables
  z_var a1(vfac["a1"], crab::INT_TYPE, 32);
  z_var a2(vfac["a2"], crab::INT_TYPE, 32);
  z_var a3(vfac["a3"], crab::INT_TYPE, 32);
  z_var a4(vfac["a4"], crab::INT_TYPE, 32);
  z_var b1(vfac["b1"], crab::INT_TYPE, 32);
  z_var b2(vfac["b2"], crab::INT_TYPE, 32);  

  // entry and exit block
   z_cfg_t *cfg = new z_cfg_t("entry", "exit",
			      function_decl<z_number, varname_t>("bar", {a1,a2,a3,a4}, {b1,b2}));
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  entry.callsite("foo", {b1,b2}, {a1,a2,a3,a4});
  exit.assertion(b1 >= 0);
  exit.assertion(b2 >= 0);      
#else
  // Definining program variables
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var i3(vfac["i3"], crab::INT_TYPE, 32);
  z_var i4(vfac["i4"], crab::INT_TYPE, 32);
  z_var o1(vfac["o1"], crab::INT_TYPE, 32);
  z_var o2(vfac["o2"], crab::INT_TYPE, 32);  

  // entry and exit block
   z_cfg_t *cfg = new z_cfg_t("entry", "exit",
  			      function_decl<z_number, varname_t>("bar", {i1,i2,i3,i4}, {o1,o2}));
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  entry >> exit;
  entry.callsite("foo", {o1,o2}, {i1,i2,i3,i4});
  exit.assertion(o1 >= 0);
  exit.assertion(o2 >= 0);      
#endif
  
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  using callgraph_t = call_graph<z_cfg_ref_t>;
  variable_factory_t vfac;

  z_cfg_t *p1 = prog(vfac);
  crab::outs() << *p1 << "\n";

  z_cfg_t *p2 = foo_cfg(vfac);
  crab::outs() << *p2 << "\n";

  z_cfg_t *p3 = bar_cfg(vfac);
  crab::outs() << *p3 << "\n";
  
  vector<z_cfg_ref_t> cfgs({*p1, *p2, *p3});
  callgraph_t cg(cfgs);
  
  using crawler_t = crab::analyzer::inter_assertion_crawler<callgraph_t>;
  crawler_t crawler(cg);
  crawler.run();

  auto print_results = [&crawler](z_cfg_t &cfg) {
    // Print results in DFS to enforce a fixed order
    std::set<crab::cfg_impl::basic_block_label_t> visited;
    std::vector<crab::cfg_impl::basic_block_label_t> worklist;
    worklist.push_back(cfg.entry());
    visited.insert(cfg.entry());
    while (!worklist.empty()) {
      auto cur_label = worklist.back();
      worklist.pop_back();
      auto results = crawler.get_results(cfg, cur_label);
      crab::outs() << crab::basic_block_traits<crab::cfg_impl::z_basic_block_t>::to_string(cur_label)
		   << "=" << results << "\n";
      auto const &cur_node = cfg.get_node(cur_label);
      for (auto const& kid_label :
         boost::make_iterator_range(cur_node.next_blocks())) {
	if (visited.insert(kid_label).second) {
	  worklist.push_back(kid_label);
	}
      }
    }};

  crab::outs() << "Assertion Crawler Analysis for main\n";
  print_results(*p1);
  crab::outs() << "Assertion Crawler Analysis for bar\n";  
  print_results(*p3);
  crab::outs() << "Assertion Crawler Analysis for foo\n";    
  print_results(*p2);
  
  //crawler.write(crab::outs());

  delete p1;
  delete p2;
  delete p3;
  return 0;
}
