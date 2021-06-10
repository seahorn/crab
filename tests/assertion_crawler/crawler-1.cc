#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/dataflow/assertion_crawler.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

using assertion_crawler_t = crab::analyzer::assertion_crawler<z_cfg_ref_t>;
using assertion_crawler_domain_t = typename assertion_crawler_t::assertion_crawler_domain_t;
using var_dom_t = typename assertion_crawler_domain_t::var_dom_t;
using summary_map_t = typename assertion_crawler_t::summary_map_t;  
using assert_map_t = typename assertion_crawler_t::assert_map_t;

std::unique_ptr<z_cfg_t> main_cfg(variable_factory_t &vfac) {
  /*
     (x1,x2) := foo(y1,y2,y3,y4);
     assert(x1 >= 100);
     assert(x2 >= 100);
   */
  // Definining program variables
  z_var y1(vfac["y1"], crab::INT_TYPE, 32);
  z_var y2(vfac["y2"], crab::INT_TYPE, 32);
  z_var y3(vfac["y3"], crab::INT_TYPE, 32);
  z_var y4(vfac["y4"], crab::INT_TYPE, 32);
  z_var x1(vfac["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac["x2"], crab::INT_TYPE, 32);  
  
  // entry and exit block
  function_decl<z_number, varname_t> decl("main", {}, {});  
  std::unique_ptr<z_cfg_t> cfg(new z_cfg_t("exit", "exit", decl));
  // adding blocks
  //z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  //entry >> exit;
  // adding statements
  //entry.havoc(y1);
  //entry.havoc(y2);
  //entry.havoc(y3);
  //entry.havoc(y4);  
  exit.callsite("foo",{x1,x2},{y1,y2,y3,y4});
  exit.assertion(x1 >= 100);
  exit.assertion(x2 >= 200);    
  return cfg;
}

z_cfg_t *foo_cfg(variable_factory_t &vfac, summary_map_t &summaries) {
  // Definining program variables
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var i3(vfac["i3"], crab::INT_TYPE, 32);
  z_var i4(vfac["i4"], crab::INT_TYPE, 32);
  z_var o1(vfac["o1"], crab::INT_TYPE, 32);
  z_var o2(vfac["o2"], crab::INT_TYPE, 32);  
  
  // entry and exit block
   z_cfg_t *cfg = new z_cfg_t("entry", "entry",
			      function_decl<z_number, varname_t>("foo", {i1,i2,i3,i4}, {o1,o2}));
  // adding blocks
  /* z_basic_block_t &entry =*/cfg->insert("entry");

  // Summary for foo
  assertion_crawler_domain_t summary_foo;
  var_dom_t s1 = var_dom_t::bottom();
  s1 += i1;
  s1 += i2;
  summary_foo.get_second().set(o1, s1);
  var_dom_t s2 = var_dom_t::bottom();
  s2 += i3;
  s2 += i4;
  summary_foo.get_second().set(o2, s2);

  typename summary_map_t::mapped_type summary_foo_info({i1,i2,i3,i4},{o1,o2},summary_foo);
  summaries.insert({&cfg->get_func_decl(), std::move(summary_foo_info)});
  
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  summary_map_t summaries;  
  assert_map_t assert_map;
  variable_factory_t vfac;

  std::unique_ptr<z_cfg_t> p1 = main_cfg(vfac);
  crab::outs() << *p1 << "\n";
  z_cfg_t *p2 = foo_cfg(vfac, summaries);
  crab::outs() << *p2;
  crab::outs() << "Summary table:\n";
  for (auto &kv: summaries) {
    crab::outs() << "\t";
    kv.first.write(crab::outs());
    crab::outs() << " -> ";
    kv.second.write(crab::outs());
    crab::outs() << "\n";
  }
  
  assertion_crawler_t assert_crawler(*p1, assert_map, summaries);
  assert_crawler.exec();
  crab::outs() << "\n";
  assert_crawler.write(crab::outs());
  crab::outs() << "\n";

  delete p2;
  return 0;
}
