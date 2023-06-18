#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>
#include <crab/transforms/lower_safe_assertions.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::checker;

z_cfg_t *prog(variable_factory_t &vfac, crab::tag_manager &as_man) {
  /*
    i := 0;
    x := 1;
    y := 0;
    p := malloc(...);
    l := p;
    while(i <= 99) {
       x+=y;
       y++;
       *p := i;
       q  := malloc(...)
       p->next = q;
       p := p->next;
       i++;
    }
    assume(x <= y);
    assert(i == 100);
    assert(i >= 200);
    assert(x >= y);
    assert(x >= 200);
   */

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var l(vfac["l"], crab::REF_TYPE);
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var p_next(vfac["p_next"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);
  // create an  memory region
  z_var mem(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;


  // adding statements
  entry.assign(i, 0);
  entry.assign(x, 1);
  entry.assign(y, 0);
  entry.make_ref(p, mem, size8, as_man.mk_tag());
  entry.gep_ref(l, mem, p, mem, z_number(0));
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(x, x, y);
  bb2.add(y, y, 1);
  // *p := i
  bb2.store_to_ref(p, mem, i);
  //  q := malloc(...)
  bb2.make_ref(q, mem, size8, as_man.mk_tag());
  //  p_next := p+4
  bb2.gep_ref(p_next, mem, p, mem, z_number(4));
  //  *p_next := q
  bb2.store_to_ref(p_next, mem, q);
  //  p := p_next
  bb2.gep_ref(p, mem, p_next, mem, z_number(0));
  bb2.add(i, i, 1);
  ret.assume(x <= y);
  ret.assertion(i == 100);
  ret.assertion(i >= 200);
  ret.assertion(x >= y);
  ret.assertion(x >= 200);

  return cfg;
}

// // Print invariants by traversing the cfg in dfs.
// template<typename Analyzer>
// static void print_invariants(z_cfg_ref_t cfg, Analyzer& analyser) {
//   std::set<crab::cfg_impl::basic_block_label_t> visited;
//   std::vector<crab::cfg_impl::basic_block_label_t> worklist;
//   worklist.push_back(cfg.entry());
//   visited.insert(cfg.entry());
//   while (!worklist.empty()) {
//     auto cur_label = worklist.back();
//     worklist.pop_back();

//     auto pre = analyser.get_pre (cur_label);
//     auto post = analyser.get_post (cur_label);
//     crab::outs() << get_label_str (cur_label) << "="
//               << pre
//               << " ==> "
//               << post << "\n";

//     auto const &cur_node = cfg.get_node (cur_label);
//     for (auto const kid_label :
//     boost::make_iterator_range(cur_node.next_blocks())) {
//       if (visited.insert(kid_label).second) {
// 	worklist.push_back(kid_label);
//       }
//     }
//   }

// }

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  using num_analyzer_t = intra_fwd_analyzer<z_cfg_ref_t, z_sdbm_domain_t>;
  using num_checker_t = intra_checker<num_analyzer_t>;
  using assert_prop_num_checker_t = assert_property_checker<num_analyzer_t>;

  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_cfg_ref_t cfg_ref(*(prog(vfac, as_man)));

  // Run numerical analysis
  z_sdbm_domain_t absval_fac;
  crab::fixpoint_parameters fixpo_params;  
  num_analyzer_t num_a(cfg_ref, absval_fac, nullptr, fixpo_params);
  num_a.run(absval_fac.make_top());
  crab::outs() << "Analysis using " << absval_fac.domain_name() << "\n";
  // print_invariants(cfg_ref, num_a);

  // Run the checker
  const int verbose = 3;
  typename num_checker_t::prop_checker_ptr prop(
      new assert_prop_num_checker_t(verbose));
  num_checker_t checker(num_a, {prop});
  checker.run();
  checker.show(crab::outs());

  // Transformation
  crab::outs() << "Before replacing safe assertions with assume statements\n"
               << cfg_ref << "\n";
  std::set<const z_cfg_ref_t::statement_t *> safe_checks;
  safe_checks.insert(prop->get_safe_checks().begin(),
                     prop->get_safe_checks().end());
  crab::transforms::lower_safe_assertions<z_cfg_ref_t> lsa(safe_checks);
  lsa.run(cfg_ref);
  crab::outs() << "After replacing safe assertions with assume statements\n"
               << cfg_ref << "\n";
  return 0;
}
