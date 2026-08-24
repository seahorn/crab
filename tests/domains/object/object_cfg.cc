// End-to-end CFG tests for the object domain: each program is analyzed with
// intra_fwd_analyzer and the assertion checker's verdicts are asserted, both
// as totals and per tagged assertion (via debug_info ids).
#define BOOST_TEST_MODULE object_cfg
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "object_dom.hpp"
#include <crab/analysis/fwd_analyzer.hpp>
#include <crab/checkers/assertion.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/checkers/checker.hpp>
#include <crab/config.h>
#include <crab/domains/object_domain.hpp>
#include <crab/fixpoint/fixpoint_params.hpp>

#include <memory>
#include <string>
#include <vector>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::object_domain_impl;

namespace {
using variable_or_constant = z_var_or_cst_t;
using variable_or_constant_vector_t = std::vector<z_var_or_cst_t>;

// Scope the global object-domain parameters: save on entry, restore on exit,
// so test cases with different configurations cannot leak into each other.
struct scoped_object_params {
  crab::domains::object_domain_params m_saved;
  scoped_object_params(
      crab::domains::object_domain_params::reduction_level_t level)
      : m_saved(crab_domain_params_man::get().reduction_level()) {
    crab::domains::object_domain_params p(level);
    crab_domain_params_man::get().update_params(p);
  }
  ~scoped_object_params() { crab_domain_params_man::get().update_params(m_saved); }
};

// The production configuration the programs below are documented against.
struct opt_reduction_scope : scoped_object_params {
  opt_reduction_scope()
      : scoped_object_params(crab::domains::object_domain_params::
                                 reduction_level_t::REDUCTION_BEFORE_CHECK) {}
};

template <typename Dom>
crab::checker::checks_db analyze_and_check(z_cfg_t &cfg,
                                           unsigned widening_delay,
                                           unsigned descending_iterations,
                                           unsigned max_thresholds) {
  using analyzer_t = crab::analyzer::intra_fwd_analyzer<z_cfg_ref_t, Dom>;
  z_cfg_ref_t cfg_ref(cfg);
  Dom init;
  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening_delay;
  fixpo_params.get_descending_iterations() = descending_iterations;
  fixpo_params.get_max_thresholds() = max_thresholds;
  analyzer_t analyzer(cfg_ref, init.make_top(), nullptr, fixpo_params);
  typename analyzer_t::assumption_map_t assumptions;
  analyzer.run(cfg.entry(), init, assumptions);
  using checker_t = crab::checker::intra_checker<analyzer_t>;
  using prop_t = crab::checker::assert_property_checker<analyzer_t>;
  typename checker_t::prop_checker_ptr prop(new prop_t(0));
  checker_t checker(analyzer, {prop});
  checker.run();
  return checker.get_all_checks();
}

// The assertion tagged with `id` was checked and every verdict is safe (or
// unreachable).
bool assertion_is_safe(const crab::checker::checks_db &db, int64_t id) {
  crab::cfg::debug_info dbg(id);
  if (!db.has_checks(dbg)) {
    return false;
  }
  for (crab::checker::check_kind k : db.get_checks(dbg)) {
    if (k != crab::checker::check_kind::CRAB_SAFE &&
        k != crab::checker::check_kind::CRAB_UNREACH) {
      return false;
    }
  }
  return true;
}

// The assertion tagged with `id` was checked and produced a warning.
bool assertion_is_warning(const crab::checker::checks_db &db, int64_t id) {
  crab::cfg::debug_info dbg(id);
  if (!db.has_checks(dbg)) {
    return false;
  }
  for (crab::checker::check_kind k : db.get_checks(dbg)) {
    if (k == crab::checker::check_kind::CRAB_WARN) {
      return true;
    }
  }
  return false;
}

/*
  Object A = {V_len, V_cap}

  void main() {
    Object A *ary[2];
    for (int i = 0; i < 1; i++) {   // runs for i == 0 only; unrolled below
      ary[i] = malloc(sizeof(Object A));
      int chunk = 10 * i;
      ary[i]->len = chunk;
      ary[i]->cap = chunk + 5;
    }
    assert(ary[1]->len + 5 == ary[1]->cap);   // id 1
  }

  Both allocations come from the same DSA node {V_len, V_cap}, so they are
  instances of the same abstract object. The first make_ref keeps the object
  a singleton and the writes are recorded precisely in its cache; the second
  make_ref promotes the (dirty) cache into the summary, so the loads through
  ary1 still see V_len == 0 and V_cap == 5 and the assertion is proved.
*/
z_cfg_t *prog1(variable_factory_t &vfac, crab::tag_manager &as_man) {
  z_var A0(vfac["ary0"], crab::REF_TYPE, 32);
  z_var A1(vfac["ary1"], crab::REF_TYPE, 32);
  z_var cap_ref0(vfac["cap_ref0"], crab::REF_TYPE, 32);
  z_var cap_ref1(vfac["cap_ref1"], crab::REF_TYPE, 32);
  z_var deref_len1(vfac["*len1"], crab::INT_TYPE, 32);
  z_var deref_cap1(vfac["*cap1"], crab::INT_TYPE, 32);
  z_var rgn_len(vfac["V_len"], crab::REG_INT_TYPE, 32);
  z_var rgn_cap(vfac["V_cap"], crab::REG_INT_TYPE, 32);

  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn_len));
  obj1.push_back(variable_or_constant(rgn_cap));

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t five32(z_number(5), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &ret = cfg->insert("ret");
  entry.add_succ(ret);

  entry.region_init(rgn_len);
  entry.region_init(rgn_cap);
  entry.intrinsic("regions_from_memory_object", {}, obj1);

  entry.make_ref(A0, rgn_len, size8, as_man.mk_tag());
  entry.gep_ref(cap_ref0, rgn_cap, A0, rgn_len, z_number(4));
  entry.store_to_ref(A0, rgn_len, zero32);
  entry.store_to_ref(cap_ref0, rgn_cap, five32);

  entry.make_ref(A1, rgn_len, size8, as_man.mk_tag());
  entry.gep_ref(cap_ref1, rgn_cap, A1, rgn_len, z_number(4));

  ret.load_from_ref(deref_len1, A1, rgn_len);
  ret.load_from_ref(deref_cap1, cap_ref1, rgn_cap);
  ret.assertion(deref_len1 + 5 == deref_cap1, debug_info(1));

  return cfg;
}

#ifdef HAVE_ELINA
/*
  Object A = {V_x, V_y}

  Object *objA = new Object A;
  int *x = &objA->x;  int *y = &objA->y;
  i := 0;  *x := 1;  *y := 0;
  while (i <= 5) {
    xt = *x + 2;  yt = *y + 1;
    *x = xt;      *y = yt;
    i = i + 1;
  }
  assert(*x >= *y);   // id 2

  The loop preserves the intra-object relation V_x - 2*V_y = 1, which with
  *y >= 0 entails *x >= *y. The relation is between two fields of the SAME
  object, exactly what the per-object summary/cache carries across widening
  and hands back to the base domain via the reduction. The relation has a
  coefficient of 2, hence the polyhedra (ELINA pk) base domain.
*/
z_cfg_t *prog2(variable_factory_t &vfac, crab::tag_manager &as_man) {
  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var y(vfac["y"], crab::REF_TYPE);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var deref_x(vfac["*x"], crab::INT_TYPE, 32);
  z_var xt(vfac["xt"], crab::INT_TYPE, 32);
  z_var deref_y(vfac["*y"], crab::INT_TYPE, 32);
  z_var yt(vfac["yt"], crab::INT_TYPE, 32);
  z_var rgn_x(vfac["V_x"], crab::REG_INT_TYPE, 32);
  z_var rgn_y(vfac["V_y"], crab::REG_INT_TYPE, 32);

  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn_x));
  obj1.push_back(variable_or_constant(rgn_y));

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &ret = cfg->insert("ret");
  entry.add_succ(bb1);
  bb1.add_succ(bb1_t);
  bb1.add_succ(bb1_f);
  bb1_t.add_succ(bb2);
  bb2.add_succ(bb1);
  bb1_f.add_succ(bb3);
  bb3.add_succ(ret);

  entry.region_init(rgn_x);
  entry.region_init(rgn_y);
  entry.intrinsic("regions_from_memory_object", {}, obj1);

  entry.make_ref(A, rgn_x, size8, as_man.mk_tag());
  entry.gep_ref(x, rgn_x, A, rgn_x, z_number(0));
  entry.gep_ref(y, rgn_y, A, rgn_x, z_number(4));
  entry.assign(i, 0);
  entry.store_to_ref(x, rgn_x, one32);
  entry.store_to_ref(y, rgn_y, zero32);
  bb1_t.assume(i <= 5);
  bb1_f.assume(i >= 6);
  bb2.load_from_ref(deref_x, x, rgn_x);
  bb2.add(xt, deref_x, 2);
  bb2.load_from_ref(deref_y, y, rgn_y);
  bb2.add(yt, deref_y, 1);
  bb2.store_to_ref(x, rgn_x, xt);
  bb2.store_to_ref(y, rgn_y, yt);
  bb2.add(i, i, 1);
  bb3.load_from_ref(deref_x, x, rgn_x);
  bb3.load_from_ref(deref_y, y, rgn_y);
  ret.assertion(deref_x >= deref_y, debug_info(2));
  return cfg;
}

/*
  Object A = {V_i}, Object B = {V_j, V_x}

  *i := 0;
  while (*i <= 3) {
    *j := 0;
    while (*j <= 3) {
      assert(*i <= *j + 3);   // id 3 -- EXPECTED WARNING, see below
      *i++;  *j++;  *x++;
    }
    *i = *i - *j + 1;
  }

  The asserted invariant relates fields of TWO DIFFERENT objects: V_i is in
  Object A, V_j in Object B. The per-object summary/cache can only retain
  relations among the fields of a single object across widening; a
  cross-object fact lives transiently in the base domain and does not
  survive the widening at the outer loop head. The warning is a structural
  precision boundary of the per-object design, not a bug -- pinned as an
  expected warning.
*/
z_cfg_t *prog3(variable_factory_t &vfac, crab::tag_manager &as_man) {
  z_cfg_t *cfg = new z_cfg_t("entry");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &l1 = cfg->insert("l1");
  z_basic_block_t &l1_entry = cfg->insert("l1_entry");
  z_basic_block_t &l1_cont = cfg->insert("l1_cont");
  z_basic_block_t &l2 = cfg->insert("l2");
  z_basic_block_t &l2_body = cfg->insert("l2_body");
  z_basic_block_t &l2_exit = cfg->insert("l2_exit");

  entry >> l1;
  l1 >> l1_entry;
  l1_entry >> l2;
  l2 >> l2_body;
  l2_body >> l2;
  l2 >> l2_exit;
  l2_exit >> l1_cont;
  l1_cont >> l1;

  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var B(vfac["objB"], crab::REF_TYPE, 32);
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var j(vfac["j"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var deref_j(vfac["*j"], crab::INT_TYPE, 32);
  z_var rgn_i(vfac["V_i"], crab::REG_INT_TYPE, 32);
  z_var rgn_x(vfac["V_x"], crab::REG_INT_TYPE, 32);
  z_var rgn_j(vfac["V_j"], crab::REG_INT_TYPE, 32);

  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn_i));
  variable_or_constant_vector_t obj2;
  obj2.push_back(variable_or_constant(rgn_j));
  obj2.push_back(variable_or_constant(rgn_x));

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  entry.region_init(rgn_i);
  entry.region_init(rgn_j);
  entry.region_init(rgn_x);
  entry.intrinsic("regions_from_memory_object", {}, obj1);
  entry.intrinsic("regions_from_memory_object", {}, obj2);

  entry.make_ref(A, rgn_i, size4, as_man.mk_tag());
  entry.make_ref(B, rgn_j, size8, as_man.mk_tag());
  entry.gep_ref(i, rgn_i, A, rgn_i, z_number(0));
  entry.gep_ref(j, rgn_j, B, rgn_j, z_number(0));
  entry.gep_ref(x, rgn_x, A, rgn_i, z_number(4));
  entry.store_to_ref(i, rgn_i, zero32);
  l1.load_from_ref(deref_i, i, rgn_i);
  l1.assume(deref_i <= 3);
  l1_entry.store_to_ref(j, rgn_j, zero32);

  l2_body.load_from_ref(deref_j, j, rgn_j);
  l2_body.load_from_ref(deref_i, i, rgn_i);
  l2_body.assume(deref_j <= 3);
  l2_body.assertion(deref_i <= deref_j + 3, debug_info(3));
  l2_body.add(deref_i, deref_i, 1);
  l2_body.add(deref_j, deref_j, 1);
  l2_body.store_to_ref(i, rgn_i, deref_i);
  l2_body.store_to_ref(j, rgn_j, deref_j);

  l2_exit.load_from_ref(deref_j, j, rgn_j);
  l2_exit.assume(deref_j >= 4);
  l1_cont.load_from_ref(deref_i, i, rgn_i);
  l1_cont.load_from_ref(deref_j, j, rgn_j);
  l1_cont.assign(deref_i, deref_i - deref_j + 1);
  l1_cont.store_to_ref(i, rgn_i, deref_i);

  return cfg;
}
#endif // HAVE_ELINA
} // namespace

BOOST_AUTO_TEST_CASE(prog1_same_object_relation_is_summarized) {
  opt_reduction_scope params;
  variable_factory_t vfac;
  crab::tag_manager as_man;
  std::unique_ptr<z_cfg_t> cfg(prog1(vfac, as_man));
  cfg->simplify();
  crab::checker::checks_db db = analyze_and_check<z_obj_zones_t>(
      *cfg, 2 /*widening delay*/, 2 /*descending*/, 20 /*thresholds*/);
  BOOST_TEST(db.get_total_safe() == 1u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_is_safe(db, 1),
             "assert(ary1->len + 5 == ary1->cap) must be proved");
}

#ifdef HAVE_ELINA
BOOST_AUTO_TEST_CASE(prog2_intra_object_relation_survives_widening) {
  opt_reduction_scope params;
  variable_factory_t vfac;
  crab::tag_manager as_man;
  std::unique_ptr<z_cfg_t> cfg(prog2(vfac, as_man));
  crab::checker::checks_db db = analyze_and_check<z_obj_pk_t>(
      *cfg, 2 /*widening delay*/, 2 /*descending*/, 20 /*thresholds*/);
  BOOST_TEST(db.get_total_safe() == 1u);
  BOOST_TEST(db.get_total_warning() == 0u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_is_safe(db, 2), "assert(*x >= *y) must be proved");
}

BOOST_AUTO_TEST_CASE(prog3_cross_object_relation_is_a_known_warning) {
  opt_reduction_scope params;
  variable_factory_t vfac;
  crab::tag_manager as_man;
  std::unique_ptr<z_cfg_t> cfg(prog3(vfac, as_man));
  cfg->simplify();
  crab::checker::checks_db db = analyze_and_check<z_obj_pk_t>(
      *cfg, 1 /*widening delay*/, 2 /*descending*/, 20 /*thresholds*/);
  // EXPECTED warning: the assertion relates fields of two different
  // objects, which the per-object summary cannot retain across widening
  // (see the comment on prog3). If this ever becomes safe, the domain
  // gained cross-object precision and this pin should be revisited.
  BOOST_TEST(db.get_total_safe() == 0u);
  BOOST_TEST(db.get_total_warning() == 1u);
  BOOST_TEST(db.get_total_error() == 0u);
  BOOST_TEST(assertion_is_warning(db, 3),
             "assert(*i <= *j + 3) is expected to be a warning");
}
#endif // HAVE_ELINA

int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int idx = 0; idx < argc; ++idx) {
    if (std::string(argv[idx]) == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    args.push_back(argv[idx]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
