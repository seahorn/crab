#include "../../common.hpp"
#include "../../program_options.hpp"
#include <crab/domains/object_domain.hpp>


using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::object_domain_impl;

namespace {
  using variable_vector_t = std::vector<z_var>;
  using variable_or_constant = z_var_or_cst_t;
  using variable_or_constant_vector_t = std::vector<z_var_or_cst_t>;
}

z_cfg_t *cfg1(variable_factory_t &vfac) {

  /*
    Object A = {V_i, V_x}
    Object B = {V_y}

    Object *objA = new Object A;
    Object *objB = new Object B;

    int *i = objA->i;
    int *x = objA->x;
    int *y = objB->y;
    *i := 0;
    *x := 1;
    *y := 0;
    while(*i <= 99) {
      *x = *x + *y;
      *y = *y + 1;
      *i = *i + 1;
    }
    assert(*x >= *y);
   */

  // === Define program variables
  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var B(vfac["objB"], crab::REF_TYPE, 32);
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var y(vfac["y"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var deref_x(vfac["*x"], crab::INT_TYPE, 32);
  z_var deref_y(vfac["*y"], crab::INT_TYPE, 32);
  // === Define memory regions
  z_var rgn_i(vfac["V_i"], crab::REG_INT_TYPE, 32);
  z_var rgn_x(vfac["V_x"], crab::REG_INT_TYPE, 32);
  z_var rgn_y(vfac["V_y"], crab::REG_INT_TYPE, 32);

  // assume dsa node
  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn_i));
  obj1.push_back(variable_or_constant(rgn_x));
  obj1.push_back(variable_or_constant(rgn_y));

  // variable_or_constant_vector_t obj2;
  // obj2.push_back(variable_or_constant(rgn_y));

  // === Create allocation sites
  crab::tag_manager as_man;
  // === Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // === Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &bb3 = cfg->insert("bb3");
  z_basic_block_t &ret = cfg->insert("ret");
  // === Adding CFG edges
  entry.add_succ(bb1);
  bb1.add_succ(bb1_t);
  bb1.add_succ(bb1_f);
  bb1_t.add_succ(bb2);
  bb2.add_succ(bb1);
  bb1_f.add_succ(bb3);
  bb3.add_succ(ret);

  // === Adding statements

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size12(z_number(12), crab::variable_type(crab::INT_TYPE, 32));

  // Intialization of memory regions
  entry.intrinsic("regions_from_memory_object", {}, obj1);
  // entry.intrinsic("regions_from_memory_object", {}, obj2);
  entry.region_init(rgn_i);
  entry.region_init(rgn_x);
  entry.region_init(rgn_y);

  //// Create references
  entry.make_ref(A, rgn_i, size12, as_man.mk_tag());
  // entry.make_ref(B, rgn_y, size4, as_man.mk_tag());
  entry.gep_ref(i, rgn_i, A, rgn_i, z_number(0));
  entry.gep_ref(x, rgn_x, A, rgn_i, z_number(4));
  entry.gep_ref(y, rgn_y, x, rgn_x, z_number(4));
  //// *i := 0;
  entry.store_to_ref(i, rgn_i, zero32);
  //// *x := 1;
  entry.store_to_ref(x, rgn_x, one32);
  //// *y := 0;
  entry.store_to_ref(y, rgn_y, zero32);
  //// assume(*i <= 99);
  bb1_t.load_from_ref(deref_i, i, rgn_i);
  bb1_t.assume(deref_i <= 99);
  //// assume(*i >= 100);
  bb1_f.load_from_ref(deref_i, i, rgn_i);
  bb1_f.assume(deref_i >= 100);
  //// *x = *x + *y
  bb2.load_from_ref(deref_x, x, rgn_x);
  bb2.load_from_ref(deref_y, y, rgn_y);
  bb2.add(deref_x, deref_x, deref_y);
  bb2.store_to_ref(x, rgn_x, deref_x);
  //// *y = *y + 1
  bb2.load_from_ref(deref_y, y, rgn_y);
  bb2.add(deref_y, deref_y, 1);
  bb2.store_to_ref(y, rgn_y, deref_y);
  //// *i = *i + 1
  bb2.load_from_ref(deref_i, i, rgn_i);
  bb2.add(deref_i, deref_i, 1);
  bb2.store_to_ref(i, rgn_i, deref_i);
  //// assert(*x >= *y)
  bb3.load_from_ref(deref_x, x, rgn_x);
  bb3.load_from_ref(deref_y, y, rgn_y);
  ret.assertion(deref_x >= deref_y);
  return cfg;
}


int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  z_obj_sdbm_t init;
  variable_factory_t vfac;
  z_cfg_t *cfg = cfg1(vfac);
  crab::outs() << *cfg << "\n";

  run(cfg, cfg->entry(), init, false, 2, 2, 20, stats_enabled);

  return 0;
}