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

z_cfg_t *prog(variable_factory_t &vfac) {
  ////
  // Building the CFG
  ////

  /*
    Object A = {V_a, V_b, V_c}
    Object B = {V_d, V_e}

    Object *objA = new Object A;
    Object *objB = new Object B;

    int *p = objA->b;
    int x = *p;
    int *r = objB->d;
    int *s = objB->e;
    x = *r;
    y = *s;
  */

  // Definining program variables
  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var B(vfac["objB"], crab::REF_TYPE, 32);
  z_var a(vfac["a"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE, 32);
  z_var r(vfac["r"], crab::REF_TYPE, 32);
  z_var s(vfac["s"], crab::REF_TYPE, 32);
  z_var rgn1(vfac["V_a"], crab::REG_INT_TYPE, 32);
  z_var rgn2(vfac["V_b"], crab::REG_INT_TYPE, 32);
  z_var rgn3(vfac["V_c"], crab::REG_INT_TYPE, 32);

  z_var rgn4(vfac["V_d"], crab::REG_INT_TYPE, 32);
  z_var rgn5(vfac["V_e"], crab::REG_INT_TYPE, 32);

  // assume dsa node
  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn1));
  obj1.push_back(variable_or_constant(rgn2));
  obj1.push_back(variable_or_constant(rgn3));

  variable_or_constant_vector_t obj2;
  obj2.push_back(variable_or_constant(rgn4));
  obj2.push_back(variable_or_constant(rgn5));

  // Defining size
  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size12(z_number(12), crab::variable_type(crab::INT_TYPE, 32));
  // Create allocation sites
  crab::tag_manager as_man;
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");

  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");

  // adding control flow
  entry.add_succ(bb1);
  bb1.add_succ(bb2);
  bb2.add_succ(ret);
  // adding statements

  // Intialization of memory regions
  entry.intrinsic("regions_from_memory_object", {}, obj1);
  entry.intrinsic("regions_from_memory_object", {}, obj2);
  entry.region_init(rgn1);
  entry.region_init(rgn2);
  entry.region_init(rgn3);
  entry.region_init(rgn4);
  entry.region_init(rgn5);
  // Create object
  entry.make_ref(A, rgn1, size12, as_man.mk_tag());
  entry.make_ref(B, rgn3, size8, as_man.mk_tag());
  bb1.gep_ref(r, rgn2, A, rgn1, z_number(4)); // V_b, p = gep_ref(V_a, &A, 4)

  bb1.load_from_ref(x, r, rgn2); // x = load_ref(V_b, p)

  bb2.gep_ref(r, rgn4, B, rgn4, z_number(0)); // V_d,r = gep_ref(V_d, &B, 0)
  bb2.gep_ref(s, rgn5, r, rgn4, z_number(4)); // V_e,s = gep_ref(V_d, r, 4)
  bb2.load_from_ref(x, r, rgn4); // x = load_ref(V_d, r)
  bb2.load_from_ref(y, s, rgn5); // y = load_ref(V_e, s)
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  z_obj_zones_t init;
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  cfg->simplify();
  crab::outs() << *cfg << "\n";

  run(cfg, cfg->entry(), init, false, 2, 2, 20, stats_enabled);

  return 0;
}