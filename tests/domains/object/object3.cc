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

  /*
    Object A = {V_i}
    Object B = {V_j, V_x}

    Object *objA = new Object A;
    Object *objB = new Object B;

    int *i = objA->i;
    int *j = objB->j;
    int *x = objB->x;
    *i := 0;

    while (i<=3) {
      *j := 0;
      while (*j <= 3) {
        assert (*i <= *j + 3);
        *i++;
        *j++;
        *x++;
      }
      *i = *i - *j + 1;
    }
   */
  // === Create allocation sites
  crab::tag_manager as_man;
  // === Create empty CFG
  z_cfg_t *cfg = new z_cfg_t("entry");
  // === Adding CFG blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &l1 = cfg->insert("l1");
  z_basic_block_t &l1_entry = cfg->insert("l1_entry");
  z_basic_block_t &l1_cont = cfg->insert("l1_cont");
  z_basic_block_t &l2 = cfg->insert("l2");
  z_basic_block_t &l2_body = cfg->insert("l2_body");
  z_basic_block_t &l2_exit = cfg->insert("l2_exit");

  entry >> l1;
  // outer loop
  l1 >> l1_entry;
  l1_entry >> l2;
  // inner loop
  l2 >> l2_body;
  l2_body >> l2;
  l2 >> l2_exit;
  // outer loop again
  l2_exit >> l1_cont;
  l1_cont >> l1;

  // === Define program variables
  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var B(vfac["objB"], crab::REF_TYPE, 32);
  z_var i(vfac["i"], crab::REF_TYPE);
  z_var x(vfac["x"], crab::REF_TYPE);
  z_var j(vfac["j"], crab::REF_TYPE);
  z_var deref_i(vfac["*i"], crab::INT_TYPE, 32);
  z_var deref_x(vfac["*x"], crab::INT_TYPE, 32);
  z_var deref_j(vfac["*j"], crab::INT_TYPE, 32);
  // === Define memory regions
  z_var rgn_i(vfac["V_i"], crab::REG_INT_TYPE, 32);
  z_var rgn_x(vfac["V_x"], crab::REG_INT_TYPE, 32);
  z_var rgn_j(vfac["V_j"], crab::REG_INT_TYPE, 32);

  // assume dsa node
  variable_or_constant_vector_t obj1;
  obj1.push_back(variable_or_constant(rgn_i));

  variable_or_constant_vector_t obj2;
  obj2.push_back(variable_or_constant(rgn_j));
  obj2.push_back(variable_or_constant(rgn_x));

  z_var_or_cst_t zero32(z_number(0), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t one32(z_number(1), crab::variable_type(crab::INT_TYPE, 32));  
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));

  // === Adding statements

  // Intialization of memory regions
  entry.intrinsic("regions_from_memory_object", {}, obj1);
  entry.intrinsic("regions_from_memory_object", {}, obj2);
  entry.region_init(rgn_i);
  entry.region_init(rgn_j);
  entry.region_init(rgn_x);

  //// Create references
  entry.make_ref(A, rgn_i, size4, as_man.mk_tag());
  entry.make_ref(B, rgn_j, size8, as_man.mk_tag());
  entry.gep_ref(i, rgn_i, A, rgn_i, z_number(0));
  entry.gep_ref(j, rgn_j, B, rgn_j, z_number(0));
  entry.gep_ref(x, rgn_x, A, rgn_i, z_number(4));
  //// *i := 0;
  entry.store_to_ref(i, rgn_i, zero32);
  l1.load_from_ref(deref_i, i, rgn_i);
  l1.assume(deref_i <= 3);
  //// *j := 0;
  l1_entry.store_to_ref(j, rgn_j, zero32);

  l2_body.load_from_ref(deref_j, j, rgn_j);
  l2_body.load_from_ref(deref_i, i, rgn_i);
  l2_body.assume(deref_j <= 3);
  l2_body.assertion(deref_i <= deref_j + 3);
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

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  cfg->simplify();
  crab::outs() << *cfg << "\n";

  {
    z_obj_zones_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }
  {
    z_obj_zones_t init;
    run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
  }

  delete cfg;
  return 0;
}