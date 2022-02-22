#include "../../common.hpp"
#include "../../program_options.hpp"
#include <crab/domains/object_domain.hpp>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::object_domain_impl;
using namespace ikos;

namespace {
  using variable_vector_t = std::vector<z_var>;
  using variable_or_constant = z_var_or_cst_t;
  using variable_or_constant_vector_t = std::vector<z_var_or_cst_t>;
}

void simulate_intrinsic(z_obj_zones_t &dom, variable_vector_t rgn_vars) {
  variable_or_constant_vector_t inputs;
  inputs.reserve(rgn_vars.size());
  for (auto rgn : rgn_vars) {
    variable_or_constant v(rgn);
    inputs.push_back(v);
  }
  dom.intrinsic("regions_from_memory_object", inputs, {});
}

int main(int argc, char **argv) {

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  
  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));

  variable_vector_t obj1;
  variable_vector_t obj2;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var ref(vfac["ref"], crab::REF_TYPE, 32);

  z_var rgn1(vfac["V_a"], crab::REG_INT_TYPE, 32);
  obj1.push_back(rgn1);
  z_var rgn2(vfac["V_b"], crab::REG_INT_TYPE, 32);
  obj1.push_back(rgn2);
  z_var rgn3(vfac["V_c"], crab::REG_INT_TYPE, 32);
  obj1.push_back(rgn3);

  z_var rgn4(vfac["V_d"], crab::REG_INT_TYPE, 32);
  obj2.push_back(rgn4);
  z_var rgn5(vfac["V_e"], crab::REG_INT_TYPE, 32);
  obj2.push_back(rgn5);
  
  { // join
    crab::outs() << "=== 1. join with two objects === \n";
    z_obj_zones_t left;
    z_obj_zones_t right;
    
    simulate_intrinsic(left, obj1);
    simulate_intrinsic(right, obj1);

    simulate_intrinsic(left, obj2);
    simulate_intrinsic(right, obj2);

    left += (x <= z_number(0));
    left.region_assign(rgn1, z_number(0));
    left.region_assign(rgn2, z_number(1));
    left.region_assign(rgn3, z_number(1));

    right.assign(x, z_number(2));
    right.add_cons_for_objects(rgn1 > rgn2);
    right.region_assign(rgn4, z_number(1));

    z_obj_zones_t l_join_r = left | right;
    crab::outs() << left << " | \n" << right << " = \n" << l_join_r << "\n";
  }

  { // join 2
    crab::outs() << "=== 2. join with two objects === \n";
    z_obj_zones_t left;
    z_obj_zones_t right;
    
    simulate_intrinsic(left, obj1);
    simulate_intrinsic(right, obj1);

    simulate_intrinsic(left, obj2);
    simulate_intrinsic(right, obj2);

    left += (y <= z_number(5));
    left += (x <= y);
    left.region_assign(rgn1, z_number(0));
    left.region_assign(rgn3, z_number(3));
    left.add_cons_for_objects(rgn2 <= rgn3);

    left.region_assign(rgn4, z_number(0));
    left.add_cons_for_objects(rgn4 <= rgn5);

    right += (y <= z_number(7));
    right.add_cons_for_objects(rgn1 > rgn2);
    right.region_assign(rgn4, z_number(4));

    z_obj_zones_t l_join_r = left | right;
    crab::outs() << left << " | \n" << right << " = \n" << l_join_r << "\n";
  }

  {// meet 1
    crab::outs() << "=== 3. meet with two objects === \n";
    z_obj_zones_t left;
    z_obj_zones_t right;

    simulate_intrinsic(left, obj1);
    simulate_intrinsic(right, obj1);

    simulate_intrinsic(left, obj2);
    simulate_intrinsic(right, obj2);

    left += (y >= z_number(5));
    left.region_assign(rgn1, z_number(0));

    right += (y <= z_number(7));
    right.region_assign(rgn4, z_number(4));

    z_obj_zones_t l_meet_r = left & right;
    crab::outs() << left << " & \n" << right << " = \n" << l_meet_r << "\n";
  }

  return 0;
}
