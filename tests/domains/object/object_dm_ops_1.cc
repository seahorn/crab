#include "../../common.hpp"
#include "../../program_options.hpp"
#include <crab/domains/object_domain.hpp>
#include <crab/domains/uf_domain.hpp>

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::object_domain_impl;
using namespace crab::domains::term;

namespace {
using variable_vector_t = std::vector<z_var>;
using variable_or_constant = z_var_or_cst_t;
using variable_or_constant_vector_t = std::vector<z_var_or_cst_t>;
} // namespace

template <class EUFDom> term_operator_t make_fresh_symbol(EUFDom &dom) {
  term_operator_t t = dom.make_uf(term_op_val_generator_t::get_next_val());
  return t;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_var A(vfac["objA"], crab::REF_TYPE, 32);
  z_var A2(vfac["objA2"], crab::REF_TYPE, 32);
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
  z_var t(vfac["t"], crab::REF_TYPE, 32);
  z_var p2(vfac["p2"], crab::REF_TYPE, 32);
  z_var r2(vfac["r2"], crab::REF_TYPE, 32);
  z_var rgn1(vfac["V_a"], crab::REG_INT_TYPE, 32);
  z_var rgn2(vfac["V_b"], crab::REG_INT_TYPE, 32);
  z_var rgn3(vfac["V_c"], crab::REG_INT_TYPE, 32);

  z_var rgn4(vfac["V_d"], crab::REG_INT_TYPE, 32);
  z_var rgn5(vfac["V_e"], crab::REG_INT_TYPE, 32);

  // Create allocation sites
  crab::tag_manager as_man;

  // mimic dsa node info
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
  z_var_or_cst_t four32(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t five32(z_number(5), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t six32(z_number(6), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t eight32(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t ten32(z_number(10), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t twelve32(z_number(12),
                          crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size4(z_number(4), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32));
  z_var_or_cst_t size12(z_number(12), crab::variable_type(crab::INT_TYPE, 32));

  z_obj_zones_t inv1;
  inv1.region_init(rgn1);
  inv1.region_init(rgn2);
  inv1.region_init(rgn3);
  inv1.region_init(rgn4);
  inv1.region_init(rgn5);
  inv1.intrinsic("regions_from_memory_object", obj1, {});
  inv1.intrinsic("regions_from_memory_object", obj2, {});

  // Create one singleton object
  inv1.ref_make(A, rgn1, size12, as_man.mk_tag());
  inv1 += (A > 0);
  // perform strong updates
  inv1.ref_gep(A, rgn1, p, rgn1, z_number(0)); // V_a, p = gep_ref(V_a, &A, 0)
  inv1.ref_gep(p, rgn1, r, rgn2, z_number(4)); // V_b, r = gep_ref(V_a, p, 4)
  inv1.ref_store(p, rgn1, zero32);             // store_ref(V_a, p, 0)
  inv1.ref_store(r, rgn2, one32);              // store_ref(V_b, r, 1)
  inv1.ref_load(p, rgn1, x);                   // x = load_ref(V_a, p)
  inv1.ref_load(r, rgn2, y);                   // y = load_ref(V_b, r)
  /* inv1:
  base: { 0 < objA && p = objA && r = p + 4 && V_a = 0 && x = 0 && V_b = 1 && y
  = 1 }, addrs: { objA == p && objA == q }, uf_regs: {}, odi_map : empty
  */
  {
    crab::outs()
        << "=== forget & project 1: state contains singleton objects ===\n";

    crab::outs() << "forget " << x << ", " << rgn2 << " at \n"
                 << inv1 << "\n=> \n";
    z_obj_zones_t inv2 = inv1;
    inv2.forget({x, rgn2});
    crab::outs() << inv2 << "\n";
    crab::outs() << "project " << x << ", " << rgn2 << " at \n"
                 << inv1 << "\n=> \n";
    inv2 = inv1;
    inv2.project({x, rgn2});
    crab::outs() << inv2 << "\n";
  }

  // Object moves into odi map
  z_obj_zones_t inv3 = inv1;
  inv3.ref_make(A2, rgn1, size12, as_man.mk_tag());
  inv3.ref_gep(A2, rgn1, s, rgn1, z_number(0)); // V_a, s = gep_ref(V_a, &A, 0)
  inv3.ref_store(s, rgn1, six32);               // store_ref(V_a, s, 6)
  inv3.ref_gep(s, rgn1, t, rgn2, z_number(4));  // V_b, t = gep_ref(V_a, s, 4)
  inv3.ref_store(t, rgn2, twelve32);            // store_ref(V_b, t, 12)
  /* inv3:
  base: { 0 < objA && p = objA && r = p + 4 && objA2 > 0
          && s = objA2 && t = s + 4 && x = 0 && y = 1 },
  addrs: { objA == p && objA == q && objA2 == s && s == t && mru_A == t },
  uf_regs: {},
  odi_map:
    objA(
      sum: { V_a = 0 && V_b = 1 },
      cache: { V_a = 6 && V_b = 12 },
      uf_flds: {}
    )
  */
  {
    crab::outs() << "=== forget & project 2: object uses odi map ===\n";

    crab::outs() << "forget " << t << ", " << rgn2 << " at \n"
                 << inv3 << "\n=> \n";
    z_obj_zones_t inv2 = inv3;
    inv2.forget({x, rgn2});
    crab::outs() << inv2 << "\n";
    crab::outs() << "project " << y << ", " << rgn1 << " at \n"
                 << inv3 << "\n=> \n";
    inv2 = inv3;
    inv2.project({y, rgn1});
    crab::outs() << inv2 << "\n";
  }

  // Use new mru object
  z_obj_zones_t inv5 = inv3;
  inv5.ref_store(p, rgn1, four32);  // store_ref(V_a, p, 4)
  inv5.ref_store(r, rgn2, eight32); // store_ref(V_b, r, 8)
  inv5.ref_load(p, rgn1, x);        // x = load_ref(V_a, p)
  inv5.ref_load(r, rgn2, y);        // y = load_ref(V_b, r)
  /* inv4:
  base: { 0 < objA && p = objA && r = p + 4 && objA2 > 0
          && s = objA2 && t = s + 4 },
  addrs: { objA == p && objA == q && objA2 == s && s == t && mru_A == p },
  uf_regs: { x = $symb0,  y = $symb1 },
  odi_map:
    objA(
      sum: { 0 <= V_a && V_a <= 6 && 1 <= V_b && V_b <= 12 && V_a < V_b },
      cache: { V_a = 4 && V_b = 8 },
      uf_flds: { V_a = $symb0,  V_b = $symb1 }
    )
  */
  {
    crab::outs() << "=== forget & project 3: object uses odi map, reg dom "
                    "contains equalities ===\n";

    crab::outs() << "forget " << t << ", " << x << ", " << rgn2 << " at \n"
                 << inv5 << "\n=> \n";
    z_obj_zones_t inv6 = inv5;
    inv6.forget({t, x, rgn2});
    crab::outs() << inv6 << "\n";
    crab::outs() << "project " << s << ", " << y << ", " << rgn1 << " at \n"
                 << inv5 << "\n=> \n";
    inv6 = inv5;
    inv6.project({s, y, rgn1});
    crab::outs() << inv6 << "\n";
  }

  {
    crab::outs() << "=== rename 1: state contains singleton objects ===\n";
    crab::outs() << "rename " << x << ", " << rgn2 << ", " << p << " -> to "
                 << z << ", " << rgn3 << ", " << s << " at \n"
                 << inv1 << "=> \n";
    z_obj_zones_t inv2 = inv1;
    inv2.rename({x, rgn2, p}, {z, rgn3, s});
    crab::outs() << inv2 << "\n";
  }

  {
    crab::outs() << "=== rename 2: object uses odi map, reg dom contains "
                    "equalities ===\n";
    crab::outs() << "rename " << x << ", " << rgn2 << ", " << p << " -> to "
                 << z << ", " << rgn3 << ", " << p2 << " at \n"
                 << inv5 << "=> \n";
    z_obj_zones_t inv6 = inv5;
    inv6.rename({x, rgn2, p}, {z, rgn3, p2});
    crab::outs() << inv6 << "\n";
  }
  return 0;
}