/* This test case does not show a bug but an imprecision of the
   widening operator with rationals */

#include "../common.hpp"
using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main(int argc, char **argv) {

#ifdef HAVE_LDD
  variable_factory_t vfac;
  {

    // build first widening argument
    q_lin_cst_sys_t csts1, csts2, csts3;

    csts1 += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(2));
    csts1 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(2));
    q_boxes_domain_t inv1_c1;
    inv1_c1 += csts1;

    csts2 += (q_var(vfac["x"], crab::REAL_TYPE) == q_number(1));
    csts2 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(1));
    q_boxes_domain_t inv1_c2;
    inv1_c2 += csts2;

    csts3 += (q_var(vfac["x"], crab::REAL_TYPE) == q_number(1));
    csts3 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(0));
    q_boxes_domain_t inv1_c3;
    inv1_c3 += csts3;

    q_boxes_domain_t inv1;
    inv1.set_to_bottom();
    inv1 |= inv1_c1;
    inv1 |= inv1_c2;
    inv1 |= inv1_c3;

    // build first widening argument
    q_lin_cst_sys_t csts4, csts5, csts6, csts7, csts8, csts9, csts10, csts11,
        csts12;

    csts4 += (q_var(vfac["x"], crab::REAL_TYPE) > q_number(3));
    csts4 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(3));
    q_boxes_domain_t inv2_c1;
    inv2_c1 += csts4;

    csts5 += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(2));
    csts5 += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(3));
    csts5 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(3));
    q_boxes_domain_t inv2_c2;
    inv2_c2 += csts5;

    csts6 += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(2));
    csts6 += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(3));
    csts6 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(2));
    q_boxes_domain_t inv2_c3;
    inv2_c3 += csts6;

    csts7 += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(2));
    csts7 += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(3));
    csts7 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(1));
    q_boxes_domain_t inv2_c4;
    inv2_c4 += csts7;

    csts8 += (q_var(vfac["x"], crab::REAL_TYPE) > q_number(1));
    csts8 += (q_var(vfac["x"], crab::REAL_TYPE) < q_number(2));
    csts8 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(2));
    q_boxes_domain_t inv2_c5;
    inv2_c5 += csts8;

    csts9 += (q_var(vfac["x"], crab::REAL_TYPE) > q_number(1));
    csts9 += (q_var(vfac["x"], crab::REAL_TYPE) < q_number(2));
    csts9 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(1));
    q_boxes_domain_t inv2_c6;
    inv2_c6 += csts9;

    csts10 += (q_var(vfac["x"], crab::REAL_TYPE) == q_number(1));
    csts10 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(2));
    q_boxes_domain_t inv2_c7;
    inv2_c7 += csts10;

    csts11 += (q_var(vfac["x"], crab::REAL_TYPE) == q_number(1));
    csts11 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(1));
    q_boxes_domain_t inv2_c8;
    inv2_c8 += csts11;

    csts12 += (q_var(vfac["x"], crab::REAL_TYPE) == q_number(1));
    csts12 += (q_var(vfac["y"], crab::REAL_TYPE) == q_number(0));
    q_boxes_domain_t inv2_c9;
    inv2_c9 += csts12;

    q_boxes_domain_t inv2;
    inv2.set_to_bottom();
    inv2 |= inv2_c1;
    inv2 |= inv2_c2;
    inv2 |= inv2_c3;
    inv2 |= inv2_c4;
    inv2 |= inv2_c5;
    inv2 |= inv2_c6;
    inv2 |= inv2_c7;
    inv2 |= inv2_c8;
    inv2 |= inv2_c9;

    {
      crab::outs() << "OP1=" << inv1 << "\n";
      crab::outs() << "OP2=" << inv2 << "\n";
      q_boxes_domain_t inv3 = inv1 || inv2;
      crab::outs() << "WIDENING(OP1,OP2)=" << inv3 << "\n";
      q_lin_cst_sys_t c3 = inv3.to_linear_constraint_system();
      crab::outs() << "[A] CONVEX(WIDENING(OP1,OP2))=" << c3 << "\n";
      crab::outs() << "--------------------------------\n";
    }

    {
      q_lin_cst_sys_t c1 = inv1.to_linear_constraint_system();
      q_lin_cst_sys_t c2 = inv2.to_linear_constraint_system();
      q_interval_domain_t intv1, intv2;
      intv1 += c1;
      intv2 += c2;
      crab::outs() << "CONVEX(OP1)=" << intv1 << "\n";
      crab::outs() << "CONVEX(OP2)=" << intv2 << "\n";
      q_interval_domain_t intv3 = intv1 || intv2;
      crab::outs() << "[B] WIDENING(CONVEX(OP1),CONVEX(OP2))=" << intv3 << "\n";
      crab::outs() << "--------------------------------\n";
    }
    crab::outs() << "[A] and [B] should be the same if CONVEX distributes over "
                    "WIDENING \n";
  }
#endif
  return 0;
}
