#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
#ifdef HAVE_APRON
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  crab::outs() << "Running APRON\n\n";
  variable_factory_t vfac;
  {
    q_pk_apron_domain_t inv;
    q_lin_cst_sys_t csts;
    csts += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(1));
    csts += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(1));
    crab::outs() << "INITIALLY=" << inv << "\n";
    crab::outs() << "ADDING CONSTRAINTS=" << csts << "\n";
    inv += csts;
    crab::outs() << "EXPECTED={x = 1}\n";
    crab::outs() << "RESULT=" << inv << "\n";
  }
  crab::outs() << "------------------------------------\n";
  {
    q_pk_apron_domain_t inv;
    q_lin_cst_sys_t csts;
    csts += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(0.5));
    csts += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(0.5));
    crab::outs() << "INITIALLY=" << inv << "\n";
    crab::outs() << "ADDING CONSTRAINTS=" << csts << "\n";
    inv += csts;
    crab::outs() << "EXPECTED={x = 1/2}\n";
    crab::outs() << "RESULT=" << inv << "\n";
  }
  crab::outs() << "------------------------------------\n";
  {
    q_oct_apron_domain_t inv;
    q_lin_cst_sys_t csts;
    csts += (q_var(vfac["x"], crab::REAL_TYPE) >= q_number(0.5));
    csts += (q_var(vfac["x"], crab::REAL_TYPE) <= q_number(0.5));
    crab::outs() << "INITIALLY=" << inv << "\n";
    crab::outs() << "ADDING CONSTRAINTS=" << csts << "\n";
    inv += csts;
    crab::outs() << "EXPECTED={x = 1/2}\n";
    crab::outs() << "RESULT=" << inv << "\n";
  }
  crab::outs() << "------------------------------------\n";

#endif


  return 0;
}
