#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main(int argc, char **argv) {
#ifdef HAVE_APRON  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  z_var M(vfac["M"], crab::ARR_INT_TYPE);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  { // backward array load
    z_aa_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.domain_name() << "\n";
    pre += (y >= 1);
    pre.backward_array_load(y, M, 4, 4, inv);
    crab::outs() << "EXPECTED: {M[4..7] >= 1} \n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";

  { // backward array store
    z_aa_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.domain_name() << "\n";
    pre += (x >= 2);
    pre += (y >= 1);
    /// pre += (z >= -10);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, z, false, inv);
    crab::outs() << "EXPECTED: {z >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";

  { // backward array store
    z_aa_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.domain_name() << "\n";
    pre += (x >= 2);
    pre += (y >= 1);
    pre += (z == -10);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, z, false, inv);
    crab::outs() << "EXPECTED: _|_\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";

  { // backward array store
    z_aa_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.domain_name() << "\n";
    pre += (x >= 2);
    pre += (y >= 1);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, y, false, inv);
    crab::outs() << "EXPECTED: {y >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
#endif


  return 0;
}
