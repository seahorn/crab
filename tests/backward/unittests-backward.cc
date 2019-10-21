#include "../program_options.hpp"
#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

int main(int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  z_var M(vfac["M"], crab::ARR_INT_TYPE);
  z_var x(vfac["x"], crab::INT_TYPE, 32);  
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  #ifdef HAVE_APRON
  { // backward array load
    z_ae_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";    
    pre += (y >= 1);
    pre.backward_array_load(y, M, 4, 4, inv);
    crab::outs() << "EXPECTED: {M[4..7] >= 1} \n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
  
  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
    pre += (x >= 2);
    pre += (y >= 1);
    ///pre += (z >= -10);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, z, false, inv);
    crab::outs() << "EXPECTED: {z >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
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
    z_ae_box_apron_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";            
    pre += (x >= 2);
    pre += (y >= 1);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, y, false, inv);
    crab::outs() << "EXPECTED: {y >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
  #endif 

  #ifdef HAVE_ELINA
  #if 0
  { // backward array load
    z_ae_zones_elina_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";    
    pre += (y >= 1);
    pre.backward_array_load(y, M, 4, 4, inv);
    crab::outs() << "EXPECTED: {M[4..7] >= 1} \n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
  
  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_zones_elina_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
    pre += (x >= 2);
    pre += (y >= 1);
    ///pre += (z >= -10);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, z, false, inv);
    crab::outs() << "EXPECTED: {z >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_zones_elina_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
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
    z_ae_zones_elina_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";            
    pre += (x >= 2);
    pre += (y >= 1);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, y, false, inv);
    crab::outs() << "EXPECTED: {y >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
  #endif 
  #endif  
  
  crab::outs() << "============================\n";
  
  { // backward array load
    z_ae_int_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";    
    pre += (y >= 1);
    pre.backward_array_load(y, M, 4, 4, inv);
    crab::outs() << "EXPECTED: {M[4..7] >= 1} \n";
    crab::outs() << "RESULT: " << pre << "\n";
  }
  
  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_int_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
    pre += (x >= 2);
    pre += (y >= 1);
    ///pre += (z >= -10);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, z, false, inv);
    crab::outs() << "EXPECTED: {z >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  crab::outs() << "============================\n";
  
  { // backward array store
    z_ae_int_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";        
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
    z_ae_int_t pre, inv;
    crab::outs() << "Test using " << pre.getDomainName() << "\n";            
    pre += (x >= 2);
    pre += (y >= 1);
    pre.array_store(M, 4, 4, x, false);
    pre.backward_array_store(M, 4, 4, y, false, inv);
    crab::outs() << "EXPECTED: {y >= 2, ...}\n";
    crab::outs() << "RESULT: " << pre << "\n";
  }

  return 0;
}
