#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/domains/value_partitioning_domain.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

using domain_t = product_value_partitioning_domain<z_sdbm_domain_t>;

#define START_PARTITION "value_partition_start"
#define END_PARTITION "value_partition_end"

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  {
    domain_t dom0;
    crab::outs() << "===== Initially ======\n" << dom0 << "\n";    
    dom0.intrinsic(START_PARTITION, {x}, {});
    crab::outs() << "Start partitioning on " << x << "\n";
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(x, z_number(5));
    dom1.assign(y, z_number(10));
    crab::outs() << "After x:=5;y:=10 " << dom1 << "\n";
    dom2.assign(x, z_number(10));
    dom2.assign(y, z_number(20));
    crab::outs() << "After x:=10;y:=20 " << dom2 << "\n";
    domain_t dom3 = dom1 | dom2;
    crab::outs() << "After join " << dom3 << "\n";
    dom3.intrinsic(START_PARTITION, {z}, {});
    crab::outs() << "Start partitioning on " << z << "\n";
    domain_t dom4(dom3);
    domain_t dom5(dom3);
    dom4.assign(z, z_number(5));
    dom4.assign(w, z_number(10));
    crab::outs() << "After z:=5;w:=10 " << dom4 << "\n";
    dom5.assign(z, z_number(10));
    dom5.assign(w, z_number(20));
    crab::outs() << "After z:=10;w:=20 " << dom5 << "\n";
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "After join " << dom6 << "\n";

    // zones can infer bottom without any partitioning
    // dom6 += (z_lin_exp_t(x) == 5);
    // dom6 += (z_lin_exp_t(y) == 20);    // bottom!

    // zones cannot infer bottom without partitioning
    // dom6 += (z_lin_exp_t(x) == 5);
    // dom6 += (z_lin_exp_t(y) == 11);    // bottom!

    dom6.intrinsic(END_PARTITION, {x}, {});
    crab::outs() << "End partitioning on " << x << "\n";
    crab::outs() << dom6 << "\n";
    dom6.intrinsic(END_PARTITION, {z}, {});
    crab::outs() << "End partitioning on " << z << "\n";
    crab::outs() << dom6 << "\n";
  }


  {
    domain_t dom0;
    crab::outs() << "===== Infer bottom ======\n" << dom0 << "\n";
    dom0.intrinsic(START_PARTITION, {x}, {});
    crab::outs() << "Start partitioning on " << x << "\n";
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(x, z_number(5));
    dom1.assign(y, z_number(10));
    crab::outs() << "After x:=5;y:=10 " << dom1 << "\n";
    dom2.assign(x, z_number(10));
    dom2.assign(y, z_number(20));
    crab::outs() << "After x:=10;y:=20 " << dom2 << "\n";
    domain_t dom3 = dom1 | dom2;
    crab::outs() << "After join " << dom3 << "\n";
    dom3.intrinsic(START_PARTITION, {z}, {});
    crab::outs() << "Start partitioning on " << z << "\n";
    domain_t dom4(dom3);
    domain_t dom5(dom3);
    dom4.assign(z, z_number(5));
    dom4.assign(w, z_number(10));
    crab::outs() << "After z:=5;w:=10 " << dom4 << "\n";
    dom5.assign(z, z_number(10));
    dom5.assign(w, z_number(20));
    crab::outs() << "After z:=10;w:=20 " << dom5 << "\n";
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "After join " << dom6 << "\n";

    // zones can infer bottom without any partitioning
    // dom6 += (z_lin_exp_t(x) == 5);
    // dom6 += (z_lin_exp_t(y) == 20);    // bottom!

    // zones cannot infer bottom without partitioning
    dom6 += (z_lin_exp_t(x) == 5);
    dom6 += (z_lin_exp_t(y) == 11);    // bottom!
    crab::outs() << dom6 << "\n";
  }

  z_var v1(vfac["v1"], crab::INT_TYPE, 32);
  z_var v2(vfac["v2"], crab::INT_TYPE, 32);
  z_var v3(vfac["v3"], crab::INT_TYPE, 32);
  z_var v4(vfac["v4"], crab::INT_TYPE, 32);
  z_var w1(vfac["w1"], crab::INT_TYPE, 32);
  z_var w2(vfac["w2"], crab::INT_TYPE, 32);
  z_var w3(vfac["w3"], crab::INT_TYPE, 32);
  z_var w4(vfac["w4"], crab::INT_TYPE, 32);
  
  {
    
    domain_t dom0;
    crab::outs() << "===== Join (EXPECTED Top) ======\n" ;
    dom0.intrinsic(START_PARTITION, {v1}, {});
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(v1, z_number(5));
    dom1.assign(v2, z_number(10));
    dom2.assign(v1, z_number(10));
    dom2.assign(v2, z_number(20));
    domain_t dom3 = dom1 | dom2;

    domain_t dom4(dom0);
    domain_t dom5(dom0);
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    dom4.assign(v3, z_number(50));
    dom4.assign(v4, z_number(60));
    dom5.assign(v3, z_number(52));
    dom5.assign(v4, z_number(62));
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "Join\n\t" << dom3 << "\n\t" << dom6 <<"\n";
    domain_t dom7 = dom3 | dom6;

    crab::outs() << "RES:"<< dom7 << "\n";
  }

  {
    
    domain_t dom0;
    crab::outs() << "===== Meet (EXPECTED non-bot/non-top ======\n" ;
    dom0.intrinsic(START_PARTITION, {v1}, {});
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(v1, z_number(5));
    dom1.assign(v2, z_number(10));
    dom2.assign(v1, z_number(10));
    dom2.assign(v2, z_number(20));
    domain_t dom3 = dom1 | dom2;

    domain_t dom4(dom0);
    domain_t dom5(dom0);
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    dom4.assign(v3, z_number(50));
    dom4.assign(v4, z_number(60));
    dom5.assign(v3, z_number(52));
    dom5.assign(v4, z_number(62));
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "Meet\n\t" << dom3 << "\n\t" << dom6 <<"\n";
    domain_t dom7 =  (dom3 & dom6);

    crab::outs() << "RES:"<< dom7 << "\n";
  }

  {
    
    domain_t dom0;
    crab::outs() << "===== Join EXPECTED non-top ======\n" ;
    dom0.intrinsic(START_PARTITION, {v1}, {});
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(v1, z_number(5));
    dom1.assign(v2, z_number(10));
    dom2.assign(v1, z_number(10));
    dom2.assign(v2, z_number(20));
    domain_t dom3 = dom1 | dom2;

    domain_t dom4(dom0);
    domain_t dom5(dom0);
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    dom4.assign(v1, z_number(50));
    dom4.assign(v2, z_number(60));
    dom5.assign(v1, z_number(52));
    dom5.assign(v2, z_number(62));
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "Join\n\t" << dom3 << "\n\t" << dom6 <<"\n";
    domain_t dom7 = dom3 | dom6;

    crab::outs() << "RES:"<<dom7 << "\n";
  }

  {
    
    domain_t dom0;
    crab::outs() << "===== Meet EXPECTED bot ======\n" ;
    dom0.intrinsic(START_PARTITION, {v1}, {});
    domain_t dom1(dom0);
    domain_t dom2(dom0);
    dom1.assign(v1, z_number(5));
    dom1.assign(v2, z_number(10));
    dom2.assign(v1, z_number(10));
    dom2.assign(v2, z_number(20));
    domain_t dom3 = dom1 | dom2;

    domain_t dom4(dom0);
    domain_t dom5(dom0);
    dom0.intrinsic(START_PARTITION, {v3}, {});    
    dom4.assign(v1, z_number(50));
    dom4.assign(v2, z_number(60));
    dom5.assign(v1, z_number(52));
    dom5.assign(v2, z_number(62));
    domain_t dom6 = dom4 | dom5;
    crab::outs() << "Meet\n\t" << dom3 << "\n\t" << dom6 <<"\n";
    domain_t dom7 = dom3 & dom6;

    crab::outs() << "RES:"<<dom7 << "\n";
  }
  
  return 0;
}
