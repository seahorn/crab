#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/domains/symbolic_variable_eq_domain.hpp>
#include <crab/domains/intervals.hpp>

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

using z_interval_domain_t = interval_domain<z_number, varname_t>;
using test_domain_t =
    crab::domains::symbolic_variable_equiality_domain<z_interval_domain_t>;
using value_domain_t = symbolic_variable_equiality_domain_impl::symbolic_var;

int i = 0;

void perfrom_domain_operations(const test_domain_t& dom1, const test_domain_t &dom2) {
  crab::outs() << "---- case "<<i<<"---- \n";
  crab::outs() << "Dom1=" << dom1 << "\n";
  crab::outs() << "Dom2=" << dom2 << "\n";
  bool r1 = dom1 <= dom2;
  crab::outs() << "Dom1 <= Dom2 = " << (r1 ? "true" : "false") << "\n";
  bool r2 = dom2 <= dom1;
  crab::outs() << "Dom2 <= Dom1 = " << (r2 ? "true" : "false") << "\n";
  test_domain_t dom3 = dom1 | dom2;
  crab::outs() << "Dom3 = Dom1 | Dom2 = " << dom3 << "\n";
  bool r3 = dom1 <= dom3;
  bool r4 = dom2 <= dom3;
  crab::outs() << "Dom1 <= Dom3 = " << (r3 ? "true" : "false") << "\n";
  crab::outs() << "Dom2 <= Dom3 = " << (r4 ? "true" : "false") << "\n";
  test_domain_t dom4 = dom1 & dom2;
  crab::outs() << "Dom4 = Dom1 & Dom2 = " << dom4 << "\n";
  bool r5 = dom4 <= dom1;
  bool r6 = dom4 <= dom2;
  crab::outs() << "Dom4 <= Dom1 = " << (r5 ? "true" : "false") << "\n";
  crab::outs() << "Dom4 <= Dom2 = " << (r6 ? "true" : "false") << "\n";
  i ++;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_var v1(vfac["v1"], crab::INT_TYPE, 32);
  z_var v2(vfac["v2"], crab::INT_TYPE, 32);
  z_var v3(vfac["v3"], crab::INT_TYPE, 32);
  z_var v4(vfac["v4"], crab::INT_TYPE, 32);
  z_var v5(vfac["v5"], crab::INT_TYPE, 32);
  z_var v6(vfac["v6"], crab::INT_TYPE, 32);
  z_var v7(vfac["v7"], crab::INT_TYPE, 32);
  z_var v8(vfac["v8"], crab::INT_TYPE, 32);
  z_var v9(vfac["v9"], crab::INT_TYPE, 32);
  z_var v10(vfac["v10"], crab::INT_TYPE, 32);
  z_var v11(vfac["v11"], crab::INT_TYPE, 32);
  z_var v12(vfac["v12"], crab::INT_TYPE, 32);

  {
    value_domain_t idom1(1),idom2(2),idom3(3); // some abstract domain values
    test_domain_t dom1, dom2;
    // dom1 : {v1,v2}=>#var1
    // dom2 : {v4,v5}=>#var2
    dom1.set(v1, idom1);
    dom1.add(v1, v2);

    dom2.set(v4, idom2);
    dom2.add(v4, v5);
    /*
    join: {}
    no relation between dom1 and dom2
    meet: {[v1,v2]=>#var0,[v4,v5]=>#var1}
    */
    perfrom_domain_operations(dom1, dom2);

    test_domain_t dom4, dom5;
    // dom4 : {v1,v2}=>#var1
    // dom5 : {v4,v2}=>#var2
    dom4.set(v1, idom1);
    dom4.add(v1, v2);

    dom5.set(v4, idom2);
    dom5.add(v4, v2);
    /*
    join: {}
    no relation between dom4 and dom5
    meet: {[v1,v2,v4]=>#var0}
    */
    perfrom_domain_operations(dom4, dom5);

    test_domain_t dom7, dom8;
    // dom7 : {v2,v3}=>#var1, {v1,v5}=>#var3
    // dom8 : {v1,v2}=>#var1

    dom7.set(v2, idom1);
    dom7.add(v2, v3);
    dom7.set(v5, idom1);
    dom7.add(v5, v1);

    dom8.set(v1, idom2);
    dom8.add(v1, v2);
    /*
    join: {}
    no relation between dom7 and dom8
    meet: {[v1,v2,v3,v5]=>#var0}
    */
    perfrom_domain_operations(dom7, dom8);
  }

  { // test all operations - level simple
    value_domain_t idom1(1),idom2(2); // some abstract domain values
    test_domain_t dom1, dom2;
    // dom1 : {v1,v2,v3,v4}=>#var1
    // dom2 : {v1,v2,v3}=>#var2

    dom1.set(v1, idom1);
    dom1.add(v1, v2);
    dom1.add(v2, v3);
    dom1.add(v3, v4);

    dom2.set(v2, idom2);
    dom2.add(v2, v1);
    dom2.add(v1, v3);
    /*
    join: {v1,v2,v3}=>#var0
    dom1 <= dom2
    meet: {[v1,v2,v3,v4]=>#var0}
    */
    perfrom_domain_operations(dom1, dom2);

    test_domain_t dom5(dom1);
    // dom 5: {[v1,v3,v4]=>#var1}
    dom5 -= v2;
    crab::outs() << "After forgetting " << v2 << " in Dom 1:" << dom5 << "\n";

    test_domain_t dom6(dom1);
    dom6.rename({v1,v2,v3,v4}, {v5,v6,v7,v8});
    crab::outs() << "After renaming {v1,v2,v3,v4} with {v5,v6,v7,v8} in Dom1:" << dom6 << "\n";

    test_domain_t dom7(dom1);
    // dom 7: {[v1,v3]=>#var1}
    dom7.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom1:" << dom7 << "\n";
  }

  { // test all operations - level moderate
    value_domain_t idom1(1),idom2(2), idom3(3); // some abstract domain values
    test_domain_t dom1, dom2, dom3;
    // dom1 : {v1,v2}=>#var2, {v3,v4}=>#var3
    // dom2 : {v2,v3}=>#var3, {v1,v4}=>#var1

    dom1.set(v1, idom2);
    dom1.add(v1, v2);

    dom1.set(v3, idom3);
    dom1.add(v3, v4);

    dom2.set(v2, idom3);
    dom2.add(v2, v3);
    dom2.set(v1, idom1);
    dom2.add(v1, v4);
    /*
    join: {}
    no relation between dom1 and dom2
    meet: {[v1,v2,v3,v4]=>#var0}
    */
    perfrom_domain_operations(dom1, dom2);
    test_domain_t dom5(dom1);
    // dom 5: {[v3,v4]=>#var3} if normalization
    // {[v1]=>#var2,[v3,v4]=>#var3} if not
    dom5 -= v2;
    crab::outs() << "After forgetting " << v2 << " in Dom 1:" << dom5 << "\n";

    test_domain_t dom6(dom1);
    dom6.rename({v1,v2,v3,v4}, {v5,v6,v7,v8});
    crab::outs() << "After renaming {v1,v2,v3,v4} with {v5,v6,v7,v8} in Dom1:" << dom6 << "\n";

    test_domain_t dom7(dom1);
    // dom 7: {} if normalization
    // {[v3]=>#var3,[v1]=>#var2}
    dom7.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom1:" << dom7 << "\n";
  }

  { // test all operations - level hard
    value_domain_t idom1(1),idom2(2), idom3(3); // some abstract domain values
    test_domain_t dom1, dom2;
    // dom1 : {v1,v2,v6,v8,v11}=>#var2, {v3,v7,v12}=>#var3, {v4,v5,v9,v10}=>#var1
    // dom2 : {v2,v4,v5,v10,v12}=>#var3, {v6,v7,v8,v9}=>#var1, {v1,v3,v11}=>#var2

    dom1.set(v6, idom2);
    dom1.add(v6, v1);
    dom1.add(v1, v2);
    dom1.add(v6, v8);
    dom1.add(v8, v11);

    dom1.set(v3, idom3);
    dom1.add(v3, v7);
    dom1.add(v7, v12);

    dom1.set(v4, idom1);
    dom1.add(v4, v9);
    dom1.add(v4, v10);
    dom1.add(v4, v5);

    dom2.set(v5, idom3);
    dom2.add(v5, v12);
    dom2.add(v5, v4);
    dom2.add(v4, v10);
    dom2.add(v12, v2);

    dom2.set(v6, idom1);
    dom2.add(v6, v8);
    dom2.add(v6, v7);
    dom2.add(v7, v9);

    dom2.set(v3, idom2);
    dom2.add(v3, v11);
    dom2.add(v11, v1);
    /*
    join: {[v6,v8]=>#var0,[v1,v11]=>#var1,[v4,v5,v10]=>#var2}
    no relation between dom1 and dom2
    meet: {[v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12]=>#var0}
    */
    perfrom_domain_operations(dom1, dom2);
    test_domain_t dom5(dom1);
    // dom 5: {[v1,v6,v8,v11]=>#var2,[v3,v7,v12]=>#var3,[v4,v5,v9,v10]=>#var1}
    dom5 -= v2;
    crab::outs() << "After forgetting " << v2 << " in Dom 1: " << dom5 << "\n";

    test_domain_t dom6(dom1);
    // dom6: {} if normalization
    // {[v3]=>#var3,[v1]=>#var2} if not
    dom6.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom1: " << dom6 << "\n";

    test_domain_t dom7(dom2);
    // dom7: {[v1,v3]=>#var2}
    dom7.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom2: " << dom7 << "\n";
  }

  {// test for all operations - level hard
    value_domain_t idom1(1),idom2(2), idom3(3); // some abstract domain values
    test_domain_t dom1, dom2;
    // dom1 : {v1,v2,v3,v4,v12}=>#var2, {v10,v11}=>#var3
    // dom2 : {v2,v3,v10}=>#var3, {v4,v12,v11}=>#var1
    dom1.set(v12, idom2);
    dom1.add(v12, v1);
    dom1.add(v12, v2);
    dom1.add(v12, v3);
    dom1.add(v12, v4);

    dom1.set(v11, idom3);
    dom1.add(v11, v10);

    dom2.set(v10, idom3);
    dom2.add(v10, v2);
    dom2.add(v10, v3);

    dom2.set(v11, idom3);
    dom2.add(v11, v4);
    dom2.add(v4, v12);

    /*
    join: {[v2,v3]=>#var0}
    no relation between dom1 and dom2
    meet: {[v1,v2,v3,v4,v10,v11,v12]=>#var0}
    */
    perfrom_domain_operations(dom1, dom2);
    test_domain_t dom5(dom1);
    dom5 -= v10;
    crab::outs() << "After forgetting " << v2 << " in Dom 1: " << dom5 << "\n";
    // dom5 = dom1.project(v10) := { {v1,v2,v3,v4,v12}=>#var2} } if normalization
    // dom5 := { {v1,v2,v3,v4,v12}=>#var2, {v11}=>#var3 }

    test_domain_t dom6(dom1);
    // dom6: {[v3,v1]=>#var0}
    dom6.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom1: " << dom6 << "\n";

    test_domain_t dom7(dom2);
    // dom7: {[v1,v3]=>#var2}
    // dom7: {} if normalization
    // {v3}=>#var3 if not
    dom7.project({v1,v3});
    crab::outs() << "After projecting on v1 and v3 in Dom2: " << dom7 << "\n";
  }

  {// test operation for object domain
    value_domain_t idom1(1),idom2(2), idom3(3); // some abstract domain values
    test_domain_t dom1, dom2;
    // dom1 : {v1, v2}=>#var1, {v3,v4}=>#var3
    // dom2 : {v1, v2}=>#var3
    dom1.set(v1, idom1);
    dom1.add(v1, v2);
    dom1.set(v3, idom3);
    dom1.add(v3, v4);

    dom2.set(v1, idom3);
    dom2.add(v1, v2);
    /*
    join: {[v1,v2]=>#var0}
    dom1 <= dom2
    meet: {[v1,v2]=>#var0,[v3,v4]=>#var1}
    */
    perfrom_domain_operations(dom1, dom2);
  }
}