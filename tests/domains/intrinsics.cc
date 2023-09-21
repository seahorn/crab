#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("x0", "ret");
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb1;
  bb1_f >> ret;
  // adding statements
  x0.assign(k, 2147483648);
  x0.intrinsic("foo1", {i}, {i, k});
  x0.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb1_t.havoc(nd);
  bb1_t.select(inc, nd, 1, 2);
  bb1_t.add(i, i, inc);
  bb1_t.intrinsic("foo2", {i}, {i, k});

  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  /* Show that a CFG can have intrinsics */
  variable_factory_t vfac;
  z_cfg_t *cfg = prog(vfac);
  crab::outs() << *cfg << "\n";
  delete cfg;

  /* For now all domains ignore intrinsics so nothing is printed */
  z_var i1(vfac["i1"], crab::INT_TYPE, 32);
  z_var i2(vfac["i2"], crab::INT_TYPE, 32);
  z_var o1(vfac["o1"], crab::INT_TYPE, 32);
  z_var o2(vfac["o2"], crab::INT_TYPE, 32);

  {
    z_interval_domain_t inv;
    inv.intrinsic("foo1", {o1, o2}, {i1, i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";
  }

  {
    z_term_domain_t inv;
    inv.intrinsic("foo1", {o1, o2}, {i1, i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";
  }

  {
    z_bool_num_domain_t inv;
    inv.intrinsic("foo1", {o1, o2}, {i1, i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";
  }

  {
    reduced_domain_product2<z_number, varname_t, z_sdbm_domain_t,
			    z_dis_interval_domain_t> inv;      
    inv.intrinsic("foo1", {o1, o2}, {i1, i2});
    inv.intrinsic("foo2", {o1}, {i1});
    inv.intrinsic("foo2", {}, {i1});
    inv.intrinsic("foo2", {}, {});
    crab::outs() << inv << "\n";
  }
// #ifdef HAVE_APRON
//   {
//     z_pk_apron_domain_t inv;
//     inv.intrinsic("foo1", {o1, o2}, {i1, i2});
//     inv.intrinsic("foo2", {o1}, {i1});
//     inv.intrinsic("foo2", {}, {i1});
//     inv.intrinsic("foo2", {}, {});
//     crab::outs() << inv << "\n";
//   }
// #endif
// #ifdef HAVE_ELINA
//   {
//     z_pk_elina_domain_t inv;
//     inv.intrinsic("foo1", {o1, o2}, {i1, i2});
//     inv.intrinsic("foo2", {o1}, {i1});
//     inv.intrinsic("foo2", {}, {i1});
//     inv.intrinsic("foo2", {}, {});
//     crab::outs() << inv << "\n";
//   }
// #endif

  return 0;
}
