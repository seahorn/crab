#include "../../common.hpp"
#include "../../program_options.hpp"

using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

/* Example of how to build a CFG */
z_cfg_t *prog1(variable_factory_t &vfac) {

  /*
    int i,x;
    int N = nd_int();
    __CRAB_assume(N > 0);
    i = 0;
    x = 0;
    y = 0;
    while (i < N) {
      // loop_counter(i);
      i++;
      x = x + 4;
      y = y + 8;
    }
    __CRAB_assert(x == 4*N ); // needs cross-coefficient saturation
    __CRAB_assert(y == 8*N ); // needs cross-coefficient saturation
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var x_next(vfac["x.next"], crab::INT_TYPE, 32);
  z_var n(vfac["N"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_header;
  loop_exit >> exit;
  // adding statements
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  entry.assign(y, 0);
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x == 4 * n);
  loop_exit.assertion(y == 8 * n);
  loop_body.intrinsic("loop_counter", {}, {i});
  loop_body.add(i, i, 1);
  loop_body.add(x, x, 4);
  loop_body.add(y, y, 8);

  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {
  /*
    int i,x;
    int N = nd_int();
    __CRAB_assume(N > 0);
    i = 0;
    x = 0;
    while (i < N) {
      i++;
      if (*) {
       x = x+2;
      } else {
       x = x+3;
      }
    }
    __CRAB_assert(x >= 2*N); // needs cross-coefficient saturation
    __CRAB_assert(x <= 3*N); // needs cross-coefficient saturation
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_body_then = cfg->insert("loop_body_then");
  z_basic_block_t &loop_body_else = cfg->insert("loop_body_else");
  z_basic_block_t &loop_body_tail = cfg->insert("loop_body_tail");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_body_then;
  loop_body >> loop_body_else;
  loop_body_then >> loop_body_tail;
  loop_body_else >> loop_body_tail;
  loop_body_tail >> loop_header;
  loop_exit >> exit;
  // adding statements
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x >= 2 * n);
  loop_exit.assertion(x <= 3 * n);
  loop_body.intrinsic("loop_counter", {}, {i});
  loop_body.add(i, i, 1);
  loop_body_then.add(x, x, 2);
  loop_body_else.add(x, x, 3);
  return cfg;
}

z_cfg_t *prog3(variable_factory_t &vfac) {
  /*
    int isz = 4;
    int len = nd_int();
    __CRAB_assume(len > 0);
    __CRAB_assume(len <= 10);
    int tsz = nd_int();
    __CRAB_assume(tsz >= len * isz);
    int i = 0;
    while (i < len) {
      int idx = i * isz;
      int offset = idx + isz;
      __CRAB_assert(offset <= tsz);
      i++;
    }
    __CRAB_assert(tsz >= len * isz);
   */
  // Defining program variables
  z_var isz(vfac["isz"], crab::INT_TYPE, 32);
  z_var len(vfac["len"], crab::INT_TYPE, 32);
  z_var tsz(vfac["tsz"], crab::INT_TYPE, 32);
  z_var tmp(vfac["tmp"], crab::INT_TYPE, 32);
  z_var tmp2(vfac["tmp2"], crab::INT_TYPE, 32);
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var idx(vfac["idx"], crab::INT_TYPE, 32);
  z_var offset(vfac["offset"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_body;
  loop_body >> loop_header;
  loop_header >> loop_exit;
  loop_exit >> exit;
  // adding statements
  entry.assign(isz, 4);
  entry.havoc(len);
  entry.assume(len >= 1);
  entry.assume(len <= 10);
  entry.havoc(tsz);
  entry.mul(tmp, len, isz);
  entry.assume(z_lin_exp_t(tsz) >= tmp); // tmp <= tsz, tmp = 4 * len
  entry.assign(i, 0);
  loop_body.assume(z_lin_exp_t(i) < len); // i - len <= -1
  loop_body.mul(idx, i, isz);             // idx - 4i <= 0 => idx - 4len <= -1
  loop_body.add(offset, idx, isz);        // offset = idx + 4
  loop_body.assertion(z_lin_exp_t(offset) <= tsz);
  loop_body.add(i, i, 1);
  loop_exit.assume(z_lin_exp_t(i) >= len);
  loop_exit.mul(tmp2, len, isz);
  loop_exit.assertion(z_lin_exp_t(tsz) >= tmp2);
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  auto &coeffs = crab_domain_params_man::get().coefficients();
  coeffs.insert(coeffs.end(), {2, 3, 4});

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    crab::outs() << "\n========== prog1 ==========\n";
    crab::outs() << *cfg << "\n";
    crab::outs() << "---------- TVPI(DBM) ----------\n";
    {
      z_tvpi_dbm_domain_t init;
      run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    }
    crab::outs() << "---------- TVPI(SplitDBM) ----------\n";
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac);
    crab::outs() << "\n========== prog2 ==========\n";
    crab::outs() << *cfg << "\n";
    crab::outs() << "---------- TVPI(DBM) ----------\n";
    {
      z_tvpi_dbm_domain_t init;
      run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    }
    crab::outs() << "---------- TVPI(SplitDBM) ----------\n";
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog3(vfac);
    crab::outs() << "\n========== prog3 ==========\n";
    crab::outs() << *cfg << "\n";
    crab::outs() << "---------- TVPI(DBM) ----------\n";
    {
      z_tvpi_dbm_domain_t init;
      run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    }
    crab::outs() << "---------- TVPI(SplitDBM) ----------\n";
    delete cfg;
  }

  return 0;
}
