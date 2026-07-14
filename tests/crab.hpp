#pragma once

#include <memory>

#include "./crab_lang.hpp"

template <typename Dom>
extern void z_intra_run(const crab::cfg_impl::z_cfg_t &cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, unsigned widening,
                        unsigned narrowing, unsigned jump_set_size,
                        bool enable_stats, bool enable_checker, bool print_invariants);

template <typename Dom>
extern void q_intra_run(const crab::cfg_impl::q_cfg_t &cfg,
                        crab::cfg_impl::basic_block_label_t entry, Dom init,
                        bool run_liveness, unsigned widening,
                        unsigned narrowing, unsigned jump_set_size,
                        bool enable_stats, bool enable_checker, bool print_invariants);

// To run abstract domains defined over integers
template <typename Dom>
void run(const crab::cfg_impl::z_cfg_t &cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  z_abs_domain_t init_wrapper(init);
  z_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#else
  z_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#endif
}

template <typename Dom>
void run_and_check(const crab::cfg_impl::z_cfg_t &cfg,
                   crab::cfg_impl::basic_block_label_t entry, Dom init,
                   bool run_liveness, unsigned widening, unsigned narrowing,
                   unsigned jump_set_size, bool enable_stats) {

#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  z_abs_domain_t init_wrapper(init);
  z_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, true, false);

#else
  z_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, true, false);

#endif
}

// To run abstract domains defined over rationals
template <typename Dom>
void run(const crab::cfg_impl::q_cfg_t &cfg,
         crab::cfg_impl::basic_block_label_t entry, Dom init, bool run_liveness,
         unsigned widening, unsigned narrowing, unsigned jump_set_size,
         bool enable_stats) {
#ifdef USE_GENERIC_WRAPPER
  using namespace crab::domain_impl;
  q_abs_domain_t init_wrapper(init);
  q_intra_run(cfg, entry, init_wrapper, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#else
  q_intra_run(cfg, entry, init, run_liveness, widening, narrowing,
              jump_set_size, enable_stats, false, true);
#endif
}

/**
 * Ergonomic convenience layer over the raw run()/run_and_check() above.
 *
 * The raw entry points take 8 positional arguments, four of which are the
 * same literals at almost every call site:
 *
 *   run(cfg, cfg.entry(), init, false, 1, 2, 20, stats_enabled);
 *
 * `run_config` bundles those fixpoint knobs with their usual defaults so the
 * common case reads as `run(cfg, init, stats)`, and lets a test override just
 * the one it cares about without respelling the rest:
 *
 *   run(cfg, init, stats);                                 // defaults
 *   run(cfg, init, stats, run_config().with_widening(10)); // one override
 *
 * These are plain additive overloads: the 8-argument run()/run_and_check()
 * calls elsewhere keep compiling and behaving exactly as before.
 */
struct run_config {
  bool run_liveness = false;
  unsigned widening = 1;
  unsigned narrowing = 2;
  unsigned jump_set_size = 20;

  // Fluent setters (C++14 has no designated initializers).
  run_config &with_liveness(bool v = true) { run_liveness = v; return *this; }
  run_config &with_widening(unsigned v) { widening = v; return *this; }
  run_config &with_narrowing(unsigned v) { narrowing = v; return *this; }
  run_config &with_jump_set_size(unsigned v) { jump_set_size = v; return *this; }
};

// (#1) Convenience run(): entry defaults to cfg.entry(), fixpoint knobs to
// run_config's defaults. Integer CFGs.
template <typename Dom>
void run(const crab::cfg_impl::z_cfg_t &cfg, Dom init, bool enable_stats,
         run_config cfg_opts = run_config()) {
  run(cfg, cfg.entry(), init, cfg_opts.run_liveness, cfg_opts.widening,
      cfg_opts.narrowing, cfg_opts.jump_set_size, enable_stats);
}

// Rational CFGs.
template <typename Dom>
void run(const crab::cfg_impl::q_cfg_t &cfg, Dom init, bool enable_stats,
         run_config cfg_opts = run_config()) {
  run(cfg, cfg.entry(), init, cfg_opts.run_liveness, cfg_opts.widening,
      cfg_opts.narrowing, cfg_opts.jump_set_size, enable_stats);
}

// Convenience run_and_check(): same idea for the checker entry point.
template <typename Dom>
void run_and_check(const crab::cfg_impl::z_cfg_t &cfg, Dom init, bool enable_stats,
                   run_config cfg_opts = run_config()) {
  run_and_check(cfg, cfg.entry(), init, cfg_opts.run_liveness,
                cfg_opts.widening, cfg_opts.narrowing, cfg_opts.jump_set_size,
                enable_stats);
}

// (#2) Sweep one CFG through a list of default-constructed domains, each with
// the same config, in the order given:
//
//   run_all<z_interval_domain_t, z_dbm_domain_t, z_sdbm_domain_t>(cfg, stats);
//
template <typename... Doms>
void run_all(const crab::cfg_impl::z_cfg_t &cfg, bool enable_stats,
             run_config cfg_opts = run_config()) {
  // C++14 has no fold expressions; expand through an initializer list, which
  // guarantees left-to-right evaluation so output order is deterministic.
  using expander = int[];
  (void)expander{0, (run(cfg, Doms(), enable_stats, cfg_opts), 0)...};
}

template <typename... Doms>
void run_all(const crab::cfg_impl::q_cfg_t &cfg, bool enable_stats,
             run_config cfg_opts = run_config()) {
  using expander = int[];
  (void)expander{0, (run(cfg, Doms(), enable_stats, cfg_opts), 0)...};
}

// (#2) Same sweep as run_all, but through the checker entry point: runs each
// default-constructed domain with run_and_check() in the order given:
//
//   check_all<z_interval_domain_t, z_sdbm_domain_t>(cfg, stats);
//
template <typename... Doms>
void check_all(const crab::cfg_impl::z_cfg_t &cfg, bool enable_stats,
               run_config cfg_opts = run_config()) {
  using expander = int[];
  (void)expander{0, (run_and_check(cfg, Doms(), enable_stats, cfg_opts), 0)...};
}

// (#3) unique_ptr twins: when a prog() builder returns std::unique_ptr<cfg>,
// these let call sites stay `run(cfg, ...)` instead of `run(cfg.get(), ...)`.
// unique_ptr does not implicitly convert to a raw pointer, so there is no
// ambiguity with the raw-pointer overloads above.
template <typename Dom>
void run(const std::unique_ptr<crab::cfg_impl::z_cfg_t> &cfg, Dom init,
         bool enable_stats, run_config cfg_opts = run_config()) {
  run(*cfg.get(), init, enable_stats, cfg_opts);
}
template <typename Dom>
void run(const std::unique_ptr<crab::cfg_impl::q_cfg_t> &cfg, Dom init,
         bool enable_stats, run_config cfg_opts = run_config()) {
  run(*cfg.get(), init, enable_stats, cfg_opts);
}
template <typename Dom>
void run_and_check(const std::unique_ptr<crab::cfg_impl::z_cfg_t> &cfg, Dom init,
                   bool enable_stats, run_config cfg_opts = run_config()) {
  run_and_check(*cfg.get(), init, enable_stats, cfg_opts);
}
template <typename... Doms>
void run_all(const std::unique_ptr<crab::cfg_impl::z_cfg_t> &cfg,
             bool enable_stats, run_config cfg_opts = run_config()) {
  run_all<Doms...>(*cfg.get(), enable_stats, cfg_opts);
}
template <typename... Doms>
void run_all(const std::unique_ptr<crab::cfg_impl::q_cfg_t> &cfg,
             bool enable_stats, run_config cfg_opts = run_config()) {
  run_all<Doms...>(*cfg.get(), enable_stats, cfg_opts);
}
template <typename... Doms>
void check_all(const std::unique_ptr<crab::cfg_impl::z_cfg_t> &cfg,
               bool enable_stats, run_config cfg_opts = run_config()) {
  check_all<Doms...>(*cfg.get(), enable_stats, cfg_opts);
}
