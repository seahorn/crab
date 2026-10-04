// Unit tests for the object domain and its building blocks (object_info,
// the copy-on-write abstract_domain_ref wrapper, odi_map_domain).
//
// Conventions: a state is never verified with the very operation under
// test -- expected facts are pinned through an independent observation
// (entails/at on registers after a load, inclusion against handcrafted
// states, ...). main() translates ctest's --disable-warnings and routes the
// Boost.Test log/report to stderr so stdout stays empty.
#define BOOST_TEST_MODULE object_domain_unittests
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_NO_MAIN
#include <boost/test/included/unit_test.hpp>

#include "object_dom.hpp"
#include <crab/domains/object_domain.hpp>

#include <string>
#include <vector>

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace crab::domains;
using namespace crab::object_domain_impl;

namespace {
using variable_or_constant = z_var_or_cst_t;
using variable_or_constant_vector_t = std::vector<z_var_or_cst_t>;

// crab types stream to crab_os, not std::ostream; stringify for Boost
// diagnostics.
template <typename T> std::string to_str(const T &x) {
  crab::crab_string_os os;
  os << x;
  return os.str();
}

// Scope the global object-domain parameters: save on entry, restore on exit,
// so suites with different configurations cannot leak into each other.
struct scoped_object_params {
  crab::domains::object_domain_params m_saved;
  scoped_object_params(
      crab::domains::object_domain_params::reduction_level_t level)
      : m_saved(crab_domain_params_man::get().reduction_level()) {
    crab::domains::object_domain_params p(level);
    crab_domain_params_man::get().update_params(p);
  }
  ~scoped_object_params() {
    crab_domain_params_man::get().update_params(m_saved);
  }
};

// Shared fixture: two abstract objects
//   objA = {V_a, V_b} (two integer fields), objB = {V_d} (one field),
// with region_init and the DSA intrinsic already applied, under the
// production configuration (OPT reduction, no singletons-in-base).
struct fixture {
  scoped_object_params params;
  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var A, A2, B, p, q, r, s, t, x, y, z;
  z_var rgn_a, rgn_b, rgn_d;
  z_var bcond;
  z_var_or_cst_t size8;

  fixture()
      : params(crab::domains::object_domain_params::reduction_level_t::
                   REDUCTION_BEFORE_CHECK),
        A(vfac["objA"], crab::REF_TYPE, 32),
        A2(vfac["objA2"], crab::REF_TYPE, 32),
        B(vfac["objB"], crab::REF_TYPE, 32), p(vfac["p"], crab::REF_TYPE, 32),
        q(vfac["q"], crab::REF_TYPE, 32), r(vfac["r"], crab::REF_TYPE, 32),
        s(vfac["s"], crab::REF_TYPE, 32), t(vfac["t"], crab::REF_TYPE, 32),
        x(vfac["x"], crab::INT_TYPE, 32), y(vfac["y"], crab::INT_TYPE, 32),
        z(vfac["z"], crab::INT_TYPE, 32),
        rgn_a(vfac["V_a"], crab::REG_INT_TYPE, 32),
        rgn_b(vfac["V_b"], crab::REG_INT_TYPE, 32),
        rgn_d(vfac["V_d"], crab::REG_INT_TYPE, 32),
        bcond(vfac["bc"], crab::BOOL_TYPE, 1),
        size8(z_number(8), crab::variable_type(crab::INT_TYPE, 32)) {}

  z_obj_zones_t init_state() {
    z_obj_zones_t inv;
    inv.region_init(rgn_a);
    inv.region_init(rgn_b);
    inv.region_init(rgn_d);
    variable_or_constant_vector_t objA_flds{variable_or_constant(rgn_a),
                                            variable_or_constant(rgn_b)};
    variable_or_constant_vector_t objB_flds{variable_or_constant(rgn_d)};
    inv.intrinsic("regions_from_memory_object", objA_flds, {});
    inv.intrinsic("regions_from_memory_object", objB_flds, {});
    return inv;
  }

  z_var_or_cst_t cst32(int n) {
    return z_var_or_cst_t(z_number(n),
                          crab::variable_type(crab::INT_TYPE, 32));
  }

  // Singleton objA: one allocation, fields written through p (V_a) and
  // r = p + 4 (V_b).
  z_obj_zones_t singleton_state(int va, int vb) {
    z_obj_zones_t inv = init_state();
    inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
    inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
    inv.ref_gep(p, rgn_a, r, rgn_b, z_number(4));
    inv.ref_store(p, rgn_a, cst32(va));
    inv.ref_store(r, rgn_b, cst32(vb));
    return inv;
  }

  // Non-singleton objA: a second allocation moves the first one's values
  // into the summary; the second allocation's values sit in the cache and
  // s/t own the MRU cache.
  z_obj_zones_t odi_state(int va1, int vb1, int va2, int vb2) {
    z_obj_zones_t inv = singleton_state(va1, vb1);
    inv.ref_make(A2, rgn_a, size8, as_man.mk_tag());
    inv.ref_gep(A2, rgn_a, s, rgn_a, z_number(0));
    inv.ref_gep(s, rgn_a, t, rgn_b, z_number(4));
    inv.ref_store(s, rgn_a, cst32(va2));
    inv.ref_store(t, rgn_b, cst32(vb2));
    return inv;
  }
};
} // namespace

// Check entails(cst) == expected, dumping the state on failure.
#define CHECK_ENTAILS(dom, cst, expected)                                      \
  BOOST_TEST(((dom).entails(cst) == (expected)),                               \
             std::string("entails(") + to_str(cst) + ") is not " +             \
                 ((expected) ? "true" : "false") + " in state " + to_str(dom))

// Check at(v) == expected, dumping the state on failure.
#define CHECK_AT(dom, v, expected)                                             \
  BOOST_TEST(((dom).at(v) == (expected)),                                      \
             to_str(v) + " = " + to_str((dom).at(v)) + ", expected " +         \
                 to_str(expected) + " in state " + to_str(dom))

BOOST_FIXTURE_TEST_SUITE(object_lattice, fixture)

BOOST_AUTO_TEST_CASE(lattice_constants) {
  z_obj_zones_t top;
  BOOST_TEST(top.is_top());
  BOOST_TEST(!top.is_bottom());
  z_obj_zones_t bot = top.make_bottom();
  BOOST_TEST(bot.is_bottom());
  z_obj_zones_t st = singleton_state(1, 2);
  BOOST_TEST((bot <= st));
  BOOST_TEST((st <= top));
  BOOST_TEST((!(top <= st)));
  BOOST_TEST((!(st <= bot)));
}

BOOST_AUTO_TEST_CASE(inclusion_on_register_facts) {
  z_obj_zones_t strong = init_state();
  strong += z_lin_cst_t(z_lin_exp_t(x) == z_number(0));
  z_obj_zones_t weak = init_state();
  weak += z_lin_cst_t(z_lin_exp_t(x) >= z_number(0));
  CHECK_ENTAILS(strong, z_lin_cst_t(z_lin_exp_t(x) <= z_number(0)), true);
  CHECK_ENTAILS(weak, z_lin_cst_t(z_lin_exp_t(x) <= z_number(0)), false);
  BOOST_TEST((strong <= weak));
  BOOST_TEST((!(weak <= strong)));
}

BOOST_AUTO_TEST_CASE(join_is_an_upper_bound) {
  z_obj_zones_t s1 = init_state();
  s1 += z_lin_cst_t(z_lin_exp_t(x) == z_number(0));
  z_obj_zones_t s2 = init_state();
  s2 += z_lin_cst_t(z_lin_exp_t(x) == z_number(1));
  z_obj_zones_t j = s1 | s2;
  BOOST_TEST((s1 <= j));
  BOOST_TEST((s2 <= j));
  CHECK_ENTAILS(j, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(j, z_lin_cst_t(z_lin_exp_t(x) <= z_number(1)), true);
  CHECK_ENTAILS(j, z_lin_cst_t(z_lin_exp_t(x) <= z_number(0)), false);
}

BOOST_AUTO_TEST_CASE(meet_is_a_lower_bound) {
  z_obj_zones_t s1 = init_state();
  s1 += z_lin_cst_t(z_lin_exp_t(x) <= z_number(5));
  z_obj_zones_t s2 = init_state();
  s2 += z_lin_cst_t(z_lin_exp_t(x) >= z_number(3));
  z_obj_zones_t m = s1 & s2;
  BOOST_TEST((m <= s1));
  BOOST_TEST((m <= s2));
  CHECK_AT(m, x, z_interval_t(z_number(3), z_number(5)));
  // contradictory register facts meet to bottom
  z_obj_zones_t s3 = init_state();
  s3 += z_lin_cst_t(z_lin_exp_t(x) == z_number(0));
  z_obj_zones_t s4 = init_state();
  s4 += z_lin_cst_t(z_lin_exp_t(x) == z_number(1));
  BOOST_TEST(((s3 & s4).is_bottom()));
}

// S5: both states hold objA in the odi map with committed summaries that
// contradict each other on V_a; their meet concretizes to the empty set and
// must be bottom.
BOOST_AUTO_TEST_CASE(meet_of_contradictory_fields_is_bottom) {
  z_obj_zones_t s1 = odi_state(0, 0, 0, 0); // all instances V_a = 0
  z_obj_zones_t s2 = odi_state(1, 1, 1, 1); // all instances V_a = 1
  z_obj_zones_t m = s1 & s2;
  BOOST_TEST(m.is_bottom());
}

// S15: widening with thresholds must clamp a growing field bound at the
// next threshold instead of losing it -- the fixpoint engine calls this
// entry point whenever max_thresholds > 0, so before the implementation
// object analyses silently degraded to plain widening. The witness uses an
// interval-based instantiation: crab's DBM/octagon domains currently
// implement widening_thresholds as a plain-widening fallback ("TODO: use
// thresholds"), so they cannot demonstrate the clamping; the object
// domain's job is only to route the thresholds faithfully to its base.
BOOST_AUTO_TEST_CASE(widening_with_thresholds_keeps_the_bound) {
  using z_obj_intervals_t =
      object_domain<crab::object_domain_impl::TestObjectParams<
          flat_boolean_numerical_domain<z_interval_domain_t>>>;
  variable_or_constant_vector_t objA_flds{variable_or_constant(rgn_a),
                                          variable_or_constant(rgn_b)};
  variable_or_constant_vector_t commit_args{variable_or_constant(rgn_a),
                                            variable_or_constant(rgn_b)};
  auto mk = [&](int va2, int vb2) {
    z_obj_intervals_t inv;
    inv.region_init(rgn_a);
    inv.region_init(rgn_b);
    inv.intrinsic("regions_from_memory_object", objA_flds, {});
    inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
    inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
    inv.ref_gep(p, rgn_a, r, rgn_b, z_number(4));
    inv.ref_store(p, rgn_a, cst32(0));
    inv.ref_store(r, rgn_b, cst32(0));
    inv.ref_make(A2, rgn_a, size8, as_man.mk_tag());
    inv.ref_gep(A2, rgn_a, s, rgn_a, z_number(0));
    inv.ref_gep(s, rgn_a, t, rgn_b, z_number(4));
    inv.ref_store(s, rgn_a, cst32(va2));
    inv.ref_store(t, rgn_b, cst32(vb2));
    inv.intrinsic("commit_cache", commit_args, {}); // summary: V_a in [0,va2]
    return inv;
  };
  z_obj_intervals_t s1 = mk(1, 1);
  z_obj_intervals_t s2 = mk(2, 2);
  z_obj_intervals_t j = s1 | s2;
  crab::thresholds<z_number> ts;
  ts.add(ikos::bound<z_number>(5));
  z_obj_intervals_t w_plain = s1 || j;
  z_obj_intervals_t w_ts = s1.widening_thresholds(j, ts);
  // plain widening loses the upper bound of V_a ...
  w_plain.ref_load(p, rgn_a, x);
  BOOST_TEST(!w_plain.entails(z_lin_cst_t(z_lin_exp_t(x) <= z_number(5))),
             "plain widening unexpectedly kept the bound");
  // ... widening with thresholds clamps it at the next threshold
  w_ts.ref_load(p, rgn_a, x);
  BOOST_TEST(w_ts.entails(z_lin_cst_t(z_lin_exp_t(x) <= z_number(5))),
             "thresholds were not routed to the base domain: " + to_str(w_ts));
}

// S11: a join with top must keep the shared side maps. The joined state
// must still know that rgn_a and rgn_b belong to the same object -- only
// the shared field->id map records that. Recreate the object THROUGH rgn_b
// only, then access V_a through a cross-field gep: if the grouping
// survived, the load goes through the object and forgets y; if the map was
// lost, rgn_a has no object and the load is skipped, leaving y's seeded
// value behind.
BOOST_AUTO_TEST_CASE(top_join_keeps_the_shared_maps) {
  z_obj_zones_t st = singleton_state(5, 7);
  z_obj_zones_t top;
  z_obj_zones_t j = st | top; // one operand top => join is top
  BOOST_TEST(j.is_top());
  j += z_lin_cst_t(z_lin_exp_t(y) == z_number(41)); // seed
  j.region_init(rgn_b);
  j.ref_make(A2, rgn_b, size8, as_man.mk_tag());
  j.ref_gep(A2, rgn_b, q, rgn_a, z_number(-4));
  j.ref_load(q, rgn_a, y);
  CHECK_AT(j, y, z_interval_t::top());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_memory_ops, fixture)

BOOST_AUTO_TEST_CASE(singleton_strong_update) {
  z_obj_zones_t inv = singleton_state(0, 1);
  inv.ref_load(p, rgn_a, x);
  inv.ref_load(r, rgn_b, y);
  CHECK_AT(inv, x, z_interval_t(z_number(0)));
  CHECK_AT(inv, y, z_interval_t(z_number(1)));
  // strong update: a second store overwrites, not accumulates
  inv.ref_store(p, rgn_a, cst32(7));
  inv.ref_load(p, rgn_a, x);
  CHECK_AT(inv, x, z_interval_t(z_number(7)));
}

BOOST_AUTO_TEST_CASE(cache_read_through_the_mru_object) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  // s/t refer to the MRU (second) allocation; loads hit the cache
  inv.ref_load(s, rgn_a, x);
  inv.ref_load(t, rgn_b, y);
  CHECK_AT(inv, x, z_interval_t(z_number(6)));
  CHECK_AT(inv, y, z_interval_t(z_number(12)));
  // p refers to the first allocation: a load through p misses the cache,
  // and only the summary (join of both writes) is available
  z_obj_zones_t inv2 = odi_state(0, 1, 6, 12);
  inv2.ref_load(p, rgn_a, x);
  CHECK_ENTAILS(inv2, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
  CHECK_ENTAILS(inv2, z_lin_cst_t(z_lin_exp_t(x) >= z_number(6)), false);
}

BOOST_AUTO_TEST_CASE(ref_assume_equality_gives_a_cache_hit) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  // q == s: q now refers to the MRU object, so a load through q hits the
  // cache and reads the exact stored value
  inv.ref_assume(z_ref_cst_t::mk_eq(q, s));
  inv.ref_load(q, rgn_a, x);
  CHECK_AT(inv, x, z_interval_t(z_number(6)));
}

// S-H witness: a same-region nonzero-offset gep (array stepping) must DROP
// the stepped reference's stale base-address equality, not merely refrain
// from adding one. With a reused variable (phi lowering does this), the
// stale equality lets the store through q strong-update the cache, wiping
// the other element's value: the load then claims exactly 7 where the
// sound answer is [5, 7].
BOOST_AUTO_TEST_CASE(array_gep_drops_stale_base_equality) {
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, q, rgn_a, z_number(0)); // q = &a[0]: q_base == A_base
  inv.ref_store(q, rgn_a, cst32(5));            // a[0] = 5; q's class owns MRU
  inv.ref_gep(A, rgn_a, q, rgn_a, z_number(4)); // q = &a[1]: array branch
  inv.ref_store(q, rgn_a, cst32(7));            // a[1] = 7: must NOT be a hit
  inv.ref_load(A, rgn_a, x);                    // x in [5,7], not exactly 7
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(7)), false);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) >= z_number(5)), true);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) <= z_number(7)), true);
}

// S-B witnesses: a load whose result is unknown must HAVOC the destination
// register, never leave its previous value in place.
BOOST_AUTO_TEST_CASE(load_from_unknown_region_havocs_res) {
  z_var rgn_u(vfac["V_unk"], crab::REG_UNKNOWN_TYPE, 32);
  z_obj_zones_t inv = init_state();
  inv.assign(x, 5);
  inv.ref_load(q, rgn_u, x);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(5)), false);
}

BOOST_AUTO_TEST_CASE(load_from_untracked_region_havocs_res) {
  // V_solo was never region_init'ed nor grouped: no object id exists
  z_var rgn_solo(vfac["V_solo"], crab::REG_INT_TYPE, 32);
  z_obj_zones_t inv = init_state();
  inv.assign(x, 5);
  inv.ref_load(q, rgn_solo, x);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(5)), false);
}

BOOST_AUTO_TEST_CASE(load_from_cast_region_havocs_res) {
  // region_cast installs an object whose refcount is still zero; the
  // zero-refcount exit must havoc like the others
  z_var rgn_u(vfac["V_unk2"], crab::REG_UNKNOWN_TYPE, 32);
  z_var rgn_c(vfac["V_cast"], crab::REG_INT_TYPE, 32);
  z_obj_zones_t inv = init_state();
  inv.region_cast(rgn_u, rgn_c);
  inv.assign(x, 5);
  inv.ref_load(q, rgn_c, x);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(5)), false);
}

// S-A witness: loading a REFERENCE from a field must re-point the loaded
// reference into the field's address class (keeping the field's binding),
// not the reverse. Discriminator: s owns objA's MRU cache (V_a = 7 exactly;
// the summary only knows the join with 1). A pointer field holding s is
// loaded into q; a read through q must be a cache HIT and see exactly 7.
BOOST_AUTO_TEST_CASE(ref_load_of_reference_repoints_the_loaded_ref) {
  z_var rgn_ptr(vfac["V_ptr"], crab::REG_REF_TYPE, 32);
  z_var h(vfac["h"], crab::REF_TYPE, 32);
  z_obj_zones_t inv = odi_state(1, 2, 7, 4);
  inv.region_init(rgn_ptr);
  variable_or_constant_vector_t objP_flds{variable_or_constant(rgn_ptr)};
  inv.intrinsic("regions_from_memory_object", objP_flds, {});
  inv.ref_make(h, rgn_ptr, size8, as_man.mk_tag());
  inv.ref_store(h, rgn_ptr, z_var_or_cst_t(s)); // field := s (the MRU owner)
  inv.ref_load(h, rgn_ptr, q);                  // q := *h
  inv.ref_load(q, rgn_a, x);                    // must HIT objA's cache
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(7)), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_var_ops, fixture)

BOOST_AUTO_TEST_CASE(forget_and_project_on_registers) {
  z_obj_zones_t inv = init_state();
  inv += z_lin_cst_t(z_lin_exp_t(x) == z_number(3));
  inv += z_lin_cst_t(z_lin_exp_t(y) == z_number(4));
  z_obj_zones_t forgotten = inv;
  forgotten.forget({x});
  CHECK_AT(forgotten, x, z_interval_t::top());
  CHECK_AT(forgotten, y, z_interval_t(z_number(4)));
  z_obj_zones_t projected = inv;
  projected.project({x});
  CHECK_AT(projected, x, z_interval_t(z_number(3)));
  CHECK_AT(projected, y, z_interval_t::top());
}

// Ported from the retired object_dm_ops.cc: forget/project on a state whose
// register-equality domain is live (x/y bound to fields through loads).
BOOST_AUTO_TEST_CASE(forget_and_project_with_live_equalities) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.ref_load(s, rgn_a, x); // x = 6, eq_regs x == #id(V_a)
  inv.ref_load(t, rgn_b, y); // y = 12, eq_regs y == #id(V_b)
  z_obj_zones_t fg = inv;
  fg.forget({t, x, rgn_b});
  CHECK_AT(fg, x, z_interval_t::top());
  CHECK_AT(fg, y, z_interval_t(z_number(12)));
  // V_a survives forgetting {t, x, V_b}: a load through s still hits the
  // cache
  fg.ref_load(s, rgn_a, y);
  CHECK_AT(fg, y, z_interval_t(z_number(6)));
  z_obj_zones_t pj = inv;
  pj.project({s, y, rgn_a});
  CHECK_AT(pj, y, z_interval_t(z_number(12)));
  CHECK_AT(pj, x, z_interval_t::top());
  pj.ref_load(s, rgn_a, x);
  CHECK_AT(pj, x, z_interval_t(z_number(6)));
}

// S2: singleton object: field values live in
// the odi map's cache, so renaming the region must rename inside the odi
// map, not (only) in the base domain.
BOOST_AUTO_TEST_CASE(rename_a_region_of_a_cached_object) {
  z_obj_zones_t inv = singleton_state(6, 1);
  z_var rgn_c(vfac["V_c"], crab::REG_INT_TYPE, 32);
  inv.rename({rgn_a}, {rgn_c});
  inv.ref_load(p, rgn_c, x);
  CHECK_AT(inv, x, z_interval_t(z_number(6)));
}

// Ported from the retired object_dm_ops.cc: rename on a state whose
// register-equality domain is live.
BOOST_AUTO_TEST_CASE(rename_with_live_equalities) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.ref_load(s, rgn_a, x); // x = 6
  z_var rgn_c(vfac["V_c"], crab::REG_INT_TYPE, 32);
  z_var p2(vfac["p2"], crab::REF_TYPE, 32);
  inv.rename({x, rgn_b, p}, {z, rgn_c, p2});
  CHECK_AT(inv, z, z_interval_t(z_number(6)));
  CHECK_AT(inv, x, z_interval_t::top());
  // the renamed region V_c keeps V_b's cached value
  inv.ref_load(t, rgn_c, y);
  CHECK_AT(inv, y, z_interval_t(z_number(12)));
}

// S14: forget(S) must behave exactly like folding operator-= over S -- the
// law region_domain follows. Before the unification, forget dropped a
// multi-field object's whole entry when ALL its fields were listed, while
// the fold partially forgets each field and keeps the entry (with its
// refcount, which is allocation info); the two disagreed on inclusion.
BOOST_AUTO_TEST_CASE(forget_is_the_fold_of_operator_minus) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.ref_load(s, rgn_a, x);
  inv.ref_load(t, rgn_b, y);
  z_obj_zones_t batched = inv;
  batched.forget({rgn_a, rgn_b, x, t});
  z_obj_zones_t folded = inv;
  folded -= rgn_a;
  folded -= rgn_b;
  folded -= x;
  folded -= t;
  BOOST_TEST((batched <= folded),
             "forget({...}) is not included in the operator-= fold");
  BOOST_TEST((folded <= batched),
             "the operator-= fold is not included in forget({...})");
  // single-field object: both forms drop the whole entry
  z_obj_zones_t inv2 = init_state();
  inv2.ref_make(B, rgn_d, size8, as_man.mk_tag());
  inv2.ref_gep(B, rgn_d, q, rgn_d, z_number(0));
  inv2.ref_store(q, rgn_d, cst32(7));
  z_obj_zones_t b1 = inv2;
  b1.forget({rgn_d});
  z_obj_zones_t b2 = inv2;
  b2 -= rgn_d;
  BOOST_TEST((b1 <= b2));
  BOOST_TEST((b2 <= b1));
}

// S3: projecting onto the live variables of an object must preserve the
// reference->MRU-cache relationship: a load through the MRU reference
// afterwards still reads the exact cached value.
BOOST_AUTO_TEST_CASE(project_keeps_the_mru_link) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.project({s, t, rgn_a, rgn_b});
  inv.ref_load(s, rgn_a, x);
  CHECK_AT(inv, x, z_interval_t(z_number(6)));
}

// T-6 baseline: storing a register links it to the field through the
// equality symbols; the value transfers back on the load (previously this
// path was only exercised through the ELINA-gated CFG tests).
BOOST_AUTO_TEST_CASE(store_register_then_load_reduces) {
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
  inv += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
  inv += z_lin_cst_t(z_lin_exp_t(x) <= z_number(9));
  inv.ref_store(p, rgn_a, z_var_or_cst_t(x));
  inv.ref_load(p, rgn_a, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) >= z_number(1)), true);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) <= z_number(9)), true);
}

// S-F witness: renaming a region must carry its derived ghosts along --
// the write-region ghost in the field-equality domain (exercised here)
// and, for reference regions, the base-address ghost. Same round trip as
// the baseline with the region renamed between store and load.
BOOST_AUTO_TEST_CASE(rename_carries_write_region_ghosts) {
  z_var rgn_a2(vfac["V_a2"], crab::REG_INT_TYPE, 32);
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
  inv += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
  inv += z_lin_cst_t(z_lin_exp_t(x) <= z_number(9));
  inv.ref_store(p, rgn_a, z_var_or_cst_t(x)); // pending: V_a_w == #s, x == #s
  inv.rename({rgn_a}, {rgn_a2});
  inv.ref_load(p, rgn_a2, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) >= z_number(1)), true);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) <= z_number(9)), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_conversions_and_observers, fixture)

// S4: weak_assign(x, 5) == join(state, state[x := 5]) => x in {0, 5}.
BOOST_AUTO_TEST_CASE(weak_assign_joins_with_the_old_value) {
  z_obj_zones_t inv = init_state();
  inv += z_lin_cst_t(z_lin_exp_t(x) == z_number(0));
  inv.weak_assign(x, z_lin_exp_t(z_number(5)));
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) <= z_number(5)), true);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) <= z_number(0)), false);
}

// S6: after loads bind registers to cached fields, entails() applies the
// pending reduction; the exported constraint system must agree with it.
BOOST_AUTO_TEST_CASE(to_linear_constraint_system_sees_reduced_facts) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.ref_load(s, rgn_a, x);
  z_lin_cst_t x_is_6(z_lin_exp_t(x) == z_number(6));
  CHECK_ENTAILS(inv, x_is_6, true);
  auto csts = inv.to_linear_constraint_system();
  // rebuild a plain numerical state from the exported constraints and ask
  // it the same question
  z_soct_domain_t exported;
  for (auto const &c : csts) {
    exported += c;
  }
  BOOST_TEST(exported.entails(x_is_6),
             "exported constraints lost x = 6: " + to_str(exported));
}

// S1: the selected result with an unknown condition must over-approximate
// both branches.
BOOST_AUTO_TEST_CASE(select_ref_over_approximates_both_branches) {
  z_obj_zones_t base = odi_state(0, 1, 6, 12);
  // Manually apply each branch of the select
  z_obj_zones_t b1 = base;
  b1.ref_gep(s, rgn_a, q, rgn_a, z_number(4)); // same-region step
  z_obj_zones_t b2 = base;
  b2.ref_gep(p, rgn_a, q, rgn_a, z_number(0));
  z_obj_zones_t sel = base;
  sel.select_ref(q, rgn_a, bcond, variable_or_constant(s),
                 boost::optional<z_var>(rgn_a), variable_or_constant(p),
                 boost::optional<z_var>(rgn_a));
  BOOST_TEST((b1 <= sel));
  BOOST_TEST((b2 <= sel));
}

// S9: overwriting a register through ref_to_int must kill its old field
// link, or a later reduction re-imposes the field's value on an address.
BOOST_AUTO_TEST_CASE(ref_to_int_clears_the_stale_register_symbol) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  // bind x to field V_a through a load (x = 6, eq_regs x == #id(V_a))
  inv.ref_load(s, rgn_a, x);
  inv.ref_to_int(rgn_a, s, x);
  // a second load through s re-uses V_a's class id (get_symbol_or_fresh)
  // and re-arms the reduction: with a stale symbol still on x, the
  // reduction would re-impose V_a's value on the address-valued x
  inv.ref_load(s, rgn_a, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) == z_number(6)), false);
}

// S9: s owns the MRU cache (V_a = 6 exactly); re-targeting s to an
// arbitrary address must drop that ownership, so a load through s no
// longer reads the exact cached value.
BOOST_AUTO_TEST_CASE(int_to_ref_drops_stale_cache_ownership) {
  z_obj_zones_t inv = odi_state(0, 1, 6, 12);
  inv.int_to_ref(z, rgn_a, s);
  inv.ref_load(s, rgn_a, x);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) >= z_number(6)), false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_region_copy, fixture)

BOOST_AUTO_TEST_CASE(region_copy_transfers_the_field_value) {
  z_obj_zones_t inv = init_state();
  // objB = {V_d} is a single-field object; write through it
  inv.ref_make(B, rgn_d, size8, as_man.mk_tag());
  inv.ref_gep(B, rgn_d, q, rgn_d, z_number(0));
  inv.ref_store(q, rgn_d, cst32(7));
  // copy into a fresh single-field region
  z_var rgn_e(vfac["V_e"], crab::REG_INT_TYPE, 32);
  inv.region_init(rgn_e);
  variable_or_constant_vector_t objE_flds{variable_or_constant(rgn_e)};
  inv.intrinsic("regions_from_memory_object", objE_flds, {});
  inv.region_copy(rgn_e, rgn_d);
  inv.ref_load(q, rgn_e, x);
  CHECK_AT(inv, x, z_interval_t(z_number(7)));
}

// S7: two runs that differ ONLY in the value stored to V_b -- a field of
// objA that V_e's object does not contain. After copying V_a into V_e and
// forgetting rgn_b, the states are semantically identical, so mutual
// inclusion must hold. Before the fix, rhs's sibling dims (V_b) stayed
// inside V_e's entry -- unreachable and unforgettable through lhs -- and
// broke the inclusion.
BOOST_AUTO_TEST_CASE(no_sibling_dims_leak_into_the_copied_object) {
  z_var rgn_e(vfac["V_e"], crab::REG_INT_TYPE, 32);
  variable_or_constant_vector_t objE_flds{variable_or_constant(rgn_e)};
  auto mk = [&](int vb) {
    z_obj_zones_t inv = odi_state(0, 1, 6, vb);
    inv.region_init(rgn_e);
    inv.intrinsic("regions_from_memory_object", objE_flds, {});
    inv.region_copy(rgn_e, rgn_a);
    inv.forget({rgn_b}); // drop the only real difference
    return inv;
  };
  z_obj_zones_t A_st = mk(12);
  z_obj_zones_t B_st = mk(99);
  BOOST_TEST((A_st <= B_st),
             "A <= B fails on semantically identical states");
  BOOST_TEST((B_st <= A_st),
             "B <= A fails on semantically identical states");
  // the copied value itself still transfers
  z_obj_zones_t C_st = A_st;
  C_st.ref_load(q, rgn_e, x); // q was never bound: miss reads the summary
  CHECK_ENTAILS(C_st, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_intrinsics, fixture)

// S13: an object's DSA field list may mix known and unknown regions; clam
// hands the same full list to every intrinsic. commit_cache must skip the
// unknown region and still commit the objects named by the known ones (it
// used to abandon the whole intrinsic on the first unknown).
BOOST_AUTO_TEST_CASE(commit_cache_tolerates_unknown_regions) {
  z_var rgn_u(vfac["V_u"], crab::REG_UNKNOWN_TYPE, 32);
  z_obj_zones_t inv = odi_state(0, 1, 6, 12); // dirty cache: V_a = 6
  variable_or_constant_vector_t commit_args{variable_or_constant(rgn_a),
                                            variable_or_constant(rgn_u)};
  inv.intrinsic("commit_cache", commit_args, {});
  // committed: the cache was folded into the summary and reset, so a load
  // through the (former) MRU reference reads the summary join [0,6], no
  // longer the exact cached value
  inv.ref_load(s, rgn_a, x);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) >= z_number(6)), false);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
}

// S12/S13: one unknown region in the DSA field list must not discard the
// whole copy -- object_copy filters unknown src/dst pairs itself, and only
// the known regions are registered as fields of the destination object.
BOOST_AUTO_TEST_CASE(copy_memory_object_tolerates_unknown_regions) {
  z_var rgn_u(vfac["V_u"], crab::REG_UNKNOWN_TYPE, 32);
  z_var rgn_u2(vfac["V_u2"], crab::REG_UNKNOWN_TYPE, 32);
  z_obj_zones_t inv = init_state();
  // objB = {V_d}: write 7 through it
  inv.ref_make(B, rgn_d, size8, as_man.mk_tag());
  inv.ref_gep(B, rgn_d, q, rgn_d, z_number(0));
  inv.ref_store(q, rgn_d, cst32(7));
  // destination object {V_e}; the copy list carries an unknown pair too
  z_var rgn_e(vfac["V_e"], crab::REG_INT_TYPE, 32);
  inv.region_init(rgn_e);
  variable_or_constant_vector_t objE_flds{variable_or_constant(rgn_e)};
  inv.intrinsic("regions_from_memory_object", objE_flds, {});
  variable_or_constant_vector_t copy_in{variable_or_constant(rgn_d),
                                        variable_or_constant(rgn_u)};
  std::vector<z_var> copy_out{rgn_e, rgn_u2};
  inv.intrinsic("copy_memory_object", copy_in, copy_out);
  // the known pair was copied despite the unknown one
  inv.ref_load(q, rgn_e, x);
  CHECK_AT(inv, x, z_interval_t(z_number(7)));
}

BOOST_AUTO_TEST_SUITE_END()

// Ported from the retired object_join.cc: the four documented join
// scenarios, with the printed expected states turned into assertions.
BOOST_FIXTURE_TEST_SUITE(object_join_scenarios, fixture)

// join 1: both states hold the same singleton object with different cached
// values; the join commits both caches, and the summary keeps the
// relational fact V_a < V_b.
BOOST_AUTO_TEST_CASE(join_of_two_singleton_states) {
  z_obj_zones_t inv1 = singleton_state(0, 1);
  z_obj_zones_t inv2 = inv1;
  inv2.ref_store(p, rgn_a, cst32(5));
  inv2.ref_store(r, rgn_b, cst32(10));
  z_obj_zones_t res = inv1 | inv2;
  BOOST_TEST((inv1 <= res));
  BOOST_TEST((inv2 <= res));
  res.ref_load(p, rgn_a, x);
  res.ref_load(r, rgn_b, y);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) <= z_number(5)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(y) >= z_number(1)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(y) <= z_number(10)), true);
  // the intra-object relation V_a < V_b survives the join
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) <= z_lin_exp_t(y) - 1), true);
}

// join 2: one state singleton, the other already in the odi map.
BOOST_AUTO_TEST_CASE(join_of_singleton_with_odi_state) {
  z_obj_zones_t inv1 = singleton_state(0, 1);
  z_obj_zones_t inv3 = inv1;
  inv3.ref_make(A2, rgn_a, size8, as_man.mk_tag());
  inv3.ref_gep(A2, rgn_a, s, rgn_a, z_number(0));
  inv3.ref_store(s, rgn_a, cst32(6));
  inv3.ref_gep(s, rgn_a, t, rgn_b, z_number(4));
  inv3.ref_store(t, rgn_b, cst32(12));
  z_obj_zones_t res = inv1 | inv3;
  BOOST_TEST((inv1 <= res));
  BOOST_TEST((inv3 <= res));
  res.ref_load(p, rgn_a, x);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
}

// join 3: two odi states sharing the SAME MRU object. The implementation
// always flushes caches when joining -- even when both sides share the MRU
// object -- so the exact cached values are not preserved. This is a known,
// documented imprecision (see object_join.cc history); if the join ever
// learns to keep a shared cache, the last check below should flip to exact.
BOOST_AUTO_TEST_CASE(join_of_odi_states_with_same_mru_flushes_the_cache) {
  z_obj_zones_t inv3 = odi_state(0, 1, 6, 12);
  z_obj_zones_t inv4 = inv3;
  inv4.ref_store(s, rgn_a, cst32(4));
  inv4.ref_store(t, rgn_b, cst32(8));
  z_obj_zones_t res = inv3 | inv4;
  BOOST_TEST((inv3 <= res));
  BOOST_TEST((inv4 <= res));
  res.ref_load(s, rgn_a, x);
  // summary join bounds hold ...
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
  // ... but the cache was flushed, so the load is NOT pinned to {4,6}
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) >= z_number(4)), false);
}

// join 4: two odi states with DIFFERENT MRU objects; caches flushed, the
// joined summary bounds all instances.
BOOST_AUTO_TEST_CASE(join_of_odi_states_with_different_mru) {
  z_obj_zones_t inv3 = odi_state(0, 1, 6, 12);
  z_obj_zones_t inv5 = inv3;
  inv5.ref_store(p, rgn_a, cst32(4)); // p's instance becomes the MRU
  inv5.ref_store(r, rgn_b, cst32(8));
  z_obj_zones_t res = inv3 | inv5;
  BOOST_TEST((inv3 <= res));
  BOOST_TEST((inv5 <= res));
  res.ref_load(p, rgn_a, x);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) >= z_number(0)), true);
  CHECK_ENTAILS(res, z_lin_cst_t(z_lin_exp_t(x) <= z_number(6)), true);
}

// S-E witness: joining with a top state must merge BOTH sides' shared DSA
// maps -- ideally all states converge on one shared map. The dropped odi
// ENTRY correctly stays dropped (the object is top), so the map's survival
// is observed through copy_memory_object: its source side requires the
// binding (get_obj_id_or_fail aborts the analysis without it).
BOOST_AUTO_TEST_CASE(top_join_merges_both_sides_shared_maps) {
  variable_or_constant_vector_t fa{variable_or_constant(rgn_a)};
  variable_or_constant_vector_t fb{variable_or_constant(rgn_d)};
  auto mk_a = [&]() {
    z_obj_zones_t st;
    st.region_init(rgn_a);
    st.intrinsic("regions_from_memory_object", fa, {});
    st.ref_make(A, rgn_a, size8, as_man.mk_tag());
    st.ref_store(A, rgn_a, cst32(1)); // non-top
    return st;
  };
  z_obj_zones_t b;
  b.region_init(rgn_d);
  b.intrinsic("regions_from_memory_object", fb, {});
  b.set_to_top(); // keeps b's shared maps

  z_var rgn_e2(vfac["V_e2"], crab::REG_INT_TYPE, 32);
  variable_or_constant_vector_t copy_in{variable_or_constant(rgn_d)};
  std::vector<z_var> copy_out{rgn_e2};

  z_obj_zones_t j = mk_a() | b; // operator|
  BOOST_TEST(j.is_top());
  j.intrinsic("copy_memory_object", copy_in, copy_out);
  BOOST_TEST(j.is_top());

  z_obj_zones_t a2 = mk_a();
  a2 |= b; // operator|=
  BOOST_TEST(a2.is_top());
  a2.intrinsic("copy_memory_object", copy_in, copy_out);
  BOOST_TEST(a2.is_top());
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_inter_ops, fixture)

// T-2: callee_entry must transfer the caller's object invariants onto the
// callee's formal regions (grouped DSA parameters -> object_copy path).
BOOST_AUTO_TEST_CASE(callee_entry_propagates_object_invariants) {
  z_var c_rgn(vfac["V_c"], crab::REG_INT_TYPE, 32);
  z_var f_rgn(vfac["V_f"], crab::REG_INT_TYPE, 32);
  z_var pf(vfac["pf"], crab::REF_TYPE, 32);
  z_obj_zones_t caller;
  caller.region_init(c_rgn);
  variable_or_constant_vector_t fc{variable_or_constant(c_rgn)};
  caller.intrinsic("regions_from_memory_object", fc, {});
  caller.ref_make(A, c_rgn, size8, as_man.mk_tag());
  caller.ref_store(A, c_rgn, cst32(4));

  std::vector<z_var> caller_in{c_rgn}, caller_out, callee_in{f_rgn}, callee_out;
  crab::domains::callsite_info<z_var> cs("f", caller_in, caller_out, callee_in,
                                         callee_out, {{c_rgn}}, {}, {{f_rgn}},
                                         {});
  z_obj_zones_t callee;
  callee.callee_entry(cs, caller);
  callee.ref_make(pf, f_rgn, size8, as_man.mk_tag());
  callee.ref_load(pf, f_rgn, y);
  CHECK_ENTAILS(callee, z_lin_cst_t(z_lin_exp_t(y) == z_number(4)), true);
}

// T-3: caller_continuation must propagate the callee's exit-state object
// invariants back onto the caller's actual regions.
BOOST_AUTO_TEST_CASE(caller_continuation_returns_object_facts) {
  z_var c_rgn(vfac["V_c"], crab::REG_INT_TYPE, 32);
  z_var f_rgn(vfac["V_f"], crab::REG_INT_TYPE, 32);
  z_var pc(vfac["pc"], crab::REF_TYPE, 32);
  // callee at exit: object over f_rgn holding exactly 7
  z_obj_zones_t callee;
  callee.region_init(f_rgn);
  variable_or_constant_vector_t ff{variable_or_constant(f_rgn)};
  callee.intrinsic("regions_from_memory_object", ff, {});
  callee.ref_make(B, f_rgn, size8, as_man.mk_tag());
  callee.ref_store(B, f_rgn, cst32(7));
  // caller at the callsite: knows its own object over c_rgn
  z_obj_zones_t caller;
  caller.region_init(c_rgn);
  variable_or_constant_vector_t fc{variable_or_constant(c_rgn)};
  caller.intrinsic("regions_from_memory_object", fc, {});

  std::vector<z_var> caller_in{c_rgn}, caller_out, callee_in{f_rgn}, callee_out;
  crab::domains::callsite_info<z_var> cs("f", caller_in, caller_out, callee_in,
                                         callee_out, {{c_rgn}}, {}, {{f_rgn}},
                                         {});
  caller.caller_continuation(cs, callee);
  caller.ref_make(pc, c_rgn, size8, as_man.mk_tag());
  caller.ref_load(pc, c_rgn, y);
  CHECK_ENTAILS(caller, z_lin_cst_t(z_lin_exp_t(y) == z_number(7)), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(object_reduction_configs, fixture)

// T-4: under NO_REDUCTION the symbolic store must NOT smuggle the
// register's bounds into the field (that transfer needs the reduction);
// soundness is preserved by dedicated forgets in ref_store.
BOOST_AUTO_TEST_CASE(no_reduction_store_is_sound) {
  scoped_object_params no_red(
      crab::domains::object_domain_params::reduction_level_t::NO_REDUCTION);
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
  inv += z_lin_cst_t(z_lin_exp_t(x) >= z_number(1));
  inv += z_lin_cst_t(z_lin_exp_t(x) <= z_number(9));
  inv.ref_store(p, rgn_a, z_var_or_cst_t(x));
  inv.ref_load(p, rgn_a, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) >= z_number(1)), false);
}

// T-5: FULL_REDUCTION must agree with the default configuration on the
// core scenarios (reduce-at-every-transfer-function composition).
BOOST_AUTO_TEST_CASE(full_reduction_matches_default_reduction) {
  scoped_object_params full(
      crab::domains::object_domain_params::reduction_level_t::FULL_REDUCTION);
  z_obj_zones_t s1 = singleton_state(1, 2);
  s1.ref_load(p, rgn_a, x);
  CHECK_ENTAILS(s1, z_lin_cst_t(z_lin_exp_t(x) == z_number(1)), true);
  z_obj_zones_t s2 = odi_state(1, 2, 7, 4);
  s2.ref_load(s, rgn_a, y);
  CHECK_ENTAILS(s2, z_lin_cst_t(z_lin_exp_t(y) == z_number(7)), true);
  s2.ref_load(t, rgn_b, z);
  CHECK_ENTAILS(s2, z_lin_cst_t(z_lin_exp_t(z) == z_number(4)), true);
}

// W-1: the second store kills the field's old cache value. x is a
// register: `x := 7` AFTER the store does not change memory, so the sound
// result is "y = x's value at store time, somewhere in [0,10]" -- never
// exactly the PRE-store 5. Pre-fix, the default configuration answered
// y = [5,5]: the pending store's symbol died with the register
// reassignment and the commit that would have dropped the stale 5 was
// skipped with it (regs.empty() -> continue).
BOOST_AUTO_TEST_CASE(pending_store_never_resurrects_old_value) {
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
  inv.ref_store(p, rgn_a, cst32(5));
  inv += z_lin_cst_t(z_lin_exp_t(x) >= z_number(0));
  inv += z_lin_cst_t(z_lin_exp_t(x) <= z_number(10));
  inv.ref_store(p, rgn_a, z_var_or_cst_t(x)); // overwrites the 5
  inv.assign(x, z_lin_exp_t(z_number(7)));    // register only; memory keeps
                                              // x's value at store time
  inv.ref_load(p, rgn_a, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) == z_number(5)), false);
}

// Same program under FULL_REDUCTION: the reduction at the assignment
// commits the pending store while the register is still linked, so this
// configuration was sound even before the fix.
BOOST_AUTO_TEST_CASE(pending_store_never_resurrects_old_value_full) {
  scoped_object_params full(
      crab::domains::object_domain_params::reduction_level_t::FULL_REDUCTION);
  z_obj_zones_t inv = init_state();
  inv.ref_make(A, rgn_a, size8, as_man.mk_tag());
  inv.ref_gep(A, rgn_a, p, rgn_a, z_number(0));
  inv.ref_store(p, rgn_a, cst32(5));
  inv += z_lin_cst_t(z_lin_exp_t(x) >= z_number(0));
  inv += z_lin_cst_t(z_lin_exp_t(x) <= z_number(10));
  inv.ref_store(p, rgn_a, z_var_or_cst_t(x));
  inv.assign(x, z_lin_exp_t(z_number(7)));
  inv.ref_load(p, rgn_a, y);
  CHECK_ENTAILS(inv, z_lin_cst_t(z_lin_exp_t(y) == z_number(5)), false);
}

BOOST_AUTO_TEST_SUITE_END()

// ---------------------------------------------------------------------------
// Subdomain tests: object_info, the COW wrapper and the odi map itself.
// ---------------------------------------------------------------------------

namespace {
using base_t =
    flat_boolean_numerical_domain<crab::domain_impl::z_soct_domain_t>;
// Force object_domain's instantiation first: odi_map_domain's nested types
// depend on it, and naming them below would otherwise trigger a recursive
// instantiation while object_domain is still incomplete.
static_assert(sizeof(z_obj_zones_t) > 0, "complete the object domain");
using odi_map_t =
    object_domain_impl::odi_map_domain<z_var, z_obj_zones_t, base_t>;
using odi_info_t = odi_map_t::odi_info_t;
using odi_value_t = odi_map_t::odi_value_t;
using map_raw_value_t = odi_map_t::map_raw_value_t;
using cow_ref_t = object_domain_impl::abstract_domain_ref<z_var, base_t>;

odi_info_t mk_info(const small_range &count, bool used, bool dirty) {
  return odi_info_t(count, /*obj_init=*/boolean_value::get_true(),
                    /*sum_presence=*/boolean_value::get_true(),
                    /*cache_used=*/used ? boolean_value::get_true()
                                        : boolean_value::get_false(),
                    /*cache_dirty=*/dirty ? boolean_value::get_true()
                                          : boolean_value::get_false(),
                    /*is_loaded=*/false, /*is_stored=*/false);
}

// An odi map entry whose SUMMARY holds the single constraint `v == n`,
// with a clean, unused cache.
map_raw_value_t mk_summary_entry(const z_var &v, int n) {
  odi_value_t val;
  (*val.first()) += z_lin_cst_t(z_lin_exp_t(v) == z_number(n));
  return map_raw_value_t(mk_info(small_range::oneOrMore(), false, false),
                         std::move(val));
}

struct subdom_fixture {
  scoped_object_params params;
  variable_factory_t vfac;
  z_var id, id2, va, vb;
  subdom_fixture()
      : params(crab::domains::object_domain_params::reduction_level_t::
                   REDUCTION_BEFORE_CHECK),
        id(vfac["objA"], crab::REG_INT_TYPE, 32),
        id2(vfac["objB"], crab::REG_INT_TYPE, 32),
        va(vfac["V_a"], crab::INT_TYPE, 32),
        vb(vfac["V_b"], crab::INT_TYPE, 32) {}
};
} // namespace

BOOST_FIXTURE_TEST_SUITE(object_subdomains, subdom_fixture)

BOOST_AUTO_TEST_CASE(object_info_accessors) {
  // ExactlyOne tracks the allocation it counts, so it is built the way
  // the domain builds it: zero() then increment(ref).
  small_range one = small_range::zero();
  one.increment(id);
  odi_info_t info = mk_info(one, true, false);
  BOOST_TEST(info.refcount_val().is_one());
  BOOST_TEST(info.cacheused_val().is_true());
  BOOST_TEST(info.cachedirty_val().is_false());
  BOOST_TEST(!info.cache_reg_loaded_val());
  info.cache_reg_loaded_val() = true;
  BOOST_TEST(info.cache_reg_loaded_val());
  odi_info_t copy = info;
  copy.cache_reg_loaded_val() = false;
  BOOST_TEST(info.cache_reg_loaded_val(),
             "object_info copies must not alias");
}

BOOST_AUTO_TEST_CASE(copy_on_write_wrapper) {
  cow_ref_t a;
  cow_ref_t b = a; // shared
  (*b) += z_lin_cst_t(z_lin_exp_t(va) == z_number(1));
  const cow_ref_t &a_view = a;
  BOOST_TEST((*a_view).is_top(),
             "mutation of a copy leaked into the original: " +
                 to_str(*a_view));
  BOOST_TEST((*b).entails(z_lin_cst_t(z_lin_exp_t(va) == z_number(1))));
  cow_ref_t c = b;
  BOOST_TEST((*static_cast<const cow_ref_t &>(c))
                 .entails(z_lin_cst_t(z_lin_exp_t(va) == z_number(1))));
}

BOOST_AUTO_TEST_CASE(odi_map_join_and_inclusion) {
  odi_map_t m1;
  m1.set(id, mk_summary_entry(va, 0));
  odi_map_t m2;
  m2.set(id, mk_summary_entry(va, 1));
  BOOST_TEST((!(m1 <= m2)));
  odi_map_t j = m1.join(m2);
  BOOST_TEST((m1 <= j));
  BOOST_TEST((m2 <= j));
  // join with an empty (top) map is top
  odi_map_t top_map;
  BOOST_TEST(m1.join(top_map).is_top());
  // an entry present on one side only is dropped by the join
  // (default_is_absorbing = true): the missing side means top
  odi_map_t m3;
  m3.set(id2, mk_summary_entry(vb, 5));
  BOOST_TEST(m1.join(m3).is_top());
}

BOOST_AUTO_TEST_CASE(odi_map_meet_keeps_both_objects) {
  odi_map_t m1;
  m1.set(id, mk_summary_entry(va, 0));
  odi_map_t m3;
  m3.set(id2, mk_summary_entry(vb, 5));
  odi_map_t m = m1.meet(m3);
  BOOST_TEST((m <= m1));
  BOOST_TEST((m <= m3));
}

// S5: a contradictory per-object meet must make the whole map bottom.
BOOST_AUTO_TEST_CASE(odi_map_meet_signals_bottom) {
  odi_map_t m1;
  m1.set(id, mk_summary_entry(va, 0));
  odi_map_t m2;
  m2.set(id, mk_summary_entry(va, 1));
  odi_map_t m = m1.meet(m2);
  BOOST_TEST(m.is_bottom());
}

// S17: the tree's equality functor decides equality in O(1) by identity:
// a copy shares the three inner COW refs, so it compares equal; detaching
// any component or changing an info field breaks the equality.
BOOST_AUTO_TEST_CASE(odi_value_equality_is_identity_based) {
  object_domain_impl::object_equal_to<odi_map_t::mapped_type> eq;
  auto p1 = std::make_shared<map_raw_value_t>(mk_summary_entry(va, 0));
  BOOST_TEST(eq(p1, p1)); // same allocation
  // fresh outer allocation, shared inner COW refs -> still equal
  auto p2 = std::make_shared<map_raw_value_t>(*p1);
  BOOST_TEST(eq(p1, p2));
  // detaching one component (a cache write) breaks the identity
  auto p3 = std::make_shared<map_raw_value_t>(*p1);
  (*p3->second().second().first()) +=
      z_lin_cst_t(z_lin_exp_t(vb) == z_number(1));
  BOOST_TEST(!eq(p1, p3));
  // same subdomain refs, but an info field differs -> not equal
  auto p4 = std::make_shared<map_raw_value_t>(*p1);
  p4->first().cachedirty_val() = boolean_value::get_true();
  BOOST_TEST(!eq(p1, p4));
}

// S17 end to end: set() with a value equal to the current one must keep the
// OLD tree leaf (find() returns the same address), preserving structural
// sharing across states instead of allocating a fresh node.
BOOST_AUTO_TEST_CASE(odi_map_noop_set_keeps_leaf) {
  odi_map_t m1;
  m1.set(id, mk_summary_entry(va, 0));
  odi_map_t m2 = m1; // shares the tree, and thus the leaf
  const map_raw_value_t *before = m2.find(id);
  map_raw_value_t copy = *before; // same info, shared inner COW refs
  m2.set(id, copy);
  BOOST_TEST((m2.find(id) == before), "no-op set() dropped leaf sharing");
  BOOST_TEST((m1.find(id) == before));
}

// Dirty cache {va = 3} over summary {va = 0}: committing must widen the
// summary to cover the cached value and reset the cache to top.
BOOST_AUTO_TEST_CASE(commit_cache_if_dirty_folds_and_resets) {
  odi_value_t val;
  (*val.first()) += z_lin_cst_t(z_lin_exp_t(va) == z_number(0));
  (*odi_map_t::object_cache_val(val)) +=
      z_lin_cst_t(z_lin_exp_t(va) == z_number(3));
  odi_info_t info = mk_info(small_range::oneOrMore(), true, true);
  odi_map_t m;
  m.commit_cache_if_dirty(info, val);
  const base_t &sum = *val.first();
  BOOST_TEST(sum.entails(z_lin_cst_t(z_lin_exp_t(va) >= z_number(0))),
             "summary lost the old value: " + to_str(sum));
  BOOST_TEST(sum.entails(z_lin_cst_t(z_lin_exp_t(va) <= z_number(3))),
             "summary lost the cached value: " + to_str(sum));
  BOOST_TEST((*odi_map_t::object_cache_val(val)).is_top(),
             "cache was not reset: " +
                 to_str(*odi_map_t::object_cache_val(val)));
}

BOOST_AUTO_TEST_SUITE_END()

int main(int argc, char **argv) {
  std::vector<char *> args;
  for (int idx = 0; idx < argc; ++idx) {
    if (std::string(argv[idx]) == "--disable-warnings") {
      crab::CrabEnableWarningMsg(false);
      continue;
    }
    args.push_back(argv[idx]);
  }
  char log_sink[] = "--log_sink=stderr";
  char report_sink[] = "--report_sink=stderr";
  args.push_back(log_sink);
  args.push_back(report_sink);
  return boost::unit_test::unit_test_main(
      &init_unit_test, static_cast<int>(args.size()), args.data());
}
