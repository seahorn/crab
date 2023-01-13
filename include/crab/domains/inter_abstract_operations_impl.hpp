#pragma once
/**==============================================================**/
/**  To be included only by inter_abstract_operations.hpp        **/
/**==============================================================**/
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/variable.hpp>

#include <algorithm> // sort, set_difference
#include <string>
#include <vector>

namespace crab {
namespace domains {

template <class Domain>
class inter_abstract_operations<Domain, false /*UseDefaultImpl*/> {
public:
  using variable_t = typename Domain::variable_t;
  static void callee_entry(const callsite_info<variable_t> &callsite,
                           const Domain &caller, Domain &callee) {
    CRAB_ERROR(caller.domain_name(), " does not implement callee_entry");
  }

  static void caller_continuation(const callsite_info<variable_t> &callsite,
                                  const Domain &callee, Domain &caller) {
    CRAB_ERROR(callee.domain_name(), " does not implement caller_continuation");
  }
};

namespace inter_transformers_impl {

template <class Domain>
inline void unify(Domain &inv, const typename Domain::variable_t &lhs,
                  const typename Domain::variable_t &rhs) {
  using reference_constraint_t = typename Domain::reference_constraint_t;
  using number_t = typename Domain::number_t;

  assert(lhs.get_type() == rhs.get_type());
  auto ty = lhs.get_type();
  if (ty.is_bool()) {
    inv.assign_bool_var(lhs, rhs, false);
  } else if (ty.is_integer() || ty.is_real()) {
    inv.assign(lhs, rhs);
  } else if (ty.is_reference()) {
    inv -= lhs;
    inv.ref_assume(reference_constraint_t::mk_eq(lhs, rhs, number_t(0)));
  } else if (ty.is_region()) {
    inv.region_copy(lhs, rhs);
  } else if (ty.is_array()) {
    inv.array_assign(lhs, rhs);
  } else {
    CRAB_ERROR("abs_transformer::unify unsupported type");
  }
}

template <class V>
inline std::vector<V> set_difference(std::vector<V> v1, std::vector<V> v2) {
  std::vector<V> out;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),
                      std::back_inserter(out));
  return out;
}

template <class V>
inline std::vector<V> set_intersection(std::vector<V> v1, std::vector<V> v2) {
  std::vector<V> out;
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                        std::back_inserter(out));
  return out;
}

} // end namespace inter_transformers_impl

template <class Domain>
class inter_abstract_operations<Domain, true /*UseDefaultImpl*/> {
public:
  using variable_t = typename Domain::variable_t;
  static void callee_entry(const callsite_info<variable_t> &callsite,
                           const Domain &caller, Domain &callee);
  static void caller_continuation(const callsite_info<variable_t> &callsite,
                                  const Domain &callee, Domain &caller);
};

template <class Domain>
void inter_abstract_operations<Domain, true>::callee_entry(
    const callsite_info<typename Domain::variable_t> &callsite,
    const Domain &_caller, Domain &callee_at_entry) {

  using variable_t = typename Domain::variable_t;

  crab::CrabStats::count(_caller.domain_name() + ".callee_entry");
  crab::ScopedCrabStats __st__(_caller.domain_name() + ".callee_entry");

  if (_caller.is_bottom()) {
    callee_at_entry.set_to_bottom();
    return;
  }

  Domain caller(_caller);

  CRAB_LOG("inter-restrict", errs() << "Inv at the caller: " << caller << "\n");

  // 1. propagate from actual to formal parameters
  CRAB_LOG("inter-restrict", errs()
                                 << "Unifying formal and actual parameters\n";);

  for (unsigned i = 0, e = callsite.get_caller_in_params().size(); i < e; ++i) {
    const variable_t &formal = callsite.get_callee_in_params()[i];
    const variable_t &actual = callsite.get_caller_in_params()[i];
    if (!(formal == actual)) {
      CRAB_LOG("inter-restrict",
               errs() << "\t" << formal << ":" << formal.get_type() << " and "
                      << actual << ":" << actual.get_type() << "\n";);
      inter_transformers_impl::unify(caller, formal, actual);
      if (::crab::CrabSanityCheckFlag) {
        if (caller.is_bottom()) {
          CRAB_ERROR("Obtained bottom after unification");
        }
      }
    }
  }
  CRAB_LOG("inter-restrict", errs() << "Inv after formal/actual unification: "
                                    << caller << "\n";);
  // 2. Meet
  callee_at_entry &= caller;
  CRAB_LOG("inter-restrict", errs() << "Inv after meet with callee  "
                                    << callee_at_entry << "\n";);

  // 3. Project onto **input** formal parameters
  callee_at_entry.project(callsite.get_callee_in_params());
  CRAB_LOG(
      "inter-restrict",
      errs() << "Inv at the callee after projecting onto formals: ";
      for (auto &v
           : callsite.get_callee_in_params()) { errs() << v << ";"; } errs()
      << "\n"
      << callee_at_entry << "\n";);
}

template <class Domain>
void inter_abstract_operations<Domain, true>::caller_continuation(
    const callsite_info<typename Domain::variable_t> &callsite,
    const Domain &_callee_at_exit, Domain &caller) {

  using variable_t = typename Domain::variable_t;

  crab::CrabStats::count(_callee_at_exit.domain_name() +
                         ".caller_continuation");
  crab::ScopedCrabStats __st__(_callee_at_exit.domain_name() +
                               ".caller_continuation");

  if (caller.is_bottom()) {
    return;
  }

  if (_callee_at_exit.is_bottom()) {
    caller.set_to_bottom();
    return;
  }

  Domain callee_at_exit(_callee_at_exit);

  // 1. make sure output parameters at the callsite are unconstrained
  caller.forget(callsite.get_caller_out_params());

  CRAB_LOG("inter-extend",
           crab::outs() << "Caller after forgetting lhs variables=" << caller
                        << "\n";);

  // 2. Wire-up outputs: propagate from callee's outputs to caller's
  // lhs of the callsite
  for (unsigned i = 0, e = callsite.get_callee_out_params().size(); i < e;
       ++i) {
    const variable_t &out_formal = callsite.get_callee_out_params()[i];
    const variable_t &out_actual = callsite.get_caller_out_params()[i];
    if (!(out_formal == out_actual)) {
      CRAB_LOG("inter-extend", crab::outs() << "Unifying output " << out_actual
                                            << ":= " << out_formal << "\n";);
      inter_transformers_impl::unify(callee_at_exit, out_actual, out_formal);
    }
  }

  // 3. Wire-up inputs (propagate from callee's inputs to caller's
  // inputs at callsite) and remove the callee variables from the
  // caller continuation. This is needed to propagate up new
  // input-output relationships and also although an input variable
  // cannot be re-assigned its value can be further constrained via
  // assume's.
  //
  // This step is a bit tricky because of two things we need to consider:
  // 1) We cannot forget a callee variable if it appears on the callsite
  // 2) We cannot propagate up if a callsite input parameter is
  //    killed at the callsite (i.e., re-defined via callsite output)
  std::set<variable_t> cs_in_args(callsite.get_caller_in_params().begin(),
                                  callsite.get_caller_in_params().end());
  std::vector<variable_t> killed_cs_in_args =
      inter_transformers_impl::set_intersection(
          callsite.get_caller_in_params(), callsite.get_caller_out_params());

  // Parameters that appear both as callsite argument and callee's
  // formal parameter so they shouldn't be forgotten.
  std::vector<variable_t> caller_and_callee_params;
  caller_and_callee_params.reserve(callsite.get_callee_in_params().size());
  for (unsigned i = 0, e = callsite.get_callee_in_params().size(); i < e; ++i) {
    const variable_t &in_formal = callsite.get_callee_in_params()[i];
    if (cs_in_args.count(in_formal) > 0) {
      caller_and_callee_params.push_back(in_formal);
    } else {
      const variable_t &in_actual = callsite.get_caller_in_params()[i];
      auto lower = std::lower_bound(killed_cs_in_args.begin(),
                                    killed_cs_in_args.end(), in_actual);
      if (lower == killed_cs_in_args.end() || in_actual < *lower) {
        // not found
        //
        // the formal parameter will be forgotten so we need to
        // unify it with the corresponding actual parameter at the
        // callsite to propagate up new relationships created in the
        // callee.
        //
        // Note that this can destroy relationships of the actual
        // parameter at the caller before the call but the meet will
        // restore them later.
        //
        CRAB_LOG("inter-extend", crab::outs() << "Unifying input " << in_actual
                                              << ":=" << in_formal << "\n";);
        inter_transformers_impl::unify(callee_at_exit, in_actual, in_formal);
      }
    }
  }

  std::vector<variable_t> callee_params;
  callsite.get_all_callee_params(callee_params);
  caller_and_callee_params.insert(caller_and_callee_params.end(),
                                  callsite.get_caller_out_params().begin(),
                                  callsite.get_caller_out_params().end());
  callee_params = inter_transformers_impl::set_difference(
      callee_params, caller_and_callee_params);

  // 4. Forget callee's parameters
  callee_at_exit.forget(callee_params);

  CRAB_LOG(
      "inter-extend", crab::outs() << "Forgotten all callee parameters {";
      for (auto const &v
           : callee_params) { crab::outs() << v << ";"; } crab::outs()
      << "}\n";);

  CRAB_LOG("inter-extend2", crab::outs()
                                << "Meet caller with callee:\n"
                                << "CALLER:" << caller << "\n"
                                << "CALLEE:" << callee_at_exit << "\n";);

  // 5. Meet with the caller
  caller &= callee_at_exit;

  CRAB_LOG("inter-extend", crab::outs() << caller << "\n";);
}

} // end namespace domains
} // end namespace crab
