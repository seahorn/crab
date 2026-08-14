#pragma once

/*
 * Ghost variable manager for the template DBM (tDBM) domain.
 *
 * A ghost variable ghost(v,k) is the auxiliary dimension representing the
 * scaled term k*v.  Its varname comes from the variable factory's cached
 * string API (get_or_insert_varname(var, suffix)), so identity — the same
 * (v,k) always yields the same varname — and lifetime are owned by the
 * factory.  The manager is stateless; it fixes the naming convention:
 *
 *     name(ghost(v,k)) = v.str() + "*" + k        e.g.  "x*3"  (= 3*x)
 */

#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>

#include <string>

namespace crab {
namespace domains {
namespace tvpi_utils {

template <typename Variable> class ghost_variable_manager {
public:
  using variable_t = Variable;
  using varname_t = typename variable_t::varname_t;
  using number_t = ikos::z_number;

  /// Return (creating if necessary) ghost(v,k).  Identity is cached by the
  /// variable factory.  @p k must be > 1 (k == 1 is the identity: callers
  /// return @p v themselves).
  static variable_t get_or_insert(const variable_t &v, const number_t &k) {
    assert(k > 1);
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
    return variable_t(vfac.get_or_insert_varname(v.name(), encode(k)),
                      v.get_type());
  }

private:
  /// Marker separating the original name from the coefficient.
  static constexpr char marker = '*';

  /// Suffix encoding of coefficient @p k, e.g. "*3".
  static std::string encode(const number_t &k) {
    std::string ks = k.get_str();
    std::string res;
    res.reserve(ks.size() + 1);
    res.push_back(marker);
    res.append(ks);
    return res;
  }
};

} // namespace tvpi_utils
} // end namespace domains
} // end namespace crab
