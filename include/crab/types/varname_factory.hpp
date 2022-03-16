#pragma once

/*
 * Factories for variable names. The key feature of a variable name is
 * that it can be indexed in a datastructure such as patricia trees by
 * providing the method index().
 */

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp>

#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

namespace crab {

/* T is the template parameter of var_factory_impl::variable_factory */
template <typename T> class variable_name_traits {
public:
  static std::string to_string(T varname);
};

namespace var_factory_impl {

template <typename T> class variable_factory;

template <class T> class indexed_varname : public indexable {
  template <typename Any> friend class variable_factory;
public:  
  using variable_factory_t = variable_factory<T>;
  
private:
  using this_type = indexed_varname<T>;
  
  boost::optional<T> m_s;
  ikos::index_t m_id;
  // optional string name associated with m_id if m_s is boost::none
  std::shared_ptr<std::string> m_name;
  variable_factory_t *m_vfac;

  indexed_varname() = delete;
  // first constructor
  indexed_varname(ikos::index_t id, variable_factory_t *vfac,
                  std::string name = "")
      : m_s(boost::none), m_id(id),
        m_name(name == "" ? nullptr : std::make_shared<std::string>(name)),
        m_vfac(vfac) {}
  // second constructor
  indexed_varname(T s, ikos::index_t id, variable_factory_t *vfac)
      : m_s(s), m_id(id), m_name(nullptr), m_vfac(vfac) {}

  std::string rename(const std::string &s) const {
    auto it = m_vfac->get_renaming_map().find(s);
    if (it != m_vfac->get_renaming_map().end()) {
      return it->second;
    } else {
      return s;
    }
  }

public:
  indexed_varname(const this_type &is) = default;
  indexed_varname(this_type &&is) = default;
  ~indexed_varname() = default;
  this_type &operator=(const this_type &is) = default;
  this_type &operator=(this_type &&is) = default;

  virtual ikos::index_t index() const override { return m_id; }

  std::string str() const {
    if (m_s) {
      return rename(crab::variable_name_traits<T>::to_string(*m_s));
    } else if (m_name && (*m_name != "")) {
      return rename(*m_name);
    } else {
      // unlikely prefix
      return rename("@V_" + std::to_string(m_id));
    }
  }

  boost::optional<T> get() const { return m_s; }

  variable_factory_t &get_var_factory() { return *m_vfac; }

  bool operator<(const this_type &s) const { return (m_id < s.m_id); }

  bool operator==(const this_type &s) const { return (m_id == s.m_id); }

  virtual void write(crab_os &o) const override { o << str(); }

  size_t hash() const {
    std::hash<ikos::index_t> hasher;
    return hasher(index());
  }
  
  friend crab_os &operator<<(crab_os &o, const this_type &s) {
    s.write(o);
    return o;
  }

};
// used by boost::hash_combine (no std::hash_combine in C+11)  
//template <typename T>
//inline std::size_t hash_value(const crab::var_factory_impl::indexed_varname<T>&v) {
//  return v.hash();
//}  
} // end namespace var_factory_impl
} // end namespace crab

namespace std {
template <typename T> struct hash<crab::var_factory_impl::indexed_varname<T>> {
  size_t operator()(
      const typename crab::var_factory_impl::indexed_varname<T> &v) const {
    return v.hash();
  }
};
} // end namespace std

namespace crab {
namespace var_factory_impl {

// This variable factory (it's actually a factory of
// indexed_varname's) creates a new indexed_variable associated to an
// element of type T if provided. It can also create indexed_varname's
// that are not associated to an element of type T. We call them
// shadow variables.
//
// The factory uses a counter of type index_t to generate variable
// id's that always increases.
template <class T> class variable_factory {
  using variable_factory_t = variable_factory<T>;
  using t_map_t = std::unordered_map<T, indexed_varname<T>>;
  using shadow_map_t =
      std::unordered_map<indexed_varname<T>,
                         std::map<std::string, indexed_varname<T>>>;
  // global counter to generate indexes
  ikos::index_t m_next_id;
  // (cached) indexed_varname's associated with a T-instance.
  t_map_t m_map;
  // (non-cached) fresh indexed_varname's.
  std::vector<indexed_varname<T>> m_shadow_vars;
  // (cached) indexed_varname's associated with another indexed_varname.
  shadow_map_t m_shadow_map;
  mutable std::unordered_map<std::string, std::string> m_renaming_map;

  ikos::index_t get_and_increment_id(void) {
    if (m_next_id == std::numeric_limits<ikos::index_t>::max()) {
      CRAB_ERROR("Reached limit of ", std::numeric_limits<ikos::index_t>::max(),
                 " variables");
    }
    ikos::index_t res = m_next_id;
    ++m_next_id;
    return res;
  }

public:
  using varname_t = indexed_varname<T>;

  variable_factory() : m_next_id(1) {}

  virtual ~variable_factory() = default;

  variable_factory(ikos::index_t start_id) : m_next_id(start_id) {}

  variable_factory(const variable_factory_t &o) = delete;

  variable_factory_t &operator=(const variable_factory_t &o) = delete;

  virtual varname_t operator[](T s) {
    auto it = m_map.find(s);
    if (it == m_map.end()) {
      varname_t iv(s, get_and_increment_id(), this);
      m_map.insert({s, iv});
      return iv;
    } else {
      return it->second;
    }
  }

  // Generate a fresh indexed_varname's without being associated with
  // a particular instance of T.
  //
  // If you are an abstract domain then do not use it unless strictly
  // necessary because it can produce an unbounded number of
  // indexed_varname objects.
  virtual varname_t get(std::string name = "") {
    varname_t iv(get_and_increment_id(), this, name);
    m_shadow_vars.push_back(iv);
    return iv;
  }

  // API for abstract domains
  //
  // Create a fresh indexed_varname associated to var.
  // Given the same var and name it always return the same indexed_varname.
  // The returned indexed_varname's name is var's name concatenated with name.
  virtual varname_t get(const varname_t &var, std::string name) {
    auto it = m_shadow_map.find(var);
    if (it == m_shadow_map.end()) {
      varname_t iv(get_and_increment_id(), this, var.str() + name);
      std::map<std::string, varname_t> named_shadows;
      named_shadows.insert({name, iv});
      m_shadow_map.insert({var, named_shadows});
      return iv;
    } else {
      std::map<std::string, varname_t> &named_shadows = it->second;
      auto nit = named_shadows.find(name);
      if (nit != named_shadows.end()) {
        return nit->second;
      } else {
        varname_t iv(get_and_increment_id(), this, var.str() + name);
        named_shadows.insert({name, iv});
        return iv;
      }
    }
  }

  // Allow temporary renaming for pretty printing
  void add_renaming_map(
      const std::unordered_map<std::string, std::string> &smap) const {
    clear_renaming_map();
    m_renaming_map.insert(smap.begin(), smap.end());
  }

  void clear_renaming_map() const { m_renaming_map.clear(); }

  const std::unordered_map<std::string, std::string> &get_renaming_map() const {
    return m_renaming_map;
  }

  // return all the non-T variables created by the factory.  
  virtual std::vector<varname_t> get_shadow_vars() const {
    std::vector<varname_t> out(m_shadow_vars.begin(), m_shadow_vars.end());
    for (auto &kv_ : m_shadow_map) {
      for (auto &kv : kv_.second) {
        out.push_back(varname_t(kv.second));
      }
    }
    return out;
  }
};

//! Specialized factory for strings
class str_variable_factory : public variable_factory<std::string> {
  using variable_factory_t = variable_factory<std::string>;

public:
  using varname_t = typename variable_factory_t::varname_t;

  str_variable_factory() : variable_factory_t() {}
};

inline int fresh_colour(int col_x, int col_y) {
  switch (col_x) {
  case 0:
    return col_y == 1 ? 2 : 1;
  case 1:
    return col_y == 0 ? 2 : 0;
  case 2:
    return col_y == 0 ? 1 : 0;
  default:
    CRAB_ERROR("Unreachable");
  }
}

//! Three-coloured variable allocation. So the number of variables
//  is bounded by 3|Tbl|, rather than always increasing.
class str_var_alloc_col {
  static const char **col_prefix;

  // Ideally, str_var_alloc_col should inherit from
  // str_variable_factory.  However, we need public copy constructors
  // for str_var_alloc_col. The hack here is to have vfac as a global
  // factory that all str_var_alloc_col's share.
  static str_variable_factory &get_vfac() {
    static str_variable_factory vfac;
    return vfac;
  }

public:
  using varname_t = str_variable_factory::varname_t;

  str_var_alloc_col() : colour(0), next_id(0) {}

  str_var_alloc_col(const str_var_alloc_col &o)
      : colour(o.colour), next_id(o.next_id) {}

  str_var_alloc_col(const str_var_alloc_col &x, const str_var_alloc_col &y)
      : colour(fresh_colour(x.colour, y.colour)), next_id(0) {
    assert(colour != x.colour);
    assert(colour != y.colour);
  }

  str_var_alloc_col &operator=(const str_var_alloc_col &x) {
    if (this != &x) {
      colour = x.colour;
      next_id = x.next_id;
    }
    return *this;
  }

  str_variable_factory::varname_t next() {
    std::string v = col_prefix[colour] + std::to_string(next_id++);
    return get_vfac()[v];
  }

  void add_renaming_map(
      const std::unordered_map<std::string, std::string> &smap) const {
    get_vfac().add_renaming_map(smap);
  }

  void clear_renaming_map() const { get_vfac().clear_renaming_map(); }

protected:
  int colour;
  ikos::index_t next_id;
};

//// The type ikos::index_t cannot be used directly as varname_t. We
//// need a wrapper that implements index() function.
// //! Specialized factory for integers
// class int_variable_factory {
// public:
//   using varname_t = ikos::index_t;
//   int_variable_factory() {}
//   int_variable_factory(const int_variable_factory &o) = delete;
//   int_variable_factory &operator=(const int_variable_factory &o) = delete;
//   varname_t operator[](ikos::index_t v) { return v; }
// };
// class int_var_alloc_col {
// public:
//   using varname_t = ikos::index_t;
//   static int_variable_factory vfac;
//   int_var_alloc_col() : colour(0), next_id(0) {}
//   int_var_alloc_col(const int_var_alloc_col &o)
//       : colour(o.colour), next_id(o.next_id) {}
//   int_var_alloc_col(const int_var_alloc_col &x, const int_var_alloc_col &y)
//       : colour(fresh_colour(x.colour, y.colour)), next_id(0) {
//     assert(colour != x.colour);
//     assert(colour != y.colour);
//   }
//   int_var_alloc_col &operator=(const int_var_alloc_col &x) {
//     colour = x.colour;
//     next_id = x.next_id;
//     return *this;
//   }
//   int_variable_factory::varname_t next() {
//     ikos::index_t id = next_id++;
//     return 3 * id + colour;
//   }
// protected:
//   int colour;
//   ikos::index_t next_id;
// };

} // end namespace var_factory_impl
} // end namespace crab
