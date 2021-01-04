#pragma once

/*
 * Factories for variable names.
 */

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp>

#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>

#include <functional>
#include <limits>
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

// This variable factory creates a new variable associated to an
// element of type T. It can also create variables that are not
// associated to an element of type T. We call them shadow variables.
//
// The factory uses a counter of type index_t to generate variable
// id's that always increases.
template <class T> class variable_factory {
public:
  class indexed_varname : public indexable {
    template <typename Any> friend class variable_factory;
    boost::optional<T> m_s;
    ikos::index_t m_id;
    // optional string name associated with m_id if m_s is boost::none
    std::shared_ptr<std::string> m_name;
    variable_factory *m_vfac;

    indexed_varname() = delete;
    // first constructor
    indexed_varname(ikos::index_t id, variable_factory *vfac,
                    std::string name = "")
        : m_s(boost::none), m_id(id),
          m_name(name == "" ? nullptr : std::make_shared<std::string>(name)),
          m_vfac(vfac) {}
    // second constructor
    indexed_varname(T s, ikos::index_t id, variable_factory *vfac)
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
    ~indexed_varname() = default;
    indexed_varname(const indexed_varname &is) = default;
    indexed_varname(indexed_varname &&is) = default;
    indexed_varname &operator=(const indexed_varname &is) = default;
    indexed_varname &operator=(indexed_varname &&is) = default;

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

    variable_factory &get_var_factory() { return *m_vfac; }

    bool operator<(const indexed_varname &s) const { return (m_id < s.m_id); }

    bool operator==(const indexed_varname &s) const { return (m_id == s.m_id); }

    virtual void write(crab_os &o) const override { o << str(); }

    friend crab_os &operator<<(crab_os &o, const indexed_varname &s) {
      o << s.str();
      return o;
    }

    friend size_t hash_value(const indexed_varname &s) {
      std::hash<ikos::index_t> hasher;
      return hasher(s.index());
    }
  };

private:
  using variable_factory_t = variable_factory<T>;
  using t_map_t = std::unordered_map<T, indexed_varname>;
  using shadow_map_t = std::unordered_map<ikos::index_t, indexed_varname>;

  ikos::index_t m_next_id;
  t_map_t m_map;
  shadow_map_t m_shadow_map;
  std::vector<indexed_varname> m_shadow_vars;
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
  using varname_t = indexed_varname;
  using var_range =
      boost::iterator_range<typename std::vector<indexed_varname>::iterator>;
  using const_var_range = boost::iterator_range<
      typename std::vector<indexed_varname>::const_iterator>;

public:
  /************************* High-level API *************************/

  variable_factory() : m_next_id(1) {}

  virtual ~variable_factory() = default;

  variable_factory(ikos::index_t start_id) : m_next_id(start_id) {}

  variable_factory(const variable_factory_t &o) = delete;

  variable_factory_t &operator=(const variable_factory_t &o) = delete;

  virtual indexed_varname operator[](T s) {
    auto it = m_map.find(s);
    if (it == m_map.end()) {
      indexed_varname is(s, get_and_increment_id(), this);
      m_map.insert({s, is});
      return is;
    } else {
      return it->second;
    }
  }

  /************************* Low-level API *************************/
  // generate indexed_varname's without being associated with a
  // particular T (w/o caching).
  //
  // XXX: do not use it unless strictly necessary because it can
  // produce an unbounded number of indexed_varname objects. This
  // commonly used by abstract domains to generate synthetic
  // variables.
  virtual indexed_varname get(std::string name = "") {
    indexed_varname is(get_and_increment_id(), this, name);
    m_shadow_vars.push_back(is);
    return is;
  }

  // generate a shadow indexed_varname's associated to some key.  Used
  // also by abstract domains to generate synthetic variables, but in
  // this case, the number is bounded by the number of T;s objects.
  virtual indexed_varname get(ikos::index_t key, std::string name = "") {
    auto it = m_shadow_map.find(key);
    if (it == m_shadow_map.end()) {
      indexed_varname is(get_and_increment_id(), this, name);
      m_shadow_map.insert(typename shadow_map_t::value_type(key, is));
      m_shadow_vars.push_back(is);
      return is;
    } else {
      return it->second;
    }
  }

  // return all the shadow variables created by the factory.
  virtual const_var_range get_shadow_vars() const {
    return boost::make_iterator_range(m_shadow_vars.begin(),
                                      m_shadow_vars.end());
  }

  // When varname_t::str() is called this map is used to do renaming
  void add_renaming_map(
      const std::unordered_map<std::string, std::string> &smap) const {
    clear_renaming_map();
    m_renaming_map.insert(smap.begin(), smap.end());
  }

  void clear_renaming_map() const { m_renaming_map.clear(); }

  const std::unordered_map<std::string, std::string> &get_renaming_map() const {
    return m_renaming_map;
  }
};

//! Specialized factory for strings
class str_variable_factory : public variable_factory<std::string> {
  using variable_factory_t = variable_factory<std::string>;

public:
  using varname_t = variable_factory_t::varname_t;
  using const_var_range = variable_factory_t::const_var_range;

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
