#pragma once

namespace crab {
namespace domains {
namespace object_domain_impl {

/* 
  This class wraps the numerical domain used in ODI map
  The purpose of this class is for handling operator == used in Patricia trees.
*/
template <class Domain>
class region_object {
  using domain_t = Domain;
  using variable_t = typename domain_t::variable_t;
  using linear_expression_t = typename domain_t::linear_expression_t;
  using linear_constraint_system_t = typename domain_t::linear_constraint_system_t;
  using interval_t = typename domain_t::interval_t;

  domain_t m_base;
  region_object(domain_t &&base) : m_base(base) {}
public:
  region_object() { m_base.set_to_top(); }
  region_object(const region_object &o) : m_base(o.m_base) {}
  region_object(region_object &&o) : m_base(std::move(o.m_base)) {}
  region_object &operator=(const region_object &o) {
    if (this != &o) {
      m_base = o.m_base;
    }
    return *this;
  }

  region_object &operator=(region_object &&o) {
    if (this != &o) {
      m_base = std::move(o.m_base);
    }
    return *this;
  }

  static region_object bottom() {
    region_object res;
    res.m_base.set_to_bottom();
    return res;
  }

  static region_object top() {
    region_object res;
    return res;
  }

  domain_t &dom() {
    return m_base;
  }

  const domain_t&dom() const {
    return m_base;
  }

  void set_dom(const domain_t &dom) {
    m_base = dom;
  }

  bool is_bottom() const { return m_base.is_bottom(); }
  bool is_top() const { return m_base.is_top(); }

  bool operator<=(const region_object &o) const {
    return m_base <= o.m_base;
  }
  bool operator==(const region_object &o) const {
    return (m_base <= o.m_base && o.m_base <= m_base);
  }
  void operator|=(const region_object &o) {
    m_base |= o.m_base;
  }
  region_object operator|(const region_object &o) const {
    return region_object(m_base | o.m_base);
  }
  region_object operator||(const region_object &o) const {
    return region_object(m_base || o.m_base);
  }
  region_object operator&(const region_object &o) const {
    return region_object(m_base & o.m_base);
  }
  region_object operator&&(const region_object &o) const {
    return region_object(m_base && o.m_base);
  }

  void assign(const variable_t &x, const linear_expression_t &e) {
    m_base.assign(x, e);
  }

  interval_t operator[](const variable_t &v) {
    return m_base[v];
  }

  void operator+=(const linear_constraint_system_t &csts) {
    m_base += csts;
  }

  void operator-=(const variable_t &v) {
    m_base -= v;
  }

  void write(crab::crab_os &o) const {
    o << m_base;
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, region_object &dom) {
    dom.write(o);
    return o;
  }
};
} // end namespace object_domain_impl
} // end namespace domains
} // end namespace crab