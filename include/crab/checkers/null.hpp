#pragma once

/*
   Property checker for null-dereference
 */

#include <crab/checkers/base_property.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/nullity.hpp>

namespace crab {

namespace checker {

namespace null_detail {

// arbitrary abstract domain without nullity information
template <typename Domain> struct get_as {
  typedef typename Domain::variable_t variable_t;
  typedef crab::domains::nullity_value nullity_value_t;
  get_as(Domain &) {}
  nullity_value_t operator[](variable_t v) { return nullity_value_t::top(); }
};

// nullity domain
template <typename Number, typename VariableName>
class get_as<crab::domains::nullity_domain<Number, VariableName>> {
  typedef crab::domains::nullity_value nullity_value_t;
  typedef crab::domains::nullity_domain<Number, VariableName> nullity_domain_t;
  typedef typename nullity_domain_t::variable_t variable_t;

  nullity_domain_t &m_inv;

public:
  get_as(nullity_domain_t &inv) : m_inv(inv) {}

  nullity_value_t operator[](variable_t v) { return m_inv.get_nullity(v); }
};

// Reduced product of an arbitrary abstract domain with nullity
template <typename Dom, typename Number, typename VariableName>
class get_as<crab::domains::domain_product2<
    Number, VariableName, Dom,
    crab::domains::nullity_domain<Number, VariableName>>> {

  typedef crab::domains::nullity_value nullity_value_t;
  typedef crab::domains::nullity_domain<Number, VariableName> nullity_domain_t;
  typedef crab::domains::domain_product2<Number, VariableName, Dom,
                                         nullity_domain_t>
      domain_product2_t;
  typedef typename nullity_domain_t::variable_t variable_t;

  domain_product2_t &m_inv;

public:
  get_as(domain_product2_t &inv) : m_inv(inv) {}

  nullity_value_t operator[](variable_t v) {
    return m_inv.second().get_nullity(v);
  }
};

// Reduced product of an arbitrary abstract domain with nullity
template <typename Dom>
class get_as<crab::domains::numerical_nullity_domain<Dom>> {

  typedef crab::domains::numerical_nullity_domain<Dom> domain_t;

  typedef typename domain_t::varname_t varname_t;
  typedef typename domain_t::number_t number_t;
  typedef typename domain_t::variable_t variable_t;

  typedef crab::domains::nullity_value nullity_value_t;
  typedef crab::domains::nullity_domain<number_t, varname_t> nullity_domain_t;

  domain_t &m_inv;

public:
  get_as(domain_t &inv) : m_inv(inv) {}

  nullity_value_t operator[](variable_t v) {
    return m_inv.second().get_nullity(v);
  }
};

} // namespace null_detail

template <typename Analyzer>
class null_property_checker : public property_checker<Analyzer> {

  typedef typename Analyzer::abs_dom_t abs_dom_t;
  typedef typename abs_dom_t::varname_t varname_t;
  typedef typename abs_dom_t::number_t number_t;
  typedef crab::domains::nullity_domain<number_t, varname_t> nullity_domain_t;
  typedef property_checker<Analyzer> base_checker_t;
  using typename base_checker_t::lin_cst_sys_t;
  using typename base_checker_t::lin_cst_t;
  using typename base_checker_t::lin_exp_t;
  using typename base_checker_t::ptr_load_t;
  using typename base_checker_t::ptr_store_t;
  using typename base_checker_t::var_t;

  std::string checked_prop_str(var_t p) {
    std::string res = p.name().str() + " != 0";
    return res;
  }

public:
  using analyzer_t = Analyzer;

  null_property_checker(int verbose = 0) : base_checker_t(verbose) {}

  std::string get_property_name() const override {
    return "null-dereference checker";
  }

  void check(ptr_store_t &s) override {
    if (!this->m_abs_tr)
      return;

    auto &inv = this->m_abs_tr->get_abs_value();
    auto ptr = s.lhs();
    null_detail::get_as<abs_dom_t> null_inv(inv);
    crab::domains::nullity_value val = null_inv[ptr];

    if (val.is_bottom()) {
      this->m_db.add(_UNREACH);
    } else if (val.is_non_null()) {
      crab::crab_string_os os;
      if (this->m_verbose >= 3) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_safe(os.str(), &s);
    } else if (val.is_null()) {
      crab::crab_string_os os;
      if (this->m_verbose >= 1) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_error(os.str(), &s);
    } else {
      crab::crab_string_os os;
      if (this->m_verbose >= 2) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_warning(os.str(), &s);
    }

    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  void check(ptr_load_t &s) override {
    if (!this->m_abs_tr)
      return;

    auto &inv = this->m_abs_tr->get_abs_value();
    auto ptr = s.rhs();
    null_detail::get_as<abs_dom_t> null_inv(inv);
    crab::domains::nullity_value val = null_inv[ptr];

    if (val.is_bottom()) {
      this->m_db.add(_UNREACH);
    } else if (val.is_non_null()) {
      crab::crab_string_os os;
      if (this->m_verbose >= 3) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_safe(os.str(), &s);
    } else if (val.is_null()) {
      crab::crab_string_os os;
      if (this->m_verbose >= 1) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_error(os.str(), &s);
    } else {
      crab::crab_string_os os;
      if (this->m_verbose >= 2) {
        os << "Property : " << checked_prop_str(ptr) << "\n";
        os << "Invariant: " << inv;
      }
      this->add_warning(os.str(), &s);
    }
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }
};
} // namespace checker
} // namespace crab
