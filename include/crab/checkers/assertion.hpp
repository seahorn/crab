#pragma once

/**
 *  User-definable assertion checker
 **/

#include <crab/checkers/base_property.hpp>

namespace crab {
namespace checker {

template <typename Analyzer>
class assert_property_checker : public property_checker<Analyzer> {

  using varname_t = typename Analyzer::varname_t;
  using number_t = typename Analyzer::number_t;

  using interval_t = ikos::interval<number_t>;
  using abs_dom_t = typename Analyzer::abs_dom_t;
  using base_checker_t = property_checker<Analyzer>;

  using typename base_checker_t::assert_ref_t;
  using typename base_checker_t::assert_t;
  using typename base_checker_t::assign_t;
  using typename base_checker_t::assume_t;
  using typename base_checker_t::basic_block_t;
  using typename base_checker_t::bin_op_t;
  using typename base_checker_t::bool_assert_t;
  using typename base_checker_t::lin_cst_sys_t;
  using typename base_checker_t::lin_cst_t;
  using typename base_checker_t::lin_exp_t;
  using typename base_checker_t::var_t;

  bool entails(const abs_dom_t &lhs, const lin_cst_t &rhs) const {
    CRAB_LOG("checker-entailment", 
	     crab::outs() << "Checking whether\n"
	                  << lhs << "\nentails " << rhs << "\n";);

     bool res = lhs.entails(rhs);

     CRAB_LOG("checker-entailment",
              if (res) { crab::outs() << "\t**entailment holds.\n"; } else {
                crab::outs() << "\t**entailment does not hold.\n";
              });
     return res;
  }
  
public:
  using analyzer_t = Analyzer;

  assert_property_checker(int verbose = 0) : base_checker_t(verbose) {}

  virtual std::string get_property_name() const override {
    return "user-defined assertion checker";
  }

  virtual bool is_interesting(const basic_block_t &bb) const override {
    for (auto &s : bb) {
      if (s.is_assert() || s.is_bool_assert() || s.is_ref_assert()) {
        return true;
      }
    }
    return false;
  }

  virtual void check(assert_t &s) override {
    if (!this->m_abs_tr)
      return;

    const lin_cst_t &cst = s.constraint();

    if (this->m_safe_assertions.count(&s) > 0) {
      crab::crab_string_os os;
      if (this->m_verbose >= 3) {
        os << "Property : " << cst << "\n";
        auto &inv = this->m_abs_tr->get_abs_value();
        os << "Invariant: " << inv << "\n";
        os << "Note: it was proven by the forward+backward analysis";
      }
      this->add_safe(os.str(), &s);
    } else {
      if (cst.is_contradiction()) {
        if (this->m_abs_tr->get_abs_value().is_bottom()) {
          crab::crab_string_os os;
          if (this->m_verbose >= 3) {
            os << "Property : " << cst << "\n";
            auto &inv = this->m_abs_tr->get_abs_value();
            os << "Invariant: " << inv;
          }
          this->add_safe(os.str(), &s);
        } else {
          crab::crab_string_os os;
          if (this->m_verbose >= 2) {
            os << "Property : " << cst << "\n";
            auto &inv = this->m_abs_tr->get_abs_value();
            os << "Invariant: " << inv;
          }
          this->add_warning(os.str(), &s);
        }
        return;
      }

      if (this->m_abs_tr->get_abs_value().is_bottom()) {
        this->m_db.add(_UNREACH, s.get_debug_info());
        return;
      }

      const abs_dom_t &inv = this->m_abs_tr->get_abs_value();
      if (entails(inv, cst)) {
        crab::crab_string_os os;
        if (this->m_verbose >= 3) {
          os << "Property : " << cst << "\n";
          os << "Invariant: " << inv;
        }
        this->add_safe(os.str(), &s);
      } else {
        crab::crab_string_os os;
        if (this->m_verbose >= 2) {
          os << "Property : " << cst << "\n";
          os << "Invariant: " << inv;
        }
        this->add_warning(os.str(), &s);
      }
      // else if (crab::domains::checker_domain_traits<abs_dom_t>::intersect(
      //                inv, cst)) {
      //   crab::crab_string_os os;
      //   if (this->m_verbose >= 2) {
      //     os << "Property : " << cst << "\n";
      //     os << "Invariant: " << inv;
      //   }
      //   this->add_warning(os.str(), &s);
      // } else {
      //   /* Instead this program:
      //      x:=0;
      //      y:=1;
      //      if (x=34) {
      //        assert(y==2);
      //      }
      //      Suppose due to some abstraction we have:
      //       havoc(x);
      //       y:=1;
      //       if (x=34) {
      //         assert(y==2);
      //       }
      //      As a result, we have inv={y=1,x=34}  and cst={y=2}
      //      Note that inv does not either entail or intersect with cst.
      //      However, the original program does not violate the assertion.
      //   */
      //   crab::crab_string_os os;
      //   if (this->m_verbose >= 2) {
      //     os << "Property : " << cst << "\n";
      //     os << "Invariant: " << inv;
      //   }
      //   this->add_warning(os.str(), &s);
      // }
    }
    s.accept(&*this->m_abs_tr); // propagate invariants to the next stmt
  }

  virtual void check(bool_assert_t &s) override {
    if (!this->m_abs_tr) {
      return;
    }

    if (this->m_safe_assertions.count(&s) > 0) {
      crab::crab_string_os os;
      if (this->m_verbose >= 3) {
        os << "Property : " << s << "\n";
        os << "Invariant: " << this->m_abs_tr->get_abs_value() << "\n";
        os << "Note: it was proven by the forward+backward analysis";
      }
      this->add_safe(os.str(), &s);
    } else {

      if (this->m_abs_tr->get_abs_value().is_bottom()) {
        this->m_db.add(_UNREACH, s.get_debug_info());
        return;
      }

      abs_dom_t inv1(this->m_abs_tr->get_abs_value());
      inv1.assume_bool(s.cond(), true /*is_negated*/);
      if (inv1.is_bottom()) {
        crab::crab_string_os os;
        if (this->m_verbose >= 3) {
          os << "Property : " << s << "\n";
          os << "Invariant: " << this->m_abs_tr->get_abs_value();
        }
        this->add_safe(os.str(), &s);
      } else {
        crab::crab_string_os os;
        if (this->m_verbose >= 2) {
          os << "Property : " << s << "\n";
          os << "Invariant: " << this->m_abs_tr->get_abs_value();
        }
        this->add_warning(os.str(), &s);
      }
    }
    s.accept(&*this->m_abs_tr); // propagate invariants to the next stmt
  }

  virtual void check(assert_ref_t &s) override {
    if (!this->m_abs_tr) {
      return;
    }

    if (this->m_safe_assertions.count(&s) > 0) {
      crab::crab_string_os os;
      if (this->m_verbose >= 3) {
        os << "Property : " << s << "\n";
        os << "Invariant: " << this->m_abs_tr->get_abs_value() << "\n";
        os << "Note: it was proven by the forward+backward analysis";
      }
      this->add_safe(os.str(), &s);
    } else {

      if (this->m_abs_tr->get_abs_value().is_bottom()) {
        this->m_db.add(_UNREACH, s.get_debug_info());
        return;
      }

      abs_dom_t inv1(this->m_abs_tr->get_abs_value());
      inv1.ref_assume(s.constraint().negate());
      if (inv1.is_bottom()) {
        crab::crab_string_os os;
        if (this->m_verbose >= 3) {
          os << "Property : " << s << "\n";
          os << "Invariant: " << this->m_abs_tr->get_abs_value();
        }
        this->add_safe(os.str(), &s);
      } else {
        crab::crab_string_os os;
        if (this->m_verbose >= 2) {
          os << "Property : " << s << "\n";
          os << "Invariant: " << this->m_abs_tr->get_abs_value();
        }
        this->add_warning(os.str(), &s);
      }
    }
    s.accept(&*this->m_abs_tr); // propagate invariants to the next stmt
  }
};
} // namespace checker
} // namespace crab
