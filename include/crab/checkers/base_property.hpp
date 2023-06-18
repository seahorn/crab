#pragma once

/**
 *  Base class for a property checker
 **/

#include <crab/cfg/cfg.hpp>
#include <crab/support/debug.hpp>

#include <map>
#include <set>
#include <vector>

namespace crab {

namespace checker {

enum class check_kind {
  CRAB_SAFE,
  CRAB_ERR,
  CRAB_WARN,
  CRAB_UNREACH
};

// Toy database to store invariants
class checks_db {
public:  
  using checks_map_t =
      std::map<crab::cfg::debug_info, std::vector<check_kind>>;
private:
  
  checks_map_t m_db;
  unsigned m_total_safe;
  unsigned m_total_err;
  unsigned m_total_unreach;
  unsigned m_total_warn;

  void insert_db(const crab::cfg::debug_info &di, check_kind check) {
    m_db[di].push_back(check);
  }

  void merge_db(const checks_map_t &o) {
    for (auto const &kv : o) {
      m_db[kv.first].insert(m_db[kv.first].end(), kv.second.begin(),
                            kv.second.end());
    }
  }

public:
  checks_db()
      : m_total_safe(0),
	m_total_err(0),
	m_total_unreach(0),
	m_total_warn(0) {}

  void clear() {
    m_db.clear();
    m_total_safe = 0;
    m_total_err = 0;
    m_total_unreach = 0;
    m_total_warn = 0;
  }

  bool has_checks(const crab::cfg::debug_info &dbg) const {
    return m_db.find(dbg) != m_db.end();
  }

  // precondition: has_checks(dbg) returns true
  const std::vector<check_kind> &
  get_checks(const crab::cfg::debug_info &dbg) const {
    auto it = m_db.find(dbg);
    if (it != m_db.end()) {
      return it->second;
    } else {
      CRAB_ERROR("Cannot find check associated with ", dbg);
    }
  }

  const checks_map_t& get_all_checks() const {
    return m_db;
  }
  
  unsigned get_total_safe() const { return m_total_safe + m_total_unreach; }

  unsigned get_total_warning() const { return m_total_warn; }

  unsigned get_total_error() const { return m_total_err; }

  // add an entry in the database
  void add(check_kind status, crab::cfg::debug_info dbg) {
    switch (status) {
    case check_kind::CRAB_SAFE:
      m_total_safe++;
      break;
    case check_kind::CRAB_ERR:
      m_total_err++;
      break;
    case check_kind::CRAB_UNREACH:
      m_total_unreach++;
      break;
    default:
      m_total_warn++;
    }
    if (dbg.has_debug()) {
      insert_db(dbg, status);
    }
  }

  // merge two databases
  void operator+=(const checks_db &other) {
    merge_db(other.m_db);
    m_total_safe += other.m_total_safe;
    m_total_err += other.m_total_err;
    m_total_warn += other.m_total_warn;
    m_total_unreach += other.m_total_unreach;
  }

  void write(crab_os &o, bool print_content = false) const {
    std::vector<unsigned> cnts = {m_total_safe, m_total_err, m_total_warn,
                                  m_total_unreach};
    unsigned MaxValLen = 0;
    for (auto c : cnts) {
      MaxValLen = std::max(MaxValLen, (unsigned)std::to_string(c).size());
    }

    o << std::string((int)MaxValLen - std::to_string(m_total_safe).size(), ' ')
      << m_total_safe << std::string(2, ' ') << "Number of total safe checks\n";
    o << std::string((int)MaxValLen - std::to_string(m_total_err).size(), ' ')
      << m_total_err << std::string(2, ' ') << "Number of total error checks\n";
    o << std::string((int)MaxValLen - std::to_string(m_total_warn).size(), ' ')
      << m_total_warn << std::string(2, ' ')
      << "Number of total warning checks\n";
    o << std::string((int)MaxValLen - std::to_string(m_total_unreach).size(),
                     ' ')
      << m_total_unreach << std::string(2, ' ')
      << "Number of total unreachable checks\n";

    if (print_content) {
      o << "\nCheck database content:\n";
      unsigned MaxFileLen = 0;
      for (auto const &kv : m_db) {
        MaxFileLen = std::max(MaxFileLen, (unsigned)kv.first.get_file().size());
      }
      for (auto const &kv : m_db) {
        o << kv.first.get_file()
          << std::string((int)MaxFileLen - kv.first.get_file().size(), ' ')
          << std::string(2, ' ') << " line " << kv.first.get_line() << " col "
          << kv.first.get_column() << ":\n"
          << "\t";
        auto const &checks = kv.second;
        for (unsigned i = 0, num_checks = checks.size(); i < num_checks;) {
          switch (checks[i]) {
          case check_kind::CRAB_SAFE:
            o << "safe";
            break;
          case check_kind::CRAB_ERR:
            o << "error";
            break;
          case check_kind::CRAB_WARN:
            o << "warning";
            break;
          case check_kind::CRAB_UNREACH:
            o << "unreachable (safe)";
            break;
          }
          ++i;
          if (i < num_checks) {
            o << " ";
          }
        }
        o << "\n";
      }
    }
  }
};

template <typename Analyzer>
class property_checker
    : public crab::cfg::statement_visitor<
          typename Analyzer::basic_block_label_t, typename Analyzer::number_t,
          typename Analyzer::varname_t> {
public:
  using analyzer_t = Analyzer;
  using abs_tr_t = typename Analyzer::abs_tr_t;
  using varname_t = typename Analyzer::varname_t;
  using number_t = typename Analyzer::number_t;
  using abs_dom_t = typename Analyzer::abs_dom_t;

  using var_t = typename abs_dom_t::variable_t;
  using lin_exp_t = typename abs_dom_t::linear_expression_t;
  using lin_cst_t = typename abs_dom_t::linear_constraint_t;
  using lin_cst_sys_t = typename abs_dom_t::linear_constraint_system_t;

  using cfg_t = typename Analyzer::cfg_t;
  using basic_block_t = typename cfg_t::basic_block_t;
  using basic_block_label_t = typename cfg_t::basic_block_label_t;

  using statement_t =
      crab::cfg::statement<basic_block_label_t, number_t, varname_t>;
  using bin_op_t =
      crab::cfg::binary_op<basic_block_label_t, number_t, varname_t>;
  using assign_t =
      crab::cfg::assignment<basic_block_label_t, number_t, varname_t>;
  using assume_t =
      crab::cfg::assume_stmt<basic_block_label_t, number_t, varname_t>;
  using assert_t =
      crab::cfg::assert_stmt<basic_block_label_t, number_t, varname_t>;
  using int_cast_t =
      crab::cfg::int_cast_stmt<basic_block_label_t, number_t, varname_t>;
  using select_t =
      crab::cfg::select_stmt<basic_block_label_t, number_t, varname_t>;
  using havoc_t =
      crab::cfg::havoc_stmt<basic_block_label_t, number_t, varname_t>;
  using unreach_t =
      crab::cfg::unreachable_stmt<basic_block_label_t, number_t, varname_t>;
  using callsite_t =
      crab::cfg::callsite_stmt<basic_block_label_t, number_t, varname_t>;
  using intrinsic_t =
      crab::cfg::intrinsic_stmt<basic_block_label_t, number_t, varname_t>;
  using arr_init_t =
      crab::cfg::array_init_stmt<basic_block_label_t, number_t, varname_t>;
  using arr_store_t =
      crab::cfg::array_store_stmt<basic_block_label_t, number_t, varname_t>;
  using arr_load_t =
      crab::cfg::array_load_stmt<basic_block_label_t, number_t, varname_t>;
  using arr_assign_t =
      crab::cfg::array_assign_stmt<basic_block_label_t, number_t, varname_t>;
  using region_init_t =
      crab::cfg::region_init_stmt<basic_block_label_t, number_t, varname_t>;
  using region_copy_t =
      crab::cfg::region_copy_stmt<basic_block_label_t, number_t, varname_t>;
  using region_cast_t =
      crab::cfg::region_cast_stmt<basic_block_label_t, number_t, varname_t>;
  using make_ref_t =
      crab::cfg::make_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using remove_ref_t =
    crab::cfg::remove_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using load_from_ref_t =
      crab::cfg::load_from_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using store_to_ref_t =
      crab::cfg::store_to_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using gep_ref_t =
      crab::cfg::gep_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using assume_ref_t =
      crab::cfg::assume_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using assert_ref_t =
      crab::cfg::assert_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using select_ref_t =
      crab::cfg::ref_select_stmt<basic_block_label_t, number_t, varname_t>;
  using int_to_ref_t =
      crab::cfg::int_to_ref_stmt<basic_block_label_t, number_t, varname_t>;
  using ref_to_int_t =
      crab::cfg::ref_to_int_stmt<basic_block_label_t, number_t, varname_t>;
  using bool_bin_op_t =
      crab::cfg::bool_binary_op<basic_block_label_t, number_t, varname_t>;
  using bool_assign_cst_t =
      crab::cfg::bool_assign_cst<basic_block_label_t, number_t, varname_t>;
  using bool_assign_var_t =
      crab::cfg::bool_assign_var<basic_block_label_t, number_t, varname_t>;
  using bool_assume_t =
      crab::cfg::bool_assume_stmt<basic_block_label_t, number_t, varname_t>;
  using bool_assert_t =
      crab::cfg::bool_assert_stmt<basic_block_label_t, number_t, varname_t>;
  using bool_select_t =
      crab::cfg::bool_select_stmt<basic_block_label_t, number_t, varname_t>;

protected:
  // The abstract transformer
  abs_tr_t *m_abs_tr;
  // Known safe assertions before start forward propagation (it can be empty)
  std::set<const statement_t *> m_safe_assertions;
  // Verbosity to print user messages
  int m_verbose;
  // Store debug information about the checks
  checks_db m_db;
  // Statements where checks occur
  std::vector<const statement_t *> m_safe_checks;
  std::vector<const statement_t *> m_warning_checks;
  std::vector<const statement_t *> m_error_checks;

  void add_safe(std::string msg, const statement_t *s) {
    m_db.add(check_kind::CRAB_SAFE, s->get_debug_info());
    m_safe_checks.push_back(s);

    if (m_verbose >= 3) {
      crab::outs() << " --- SAFE --------------------\n";
      if (s->get_debug_info().has_debug()) {
        crab::outs() << s->get_debug_info() << "\n";
      }
      crab::outs() << msg << "\n";
      crab::outs() << " -----------------------------\n";
    }
  }

  void add_warning(std::string msg, const statement_t *s) {
    m_db.add(check_kind::CRAB_WARN, s->get_debug_info());
    m_warning_checks.push_back(s);

    if (m_verbose >= 2) {
      crab::outs() << " --- WARNING -----------------\n";
      if (s->get_debug_info().has_debug()) {
        crab::outs() << s->get_debug_info() << "\n";
      }
      crab::outs() << msg << "\n";
      crab::outs() << " -----------------------------\n";
    }
  }

  void add_error(std::string msg, const statement_t *s) {
    m_db.add(check_kind::CRAB_ERR, s->get_debug_info());
    m_error_checks.push_back(s);

    if (m_verbose >= 1) {
      crab::outs() << " --- ERROR -------------------\n";
      if (s->get_debug_info().has_debug()) {
        crab::outs() << s->get_debug_info() << "\n";
      }
      crab::outs() << msg << "\n";
      crab::outs() << " -----------------------------\n";
    }
  }

  virtual void check(assert_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bin_op_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(assign_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(assume_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(select_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(int_cast_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(havoc_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(unreach_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(callsite_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(intrinsic_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(arr_init_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(arr_assign_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(arr_store_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(arr_load_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(region_init_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(region_copy_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(region_cast_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }
  
  virtual void check(make_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(remove_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }
  
  virtual void check(load_from_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(store_to_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(gep_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(assume_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(assert_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(select_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(int_to_ref_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(ref_to_int_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_assert_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_bin_op_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_assign_cst_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_assign_var_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_assume_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

  virtual void check(bool_select_t &s) {
    if (!this->m_abs_tr)
      return;
    s.accept(&*this->m_abs_tr); // propagate m_inv to the next stmt
  }

public:
  /* Visitor API */
  virtual void visit(bin_op_t &s) override { check(s); }
  virtual void visit(assign_t &s) override { check(s); }
  virtual void visit(assume_t &s) override { check(s); }
  virtual void visit(select_t &s) override { check(s); }
  virtual void visit(assert_t &s) override { check(s); }
  virtual void visit(int_cast_t &s) override { check(s); }
  virtual void visit(havoc_t &s) override { check(s); }
  virtual void visit(unreach_t &s) override { check(s); }
  virtual void visit(callsite_t &s) override { check(s); }
  virtual void visit(intrinsic_t &s) override { check(s); }
  virtual void visit(arr_init_t &s) override { check(s); }
  virtual void visit(arr_assign_t &s) override { check(s); }
  virtual void visit(arr_store_t &s) override { check(s); }
  virtual void visit(arr_load_t &s) override { check(s); }
  virtual void visit(region_init_t &s) override { check(s); }
  virtual void visit(region_copy_t &s) override { check(s); }
  virtual void visit(region_cast_t &s) override { check(s); }  
  virtual void visit(make_ref_t &s) override { check(s); }
  virtual void visit(remove_ref_t &s) override { check(s); }  
  virtual void visit(load_from_ref_t &s) override { check(s); }
  virtual void visit(store_to_ref_t &s) override { check(s); }
  virtual void visit(gep_ref_t &s) override { check(s); }
  virtual void visit(assume_ref_t &s) override { check(s); }
  virtual void visit(assert_ref_t &s) override { check(s); }
  virtual void visit(select_ref_t &s) override { check(s); }
  virtual void visit(int_to_ref_t &s) override { check(s); }
  virtual void visit(ref_to_int_t &s) override { check(s); }
  virtual void visit(bool_bin_op_t &s) override { check(s); }
  virtual void visit(bool_assign_cst_t &s) override { check(s); }
  virtual void visit(bool_assign_var_t &s) override { check(s); }
  virtual void visit(bool_assume_t &s) override { check(s); }
  virtual void visit(bool_select_t &s) override { check(s); }
  virtual void visit(bool_assert_t &s) override { check(s); }

  property_checker(int verbose) : m_abs_tr(nullptr), m_verbose(verbose) {}

  // whether the basic block is of interest for the checker
  virtual bool is_interesting(const basic_block_t &b) const { return true; }

  // set internal state for the checker
  void set(abs_tr_t *abs_tr,
           const std::set<const statement_t *> &safe_assertions) {
    m_abs_tr = abs_tr;
    m_safe_assertions.insert(safe_assertions.begin(), safe_assertions.end());
  }

  const checks_db &get_db() const { return m_db; }

  checks_db &get_db() { return m_db; }

  const std::vector<const statement_t *> &get_safe_checks() const {
    return m_safe_checks;
  }

  const std::vector<const statement_t *> &get_warning_checks() const {
    return m_warning_checks;
  }

  const std::vector<const statement_t *> &get_error_checks() const {
    return m_error_checks;
  }

  virtual std::string get_property_name() const { return "dummy property"; }

  void write(crab_os &o) const {
    o << get_property_name() << "\n";
    m_db.write(o);
  }
};
} // namespace checker
} // namespace crab
