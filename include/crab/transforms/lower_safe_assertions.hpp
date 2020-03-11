#pragma once

#include <crab/cfg/cfg.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/transforms/transform.hpp>

#include <set>

namespace crab {
namespace transforms {

/**
 * Replace safe assertions with assume statements.
 **/
template <typename CFG> class lower_safe_assertions : public transform<CFG> {

  // some typedefs here
  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;

  using statement_t = typename CFG::statement_t;

  using assert_t = crab::cfg::assert_stmt<number_t, varname_t>;
  using ptr_assert_t = crab::cfg::ptr_assert_stmt<number_t, varname_t>;
  using bool_assert_t = crab::cfg::bool_assert_stmt<number_t, varname_t>;

  using assume_t = crab::cfg::assume_stmt<number_t, varname_t>;
  using ptr_assume_t = crab::cfg::ptr_assume_stmt<number_t, varname_t>;
  using bool_assume_t = crab::cfg::bool_assume_stmt<number_t, varname_t>;

  using basic_block_t = typename CFG::basic_block_t;

  bool is_assert(const statement_t *s) const {
    return s->is_assert() || s->is_ptr_assert() || s->is_bool_assert();
  }

public:
  lower_safe_assertions(const std::set<const statement_t *> &safe_checks)
      : transform<CFG>(), m_safe_checks(safe_checks) {
    CRAB_LOG("lsa", crab::outs()
                        << "Found " << safe_checks.size() << " assertions\n";);
  }

  virtual bool run(CFG &cfg) override {
    bool change = false;
    std::vector<std::pair<statement_t *, basic_block_t *>> to_replace;
    for (auto b_it = cfg.begin(), b_et = cfg.end(); b_it != b_et; ++b_it) {
      basic_block_t &bb = *b_it;
      for (auto s_it = bb.begin(), s_et = bb.end(); s_it != s_et; ++s_it) {
        statement_t &s = *s_it;
        if (is_assert(&s) && m_safe_checks.count(&s) > 0) {
          CRAB_LOG("lsa", crab::outs() << "Found assertion to be lowered: " << s
                                       << "\n";);
          to_replace.push_back({&s, &bb});
          change = true;
        }
      }
    }

    while (!to_replace.empty()) {
      auto p = to_replace.back();
      statement_t *s = p.first;
      basic_block_t *parent = p.second;
      to_replace.pop_back();
      if (s->is_assert()) {
        statement_t *new_s =
            new assume_t((static_cast<assert_t *>(s))->constraint());
        CRAB_LOG("lsa", crab::outs() << "Replacing " << *s << " with " << *new_s
                                     << "\n";);
        parent->replace(s, new_s);
      } else if (s->is_ptr_assert()) {
        statement_t *new_s =
            new ptr_assume_t((static_cast<ptr_assert_t *>(s))->constraint());
        CRAB_LOG("lsa", crab::outs() << "Replacing " << *s << " with " << *new_s
                                     << "\n";);
        parent->replace(s, new_s);
      } else if (s->is_bool_assert()) {
        statement_t *new_s = new bool_assume_t(
            (static_cast<bool_assert_t *>(s))->cond(), false /* not negated*/);
        CRAB_LOG("lsa", crab::outs() << "Replacing " << *s << " with " << *new_s
                                     << "\n";);
        parent->replace(s, new_s);
      } else {
        // unreachable
      }
    }
    return change;
  }

  virtual std::string get_name() const override {
    return "Replace safe assertions with assume statements";
  }

private:
  const std::set<const statement_t *> &m_safe_checks;
};

} // namespace transforms
} // namespace crab
