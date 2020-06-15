#pragma once

#include <crab/common/debug.hpp>
#include <crab/common/os.hpp>
#include <crab/common/types.hpp>

#include <boost/optional.hpp>

namespace crab {

// Very simple language for expressing constraints over reference variables.
template <typename Number, typename VariableName>
class reference_constraint {
public:
  using number_t = Number;
  using variable_t = ikos::variable<number_t, VariableName>;
  using reference_constraint_t = reference_constraint<number_t, VariableName> ;

private:
  using opt_var_t = boost::optional<variable_t>;
  typedef enum { REF_EQ,
		 REF_LESS,
		 REF_LESS_OR_EQ,
		 REF_DISEQ } cst_kind_t;

  
  opt_var_t m_lhs;  // if !m_lhs then NULL
  opt_var_t m_rhs;  // if !m_rhs then NULL
  number_t m_offset; 
  cst_kind_t m_kind;

  reference_constraint(opt_var_t lhs, opt_var_t rhs, cst_kind_t kind)
    : m_lhs(lhs), m_rhs(rhs), m_offset(0), m_kind(kind) {
    if (!m_lhs && m_rhs) { // normalize
      std::swap(m_lhs, m_rhs);
    }
  }
  
  reference_constraint(opt_var_t lhs, opt_var_t rhs, number_t offset, cst_kind_t kind)
    : m_lhs(lhs), m_rhs(rhs), m_offset(offset), m_kind(kind) {
    if (!m_lhs && m_rhs) { // normalize
      std::swap(m_lhs, m_rhs);
    }
  }

public:
  
  reference_constraint() : m_kind(REF_EQ) {}

  // return true iff null != null or null < null
  bool is_contradiction() const {
    return (!m_lhs && !m_rhs &&
	    (m_kind == REF_DISEQ || m_kind == REF_LESS));
  }

  // return true iff null == null or null <= null
  bool is_tautology() const {
    return (!m_lhs && !m_rhs &&
	    (m_kind == REF_EQ || m_kind == REF_LESS_OR_EQ)) ;
  }

  bool is_equality() const { return m_kind == REF_EQ; }

  bool is_disequality() const { return m_kind == REF_DISEQ; }

  bool is_less_or_equal_than() const { return m_kind == REF_LESS_OR_EQ; }

  bool is_less_than() const { return m_kind == REF_LESS; }  

  // return true iff  p REL_OP null
  bool is_unary() const { return (m_lhs && !m_rhs); }

  // return true iff p REL_OP q + offset
  bool is_binary() const { return (m_lhs && m_rhs); }

  variable_t lhs() const {
    if (!m_lhs)
      CRAB_ERROR("reference_constraint lhs is null");
    return *m_lhs;
  }

  variable_t rhs() const {
    if (!m_rhs)
      CRAB_ERROR("reference_constraint rhs is null");
    return *m_rhs;
  }

  number_t offset() const {
    return m_offset;
  }

  reference_constraint_t negate() const {
    if (is_contradiction()) {
      return mk_true();
    } else if (is_tautology()) {
      return mk_false();
    } else if (is_unary() && is_equality()) {
      // p == NULL --> p != NULL 
      return mk_not_null(lhs());
    } else if (is_unary() && is_disequality()) {
      // p != NULL --> q == NULL
      return mk_null(lhs());
    } else if (is_binary()) {
      if (is_equality()) {
	// p == q + k --> p != q + k
 	return mk_not_eq(lhs(), rhs(), offset());
      } else if (is_disequality()) {
	// p != q + k --> p = q + k
	return mk_eq(lhs(), rhs(), offset());
      } else if (is_less_or_equal_than()) {
	// p <= q + k --> q < p - k
	return mk_lt(rhs(), lhs(), -(offset()));
      } else if (is_less_than()) {
	// p <= q + k --> q <= p - k
 	return mk_le(rhs(), lhs(), -(offset()));
      }
    }
    CRAB_ERROR("reference_constraints::negate: unsupported case ", *this);
  }
  static reference_constraint_t mk_true() {
    return reference_constraint_t();
  }

  static reference_constraint_t mk_false() {
    return reference_constraint_t(opt_var_t(), opt_var_t(), REF_DISEQ);
  }

  static reference_constraint_t mk_null(variable_t v) {
    if (!v.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_null:", v, " must be a reference");
    }
    return reference_constraint_t(opt_var_t(v), opt_var_t(), REF_EQ);
  }

  static reference_constraint_t mk_not_null(variable_t v) {
    if (!v.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_not_null:", v, " must be a reference");
    }    
    return reference_constraint_t(opt_var_t(v), opt_var_t(), REF_DISEQ);
  }

  static reference_constraint_t mk_eq(variable_t v1, variable_t v2, number_t offset = number_t(0)) {
    if (!v1.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_eq:", v1, " must be a reference");
    }
    if (!v2.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_eq:", v2, " must be a reference");
    }
    return reference_constraint_t(opt_var_t(v1), opt_var_t(v2), offset, REF_EQ);
  }

  static reference_constraint_t mk_not_eq(variable_t v1, variable_t v2, number_t offset = number_t(0)) {
    if (!v1.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_not_eq:", v1, " must be a reference");
    }
    if (!v2.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_not_eq:", v2, " must be a reference");
    }    
    return reference_constraint_t(opt_var_t(v1), opt_var_t(v2), offset, REF_DISEQ);
  }

  static reference_constraint_t mk_lt(variable_t v1, variable_t v2, number_t offset = number_t(0)) {
    if (!v1.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_lt:", v1, " must be a reference");
    }
    if (!v2.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_lt:", v2, " must be a reference");
    }    
    return reference_constraint_t(opt_var_t(v1), opt_var_t(v2), offset, REF_LESS);
  }

  static reference_constraint_t mk_le(variable_t v1, variable_t v2, number_t offset = number_t(0)) {
    if (!v1.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_leq:", v1, " must be a reference");
    }
    if (!v2.is_ref_type()) {
      CRAB_ERROR("reference_constraint::mk_leq:", v2, " must be a reference");
    }    
    return reference_constraint_t(opt_var_t(v1), opt_var_t(v2), offset, REF_LESS_OR_EQ);
  }
  
  void write(crab_os &o) const {
    if (is_contradiction()) {
      o << "false";
    } else if (is_tautology()) {
      o << "true";
    } else {
      assert(m_lhs);
      o << lhs();
      switch (m_kind) {
      case REF_EQ:
        o << " == ";
	break;
      case REF_DISEQ:
	o << " != ";
	break;
      case REF_LESS_OR_EQ:
	o << " <= ";
	break;
      case REF_LESS:
	o << " < ";
	break;
      default:
	// unreachable
	;;;
      }
      if (!m_rhs)
        o << "NULL";
      else
        o << rhs();

      if (m_offset != 0) {
	o << " + " << m_offset;
      }
    }
  }
};

template <typename Number, typename VariableName>
inline crab_os &operator<<(crab_os &o,
			   const reference_constraint<Number, VariableName> &cst) {
  cst.write(o);
  return o;
}
} //end namespace crab
