#pragma once

#include <crab/support/debug.hpp>
#include <string>
#include <vector>

namespace crab {
namespace domains {

/**
 * Contains all the information necessary to execute abstractly a callsite.
 **/
template <class Variable> class callsite_info {
  void check_consistency() const;
  const std::string &m_function;
  const std::vector<Variable> &m_caller_in_params;
  const std::vector<Variable> &m_caller_out_params;
  const std::vector<Variable> &m_callee_in_params;
  const std::vector<Variable> &m_callee_out_params;

public:
  callsite_info(const std::string &function,
                const std::vector<Variable> &caller_in_params,
                const std::vector<Variable> &caller_out_params,
                const std::vector<Variable> &callee_in_params,
                const std::vector<Variable> &callee_out_params);

  const std::vector<Variable> &get_caller_in_params() const;
  const std::vector<Variable> &get_caller_out_params() const;
  const std::vector<Variable> &get_callee_in_params() const;
  const std::vector<Variable> &get_callee_out_params() const;
  const std::string &get_function() const;

  void get_all_caller_params(std::vector<Variable> &out) const;
  void get_all_callee_params(std::vector<Variable> &out) const;
};

template <class Variable>
callsite_info<Variable>::callsite_info(
    const std::string &function, const std::vector<Variable> &caller_in_params,
    const std::vector<Variable> &caller_out_params,
    const std::vector<Variable> &callee_in_params,
    const std::vector<Variable> &callee_out_params)
    : m_function(function), m_caller_in_params(caller_in_params),
      m_caller_out_params(caller_out_params),
      m_callee_in_params(callee_in_params),
      m_callee_out_params(callee_out_params) {
  check_consistency();
}

template <class Variable>
const std::vector<Variable> &
callsite_info<Variable>::get_caller_in_params() const {
  return m_caller_in_params;
}
template <class Variable>
const std::vector<Variable> &
callsite_info<Variable>::get_caller_out_params() const {
  return m_caller_out_params;
}
template <class Variable>
const std::vector<Variable> &
callsite_info<Variable>::get_callee_in_params() const {
  return m_callee_in_params;
}
template <class Variable>
const std::vector<Variable> &
callsite_info<Variable>::get_callee_out_params() const {
  return m_callee_out_params;
}
template <class Variable>
const std::string &callsite_info<Variable>::get_function() const {
  return m_function;
}

template <class Variable>
void callsite_info<Variable>::get_all_caller_params(
    std::vector<Variable> &out) const {
  out.reserve(m_caller_in_params.size() + m_caller_out_params.size());
  out.insert(out.end(), m_caller_in_params.begin(), m_caller_in_params.end());
  out.insert(out.end(), m_caller_out_params.begin(), m_caller_out_params.end());
}

template <class Variable>
void callsite_info<Variable>::get_all_callee_params(
    std::vector<Variable> &out) const {
  out.reserve(m_callee_in_params.size() + m_callee_out_params.size());
  out.insert(out.end(), m_callee_in_params.begin(), m_callee_in_params.end());
  out.insert(out.end(), m_callee_out_params.begin(), m_callee_out_params.end());
}

template <class Variable>
void callsite_info<Variable>::check_consistency() const {
  if (m_caller_in_params.size() != m_callee_in_params.size()) {
    CRAB_ERROR("caller and callee have different number of input parameters in "
               "callsite_info");
  }
  if (m_caller_out_params.size() != m_callee_out_params.size()) {
    CRAB_ERROR("caller and callee have different number of output parameters "
               "in callsite_info");
  }
}
} // end namespace domains
} // end namespace crab
