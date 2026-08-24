#include <crab/support/debug.hpp>

#ifndef NCRABLOG
namespace crab {
bool CrabLogFlag = false;
std::set<std::string> CrabLog;

void CrabEnableLog(std::string x) {
  if (x.empty())
    return;
  CrabLogFlag = true;
  CrabLog.insert(x);
}
} // namespace crab

#else
namespace crab {
void CrabEnableLog(std::string x) {}
} // namespace crab
#endif

namespace crab {
unsigned CrabVerbosity = 0;
void CrabEnableVerbosity(unsigned v) { CrabVerbosity = v; }

bool CrabWarningFlag = true;
void CrabEnableWarningMsg(bool v) { CrabWarningFlag = v; }

bool CrabSanityCheckFlag = false;
void CrabEnableSanityChecks(bool v) { CrabSanityCheckFlag = v; }

crab_os &get_msg_stream(bool timestamp) {
  crab::crab_os *result = &crab::outs();
  if (timestamp) {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[%Y-%m-%d.%X] ", &tstruct);
    *result << buf;
  }
  return *result;
}

source_location::source_location()
    : m_filename(""), m_line(0), m_column(0), m_id(0), m_valid(false) {}

source_location::source_location(const std::string &filename, unsigned line,
                                 unsigned column, unsigned id)
    : m_filename(filename), m_line(line), m_column(column), m_id(id),
      m_valid(true) {}

void source_location::write(crab_os &o) const {
  if (!m_valid) {
    o << "<unknown>";
    return;
  }

  o << m_filename;
  if (m_line > 0) {
    o << ":" << m_line;
    if (m_column > 0) {
      o << ":" << m_column;
    }
  }
}
} // namespace crab
