#include <crab/common/types.hpp>
#include <iostream>

namespace crab 
{
  crab_os& outs () { return *crab_os::cout();}
  crab_os& errs () { return *crab_os::cerr();}
}

namespace crab {

  crab_os* crab_os::m_cout = nullptr;
  crab_os* crab_os::m_cerr = nullptr;

  crab_os::crab_os(std::ostream& os) {
    m_os = &os;
  }

  crab_os* crab_os::cout () {
    if (!m_cout) m_cout = new crab_os(std::cout);
    return m_cout;
  }

  crab_os* crab_os::cerr () {
    if (!m_cerr) m_cerr = new crab_os(std::cerr);
    return m_cerr;
  }

  crab_os::~crab_os() { 
    if (m_cout) delete m_cout;
    if (m_cerr) delete m_cerr;
  }

  crab_os & crab_os::operator<<(char C) {
    *m_os << C;
    return *this;
  }

  crab_os & crab_os::operator<<(unsigned char C) {
    *m_os << C;
    return *this;
  }

  crab_os & crab_os::operator<<(signed char C) {
    *m_os << C;
    return *this;
  }

  crab_os & crab_os::operator<<(const char* C) {
    *m_os << C;
    return *this;
  }

  crab_os & crab_os::operator<<(const std::string& Str) {
    *m_os << Str;
    return *this;
  }

  crab_os & crab_os::operator<<(unsigned long N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(long N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(unsigned long long N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(long long N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(const void *P) {
    *m_os << P;
    return *this;
  }

  crab_os & crab_os::operator<<(unsigned int N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(int N) {
    *m_os << N;
    return *this;
  }

  crab_os & crab_os::operator<<(double N) {
    *m_os << N;
    return *this;
  }

} // end namespace
