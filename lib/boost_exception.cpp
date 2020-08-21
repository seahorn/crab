#include <cstdlib>
#include <exception>
#include <boost/assert/source_location.hpp>
#ifdef BOOST_NO_EXCEPTIONS
namespace boost {
void throw_exception(std::exception const &e) {
  // TBD: print error message
  std::exit(1);
}

// Starting with boost 1.73  
void throw_exception(std::exception const &e, boost::source_location const &loc) {
  // TBD: print error message
  std::exit(1);
}
  
} // namespace boost
#endif
