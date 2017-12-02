#include <crab/config.h>
#include <crab/common/types.hpp>
#include <crab/common/wrapint.hpp>
#include <crab/domains/wrapped_interval_domain.hpp>

using namespace std;
using namespace crab::domains;

int main(int argc, char *argv[]) {

  typedef wrapped_interval<ikos::z_number> wrapped_interval_t;
  
  wrapped_interval_t i1(crab::wrapint(-10, 8), crab::wrapint(10, 8));
  wrapped_interval_t i2(crab::wrapint(0, 8), crab::wrapint(20, 8));  
  
  crab::outs () << "i1=" << i1 << "\n";
  crab::outs () << "i2=" << i2 << "\n";
  crab::outs () << "-i2=" << (-i2) << "\n";
  crab::outs() << "i1 & i2: " << (i1 & i2) << "\n";
  crab::outs() << "i1 | i2: " << (i1 | i2) << "\n";
  
  wrapped_interval_t i3(crab::wrapint(0,8), crab::wrapint(-10,8));
  wrapped_interval_t i4(crab::wrapint(-15,8),crab::wrapint(5,8));
  crab::outs() << i3 << " & " << i4 << ": " << (i3 & i4) << "\n";
  crab::outs() << i3 << " | " << i4 << ": " << (i3 | i4) << "\n";
  
  wrapped_interval_t urk(crab::wrapint(1642571628,32), crab::wrapint(700177772,32));
  crab::outs() << "urk: " << urk << "\n";
  crab::outs() << "urk.is_top(): " << urk.is_top() << "\n";
  crab::outs() << "urk == Top: " << (urk == wrapped_interval_t::top()) << "\n";
  
  return 0;
} 
